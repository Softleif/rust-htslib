// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::borrow::{Borrow, BorrowMut};
use std::fmt;
use std::marker::PhantomData;
use std::ops::Deref;
use std::os::raw::c_char;
use std::ptr;
use std::slice;
use std::str;
use std::sync::Arc;
use std::{ffi, iter};

use bio_types::genome;
use cstr8::{cstr8, CStr8, CString8};
use derive_new::new;
use ieee754::Ieee754;
use lazy_static::lazy_static;

use crate::bcf::header::{HeaderView, Id};
use crate::bcf::BcfError as Error;

type Result<T> = std::result::Result<T, Error>;
use crate::htslib;

const MISSING_INTEGER: i32 = i32::MIN;
const VECTOR_END_INTEGER: i32 = i32::MIN + 1;

lazy_static! {
    static ref MISSING_FLOAT: f32 = Ieee754::from_bits(0x7F80_0001);
    static ref VECTOR_END_FLOAT: f32 = Ieee754::from_bits(0x7F80_0002);
}

// ---------------------------------------------------------------------------
// Pure Rust BCF binary decoding helpers (replacement for htslib inline fns)
// ---------------------------------------------------------------------------

/// Unofficial BCF type code for 64-bit integers (not exported by hts-sys).
const BCF_BT_INT64: u32 = 4;

/// Maps BCF type codes (0–15) to log2(element size in bytes).
/// Matches the C `bcf_type_shift[16]` table. Unused/invalid codes map to 0.
const BCF_TYPE_SHIFT: [u8; 16] = [
    0, // 0  BCF_BT_NULL  → 1 byte
    0, // 1  BCF_BT_INT8  → 1 byte
    1, // 2  BCF_BT_INT16 → 2 bytes
    2, // 3  BCF_BT_INT32 → 4 bytes
    3, // 4  BCF_BT_INT64 → 8 bytes  (unofficial)
    2, // 5  BCF_BT_FLOAT → 4 bytes
    0, // 6  (unused)
    0, // 7  BCF_BT_CHAR  → 1 byte
    0, 0, 0, 0, 0, 0, 0, 0, // 8–15 unused
];

/// Compute the byte length of `count` elements of BCF type `type_code`.
#[inline]
fn bcf_type_bytes(count: usize, type_code: u32) -> usize {
    count << BCF_TYPE_SHIFT[type_code as usize] as usize
}

/// Find the length of a BCF_BT_CHAR field, stopping at the first NUL byte.
/// Returns at most `max_len` (the declared element count).
#[inline]
fn char_field_len(data: &[u8], max_len: usize) -> usize {
    memchr::memchr(0, &data[..max_len]).unwrap_or(max_len)
}

/// Decode one typed integer from a BCF byte stream.
///
/// Returns `Some((value, bytes_consumed))` on success, or `None` if
/// `type_code` is unrecognised or `data` is too short. Returning `None`
/// prevents callers from advancing by zero bytes (infinite loop on corrupt data).
fn bcf_dec_int1(data: &[u8], type_code: u32) -> Option<(i64, usize)> {
    match type_code {
        htslib::BCF_BT_INT8 if !data.is_empty() => Some((data[0] as i8 as i64, 1)),
        htslib::BCF_BT_INT16 if data.len() >= 2 => {
            Some((i16::from_le_bytes([data[0], data[1]]) as i64, 2))
        }
        htslib::BCF_BT_INT32 if data.len() >= 4 => Some((
            i32::from_le_bytes([data[0], data[1], data[2], data[3]]) as i64,
            4,
        )),
        BCF_BT_INT64 if data.len() >= 8 => Some((
            i64::from_le_bytes([
                data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7],
            ]),
            8,
        )),
        _ => None,
    }
}

/// Decode a typed integer where the first byte encodes the type.
/// Returns `Some((value, bytes_consumed))` including the type byte,
/// or `None` on corrupt data.
fn bcf_dec_typed_int1(data: &[u8]) -> Option<(i64, usize)> {
    let type_code = (*data.first()? & 0xf) as u32;
    let (val, consumed) = bcf_dec_int1(data.get(1..)?, type_code)?;
    Some((val, 1 + consumed))
}

/// Decode the size+type header byte(s) of a BCF typed value.
/// Returns `Some((count, type_code, bytes_consumed))`, or `None` on corrupt data.
fn bcf_dec_size(data: &[u8]) -> Option<(usize, u32, usize)> {
    let first = *data.first()?;
    let type_code = (first & 0xf) as u32;
    let count_nibble = first >> 4;
    if count_nibble != 15 {
        Some((count_nibble as usize, type_code, 1))
    } else {
        let (val, consumed) = bcf_dec_typed_int1(data.get(1..)?)?;
        Some((val as usize, type_code, 1 + consumed))
    }
}

/// Grow a C-allocated buffer using `libc::realloc`. Aborts on allocation failure.
///
/// # Safety
/// `ptr` must be null or point to memory previously allocated by `libc::malloc`/`realloc`.
#[inline]
unsafe fn crealloc<T>(ptr: *mut T, count: usize) -> *mut T {
    let new = libc::realloc(ptr as *mut libc::c_void, count * std::mem::size_of::<T>());
    if new.is_null() && count > 0 {
        std::process::abort();
    }
    new as *mut T
}

/// RAII wrapper for a C-compatible `kbitset_t`, replacing `kbs_init`/`kbs_insert`/`kbs_destroy`.
///
/// Allocates via `libc::malloc` so the pointer can be passed to C functions like
/// `bcf_remove_allele_set`. Freed on drop via `libc::free`.
struct KBitSet {
    ptr: *mut htslib::kbitset_t,
}

impl KBitSet {
    /// Create a bitset from a boolean slice. `bits[i] == true` means bit `i` is set.
    fn from_bools(bits: &[bool]) -> Self {
        let ni = bits.len();
        let elt_bits = std::mem::size_of::<libc::c_ulong>() * 8;
        // Match C: n = (ni + KBS_ELTBITS - 1) / KBS_ELTBITS  (ceiling division)
        let n_slots = ni.div_ceil(elt_bits);
        // Alloc: sizeof(kbitset_t) already includes b[1], plus n_slots extra ulongs
        let alloc_size = std::mem::size_of::<htslib::kbitset_t>()
            + n_slots * std::mem::size_of::<libc::c_ulong>();
        let ptr = unsafe { libc::malloc(alloc_size) as *mut htslib::kbitset_t };
        if ptr.is_null() {
            std::process::abort();
        }
        unsafe {
            // Match C: bs->n = bs->n_max = n  (both store the slot count)
            (*ptr).n = n_slots;
            (*ptr).n_max = n_slots;
            // Zero all bit slots (n_slots data + 1 sentinel)
            std::ptr::write_bytes((*ptr).b.as_mut_ptr(), 0, n_slots + 1);
            // Sentinel at b[n_slots]: matches kbs_last_mask = KBS_MASK(ni) - 1
            let last_mask = {
                let m = (1usize << (ni % elt_bits)) as libc::c_ulong;
                let mask = m.wrapping_sub(1);
                if mask == 0 {
                    !0 as libc::c_ulong
                } else {
                    mask
                }
            };
            *(*ptr).b.as_mut_ptr().add(n_slots) = last_mask;
            // Set individual bits
            for (i, &set) in bits.iter().enumerate() {
                if set {
                    let elt = i / elt_bits;
                    let mask = (1usize << (i % elt_bits)) as libc::c_ulong;
                    *(*ptr).b.as_mut_ptr().add(elt) |= mask;
                }
            }
        }
        KBitSet { ptr }
    }

    fn as_ptr(&self) -> *const htslib::kbitset_t {
        self.ptr
    }
}

impl Drop for KBitSet {
    fn drop(&mut self) {
        unsafe { libc::free(self.ptr as *mut libc::c_void) };
    }
}

/// Pure Rust replacement for `htslib::bcf_unpack`.
///
/// Parses the binary BCF data in `b->shared` and `b->indiv` and populates
/// the decoded fields in `b->d`. Uses `libc::realloc` for C-compatible
/// memory management (the buffers will be freed by `bcf_destroy`).
///
/// # Safety
/// `b` must be a valid, non-null pointer to a `bcf1_t` whose `shared` and
/// `indiv` kstrings contain valid BCF binary data.
pub(crate) unsafe fn bcf_unpack_rs(b: *mut htslib::bcf1_t, which: i32) -> i32 {
    if (*b).shared.l == 0 {
        return 0;
    }

    let mut which = which;
    if which & htslib::BCF_UN_FLT as i32 != 0 {
        which |= htslib::BCF_UN_STR as i32;
    }
    if which & htslib::BCF_UN_INFO as i32 != 0 {
        which |= htslib::BCF_UN_SHR as i32;
    }

    let shared = slice::from_raw_parts((*b).shared.s as *const u8, (*b).shared.l);
    let d = std::ptr::addr_of_mut!((*b).d);

    // Dispatch to per-phase helpers. Each returns 0 on success, -1 on corrupt data.

    if (which & htslib::BCF_UN_STR as i32 != 0)
        && ((*b).unpacked & htslib::BCF_UN_STR as i32 == 0)
        && unpack_str(shared, b, d) < 0
    {
        return -1;
    }

    if (which & htslib::BCF_UN_FLT as i32 != 0)
        && ((*b).unpacked & htslib::BCF_UN_FLT as i32 == 0)
        && unpack_flt(shared, b, d) < 0
    {
        return -1;
    }

    if (which & htslib::BCF_UN_INFO as i32 != 0)
        && ((*b).unpacked & htslib::BCF_UN_INFO as i32 == 0)
        && unpack_info(shared, b, d) < 0
    {
        return -1;
    }

    if (which & htslib::BCF_UN_FMT as i32 != 0)
        && (*b).n_sample() > 0
        && ((*b).unpacked & htslib::BCF_UN_FMT as i32 == 0)
        && unpack_fmt(b, d) < 0
    {
        return -1;
    }

    0
}

/// Write a BCF_BT_CHAR field (or "." for count==0) into a C buffer.
/// Returns the number of bytes written (excluding the NUL terminator, which is also written).
///
/// # Safety
/// `dst` must have room for the string plus a NUL terminator.
unsafe fn write_char_field(dst: *mut u8, src: &[u8], count: usize) -> usize {
    if count == 0 {
        *dst = b'.';
        *dst.add(1) = 0;
        1
    } else {
        let str_len = char_field_len(src, count);
        std::ptr::copy_nonoverlapping(src.as_ptr(), dst, str_len);
        *dst.add(str_len) = 0;
        str_len
    }
}

/// Phase 1: Unpack ID + REF/ALT alleles from the shared block.
///
/// # Safety
/// `b` and `d` must be valid pointers to the record and its decoded-data struct.
unsafe fn unpack_str(shared: &[u8], b: *mut htslib::bcf1_t, d: *mut htslib::bcf_dec_t) -> i32 {
    let mut pos: usize = 0;

    // --- ID ---
    let id_start = pos;
    let (count, type_code, hdr_size) = match bcf_dec_size(shared.get(pos..).unwrap_or(&[])) {
        Some(v) => v,
        None => return -1,
    };
    pos += hdr_size;
    let id_total_bytes = bcf_type_bytes(count, type_code);

    let id_display_len = if count == 0 {
        1 // "."
    } else if type_code == htslib::BCF_BT_CHAR {
        char_field_len(&shared[pos..], count)
    } else {
        id_total_bytes
    };

    let needed = id_display_len + 1;
    if needed > (*d).m_id as usize {
        (*d).m_id = needed as i32;
        (*d).id = crealloc((*d).id, needed);
    }
    let src_end = (pos + id_total_bytes).min(shared.len());
    write_char_field((*d).id as *mut u8, &shared[pos..src_end], count);
    pos += id_total_bytes;
    (*b).unpack_size[0] = (pos - id_start) as i32;

    // --- REF + ALT alleles ---
    let allele_start = pos;
    let n_allele = (*b).n_allele() as usize;

    if n_allele > (*d).m_allele as usize {
        (*d).m_allele = n_allele as i32;
        (*d).allele = crealloc((*d).allele, n_allele);
    }

    // First pass: compute total als buffer size
    let mut scan_pos = pos;
    let mut total_als_len: usize = 0;
    for _ in 0..n_allele {
        let (cnt, tc, hs) = match bcf_dec_size(shared.get(scan_pos..).unwrap_or(&[])) {
            Some(v) => v,
            None => return -1,
        };
        scan_pos += hs;
        let data_bytes = bcf_type_bytes(cnt, tc);
        let str_len = if cnt == 0 {
            1
        } else if tc == htslib::BCF_BT_CHAR {
            char_field_len(&shared[scan_pos..], cnt)
        } else {
            data_bytes
        };
        total_als_len += str_len + 1;
        scan_pos += data_bytes;
    }

    if total_als_len > (*d).m_als as usize {
        (*d).m_als = total_als_len as i32;
        (*d).als = crealloc((*d).als, total_als_len);
    }

    // Second pass: copy allele strings
    let mut als_offset: usize = 0;
    for i in 0..n_allele {
        let (cnt, tc, hs) = match bcf_dec_size(shared.get(pos..).unwrap_or(&[])) {
            Some(v) => v,
            None => return -1,
        };
        pos += hs;
        let data_bytes = bcf_type_bytes(cnt, tc);

        *(*d).allele.add(i) = (*d).als.add(als_offset);
        let src_end = (pos + data_bytes).min(shared.len());
        let written = write_char_field(
            ((*d).als as *mut u8).add(als_offset),
            &shared[pos..src_end],
            cnt,
        );
        als_offset += written + 1; // written bytes + NUL
        pos += data_bytes;
    }
    (*b).unpack_size[1] = (pos - allele_start) as i32;
    (*b).unpacked |= htslib::BCF_UN_STR as i32;
    0
}

/// Phase 2: Unpack FILTER field.
///
/// # Safety
/// `b` and `d` must be valid pointers. `shared` must cover the full shared block.
unsafe fn unpack_flt(shared: &[u8], b: *mut htslib::bcf1_t, d: *mut htslib::bcf_dec_t) -> i32 {
    let flt_start = ((*b).unpack_size[0] + (*b).unpack_size[1]) as usize;
    let mut pos = flt_start;

    if pos >= shared.len() {
        return -1;
    }

    if shared[pos] >> 4 != 0 {
        let (count, type_code, hdr_size) = match bcf_dec_size(shared.get(pos..).unwrap_or(&[])) {
            Some(v) => v,
            None => return -1,
        };
        pos += hdr_size;
        (*d).n_flt = count as i32;

        if count > (*d).m_flt as usize {
            (*d).m_flt = count as i32;
            (*d).flt = crealloc((*d).flt, count);
        }

        for i in 0..count {
            let (val, consumed) = match bcf_dec_int1(shared.get(pos..).unwrap_or(&[]), type_code) {
                Some(v) => v,
                None => return -1,
            };
            *(*d).flt.add(i) = val as i32;
            pos += consumed;
        }
    } else {
        pos += 1;
        (*d).n_flt = 0;
    }
    (*b).unpack_size[2] = (pos - flt_start) as i32;
    (*b).unpacked |= htslib::BCF_UN_FLT as i32;
    0
}

/// Phase 3: Unpack INFO fields.
///
/// # Safety
/// `b` and `d` must be valid pointers. `shared` must cover the full shared block.
unsafe fn unpack_info(shared: &[u8], b: *mut htslib::bcf1_t, d: *mut htslib::bcf_dec_t) -> i32 {
    let info_start = ((*b).unpack_size[0] + (*b).unpack_size[1] + (*b).unpack_size[2]) as usize;
    let mut pos = info_start;
    let n_info = (*b).n_info() as usize;

    if n_info > (*d).m_info as usize {
        (*d).m_info = n_info as i32;
        (*d).info = crealloc((*d).info, n_info);
    }
    for i in 0..(*d).m_info as usize {
        (*(*d).info.add(i)).set_vptr_free(0);
    }

    for i in 0..n_info {
        let info = &mut *(*d).info.add(i);
        let ptr_start = pos;

        let (key, key_consumed) = match bcf_dec_typed_int1(shared.get(pos..).unwrap_or(&[])) {
            Some(v) => v,
            None => return -1,
        };
        info.key = key as i32;
        pos += key_consumed;

        let (len, type_code, hdr_size) = match bcf_dec_size(shared.get(pos..).unwrap_or(&[])) {
            Some(v) => v,
            None => return -1,
        };
        info.len = len as i32;
        info.type_ = type_code as i32;
        pos += hdr_size;

        info.vptr = shared.as_ptr().add(pos) as *mut u8;
        info.set_vptr_off((pos - ptr_start) as u32);
        info.set_vptr_free(0);
        info.v1.i = 0;

        let data_len;
        if len == 1 {
            // Pre-extract scalar value into v1 union
            match type_code {
                htslib::BCF_BT_INT8 | htslib::BCF_BT_CHAR => {
                    if pos < shared.len() {
                        info.v1.i = shared[pos] as i8 as i64;
                    }
                    data_len = len;
                }
                htslib::BCF_BT_FLOAT => {
                    if pos + 4 <= shared.len() {
                        info.v1.f = f32::from_le_bytes([
                            shared[pos],
                            shared[pos + 1],
                            shared[pos + 2],
                            shared[pos + 3],
                        ]);
                    }
                    data_len = bcf_type_bytes(len, type_code);
                }
                _ => {
                    // INT16, INT32, INT64
                    if let Some((v, _)) = bcf_dec_int1(shared.get(pos..).unwrap_or(&[]), type_code)
                    {
                        info.v1.i = v;
                    }
                    data_len = bcf_type_bytes(len, type_code);
                }
            }
        } else {
            data_len = bcf_type_bytes(len, type_code);
        }
        pos += data_len;
        info.vptr_len = (pos - ptr_start - info.vptr_off() as usize) as u32;
    }
    (*b).unpacked |= htslib::BCF_UN_INFO as i32;
    0
}

/// Phase 4: Unpack FORMAT + per-sample data from the indiv block.
///
/// # Safety
/// `b` and `d` must be valid pointers.
unsafe fn unpack_fmt(b: *mut htslib::bcf1_t, d: *mut htslib::bcf_dec_t) -> i32 {
    let n_fmt = (*b).n_fmt() as usize;
    let n_sample = (*b).n_sample() as usize;

    if n_fmt > (*d).m_fmt as usize {
        (*d).m_fmt = n_fmt as i32;
        (*d).fmt = crealloc((*d).fmt, n_fmt);
    }
    for i in 0..(*d).m_fmt as usize {
        (*(*d).fmt.add(i)).set_p_free(0);
    }

    let indiv = slice::from_raw_parts((*b).indiv.s as *const u8, (*b).indiv.l);
    let mut pos: usize = 0;

    for i in 0..n_fmt {
        let fmt = &mut *(*d).fmt.add(i);
        let ptr_start = pos;

        let (id, consumed) = match bcf_dec_typed_int1(indiv.get(pos..).unwrap_or(&[])) {
            Some(v) => v,
            None => return -1,
        };
        fmt.id = id as i32;
        pos += consumed;

        let (n, type_code, hdr_size) = match bcf_dec_size(indiv.get(pos..).unwrap_or(&[])) {
            Some(v) => v,
            None => return -1,
        };
        fmt.n = n as i32;
        fmt.type_ = type_code as i32;
        fmt.size = bcf_type_bytes(n, type_code) as i32;
        pos += hdr_size;

        fmt.p = indiv.as_ptr().add(pos) as *mut u8;
        fmt.set_p_off((pos - ptr_start) as u32);
        fmt.set_p_free(0);

        let data_len = n_sample * fmt.size as usize;
        pos += data_len;
        fmt.p_len = data_len as u32;
    }
    (*b).unpacked |= htslib::BCF_UN_FMT as i32;
    0
}

/// Common methods for numeric INFO and FORMAT entries
pub trait Numeric {
    /// Return true if entry is a missing value
    fn is_missing(&self) -> bool;

    /// Return missing value for storage in BCF record.
    fn missing() -> Self;
}

impl Numeric for f32 {
    fn is_missing(&self) -> bool {
        self.bits() == MISSING_FLOAT.bits()
    }

    fn missing() -> f32 {
        *MISSING_FLOAT
    }
}

impl Numeric for i32 {
    fn is_missing(&self) -> bool {
        *self == MISSING_INTEGER
    }

    fn missing() -> i32 {
        MISSING_INTEGER
    }
}

trait NumericUtils {
    /// Return true if entry marks the end of the record.
    fn is_vector_end(&self) -> bool;
}

impl NumericUtils for f32 {
    fn is_vector_end(&self) -> bool {
        self.bits() == VECTOR_END_FLOAT.bits()
    }
}

impl NumericUtils for i32 {
    fn is_vector_end(&self) -> bool {
        *self == VECTOR_END_INTEGER
    }
}

/// A trait to allow for seamless use of bytes or integer identifiers for filters
pub trait FilterId {
    fn id_from_header(&self, header: &HeaderView) -> Result<Id>;
    fn is_pass(&self) -> bool;
}

impl FilterId for [u8] {
    fn id_from_header(&self, header: &HeaderView) -> Result<Id> {
        let id = CString8::new(std::str::from_utf8(self).map_err(|_| Error::InvalidRecord)?)
            .map_err(|_| Error::InvalidRecord)?;
        header.name_to_id(&id)
    }
    fn is_pass(&self) -> bool {
        matches!(self, b"PASS" | b".")
    }
}

impl FilterId for &CStr8 {
    fn id_from_header(&self, header: &HeaderView) -> Result<Id> {
        header.name_to_id(self)
    }

    fn is_pass(&self) -> bool {
        matches!(self.as_bytes(), b"PASS" | b".")
    }
}

impl FilterId for Id {
    fn id_from_header(&self, _header: &HeaderView) -> Result<Id> {
        Ok(*self)
    }
    fn is_pass(&self) -> bool {
        *self == Id(0)
    }
}

/// A buffer for info or format data.
#[derive(Debug)]
pub struct Buffer {
    inner: *mut ::std::os::raw::c_void,
    len: i32,
}

impl Buffer {
    pub fn new() -> Self {
        Buffer {
            inner: ptr::null_mut(),
            len: 0,
        }
    }
}

impl Default for Buffer {
    fn default() -> Self {
        Self::new()
    }
}

impl Drop for Buffer {
    fn drop(&mut self) {
        // SAFETY: self.inner was allocated by htslib via bcf_get_info_values/bcf_get_format_values; free is symmetric.
        unsafe {
            ::libc::free(self.inner);
        }
    }
}

#[derive(new, Debug)]
pub struct BufferBacked<'a, T: 'a + fmt::Debug, B: Borrow<Buffer> + 'a> {
    value: T,
    _buffer: B,
    #[new(default)]
    phantom: PhantomData<&'a B>,
}

impl<'a, T: 'a + fmt::Debug, B: Borrow<Buffer> + 'a> Deref for BufferBacked<'a, T, B> {
    type Target = T;

    fn deref(&self) -> &T {
        &self.value
    }
}

impl<'a, T: 'a + fmt::Debug + fmt::Display, B: Borrow<Buffer> + 'a> fmt::Display
    for BufferBacked<'a, T, B>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.value, f)
    }
}

/// A VCF/BCF record.
/// New records can be created by the `empty_record` methods of [`bcf::Reader`](crate::bcf::Reader)
/// and [`bcf::Writer`](crate::bcf::Writer).
/// # Example
/// ```rust
/// use rust_htslib::bcf::{Format, Writer};
/// use rust_htslib::bcf::header::Header;
///
/// // Create minimal VCF header with a single sample
/// let mut header = Header::new();
/// header.push_sample("sample".as_bytes());
///
/// // Write uncompressed VCF to stdout with above header and get an empty record
/// let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
/// let mut record = vcf.empty_record();
/// ```
#[derive(Debug)]
pub struct Record {
    pub inner: *mut htslib::bcf1_t,
    header: Arc<HeaderView>,
}

impl Record {
    /// Construct record with reference to header `HeaderView`, for create-internal use.
    pub fn new(header: Arc<HeaderView>) -> Self {
        // SAFETY: bcf_init allocates a new record; bcf_unpack initializes it.
        let inner = unsafe {
            let inner = htslib::bcf_init();
            // Always unpack record.
            bcf_unpack_rs(inner, htslib::BCF_UN_ALL as i32);
            inner
        };
        Record { inner, header }
    }

    /// Force unpacking of internal record values.
    pub fn unpack(&mut self) {
        // SAFETY: self.inner is non-null (from constructor).
        unsafe { bcf_unpack_rs(self.inner, htslib::BCF_UN_ALL as i32) };
    }

    /// Return associated header.
    pub fn header(&self) -> &HeaderView {
        self.header.as_ref()
    }

    /// Translate the record to the given header.
    pub fn translate(&mut self, dst_header: &mut Arc<HeaderView>) -> Result<()> {
        // SAFETY: dst_header.inner, self.header().inner, and self.inner are non-null (from constructors).
        if unsafe { htslib::bcf_translate(dst_header.inner, self.header().inner, self.inner) } == 0
        {
            self.set_header(Arc::clone(dst_header));
            Ok(())
        } else {
            Err(Error::Translate)
        }
    }

    /// Set the record header.
    pub(crate) fn set_header(&mut self, header: Arc<HeaderView>) {
        self.header = header;
    }

    /// Return reference to the inner C struct.
    ///
    /// # Remarks
    ///
    /// Note that this function is only required as long as Rust-Htslib does not provide full
    /// access to all aspects of Htslib.
    pub fn inner(&self) -> &htslib::bcf1_t {
        // SAFETY: self.inner is non-null (from constructor or bcf_dup).
        unsafe { &*self.inner }
    }

    /// Return mutable reference to inner C struct.
    ///
    /// # Remarks
    ///
    /// Note that this function is only required as long as Rust-Htslib does not provide full
    /// access to all aspects of Htslib.
    pub fn inner_mut(&mut self) -> &mut htslib::bcf1_t {
        // SAFETY: self.inner is non-null (from constructor or bcf_dup); we have &mut self.
        unsafe { &mut *self.inner }
    }

    /// Get the reference id of the record.
    ///
    /// To look up the contig name,
    /// use [`HeaderView::rid2name`](../header/struct.HeaderView.html#method.rid2name).
    ///
    /// # Returns
    ///
    /// - `Some(rid)` if the internal `rid` is set to a value that is not `-1`
    /// - `None` if the internal `rid` is set to `-1`
    pub fn rid(&self) -> Option<u32> {
        match self.inner().rid {
            -1 => None,
            rid => Some(rid as u32),
        }
    }

    /// Update the reference id of the record.
    ///
    /// To look up reference id for a contig name,
    /// use [`HeaderView::name2rid`](../header/struct.HeaderView.html#method.name2rid).
    ///
    /// # Example
    ///
    /// Example assumes we have a Record `record` from a VCF with a header containing region
    /// named `1`. See [module documentation](../index.html#example-writing) for how to set
    /// up VCF, header, and record.
    ///
    /// ```
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// # let mut header = Header::new();
    /// # let header_contig_line = r#"##contig=<ID=1,length=10>"#;
    /// # header.push_record(header_contig_line.as_bytes());
    /// # header.push_sample("test_sample".as_bytes());
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let rid = record.header().name2rid(b"1").ok();
    /// record.set_rid(rid);
    /// assert_eq!(record.rid(), rid);
    /// let name = record.header().rid2name(record.rid().unwrap()).ok();
    /// assert_eq!(Some("1".as_bytes()), name);
    /// ```
    pub fn set_rid(&mut self, rid: Option<u32>) {
        match rid {
            Some(rid) => self.inner_mut().rid = rid as i32,
            None => self.inner_mut().rid = -1,
        }
    }

    /// Return **0-based** position
    pub fn pos(&self) -> i64 {
        self.inner().pos
    }

    /// Set **0-based** position
    pub fn set_pos(&mut self, pos: i64) {
        self.inner_mut().pos = pos;
    }

    /// Return the **0-based, exclusive** end position
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Header, Writer};
    /// # use tempfile::NamedTempFile;
    /// # let tmp = NamedTempFile::new().unwrap();
    /// # let path = tmp.path();
    /// # let header = Header::new();
    /// # let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let alleles: &[&[u8]] = &[b"AGG", b"TG"];
    /// record.set_alleles(alleles).expect("Failed to set alleles");
    /// record.set_pos(5);
    ///
    /// assert_eq!(record.end(), 8)
    /// ```
    pub fn end(&self) -> i64 {
        self.pos() + self.rlen()
    }

    /// Return the value of the ID column.
    ///
    /// When empty, returns `b".".to_vec()`.
    pub fn id(&self) -> Vec<u8> {
        if self.inner().d.id.is_null() {
            b".".to_vec()
        } else {
            // SAFETY: d.id is checked for null above; it is a valid NUL-terminated C string from htslib.
            let id = unsafe { ffi::CStr::from_ptr(self.inner().d.id) };
            id.to_bytes().to_vec()
        }
    }

    /// Update the ID string to the given value.
    pub fn set_id(&mut self, id: &[u8]) -> Result<()> {
        let c_str = ffi::CString::new(id).unwrap();
        // SAFETY: self.header().inner and self.inner are non-null (from constructor); c_str is valid.
        if unsafe {
            htslib::bcf_update_id(
                self.header().inner,
                self.inner,
                c_str.as_ptr() as *mut c_char,
            )
        } == 0
        {
            Ok(())
        } else {
            Err(Error::SetValues)
        }
    }

    /// Clear the ID column (set it to `"."`).
    pub fn clear_id(&mut self) -> Result<()> {
        let c_str = ffi::CString::new(&b"."[..]).unwrap();
        // SAFETY: self.header().inner and self.inner are non-null (from constructor); c_str is valid.
        if unsafe {
            htslib::bcf_update_id(
                self.header().inner,
                self.inner,
                c_str.as_ptr() as *mut c_char,
            )
        } == 0
        {
            Ok(())
        } else {
            Err(Error::SetValues)
        }
    }

    /// Add the ID string (the ID field is semicolon-separated), checking for duplicates.
    pub fn push_id(&mut self, id: &[u8]) -> Result<()> {
        let c_str = ffi::CString::new(id).unwrap();
        // SAFETY: self.header().inner and self.inner are non-null (from constructor); c_str is valid.
        if unsafe {
            htslib::bcf_add_id(
                self.header().inner,
                self.inner,
                c_str.as_ptr() as *mut c_char,
            )
        } == 0
        {
            Ok(())
        } else {
            Err(Error::SetValues)
        }
    }

    /// Return `Filters` iterator for enumerating all filters that have been set.
    ///
    /// A record having the `PASS` filter will return an empty `Filter` here.
    pub fn filters(&self) -> Filters<'_> {
        Filters::new(self)
    }

    /// Query whether the filter with the given ID has been set.
    ///
    /// This method can be used to check if a record passes filtering by using either `Id(0)`,
    /// `PASS` or `.`
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Header, Writer};
    /// # use rust_htslib::bcf::header::Id;
    /// # use tempfile::NamedTempFile;
    /// # let tmp = tempfile::NamedTempFile::new().unwrap();
    /// # let path = tmp.path();
    /// let mut header = Header::new();
    /// header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
    /// # let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// assert!(record.has_filter("PASS".as_bytes()));
    /// assert!(record.has_filter(".".as_bytes()));
    /// assert!(record.has_filter(&Id(0)));
    ///
    /// record.push_filter("foo".as_bytes()).unwrap();
    /// assert!(record.has_filter("foo".as_bytes()));
    /// assert!(!record.has_filter("PASS".as_bytes()))
    /// ```
    pub fn has_filter<T: FilterId + ?Sized>(&self, flt_id: &T) -> bool {
        if flt_id.is_pass() && self.inner().d.n_flt == 0 {
            return true;
        }
        let id = match flt_id.id_from_header(self.header()) {
            Ok(i) => *i,
            Err(_) => return false,
        };
        for i in 0..(self.inner().d.n_flt as isize) {
            // SAFETY: i is within [0, n_flt); d.flt is a valid pointer with n_flt elements.
            if unsafe { *self.inner().d.flt.offset(i) } == id as i32 {
                return true;
            }
        }
        false
    }

    /// Set the given filter IDs to the FILTER column.
    ///
    /// Setting an empty slice removes all filters and sets `PASS`.
    ///
    /// # Example
    /// ```rust
    /// # use cstr8::cstr8;
    /// # use rust_htslib::bcf::{Format, Header, Writer};
    /// # use rust_htslib::bcf::header::Id;
    /// # use tempfile::NamedTempFile;
    /// # let tmp = tempfile::NamedTempFile::new().unwrap();
    /// # let path = tmp.path();
    /// let mut header = Header::new();
    /// header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
    /// header.push_record(br#"##FILTER=<ID=bar,Description="a horse walks into...">"#);
    /// # let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let foo = record.header().name_to_id(cstr8!("foo")).unwrap();
    /// let bar = record.header().name_to_id(cstr8!("bar")).unwrap();
    /// assert!(record.has_filter("PASS".as_bytes()));
    /// let mut filters = vec![&foo, &bar];
    /// record.set_filters(&filters).unwrap();
    /// assert!(record.has_filter(&foo));
    /// assert!(record.has_filter(&bar));
    /// assert!(!record.has_filter("PASS".as_bytes()));
    /// filters.clear();
    /// record.set_filters(&filters).unwrap();
    /// assert!(record.has_filter("PASS".as_bytes()));
    /// assert!(!record.has_filter("foo".as_bytes()));
    /// // 'baz' isn't in the header
    /// assert!(record.set_filters(&["baz".as_bytes()]).is_err())
    /// ```
    ///
    /// # Errors
    /// If any of the filter IDs do not exist in the header, an [`BcfError::UnknownID`] is returned.
    ///
    pub fn set_filters<T: FilterId + ?Sized>(&mut self, flt_ids: &[&T]) -> Result<()> {
        let mut ids: Vec<i32> = flt_ids
            .iter()
            .map(|id| id.id_from_header(self.header()).map(|id| *id as i32))
            .collect::<Result<Vec<i32>>>()?;
        // SAFETY: self.header().inner and self.inner are non-null (from constructor); ids is a valid slice.
        unsafe {
            htslib::bcf_update_filter(
                self.header().inner,
                self.inner,
                ids.as_mut_ptr(),
                ids.len() as i32,
            );
        };
        Ok(())
    }

    /// Add the given filter to the FILTER column.
    ///
    /// If `flt_id` is `PASS` or `.` then all existing filters are removed first. Otherwise,
    /// any existing `PASS` filter is removed.
    ///
    /// # Example
    /// ```rust
    /// # use cstr8::cstr8;
    /// # use rust_htslib::bcf::{Format, Header, Writer};
    /// # use tempfile::NamedTempFile;
    /// # let tmp = tempfile::NamedTempFile::new().unwrap();
    /// # let path = tmp.path();
    /// let mut header = Header::new();
    /// header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
    /// header.push_record(br#"##FILTER=<ID=bar,Description="dranks">"#);
    /// # let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let foo = "foo".as_bytes();
    /// let bar = record.header().name_to_id(cstr8!("bar")).unwrap();
    /// assert!(record.has_filter("PASS".as_bytes()));
    ///
    /// record.push_filter(foo).unwrap();
    /// record.push_filter(&bar).unwrap();
    /// assert!(record.has_filter(foo));
    /// assert!(record.has_filter(&bar));
    /// // filter must exist in the header
    /// assert!(record.push_filter("baz".as_bytes()).is_err())
    /// ```
    ///
    /// # Errors
    /// If the `flt_id` does not exist in the header, an [`BcfError::UnknownID`] is returned.
    ///
    pub fn push_filter<T: FilterId + ?Sized>(&mut self, flt_id: &T) -> Result<()> {
        let id = flt_id.id_from_header(self.header())?;
        // SAFETY: self.header().inner and self.inner are non-null (from constructor); id is valid.
        unsafe {
            htslib::bcf_add_filter(self.header().inner, self.inner, *id as i32);
        };
        Ok(())
    }

    /// Remove the given filter from the FILTER column.
    ///
    /// # Arguments
    ///
    /// - `flt_id` - The corresponding filter ID to remove.
    /// - `pass_on_empty` - Set to `PASS` when removing the last filter.
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Header, Writer};
    /// # use tempfile::NamedTempFile;
    /// # let tmp = tempfile::NamedTempFile::new().unwrap();
    /// # let path = tmp.path();
    /// let mut header = Header::new();
    /// header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
    /// header.push_record(br#"##FILTER=<ID=bar,Description="a horse walks into...">"#);
    /// # let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let foo = "foo".as_bytes();
    /// let bar = "bar".as_bytes();
    /// record.set_filters(&[foo, bar]).unwrap();
    /// assert!(record.has_filter(foo));
    /// assert!(record.has_filter(bar));
    ///
    /// record.remove_filter(foo, true).unwrap();
    /// assert!(!record.has_filter(foo));
    /// assert!(record.has_filter(bar));
    /// // 'baz' is not in the header
    /// assert!(record.remove_filter("baz".as_bytes(), true).is_err());
    ///
    /// record.remove_filter(bar, true).unwrap();
    /// assert!(!record.has_filter(bar));
    /// assert!(record.has_filter("PASS".as_bytes()));
    /// ```
    ///
    /// # Errors
    /// If the `flt_id` does not exist in the header, an [`BcfError::UnknownID`] is returned.
    ///
    pub fn remove_filter<T: FilterId + ?Sized>(
        &mut self,
        flt_id: &T,
        pass_on_empty: bool,
    ) -> Result<()> {
        let id = flt_id.id_from_header(self.header())?;
        // SAFETY: self.header().inner and self.inner are non-null (from constructor); id is valid.
        unsafe {
            htslib::bcf_remove_filter(
                self.header().inner,
                self.inner,
                *id as i32,
                pass_on_empty as i32,
            )
        };
        Ok(())
    }

    /// Get alleles strings.
    ///
    /// The first allele is the reference allele.
    pub fn alleles(&self) -> Vec<&[u8]> {
        // SAFETY: self.inner is non-null (from constructor).
        unsafe { bcf_unpack_rs(self.inner, htslib::BCF_UN_ALL as i32) };
        let n = self.inner().n_allele() as usize;
        let dec = self.inner().d;
        // SAFETY: dec.allele points to n valid C string pointers after bcf_unpack.
        let alleles = unsafe { slice::from_raw_parts(dec.allele, n) };
        (0..n)
            // SAFETY: each alleles[i] is a valid NUL-terminated C string from htslib.
            .map(|i| unsafe { ffi::CStr::from_ptr(alleles[i]).to_bytes() })
            .collect()
    }

    /// Set alleles. The first allele is the reference allele.
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// #
    /// # // Create minimal VCF header with a single sample
    /// # let mut header = Header::new();
    /// # header.push_sample("sample".as_bytes());
    /// #
    /// # // Write uncompressed VCF to stdout with above header and get an empty record
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// assert_eq!(record.allele_count(), 0);
    ///
    /// let alleles: &[&[u8]] = &[b"A", b"TG"];
    /// record.set_alleles(alleles).expect("Failed to set alleles");
    /// assert_eq!(record.allele_count(), 2)
    /// ```
    pub fn set_alleles(&mut self, alleles: &[&[u8]]) -> Result<()> {
        let cstrings: Vec<ffi::CString> = alleles
            .iter()
            .map(|vec| ffi::CString::new(*vec).unwrap())
            .collect();
        let mut ptrs: Vec<*const c_char> = cstrings
            .iter()
            .map(|cstr| cstr.as_ptr() as *const c_char)
            .collect();
        // SAFETY: self.header().inner and self.inner are non-null; ptrs contains valid CString pointers.
        if unsafe {
            htslib::bcf_update_alleles(
                self.header().inner,
                self.inner,
                ptrs.as_mut_ptr(),
                alleles.len() as i32,
            )
        } == 0
        {
            Ok(())
        } else {
            Err(Error::SetValues)
        }
    }

    /// Get variant quality.
    pub fn qual(&self) -> f32 {
        self.inner().qual
    }

    /// Set variant quality.
    pub fn set_qual(&mut self, qual: f32) {
        self.inner_mut().qual = qual;
    }

    pub fn info<'a>(&'a self, tag: &'a [u8]) -> Info<'a, Buffer> {
        self.info_shared_buffer(tag, Buffer::new())
    }

    /// Get the value of the given info tag.
    pub fn info_shared_buffer<'a, 'b, B: BorrowMut<Buffer> + Borrow<Buffer> + 'b>(
        &'a self,
        tag: &'a [u8],
        buffer: B,
    ) -> Info<'a, B> {
        Info {
            record: self,
            tag,
            buffer,
        }
    }

    /// Get the number of samples in the record.
    pub fn sample_count(&self) -> u32 {
        self.inner().n_sample()
    }

    /// Get the number of alleles, including reference allele.
    pub fn allele_count(&self) -> u32 {
        self.inner().n_allele()
    }

    /// Add/replace genotypes in FORMAT GT tag.
    ///
    /// # Arguments
    ///
    /// - `genotypes` - a flattened, two-dimensional array of GenotypeAllele,
    ///   the first dimension contains one array for each sample.
    ///
    /// # Errors
    ///
    /// Returns error if GT tag is not present in header.
    ///
    /// # Example
    ///
    /// Example assumes we have a Record `record` from a VCF with a `GT` `FORMAT` tag.
    /// See [module documentation](../index.html#example-writing) for how to set up
    /// VCF, header, and record.
    ///
    /// ```
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// # use rust_htslib::bcf::record::GenotypeAllele;
    /// # let mut header = Header::new();
    /// # let header_contig_line = r#"##contig=<ID=1,length=10>"#;
    /// # header.push_record(header_contig_line.as_bytes());
    /// # let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
    /// # header.push_record(header_gt_line.as_bytes());
    /// # header.push_sample("test_sample".as_bytes());
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let alleles = &[GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)];
    /// record.push_genotypes(alleles);
    /// assert_eq!("1/1", &format!("{}", record.genotypes().unwrap().get(0)));
    /// ```
    pub fn push_genotypes(&mut self, genotypes: &[GenotypeAllele]) -> Result<()> {
        let encoded: Vec<i32> = genotypes.iter().map(|gt| i32::from(*gt)).collect();
        self.push_format_integer(cstr8!("GT"), &encoded)
    }

    /// Add/replace genotypes in FORMAT GT tag by providing a list of genotypes.
    ///
    /// # Arguments
    ///
    /// - `genotypes` - a two-dimensional array of GenotypeAllele
    /// - `max_ploidy` - the maximum number of alleles allowed for any genotype on any sample.
    ///
    /// # Errors
    ///
    /// Returns an error if any genotype has more allelles than `max_ploidy` or if the GT tag is not present in the header.
    ///
    /// # Example
    ///
    /// Example assumes we have a Record `record` from a VCF with a `GT` `FORMAT` tag and three samples.
    /// See [module documentation](../index.html#example-writing) for how to set up
    /// VCF, header, and record.
    ///
    /// ```
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// # use rust_htslib::bcf::record::GenotypeAllele;
    /// # use std::iter;
    /// # let mut header = Header::new();
    /// # let header_contig_line = r#"##contig=<ID=1,length=10>"#;
    /// # header.push_record(header_contig_line.as_bytes());
    /// # let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
    /// # header.push_record(header_gt_line.as_bytes());
    /// # header.push_sample("first_sample".as_bytes());
    /// # header.push_sample("second_sample".as_bytes());
    /// # header.push_sample("third_sample".as_bytes());
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf)?;
    /// # let mut record = vcf.empty_record();
    /// let alleles = vec![
    ///     vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)],
    ///     vec![GenotypeAllele::Unphased(0), GenotypeAllele::Phased(1)],
    ///     vec![GenotypeAllele::Unphased(0)],
    /// ];
    /// record.push_genotype_structured(&alleles, 2);
    /// let gts = record.genotypes()?;
    /// assert_eq!("1/1", &format!("{}", gts.get(0)));
    /// assert_eq!("0|1", &format!("{}", gts.get(1)));
    /// assert_eq!("0", &format!("{}", gts.get(2)));
    /// # Ok::<(), rust_htslib::bcf::BcfError>(())
    /// ```
    pub fn push_genotype_structured<GT>(
        &mut self,
        genotypes: &[GT],
        max_ploidy: usize,
    ) -> Result<()>
    where
        GT: AsRef<[GenotypeAllele]>,
    {
        let mut data = Vec::with_capacity(max_ploidy * genotypes.len());
        for gt in genotypes {
            if gt.as_ref().len() > max_ploidy {
                return Err(Error::SetValues);
            }
            data.extend(
                gt.as_ref()
                    .iter()
                    .map(|gta| i32::from(*gta))
                    .chain(iter::repeat_n(
                        VECTOR_END_INTEGER,
                        max_ploidy - gt.as_ref().len(),
                    )),
            );
        }
        self.push_format_integer(cstr8!("GT"), &data)
    }

    /// Get genotypes as vector of one `Genotype` per sample.
    ///
    /// # Example
    /// Parsing genotype field (`GT` tag) from a VCF record:
    /// ```
    /// use crate::rust_htslib::bcf::{Reader, Read};
    /// let mut vcf = Reader::from_path(&"test/test_string.vcf").expect("Error opening file.");
    /// let expected = ["./1", "1|1", "0/1", "0|1", "1|.", "1/1"];
    /// for (rec, exp_gt) in vcf.records().zip(expected.iter()) {
    ///     let mut rec = rec.expect("Error reading record.");
    ///     let genotypes = rec.genotypes().expect("Error reading genotypes");
    ///     assert_eq!(&format!("{}", genotypes.get(0)), exp_gt);
    /// }
    /// ```
    pub fn genotypes(&self) -> Result<Genotypes<'_, Buffer>> {
        self.genotypes_shared_buffer(Buffer::new())
    }

    /// Get genotypes as vector of one `Genotype` per sample, using a given shared buffer
    /// to avoid unnecessary allocations.
    pub fn genotypes_shared_buffer<'a, B>(&self, buffer: B) -> Result<Genotypes<'a, B>>
    where
        B: BorrowMut<Buffer> + Borrow<Buffer> + 'a,
    {
        Ok(Genotypes {
            encoded: self.format_shared_buffer(b"GT", buffer).integer()?,
        })
    }

    /// Retrieve data for a `FORMAT` field
    ///
    /// # Example
    /// *Note: some boilerplate for the example is hidden for clarity. See [module documentation](../index.html#example-writing)
    /// for an example of the setup used here.*
    ///
    /// ```rust
    /// # use cstr8::cstr8;
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// #
    /// # // Create minimal VCF header with a single sample
    /// # let mut header = Header::new();
    /// header.push_sample(b"sample1").push_sample(b"sample2").push_record(br#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">"#);
    /// #
    /// # // Write uncompressed VCF to stdout with above header and get an empty record
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// record.push_format_integer(cstr8!("DP"), &[20, 12]).expect("Failed to set DP format field");
    ///
    /// let read_depths = record.format(b"DP").integer().expect("Couldn't retrieve DP field");
    /// let sample1_depth = read_depths[0];
    /// assert_eq!(sample1_depth, &[20]);
    /// let sample2_depth = read_depths[1];
    /// assert_eq!(sample2_depth, &[12])
    /// ```
    ///
    /// # Errors
    /// **Attention:** the returned [`BufferBacked`] from [`integer()`](Format::integer)
    /// (`read_depths`), which holds the data, has to be kept in scope as long as the data is
    /// accessed. If parts of the data are accessed after the `BufferBacked` object is been
    /// dropped, you will access unallocated memory.
    pub fn format<'a>(&'a self, tag: &'a [u8]) -> Format<'a, Buffer> {
        self.format_shared_buffer(tag, Buffer::new())
    }

    /// Get the value of the given format tag for each sample.
    pub fn format_shared_buffer<'a, 'b, B: BorrowMut<Buffer> + Borrow<Buffer> + 'b>(
        &'a self,
        tag: &'a [u8],
        buffer: B,
    ) -> Format<'a, B> {
        Format::new(self, tag, buffer)
    }

    /// Add/replace an integer-typed FORMAT tag.
    ///
    /// # Arguments
    ///
    /// - `tag` - The tag's string.
    /// - `data` - a flattened, two-dimensional array, the first dimension contains one array
    ///   for each sample.
    ///
    /// # Errors
    ///
    /// Returns error if tag is not present in header.
    pub fn push_format_integer(&mut self, tag: &CStr8, data: &[i32]) -> Result<()> {
        self.push_format(tag, data, htslib::BCF_HT_INT)
    }

    /// Add/replace a float-typed FORMAT tag.
    ///
    /// # Arguments
    ///
    /// - `tag` - The tag's string.
    /// - `data` - a flattened, two-dimensional array, the first dimension contains one array
    ///   for each sample.
    ///
    /// # Errors
    ///
    /// Returns error if tag is not present in header.
    ///
    /// # Example
    ///
    /// Example assumes we have a Record `record` from a VCF with an `AF` `FORMAT` tag.
    /// See [module documentation](../index.html#example-writing) for how to set up
    /// VCF, header, and record.
    ///
    /// ```
    /// # use cstr8::cstr8;
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// # use rust_htslib::bcf::record::GenotypeAllele;
    /// # let mut header = Header::new();
    /// # let header_contig_line = r#"##contig=<ID=1,length=10>"#;
    /// # header.push_record(header_contig_line.as_bytes());
    /// # let header_af_line = r#"##FORMAT=<ID=AF,Number=1,Type=Float,Description="Frequency">"#;
    /// # header.push_record(header_af_line.as_bytes());
    /// # header.push_sample("test_sample".as_bytes());
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// record.push_format_float(cstr8!("AF"), &[0.5]);
    /// assert_eq!(0.5, record.format(b"AF").float().unwrap()[0][0]);
    /// ```
    pub fn push_format_float(&mut self, tag: &CStr8, data: &[f32]) -> Result<()> {
        self.push_format(tag, data, htslib::BCF_HT_REAL)
    }

    /// Add/replace a single-char-typed FORMAT tag.
    ///
    /// # Arguments
    ///
    /// - `tag` - The tag's string.
    /// - `data` - a flattened, two-dimensional array, the first dimension contains one array
    ///   for each sample.
    ///
    /// # Errors
    ///
    /// Returns error if tag is not present in header.
    pub fn push_format_char(&mut self, tag: &CStr8, data: &[u8]) -> Result<()> {
        self.push_format(tag, data, htslib::BCF_HT_STR)
    }

    /// Add a format tag. Data is a flattened two-dimensional array.
    /// The first dimension contains one array for each sample.
    fn push_format<T>(&mut self, tag: &CStr8, data: &[T], ht: u32) -> Result<()> {
        // SAFETY: self.header().inner and self.inner are non-null; tag is a valid CStr8; data is a valid slice.
        unsafe {
            if htslib::bcf_update_format(
                self.header().inner,
                self.inner,
                tag.as_ptr() as *mut c_char,
                data.as_ptr() as *const ::std::os::raw::c_void,
                data.len() as i32,
                ht as i32,
            ) == 0
            {
                Ok(())
            } else {
                Err(Error::SetTag { tag: tag.into() })
            }
        }
    }

    // TODO: should we add convenience methods clear_format_*?

    /// Add a string-typed FORMAT tag. Note that genotypes are treated as a special case
    /// and cannot be added with this method. See instead [push_genotypes](#method.push_genotypes).
    ///
    /// # Arguments
    ///
    /// - `tag` - The tag's string.
    /// - `data` - a two-dimensional array, the first dimension contains one array
    ///   for each sample. Must be non-empty.
    ///
    /// # Errors
    ///
    /// Returns error if tag is not present in header.
    pub fn push_format_string<D: Borrow<[u8]>>(&mut self, tag: &CStr8, data: &[D]) -> Result<()> {
        assert!(
            !data.is_empty(),
            "given string data must have at least 1 element"
        );
        let c_data = data
            .iter()
            .map(|s| ffi::CString::new(s.borrow()).unwrap())
            .collect::<Vec<ffi::CString>>();
        let c_ptrs = c_data
            .iter()
            .map(|s| s.as_ptr() as *mut i8)
            .collect::<Vec<*mut i8>>();
        // SAFETY: self.header().inner and self.inner are non-null; tag is valid; c_ptrs contains valid CString pointers.
        unsafe {
            if htslib::bcf_update_format_string(
                self.header().inner,
                self.inner,
                tag.as_ptr() as *mut c_char,
                c_ptrs.as_slice().as_ptr() as *mut *const c_char,
                data.len() as i32,
            ) == 0
            {
                Ok(())
            } else {
                Err(Error::SetTag { tag: tag.into() })
            }
        }
    }

    /// Add/replace an integer-typed INFO entry.
    pub fn push_info_integer(&mut self, tag: &CStr8, data: &[i32]) -> Result<()> {
        self.push_info(tag, data, htslib::BCF_HT_INT)
    }

    /// Remove the integer-typed INFO entry.
    pub fn clear_info_integer(&mut self, tag: &CStr8) -> Result<()> {
        self.push_info::<i32>(tag, &[], htslib::BCF_HT_INT)
    }

    /// Add/replace a float-typed INFO entry.
    pub fn push_info_float(&mut self, tag: &CStr8, data: &[f32]) -> Result<()> {
        self.push_info(tag, data, htslib::BCF_HT_REAL)
    }

    /// Remove the float-typed INFO entry.
    pub fn clear_info_float(&mut self, tag: &CStr8) -> Result<()> {
        self.push_info::<u8>(tag, &[], htslib::BCF_HT_REAL)
    }

    /// Add/replace an INFO tag.
    ///
    /// # Arguments
    /// * `tag` - the tag to add/replace
    /// * `data` - the data to set
    /// * `ht` - the HTSLib type to use
    fn push_info<T>(&mut self, tag: &CStr8, data: &[T], ht: u32) -> Result<()> {
        // SAFETY: self.header().inner and self.inner are non-null; tag is a valid CStr8; data is a valid slice.
        unsafe {
            if htslib::bcf_update_info(
                self.header().inner,
                self.inner,
                tag.as_ptr() as *mut c_char,
                data.as_ptr() as *const ::std::os::raw::c_void,
                data.len() as i32,
                ht as i32,
            ) == 0
            {
                Ok(())
            } else {
                Err(Error::SetTag { tag: tag.into() })
            }
        }
    }

    /// Set flag into the INFO column.
    pub fn push_info_flag(&mut self, tag: &CStr8) -> Result<()> {
        self.push_info_string_impl(tag, &[b""], htslib::BCF_HT_FLAG)
    }

    /// Remove the flag from the INFO column.
    pub fn clear_info_flag(&mut self, tag: &CStr8) -> Result<()> {
        self.push_info_string_impl(tag, &[], htslib::BCF_HT_FLAG)
    }

    /// Add/replace a string-typed INFO entry.
    pub fn push_info_string(&mut self, tag: &CStr8, data: &[&[u8]]) -> Result<()> {
        self.push_info_string_impl(tag, data, htslib::BCF_HT_STR)
    }

    /// Remove the string field from the INFO column.
    pub fn clear_info_string(&mut self, tag: &CStr8) -> Result<()> {
        self.push_info_string_impl(tag, &[], htslib::BCF_HT_STR)
    }

    /// Call `bcf_update_info` with a pre-validated pointer and length.
    ///
    /// # Safety
    /// The caller must ensure that `ptr` points to a valid memory region of at least `len`
    /// elements of the type described by `ht`, and that `ht` is a valid BCF_HT_* constant.
    fn update_info_unchecked(
        &mut self,
        tag: &CStr8,
        ptr: *const ::std::os::raw::c_void,
        len: i32,
        ht: u32,
    ) -> Result<()> {
        // SAFETY: caller guarantees ptr and len are valid; self.header().inner and self.inner are non-null.
        unsafe {
            if htslib::bcf_update_info(
                self.header().inner,
                self.inner,
                tag.as_ptr() as *mut c_char,
                ptr,
                len,
                ht as i32,
            ) == 0
            {
                Ok(())
            } else {
                Err(Error::SetTag { tag: tag.into() })
            }
        }
    }

    /// Add an string-valued INFO tag.
    fn push_info_string_impl(&mut self, tag: &CStr8, data: &[&[u8]], ht: u32) -> Result<()> {
        if data.is_empty() {
            // Clear the tag
            // SAFETY: empty_str is a valid NUL-terminated string; len=0 signals deletion.
            let empty_str = unsafe { CStr8::from_utf8_with_nul_unchecked(b"\0") };
            return self.update_info_unchecked(
                tag,
                empty_str.as_ptr() as *const ::std::os::raw::c_void,
                0,
                ht,
            );
        }

        if data == [b""] {
            // This is a flag
            // SAFETY: empty_str is a valid NUL-terminated string; len=1 signals flag presence.
            let empty_str = unsafe { CStr8::from_utf8_with_nul_unchecked(b"\0") };
            return self.update_info_unchecked(
                tag,
                empty_str.as_ptr() as *const ::std::os::raw::c_void,
                1,
                ht,
            );
        }

        let data_bytes: usize = data.iter().map(|x| x.len()).sum::<usize>() + data.len();
        let mut buf: Vec<u8> = Vec::with_capacity(data_bytes);
        for (i, &s) in data.iter().enumerate() {
            if i > 0 {
                buf.extend(b",");
            }
            buf.extend(s);
        }
        let c_str = ffi::CString::new(buf).unwrap();
        let len = if ht == htslib::BCF_HT_FLAG {
            data.len()
        } else {
            c_str.to_bytes().len()
        };
        // SAFETY: self.header().inner and self.inner are non-null; c_str is a valid CString; len matches data.
        unsafe {
            if htslib::bcf_update_info(
                self.header().inner,
                self.inner,
                tag.as_ptr() as *mut c_char,
                c_str.as_ptr() as *const ::std::os::raw::c_void,
                len as i32,
                ht as i32,
            ) == 0
            {
                Ok(())
            } else {
                Err(Error::SetTag { tag: tag.into() })
            }
        }
    }

    /// Remove unused alleles.
    pub fn trim_alleles(&mut self) -> Result<()> {
        // SAFETY: self.header().inner and self.inner are non-null (from constructor).
        match unsafe { htslib::bcf_trim_alleles(self.header().inner, self.inner) } {
            -1 => Err(Error::RemoveAlleles),
            _ => Ok(()),
        }
    }

    pub fn remove_alleles(&mut self, remove: &[bool]) -> Result<()> {
        let rm_set = KBitSet::from_bools(remove);
        // SAFETY: self.header().inner and self.inner are non-null; rm_set is a valid C-compatible kbitset_t.
        let ret = unsafe {
            htslib::bcf_remove_allele_set(self.header().inner, self.inner, rm_set.as_ptr())
        };
        // rm_set is freed on drop
        match ret {
            -1 => Err(Error::RemoveAlleles),
            _ => Ok(()),
        }
    }

    /// Get the length of the reference allele. If the record has no reference allele, then the
    /// result will be `0`.
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// #
    /// # // Create minimal VCF header with a single sample
    /// # let mut header = Header::new();
    /// # header.push_sample("sample".as_bytes());
    /// #
    /// # // Write uncompressed VCF to stdout with above header and get an empty record
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// # assert_eq!(record.rlen(), 0);
    /// let alleles: &[&[u8]] = &[b"AGG", b"TG"];
    /// record.set_alleles(alleles).expect("Failed to set alleles");
    /// assert_eq!(record.rlen(), 3)
    /// ```
    pub fn rlen(&self) -> i64 {
        self.inner().rlen
    }

    /// Clear all parts of the record. Useful if you plan to reuse a record object multiple times.
    ///
    /// # Example
    /// ```rust
    /// # use rust_htslib::bcf::{Format, Writer};
    /// # use rust_htslib::bcf::header::Header;
    /// #
    /// # // Create minimal VCF header with a single sample
    /// # let mut header = Header::new();
    /// # header.push_sample("sample".as_bytes());
    /// #
    /// # // Write uncompressed VCF to stdout with above header and get an empty record
    /// # let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// # let mut record = vcf.empty_record();
    /// let alleles: &[&[u8]] = &[b"AGG", b"TG"];
    /// record.set_alleles(alleles).expect("Failed to set alleles");
    /// record.set_pos(6);
    /// record.clear();
    /// assert_eq!(record.rlen(), 0);
    /// assert_eq!(record.pos(), 0)
    /// ```
    pub fn clear(&self) {
        // SAFETY: self.inner is non-null (from constructor).
        unsafe { htslib::bcf_clear(self.inner) }
    }

    /// Provide short description of record for locating it in the BCF/VCF file.
    pub fn desc(&self) -> String {
        if let Some(rid) = self.rid() {
            if let Ok(contig) = self.header.rid2name(rid) {
                return format!("{}:{}", str::from_utf8(contig).unwrap(), self.pos());
            }
        }
        "".to_owned()
    }

    /// Convert to VCF String
    ///
    /// Intended for debug only. Use Writer for efficient VCF output.
    ///
    pub fn to_vcf_string(&self) -> Result<String> {
        let mut buf = htslib::kstring_t {
            l: 0,
            m: 0,
            s: ptr::null_mut(),
        };
        // SAFETY: self.header().inner and self.inner are non-null; buf is a valid kstring_t.
        let ret = unsafe { htslib::vcf_format(self.header().inner, self.inner, &mut buf) };

        if ret < 0 {
            if !buf.s.is_null() {
                // SAFETY: buf.s was allocated by vcf_format; free is symmetric.
                unsafe {
                    libc::free(buf.s as *mut libc::c_void);
                }
            }
            return Err(Error::ToString);
        }

        // SAFETY: buf.s is non-null after successful vcf_format; it is a valid NUL-terminated C string.
        let vcf_str = unsafe {
            let vcf_str = String::from(ffi::CStr::from_ptr(buf.s).to_str().unwrap());
            if !buf.s.is_null() {
                libc::free(buf.s as *mut libc::c_void);
            }
            vcf_str
        };

        Ok(vcf_str)
    }
}

impl Clone for Record {
    fn clone(&self) -> Self {
        // SAFETY: self.inner is non-null (from constructor); bcf_dup returns a new copy.
        let inner = unsafe { htslib::bcf_dup(self.inner) };
        Record {
            inner,
            header: self.header.clone(),
        }
    }
}

impl genome::AbstractLocus for Record {
    fn contig(&self) -> &str {
        str::from_utf8(
            self.header()
                .rid2name(self.rid().expect("rid not set"))
                .expect("unable to find rid in header"),
        )
        .expect("unable to interpret contig name as UTF-8")
    }

    fn pos(&self) -> u64 {
        self.pos() as u64
    }
}

/// Phased or unphased alleles, represented as indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GenotypeAllele {
    Unphased(i32),
    Phased(i32),
    UnphasedMissing,
    PhasedMissing,
}

impl GenotypeAllele {
    /// Decode given integer according to BCF standard.
    #[deprecated(
        since = "0.36.0",
        note = "Please use the conversion trait From<i32> for GenotypeAllele instead."
    )]
    pub fn from_encoded(encoded: i32) -> Self {
        match (encoded, encoded & 1) {
            (0, 0) => GenotypeAllele::UnphasedMissing,
            (1, 1) => GenotypeAllele::PhasedMissing,
            (e, 1) => GenotypeAllele::Phased((e >> 1) - 1),
            (e, 0) => GenotypeAllele::Unphased((e >> 1) - 1),
            _ => panic!("unexpected phasing type"),
        }
    }

    /// Get the index into the list of alleles.
    pub fn index(self) -> Option<u32> {
        match self {
            GenotypeAllele::Unphased(i) | GenotypeAllele::Phased(i) => Some(i as u32),
            GenotypeAllele::UnphasedMissing | GenotypeAllele::PhasedMissing => None,
        }
    }
}

impl fmt::Display for GenotypeAllele {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.index() {
            Some(a) => write!(f, "{}", a),
            None => write!(f, "."),
        }
    }
}

impl From<GenotypeAllele> for i32 {
    fn from(allele: GenotypeAllele) -> i32 {
        let (allele, phased) = match allele {
            GenotypeAllele::UnphasedMissing => (-1, 0),
            GenotypeAllele::PhasedMissing => (-1, 1),
            GenotypeAllele::Unphased(a) => (a, 0),
            GenotypeAllele::Phased(a) => (a, 1),
        };
        ((allele + 1) << 1) | phased
    }
}

impl From<i32> for GenotypeAllele {
    fn from(encoded: i32) -> GenotypeAllele {
        match (encoded, encoded & 1) {
            (0, 0) => GenotypeAllele::UnphasedMissing,
            (1, 1) => GenotypeAllele::PhasedMissing,
            (e, 1) => GenotypeAllele::Phased((e >> 1) - 1),
            (e, 0) => GenotypeAllele::Unphased((e >> 1) - 1),
            _ => panic!("unexpected phasing type"),
        }
    }
}

custom_derive! {
    /// Genotype representation as a vector of `GenotypeAllele`.
    #[derive(NewtypeDeref, Debug, Clone, PartialEq, Eq, Hash)]
    pub struct Genotype(Vec<GenotypeAllele>);
}

impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Genotype(alleles) = self;
        write!(f, "{}", alleles[0])?;
        for a in &alleles[1..] {
            let sep = match a {
                GenotypeAllele::Phased(_) | GenotypeAllele::PhasedMissing => '|',
                GenotypeAllele::Unphased(_) | GenotypeAllele::UnphasedMissing => '/',
            };
            write!(f, "{}{}", sep, a)?;
        }
        Ok(())
    }
}

/// Lazy representation of genotypes, that does no computation until a particular genotype is queried.
#[derive(Debug)]
pub struct Genotypes<'a, B>
where
    B: Borrow<Buffer> + 'a,
{
    encoded: BufferBacked<'a, Vec<&'a [i32]>, B>,
}

impl<'a, B: Borrow<Buffer> + 'a> Genotypes<'a, B> {
    /// Get genotype of ith sample.
    ///
    /// Note that the result complies with the BCF spec. This means that the
    /// first allele will always be marked as `Unphased`. That is, if you have 1|1 in the VCF,
    /// this method will return `[Unphased(1), Phased(1)]`.
    pub fn get(&self, i: usize) -> Genotype {
        let igt = self.encoded[i];
        let allelles = igt
            .iter()
            .take_while(|&&i| i != VECTOR_END_INTEGER)
            .map(|&i| GenotypeAllele::from(i))
            .collect();
        Genotype(allelles)
    }
}

impl Drop for Record {
    fn drop(&mut self) {
        // SAFETY: self.inner was allocated by bcf_init or bcf_dup; bcf_destroy is symmetric.
        unsafe { htslib::bcf_destroy(self.inner) };
    }
}

// SAFETY: Record owns its inner bcf1_t exclusively; header is Arc<HeaderView> which is Send+Sync.
unsafe impl Send for Record {}

// SAFETY: Record owns its inner bcf1_t exclusively; header is Arc<HeaderView> which is Send+Sync.
unsafe impl Sync for Record {}

/// Info tag representation.
#[derive(Debug)]
pub struct Info<'a, B: BorrowMut<Buffer> + Borrow<Buffer>> {
    record: &'a Record,
    tag: &'a [u8],
    buffer: B,
}

pub type BufferBackedOption<'b, B> = Option<BufferBacked<'b, Vec<&'b [u8]>, B>>;

impl<'b, B: BorrowMut<Buffer> + Borrow<Buffer> + 'b> Info<'_, B> {
    /// Short description of info tag.
    pub fn desc(&self) -> String {
        str::from_utf8(self.tag).unwrap().to_owned()
    }

    fn data(&mut self, data_type: u32) -> Result<Option<i32>> {
        let mut n: i32 = self.buffer.borrow().len;
        let c_str = ffi::CString::new(self.tag).unwrap();
        // SAFETY: record header and inner are non-null; c_str is valid; buffer.inner is managed by htslib.
        let ret = unsafe {
            htslib::bcf_get_info_values(
                self.record.header().inner,
                self.record.inner,
                c_str.as_ptr() as *mut c_char,
                &mut self.buffer.borrow_mut().inner,
                &mut n,
                data_type as i32,
            )
        };
        self.buffer.borrow_mut().len = n;

        match ret {
            -1 => Err(Error::UndefinedTag { tag: self.desc() }),
            -2 => Err(Error::UnexpectedType { tag: self.desc() }),
            -3 => Ok(None),
            ret => Ok(Some(ret)),
        }
    }

    /// Get integers from tag. `None` if tag not present in record.
    ///
    /// Import `bcf::record::Numeric` for missing value handling.
    ///
    /// **Attention:** the returned BufferBacked which holds the data has to be kept in scope
    /// as along as the data is accessed. If parts of the data are accessed while
    /// the BufferBacked object is already dropped, you will access unallocated
    /// memory.
    pub fn integer(mut self) -> Result<Option<BufferBacked<'b, &'b [i32], B>>> {
        self.data(htslib::BCF_HT_INT).map(|data| {
            data.map(|ret| {
                // SAFETY: buffer.inner was filled by bcf_get_info_values; ret is the element count.
                let values = unsafe {
                    slice::from_raw_parts(self.buffer.borrow().inner as *const i32, ret as usize)
                };
                BufferBacked::new(&values[..ret as usize], self.buffer)
            })
        })
    }

    /// Get floats from tag. `None` if tag not present in record.
    ///
    /// Import `bcf::record::Numeric` for missing value handling.
    ///
    /// **Attention:** the returned BufferBacked which holds the data has to be kept in scope
    /// as along as the data is accessed. If parts of the data are accessed while
    /// the BufferBacked object is already dropped, you will access unallocated
    /// memory.
    pub fn float(mut self) -> Result<Option<BufferBacked<'b, &'b [f32], B>>> {
        self.data(htslib::BCF_HT_REAL).map(|data| {
            data.map(|ret| {
                // SAFETY: buffer.inner was filled by bcf_get_info_values; ret is the element count.
                let values = unsafe {
                    slice::from_raw_parts(self.buffer.borrow().inner as *const f32, ret as usize)
                };
                BufferBacked::new(&values[..ret as usize], self.buffer)
            })
        })
    }

    /// Get flags from tag. `false` if not set.
    pub fn flag(&mut self) -> Result<bool> {
        self.data(htslib::BCF_HT_FLAG).map(|data| match data {
            Some(ret) => ret == 1,
            None => false,
        })
    }

    /// Get strings from tag. `None` if tag not present in record.
    ///
    /// **Attention:** the returned BufferBacked which holds the data has to be kept in scope
    /// as along as the data is accessed. If parts of the data are accessed while
    /// the BufferBacked object is already dropped, you will access unallocated
    /// memory.
    pub fn string(mut self) -> Result<BufferBackedOption<'b, B>> {
        self.data(htslib::BCF_HT_STR).map(|data| {
            data.map(|ret| {
                BufferBacked::new(
                    // SAFETY: buffer.inner was filled by bcf_get_info_values; ret is the byte count.
                    unsafe {
                        slice::from_raw_parts(self.buffer.borrow().inner as *const u8, ret as usize)
                    }
                    .split(|c| *c == b',')
                    .map(|s| {
                        // stop at zero character
                        s.split(|c| *c == 0u8)
                            .next()
                            .expect("Bug: returned string should not be empty.")
                    })
                    .collect(),
                    self.buffer,
                )
            })
        })
    }
}

// SAFETY: Info borrows a Record (Send+Sync) and owns a Buffer; no shared mutable state across threads.
unsafe impl<B: BorrowMut<Buffer> + Borrow<Buffer>> Send for Info<'_, B> {}

// SAFETY: Info borrows a Record (Send+Sync) and owns a Buffer; no shared mutable state across threads.
unsafe impl<B: BorrowMut<Buffer> + Borrow<Buffer>> Sync for Info<'_, B> {}

fn trim_slice<T: PartialEq + NumericUtils>(s: &[T]) -> &[T] {
    s.split(|v| v.is_vector_end())
        .next()
        .expect("Bug: returned slice should not be empty.")
}

// Representation of per-sample data.
#[derive(Debug)]
pub struct Format<'a, B: BorrowMut<Buffer> + Borrow<Buffer>> {
    record: &'a Record,
    tag: &'a [u8],
    inner: *mut htslib::bcf_fmt_t,
    buffer: B,
}

impl<'a, 'b, B: BorrowMut<Buffer> + Borrow<Buffer> + 'b> Format<'a, B> {
    /// Create new format data in a given record.
    fn new(record: &'a Record, tag: &'a [u8], buffer: B) -> Format<'a, B> {
        let c_str = ffi::CString::new(tag).unwrap();
        // SAFETY: record header and inner are non-null (from constructor); c_str is a valid CString.
        let inner = unsafe {
            htslib::bcf_get_fmt(
                record.header().inner,
                record.inner,
                c_str.as_ptr() as *mut c_char,
            )
        };
        Format {
            record,
            tag,
            inner,
            buffer,
        }
    }

    /// Provide short description of format entry (just the tag name).
    pub fn desc(&self) -> String {
        str::from_utf8(self.tag).unwrap().to_owned()
    }

    pub fn inner(&self) -> &htslib::bcf_fmt_t {
        // SAFETY: self.inner is set by bcf_get_fmt which returns a valid pointer (or null handled by caller).
        unsafe { &*self.inner }
    }

    pub fn inner_mut(&mut self) -> &mut htslib::bcf_fmt_t {
        // SAFETY: self.inner is set by bcf_get_fmt; we have &mut self.
        unsafe { &mut *self.inner }
    }

    fn values_per_sample(&self) -> usize {
        self.inner().n as usize
    }

    /// Read and decode format data into a given type.
    fn data(&mut self, data_type: u32) -> Result<i32> {
        let mut n: i32 = self.buffer.borrow().len;
        let c_str = ffi::CString::new(self.tag).unwrap();
        // SAFETY: record header and inner are non-null; c_str is valid; buffer.inner is managed by htslib.
        let ret = unsafe {
            htslib::bcf_get_format_values(
                self.record.header().inner,
                self.record.inner,
                c_str.as_ptr() as *mut c_char,
                &mut self.buffer.borrow_mut().inner,
                &mut n,
                data_type as i32,
            )
        };
        self.buffer.borrow_mut().len = n;
        match ret {
            -1 => Err(Error::UndefinedTag { tag: self.desc() }),
            -2 => Err(Error::UnexpectedType { tag: self.desc() }),
            -3 => Err(Error::MissingTag {
                tag: self.desc(),
                record: self.record.desc(),
            }),
            ret => Ok(ret),
        }
    }

    /// Get format data as integers.
    ///
    /// **Attention:** the returned BufferBacked which holds the data has to be kept in scope
    /// as long as the data is accessed. If parts of the data are accessed while
    /// the BufferBacked object is already dropped, you will access unallocated
    /// memory.
    pub fn integer(mut self) -> Result<BufferBacked<'b, Vec<&'b [i32]>, B>> {
        self.data(htslib::BCF_HT_INT).map(|ret| {
            BufferBacked::new(
                // SAFETY: buffer.inner was filled by bcf_get_format_values; ret is the element count.
                unsafe {
                    slice::from_raw_parts(
                        self.buffer.borrow_mut().inner as *const i32,
                        ret as usize,
                    )
                }
                .chunks(self.values_per_sample())
                .map(trim_slice)
                .collect(),
                self.buffer,
            )
        })
    }

    /// Get format data as floats.
    ///
    /// **Attention:** the returned BufferBacked which holds the data has to be kept in scope
    /// as along as the data is accessed. If parts of the data are accessed while
    /// the BufferBacked object is already dropped, you will access unallocated
    /// memory.
    pub fn float(mut self) -> Result<BufferBacked<'b, Vec<&'b [f32]>, B>> {
        self.data(htslib::BCF_HT_REAL).map(|ret| {
            BufferBacked::new(
                // SAFETY: buffer.inner was filled by bcf_get_format_values; ret is the element count.
                unsafe {
                    slice::from_raw_parts(
                        self.buffer.borrow_mut().inner as *const f32,
                        ret as usize,
                    )
                }
                .chunks(self.values_per_sample())
                .map(trim_slice)
                .collect(),
                self.buffer,
            )
        })
    }

    /// Get format data as byte slices. To obtain the values strings, use `std::str::from_utf8`.
    ///
    /// **Attention:** the returned BufferBacked which holds the data has to be kept in scope
    /// as along as the data is accessed. If parts of the data are accessed while
    /// the BufferBacked object is already dropped, you will access unallocated
    /// memory.
    pub fn string(mut self) -> Result<BufferBacked<'b, Vec<&'b [u8]>, B>> {
        self.data(htslib::BCF_HT_STR).map(|ret| {
            if ret == 0 {
                return BufferBacked::new(Vec::new(), self.buffer);
            }
            BufferBacked::new(
                // SAFETY: buffer.inner was filled by bcf_get_format_values; ret is the byte count.
                unsafe {
                    slice::from_raw_parts(self.buffer.borrow_mut().inner as *const u8, ret as usize)
                }
                .chunks(self.values_per_sample())
                .map(|s| {
                    // stop at zero character
                    s.split(|c| *c == 0u8)
                        .next()
                        .expect("Bug: returned string should not be empty.")
                })
                .collect(),
                self.buffer,
            )
        })
    }
}

// SAFETY: Format borrows a Record (Send+Sync) and owns a Buffer; no shared mutable state across threads.
unsafe impl<B: BorrowMut<Buffer> + Borrow<Buffer>> Send for Format<'_, B> {}

// SAFETY: Format borrows a Record (Send+Sync) and owns a Buffer; no shared mutable state across threads.
unsafe impl<B: BorrowMut<Buffer> + Borrow<Buffer>> Sync for Format<'_, B> {}

#[derive(Debug)]
pub struct Filters<'a> {
    /// Reference to the `Record` to enumerate records for.
    record: &'a Record,
    /// Index of the next filter to return, if not at end.
    idx: i32,
}

impl<'a> Filters<'a> {
    pub fn new(record: &'a Record) -> Self {
        Filters { record, idx: 0 }
    }
}

impl Iterator for Filters<'_> {
    type Item = Id;

    fn next(&mut self) -> Option<Id> {
        if self.record.inner().d.n_flt <= self.idx {
            None
        } else {
            let i = self.idx as isize;
            self.idx += 1;
            // SAFETY: i is within [0, n_flt); d.flt is a valid pointer with n_flt elements.
            Some(Id(unsafe { *self.record.inner().d.flt.offset(i) } as u32))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bcf::{Format, Header, Writer};
    use tempfile::NamedTempFile;

    #[test]
    fn test_missing_float() {
        let expected: u32 = 0x7F80_0001;
        assert_eq!(MISSING_FLOAT.bits(), expected);
    }

    #[test]
    fn test_vector_end_float() {
        let expected: u32 = 0x7F80_0002;
        assert_eq!(VECTOR_END_FLOAT.bits(), expected);
    }

    #[test]
    fn test_record_rlen() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        assert_eq!(record.rlen(), 0);
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        assert_eq!(record.rlen(), 3)
    }

    #[test]
    fn test_record_end() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record.set_pos(5);

        assert_eq!(record.end(), 8)
    }

    #[test]
    fn test_record_clear() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        header.push_sample("sample".as_bytes());
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record.set_pos(6);
        record.clear();

        assert_eq!(record.rlen(), 0);
        assert_eq!(record.sample_count(), 0);
        assert_eq!(record.pos(), 0)
    }

    #[test]
    fn test_record_clone() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record.set_pos(6);

        let mut cloned_record = record.clone();
        cloned_record.set_pos(5);

        assert_eq!(record.pos(), 6);
        assert_eq!(record.allele_count(), 2);
        assert_eq!(cloned_record.pos(), 5);
        assert_eq!(cloned_record.allele_count(), 2);
    }

    #[test]
    fn test_record_has_filter_pass_is_default() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let record = vcf.empty_record();

        assert!(record.has_filter("PASS".as_bytes()));
        assert!(record.has_filter(".".as_bytes()));
        assert!(record.has_filter(&Id(0)));
        assert!(!record.has_filter("foo".as_bytes()));
        assert!(!record.has_filter(&Id(2)));
    }

    #[test]
    fn test_record_has_filter_custom() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        record.push_filter("foo".as_bytes()).unwrap();

        assert!(record.has_filter("foo".as_bytes()));
        assert!(!record.has_filter("PASS".as_bytes()))
    }

    #[test]
    fn test_record_push_filter() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
        header.push_record(br#"##FILTER=<ID=bar,Description="dranks">"#);
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        assert!(record.has_filter("PASS".as_bytes()));
        record.push_filter("foo".as_bytes()).unwrap();
        let bar = record.header().name_to_id(cstr8!("bar")).unwrap();
        record.push_filter(&bar).unwrap();
        assert!(record.has_filter("foo".as_bytes()));
        assert!(record.has_filter(&bar));
        assert!(!record.has_filter("PASS".as_bytes()));
        assert!(record.push_filter("baz".as_bytes()).is_err())
    }

    #[test]
    fn test_record_set_filters() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
        header.push_record(br#"##FILTER=<ID=bar,Description="a horse walks into...">"#);
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        assert!(record.has_filter("PASS".as_bytes()));
        record
            .set_filters(&["foo".as_bytes(), "bar".as_bytes()])
            .unwrap();
        assert!(record.has_filter("foo".as_bytes()));
        assert!(record.has_filter("bar".as_bytes()));
        assert!(!record.has_filter("PASS".as_bytes()));
        let filters: &[&Id] = &[];
        record.set_filters(filters).unwrap();
        assert!(record.has_filter("PASS".as_bytes()));
        assert!(!record.has_filter("foo".as_bytes()));
        assert!(record
            .set_filters(&["foo".as_bytes(), "baz".as_bytes()])
            .is_err())
    }

    #[test]
    fn test_record_remove_filter() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
        header.push_record(br#"##FILTER=<ID=bar,Description="a horse walks into...">"#);
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let foo = record.header().name_to_id(cstr8!("foo")).unwrap();
        let bar = record.header().name_to_id(cstr8!("bar")).unwrap();
        record.set_filters(&[&foo, &bar]).unwrap();
        assert!(record.has_filter(&foo));
        assert!(record.has_filter(&bar));
        record.remove_filter(&foo, true).unwrap();
        assert!(!record.has_filter(&foo));
        assert!(record.has_filter(&bar));
        assert!(record.remove_filter("baz".as_bytes(), true).is_err());
        record.remove_filter(&bar, true).unwrap();
        assert!(!record.has_filter(&bar));
        assert!(record.has_filter("PASS".as_bytes()));
    }

    #[test]
    fn test_record_to_vcf_string_err() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let record = vcf.empty_record();
        assert!(record.to_vcf_string().is_err());
    }

    #[test]
    fn test_record_to_vcf_string() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        header.push_record(b"##contig=<ID=chr1,length=1000>");
        header.push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        record.push_filter("foo".as_bytes()).unwrap();
        assert_eq!(
            record.to_vcf_string().unwrap(),
            "chr1\t1\t.\t.\t.\t0\tfoo\t.\n"
        );
    }
}

#[cfg(test)]
mod bcf_unpack_tests {
    use super::*;
    use crate::bcf::{Format, Header, Writer};
    use tempfile::NamedTempFile;

    /// Helper: create a BCF writer with a rich header (contigs, filters, INFO, FORMAT, samples).
    fn make_writer_and_header() -> (NamedTempFile, Writer) {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        header.push_record(b"##contig=<ID=chr1,length=10000>");
        header.push_record(b"##contig=<ID=chr2,length=20000>");
        header.push_record(br#"##FILTER=<ID=q10,Description="Quality below 10">"#);
        header.push_record(br#"##FILTER=<ID=s50,Description="SNP cluster">"#);
        header.push_record(br#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">"#);
        header.push_record(br#"##INFO=<ID=AF,Number=A,Type=Float,Description="AF">"#);
        header.push_record(br#"##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">"#);
        header.push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
        header.push_record(br#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">"#);
        header.push_sample(b"sample1");
        let writer = Writer::from_path(path, &header, true, Format::Bcf).unwrap();
        (tmp, writer)
    }

    /// Compare the decoded state produced by C bcf_unpack vs Rust bcf_unpack_rs.
    ///
    /// Strategy:
    /// 1. Record was already unpacked by C (via writer.empty_record() + set_* methods)
    /// 2. Read out the C-decoded values as the oracle
    /// 3. Reset the unpacked flag to 0
    /// 4. Call bcf_unpack_rs
    /// 5. Read out decoded values again and compare
    unsafe fn assert_unpack_equivalence(rec: &mut Record) {
        // Step 1: Record is already C-unpacked. Read oracle values.
        let c_id = rec.id();
        let c_alleles: Vec<Vec<u8>> = rec.alleles().iter().map(|a| a.to_vec()).collect();
        let c_filters: Vec<u32> = rec.filters().map(|f| f.0).collect();

        // Also capture INFO state
        let inner = &*rec.inner;
        let n_info = inner.n_info() as usize;
        let mut c_info_keys = Vec::new();
        let mut c_info_types = Vec::new();
        let mut c_info_lens = Vec::new();
        let mut c_info_v1 = Vec::new();
        for i in 0..n_info {
            let info = &*inner.d.info.add(i);
            c_info_keys.push(info.key);
            c_info_types.push(info.type_);
            c_info_lens.push(info.len);
            // Capture scalar v1 value (i64 covers both int and float via bits)
            c_info_v1.push(info.v1.i);
        }

        // Capture FORMAT state
        let n_fmt = inner.n_fmt() as usize;
        let mut c_fmt_ids = Vec::new();
        let mut c_fmt_ns = Vec::new();
        let mut c_fmt_types = Vec::new();
        let mut c_fmt_sizes = Vec::new();
        for i in 0..n_fmt {
            let fmt = &*inner.d.fmt.add(i);
            c_fmt_ids.push(fmt.id);
            c_fmt_ns.push(fmt.n);
            c_fmt_types.push(fmt.type_);
            c_fmt_sizes.push(fmt.size);
        }

        // Step 2: Reset unpacked flag and unpack_size to force re-unpacking
        (*rec.inner).unpacked = 0;
        (*rec.inner).unpack_size = [0; 3];

        // Step 3: Call Rust implementation
        bcf_unpack_rs(rec.inner, htslib::BCF_UN_ALL as i32);

        // Step 4: Compare decoded values
        let rs_id = rec.id();
        assert_eq!(c_id, rs_id, "ID mismatch");

        let rs_alleles: Vec<Vec<u8>> = rec.alleles().iter().map(|a| a.to_vec()).collect();
        assert_eq!(c_alleles, rs_alleles, "Alleles mismatch");

        let rs_filters: Vec<u32> = rec.filters().map(|f| f.0).collect();
        assert_eq!(c_filters, rs_filters, "Filters mismatch");

        // Compare INFO
        let inner = &*rec.inner;
        for i in 0..n_info {
            let info = &*inner.d.info.add(i);
            assert_eq!(c_info_keys[i], info.key, "INFO key mismatch at {i}");
            assert_eq!(c_info_types[i], info.type_, "INFO type mismatch at {i}");
            assert_eq!(c_info_lens[i], info.len, "INFO len mismatch at {i}");
            assert_eq!(c_info_v1[i], info.v1.i, "INFO v1 scalar mismatch at {i}");
        }

        // Compare FORMAT
        for i in 0..n_fmt {
            let fmt = &*inner.d.fmt.add(i);
            assert_eq!(c_fmt_ids[i], fmt.id, "FMT id mismatch at {i}");
            assert_eq!(c_fmt_ns[i], fmt.n, "FMT n mismatch at {i}");
            assert_eq!(c_fmt_types[i], fmt.type_, "FMT type mismatch at {i}");
            assert_eq!(c_fmt_sizes[i], fmt.size, "FMT size mismatch at {i}");
        }
    }

    #[test]
    fn test_unpack_minimal_record() {
        let (_tmp, writer) = make_writer_and_header();
        let mut record = writer.empty_record();
        record.set_rid(Some(0));
        record.set_pos(100);
        record.set_alleles(&[b"A", b"T"]).unwrap();
        unsafe { assert_unpack_equivalence(&mut record) };
    }

    #[test]
    fn test_unpack_with_id() {
        let (_tmp, writer) = make_writer_and_header();
        let mut record = writer.empty_record();
        record.set_rid(Some(0));
        record.set_pos(42);
        record.set_id(b"rs12345").unwrap();
        record.set_alleles(&[b"C", b"G"]).unwrap();
        unsafe { assert_unpack_equivalence(&mut record) };
    }

    #[test]
    fn test_unpack_with_filters() {
        let (_tmp, writer) = make_writer_and_header();
        let mut record = writer.empty_record();
        record.set_rid(Some(0));
        record.set_pos(200);
        record.set_alleles(&[b"A", b"T"]).unwrap();
        record.push_filter("q10".as_bytes()).unwrap();
        record.push_filter("s50".as_bytes()).unwrap();
        unsafe { assert_unpack_equivalence(&mut record) };
    }

    #[test]
    fn test_unpack_multiallelic() {
        let (_tmp, writer) = make_writer_and_header();
        let mut record = writer.empty_record();
        record.set_rid(Some(1));
        record.set_pos(500);
        record
            .set_alleles(&[b"ACGT", b"A", b"AC", b"ACGTACGT"])
            .unwrap();
        unsafe { assert_unpack_equivalence(&mut record) };
    }

    #[test]
    fn test_unpack_with_info_fields() {
        let (_tmp, writer) = make_writer_and_header();
        let mut record = writer.empty_record();
        record.set_rid(Some(0));
        record.set_pos(300);
        record.set_alleles(&[b"G", b"A"]).unwrap();
        record.push_info_integer(cstr8!("DP"), &[42]).unwrap();
        record.push_info_float(cstr8!("AF"), &[0.5]).unwrap();
        record.push_info_flag(cstr8!("DB")).unwrap();
        unsafe { assert_unpack_equivalence(&mut record) };
    }

    #[test]
    fn test_unpack_with_format_fields() {
        let (_tmp, writer) = make_writer_and_header();
        let mut record = writer.empty_record();
        record.set_rid(Some(0));
        record.set_pos(400);
        record.set_alleles(&[b"T", b"C"]).unwrap();
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
            .unwrap();
        record.push_format_integer(cstr8!("GQ"), &[30]).unwrap();
        unsafe { assert_unpack_equivalence(&mut record) };
    }

    #[test]
    fn test_unpack_kitchen_sink() {
        let (_tmp, writer) = make_writer_and_header();
        let mut record = writer.empty_record();
        record.set_rid(Some(1));
        record.set_pos(999);
        record.set_id(b"rs99999").unwrap();
        record.set_alleles(&[b"AAAA", b"T", b"GGGG"]).unwrap();
        record.set_qual(99.9);
        record.push_filter("q10".as_bytes()).unwrap();
        record.push_info_integer(cstr8!("DP"), &[1000]).unwrap();
        record.push_info_float(cstr8!("AF"), &[0.25]).unwrap();
        record.push_info_flag(cstr8!("DB")).unwrap();
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1), GenotypeAllele::Phased(2)])
            .unwrap();
        record.push_format_integer(cstr8!("GQ"), &[99]).unwrap();
        unsafe { assert_unpack_equivalence(&mut record) };
    }

    #[test]
    fn test_unpack_from_file() {
        // Test with real BCF data read from disk
        use crate::bcf::Read;
        let mut reader =
            crate::bcf::Reader::from_path("test/test_multi.bcf").expect("Error opening file.");
        let mut record = reader.empty_record();
        while let Some(Ok(())) = reader.read(&mut record) {
            // Record was unpacked by C via read(). Now test Rust re-unpack.
            unsafe { assert_unpack_equivalence(&mut record) };
        }
    }
}

#[cfg(test)]
mod kbitset_tests {
    use super::*;

    /// Build a kbitset_t using the C FFI (kbs_init + kbs_insert) as oracle.
    /// Returns a raw pointer that must be freed with kbs_destroy.
    unsafe fn kbitset_from_bools_c(bits: &[bool]) -> *mut htslib::kbitset_t {
        let bs = htslib::kbs_init(bits.len());
        for (i, &set) in bits.iter().enumerate() {
            if set {
                htslib::kbs_insert(bs, i as i32);
            }
        }
        bs
    }

    /// Compare C and Rust kbitset_t structs for identical bit patterns.
    unsafe fn assert_kbitset_eq(c_bs: *const htslib::kbitset_t, rs_bs: &KBitSet, label: &str) {
        let c = &*c_bs;
        let r = &*rs_bs.as_ptr();
        assert_eq!(c.n, r.n, "{label}: n mismatch");
        assert_eq!(c.n_max, r.n_max, "{label}: n_max mismatch");
        // Compare all bit slots including sentinel
        for i in 0..=c.n {
            let c_val = *c.b.as_ptr().add(i);
            let r_val = *r.b.as_ptr().add(i);
            assert_eq!(c_val, r_val, "{label}: b[{i}] mismatch");
        }
    }

    #[test]
    fn kbitset_matches_c_empty() {
        let bits = vec![false; 10];
        let rs = KBitSet::from_bools(&bits);
        unsafe {
            let c = kbitset_from_bools_c(&bits);
            assert_kbitset_eq(c, &rs, "empty 10-bit");
            htslib::kbs_destroy(c);
        }
    }

    #[test]
    fn kbitset_matches_c_all_set() {
        let bits = vec![true; 10];
        let rs = KBitSet::from_bools(&bits);
        unsafe {
            let c = kbitset_from_bools_c(&bits);
            assert_kbitset_eq(c, &rs, "all-set 10-bit");
            htslib::kbs_destroy(c);
        }
    }

    #[test]
    fn kbitset_matches_c_various_sizes() {
        // Test sizes that cross ulong boundaries (64-bit)
        for size in [0, 1, 2, 7, 8, 31, 32, 33, 63, 64, 65, 127, 128, 129, 256] {
            let bits: Vec<bool> = (0..size).map(|i| i % 3 == 0 || i % 7 == 0).collect();
            let rs = KBitSet::from_bools(&bits);
            unsafe {
                let c = kbitset_from_bools_c(&bits);
                assert_kbitset_eq(c, &rs, &format!("pattern size={size}"));
                htslib::kbs_destroy(c);
            }
        }
    }

    #[test]
    fn kbitset_matches_c_single_bits() {
        for size in [1, 64, 65, 128] {
            for bit in [0, size / 2, size - 1] {
                let mut bits = vec![false; size];
                bits[bit] = true;
                let rs = KBitSet::from_bools(&bits);
                unsafe {
                    let c = kbitset_from_bools_c(&bits);
                    assert_kbitset_eq(c, &rs, &format!("single bit={bit} size={size}"));
                    htslib::kbs_destroy(c);
                }
            }
        }
    }

    #[test]
    fn remove_alleles_matches_existing_test() {
        // Reproduce the existing test_remove_alleles from mod.rs
        use crate::bcf::Read;
        let mut bcf = crate::bcf::Reader::from_path("test/test_multi.bcf").unwrap();
        for res in bcf.records() {
            let mut record = res.unwrap();
            if record.pos() == 10080 {
                record.remove_alleles(&[false, false, true]).unwrap();
                assert_eq!(record.alleles(), [b"A", b"C"]);
            }
        }
    }
}
