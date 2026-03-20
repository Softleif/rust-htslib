// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::convert::TryFrom;
use std::convert::TryInto;
use std::ffi;
use std::fmt;
use std::marker::PhantomData;
use std::mem::{size_of, MaybeUninit};
use std::ops;
use std::os::raw::c_char;
use std::slice;
use std::str;
use std::sync::Arc;

use byteorder::{LittleEndian, ReadBytesExt};

use crate::bam::BamError as Error;
use crate::bam::HeaderView;

type Result<T> = std::result::Result<T, Error>;
use crate::htslib;
use crate::utils;
#[cfg(feature = "serde_feature")]
use serde::{self, Deserialize, Serialize};

use bio_types::alignment::{Alignment, AlignmentMode, AlignmentOperation};
use bio_types::genome;
use bio_types::sequence::SequenceRead;
use bio_types::sequence::SequenceReadPairOrientation;
use bio_types::strand::ReqStrand;

/// A macro creating methods for flag access.
macro_rules! flag {
    ($get:ident, $set:ident, $unset:ident, $bit:expr) => {
        pub fn $get(&self) -> bool {
            self.inner().core.flag & $bit != 0
        }

        pub fn $set(&mut self) {
            self.inner_mut().core.flag |= $bit;
        }

        pub fn $unset(&mut self) {
            self.inner_mut().core.flag &= !$bit;
        }
    };
}

/// A BAM record.
pub struct Record {
    pub inner: htslib::bam1_t,
    own: bool,
    cigar: Option<CigarStringView>,
    header: Option<Arc<HeaderView>>,
}

// SAFETY: Record owns its bam1_t data buffer; no interior mutability or shared references to FFI state.
unsafe impl Send for Record {}
unsafe impl Sync for Record {}

impl Clone for Record {
    fn clone(&self) -> Self {
        let mut copy = Record::new();
        // SAFETY: both pointers are valid bam1_t structs; bam_copy1 deep-copies all data.
        unsafe { htslib::bam_copy1(copy.inner_ptr_mut(), self.inner_ptr()) };
        copy
    }
}

impl PartialEq for Record {
    fn eq(&self, other: &Record) -> bool {
        self.tid() == other.tid()
            && self.pos() == other.pos()
            && self.bin() == other.bin()
            && self.mapq() == other.mapq()
            && self.flags() == other.flags()
            && self.mtid() == other.mtid()
            && self.mpos() == other.mpos()
            && self.insert_size() == other.insert_size()
            && self.data() == other.data()
            && self.inner().core.l_extranul == other.inner().core.l_extranul
    }
}

impl Eq for Record {}

impl fmt::Debug for Record {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.write_fmt(format_args!(
            "Record(tid: {}, pos: {})",
            self.tid(),
            self.pos()
        ))
    }
}

impl Default for Record {
    fn default() -> Self {
        Self::new()
    }
}

#[inline]
fn extranul_from_qname(qname: &[u8]) -> usize {
    let qlen = qname.len() + 1;
    if !qlen.is_multiple_of(4) {
        4 - qlen % 4
    } else {
        0
    }
}

impl Record {
    /// Create an empty BAM record.
    pub fn new() -> Self {
        let mut record = Record {
            // SAFETY: bam1_t is a C struct where all-zeros is a valid initial state (null data pointer, zero lengths).
            inner: unsafe { MaybeUninit::zeroed().assume_init() },
            own: true,
            cigar: None,
            header: None,
        };
        // The read/query name needs to be set as empty to properly initialize
        // the record
        record.set_qname(b"");
        // Developer note: these are needed so the returned record is properly
        // initialized as unmapped.
        record.set_unmapped();
        record.set_tid(-1);
        record.set_pos(-1);
        record.set_mpos(-1);
        record.set_mtid(-1);
        record
    }

    pub fn from_inner(from: *mut htslib::bam1_t) -> Self {
        Record {
            inner: {
                // SAFETY: inner is immediately overwritten by memcpy below; never read in uninitialized state.
                #[allow(clippy::uninit_assumed_init, invalid_value)]
                let mut inner = unsafe { MaybeUninit::uninit().assume_init() };
                // SAFETY: from must be a valid bam1_t pointer (caller invariant); memcpy copies the full struct.
                unsafe {
                    ::libc::memcpy(
                        &mut inner as *mut htslib::bam1_t as *mut ::libc::c_void,
                        from as *const ::libc::c_void,
                        size_of::<htslib::bam1_t>(),
                    );
                }
                inner
            },
            own: false,
            cigar: None,
            header: None,
        }
    }

    // Create a BAM record from a line SAM text. SAM slice need not be 0-terminated.
    pub fn from_sam(header_view: &HeaderView, sam: &[u8]) -> Result<Record> {
        let mut record = Self::new();

        let mut sam_copy = Vec::with_capacity(sam.len() + 1);
        sam_copy.extend(sam);
        sam_copy.push(0);

        let mut sam_string = htslib::kstring_t {
            s: sam_copy.as_ptr() as *mut c_char,
            l: sam_copy.len(),
            m: sam_copy.len(),
        };

        // SAFETY: sam_string points to a valid NUL-terminated copy; header and record pointers are valid.
        let succ = unsafe {
            htslib::sam_parse1(
                &mut sam_string,
                header_view.inner_ptr() as *mut htslib::bam_hdr_t,
                record.inner_ptr_mut(),
            )
        };

        if succ == 0 {
            Ok(record)
        } else {
            Err(Error::ParseSAM {
                rec: str::from_utf8(&sam_copy).unwrap().to_owned(),
            })
        }
    }

    pub fn set_header(&mut self, header: Arc<HeaderView>) {
        self.header = Some(header);
    }

    pub(super) fn data(&self) -> &[u8] {
        // SAFETY: inner.data is valid for inner.l_data bytes (maintained by htslib and our set_data methods).
        unsafe { slice::from_raw_parts(self.inner().data, self.inner().l_data as usize) }
    }

    #[inline]
    pub fn inner_mut(&mut self) -> &mut htslib::bam1_t {
        &mut self.inner
    }

    #[inline]
    pub(super) fn inner_ptr_mut(&mut self) -> *mut htslib::bam1_t {
        &mut self.inner as *mut htslib::bam1_t
    }

    #[inline]
    pub fn inner(&self) -> &htslib::bam1_t {
        &self.inner
    }

    #[inline]
    pub(super) fn inner_ptr(&self) -> *const htslib::bam1_t {
        &self.inner as *const htslib::bam1_t
    }

    /// Get target id.
    pub fn tid(&self) -> i32 {
        self.inner().core.tid
    }

    /// Set target id.
    pub fn set_tid(&mut self, tid: i32) {
        self.inner_mut().core.tid = tid;
    }

    /// Get position (0-based).
    pub fn pos(&self) -> i64 {
        self.inner().core.pos
    }

    /// Set position (0-based).
    pub fn set_pos(&mut self, pos: i64) {
        self.inner_mut().core.pos = pos;
    }

    pub fn bin(&self) -> u16 {
        self.inner().core.bin
    }

    pub fn set_bin(&mut self, bin: u16) {
        self.inner_mut().core.bin = bin;
    }

    /// Get MAPQ.
    pub fn mapq(&self) -> u8 {
        self.inner().core.qual
    }

    /// Set MAPQ.
    pub fn set_mapq(&mut self, mapq: u8) {
        self.inner_mut().core.qual = mapq;
    }

    /// Get strand information from record flags.
    pub fn strand(&self) -> ReqStrand {
        let reverse = self.flags() & 0x10 != 0;
        if reverse {
            ReqStrand::Reverse
        } else {
            ReqStrand::Forward
        }
    }

    /// Get raw flags.
    pub fn flags(&self) -> u16 {
        self.inner().core.flag
    }

    /// Set raw flags.
    pub fn set_flags(&mut self, flags: u16) {
        self.inner_mut().core.flag = flags;
    }

    /// Unset all flags.
    pub fn unset_flags(&mut self) {
        self.inner_mut().core.flag = 0;
    }

    /// Get target id of mate.
    pub fn mtid(&self) -> i32 {
        self.inner().core.mtid
    }

    /// Set target id of mate.
    pub fn set_mtid(&mut self, mtid: i32) {
        self.inner_mut().core.mtid = mtid;
    }

    /// Get mate position.
    pub fn mpos(&self) -> i64 {
        self.inner().core.mpos
    }

    /// Set mate position.
    pub fn set_mpos(&mut self, mpos: i64) {
        self.inner_mut().core.mpos = mpos;
    }

    /// Get insert size.
    pub fn insert_size(&self) -> i64 {
        self.inner().core.isize_
    }

    /// Set insert size.
    pub fn set_insert_size(&mut self, insert_size: i64) {
        self.inner_mut().core.isize_ = insert_size;
    }

    fn qname_capacity(&self) -> usize {
        self.inner().core.l_qname as usize
    }

    fn qname_len(&self) -> usize {
        // discount all trailing zeros (the default one and extra nulls)
        self.qname_capacity() - 1 - self.inner().core.l_extranul as usize
    }

    /// Get qname (read name). Complexity: O(1).
    pub fn qname(&self) -> &[u8] {
        &self.data()[..self.qname_len()]
    }

    /// Set the variable length data buffer
    pub fn set_data(&mut self, new_data: &[u8]) {
        self.cigar = None;

        self.inner_mut().l_data = new_data.len() as i32;
        if (self.inner().m_data as i32) < self.inner().l_data {
            // Verbosity due to lexical borrowing
            let l_data = self.inner().l_data;
            self.realloc_var_data(l_data as usize);
        }

        // Copy new data into buffer
        // SAFETY: inner.data is valid for at least l_data bytes (realloc'd above if needed).
        let data =
            unsafe { slice::from_raw_parts_mut(self.inner.data, self.inner().l_data as usize) };
        utils::copy_memory(new_data, data);
    }

    /// Set variable length data (qname, cigar, seq, qual).
    /// The aux data is left unchanged.
    /// `qual` is Phred-scaled quality values, without any offset.
    /// NOTE: seq.len() must equal qual.len() or this method
    /// will panic. If you don't have quality values use
    /// `let quals = vec![ 255 as u8; seq.len()];` as a placeholder that will
    /// be recognized as missing QVs by `samtools`.
    pub fn set(&mut self, qname: &[u8], cigar: Option<&CigarString>, seq: &[u8], qual: &[u8]) {
        assert!(qname.len() < 255);
        assert_eq!(seq.len(), qual.len(), "seq.len() must equal qual.len()");

        self.cigar = None;

        let cigar_width = if let Some(cigar_string) = cigar {
            cigar_string.len()
        } else {
            0
        } * 4;
        let q_len = qname.len() + 1;
        let extranul = extranul_from_qname(qname);

        let orig_aux_offset = self.qname_capacity()
            + 4 * self.cigar_len()
            + self.seq_len().div_ceil(2)
            + self.seq_len();
        let new_aux_offset = q_len + extranul + cigar_width + seq.len().div_ceil(2) + qual.len();
        assert!(orig_aux_offset <= self.inner.l_data as usize);
        let aux_len = self.inner.l_data as usize - orig_aux_offset;
        self.inner_mut().l_data = (new_aux_offset + aux_len) as i32;
        if (self.inner().m_data as i32) < self.inner().l_data {
            // Verbosity due to lexical borrowing
            let l_data = self.inner().l_data;
            self.realloc_var_data(l_data as usize);
        }

        // Copy the aux data.
        if aux_len > 0 && orig_aux_offset != new_aux_offset {
            // SAFETY: inner.data is valid for m_data bytes; offsets are within bounds.
            let data =
                unsafe { slice::from_raw_parts_mut(self.inner.data, self.inner().m_data as usize) };
            data.copy_within(orig_aux_offset..orig_aux_offset + aux_len, new_aux_offset);
        }

        // SAFETY: inner.data is valid for l_data bytes (realloc'd above if needed).
        let data =
            unsafe { slice::from_raw_parts_mut(self.inner.data, self.inner().l_data as usize) };

        // qname
        utils::copy_memory(qname, data);
        for i in 0..=extranul {
            data[qname.len() + i] = b'\0';
        }
        let mut i = q_len + extranul;
        self.inner_mut().core.l_qname = i as u16;
        self.inner_mut().core.l_extranul = extranul as u8;

        // cigar
        if let Some(cigar_string) = cigar {
            // SAFETY: cigar offset is always 4-byte aligned (qname padded with extranul); length from cigar_string.
            let cigar_data = unsafe {
                //cigar is always aligned to 4 bytes (see extranul above) - so this is safe
                #[allow(clippy::cast_ptr_alignment)]
                slice::from_raw_parts_mut(data[i..].as_ptr() as *mut u32, cigar_string.len())
            };
            for (i, c) in cigar_string.iter().enumerate() {
                cigar_data[i] = c.encode();
            }
            self.inner_mut().core.n_cigar = cigar_string.len() as u32;
            i += cigar_string.len() * 4;
        } else {
            self.inner_mut().core.n_cigar = 0;
        };

        // seq
        {
            for j in (0..seq.len()).step_by(2) {
                data[i + j / 2] = (ENCODE_BASE[seq[j] as usize] << 4)
                    | (if j + 1 < seq.len() {
                        ENCODE_BASE[seq[j + 1] as usize]
                    } else {
                        0
                    });
            }
            self.inner_mut().core.l_qseq = seq.len() as i32;
            i += seq.len().div_ceil(2);
        }

        // qual
        utils::copy_memory(qual, &mut data[i..]);
    }

    /// Replace the sequence with a new one of the same length.
    pub fn set_seq(&mut self, new_seq: &[u8]) {
        assert_eq!(
            new_seq.len(),
            self.seq_len(),
            "new_seq.len() must equal current seq.len()"
        );

        let seq_offset = self.qname_capacity() + self.cigar_len() * 4;
        // SAFETY: inner.data is valid for l_data bytes; seq_offset is within bounds (same length seq).
        let data =
            unsafe { slice::from_raw_parts_mut(self.inner.data, self.inner().l_data as usize) };
        for j in (0..new_seq.len()).step_by(2) {
            data[seq_offset + j / 2] = (ENCODE_BASE[new_seq[j] as usize] << 4)
                | (if j + 1 < new_seq.len() {
                    ENCODE_BASE[new_seq[j + 1] as usize]
                } else {
                    0
                });
        }
    }

    /// Replace current qname with a new one.
    pub fn set_qname(&mut self, new_qname: &[u8]) {
        // 251 + 1NUL is the max 32-bit aligned value that fits in u8
        assert!(new_qname.len() < 252);

        let old_q_len = self.qname_capacity();
        // We're going to add a terminal NUL
        let extranul = extranul_from_qname(new_qname);
        let new_q_len = new_qname.len() + 1 + extranul;

        // Length of data after qname
        let other_len = self.inner_mut().l_data - old_q_len as i32;

        if new_q_len < old_q_len && self.inner().l_data > (old_q_len as i32) {
            self.inner_mut().l_data -= (old_q_len - new_q_len) as i32;
        } else if new_q_len > old_q_len {
            self.inner_mut().l_data += (new_q_len - old_q_len) as i32;

            // Reallocate if necessary
            if (self.inner().m_data as i32) < self.inner().l_data {
                // Verbosity due to lexical borrowing
                let l_data = self.inner().l_data;
                self.realloc_var_data(l_data as usize);
            }
        }

        if new_q_len != old_q_len {
            // Move other data to new location
            // SAFETY: inner.data is valid for l_data bytes (realloc'd if needed above); memmove handles overlap.
            unsafe {
                let data = slice::from_raw_parts_mut(self.inner.data, self.inner().l_data as usize);

                ::libc::memmove(
                    data.as_mut_ptr().add(new_q_len) as *mut ::libc::c_void,
                    data.as_mut_ptr().add(old_q_len) as *mut ::libc::c_void,
                    other_len as usize,
                );
            }
        }

        // Copy qname data
        // SAFETY: inner.data is valid for l_data bytes.
        let data =
            unsafe { slice::from_raw_parts_mut(self.inner.data, self.inner().l_data as usize) };
        utils::copy_memory(new_qname, data);
        for i in 0..=extranul {
            data[new_q_len - i - 1] = b'\0';
        }
        self.inner_mut().core.l_qname = new_q_len as u16;
        self.inner_mut().core.l_extranul = extranul as u8;
    }

    /// Replace current cigar with a new one.
    pub fn set_cigar(&mut self, new_cigar: Option<&CigarString>) {
        self.cigar = None;

        let qname_data_len = self.qname_capacity();
        let old_cigar_data_len = self.cigar_len() * 4;

        // Length of data after cigar
        let other_data_len = self.inner_mut().l_data - (qname_data_len + old_cigar_data_len) as i32;

        let new_cigar_len = match new_cigar {
            Some(x) => x.len(),
            None => 0,
        };
        let new_cigar_data_len = new_cigar_len * 4;

        if new_cigar_data_len < old_cigar_data_len {
            self.inner_mut().l_data -= (old_cigar_data_len - new_cigar_data_len) as i32;
        } else if new_cigar_data_len > old_cigar_data_len {
            self.inner_mut().l_data += (new_cigar_data_len - old_cigar_data_len) as i32;

            // Reallocate if necessary
            if (self.inner().m_data as i32) < self.inner().l_data {
                // Verbosity due to lexical borrowing
                let l_data = self.inner().l_data;
                self.realloc_var_data(l_data as usize);
            }
        }

        if new_cigar_data_len != old_cigar_data_len {
            // Move other data to new location
            // SAFETY: inner.data is valid for l_data bytes (realloc'd if needed above); memmove handles overlap.
            unsafe {
                ::libc::memmove(
                    self.inner.data.add(qname_data_len + new_cigar_data_len) as *mut ::libc::c_void,
                    self.inner.data.add(qname_data_len + old_cigar_data_len) as *mut ::libc::c_void,
                    other_data_len as usize,
                );
            }
        }

        // Copy cigar data
        if let Some(cigar_string) = new_cigar {
            // SAFETY: cigar offset is 4-byte aligned (qname padded); inner.data valid for l_data bytes.
            let cigar_data = unsafe {
                #[allow(clippy::cast_ptr_alignment)]
                slice::from_raw_parts_mut(
                    self.inner.data.add(qname_data_len) as *mut u32,
                    cigar_string.len(),
                )
            };
            for (i, c) in cigar_string.iter().enumerate() {
                cigar_data[i] = c.encode();
            }
        }
        self.inner_mut().core.n_cigar = new_cigar_len as u32;
    }

    fn realloc_var_data(&mut self, new_len: usize) {
        // pad request
        let new_len = new_len as u32;
        let new_request = new_len + 32 - (new_len % 32);

        // SAFETY: inner.data is either null (initial state) or was previously allocated by malloc/realloc.
        let ptr = unsafe {
            ::libc::realloc(
                self.inner().data as *mut ::libc::c_void,
                new_request as usize,
            ) as *mut u8
        };

        if ptr.is_null() {
            panic!("ran out of memory in rust_htslib trying to realloc");
        }

        // don't update m_data until we know we have
        // a successful allocation.
        self.inner_mut().m_data = new_request;
        self.inner_mut().data = ptr;

        // we now own inner.data
        self.own = true;
    }

    pub fn cigar_len(&self) -> usize {
        self.inner().core.n_cigar as usize
    }

    /// Get reference to raw cigar string representation (as stored in BAM file).
    /// Usually, the method `Record::cigar` should be used instead.
    pub fn raw_cigar(&self) -> &[u32] {
        // SAFETY: cigar data starts at a 4-byte-aligned offset (qname is padded); length from n_cigar.
        #[allow(clippy::cast_ptr_alignment)]
        unsafe {
            slice::from_raw_parts(
                self.data()[self.qname_capacity()..].as_ptr() as *const u32,
                self.cigar_len(),
            )
        }
    }

    /// Return unpacked cigar string. If the cigar has been cached via
    /// [`cache_cigar`](Self::cache_cigar), this is a cheap clone (Arc refcount bump).
    /// Otherwise, the cigar is decoded from the raw BAM data.
    pub fn cigar(&self) -> CigarStringView {
        match self.cigar {
            Some(ref c) => c.clone(),
            None => self.unpack_cigar(),
        }
    }

    // Return unpacked cigar string. This returns None unless you have first called `bam::Record::cache_cigar`.
    pub fn cigar_cached(&self) -> Option<&CigarStringView> {
        self.cigar.as_ref()
    }

    /// Decode the cigar string and cache it inside the `Record`
    pub fn cache_cigar(&mut self) {
        self.cigar = Some(self.unpack_cigar())
    }

    /// Unpack cigar string. Complexity: O(k) with k being the length of the cigar string.
    fn unpack_cigar(&self) -> CigarStringView {
        CigarString(
            self.raw_cigar()
                .iter()
                .map(|&c| Cigar::from_raw(c))
                .collect(),
        )
        .into_view(self.pos())
    }

    pub fn seq_len(&self) -> usize {
        self.inner().core.l_qseq as usize
    }

    fn seq_data(&self) -> &[u8] {
        let offset = self.qname_capacity() + self.cigar_len() * 4;
        &self.data()[offset..][..self.seq_len().div_ceil(2)]
    }

    /// Get read sequence. Complexity: O(1).
    pub fn seq(&self) -> Seq<'_> {
        Seq {
            encoded: self.seq_data(),
            len: self.seq_len(),
        }
    }

    /// Get base qualities (PHRED-scaled probability that base is wrong).
    /// This does not entail any offsets, hence the qualities can be used directly without
    /// e.g. subtracting 33. Complexity: O(1).
    pub fn qual(&self) -> &[u8] {
        &self.data()[self.qname_capacity() + self.cigar_len() * 4 + self.seq_len().div_ceil(2)..]
            [..self.seq_len()]
    }

    /// Get the raw auxiliary data as a byte slice.
    fn raw_aux_data(&self) -> &[u8] {
        &self.data()[self.qname_capacity()
            + self.cigar_len() * std::mem::size_of::<u32>()
            + self.seq_len().div_ceil(2)
            + self.seq_len()..]
    }

    /// Look up an auxiliary field by its tag.
    ///
    /// Only the first two bytes of a given tag are used for the look-up of a field.
    /// See [`Aux`] for more details.
    pub fn aux(&self, tag: &[u8]) -> Result<Aux<'_>> {
        if tag.len() < 2 {
            return Err(Error::AuxStringError);
        }
        let raw = self.raw_aux_data();
        match aux_tag_search(raw, tag) {
            Some(offset) => {
                // SAFETY: offset is within raw (returned by aux_tag_search); raw is valid for record's lifetime.
                unsafe { parse_aux_field(raw.as_ptr().add(offset)).map(|(v, _)| v) }
            }
            None => Err(Error::AuxTagNotFound),
        }
    }

    /// Returns an iterator over the auxiliary fields of the record.
    ///
    /// When an error occurs, the `Err` variant will be returned
    /// and the iterator will not be able to advance anymore.
    pub fn aux_iter(&self) -> AuxIter<'_> {
        AuxIter {
            // In order to get to the aux data section of a `bam::Record`
            // we need to skip fields in front of it
            aux: &self.data()[
                // NUL terminated read name:
                self.qname_capacity()
                // CIGAR (uint32_t):
                + self.cigar_len() * std::mem::size_of::<u32>()
                // Read sequence (4-bit encoded):
                + self.seq_len().div_ceil(2)
                // Base qualities (char):
                + self.seq_len()..],
        }
    }

    /// Add auxiliary data.
    pub fn push_aux(&mut self, tag: &[u8], value: Aux<'_>) -> Result<()> {
        // Don't allow pushing aux data when the given tag is already present in the record.
        // `htslib` seems to allow this (for non-array values), which can lead to problems
        // since retrieving aux fields consumes &[u8; 2] and yields one field only.
        if self.aux(tag).is_ok() {
            return Err(Error::AuxTagAlreadyPresent);
        }
        self.push_aux_unchecked(tag, value)
    }

    /// Add auxiliary data, without checking if the tag is present.
    ///
    /// The caller should ensure that the same tag is not pushed more than once.
    /// This is provided as a performance optimization.
    pub fn push_aux_unchecked(&mut self, tag: &[u8], value: Aux<'_>) -> Result<()> {
        let ctag = tag.as_ptr() as *mut c_char;
        // SAFETY: self.inner_ptr_mut() is valid; ctag points to at least 2 bytes; value data is valid for the call.
        let ret = unsafe {
            match value {
                Aux::Char(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'A' as c_char,
                    size_of::<u8>() as i32,
                    [v].as_mut_ptr(),
                ),
                Aux::I8(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'c' as c_char,
                    size_of::<i8>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::U8(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'C' as c_char,
                    size_of::<u8>() as i32,
                    [v].as_mut_ptr(),
                ),
                Aux::I16(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b's' as c_char,
                    size_of::<i16>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::U16(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'S' as c_char,
                    size_of::<u16>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::I32(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'i' as c_char,
                    size_of::<i32>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::U32(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'I' as c_char,
                    size_of::<u32>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::Float(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'f' as c_char,
                    size_of::<f32>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                // Not part of specs but implemented in `htslib`:
                Aux::Double(v) => htslib::bam_aux_append(
                    self.inner_ptr_mut(),
                    ctag,
                    b'd' as c_char,
                    size_of::<f64>() as i32,
                    [v].as_mut_ptr() as *mut u8,
                ),
                Aux::String(v) => {
                    let c_str = ffi::CString::new(v).map_err(|_| Error::AuxStringError)?;
                    htslib::bam_aux_append(
                        self.inner_ptr_mut(),
                        ctag,
                        b'Z' as c_char,
                        (v.len() + 1) as i32,
                        c_str.as_ptr() as *mut u8,
                    )
                }
                Aux::HexByteArray(v) => {
                    let c_str = ffi::CString::new(v).map_err(|_| Error::AuxStringError)?;
                    htslib::bam_aux_append(
                        self.inner_ptr_mut(),
                        ctag,
                        b'H' as c_char,
                        (v.len() + 1) as i32,
                        c_str.as_ptr() as *mut u8,
                    )
                }
                // Not sure it's safe to cast an immutable slice to a mutable pointer in the following branches
                Aux::ArrayI8(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'c',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'c',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayU8(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'C',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'C',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayI16(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b's',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b's',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayU16(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'S',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'S',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayI32(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'i',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'i',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayU32(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'I',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'I',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayFloat(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'f',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'f',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
            }
        };

        if ret < 0 {
            Err(Error::Aux)
        } else {
            Ok(())
        }
    }

    /// Update or add auxiliary data.
    pub fn update_aux(&mut self, tag: &[u8], value: Aux<'_>) -> Result<()> {
        // Update existing aux data for the given tag if already present in the record
        // without changing the ordering of tags in the record or append aux data at
        // the end of the existing aux records if it is a new tag.

        let ctag = tag.as_ptr() as *mut c_char;
        // SAFETY: self.inner_ptr_mut() is valid; ctag points to at least 2 bytes; value data is valid for the call.
        let ret = unsafe {
            match value {
                Aux::Char(_v) => return Err(Error::AuxTagUpdatingNotSupported),
                Aux::I8(v) => htslib::bam_aux_update_int(self.inner_ptr_mut(), ctag, v as i64),
                Aux::U8(v) => htslib::bam_aux_update_int(self.inner_ptr_mut(), ctag, v as i64),
                Aux::I16(v) => htslib::bam_aux_update_int(self.inner_ptr_mut(), ctag, v as i64),
                Aux::U16(v) => htslib::bam_aux_update_int(self.inner_ptr_mut(), ctag, v as i64),
                Aux::I32(v) => htslib::bam_aux_update_int(self.inner_ptr_mut(), ctag, v as i64),
                Aux::U32(v) => htslib::bam_aux_update_int(self.inner_ptr_mut(), ctag, v as i64),
                Aux::Float(v) => htslib::bam_aux_update_float(self.inner_ptr_mut(), ctag, v),
                // Not part of specs but implemented in `htslib`:
                Aux::Double(v) => {
                    htslib::bam_aux_update_float(self.inner_ptr_mut(), ctag, v as f32)
                }
                Aux::String(v) => {
                    let c_str = ffi::CString::new(v).map_err(|_| Error::AuxStringError)?;
                    htslib::bam_aux_update_str(
                        self.inner_ptr_mut(),
                        ctag,
                        (v.len() + 1) as i32,
                        c_str.as_ptr() as *const c_char,
                    )
                }
                Aux::HexByteArray(_v) => return Err(Error::AuxTagUpdatingNotSupported),
                // Not sure it's safe to cast an immutable slice to a mutable pointer in the following branches
                Aux::ArrayI8(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'c',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'c',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayU8(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'C',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'C',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayI16(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b's',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b's',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayU16(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'S',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'S',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayI32(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'i',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'i',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayU32(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'I',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'I',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
                Aux::ArrayFloat(aux_array) => match aux_array {
                    AuxArray::TargetType(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'f',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                    AuxArray::RawLeBytes(inner) => htslib::bam_aux_update_array(
                        self.inner_ptr_mut(),
                        ctag,
                        b'f',
                        inner.len() as u32,
                        inner.slice.as_ptr() as *mut ::libc::c_void,
                    ),
                },
            }
        };

        if ret < 0 {
            Err(Error::Aux)
        } else {
            Ok(())
        }
    }

    // Delete auxiliary tag.
    pub fn remove_aux(&mut self, tag: &[u8]) -> Result<()> {
        if tag.len() < 2 {
            return Err(Error::AuxStringError);
        }
        let raw = self.raw_aux_data();
        match aux_tag_search(raw, tag) {
            Some(offset) => {
                // SAFETY: offset is within raw_aux_data; we own the record mutably;
                // bam_aux_del modifies data in place.
                let ptr = unsafe { raw.as_ptr().add(offset) as *mut u8 };
                unsafe { htslib::bam_aux_del(self.inner_ptr_mut(), ptr) };
                Ok(())
            }
            None => Err(Error::AuxTagNotFound),
        }
    }

    /// Access the base modifications associated with this Record through the MM tag.
    /// Example:
    /// ```
    ///    use rust_htslib::bam::{Read, Reader, Record};
    ///    let mut bam = Reader::from_path("test/base_mods/MM-orient.sam").unwrap();
    ///    let mut mod_count = 0;
    ///    for r in bam.records() {
    ///        let record = r.unwrap();
    ///        if let Ok(mods) = record.basemods_iter() {
    ///            // print metadata for the modifications present in this record
    ///            for mod_code in mods.recorded() {
    ///                if let Ok(mod_metadata) = mods.query_type(*mod_code) {
    ///                    println!("mod found with code {}/{} flags: [{} {} {}]",
    ///                              mod_code, *mod_code as u8 as char,
    ///                              mod_metadata.strand, mod_metadata.implicit, mod_metadata.canonical as u8 as char);
    ///                }
    ///            }
    ///
    ///            // iterate over the modifications in this record
    ///            // the modifications are returned as a tuple with the
    ///            // position within SEQ and an hts_base_mod struct
    ///            for res in mods {
    ///                if let Ok( (position, m) ) = res {
    ///                    println!("{} {},{}", position, m.modified_base as u8 as char, m.qual);
    ///                    mod_count += 1;
    ///                }
    ///            }
    ///        };
    ///    }
    ///    assert_eq!(mod_count, 14);
    /// ```
    pub fn basemods_iter(&self) -> Result<BaseModificationsIter<'_>> {
        BaseModificationsIter::new(self)
    }

    /// An iterator that returns all of the modifications for each position as a vector.
    /// This is useful for the case where multiple possible modifications can be annotated
    /// at a single position (for example a C could be 5-mC or 5-hmC)
    pub fn basemods_position_iter(&self) -> Result<BaseModificationsPositionIter<'_>> {
        BaseModificationsPositionIter::new(self)
    }

    /// Infer read pair orientation from record. Returns `SequenceReadPairOrientation::None` if record
    /// is not paired, mates are not mapping to the same contig, or mates start at the
    /// same position.
    pub fn read_pair_orientation(&self) -> SequenceReadPairOrientation {
        if self.is_paired()
            && !self.is_unmapped()
            && !self.is_mate_unmapped()
            && self.tid() == self.mtid()
        {
            if self.pos() == self.mpos() {
                // both reads start at the same position, we cannot decide on the orientation.
                return SequenceReadPairOrientation::None;
            }

            let (pos_1, pos_2, fwd_1, fwd_2) = if self.is_first_in_template() {
                (
                    self.pos(),
                    self.mpos(),
                    !self.is_reverse(),
                    !self.is_mate_reverse(),
                )
            } else {
                (
                    self.mpos(),
                    self.pos(),
                    !self.is_mate_reverse(),
                    !self.is_reverse(),
                )
            };

            if pos_1 < pos_2 {
                match (fwd_1, fwd_2) {
                    (true, true) => SequenceReadPairOrientation::F1F2,
                    (true, false) => SequenceReadPairOrientation::F1R2,
                    (false, true) => SequenceReadPairOrientation::R1F2,
                    (false, false) => SequenceReadPairOrientation::R1R2,
                }
            } else {
                match (fwd_2, fwd_1) {
                    (true, true) => SequenceReadPairOrientation::F2F1,
                    (true, false) => SequenceReadPairOrientation::F2R1,
                    (false, true) => SequenceReadPairOrientation::R2F1,
                    (false, false) => SequenceReadPairOrientation::R2R1,
                }
            }
        } else {
            SequenceReadPairOrientation::None
        }
    }

    flag!(is_paired, set_paired, unset_paired, 1u16);
    flag!(is_proper_pair, set_proper_pair, unset_proper_pair, 2u16);
    flag!(is_unmapped, set_unmapped, unset_unmapped, 4u16);
    flag!(
        is_mate_unmapped,
        set_mate_unmapped,
        unset_mate_unmapped,
        8u16
    );
    flag!(is_reverse, set_reverse, unset_reverse, 16u16);
    flag!(is_mate_reverse, set_mate_reverse, unset_mate_reverse, 32u16);
    flag!(
        is_first_in_template,
        set_first_in_template,
        unset_first_in_template,
        64u16
    );
    flag!(
        is_last_in_template,
        set_last_in_template,
        unset_last_in_template,
        128u16
    );
    flag!(is_secondary, set_secondary, unset_secondary, 256u16);
    flag!(
        is_quality_check_failed,
        set_quality_check_failed,
        unset_quality_check_failed,
        512u16
    );
    flag!(is_duplicate, set_duplicate, unset_duplicate, 1024u16);
    flag!(
        is_supplementary,
        set_supplementary,
        unset_supplementary,
        2048u16
    );
}

impl Drop for Record {
    fn drop(&mut self) {
        if self.own {
            // SAFETY: inner.data was allocated by malloc/realloc (tracked by self.own); free is symmetric.
            unsafe { ::libc::free(self.inner.data as *mut ::libc::c_void) }
        }
    }
}

impl SequenceRead for Record {
    fn name(&self) -> &[u8] {
        self.qname()
    }

    fn base(&self, i: usize) -> u8 {
        *decode_base_unchecked(encoded_base(self.seq_data(), i))
    }

    fn base_qual(&self, i: usize) -> u8 {
        self.qual()[i]
    }

    fn len(&self) -> usize {
        self.seq_len()
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl genome::AbstractInterval for Record {
    /// Return contig name. Panics if record does not know its header (which happens if it has not been read from a file).
    fn contig(&self) -> &str {
        let tid = self.tid();
        if tid < 0 {
            panic!("invalid tid, must be at least zero");
        }
        str::from_utf8(
            self.header
                .as_ref()
                .expect(
                    "header must be set (this is the case if the record has been read from a file)",
                )
                .tid2name(tid as u32)
                .expect("tid out of bounds"),
        )
        .expect("unable to interpret contig name as UTF-8")
    }

    /// Return genomic range covered by alignment. Panics if `Record::cache_cigar()` has not been called first or `Record::pos()` is less than zero.
    fn range(&self) -> ops::Range<genome::Position> {
        let end_pos = self
            .cigar_cached()
            .expect("cigar has not been cached yet, call cache_cigar() first")
            .end_pos() as u64;

        if self.pos() < 0 {
            panic!("invalid position, must be positive")
        }

        self.pos() as u64..end_pos
    }
}

/// Auxiliary record data
///
/// The specification allows a wide range of types to be stored as an auxiliary data field of a BAM record.
///
/// Please note that the [`Aux::Double`] variant is _not_ part of the specification, but it is supported by `htslib`.
///
/// # Examples
///
/// ```
/// use rust_htslib::{
///     bam,
///     bam::record::{Aux, AuxArray},
///     errors::Error,
/// };
///
/// //Set up BAM record
/// let bam_header = bam::Header::new();
/// let mut record = bam::Record::from_sam(
///     &mut bam::HeaderView::from_header(&bam_header),
///     "ali1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF".as_bytes(),
/// )
/// .unwrap();
///
/// // Add an integer field
/// let aux_integer_field = Aux::I32(1234);
/// record.push_aux(b"XI", aux_integer_field).unwrap();
///
/// match record.aux(b"XI") {
///     Ok(value) => {
///         // Typically, callers expect an aux field to be of a certain type.
///         // If that's not the case, the value can be `match`ed exhaustively.
///         if let Aux::I32(v) = value {
///             assert_eq!(v, 1234);
///         }
///     }
///     Err(e) => {
///         panic!("Error reading aux field: {}", e);
///     }
/// }
///
/// // Add an array field
/// let array_like_data = vec![0.4, 0.3, 0.2, 0.1];
/// let slice_of_data = &array_like_data;
/// let aux_array: AuxArray<f32> = slice_of_data.into();
/// let aux_array_field = Aux::ArrayFloat(aux_array);
/// record.push_aux(b"XA", aux_array_field).unwrap();
///
/// if let Ok(Aux::ArrayFloat(array)) = record.aux(b"XA") {
///     let read_array = array.iter().collect::<Vec<_>>();
///     assert_eq!(read_array, array_like_data);
/// } else {
///     panic!("Could not read array data");
/// }
/// ```
#[derive(Debug, PartialEq)]
pub enum Aux<'a> {
    Char(u8),
    I8(i8),
    U8(u8),
    I16(i16),
    U16(u16),
    I32(i32),
    U32(u32),
    Float(f32),
    Double(f64), // Not part of specs but implemented in `htslib`
    String(&'a str),
    HexByteArray(&'a str),
    ArrayI8(AuxArray<'a, i8>),
    ArrayU8(AuxArray<'a, u8>),
    ArrayI16(AuxArray<'a, i16>),
    ArrayU16(AuxArray<'a, u16>),
    ArrayI32(AuxArray<'a, i32>),
    ArrayU32(AuxArray<'a, u32>),
    ArrayFloat(AuxArray<'a, f32>),
}

// SAFETY: Aux borrows only immutable data (&str, &[u8]) from the record's data buffer.
unsafe impl Send for Aux<'_> {}
unsafe impl Sync for Aux<'_> {}

/// Types that can be used in aux arrays.
pub trait AuxArrayElement: Copy {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self>;
}

impl AuxArrayElement for i8 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_i8().ok()
    }
}
impl AuxArrayElement for u8 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_u8().ok()
    }
}
impl AuxArrayElement for i16 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_i16::<LittleEndian>().ok()
    }
}
impl AuxArrayElement for u16 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_u16::<LittleEndian>().ok()
    }
}
impl AuxArrayElement for i32 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_i32::<LittleEndian>().ok()
    }
}
impl AuxArrayElement for u32 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_u32::<LittleEndian>().ok()
    }
}
impl AuxArrayElement for f32 {
    fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        std::io::Cursor::new(bytes).read_f32::<LittleEndian>().ok()
    }
}

/// Provides access to aux arrays.
///
/// Provides methods to either retrieve single elements or an iterator over the
/// array.
///
/// This type is used for wrapping both, array data that was read from a
/// BAM record and slices of data that are going to be stored in one.
///
/// In order to be able to add an `AuxArray` field to a BAM record, `AuxArray`s
/// can be constructed via the `From` trait which is implemented for all
/// supported types (see [`AuxArrayElement`] for a list).
///
/// # Examples
///
/// ```
/// use rust_htslib::{
///     bam,
///     bam::record::{Aux, AuxArray},
/// };
///
/// //Set up BAM record
/// let bam_header = bam::Header::new();
/// let mut record = bam::Record::from_sam(
///     &mut bam::HeaderView::from_header(&bam_header),
///     "ali1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF".as_bytes(),
/// ).unwrap();
///
/// let data = vec![0.4, 0.3, 0.2, 0.1];
/// let slice_of_data = &data;
/// let aux_array: AuxArray<f32> = slice_of_data.into();
/// let aux_field = Aux::ArrayFloat(aux_array);
/// record.push_aux(b"XA", aux_field);
///
/// if let Ok(Aux::ArrayFloat(array)) = record.aux(b"XA") {
///     // Retrieve the second element from the array
///     assert_eq!(array.get(1).unwrap(), 0.3);
///     // Iterate over the array and collect it into a `Vec`
///     let read_array = array.iter().collect::<Vec<_>>();
///     assert_eq!(read_array, data);
/// } else {
///     panic!("Could not read array data");
/// }
/// ```
#[derive(Debug)]
pub enum AuxArray<'a, T> {
    TargetType(AuxArrayTargetType<'a, T>),
    RawLeBytes(AuxArrayRawLeBytes<'a, T>),
}

impl<T> PartialEq<AuxArray<'_, T>> for AuxArray<'_, T>
where
    T: AuxArrayElement + PartialEq,
{
    fn eq(&self, other: &AuxArray<'_, T>) -> bool {
        use AuxArray::*;
        match (self, other) {
            (TargetType(v), TargetType(v_other)) => v == v_other,
            (RawLeBytes(v), RawLeBytes(v_other)) => v == v_other,
            (TargetType(_), RawLeBytes(_)) => self.iter().eq(other.iter()),
            (RawLeBytes(_), TargetType(_)) => self.iter().eq(other.iter()),
        }
    }
}

/// Create AuxArrays from slices of allowed target types.
impl<'a, I, T> From<&'a T> for AuxArray<'a, I>
where
    I: AuxArrayElement,
    T: AsRef<[I]> + ?Sized,
{
    fn from(src: &'a T) -> Self {
        AuxArray::TargetType(AuxArrayTargetType {
            slice: src.as_ref(),
        })
    }
}

impl<'a, T> AuxArray<'a, T>
where
    T: AuxArrayElement,
{
    /// Returns the element at a position or None if out of bounds.
    pub fn get(&self, index: usize) -> Option<T> {
        match self {
            AuxArray::TargetType(v) => v.get(index),
            AuxArray::RawLeBytes(v) => v.get(index),
        }
    }

    /// Returns the number of elements in the array.
    pub fn len(&self) -> usize {
        match self {
            AuxArray::TargetType(a) => a.len(),
            AuxArray::RawLeBytes(a) => a.len(),
        }
    }

    /// Returns true if the array contains no elements.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns an iterator over the array.
    pub fn iter(&self) -> AuxArrayIter<'_, T> {
        AuxArrayIter {
            index: 0,
            array: self,
        }
    }

    /// Create AuxArrays from raw byte slices borrowed from `bam::Record`.
    fn from_bytes(bytes: &'a [u8]) -> Self {
        Self::RawLeBytes(AuxArrayRawLeBytes {
            slice: bytes,
            phantom_data: PhantomData,
        })
    }
}

/// Encapsulates slice of target type.
#[doc(hidden)]
#[derive(Debug, PartialEq)]
pub struct AuxArrayTargetType<'a, T> {
    slice: &'a [T],
}

impl<T> AuxArrayTargetType<'_, T>
where
    T: AuxArrayElement,
{
    fn get(&self, index: usize) -> Option<T> {
        self.slice.get(index).copied()
    }

    fn len(&self) -> usize {
        self.slice.len()
    }
}

/// Encapsulates slice of raw bytes to prevent it from being accidentally accessed.
#[doc(hidden)]
#[derive(Debug, PartialEq)]
pub struct AuxArrayRawLeBytes<'a, T> {
    slice: &'a [u8],
    phantom_data: PhantomData<T>,
}

impl<T> AuxArrayRawLeBytes<'_, T>
where
    T: AuxArrayElement,
{
    fn get(&self, index: usize) -> Option<T> {
        let type_size = std::mem::size_of::<T>();
        if index * type_size + type_size > self.slice.len() {
            return None;
        }
        T::from_le_bytes(&self.slice[index * type_size..][..type_size])
    }

    fn len(&self) -> usize {
        self.slice.len() / std::mem::size_of::<T>()
    }
}

/// Aux array iterator
///
/// This struct is created by the [`AuxArray::iter`] method.
pub struct AuxArrayIter<'a, T> {
    index: usize,
    array: &'a AuxArray<'a, T>,
}

impl<T> Iterator for AuxArrayIter<'_, T>
where
    T: AuxArrayElement,
{
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        let value = self.array.get(self.index);
        self.index += 1;
        value
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let array_length = self.array.len() - self.index;
        (array_length, Some(array_length))
    }
}

/// Pure Rust replacement for `htslib::bam_aux_get`.
///
/// Performs a linear scan through the packed BAM auxiliary data looking for a
/// 2-byte tag. Returns the byte offset of the type byte within `aux_data` for
/// the matching field, or `None` if the tag is not found or the data is
/// malformed.
///
/// The aux data format is a packed sequence of fields:
///   `[tag0, tag1, type, payload...]`
/// where payload length depends on the type byte.
fn aux_tag_search(aux_data: &[u8], tag: &[u8]) -> Option<usize> {
    debug_assert!(tag.len() >= 2);
    let tag0 = tag[0];
    let tag1 = tag[1];
    let mut pos = 0;
    let len = aux_data.len();

    while pos + 3 <= len {
        if aux_data[pos] == tag0 && aux_data[pos + 1] == tag1 {
            return Some(pos + 2);
        }

        let type_byte = aux_data[pos + 2];
        pos += 3;
        match type_byte {
            b'A' | b'c' | b'C' => pos += 1,
            b's' | b'S' => pos += 2,
            b'i' | b'I' | b'f' => pos += 4,
            b'd' => pos += 8,
            b'Z' | b'H' => match memchr::memchr(0, &aux_data[pos..]) {
                Some(nul_offset) => pos += nul_offset + 1,
                None => return None,
            },
            b'B' => {
                if pos + 5 > len {
                    return None;
                }
                let sub_type = aux_data[pos];
                let count = u32::from_le_bytes([
                    aux_data[pos + 1],
                    aux_data[pos + 2],
                    aux_data[pos + 3],
                    aux_data[pos + 4],
                ]) as usize;
                let elem_size = match sub_type {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return None,
                };
                pos += 5 + count * elem_size;
            }
            _ => return None,
        }
    }
    None
}

/// Parse a single aux field from a raw pointer to the type byte.
///
/// `aux` must point to the type-identifier byte of the field (the byte
/// returned by `bam_aux_get`, or the third byte of an on-disk aux entry).
/// Returns the parsed value and the total number of bytes consumed
/// (tag bytes + type byte + payload), so that callers can advance a
/// cursor through a packed aux buffer.
///
/// # Safety
/// `aux` must be non-null and point into a valid, correctly formatted BAM
/// aux buffer that lives at least as long as `'a`.
unsafe fn parse_aux_field<'a>(aux: *const u8) -> Result<(Aux<'a>, usize)> {
    const TAG_LEN: isize = 2;
    // Used for skipping type identifier
    const TYPE_ID_LEN: isize = 1;

    if aux.is_null() {
        return Err(Error::AuxTagNotFound);
    }

    let (data, type_size) = match *aux {
        b'A' => {
            let type_size = size_of::<u8>();
            (Aux::Char(*aux.offset(TYPE_ID_LEN)), type_size)
        }
        b'c' => {
            let type_size = size_of::<i8>();
            (Aux::I8(*aux.offset(TYPE_ID_LEN).cast::<i8>()), type_size)
        }
        b'C' => {
            let type_size = size_of::<u8>();
            (Aux::U8(*aux.offset(TYPE_ID_LEN)), type_size)
        }
        b's' => {
            let type_size = size_of::<i16>();
            (
                Aux::I16(
                    slice::from_raw_parts(aux.offset(TYPE_ID_LEN), type_size)
                        .read_i16::<LittleEndian>()
                        .map_err(|_| Error::AuxParsingError)?,
                ),
                type_size,
            )
        }
        b'S' => {
            let type_size = size_of::<u16>();
            (
                Aux::U16(
                    slice::from_raw_parts(aux.offset(TYPE_ID_LEN), type_size)
                        .read_u16::<LittleEndian>()
                        .map_err(|_| Error::AuxParsingError)?,
                ),
                type_size,
            )
        }
        b'i' => {
            let type_size = size_of::<i32>();
            (
                Aux::I32(
                    slice::from_raw_parts(aux.offset(TYPE_ID_LEN), type_size)
                        .read_i32::<LittleEndian>()
                        .map_err(|_| Error::AuxParsingError)?,
                ),
                type_size,
            )
        }
        b'I' => {
            let type_size = size_of::<u32>();
            (
                Aux::U32(
                    slice::from_raw_parts(aux.offset(TYPE_ID_LEN), type_size)
                        .read_u32::<LittleEndian>()
                        .map_err(|_| Error::AuxParsingError)?,
                ),
                type_size,
            )
        }
        b'f' => {
            let type_size = size_of::<f32>();
            (
                Aux::Float(
                    slice::from_raw_parts(aux.offset(TYPE_ID_LEN), type_size)
                        .read_f32::<LittleEndian>()
                        .map_err(|_| Error::AuxParsingError)?,
                ),
                type_size,
            )
        }
        b'd' => {
            let type_size = size_of::<f64>();
            (
                Aux::Double(
                    slice::from_raw_parts(aux.offset(TYPE_ID_LEN), type_size)
                        .read_f64::<LittleEndian>()
                        .map_err(|_| Error::AuxParsingError)?,
                ),
                type_size,
            )
        }
        b'Z' | b'H' => {
            let c_str = ffi::CStr::from_ptr(aux.offset(TYPE_ID_LEN).cast::<c_char>());
            let rust_str = c_str.to_str().map_err(|_| Error::AuxParsingError)?;
            (Aux::String(rust_str), c_str.to_bytes_with_nul().len())
        }
        b'B' => {
            const ARRAY_INNER_TYPE_LEN: isize = 1;
            const ARRAY_COUNT_LEN: isize = 4;

            // Used for skipping metadata
            let array_data_offset = TYPE_ID_LEN + ARRAY_INNER_TYPE_LEN + ARRAY_COUNT_LEN;

            let length = slice::from_raw_parts(aux.offset(TYPE_ID_LEN + ARRAY_INNER_TYPE_LEN), 4)
                .read_u32::<LittleEndian>()
                .map_err(|_| Error::AuxParsingError)? as usize;

            // Return tuples of an `Aux` enum and the length of data + metadata in bytes
            let (array_data, array_size) = match *aux.offset(TYPE_ID_LEN) {
                b'c' => (
                    Aux::ArrayI8(AuxArray::<'a, i8>::from_bytes(slice::from_raw_parts(
                        aux.offset(array_data_offset),
                        length,
                    ))),
                    length,
                ),
                b'C' => (
                    Aux::ArrayU8(AuxArray::<'a, u8>::from_bytes(slice::from_raw_parts(
                        aux.offset(array_data_offset),
                        length,
                    ))),
                    length,
                ),
                b's' => (
                    Aux::ArrayI16(AuxArray::<'a, i16>::from_bytes(slice::from_raw_parts(
                        aux.offset(array_data_offset),
                        length * size_of::<i16>(),
                    ))),
                    length * std::mem::size_of::<i16>(),
                ),
                b'S' => (
                    Aux::ArrayU16(AuxArray::<'a, u16>::from_bytes(slice::from_raw_parts(
                        aux.offset(array_data_offset),
                        length * size_of::<u16>(),
                    ))),
                    length * std::mem::size_of::<u16>(),
                ),
                b'i' => (
                    Aux::ArrayI32(AuxArray::<'a, i32>::from_bytes(slice::from_raw_parts(
                        aux.offset(array_data_offset),
                        length * size_of::<i32>(),
                    ))),
                    length * std::mem::size_of::<i32>(),
                ),
                b'I' => (
                    Aux::ArrayU32(AuxArray::<'a, u32>::from_bytes(slice::from_raw_parts(
                        aux.offset(array_data_offset),
                        length * size_of::<u32>(),
                    ))),
                    length * std::mem::size_of::<u32>(),
                ),
                b'f' => (
                    Aux::ArrayFloat(AuxArray::<f32>::from_bytes(slice::from_raw_parts(
                        aux.offset(array_data_offset),
                        length * size_of::<f32>(),
                    ))),
                    length * std::mem::size_of::<f32>(),
                ),
                _ => {
                    return Err(Error::AuxUnknownType);
                }
            };
            (
                array_data,
                // Offset: array-specific metadata + array size
                ARRAY_INNER_TYPE_LEN as usize + ARRAY_COUNT_LEN as usize + array_size,
            )
        }
        _ => {
            return Err(Error::AuxUnknownType);
        }
    };

    // Offset: metadata + type size
    Ok((data, TAG_LEN as usize + TYPE_ID_LEN as usize + type_size))
}

/// Auxiliary data iterator
///
/// This struct is created by the [`Record::aux_iter`] and
/// [`RecordView::aux_iter`] methods.
///
/// This iterator returns `Result`s that wrap tuples containing
/// a slice which represents the two-byte tag (`&[u8; 2]`) as
/// well as an `Aux` enum that wraps the associated value.
///
/// When an error occurs, the `Err` variant will be returned
/// and the iterator will not be able to advance anymore.
pub struct AuxIter<'a> {
    aux: &'a [u8],
}

impl<'a> Iterator for AuxIter<'a> {
    type Item = Result<(&'a [u8], Aux<'a>)>;

    fn next(&mut self) -> Option<Self::Item> {
        // We're finished
        if self.aux.is_empty() {
            return None;
        }
        // Incomplete aux data
        if (1..=3).contains(&self.aux.len()) {
            // In the case of an error, we can not safely advance in the aux data, so we terminate the Iteration
            self.aux = &[];
            return Some(Err(Error::AuxParsingError));
        }
        let tag = &self.aux[..2];
        // SAFETY: data_ptr points into the record's aux data buffer; parse_aux_field validates the format.
        Some(unsafe {
            let data_ptr = self.aux[2..].as_ptr();
            parse_aux_field(data_ptr)
                .map(|(aux, offset)| {
                    self.aux = &self.aux[offset..];
                    (tag, aux)
                })
                .inspect_err(|_e| {
                    // In the case of an error, we can not safely advance in the aux data, so we terminate the Iteration
                    self.aux = &[];
                })
        })
    }
}

static DECODE_BASE: &[u8] = b"=ACMGRSVTWYHKDBN";

/// Lookup table: maps each packed byte to its two decoded ASCII bases.
/// `DECODE_PAIR[byte] = [high_nibble_base, low_nibble_base]`
static DECODE_PAIR: [[u8; 2]; 256] = {
    const BASE: [u8; 16] = *b"=ACMGRSVTWYHKDBN";
    let mut table = [[0u8; 2]; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = [BASE[i >> 4], BASE[i & 0xf]];
        i += 1;
    }
    table
};
static ENCODE_BASE: [u8; 256] = [
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    1, 2, 4, 8, 15, 15, 15, 15, 15, 15, 15, 15, 15, 0, 15, 15, 15, 1, 14, 2, 13, 15, 15, 4, 11, 15,
    15, 12, 15, 3, 15, 15, 15, 15, 5, 6, 8, 15, 7, 9, 15, 10, 15, 15, 15, 15, 15, 15, 15, 1, 14, 2,
    13, 15, 15, 4, 11, 15, 15, 12, 15, 3, 15, 15, 15, 15, 5, 6, 8, 15, 7, 9, 15, 10, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
];

/// Scalar sequence decoder using the DECODE_PAIR lookup table.
fn decode_seq_scalar(encoded: &[u8], len: usize) -> Vec<u8> {
    let full_bytes = len / 2;
    assert!(encoded.len() >= len.div_ceil(2));
    let mut result = vec![0u8; len];

    for (chunk, &byte) in result[..full_bytes * 2]
        .chunks_exact_mut(2)
        .zip(&encoded[..full_bytes])
    {
        let pair = DECODE_PAIR[byte as usize];
        chunk[0] = pair[0];
        chunk[1] = pair[1];
    }

    if len % 2 == 1 {
        result[len - 1] = DECODE_PAIR[encoded[full_bytes] as usize][0];
    }

    result
}

/// SSSE3 sequence decoder: uses `pshufb` as a 16-entry SIMD LUT to decode
/// 32 bases (16 packed bytes) per iteration.
///
/// # Safety
///
/// Caller must ensure SSSE3 is available on the current CPU.
/// `encoded` must contain at least `(len + 1) / 2` bytes (the standard
/// BAM packing invariant: two bases per byte, with a possible trailing
/// half-byte for odd `len`).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "ssse3")]
unsafe fn decode_seq_ssse3(encoded: &[u8], len: usize) -> Vec<u8> {
    use std::arch::x86_64::*;

    let full_bytes = len / 2;
    assert!(encoded.len() >= len.div_ceil(2));
    let mut result = vec![0u8; len];

    // Load the 16-byte DECODE_BASE table into a SIMD register
    let lut = _mm_loadu_si128(DECODE_BASE.as_ptr() as *const __m128i);
    let mask_lo = _mm_set1_epi8(0x0F);

    let mut i = 0; // index into encoded bytes
    let mut o = 0; // index into result

    // Main SIMD loop: 16 packed bytes → 32 decoded bases per iteration.
    //
    // Bounds reasoning:
    //   - Read:  encoded[i..i+16], guarded by i + 16 <= full_bytes <= encoded.len()
    //   - Write: result[o..o+32],  guarded by o + 32 = 2*(i+16) <= 2*full_bytes <= len
    while i + 16 <= full_bytes {
        let packed = _mm_loadu_si128(encoded.as_ptr().add(i) as *const __m128i);

        // Extract high and low nibbles of each byte.
        // _mm_srli_epi16 shifts 16-bit lanes, but the & 0x0F mask discards
        // any cross-byte bleed, giving correct per-byte high nibbles.
        let hi = _mm_and_si128(_mm_srli_epi16(packed, 4), mask_lo);
        let lo = _mm_and_si128(packed, mask_lo);

        // Decode nibbles → ASCII bases via shuffle LUT
        let decoded_hi = _mm_shuffle_epi8(lut, hi);
        let decoded_lo = _mm_shuffle_epi8(lut, lo);

        // Interleave: [hi0,lo0, hi1,lo1, ..., hi7,lo7] and [hi8,lo8, ..., hi15,lo15]
        let out_a = _mm_unpacklo_epi8(decoded_hi, decoded_lo);
        let out_b = _mm_unpackhi_epi8(decoded_hi, decoded_lo);

        _mm_storeu_si128(result.as_mut_ptr().add(o) as *mut __m128i, out_a);
        _mm_storeu_si128(result.as_mut_ptr().add(o + 16) as *mut __m128i, out_b);

        i += 16;
        o += 32;
    }

    // Scalar tail for remaining full bytes (at most 15 iterations)
    while i < full_bytes {
        let pair = DECODE_PAIR[encoded[i] as usize];
        result[o] = pair[0];
        result[o + 1] = pair[1];
        i += 1;
        o += 2;
    }

    // Trailing odd base
    if len % 2 == 1 {
        result[o] = DECODE_PAIR[encoded[i] as usize][0];
    }

    result
}

/// NEON sequence decoder: uses `vqtbl1q_u8` as a 16-entry SIMD LUT to
/// decode 32 bases (16 packed bytes) per iteration.
///
/// # Safety
///
/// `encoded` must contain at least `(len + 1) / 2` bytes (the standard
/// BAM packing invariant: two bases per byte, with a possible trailing
/// half-byte for odd `len`).
#[cfg(target_arch = "aarch64")]
unsafe fn decode_seq_neon(encoded: &[u8], len: usize) -> Vec<u8> {
    use std::arch::aarch64::*;

    let full_bytes = len / 2;
    assert!(encoded.len() >= len.div_ceil(2));
    let mut result = vec![0u8; len];

    // Load the 16-byte DECODE_BASE table into a NEON register
    let lut = vld1q_u8(DECODE_BASE.as_ptr());
    let mask_lo = vdupq_n_u8(0x0F);

    let mut i = 0;
    let mut o = 0;

    // Main SIMD loop: 16 packed bytes → 32 decoded bases per iteration.
    //
    // Bounds reasoning: same as SSSE3 path above.
    while i + 16 <= full_bytes {
        let packed = vld1q_u8(encoded.as_ptr().add(i));

        // Extract nibbles (vshrq_n_u8 is a true per-byte shift)
        let hi = vshrq_n_u8(packed, 4);
        let lo = vandq_u8(packed, mask_lo);

        // Decode via table lookup
        let decoded_hi = vqtbl1q_u8(lut, hi);
        let decoded_lo = vqtbl1q_u8(lut, lo);

        // Interleave into pairs: [hi0,lo0, hi1,lo1, ...]
        let out_a = vzip1q_u8(decoded_hi, decoded_lo);
        let out_b = vzip2q_u8(decoded_hi, decoded_lo);

        vst1q_u8(result.as_mut_ptr().add(o), out_a);
        vst1q_u8(result.as_mut_ptr().add(o + 16), out_b);

        i += 16;
        o += 32;
    }

    // Scalar tail (at most 15 iterations)
    while i < full_bytes {
        let pair = DECODE_PAIR[encoded[i] as usize];
        result[o] = pair[0];
        result[o + 1] = pair[1];
        i += 1;
        o += 2;
    }

    if len % 2 == 1 {
        result[o] = DECODE_PAIR[encoded[i] as usize][0];
    }

    result
}

#[inline]
fn encoded_base(encoded_seq: &[u8], i: usize) -> u8 {
    (encoded_seq[i / 2] >> ((!i & 1) << 2)) & 0b1111
}

#[inline]
/// # Safety
/// Caller must ensure `i / 2 < encoded_seq.len()`.
unsafe fn encoded_base_unchecked(encoded_seq: &[u8], i: usize) -> u8 {
    // SAFETY: caller guarantees i/2 is in bounds.
    (encoded_seq.get_unchecked(i / 2) >> ((!i & 1) << 2)) & 0b1111
}

#[inline]
fn decode_base_unchecked(base: u8) -> &'static u8 {
    // SAFETY: base is a 4-bit value (0..15) from encoded BAM data; DECODE_BASE has 16 entries.
    unsafe { DECODE_BASE.get_unchecked(base as usize) }
}

/// The sequence of a record.
#[derive(Debug, Copy, Clone)]
pub struct Seq<'a> {
    pub encoded: &'a [u8],
    len: usize,
}

impl Seq<'_> {
    /// Return encoded base. Complexity: O(1).
    #[inline]
    pub fn encoded_base(&self, i: usize) -> u8 {
        encoded_base(self.encoded, i)
    }

    /// Return encoded base. Complexity: O(1).
    ///
    /// # Safety
    ///
    /// TODO
    #[inline]
    pub unsafe fn encoded_base_unchecked(&self, i: usize) -> u8 {
        encoded_base_unchecked(self.encoded, i)
    }

    /// Obtain decoded base without performing bounds checking.
    /// Use index based access seq()[i], for checked, safe access.
    /// Complexity: O(1).
    ///
    /// # Safety
    ///
    /// TODO
    #[inline]
    pub unsafe fn decoded_base_unchecked(&self, i: usize) -> u8 {
        *decode_base_unchecked(self.encoded_base_unchecked(i))
    }

    /// Return decoded sequence. Complexity: O(m) with m being the read length.
    ///
    /// Uses SIMD acceleration when available (SSSE3 on x86_64, NEON on aarch64).
    pub fn as_bytes(&self) -> Vec<u8> {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("ssse3") {
                // SAFETY: we just verified SSSE3 is available.
                return unsafe { decode_seq_ssse3(self.encoded, self.len) };
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
            // NEON is always available on aarch64.
            return unsafe { decode_seq_neon(self.encoded, self.len) };
        }

        #[allow(unreachable_code)]
        decode_seq_scalar(self.encoded, self.len)
    }

    /// Return length (in bases) of the sequence.
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Return decoded base at given position, or `None` if out of bounds.
    /// Complexity: O(1).
    #[inline]
    pub fn get(&self, index: usize) -> Option<u8> {
        if index < self.len {
            Some(*decode_base_unchecked(self.encoded_base(index)))
        } else {
            None
        }
    }
}

impl ops::Index<usize> for Seq<'_> {
    type Output = u8;

    /// Return decoded base at given position within read. Complexity: O(1).
    fn index(&self, index: usize) -> &u8 {
        decode_base_unchecked(self.encoded_base(index))
    }
}

// SAFETY: Seq borrows only an immutable &[u8] slice from the record's data buffer.
unsafe impl Send for Seq<'_> {}
unsafe impl Sync for Seq<'_> {}

/// A borrowed, read-only view of a BAM record.
///
/// Zero-copy alternative to cloning via `Record::from_inner`. Exposes the most
/// commonly-needed fields from a pileup alignment without heap allocation.
///
/// For fields not available here (chromosome ID, position, insert size, and most
/// flag predicates), use [`crate::bam::pileup::Alignment::record`] instead, which
/// returns a full owned [`Record`].
pub struct RecordView<'a> {
    inner: &'a htslib::bam1_t,
}

impl<'a> std::fmt::Debug for RecordView<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RecordView")
            .field("qname", &self.qname())
            .field("seq_len", &self.seq_len())
            .field("flags", &self.flags())
            .finish()
    }
}

impl<'a> RecordView<'a> {
    const BAM_FREVERSE: u16 = 16;
    const BAM_FMREVERSE: u16 = 32;
    const BAM_FREAD1: u16 = 64;
    const BAM_FREAD2: u16 = 128;

    /// Create a `RecordView` from a raw `bam1_t` pointer.
    ///
    /// # Safety
    /// The pointer must be valid and the referenced `bam1_t` (including its `data` buffer)
    /// must live for at least `'a`.
    pub unsafe fn from_raw(ptr: *const htslib::bam1_t) -> RecordView<'a> {
        RecordView { inner: &*ptr }
    }

    fn data(&self) -> &'a [u8] {
        // SAFETY: inner.data is valid for inner.l_data bytes (maintained by htslib); lifetime tied to 'a.
        unsafe { slice::from_raw_parts(self.inner.data, self.inner.l_data as usize) }
    }

    fn qname_capacity(&self) -> usize {
        self.inner.core.l_qname as usize
    }

    fn qname_len(&self) -> usize {
        self.qname_capacity() - 1 - self.inner.core.l_extranul as usize
    }

    fn cigar_len(&self) -> usize {
        self.inner.core.n_cigar as usize
    }

    fn seq_data(&self) -> &'a [u8] {
        let offset = self.qname_capacity() + self.cigar_len() * 4;
        &self.data()[offset..][..self.seq_len().div_ceil(2)]
    }

    /// Get qname (read name). Complexity: O(1).
    pub fn qname(&self) -> &'a [u8] {
        &self.data()[..self.qname_len()]
    }

    /// Get reference to raw cigar string representation.
    pub fn raw_cigar(&self) -> &'a [u32] {
        // SAFETY: cigar data starts at a 4-byte-aligned offset (qname is padded); length from n_cigar.
        #[allow(clippy::cast_ptr_alignment)]
        unsafe {
            slice::from_raw_parts(
                self.data()[self.qname_capacity()..].as_ptr() as *const u32,
                self.cigar_len(),
            )
        }
    }

    /// Get read sequence. Complexity: O(1).
    pub fn seq(&self) -> Seq<'a> {
        Seq {
            encoded: self.seq_data(),
            len: self.seq_len(),
        }
    }

    /// Get base qualities. Complexity: O(1).
    pub fn qual(&self) -> &'a [u8] {
        &self.data()[self.qname_capacity() + self.cigar_len() * 4 + self.seq_len().div_ceil(2)..]
            [..self.seq_len()]
    }

    /// Get read sequence and base qualities in a single pass over the record layout.
    ///
    /// Equivalent to calling `seq()` and `qual()` separately, but computes the
    /// internal byte offsets only once. Prefer this when you need both.
    pub fn seq_and_qual(&self) -> (Seq<'a>, &'a [u8]) {
        let seq_off = self.qname_capacity() + self.cigar_len() * 4;
        let seq_len = self.seq_len();
        let data = self.data();
        let seq = Seq {
            encoded: &data[seq_off..][..seq_len.div_ceil(2)],
            len: seq_len,
        };
        let qual = &data[seq_off + seq_len.div_ceil(2)..][..seq_len];
        (seq, qual)
    }

    /// Get MAPQ.
    pub fn mapq(&self) -> u8 {
        self.inner.core.qual
    }

    /// Get raw flags.
    pub fn flags(&self) -> u16 {
        self.inner.core.flag
    }

    /// Get sequence length.
    pub fn seq_len(&self) -> usize {
        self.inner.core.l_qseq as usize
    }

    /// Returns true if the record is on the reverse strand.
    pub fn is_reverse(&self) -> bool {
        self.flags() & Self::BAM_FREVERSE != 0
    }

    /// Returns true if this is the first segment in the template.
    pub fn is_first_in_template(&self) -> bool {
        self.flags() & Self::BAM_FREAD1 != 0
    }

    /// Returns true if this is the last segment in the template.
    pub fn is_last_in_template(&self) -> bool {
        self.flags() & Self::BAM_FREAD2 != 0
    }

    /// Returns true if the mate is on the reverse strand.
    pub fn is_mate_reverse(&self) -> bool {
        self.flags() & Self::BAM_FMREVERSE != 0
    }

    /// Get the raw auxiliary data as a byte slice.
    ///
    /// Returns `None` if the computed offset exceeds the record's data length,
    /// which can happen with malformed BAM records. For well-formed records
    /// this always returns `Some`.
    ///
    /// The returned slice is the unparsed aux segment of the BAM record,
    /// suitable for fast, custom tag scanning without FFI or full type dispatch
    /// overhead.
    pub fn raw_aux_data(&self) -> Option<&'a [u8]> {
        let seq_len = self.seq_len();
        let offset = self.qname_capacity() + self.cigar_len() * 4 + seq_len.div_ceil(2) + seq_len;
        self.data().get(offset..)
    }

    /// Look up an auxiliary field by its tag.
    ///
    /// Only the first two bytes of a given tag are used for the look-up of a field.
    /// See [`Aux`] for more details.
    pub fn aux(&self, tag: &[u8]) -> Result<Aux<'a>> {
        if tag.len() < 2 {
            return Err(Error::AuxStringError);
        }
        let raw = self.raw_aux_data().unwrap_or(&[]);
        match aux_tag_search(raw, tag) {
            Some(offset) => {
                // SAFETY: offset is within raw (returned by aux_tag_search); raw is valid for 'a.
                unsafe { parse_aux_field(raw.as_ptr().add(offset)).map(|(v, _)| v) }
            }
            None => Err(Error::AuxTagNotFound),
        }
    }

    /// Returns an iterator over the auxiliary fields of the record.
    ///
    /// When an error occurs, the `Err` variant will be returned
    /// and the iterator will not be able to advance anymore.
    pub fn aux_iter(&self) -> AuxIter<'a> {
        AuxIter {
            aux: self.raw_aux_data().unwrap_or(&[]),
        }
    }
}

#[cfg_attr(feature = "serde_feature", derive(Serialize, Deserialize))]
#[derive(PartialEq, PartialOrd, Eq, Debug, Clone, Copy, Hash)]
pub enum Cigar {
    Match(u32),    // M
    Ins(u32),      // I
    Del(u32),      // D
    RefSkip(u32),  // N
    SoftClip(u32), // S
    HardClip(u32), // H
    Pad(u32),      // P
    Equal(u32),    // =
    Diff(u32),     // X
}

/// BAM CIGAR operation codes and utilities for working with raw-encoded
/// CIGAR values (u32 where low 4 bits = op, upper 28 bits = length).
pub mod cigar_op {
    /// M — alignment match (can be a sequence match or mismatch)
    pub const MATCH: u32 = 0;
    /// I — insertion to the reference
    pub const INS: u32 = 1;
    /// D — deletion from the reference
    pub const DEL: u32 = 2;
    /// N — skipped region from the reference (e.g. intron)
    pub const REF_SKIP: u32 = 3;
    /// S — soft clipping (clipped sequences present in SEQ)
    pub const SOFT_CLIP: u32 = 4;
    /// H — hard clipping (clipped sequences NOT present in SEQ)
    pub const HARD_CLIP: u32 = 5;
    /// P — padding (silent deletion from padded reference)
    pub const PAD: u32 = 6;
    /// = — sequence match
    pub const EQUAL: u32 = 7;
    /// X — sequence mismatch
    pub const DIFF: u32 = 8;

    /// Bitmask of ops that consume query/read bases: M(0), I(1), S(4), =(7), X(8).
    const QUERY_MASK: u16 =
        (1 << MATCH) | (1 << INS) | (1 << SOFT_CLIP) | (1 << EQUAL) | (1 << DIFF);
    /// Bitmask of ops that consume reference bases: M(0), D(2), N(3), =(7), X(8).
    const REF_MASK: u16 = (1 << MATCH) | (1 << DEL) | (1 << REF_SKIP) | (1 << EQUAL) | (1 << DIFF);

    /// Whether a raw CIGAR operation code consumes query/read bases.
    #[inline(always)]
    pub fn consumes_query(op: u32) -> bool {
        (QUERY_MASK >> op) & 1 != 0
    }

    /// Whether a raw CIGAR operation code consumes reference bases.
    #[inline(always)]
    pub fn consumes_ref(op: u32) -> bool {
        (REF_MASK >> op) & 1 != 0
    }

    /// Extract the operation code (low 4 bits) from a raw CIGAR u32.
    #[inline(always)]
    pub fn op(raw: u32) -> u32 {
        raw & 0xF
    }

    /// Extract the length (upper 28 bits) from a raw CIGAR u32.
    #[inline(always)]
    pub fn len(raw: u32) -> u32 {
        raw >> 4
    }
}

impl Cigar {
    /// Decode a raw BAM CIGAR u32 into a `Cigar` variant.
    ///
    /// In BAM format, each CIGAR operation is stored as a u32 where the
    /// low 4 bits encode the operation type and the upper 28 bits encode
    /// the length.
    ///
    /// # Panics
    ///
    /// Panics if the operation code (low 4 bits) is not in the range 0..=8.
    #[inline]
    pub fn from_raw(raw: u32) -> Self {
        let len = raw >> 4;
        match raw & 0xF {
            0 => Cigar::Match(len),
            1 => Cigar::Ins(len),
            2 => Cigar::Del(len),
            3 => Cigar::RefSkip(len),
            4 => Cigar::SoftClip(len),
            5 => Cigar::HardClip(len),
            6 => Cigar::Pad(len),
            7 => Cigar::Equal(len),
            8 => Cigar::Diff(len),
            op => panic!("unexpected CIGAR operation code: {}", op),
        }
    }

    fn encode(self) -> u32 {
        match self {
            Cigar::Match(len) => len << 4, // | 0,
            Cigar::Ins(len) => (len << 4) | 1,
            Cigar::Del(len) => (len << 4) | 2,
            Cigar::RefSkip(len) => (len << 4) | 3,
            Cigar::SoftClip(len) => (len << 4) | 4,
            Cigar::HardClip(len) => (len << 4) | 5,
            Cigar::Pad(len) => (len << 4) | 6,
            Cigar::Equal(len) => (len << 4) | 7,
            Cigar::Diff(len) => (len << 4) | 8,
        }
    }

    /// Return the length of the CIGAR.
    pub fn len(self) -> u32 {
        match self {
            Cigar::Match(len) => len,
            Cigar::Ins(len) => len,
            Cigar::Del(len) => len,
            Cigar::RefSkip(len) => len,
            Cigar::SoftClip(len) => len,
            Cigar::HardClip(len) => len,
            Cigar::Pad(len) => len,
            Cigar::Equal(len) => len,
            Cigar::Diff(len) => len,
        }
    }

    pub fn is_empty(self) -> bool {
        self.len() == 0
    }

    /// Return the character representing the CIGAR.
    pub fn char(self) -> char {
        match self {
            Cigar::Match(_) => 'M',
            Cigar::Ins(_) => 'I',
            Cigar::Del(_) => 'D',
            Cigar::RefSkip(_) => 'N',
            Cigar::SoftClip(_) => 'S',
            Cigar::HardClip(_) => 'H',
            Cigar::Pad(_) => 'P',
            Cigar::Equal(_) => '=',
            Cigar::Diff(_) => 'X',
        }
    }
}

impl fmt::Display for Cigar {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.write_fmt(format_args!("{}{}", self.len(), self.char()))
    }
}

// SAFETY: Cigar is a simple Copy enum with no pointers or interior mutability.
unsafe impl Send for Cigar {}
unsafe impl Sync for Cigar {}

/// A CIGAR string
///
/// Backed by an `Arc<[Cigar]>` for cheap cloning.
///
/// # Example
///
/// ```
/// use rust_htslib::bam::record::{Cigar, CigarString};
///
/// let cigar = CigarString::from(vec![Cigar::Match(100), Cigar::SoftClip(10)]);
///
/// // access by index
/// assert_eq!(cigar[0], Cigar::Match(100));
/// // format into classical string representation
/// assert_eq!(format!("{}", cigar), "100M10S");
/// // iterate
/// for op in &cigar {
///    println!("{}", op);
/// }
/// // cheap clone (just a reference count bump)
/// let cigar2 = cigar.clone();
/// ```
#[derive(PartialEq, PartialOrd, Eq, Clone, Hash, Debug)]
pub struct CigarString(pub Arc<[Cigar]>);

#[cfg(feature = "serde_feature")]
impl Serialize for CigarString {
    fn serialize<S: serde::Serializer>(
        &self,
        serializer: S,
    ) -> std::result::Result<S::Ok, S::Error> {
        self.0.as_ref().serialize(serializer)
    }
}

#[cfg(feature = "serde_feature")]
impl<'de> Deserialize<'de> for CigarString {
    fn deserialize<D: serde::Deserializer<'de>>(
        deserializer: D,
    ) -> std::result::Result<Self, D::Error> {
        let v = Vec::<Cigar>::deserialize(deserializer)?;
        Ok(CigarString::from(v))
    }
}

impl ops::Deref for CigarString {
    type Target = [Cigar];

    fn deref(&self) -> &[Cigar] {
        &self.0
    }
}

impl ops::Index<usize> for CigarString {
    type Output = Cigar;

    fn index(&self, index: usize) -> &Cigar {
        &self.0[index]
    }
}

impl From<Vec<Cigar>> for CigarString {
    fn from(v: Vec<Cigar>) -> Self {
        CigarString(v.into())
    }
}

impl<const N: usize> From<[Cigar; N]> for CigarString {
    fn from(a: [Cigar; N]) -> Self {
        CigarString(Arc::from(a))
    }
}

impl std::iter::FromIterator<Cigar> for CigarString {
    fn from_iter<I: IntoIterator<Item = Cigar>>(iter: I) -> Self {
        CigarString(iter.into_iter().collect())
    }
}

impl CigarString {
    /// Create a `CigarStringView` from this CigarString at position `pos`
    pub fn into_view(self, pos: i64) -> CigarStringView {
        CigarStringView::new(self, pos)
    }

    /// Calculate the bam cigar from the alignment struct. x is the target string
    /// and y is the reference. `hard_clip` controls how unaligned read bases are encoded in the
    /// cigar string. Set to true to use the hard clip (`H`) code, or false to use soft clip
    /// (`S`) code. See the [SAM spec](https://samtools.github.io/hts-specs/SAMv1.pdf) for more details.
    pub fn from_alignment(alignment: &Alignment, hard_clip: bool) -> Self {
        match alignment.mode {
            AlignmentMode::Global => {
                panic!(" Bam cigar fn not supported for Global Alignment mode")
            }
            AlignmentMode::Local => panic!(" Bam cigar fn not supported for Local Alignment mode"),
            _ => {}
        }

        let mut cigar = Vec::new();
        if alignment.operations.is_empty() {
            return CigarString::from(cigar);
        }

        let add_op = |op: AlignmentOperation, length: u32, cigar: &mut Vec<Cigar>| match op {
            AlignmentOperation::Del => cigar.push(Cigar::Del(length)),
            AlignmentOperation::Ins => cigar.push(Cigar::Ins(length)),
            AlignmentOperation::Subst => cigar.push(Cigar::Diff(length)),
            AlignmentOperation::Match => cigar.push(Cigar::Equal(length)),
            _ => {}
        };

        if alignment.xstart > 0 {
            cigar.push(if hard_clip {
                Cigar::HardClip(alignment.xstart as u32)
            } else {
                Cigar::SoftClip(alignment.xstart as u32)
            });
        }

        let mut last = alignment.operations[0];
        let mut k = 1u32;
        for &op in alignment.operations[1..].iter() {
            if op == last {
                k += 1;
            } else {
                add_op(last, k, &mut cigar);
                k = 1;
            }
            last = op;
        }
        add_op(last, k, &mut cigar);
        if alignment.xlen > alignment.xend {
            cigar.push(if hard_clip {
                Cigar::HardClip((alignment.xlen - alignment.xend) as u32)
            } else {
                Cigar::SoftClip((alignment.xlen - alignment.xend) as u32)
            });
        }

        CigarString::from(cigar)
    }
}

impl TryFrom<&[u8]> for CigarString {
    type Error = Error;

    /// Create a CigarString from given &[u8].
    /// # Example
    /// ```
    /// use rust_htslib::bam::record::*;
    /// use rust_htslib::bam::record::CigarString;
    /// use rust_htslib::bam::record::Cigar::*;
    /// use std::convert::TryFrom;
    ///
    /// let cigar_str = "2H10M5X3=2H".as_bytes();
    /// let cigar = CigarString::try_from(cigar_str)
    ///     .expect("Unable to parse cigar string.");
    /// let expected_cigar = CigarString::from(vec![
    ///     HardClip(2),
    ///     Match(10),
    ///     Diff(5),
    ///     Equal(3),
    ///     HardClip(2),
    /// ]);
    /// assert_eq!(cigar, expected_cigar);
    /// ```
    fn try_from(bytes: &[u8]) -> Result<Self> {
        let mut inner = Vec::new();
        let mut i = 0;
        let text_len = bytes.len();
        while i < text_len {
            let mut j = i;
            while j < text_len && bytes[j].is_ascii_digit() {
                j += 1;
            }
            // check that length is provided
            if i == j {
                return Err(Error::ParseCigar {
                    msg: "Expected length before cigar operation [0-9]+[MIDNSHP=X]".to_owned(),
                });
            }
            // get the length of the operation
            let s = str::from_utf8(&bytes[i..j]).map_err(|_| Error::ParseCigar {
                msg: format!("Invalid utf-8 bytes '{:?}'.", &bytes[i..j]),
            })?;
            let n = s.parse().map_err(|_| Error::ParseCigar {
                msg: format!("Unable to parse &str '{:?}' to u32.", s),
            })?;
            // get the operation
            let op = &bytes[j];
            inner.push(match op {
                b'M' => Cigar::Match(n),
                b'I' => Cigar::Ins(n),
                b'D' => Cigar::Del(n),
                b'N' => Cigar::RefSkip(n),
                b'H' => {
                    if i == 0 || j + 1 == text_len {
                        Cigar::HardClip(n)
                    } else {
                        return Err(Error::ParseCigar {
                            msg: "Hard clipping ('H') is only valid at the start or end of a cigar."
                                .to_owned(),
                        });
                    }
                }
                b'S' => {
                    if i == 0
                        || j + 1 == text_len
                        || bytes[i-1] == b'H'
                        || bytes[j+1..].iter().all(|c| c.is_ascii_digit() || *c == b'H') {
                        Cigar::SoftClip(n)
                    } else {
                        return Err(Error::ParseCigar {
                        msg: "Soft clips ('S') can only have hard clips ('H') between them and the end of the CIGAR string."
                            .to_owned(),
                        });
                    }
                },
                b'P' => Cigar::Pad(n),
                b'=' => Cigar::Equal(n),
                b'X' => Cigar::Diff(n),
                op => {
                    return Err(Error::ParseCigar {
                        msg: format!("Expected cigar operation [MIDNSHP=X] but got [{}]", op),
                    })
                }
            });
            i = j + 1;
        }
        Ok(CigarString::from(inner))
    }
}

impl TryFrom<&str> for CigarString {
    type Error = Error;

    /// Create a CigarString from given &str.
    /// # Example
    /// ```
    /// use rust_htslib::bam::record::*;
    /// use rust_htslib::bam::record::CigarString;
    /// use rust_htslib::bam::record::Cigar::*;
    /// use std::convert::TryFrom;
    ///
    /// let cigar_str = "2H10M5X3=2H";
    /// let cigar = CigarString::try_from(cigar_str)
    ///     .expect("Unable to parse cigar string.");
    /// let expected_cigar = CigarString::from(vec![
    ///     HardClip(2),
    ///     Match(10),
    ///     Diff(5),
    ///     Equal(3),
    ///     HardClip(2),
    /// ]);
    /// assert_eq!(cigar, expected_cigar);
    /// ```
    fn try_from(text: &str) -> Result<Self> {
        let bytes = text.as_bytes();
        if !text.is_ascii() {
            return Err(Error::ParseCigar {
                msg: "CIGAR string contained non-ASCII characters, which are not valid. Valid are [0-9MIDNSHP=X].".to_owned(),
            });
        }
        CigarString::try_from(bytes)
    }
}

impl<'a> IntoIterator for &'a CigarString {
    type Item = &'a Cigar;
    type IntoIter = std::slice::Iter<'a, Cigar>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl fmt::Display for CigarString {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        for op in self.iter() {
            fmt.write_fmt(format_args!("{}{}", op.len(), op.char()))?;
        }
        Ok(())
    }
}

// Get number of leading/trailing softclips if a CigarString taking hardclips into account
fn calc_softclips<'a>(mut cigar: impl DoubleEndedIterator<Item = &'a Cigar>) -> i64 {
    match (cigar.next(), cigar.next()) {
        (Some(Cigar::HardClip(_)), Some(Cigar::SoftClip(s))) | (Some(Cigar::SoftClip(s)), _) => {
            *s as i64
        }
        _ => 0,
    }
}

#[derive(Eq, PartialEq, Clone, Debug)]
pub struct CigarStringView {
    inner: CigarString,
    pos: i64,
}

impl CigarStringView {
    /// Construct a new CigarStringView from a CigarString at a position
    pub fn new(c: CigarString, pos: i64) -> CigarStringView {
        CigarStringView { inner: c, pos }
    }

    /// Get (exclusive) end position of alignment.
    pub fn end_pos(&self) -> i64 {
        let mut pos = self.pos;
        for c in self {
            match c {
                Cigar::Match(l)
                | Cigar::RefSkip(l)
                | Cigar::Del(l)
                | Cigar::Equal(l)
                | Cigar::Diff(l) => pos += *l as i64,
                // these don't add to end_pos on reference
                Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => (),
            }
        }
        pos
    }

    /// Get the start position of the alignment (0-based).
    pub fn pos(&self) -> i64 {
        self.pos
    }

    /// Get number of bases softclipped at the beginning of the alignment.
    pub fn leading_softclips(&self) -> i64 {
        calc_softclips(self.iter())
    }

    /// Get number of bases softclipped at the end of the alignment.
    pub fn trailing_softclips(&self) -> i64 {
        calc_softclips(self.iter().rev())
    }

    /// Get number of bases hardclipped at the beginning of the alignment.
    pub fn leading_hardclips(&self) -> i64 {
        self.first().map_or(0, |cigar| {
            if let Cigar::HardClip(s) = cigar {
                *s as i64
            } else {
                0
            }
        })
    }

    /// Get number of bases hardclipped at the end of the alignment.
    pub fn trailing_hardclips(&self) -> i64 {
        self.last().map_or(0, |cigar| {
            if let Cigar::HardClip(s) = cigar {
                *s as i64
            } else {
                0
            }
        })
    }

    /// For a given position in the reference, get corresponding position within read.
    /// If reference position is outside of the read alignment, return None.
    ///
    /// # Arguments
    ///
    /// * `ref_pos` - the reference position
    /// * `include_softclips` - if true, softclips will be considered as matches or mismatches
    /// * `include_dels` - if true, positions within deletions will be considered (first reference matching read position after deletion will be returned)
    ///
    pub fn read_pos(
        &self,
        ref_pos: u32,
        include_softclips: bool,
        include_dels: bool,
    ) -> Result<Option<u32>> {
        let mut rpos = self.pos as u32; // reference position
        let mut qpos = 0u32; // position within read
        let mut j = 0; // index into cigar operation vector

        // find first cigar operation referring to qpos = 0 (and thus bases in record.seq()),
        // because all augmentations of qpos and rpos before that are invalid
        for (i, c) in self.iter().enumerate() {
            match c {
                Cigar::Match(_) |
                Cigar::Diff(_)  |
                Cigar::Equal(_) |
                // this is unexpected, but bwa + GATK indel realignment can produce insertions
                // before matching positions
                Cigar::Ins(_) => {
                    j = i;
                    break;
                },
                Cigar::SoftClip(l) => {
                    j = i;
                    if include_softclips {
                        // Alignment starts with softclip and we want to include it in the
                        // projection of the reference position. However, the POS field does not
                        // include the softclip. Hence we have to subtract its length.
                        rpos = rpos.saturating_sub(*l);
                    }
                    break;
                },
                Cigar::Del(l) => {
                    // METHOD: leading deletions can happen in case of trimmed reads where
                    // a primer has been removed AFTER read mapping.
                    // Example: 24M8I8D18M9S before trimming, 32H8D18M9S after trimming
                    // with fgbio. While leading deletions should be impossible with
                    // normal read mapping, they make perfect sense with primer trimming
                    // because the mapper still had the evidence to decide in favor of
                    // the deletion via the primer sequence.
                    rpos += l;
                },
                Cigar::RefSkip(_) => {
                    return Err(Error::UnexpectedCigarOperation {
                        msg: "'reference skip' (N) found before any operation describing read sequence".to_owned()
                    });
                },
                Cigar::HardClip(_) if i > 0 && i < self.len()-1 => {
                    return Err(Error::UnexpectedCigarOperation{
                        msg: "'hard clip' (H) found in between operations, contradicting SAMv1 spec that hard clips can only be at the ends of reads".to_owned()
                    });
                },
                // if we have reached the end of the CigarString with only pads and hard clips, we have no read position matching the variant
                Cigar::Pad(_) | Cigar::HardClip(_) if i == self.len()-1 => return Ok(None),
                // skip leading HardClips and Pads, as they consume neither read sequence nor reference sequence
                Cigar::Pad(_) | Cigar::HardClip(_) => ()
            }
        }

        let contains_ref_pos = |cigar_op_start: u32, cigar_op_length: u32| {
            cigar_op_start <= ref_pos && cigar_op_start + cigar_op_length > ref_pos
        };

        while rpos <= ref_pos && j < self.len() {
            match self[j] {
                // potential SNV evidence
                Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) if contains_ref_pos(rpos, l) => {
                    // difference between desired position and first position of current cigar
                    // operation
                    qpos += ref_pos - rpos;
                    return Ok(Some(qpos));
                }
                Cigar::SoftClip(l) if include_softclips && contains_ref_pos(rpos, l) => {
                    qpos += ref_pos - rpos;
                    return Ok(Some(qpos));
                }
                Cigar::Del(l) if include_dels && contains_ref_pos(rpos, l) => {
                    // qpos shall resemble the start of the deletion
                    return Ok(Some(qpos));
                }
                // for others, just increase pos and qpos as needed
                Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) => {
                    rpos += l;
                    qpos += l;
                    j += 1;
                }
                Cigar::SoftClip(l) => {
                    qpos += l;
                    j += 1;
                    if include_softclips {
                        rpos += l;
                    }
                }
                Cigar::Ins(l) => {
                    qpos += l;
                    j += 1;
                }
                Cigar::RefSkip(l) | Cigar::Del(l) => {
                    rpos += l;
                    j += 1;
                }
                Cigar::Pad(_) => {
                    j += 1;
                }
                Cigar::HardClip(_) if j < self.len() - 1 => {
                    return Err(Error::UnexpectedCigarOperation{
                        msg: "'hard clip' (H) found in between operations, contradicting SAMv1 spec that hard clips can only be at the ends of reads".to_owned()
                    });
                }
                Cigar::HardClip(_) => return Ok(None),
            }
        }

        Ok(None)
    }

    /// transfer ownership of the Cigar out of the CigarView
    pub fn take(self) -> CigarString {
        self.inner
    }
}

impl ops::Deref for CigarStringView {
    type Target = CigarString;

    fn deref(&self) -> &CigarString {
        &self.inner
    }
}

impl ops::Index<usize> for CigarStringView {
    type Output = Cigar;

    fn index(&self, index: usize) -> &Cigar {
        &self.inner[index]
    }
}

impl<'a> IntoIterator for &'a CigarStringView {
    type Item = &'a Cigar;
    type IntoIter = std::slice::Iter<'a, Cigar>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.iter()
    }
}

impl fmt::Display for CigarStringView {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.inner.fmt(fmt)
    }
}

pub struct BaseModificationMetadata {
    pub strand: i32,
    pub implicit: i32,
    pub canonical: u8,
}

/// struct containing the internal state required to access
/// the base modifications for a bam::Record
pub struct BaseModificationState<'a> {
    record: &'a Record,
    state: *mut htslib::hts_base_mod_state,
    buffer: Vec<htslib::hts_base_mod>,
    buffer_pos: i32,
}

impl BaseModificationState<'_> {
    /// Initialize a new BaseModification struct from a bam::Record
    /// This function allocates memory for the state structure
    /// and initializes the iterator to the start of the modification
    /// records.
    fn new(r: &Record) -> Result<BaseModificationState<'_>> {
        // SAFETY: hts_base_mod_state_alloc returns a valid pointer or null (checked below).
        let mut bm = unsafe {
            BaseModificationState {
                record: r,
                state: hts_sys::hts_base_mod_state_alloc(),
                buffer: Vec::new(),
                buffer_pos: -1,
            }
        };

        if bm.state.is_null() {
            panic!("Unable to allocate memory for hts_base_mod_state");
        }

        // parse the MM tag to initialize the state
        // SAFETY: bm.record.inner_ptr() is valid; bm.state is non-null (checked above).
        unsafe {
            let ret = hts_sys::bam_parse_basemod(bm.record.inner_ptr(), bm.state);
            if ret != 0 {
                return Err(Error::BaseModificationTagNotFound);
            }
        }

        let types = bm.recorded();
        bm.buffer.reserve(types.len());
        Ok(bm)
    }

    pub fn buffer_next_mods(&mut self) -> Result<usize> {
        // SAFETY: record and state pointers are valid (from constructor); buffer has sufficient capacity.
        unsafe {
            let ret = hts_sys::bam_next_basemod(
                self.record.inner_ptr(),
                self.state,
                self.buffer.as_mut_ptr(),
                self.buffer.capacity() as i32,
                &mut self.buffer_pos,
            );

            if ret < 0 {
                return Err(Error::BaseModificationIterationFailed);
            }

            // the htslib API won't write more than buffer.capacity() mods to the output array but it will
            // return the actual number of modifications found. We return an error to the caller
            // in the case where there was insufficient storage to return all mods.
            if ret as usize > self.buffer.capacity() {
                return Err(Error::BaseModificationTooManyMods);
            }

            // we read the modifications directly into the vector, which does
            // not update the length so needs to be manually set
            self.buffer.set_len(ret as usize);

            Ok(ret as usize)
        }
    }

    /// Return an array containing the modification codes listed for this record.
    /// Positive values are ascii character codes (eg m), negative values are chEBI codes.
    pub fn recorded<'a>(&self) -> &'a [i32] {
        // SAFETY: self.state is non-null (from constructor); bam_mods_recorded returns pointer into state's data.
        unsafe {
            let mut n: i32 = 0;
            let data_ptr: *const i32 = hts_sys::bam_mods_recorded(self.state, &mut n);

            // htslib should not return a null pointer, even when there are no base mods
            if data_ptr.is_null() {
                panic!("Unable to obtain pointer to base modifications");
            }
            assert!(n >= 0);
            // SAFETY: data_ptr is non-null (checked above); n is the count returned by htslib.
            slice::from_raw_parts(data_ptr, n as usize)
        }
    }

    /// Return metadata for the specified character code indicating the strand
    /// the base modification was called on, whether the tag uses implicit mode
    /// and the ascii code for the canonical base.
    /// If there are multiple modifications with the same code this will return the data
    /// for the first mod.  See https://github.com/samtools/htslib/issues/1635
    pub fn query_type(&self, code: i32) -> Result<BaseModificationMetadata> {
        // SAFETY: self.state is non-null (from constructor); output pointers are valid stack references.
        unsafe {
            let mut strand: i32 = 0;
            let mut implicit: i32 = 0;
            // This may be i8 or u8 in hts_sys.
            let mut canonical: c_char = 0;

            let ret = hts_sys::bam_mods_query_type(
                self.state,
                code,
                &mut strand,
                &mut implicit,
                &mut canonical,
            );
            if ret == -1 {
                Err(Error::BaseModificationTypeNotFound)
            } else {
                Ok(BaseModificationMetadata {
                    strand,
                    implicit,
                    canonical: canonical.try_into().unwrap(),
                })
            }
        }
    }
}

impl Drop for BaseModificationState<'_> {
    fn drop<'a>(&mut self) {
        // SAFETY: self.state was allocated by hts_base_mod_state_alloc; free is symmetric.
        unsafe {
            hts_sys::hts_base_mod_state_free(self.state);
        }
    }
}

/// Iterator over the base modifications that returns
/// a vector for all of the mods at each position
pub struct BaseModificationsPositionIter<'a> {
    mod_state: BaseModificationState<'a>,
}

impl BaseModificationsPositionIter<'_> {
    fn new(r: &Record) -> Result<BaseModificationsPositionIter<'_>> {
        let state = BaseModificationState::new(r)?;
        Ok(BaseModificationsPositionIter { mod_state: state })
    }

    pub fn recorded<'a>(&self) -> &'a [i32] {
        self.mod_state.recorded()
    }

    pub fn query_type(&self, code: i32) -> Result<BaseModificationMetadata> {
        self.mod_state.query_type(code)
    }
}

impl Iterator for BaseModificationsPositionIter<'_> {
    type Item = Result<(i32, Vec<hts_sys::hts_base_mod>)>;

    fn next(&mut self) -> Option<Self::Item> {
        let ret = self.mod_state.buffer_next_mods();

        // Three possible things happened in buffer_next_mods:
        // 1. the htslib API call was successful but there are no more mods
        // 2. ths htslib API call was successful and we read some mods
        // 3. the htslib API call failed, we propogate the error wrapped in an option
        match ret {
            Ok(num_mods) => {
                if num_mods == 0 {
                    None
                } else {
                    let data = (self.mod_state.buffer_pos, self.mod_state.buffer.clone());
                    Some(Ok(data))
                }
            }
            Err(e) => Some(Err(e)),
        }
    }
}

/// Iterator over the base modifications that returns
/// the next modification found, one by one
pub struct BaseModificationsIter<'a> {
    mod_state: BaseModificationState<'a>,
    buffer_idx: usize,
}

impl BaseModificationsIter<'_> {
    fn new(r: &Record) -> Result<BaseModificationsIter<'_>> {
        let state = BaseModificationState::new(r)?;
        Ok(BaseModificationsIter {
            mod_state: state,
            buffer_idx: 0,
        })
    }

    pub fn recorded<'a>(&self) -> &'a [i32] {
        self.mod_state.recorded()
    }

    pub fn query_type(&self, code: i32) -> Result<BaseModificationMetadata> {
        self.mod_state.query_type(code)
    }
}

impl Iterator for BaseModificationsIter<'_> {
    type Item = Result<(i32, hts_sys::hts_base_mod)>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.buffer_idx == self.mod_state.buffer.len() {
            // need to use the internal state to read the next
            // set of modifications into the buffer
            let ret = self.mod_state.buffer_next_mods();

            match ret {
                Ok(num_mods) => {
                    if num_mods == 0 {
                        // done iterating
                        return None;
                    } else {
                        // we read some mods, reset the position in the buffer then fall through
                        self.buffer_idx = 0;
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        // if we got here when there are mods buffered that we haven't emitted yet
        assert!(self.buffer_idx < self.mod_state.buffer.len());
        let data = (
            self.mod_state.buffer_pos,
            self.mod_state.buffer[self.buffer_idx],
        );
        self.buffer_idx += 1;
        Some(Ok(data))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_string() {
        let cigar = CigarString::from(vec![Cigar::Match(100), Cigar::SoftClip(10)]);

        assert_eq!(cigar[0], Cigar::Match(100));
        assert_eq!(format!("{}", cigar), "100M10S");
        for op in &cigar {
            println!("{}", op);
        }
    }

    #[test]
    fn test_cigar_string_view_pos() {
        let cigar = CigarString::from(vec![Cigar::Match(100), Cigar::SoftClip(10)]).into_view(5);
        assert_eq!(cigar.pos(), 5);
    }

    #[test]
    fn test_cigar_string_leading_softclips() {
        let cigar = CigarString::from(vec![Cigar::SoftClip(10), Cigar::Match(100)]).into_view(0);
        assert_eq!(cigar.leading_softclips(), 10);
        let cigar2 = CigarString::from(vec![
            Cigar::HardClip(5),
            Cigar::SoftClip(10),
            Cigar::Match(100),
        ])
        .into_view(0);
        assert_eq!(cigar2.leading_softclips(), 10);
    }

    #[test]
    fn test_cigar_string_trailing_softclips() {
        let cigar = CigarString::from(vec![Cigar::Match(100), Cigar::SoftClip(10)]).into_view(0);
        assert_eq!(cigar.trailing_softclips(), 10);
        let cigar2 = CigarString::from(vec![
            Cigar::Match(100),
            Cigar::SoftClip(10),
            Cigar::HardClip(5),
        ])
        .into_view(0);
        assert_eq!(cigar2.trailing_softclips(), 10);
    }

    #[test]
    fn test_cigar_read_pos() {
        let vpos = 5; // variant position

        // Ignore leading HardClip
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c01: 7H                 M  M
        // qpos:                  00 01
        let c01 = CigarString::from(vec![Cigar::HardClip(7), Cigar::Match(2)]).into_view(4);
        assert_eq!(c01.read_pos(vpos, false, false).unwrap(), Some(1));

        // Skip leading SoftClip or use as pre-POS matches
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c02: 5H2S         M  M  M  M  M  M
        // qpos:  00        02 03 04 05 06 07
        // c02: 5H     S  S  M  M  M  M  M  M
        // qpos:      00 01 02 03 04 05 06 07
        let c02 = CigarString::from(vec![Cigar::SoftClip(2), Cigar::Match(6)]).into_view(2);
        assert_eq!(c02.read_pos(vpos, false, false).unwrap(), Some(5));
        assert_eq!(c02.read_pos(vpos, true, false).unwrap(), Some(5));

        // Skip leading SoftClip returning None for unmatched reference positiong or use as
        // pre-POS matches
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c03:  3S                      M  M
        // qpos: 00                     03 04
        // c03:                 S  S  S  M  M
        // qpos:               00 01 02 03 04
        let c03 = CigarString::from(vec![Cigar::SoftClip(3), Cigar::Match(6)]).into_view(6);
        assert_eq!(c03.read_pos(vpos, false, false).unwrap(), None);
        assert_eq!(c03.read_pos(vpos, true, false).unwrap(), Some(2));

        // Skip leading Insertion before variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c04:  3I                X  X  X
        // qpos: 00               03 04 05
        let c04 = CigarString::from(vec![Cigar::Ins(3), Cigar::Diff(3)]).into_view(4);
        assert_eq!(c04.read_pos(vpos, true, false).unwrap(), Some(4));

        // Matches and deletion before variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c05:        =  =  D  D  X  =  =
        // qpos:      00 01       02 03 04 05
        let c05 = CigarString::from(vec![
            Cigar::Equal(2),
            Cigar::Del(2),
            Cigar::Diff(1),
            Cigar::Equal(2),
        ])
        .into_view(0);
        assert_eq!(c05.read_pos(vpos, true, false).unwrap(), Some(3));

        // single nucleotide Deletion covering variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c06:                 =  =  D  X  X
        // qpos:               00 01    02 03
        let c06 =
            CigarString::from(vec![Cigar::Equal(2), Cigar::Del(1), Cigar::Diff(2)]).into_view(3);
        assert_eq!(c06.read_pos(vpos, false, true).unwrap(), Some(2));
        assert_eq!(c06.read_pos(vpos, false, false).unwrap(), None);

        // three nucleotide Deletion covering variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c07:              =  =  D  D  D  M  M
        // qpos:            00 01          02 03
        let c07 =
            CigarString::from(vec![Cigar::Equal(2), Cigar::Del(3), Cigar::Match(2)]).into_view(2);
        assert_eq!(c07.read_pos(vpos, false, true).unwrap(), Some(2));
        assert_eq!(c07.read_pos(vpos, false, false).unwrap(), None);

        // three nucleotide RefSkip covering variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c08:              =  X  N  N  N  M  M
        // qpos:            00 01          02 03
        let c08 = CigarString::from(vec![
            Cigar::Equal(1),
            Cigar::Diff(1),
            Cigar::RefSkip(3),
            Cigar::Match(2),
        ])
        .into_view(2);
        assert_eq!(c08.read_pos(vpos, false, true).unwrap(), None);
        assert_eq!(c08.read_pos(vpos, false, false).unwrap(), None);

        // internal hard clip before variant pos
        // ref:       00 01 02 03    04 05 06 07 08 09 10 11 12 13 14 15
        // var:                          V
        // c09: 3H           =  = 3H  =  =
        // qpos:            00 01    02 03
        let c09 = CigarString::from(vec![
            Cigar::HardClip(3),
            Cigar::Equal(2),
            Cigar::HardClip(3),
            Cigar::Equal(2),
        ])
        .into_view(2);
        assert_eq!(c09.read_pos(vpos, false, true).is_err(), true);

        // Deletion right before variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c10:           M  M  D  D  M  M
        // qpos:         00 01       02 03
        let c10 =
            CigarString::from(vec![Cigar::Match(2), Cigar::Del(2), Cigar::Match(2)]).into_view(1);
        assert_eq!(c10.read_pos(vpos, false, false).unwrap(), Some(2));

        // Insertion right before variant position
        // ref:       00 01 02 03 04    05 06 07 08 09 10 11 12 13 14 15
        // var:                          V
        // c11:                 M  M 3I  M
        // qpos:               00 01 02 05 06
        let c11 =
            CigarString::from(vec![Cigar::Match(2), Cigar::Ins(3), Cigar::Match(2)]).into_view(3);
        assert_eq!(c11.read_pos(vpos, false, false).unwrap(), Some(5));

        // Insertion right after variant position
        // ref:       00 01 02 03 04 05    06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c12:                 M  M  M 2I  =
        // qpos:               00 01 02 03 05
        let c12 =
            CigarString::from(vec![Cigar::Match(3), Cigar::Ins(2), Cigar::Equal(1)]).into_view(3);
        assert_eq!(c12.read_pos(vpos, false, false).unwrap(), Some(2));

        // Deletion right after variant position
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c13:                 M  M  M  D  =
        // qpos:               00 01 02    03
        let c13 =
            CigarString::from(vec![Cigar::Match(3), Cigar::Del(1), Cigar::Equal(1)]).into_view(3);
        assert_eq!(c13.read_pos(vpos, false, false).unwrap(), Some(2));

        // A messy and complicated example, including a Pad operation
        let vpos2 = 15;
        // ref:       00    01 02    03 04 05    06 07 08 09 10 11 12 13 14 15
        // var:                                                           V
        // c14: 5H3S   = 2P  M  X 3I  M  M  D 2I  =  =  N  N  N  M  M  M  =  =  5S2H
        // qpos:  00  03    04 05 06 09 10    11 13 14          15 16 17 18 19
        let c14 = CigarString::from(vec![
            Cigar::HardClip(5),
            Cigar::SoftClip(3),
            Cigar::Equal(1),
            Cigar::Pad(2),
            Cigar::Match(1),
            Cigar::Diff(1),
            Cigar::Ins(3),
            Cigar::Match(2),
            Cigar::Del(1),
            Cigar::Ins(2),
            Cigar::Equal(2),
            Cigar::RefSkip(3),
            Cigar::Match(3),
            Cigar::Equal(2),
            Cigar::SoftClip(5),
            Cigar::HardClip(2),
        ])
        .into_view(0);
        assert_eq!(c14.read_pos(vpos2, false, false).unwrap(), Some(19));

        // HardClip after Pad
        // ref:       00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
        // var:                       V
        // c15: 5P1H            =  =  =
        // qpos:               00 01 02
        let c15 = CigarString::from(vec![Cigar::Pad(5), Cigar::HardClip(1), Cigar::Equal(3)])
            .into_view(3);
        assert_eq!(c15.read_pos(vpos, false, false).is_err(), true);

        // only HardClip and Pad operations
        // c16: 7H5P2H
        let c16 = CigarString::from(vec![Cigar::HardClip(7), Cigar::Pad(5), Cigar::HardClip(2)])
            .into_view(3);
        assert_eq!(c16.read_pos(vpos, false, false).unwrap(), None);
    }

    #[test]
    fn test_clone() {
        let mut rec = Record::new();
        rec.set_pos(300);
        rec.set_qname(b"read1");
        let clone = rec.clone();
        assert_eq!(rec, clone);
    }

    #[test]
    fn test_flags() {
        let mut rec = Record::new();

        rec.set_paired();
        assert_eq!(rec.is_paired(), true);

        rec.set_supplementary();
        assert_eq!(rec.is_supplementary(), true);
        assert_eq!(rec.is_supplementary(), true);

        rec.unset_paired();
        assert_eq!(rec.is_paired(), false);
        assert_eq!(rec.is_supplementary(), true);

        rec.unset_supplementary();
        assert_eq!(rec.is_paired(), false);
        assert_eq!(rec.is_supplementary(), false);
    }

    #[test]
    fn test_cigar_parse() {
        let cigar = "1S20M1D2I3X1=2H";
        let parsed = CigarString::try_from(cigar).unwrap();
        assert_eq!(parsed.to_string(), cigar);
    }

    // Helper: build a Record with a given sequence (and matching dummy quals).
    fn make_record_with_seq(seq: &[u8]) -> Record {
        let mut rec = Record::new();
        let cigar = CigarString::from(vec![Cigar::Match(seq.len() as u32)]);
        let qual: Vec<u8> = vec![30u8; seq.len()];
        rec.set(b"read1", Some(&cigar), seq, &qual);
        rec
    }

    #[test]
    fn test_set_seq_roundtrip() {
        let seq = b"ACGTN";
        let mut rec = make_record_with_seq(seq);
        // Overwrite with the same content via set_seq to exercise the code path.
        rec.set_seq(seq);
        let read_back = rec.seq();
        assert_eq!(read_back.len(), seq.len());
        for i in 0..seq.len() {
            assert_eq!(read_back[i], seq[i], "mismatch at position {}", i);
        }
    }

    #[test]
    fn test_set_seq_odd_length() {
        let original = b"ACGNN";
        let replacement = b"TTTNN";
        let mut rec = make_record_with_seq(original);
        // Replace last 3 bases using a record that starts as length 3
        // We need a record of length 3 for the odd-length test.
        let mut rec3 = make_record_with_seq(b"ACG");
        rec3.set_seq(b"TGC");
        let read_back = rec3.seq();
        assert_eq!(read_back.len(), 3);
        assert_eq!(read_back[0], b'T');
        assert_eq!(read_back[1], b'G');
        assert_eq!(read_back[2], b'C');
        // The 5-base record replacement still works correctly.
        rec.set_seq(replacement);
        let read_back5 = rec.seq();
        assert_eq!(read_back5.len(), 5);
        assert_eq!(read_back5[0], b'T');
        assert_eq!(read_back5[1], b'T');
        assert_eq!(read_back5[2], b'T');
    }

    #[test]
    fn test_set_seq_all_bases() {
        let seq = b"ACGTN";
        let mut rec = make_record_with_seq(seq);
        rec.set_seq(seq);
        let read_back = rec.seq();
        assert_eq!(read_back[0], b'A');
        assert_eq!(read_back[1], b'C');
        assert_eq!(read_back[2], b'G');
        assert_eq!(read_back[3], b'T');
        assert_eq!(read_back[4], b'N');
    }

    #[test]
    #[should_panic(expected = "new_seq.len() must equal current seq.len()")]
    fn test_set_seq_wrong_length_panics() {
        let mut rec = make_record_with_seq(b"ACGT");
        // Attempting to set a sequence of different length must panic.
        rec.set_seq(b"ACGTT");
    }

    #[test]
    fn test_record_view_matches_record() {
        let seq = b"ACGTN";
        let qual: Vec<u8> = vec![40u8; seq.len()];
        let cigar = CigarString::from(vec![Cigar::Match(seq.len() as u32)]);
        let mut rec = Record::new();
        rec.set(b"testread", Some(&cigar), seq, &qual);
        rec.set_mapq(60);
        rec.set_flags(0x10); // reverse strand

        // SAFETY: rec is alive for the duration of this test; inner_ptr() returns
        // a pointer to the bam1_t embedded in rec, which remains valid as long as
        // rec is not moved or dropped.
        let view = unsafe { RecordView::from_raw(rec.inner_ptr()) };

        assert_eq!(view.qname(), rec.qname());
        assert_eq!(view.seq_len(), rec.seq_len());
        assert_eq!(view.mapq(), rec.mapq());
        assert_eq!(view.flags(), rec.flags());
        assert_eq!(view.is_reverse(), rec.is_reverse());
    }

    // Build a record with a small set of aux fields covering several types.
    fn make_record_with_aux() -> Record {
        let mut rec = make_record_with_seq(b"ACGT");
        rec.push_aux(b"XI", Aux::I32(42)).unwrap();
        rec.push_aux(b"XS", Aux::String("hello")).unwrap();
        rec.push_aux(b"XF", Aux::Float(1.5)).unwrap();
        rec
    }

    #[test]
    fn test_record_view_aux_roundtrip() {
        let rec = make_record_with_aux();
        // SAFETY: rec lives for the duration of this test.
        let view = unsafe { RecordView::from_raw(rec.inner_ptr()) };

        assert_eq!(view.aux(b"XI").unwrap(), rec.aux(b"XI").unwrap());
        assert_eq!(view.aux(b"XS").unwrap(), rec.aux(b"XS").unwrap());
        assert_eq!(view.aux(b"XF").unwrap(), rec.aux(b"XF").unwrap());
    }

    #[test]
    fn test_record_view_aux_not_found() {
        let rec = make_record_with_aux();
        // SAFETY: rec lives for the duration of this test.
        let view = unsafe { RecordView::from_raw(rec.inner_ptr()) };

        assert!(matches!(view.aux(b"ZZ"), Err(Error::AuxTagNotFound)));
    }

    #[test]
    fn test_record_view_aux_short_tag() {
        let rec = make_record_with_seq(b"ACGT");
        // SAFETY: rec lives for the duration of this test.
        let view = unsafe { RecordView::from_raw(rec.inner_ptr()) };

        assert!(matches!(view.aux(b"X"), Err(Error::AuxStringError)));
    }

    #[test]
    fn test_record_view_raw_aux_data_no_aux() {
        let rec = make_record_with_seq(b"ACGT");
        // SAFETY: rec lives for the duration of this test.
        let view = unsafe { RecordView::from_raw(rec.inner_ptr()) };

        // A record with no aux fields must return Some of an empty slice.
        assert_eq!(view.raw_aux_data(), Some(&[][..]));
    }

    #[test]
    fn test_record_view_raw_aux_data_matches_aux_iter() {
        let rec = make_record_with_aux();
        // SAFETY: rec lives for the duration of this test.
        let view = unsafe { RecordView::from_raw(rec.inner_ptr()) };

        // raw_aux_data must cover at least all bytes reported by aux_iter.
        let raw = view.raw_aux_data().expect("raw_aux_data returned None");
        assert!(!raw.is_empty(), "expected non-empty aux section");

        // Every tag yielded by aux_iter must also be found via aux().
        for item in view.aux_iter() {
            let (tag, _val) = item.unwrap();
            view.aux(tag).expect("aux() must find tag seen in aux_iter");
        }
    }

    #[test]
    fn test_record_view_aux_iter_matches_record() {
        let rec = make_record_with_aux();
        // SAFETY: rec lives for the duration of this test.
        let view = unsafe { RecordView::from_raw(rec.inner_ptr()) };

        let record_tags: Vec<Vec<u8>> = rec.aux_iter().map(|r| r.unwrap().0.to_vec()).collect();
        let view_tags: Vec<Vec<u8>> = view.aux_iter().map(|r| r.unwrap().0.to_vec()).collect();

        assert_eq!(record_tags, view_tags);
    }
}

#[cfg(test)]
mod aux_tag_search_tests {
    use super::*;
    use proptest::prelude::*;
    use std::os::raw::c_char;

    /// Helper: build a Record with a given sequence (and matching dummy quals).
    fn make_record(seq: &[u8]) -> Record {
        let mut rec = Record::new();
        let cigar = CigarString::from(vec![Cigar::Match(seq.len() as u32)]);
        let qual: Vec<u8> = vec![30u8; seq.len()];
        rec.set(b"read1", Some(&cigar), seq, &qual);
        rec
    }

    /// Call the C `bam_aux_get` on a record and parse the result via `parse_aux_field`,
    /// returning the same `Result<Aux>` that `Record::aux()` returns.
    /// This is the **oracle** implementation — whatever C returns is correct.
    fn aux_get_c<'a>(rec: &'a Record, tag: &[u8]) -> Result<Aux<'a>> {
        let aux = unsafe {
            htslib::bam_aux_get(
                &rec.inner as *const htslib::bam1_t,
                tag.as_ptr() as *const c_char,
            )
        };
        unsafe { parse_aux_field(aux).map(|(aux_field, _length)| aux_field) }
    }

    /// All valid 2-byte tag characters for BAM aux fields.
    /// Tags are [A-Za-z][A-Za-z0-9] per the SAM spec.
    fn tag_strategy() -> impl Strategy<Value = [u8; 2]> {
        let first = prop::sample::select((b'A'..=b'Z').chain(b'a'..=b'z').collect::<Vec<u8>>());
        let second = prop::sample::select(
            (b'A'..=b'Z')
                .chain(b'a'..=b'z')
                .chain(b'0'..=b'9')
                .collect::<Vec<u8>>(),
        );
        (first, second).prop_map(|(a, b)| [a, b])
    }

    /// Strategy to generate an arbitrary scalar Aux value.
    fn scalar_aux_strategy() -> impl Strategy<Value = Aux<'static>> {
        prop_oneof![
            any::<u8>().prop_map(Aux::Char),
            any::<i8>().prop_map(Aux::I8),
            any::<u8>().prop_map(Aux::U8),
            any::<i16>().prop_map(Aux::I16),
            any::<u16>().prop_map(Aux::U16),
            any::<i32>().prop_map(Aux::I32),
            any::<u32>().prop_map(Aux::U32),
            any::<f32>()
                .prop_filter("must be finite", |f| f.is_finite())
                .prop_map(Aux::Float),
            any::<f64>()
                .prop_filter("must be finite", |f| f.is_finite())
                .prop_map(Aux::Double),
        ]
    }

    /// Strategy to generate a list of (unique tag, scalar aux value) pairs.
    /// Tags are deduplicated because BAM records cannot have duplicate tags.
    fn aux_fields_strategy() -> impl Strategy<Value = Vec<([u8; 2], Aux<'static>)>> {
        proptest::collection::vec((tag_strategy(), scalar_aux_strategy()), 0..=8).prop_map(
            |fields| {
                let mut seen = std::collections::HashSet::new();
                fields
                    .into_iter()
                    .filter(|(tag, _)| seen.insert(*tag))
                    .collect()
            },
        )
    }

    /// Copy a scalar Aux value. Panics on borrowed variants (String/Array).
    fn copy_scalar_aux(aux: &Aux<'_>) -> Aux<'static> {
        match *aux {
            Aux::Char(v) => Aux::Char(v),
            Aux::I8(v) => Aux::I8(v),
            Aux::U8(v) => Aux::U8(v),
            Aux::I16(v) => Aux::I16(v),
            Aux::U16(v) => Aux::U16(v),
            Aux::I32(v) => Aux::I32(v),
            Aux::U32(v) => Aux::U32(v),
            Aux::Float(v) => Aux::Float(v),
            Aux::Double(v) => Aux::Double(v),
            _ => panic!("copy_scalar_aux called on non-scalar Aux"),
        }
    }

    /// Build a record with the given scalar aux fields pushed.
    fn build_record_with_aux(fields: &[([u8; 2], Aux<'static>)]) -> Record {
        let mut rec = make_record(b"ACGT");
        for (tag, val) in fields {
            rec.push_aux(tag, copy_scalar_aux(val)).unwrap();
        }
        rec
    }

    // ---- Differential proptests: Rust aux_tag_search vs C bam_aux_get ----

    proptest! {
        /// For every tag present in the record, the Rust implementation must
        /// return the same parsed Aux value as the C implementation.
        #[test]
        fn aux_tag_search_matches_c_for_present_tags(
            fields in aux_fields_strategy(),
        ) {
            let rec = build_record_with_aux(&fields);
            for (tag, _expected) in &fields {
                let c_result = aux_get_c(&rec, tag);
                let rs_result = rec.aux(tag);
                prop_assert_eq!(
                    c_result, rs_result,
                    "Mismatch for tag {:?}", std::str::from_utf8(tag)
                );
            }
        }

        /// For random tags NOT in the record, both implementations must return
        /// the same error (AuxTagNotFound).
        #[test]
        fn aux_tag_search_matches_c_for_absent_tags(
            fields in aux_fields_strategy(),
            query_tag in tag_strategy(),
        ) {
            let rec = build_record_with_aux(&fields);
            let tag_present = fields.iter().any(|(t, _)| *t == query_tag);
            if !tag_present {
                let c_result = aux_get_c(&rec, &query_tag);
                let rs_result = rec.aux(&query_tag);
                prop_assert_eq!(
                    c_result, rs_result,
                    "Mismatch for absent tag {:?}", std::str::from_utf8(&query_tag)
                );
            }
        }

        /// The Rust implementation must agree with C on records with no aux data.
        #[test]
        fn aux_tag_search_matches_c_on_empty_aux(
            query_tag in tag_strategy(),
        ) {
            let rec = make_record(b"ACGT");
            let c_result = aux_get_c(&rec, &query_tag);
            let rs_result = rec.aux(&query_tag);
            prop_assert_eq!(c_result, rs_result);
        }

        /// For records with multiple fields, the Rust implementation must find
        /// each field in the correct position (first, middle, last).
        #[test]
        fn aux_tag_search_matches_c_positional(
            fields in aux_fields_strategy().prop_filter(
                "need at least 3 fields",
                |f| f.len() >= 3
            ),
        ) {
            let rec = build_record_with_aux(&fields);

            // Check first, last, and a middle tag
            let first_tag = &fields[0].0;
            let last_tag = &fields[fields.len() - 1].0;
            let mid_tag = &fields[fields.len() / 2].0;

            for tag in [first_tag, mid_tag, last_tag] {
                let c_result = aux_get_c(&rec, tag);
                let rs_result = rec.aux(tag);
                prop_assert_eq!(
                    c_result, rs_result,
                    "Positional mismatch for tag {:?}", std::str::from_utf8(tag)
                );
            }
        }
    }

    // ---- Deterministic tests for string and array types ----
    // (These can't easily be generated with proptest due to lifetime constraints
    //  on Aux::String and Aux::Array* — so we test them deterministically.)

    #[test]
    fn aux_tag_search_matches_c_for_string_fields() {
        let mut rec = make_record(b"ACGT");
        rec.push_aux(b"XI", Aux::I32(42)).unwrap();
        rec.push_aux(b"XS", Aux::String("hello world")).unwrap();
        rec.push_aux(b"XF", Aux::Float(1.5)).unwrap();

        // String tag: both impls must agree
        assert_eq!(aux_get_c(&rec, b"XS"), rec.aux(b"XS"));
        // Tags after the string: both impls must agree
        assert_eq!(aux_get_c(&rec, b"XF"), rec.aux(b"XF"));
        // Tag before the string: both impls must agree
        assert_eq!(aux_get_c(&rec, b"XI"), rec.aux(b"XI"));
        // Absent tag
        assert_eq!(aux_get_c(&rec, b"ZZ"), rec.aux(b"ZZ"));
    }

    #[test]
    fn aux_tag_search_matches_c_for_array_fields() {
        let mut rec = make_record(b"ACGT");
        let arr_i32: Vec<i32> = vec![1, 2, 3, 4, 5];
        let arr_f32: Vec<f32> = vec![1.0, 2.0, 3.0];
        rec.push_aux(b"XI", Aux::I32(42)).unwrap();
        rec.push_aux(b"XE", Aux::ArrayI32((&arr_i32).into()))
            .unwrap();
        rec.push_aux(b"XG", Aux::ArrayFloat((&arr_f32).into()))
            .unwrap();
        rec.push_aux(b"XS", Aux::String("after_array")).unwrap();

        for tag in [b"XI", b"XE", b"XG", b"XS"] {
            assert_eq!(
                aux_get_c(&rec, tag.as_slice()),
                rec.aux(tag.as_slice()),
                "Mismatch for tag {:?}",
                std::str::from_utf8(tag.as_slice())
            );
        }
        assert_eq!(aux_get_c(&rec, b"ZZ"), rec.aux(b"ZZ"));
    }

    #[test]
    fn aux_tag_search_matches_c_for_all_array_subtypes() {
        let mut rec = make_record(b"ACGT");
        let arr_i8: Vec<i8> = vec![-1, 0, 1];
        let arr_u8: Vec<u8> = vec![0, 128, 255];
        let arr_i16: Vec<i16> = vec![-1000, 0, 1000];
        let arr_u16: Vec<u16> = vec![0, 30000, 65535];
        let arr_i32: Vec<i32> = vec![-100000, 0, 100000];
        let arr_u32: Vec<u32> = vec![0, 2000000000, 4000000000];
        let arr_f32: Vec<f32> = vec![-1.5, 0.0, 1.5];
        rec.push_aux(b"Xa", Aux::ArrayI8((&arr_i8).into())).unwrap();
        rec.push_aux(b"Xb", Aux::ArrayU8((&arr_u8).into())).unwrap();
        rec.push_aux(b"Xc", Aux::ArrayI16((&arr_i16).into()))
            .unwrap();
        rec.push_aux(b"Xd", Aux::ArrayU16((&arr_u16).into()))
            .unwrap();
        rec.push_aux(b"Xe", Aux::ArrayI32((&arr_i32).into()))
            .unwrap();
        rec.push_aux(b"Xf", Aux::ArrayU32((&arr_u32).into()))
            .unwrap();
        rec.push_aux(b"Xg", Aux::ArrayFloat((&arr_f32).into()))
            .unwrap();
        // A scalar after all arrays — tests correct size skipping
        rec.push_aux(b"ZZ", Aux::I32(99)).unwrap();

        for tag in [b"Xa", b"Xb", b"Xc", b"Xd", b"Xe", b"Xf", b"Xg", b"ZZ"] {
            assert_eq!(
                aux_get_c(&rec, tag.as_slice()),
                rec.aux(tag.as_slice()),
                "Mismatch for tag {:?}",
                std::str::from_utf8(tag.as_slice())
            );
        }
    }

    #[test]
    fn aux_tag_search_matches_c_after_remove() {
        let mut rec = make_record(b"ACGT");
        rec.push_aux(b"XI", Aux::I32(42)).unwrap();
        rec.push_aux(b"XS", Aux::String("hello")).unwrap();
        rec.push_aux(b"XF", Aux::Float(1.5)).unwrap();

        rec.remove_aux(b"XS").unwrap();

        for tag in [b"XI", b"XS", b"XF"] {
            assert_eq!(
                aux_get_c(&rec, tag.as_slice()),
                rec.aux(tag.as_slice()),
                "Mismatch for tag {:?} after remove",
                std::str::from_utf8(tag.as_slice())
            );
        }
    }
}

#[cfg(test)]
mod alignment_cigar_tests {
    use super::*;
    use crate::bam::{Read, Reader};
    use bio_types::alignment::AlignmentOperation::{Del, Ins, Match, Subst, Xclip, Yclip};
    use bio_types::alignment::{Alignment, AlignmentMode};

    #[test]
    fn test_cigar() {
        let alignment = Alignment {
            score: 5,
            xstart: 3,
            ystart: 0,
            xend: 9,
            yend: 10,
            ylen: 10,
            xlen: 10,
            operations: vec![Match, Match, Match, Subst, Ins, Ins, Del, Del],
            mode: AlignmentMode::Semiglobal,
        };
        assert_eq!(alignment.cigar(false), "3S3=1X2I2D1S");
        assert_eq!(
            CigarString::from_alignment(&alignment, false),
            CigarString::from([
                Cigar::SoftClip(3),
                Cigar::Equal(3),
                Cigar::Diff(1),
                Cigar::Ins(2),
                Cigar::Del(2),
                Cigar::SoftClip(1),
            ])
        );

        let alignment = Alignment {
            score: 5,
            xstart: 0,
            ystart: 5,
            xend: 4,
            yend: 10,
            ylen: 10,
            xlen: 5,
            operations: vec![Yclip(5), Match, Subst, Subst, Ins, Del, Del, Xclip(1)],
            mode: AlignmentMode::Custom,
        };
        assert_eq!(alignment.cigar(false), "1=2X1I2D1S");
        assert_eq!(alignment.cigar(true), "1=2X1I2D1H");
        assert_eq!(
            CigarString::from_alignment(&alignment, false),
            CigarString::from([
                Cigar::Equal(1),
                Cigar::Diff(2),
                Cigar::Ins(1),
                Cigar::Del(2),
                Cigar::SoftClip(1),
            ])
        );
        assert_eq!(
            CigarString::from_alignment(&alignment, true),
            CigarString::from([
                Cigar::Equal(1),
                Cigar::Diff(2),
                Cigar::Ins(1),
                Cigar::Del(2),
                Cigar::HardClip(1),
            ])
        );

        let alignment = Alignment {
            score: 5,
            xstart: 0,
            ystart: 5,
            xend: 3,
            yend: 8,
            ylen: 10,
            xlen: 3,
            operations: vec![Yclip(5), Subst, Match, Subst, Yclip(2)],
            mode: AlignmentMode::Custom,
        };
        assert_eq!(alignment.cigar(false), "1X1=1X");
        assert_eq!(
            CigarString::from_alignment(&alignment, false),
            CigarString::from([Cigar::Diff(1), Cigar::Equal(1), Cigar::Diff(1)])
        );

        let alignment = Alignment {
            score: 5,
            xstart: 0,
            ystart: 5,
            xend: 3,
            yend: 8,
            ylen: 10,
            xlen: 3,
            operations: vec![Subst, Match, Subst],
            mode: AlignmentMode::Semiglobal,
        };
        assert_eq!(alignment.cigar(false), "1X1=1X");
        assert_eq!(
            CigarString::from_alignment(&alignment, false),
            CigarString::from([Cigar::Diff(1), Cigar::Equal(1), Cigar::Diff(1)])
        );
    }

    #[test]
    fn test_read_orientation_f1r2() {
        let mut bam = Reader::from_path("test/test_paired.sam").unwrap();

        for res in bam.records() {
            let record = res.unwrap();
            assert_eq!(
                record.read_pair_orientation(),
                SequenceReadPairOrientation::F1R2
            );
        }
    }

    #[test]
    fn test_read_orientation_f2r1() {
        let mut bam = Reader::from_path("test/test_nonstandard_orientation.sam").unwrap();

        for res in bam.records() {
            let record = res.unwrap();
            assert_eq!(
                record.read_pair_orientation(),
                SequenceReadPairOrientation::F2R1
            );
        }
    }

    #[test]
    fn test_read_orientation_supplementary() {
        let mut bam = Reader::from_path("test/test_orientation_supplementary.sam").unwrap();

        for res in bam.records() {
            let record = res.unwrap();
            assert_eq!(
                record.read_pair_orientation(),
                SequenceReadPairOrientation::F2R1
            );
        }
    }

    #[test]
    pub fn test_cigar_parsing_non_ascii_error() {
        let cigar_str = "43ጷ";
        let expected_error = Err(Error::ParseCigar {
                msg: "CIGAR string contained non-ASCII characters, which are not valid. Valid are [0-9MIDNSHP=X].".to_owned(),
            });

        let result = CigarString::try_from(cigar_str);
        assert_eq!(expected_error, result);
    }

    #[test]
    pub fn test_cigar_parsing() {
        // parsing test cases
        let cigar_strs = [
            "1H10M4D100I300N1102=10P25X11S", // test every cigar opt
            "100M",                          // test a single op
            "",                              // test empty input
            "1H1=1H",                        // test simple hardclip
            "1S1=1S",                        // test simple softclip
            "11H11S11=11S11H",               // test complex softclip
            "10H",
            "10S",
        ];
        // expected results
        let cigars = [
            CigarString::from(vec![
                Cigar::HardClip(1),
                Cigar::Match(10),
                Cigar::Del(4),
                Cigar::Ins(100),
                Cigar::RefSkip(300),
                Cigar::Equal(1102),
                Cigar::Pad(10),
                Cigar::Diff(25),
                Cigar::SoftClip(11),
            ]),
            CigarString::from(vec![Cigar::Match(100)]),
            CigarString::from(vec![]),
            CigarString::from(vec![
                Cigar::HardClip(1),
                Cigar::Equal(1),
                Cigar::HardClip(1),
            ]),
            CigarString::from(vec![
                Cigar::SoftClip(1),
                Cigar::Equal(1),
                Cigar::SoftClip(1),
            ]),
            CigarString::from(vec![
                Cigar::HardClip(11),
                Cigar::SoftClip(11),
                Cigar::Equal(11),
                Cigar::SoftClip(11),
                Cigar::HardClip(11),
            ]),
            CigarString::from(vec![Cigar::HardClip(10)]),
            CigarString::from(vec![Cigar::SoftClip(10)]),
        ];
        // compare
        for (&cigar_str, truth) in cigar_strs.iter().zip(cigars.iter()) {
            let cigar_parse = CigarString::try_from(cigar_str)
                .unwrap_or_else(|_| panic!("Unable to parse cigar: {}", cigar_str));
            assert_eq!(&cigar_parse, truth);
        }
    }
}

#[cfg(test)]
mod basemod_tests {
    use crate::bam::{Read, Reader};

    #[test]
    pub fn test_count_recorded() {
        let mut bam = Reader::from_path("test/base_mods/MM-double.sam").unwrap();

        for r in bam.records() {
            let record = r.unwrap();
            if let Ok(mods) = record.basemods_iter() {
                let n = mods.recorded().len();
                assert_eq!(n, 3);
            };
        }
    }

    #[test]
    pub fn test_query_type() {
        let mut bam = Reader::from_path("test/base_mods/MM-orient.sam").unwrap();

        let mut n_fwd = 0;
        let mut n_rev = 0;

        for r in bam.records() {
            let record = r.unwrap();
            if let Ok(mods) = record.basemods_iter() {
                for mod_code in mods.recorded() {
                    if let Ok(mod_metadata) = mods.query_type(*mod_code) {
                        if mod_metadata.strand == 0 {
                            n_fwd += 1;
                        }
                        if mod_metadata.strand == 1 {
                            n_rev += 1;
                        }
                    }
                }
            };
        }
        assert_eq!(n_fwd, 2);
        assert_eq!(n_rev, 2);
    }

    #[test]
    pub fn test_mod_iter() {
        let mut bam = Reader::from_path("test/base_mods/MM-double.sam").unwrap();
        let expected_positions = [1, 7, 12, 13, 13, 22, 30, 31];
        let mut i = 0;

        for r in bam.records() {
            let record = r.unwrap();
            for res in record.basemods_iter().unwrap().flatten() {
                let (position, _m) = res;
                assert_eq!(position, expected_positions[i]);
                i += 1;
            }
        }
    }

    #[test]
    pub fn test_position_iter() {
        let mut bam = Reader::from_path("test/base_mods/MM-double.sam").unwrap();
        let expected_positions = [1, 7, 12, 13, 22, 30, 31];
        let expected_counts = [1, 1, 1, 2, 1, 1, 1];
        let mut i = 0;

        for r in bam.records() {
            let record = r.unwrap();
            for res in record.basemods_position_iter().unwrap().flatten() {
                let (position, elements) = res;
                assert_eq!(position, expected_positions[i]);
                assert_eq!(elements.len(), expected_counts[i]);
                i += 1;
            }
        }
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    /// The 16 valid BAM-encodable bases (IUPAC + '=').
    const VALID_BASES: &[u8] = b"=ACMGRSVTWYHKDBN";

    /// Reference implementation: decode one base at a time using the original
    /// per-base nibble-extraction logic.
    fn as_bytes_naive(seq: &Seq<'_>) -> Vec<u8> {
        (0..seq.len())
            .map(|i| *decode_base_unchecked(encoded_base(seq.encoded, i)))
            .collect()
    }

    /// Pack ASCII bases into BAM 4-bit encoding (same logic as Record::set).
    fn encode_bases(bases: &[u8]) -> Vec<u8> {
        let mut encoded = vec![0u8; bases.len().div_ceil(2)];
        for (j, chunk) in bases.chunks(2).enumerate() {
            encoded[j] = ENCODE_BASE[chunk[0] as usize] << 4
                | if chunk.len() == 2 {
                    ENCODE_BASE[chunk[1] as usize]
                } else {
                    0
                };
        }
        encoded
    }

    // -- Strategies --

    /// Random packed bytes with a random base count. Exercises the full
    /// nibble space (0..16) including ambiguity codes that real BAM data
    /// might not contain.
    fn raw_seq_strategy() -> impl Strategy<Value = (Vec<u8>, usize)> {
        (0usize..=512).prop_flat_map(|len| {
            let n_bytes = len.div_ceil(2);
            (proptest::collection::vec(any::<u8>(), n_bytes), Just(len))
        })
    }

    /// Random sequence of valid ASCII bases (like "ACGTNN=").
    /// More realistic than raw bytes — these are the strings a caller
    /// would pass to Record::set.
    fn ascii_bases_strategy() -> impl Strategy<Value = Vec<u8>> {
        proptest::collection::vec(proptest::sample::select(VALID_BASES), 0..=512)
    }

    proptest! {
        /// The optimized `as_bytes()` must produce identical output to
        /// the naive per-base decoder for arbitrary packed data.
        #[test]
        fn as_bytes_matches_naive((encoded, len) in raw_seq_strategy()) {
            let seq = Seq { encoded: &encoded, len };
            prop_assert_eq!(seq.as_bytes(), as_bytes_naive(&seq));
        }

        /// `Seq::get(i)` must agree with `seq[i]` for valid indices and
        /// return `None` for out-of-bounds.
        #[test]
        fn get_matches_index((encoded, len) in raw_seq_strategy()) {
            let seq = Seq { encoded: &encoded, len };
            for i in 0..len {
                prop_assert_eq!(seq.get(i), Some(seq[i]));
            }
            prop_assert_eq!(seq.get(len), None);
            prop_assert_eq!(seq.get(usize::MAX), None);
        }

        /// Every byte produced by `as_bytes()` must be one of the 16
        /// valid BAM base characters.
        #[test]
        fn as_bytes_only_valid_bases((encoded, len) in raw_seq_strategy()) {
            let seq = Seq { encoded: &encoded, len };
            for (i, &b) in seq.as_bytes().iter().enumerate() {
                prop_assert!(
                    VALID_BASES.contains(&b),
                    "base {} at position {} is not a valid BAM base", b as char, i
                );
            }
        }

        /// Roundtrip: ASCII bases → encode → decode → must recover the
        /// original bases exactly. This is the fundamental correctness
        /// invariant of the codec.
        #[test]
        fn roundtrip_encode_decode(bases in ascii_bases_strategy()) {
            let encoded = encode_bases(&bases);
            let seq = Seq { encoded: &encoded, len: bases.len() };
            prop_assert_eq!(seq.as_bytes(), bases);
        }

        /// Roundtrip in reverse: packed bytes → decode → re-encode →
        /// must recover the same packed bytes. The trailing nibble of an
        /// odd-length sequence is zeroed by encode, so we mask it before
        /// comparing.
        #[test]
        fn roundtrip_decode_encode((encoded, len) in raw_seq_strategy()) {
            let seq = Seq { encoded: &encoded, len };
            let decoded = seq.as_bytes();
            let re_encoded = encode_bases(&decoded);

            // Mask the unused low nibble of the last byte for odd lengths
            let mut expected = encoded.clone();
            if len % 2 == 1 {
                if let Some(last) = expected.last_mut() {
                    *last &= 0xF0;
                }
            }
            prop_assert_eq!(re_encoded, expected);
        }

        /// Changing one nibble must affect exactly one decoded base and
        /// leave all others unchanged.
        #[test]
        fn nibble_independence(
            (encoded, len) in raw_seq_strategy().prop_filter(
                "need at least 1 base", |(_, l)| *l >= 1
            ),
            target_base in 0usize..512,
            new_nibble in 0u8..16,
        ) {
            let target_base = target_base % len;
            let byte_idx = target_base / 2;
            let is_high = target_base % 2 == 0;

            let seq_before = Seq { encoded: &encoded, len };
            let decoded_before = seq_before.as_bytes();

            let mut modified = encoded.clone();
            if is_high {
                modified[byte_idx] = (new_nibble << 4) | (modified[byte_idx] & 0x0F);
            } else {
                modified[byte_idx] = (modified[byte_idx] & 0xF0) | new_nibble;
            }

            let seq_after = Seq { encoded: &modified, len };
            let decoded_after = seq_after.as_bytes();

            for i in 0..len {
                if i == target_base {
                    // This base may have changed (or stayed same if nibble matches)
                    let expected_nibble = if is_high {
                        encoded[byte_idx] >> 4
                    } else {
                        encoded[byte_idx] & 0x0F
                    };
                    if new_nibble == expected_nibble {
                        prop_assert_eq!(decoded_after[i], decoded_before[i]);
                    }
                } else {
                    prop_assert_eq!(
                        decoded_after[i], decoded_before[i],
                        "base at position {} changed when only position {} was modified",
                        i, target_base
                    );
                }
            }
        }

        /// DECODE_PAIR table must be consistent with DECODE_BASE for
        /// every possible byte value.
        #[test]
        fn decode_pair_consistent_with_decode_base(byte in any::<u8>()) {
            let pair = DECODE_PAIR[byte as usize];
            prop_assert_eq!(pair[0], DECODE_BASE[(byte >> 4) as usize]);
            prop_assert_eq!(pair[1], DECODE_BASE[(byte & 0xf) as usize]);
        }

        /// The SIMD/optimized path must produce identical output to the
        /// scalar path for all inputs. This catches SIMD-specific bugs
        /// (wrong interleave order, nibble extraction errors, etc.).
        #[test]
        fn simd_matches_scalar((encoded, len) in raw_seq_strategy()) {
            let scalar = decode_seq_scalar(&encoded, len);
            let seq = Seq { encoded: &encoded, len };
            prop_assert_eq!(seq.as_bytes(), scalar);
        }

        /// CigarString construction methods must all produce equivalent
        /// results: From<Vec>, FromIterator, and collect().
        #[test]
        fn cigar_string_construction_equivalence(
            ops in proptest::collection::vec(
                prop_oneof![
                    (1u32..1000).prop_map(Cigar::Match),
                    (1u32..1000).prop_map(Cigar::Ins),
                    (1u32..1000).prop_map(Cigar::Del),
                    (1u32..1000).prop_map(Cigar::Equal),
                    (1u32..1000).prop_map(Cigar::Diff),
                    (1u32..1000).prop_map(Cigar::SoftClip),
                ],
                0..=20,
            )
        ) {
            let from_vec = CigarString::from(ops.clone());
            let from_iter: CigarString = ops.iter().copied().collect();

            prop_assert_eq!(&from_vec, &from_iter);
            prop_assert_eq!(from_vec.len(), ops.len());

            // Deref, Index, and iter() must agree
            for (i, expected) in ops.iter().enumerate() {
                prop_assert_eq!(&from_vec[i], expected);
            }
            let via_iter: Vec<&Cigar> = from_vec.iter().collect();
            let via_into_iter: Vec<&Cigar> = (&from_vec).into_iter().collect();
            prop_assert_eq!(via_iter, via_into_iter);
        }

        /// Display → TryFrom roundtrip for SAM-valid CIGAR strings
        /// (no clips, since clip placement rules make arbitrary
        /// op sequences invalid for parsing).
        #[test]
        fn cigar_display_roundtrip(
            ops in proptest::collection::vec(
                prop_oneof![
                    (1u32..1000).prop_map(Cigar::Match),
                    (1u32..1000).prop_map(Cigar::Ins),
                    (1u32..1000).prop_map(Cigar::Del),
                    (1u32..1000).prop_map(Cigar::Equal),
                    (1u32..1000).prop_map(Cigar::Diff),
                    (1u32..1000).prop_map(Cigar::RefSkip),
                    (1u32..1000).prop_map(Cigar::Pad),
                ],
                1..=20,
            )
        ) {
            let cigar = CigarString::from(ops);
            let text = format!("{}", cigar);
            let parsed = CigarString::try_from(text.as_str()).unwrap();
            prop_assert_eq!(cigar, parsed);
        }
    }

    /// Deterministic edge-case tests for lengths around SIMD boundaries.
    /// The SIMD loop processes 16 packed bytes (32 bases) per iteration,
    /// so we test lengths right around that stride.
    #[test]
    fn seq_decode_simd_boundary_lengths() {
        // Lengths that exercise: no SIMD (0..31), exact boundary (32, 64),
        // one past (33, 65), one before (31, 63), and odd variants.
        for len in [
            0usize, 1, 2, 15, 16, 30, 31, 32, 33, 34, 63, 64, 65, 127, 128, 129,
        ] {
            let n_bytes = len.div_ceil(2);
            // Use a recognizable pattern: byte i = i as u8
            let encoded: Vec<u8> = (0..n_bytes).map(|i| i as u8).collect();
            let seq = Seq {
                encoded: &encoded,
                len,
            };
            let result = seq.as_bytes();
            let scalar = decode_seq_scalar(&encoded, len);

            assert_eq!(result.len(), len, "wrong length for len={len}");
            assert_eq!(result, scalar, "SIMD/scalar mismatch for len={len}");

            // Verify each base individually
            for (i, &expected) in result.iter().enumerate() {
                assert_eq!(
                    seq.get(i),
                    Some(expected),
                    "get({i}) mismatch for len={len}"
                );
            }
            assert_eq!(seq.get(len), None);
        }
    }

    /// Test with uniform byte patterns to catch nibble-swap bugs.
    #[test]
    fn seq_decode_uniform_patterns() {
        let len = 128;
        let n_bytes = len / 2;

        for byte_val in [0x00, 0x11, 0x24, 0x42, 0x88, 0xFF] {
            let encoded = vec![byte_val; n_bytes];
            let seq = Seq {
                encoded: &encoded,
                len,
            };
            let result = seq.as_bytes();
            let scalar = decode_seq_scalar(&encoded, len);
            assert_eq!(result, scalar, "mismatch for uniform byte 0x{byte_val:02X}");

            // For uniform bytes, even positions should all be the same
            // and odd positions should all be the same.
            let hi_base = DECODE_BASE[(byte_val >> 4) as usize];
            let lo_base = DECODE_BASE[(byte_val & 0xF) as usize];
            for i in (0..len).step_by(2) {
                assert_eq!(result[i], hi_base);
                assert_eq!(result[i + 1], lo_base);
            }
        }

        // Asymmetric: 0x12 vs 0x21 — catches interleave-order bugs
        let enc_12 = vec![0x12u8; n_bytes];
        let enc_21 = vec![0x21u8; n_bytes];
        let r12 = Seq {
            encoded: &enc_12,
            len,
        }
        .as_bytes();
        let r21 = Seq {
            encoded: &enc_21,
            len,
        }
        .as_bytes();
        // High/low nibbles are swapped, so decoded bases should swap
        for i in 0..len / 2 {
            assert_eq!(r12[2 * i], r21[2 * i + 1]);
            assert_eq!(r12[2 * i + 1], r21[2 * i]);
        }
    }
}
