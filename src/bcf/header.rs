// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.
//! Module for working with VCF or BCF headers.
//!
//! # Examples
//! From the header of a VCF file we can
//!   - Output sample count of a VCF file
//!   - Output sample names of a VCF file
//!   - Output sample index given a sample name of a VCF file.
//! ```
//! use crate::rust_htslib::bcf::{Reader, Read};
//! use cstr8::cstr8;
//! use std::io::Read as IoRead;
//!
//! let path = &"test/test_string.vcf";
//! let mut bcf = Reader::from_path(path).expect("Error opening file.");
//! let header = bcf.header();
//! assert_eq!(header.sample_count(), 2);  // Sample count
//! let mut s = String::new();
//! for (i, mut x) in header.samples().into_iter().enumerate() {
//!     x.read_to_string(&mut s);  // Read sample name in to `s`
//!     println!("{}", s);  // output sample name
//! }
//! assert_eq!(header.sample_id(b"one").unwrap(), 0);  // Sample index wrapped in Option<usize>
//! assert_eq!(header.sample_id(b"two").unwrap(), 1);  // Sample index wrapped in Option<usize>
//! assert!(header.sample_id(b"non existent sample").is_none());  // Return none if not found
//!
//! assert_eq!(header.contig_count(), 1); // Number of contig in header.
//! // obtain the data type of an INFO field
//! let (tag_type, tag_length) = header.info_type(cstr8!("S1")).unwrap();
//! let (fmt_type, fmt_length) = header.format_type(cstr8!("GT")).unwrap();
//! ```

use std::collections::HashMap;
use std::ffi;
use std::os::raw::c_char;
use std::slice;
use std::sync::Arc;

use crate::htslib;

use cstr8::{CStr8, CompactCStr8};
use linear_map::LinearMap;

use crate::bcf::BcfError as Error;

type Result<T> = std::result::Result<T, Error>;

pub type SampleSubset = Vec<i32>;

custom_derive! {
    /// A newtype for IDs from BCF headers.
    #[derive(
        NewtypeFrom,
        NewtypeDeref,
        PartialEq,
        PartialOrd,
        Eq,
        Ord,
        Copy,
        Clone,
        Debug
    )]
    pub struct Id(pub u32);
}

/// A BCF header.
#[derive(Debug)]
pub struct Header {
    pub(crate) inner: *mut htslib::bcf_hdr_t,
    pub subset: Option<SampleSubset>,
}

// SAFETY: Header owns its inner pointer exclusively; no shared mutable state.
unsafe impl Send for Header {}
// SAFETY: Header owns its inner pointer exclusively; no shared mutable state.
unsafe impl Sync for Header {}

impl Default for Header {
    fn default() -> Self {
        Self::new()
    }
}

impl Header {
    /// Create a new (empty) `Header`.
    pub fn new() -> Self {
        let c_str = ffi::CString::new(&b"w"[..]).unwrap();
        Header {
            // SAFETY: c_str is a valid NUL-terminated string; bcf_hdr_init returns a valid pointer.
            inner: unsafe { htslib::bcf_hdr_init(c_str.as_ptr()) },
            subset: None,
        }
    }

    /// Get a pointer to the raw header.
    ///
    /// # Safety
    /// The caller must ensure that the pointer is not used after this `Header`
    /// is dropped
    pub unsafe fn inner_ptr(&self) -> *mut htslib::bcf_hdr_t {
        self.inner
    }

    /// Create a new `Header` using the given `HeaderView` as the template.
    ///
    /// After construction, you can modify the header independently from the template `header`.
    ///
    /// # Arguments
    ///
    /// - `header` - The `HeaderView` to use as the template.
    pub fn from_template(header: &HeaderView) -> Self {
        Header {
            // SAFETY: header.inner is a valid bcf_hdr_t pointer; bcf_hdr_dup returns a new copy.
            inner: unsafe { htslib::bcf_hdr_dup(header.inner) },
            subset: None,
        }
    }

    /// Create a new `Header` using the given `HeaderView` as as template, but subsetting to the
    /// given `samples`.
    ///
    /// # Arguments
    ///
    /// - `header` - The `HeaderView` to use for the template.
    /// - `samples` - A slice of byte-encoded (`[u8]`) sample names.
    pub fn from_template_subset(header: &HeaderView, samples: &[&[u8]]) -> Result<Self> {
        let mut imap = vec![0; samples.len()];
        let names: Vec<_> = samples
            .iter()
            .map(|&s| ffi::CString::new(s).unwrap())
            .collect();
        let name_pointers: Vec<_> = names.iter().map(|s| s.as_ptr() as *mut i8).collect();
        #[allow(clippy::unnecessary_cast)]
        let name_pointers_ptr = name_pointers.as_ptr() as *const *mut c_char;
        // SAFETY: header.inner is valid; name_pointers and imap are correctly sized.
        let inner = unsafe {
            htslib::bcf_hdr_subset(
                header.inner,
                samples.len() as i32,
                name_pointers_ptr,
                imap.as_mut_ptr(),
            )
        };
        if inner.is_null() {
            Err(Error::DuplicateSampleNames)
        } else {
            Ok(Header {
                inner,
                subset: Some(imap),
            })
        }
    }

    /// Add a `sample` to the header.
    ///
    /// # Arguments
    ///
    /// - `sample` - Name of the sample to add (to the end of the sample list).
    pub fn push_sample(&mut self, sample: &[u8]) -> &mut Self {
        let c_str = ffi::CString::new(sample).unwrap();
        // SAFETY: self.inner is non-null (from constructor); c_str is a valid CString.
        unsafe { htslib::bcf_hdr_add_sample(self.inner, c_str.as_ptr()) };
        self
    }

    /// Add a record to the header.
    ///
    /// # Arguments
    ///
    /// - `record` - String representation of the header line
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// header.push_record(format!("##contig=<ID={},length={}>", "chrX", 155270560).as_bytes());
    /// ```
    pub fn push_record(&mut self, record: &[u8]) -> &mut Self {
        let c_str = ffi::CString::new(record).unwrap();
        // SAFETY: self.inner is non-null (from constructor); c_str is a valid CString.
        unsafe { htslib::bcf_hdr_append(self.inner, c_str.as_ptr()) };
        self
    }

    /// Remove a `FILTER` entry from the header.
    ///
    /// # Arguments
    ///
    /// - `tag` - Name of the `FLT` tag to remove.
    pub fn remove_filter(&mut self, tag: &[u8]) -> &mut Self {
        self.remove_impl(tag, htslib::BCF_HL_FLT)
    }

    /// Remove an `INFO` entry from the header.
    ///
    /// # Arguments
    ///
    /// - `tag` - Name of the `INFO` tag to remove.
    pub fn remove_info(&mut self, tag: &[u8]) -> &mut Self {
        self.remove_impl(tag, htslib::BCF_HL_INFO)
    }

    /// Remove a `FORMAT` entry from the header.
    ///
    /// # Arguments
    ///
    /// - `tag` - Name of the `FORMAT` tag to remove.
    pub fn remove_format(&mut self, tag: &[u8]) -> &mut Self {
        self.remove_impl(tag, htslib::BCF_HL_FMT)
    }

    /// Remove a contig entry from the header.
    ///
    /// # Arguments
    ///
    /// - `tag` - Name of the `FORMAT` tag to remove.
    pub fn remove_contig(&mut self, tag: &[u8]) -> &mut Self {
        self.remove_impl(tag, htslib::BCF_HL_CTG)
    }

    /// Remove a structured entry from the header.
    ///
    /// # Arguments
    ///
    /// - `tag` - Name of the structured tag to remove.
    pub fn remove_structured(&mut self, tag: &[u8]) -> &mut Self {
        self.remove_impl(tag, htslib::BCF_HL_STR)
    }

    /// Remove a generic entry from the header.
    ///
    /// # Arguments
    ///
    /// - `tag` - Name of the generic tag to remove.
    pub fn remove_generic(&mut self, tag: &[u8]) -> &mut Self {
        self.remove_impl(tag, htslib::BCF_HL_GEN)
    }

    /// Implementation of removing header tags.
    fn remove_impl(&mut self, tag: &[u8], type_: u32) -> &mut Self {
        // SAFETY: self.inner is non-null (from constructor); c_str is a valid CString.
        unsafe {
            let v = tag.to_vec();
            let c_str = ffi::CString::new(v).unwrap();
            htslib::bcf_hdr_remove(self.inner, type_ as i32, c_str.as_ptr());
        }
        self
    }
}

impl Drop for Header {
    fn drop(&mut self) {
        // SAFETY: self.inner was allocated by bcf_hdr_init or bcf_hdr_dup; bcf_hdr_destroy is symmetric.
        unsafe { htslib::bcf_hdr_destroy(self.inner) };
    }
}

/// A header record.
#[derive(Debug)]
pub enum HeaderRecord {
    /// A `FILTER` header record.
    Filter {
        key: String,
        values: LinearMap<String, String>,
    },
    /// An `INFO` header record.
    Info {
        key: String,
        values: LinearMap<String, String>,
    },
    /// A `FORMAT` header record.
    Format {
        key: String,
        values: LinearMap<String, String>,
    },
    /// A `contig` header record.
    Contig {
        key: String,
        values: LinearMap<String, String>,
    },
    /// A structured header record.
    Structured {
        key: String,
        values: LinearMap<String, String>,
    },
    /// A generic, unstructured header record.
    Generic { key: String, value: String },
}

/// Number of dictionary types in BCF headers (BCF_DT_ID, BCF_DT_CTG, BCF_DT_SAMPLE).
const BCF_DICT_COUNT: usize = 3;

#[derive(Debug)]
pub struct HeaderView {
    pub(crate) inner: *mut htslib::bcf_hdr_t,
    /// Pre-built name→id caches for O(1) lookups, one per BCF_DT_* dictionary.
    /// Replaces the opaque C khash tables that are inaccessible from Rust.
    id_cache: [HashMap<CompactCStr8, i32>; BCF_DICT_COUNT],
}

// SAFETY: HeaderView owns its inner pointer exclusively; no shared mutable state.
unsafe impl Send for HeaderView {}
// SAFETY: HeaderView owns its inner pointer exclusively; no shared mutable state.
unsafe impl Sync for HeaderView {}

/// Build the name→id HashMap caches from the C header's `id[]` arrays.
///
/// # Safety
/// `inner` must be a valid, non-null pointer to an initialized `bcf_hdr_t`.
unsafe fn build_id_caches(
    inner: *const htslib::bcf_hdr_t,
) -> [HashMap<CompactCStr8, i32>; BCF_DICT_COUNT] {
    let hdr = &*inner;
    std::array::from_fn(|which| {
        let n = hdr.n[which] as usize;
        if n == 0 || hdr.id[which].is_null() {
            return HashMap::new();
        }
        let entries = slice::from_raw_parts(hdr.id[which], n);
        let mut map = HashMap::with_capacity(n);
        for entry in entries {
            if entry.key.is_null() {
                continue;
            }
            let key = match CompactCStr8::from_ptr(entry.key as *const u8) {
                Ok(k) => k,
                Err(_) => continue,
            };
            let id = (*entry.val).id;
            map.insert(key, id);
        }
        map
    })
}

impl HeaderView {
    /// Create a view from a raw pointer to a header.
    ///
    /// # Safety
    /// The caller must ensure that the header is initialized.
    pub unsafe fn from_ptr(inner: *mut htslib::bcf_hdr_t) -> Self {
        let id_cache = build_id_caches(inner);
        HeaderView { inner, id_cache }
    }

    /// Look up a name in one of the three BCF dictionaries.
    /// Returns the numeric ID, or `None` if the name is not found.
    fn dict_lookup(&self, which: usize, name: &CStr8) -> Option<i32> {
        self.id_cache[which].get(name).copied()
    }

    /// Get a pointer to the underlying raw header.
    ///
    /// # Safety
    /// The caller must ensure that the pointer is not used after this
    /// `HeaderView` is dropped
    pub unsafe fn as_ptr(&self) -> *mut htslib::bcf_hdr_t {
        self.inner
    }

    #[inline]
    fn inner(&self) -> htslib::bcf_hdr_t {
        // SAFETY: self.inner is non-null, initialized by constructor or from_ptr.
        unsafe { *self.inner }
    }

    /// Get the number of samples defined in the header.
    pub fn sample_count(&self) -> u32 {
        self.inner().n[htslib::BCF_DT_SAMPLE as usize] as u32
    }

    /// Get vector of sample names defined in the header.
    pub fn samples(&self) -> Vec<&[u8]> {
        // SAFETY: inner().samples is a valid pointer to sample_count() C strings, set by htslib.
        let names =
            unsafe { slice::from_raw_parts(self.inner().samples, self.sample_count() as usize) };
        names
            .iter()
            // SAFETY: each element is a non-null, NUL-terminated C string from htslib.
            .map(|name| unsafe { ffi::CStr::from_ptr(*name).to_bytes() })
            .collect()
    }

    /// Obtain id (column index) of given sample.
    /// Returns `None` if sample is not present in header.
    pub fn sample_id(&self, sample: &[u8]) -> Option<usize> {
        self.samples().iter().position(|s| *s == sample)
    }

    /// Get the number of contigs defined in the header.
    pub fn contig_count(&self) -> u32 {
        self.inner().n[htslib::BCF_DT_CTG as usize] as u32
    }

    pub fn rid2name(&self, rid: u32) -> Result<&[u8]> {
        if rid < self.contig_count() {
            // SAFETY: rid is bounds-checked above; dict entries and their keys are valid htslib pointers.
            unsafe {
                let dict = self.inner().id[htslib::BCF_DT_CTG as usize];
                let ptr = (*dict.offset(rid as isize)).key;
                Ok(ffi::CStr::from_ptr(ptr).to_bytes())
            }
        } else {
            Err(Error::UnknownRID { rid })
        }
    }

    /// Retrieve the (internal) chromosome identifier
    /// # Examples
    /// ```rust
    /// use cstr8::cstr8;
    /// use rust_htslib::bcf::header::Header;
    /// use rust_htslib::bcf::{Format, Writer};
    ///
    /// let mut header = Header::new();
    /// let contig_field = br#"##contig=<ID=foo,length=10>"#;
    /// header.push_record(contig_field);
    /// let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    /// let header_view = vcf.header();
    /// let rid = header_view.name2rid(cstr8!("foo")).unwrap();
    /// assert_eq!(rid, 0);
    /// // try and retrieve a contig not in the header
    /// let result = header_view.name2rid(cstr8!("bar"));
    /// assert!(result.is_err())
    /// ```
    /// # Errors
    /// If `name` does not match a chromosome currently in the VCF header, returns [`BcfError::UnknownContig`]
    pub fn name2rid(&self, name: &CStr8) -> Result<u32> {
        match self.dict_lookup(htslib::BCF_DT_CTG as usize, name) {
            Some(id) => Ok(id as u32),
            None => Err(Error::UnknownContig {
                contig: name.as_str().to_owned(),
            }),
        }
    }

    pub fn info_type(&self, tag: &CStr8) -> Result<(TagType, TagLength)> {
        self.tag_type(tag, htslib::BCF_HL_INFO)
    }

    pub fn format_type(&self, tag: &CStr8) -> Result<(TagType, TagLength)> {
        self.tag_type(tag, htslib::BCF_HL_FMT)
    }

    fn tag_type(&self, tag: &CStr8, hdr_type: ::libc::c_uint) -> Result<(TagType, TagLength)> {
        let tag_desc = || tag.as_str().to_owned();
        let id = self
            .dict_lookup(htslib::BCF_DT_ID as usize, tag)
            .ok_or_else(|| Error::UndefinedTag { tag: tag_desc() })?;

        // SAFETY: id is a valid index (from the header's own dictionary);
        // entry/val pointers are valid htslib internals.
        let (_type, length, num_values) = unsafe {
            let n = (*self.inner).n[htslib::BCF_DT_ID as usize] as usize;
            let entry = slice::from_raw_parts((*self.inner).id[htslib::BCF_DT_ID as usize], n);
            let d = (*entry[id as usize].val).info[hdr_type as usize];
            ((d >> 4) & 0xf, (d >> 8) & 0xf, d >> 12)
        };
        let _type = match _type as ::libc::c_uint {
            htslib::BCF_HT_FLAG => TagType::Flag,
            htslib::BCF_HT_INT => TagType::Integer,
            htslib::BCF_HT_REAL => TagType::Float,
            htslib::BCF_HT_STR => TagType::String,
            _ => return Err(Error::UnexpectedType { tag: tag_desc() }),
        };
        let length = match length as ::libc::c_uint {
            htslib::BCF_VL_FIXED => TagLength::Fixed(num_values as u32),
            htslib::BCF_VL_VAR => TagLength::Variable,
            htslib::BCF_VL_A => TagLength::AltAlleles,
            htslib::BCF_VL_R => TagLength::Alleles,
            htslib::BCF_VL_G => TagLength::Genotypes,
            _ => return Err(Error::UnexpectedType { tag: tag_desc() }),
        };

        Ok((_type, length))
    }

    /// Convert string ID (e.g., for a `FILTER` value) to its numeric identifier.
    pub fn name_to_id(&self, id: &CStr8) -> Result<Id> {
        match self.dict_lookup(htslib::BCF_DT_ID as usize, id) {
            Some(i) => Ok(Id(i as u32)),
            None => Err(Error::UnknownID { id: id.into() }),
        }
    }

    /// Convert integer representing an identifier (e.g., a `FILTER` value) to its string
    /// name.
    ///
    pub fn id_to_name(&self, id: Id) -> Result<Vec<u8>> {
        // SAFETY: self.inner is non-null (from constructor).
        let n = unsafe { (*self.inner).n[htslib::BCF_DT_ID as usize] } as u32;
        if *id >= n {
            return Err(Error::UnknownID {
                id: format!("{}", *id),
            });
        }
        // SAFETY: id is bounds-checked above; dict entry key is a valid NUL-terminated C string from htslib.
        let key = unsafe {
            ffi::CStr::from_ptr(
                (*(*self.inner).id[htslib::BCF_DT_ID as usize].offset(*id as isize)).key,
            )
        };
        Ok(key.to_bytes().to_vec())
    }

    /// Convert string sample name to its numeric identifier.
    pub fn sample_to_id(&self, id: &CStr8) -> Result<Id> {
        match self.dict_lookup(htslib::BCF_DT_SAMPLE as usize, id) {
            Some(i) => Ok(Id(i as u32)),
            None => Err(Error::UnknownSample {
                name: id.as_str().to_owned(),
            }),
        }
    }

    /// Convert integer representing a sample to its name.
    pub fn id_to_sample(&self, id: Id) -> Result<Vec<u8>> {
        if *id >= self.sample_count() {
            return Err(Error::UnknownSample {
                name: format!("{}", *id),
            });
        }
        // SAFETY: id is bounds-checked above; dict entry key is a valid NUL-terminated C string from htslib.
        let key = unsafe {
            ffi::CStr::from_ptr(
                (*(*self.inner).id[htslib::BCF_DT_SAMPLE as usize].offset(*id as isize)).key,
            )
        };
        Ok(key.to_bytes().to_vec())
    }

    /// Return structured `HeaderRecord`s.
    pub fn header_records(&self) -> Vec<HeaderRecord> {
        fn parse_kv(rec: &htslib::bcf_hrec_t) -> LinearMap<String, String> {
            let mut result: LinearMap<String, String> = LinearMap::new();
            for i in 0_i32..(rec.nkeys) {
                // SAFETY: i is within [0, nkeys); keys[i] is a valid NUL-terminated C string from htslib.
                let key = unsafe {
                    ffi::CStr::from_ptr(*rec.keys.offset(i as isize))
                        .to_str()
                        .unwrap()
                        .to_string()
                };
                // SAFETY: i is within [0, nkeys); vals[i] is a valid NUL-terminated C string from htslib.
                let value = unsafe {
                    ffi::CStr::from_ptr(*rec.vals.offset(i as isize))
                        .to_str()
                        .unwrap()
                        .to_string()
                };
                result.insert(key, value);
            }
            result
        }

        let mut result: Vec<HeaderRecord> = Vec::new();
        // SAFETY: self.inner is non-null; nhrec and hrec are valid htslib header internals.
        for i in 0_i32..unsafe { (*self.inner).nhrec } {
            // SAFETY: i is within [0, nhrec); hrec[i] is a valid pointer to a bcf_hrec_t.
            let rec = unsafe { &(**(*self.inner).hrec.offset(i as isize)) };
            // SAFETY: rec.key is a valid NUL-terminated C string from htslib.
            let key = unsafe { ffi::CStr::from_ptr(rec.key).to_str().unwrap().to_string() };
            let record = match rec.type_ {
                0 => HeaderRecord::Filter {
                    key,
                    values: parse_kv(rec),
                },
                1 => HeaderRecord::Info {
                    key,
                    values: parse_kv(rec),
                },
                2 => HeaderRecord::Format {
                    key,
                    values: parse_kv(rec),
                },
                3 => HeaderRecord::Contig {
                    key,
                    values: parse_kv(rec),
                },
                4 => HeaderRecord::Structured {
                    key,
                    values: parse_kv(rec),
                },
                5 => HeaderRecord::Generic {
                    key,
                    // SAFETY: rec.value is a valid NUL-terminated C string for generic header records.
                    value: unsafe { ffi::CStr::from_ptr(rec.value).to_str().unwrap().to_string() },
                },
                _ => panic!("Unknown type: {}", rec.type_),
            };
            result.push(record);
        }
        result
    }

    /// Create an empty record using this header view.
    ///
    /// The record can be reused multiple times.
    pub fn empty_record(self: &Arc<Self>) -> crate::bcf::Record {
        crate::bcf::Record::new(self.clone())
    }
}

impl Clone for HeaderView {
    fn clone(&self) -> Self {
        // SAFETY: self.inner is non-null (from constructor); bcf_hdr_dup returns a new copy.
        let inner = unsafe { htslib::bcf_hdr_dup(self.inner) };
        // Rebuild caches from the new (duplicated) header pointer.
        let id_cache = unsafe { build_id_caches(inner) };
        HeaderView { inner, id_cache }
    }
}

impl Drop for HeaderView {
    fn drop(&mut self) {
        // SAFETY: self.inner was allocated by bcf_hdr_dup or htslib; bcf_hdr_destroy is symmetric.
        unsafe {
            htslib::bcf_hdr_destroy(self.inner);
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum TagType {
    Flag,
    Integer,
    Float,
    String,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum TagLength {
    Fixed(u32),
    AltAlleles,
    Alleles,
    Genotypes,
    Variable,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bcf::{Read, Reader};
    use crate::htslib;

    #[test]
    fn test_header_view_empty_record() {
        // Open a VCF file to get a HeaderView
        let vcf = Reader::from_path("test/test_string.vcf").expect("Error opening file");
        let header_view = vcf.header.clone();

        // Create an empty record from the HeaderView
        let record = header_view.empty_record();
        eprintln!("{:?}", record.rid());

        // Verify the record is properly initialized with default/empty values
        assert_eq!(record.rid(), Some(0)); // No chromosome/contig set
        assert_eq!(record.pos(), 0); // No position set
        assert_eq!(record.qual(), 0.0); // No quality score set
    }

    #[test]
    fn test_header_add_sample_via_raw_pointer() {
        let sample_name = b"test-sample";

        let header = Header::new();
        let sample = std::ffi::CString::new(sample_name).unwrap();

        let view = unsafe {
            let ptr = header.inner_ptr();
            // to avoid double free, as we will wrap this later in a HeaderView
            std::mem::forget(header);
            htslib::bcf_hdr_add_sample(ptr, sample.as_ptr());
            htslib::bcf_hdr_sync(ptr);
            // When the HeaderView is dropped, the bcf_hdr is freed
            HeaderView::from_ptr(ptr)
        };

        assert_eq!(view.samples(), vec![sample_name]);
    }

    #[test]
    fn test_header_view_version_via_raw_pointer() {
        let vcf = Reader::from_path("test/test_string.vcf").expect("Error opening file");
        let hv = vcf.header.clone();

        let version = unsafe {
            // the header view will outlive this pointer
            let ptr = hv.as_ptr();
            let version_charptr = htslib::bcf_hdr_get_version(ptr);
            std::ffi::CStr::from_ptr(version_charptr).to_str().unwrap()
        };

        assert_eq!(version, "VCFv4.1");
    }

    #[test]
    fn test_rid2name_valid() {
        let vcf = Reader::from_path("test/test_string.vcf").expect("Error opening file");
        let header = vcf.header();
        assert!(header.contig_count() > 0);
        for rid in 0..header.contig_count() {
            assert!(header.rid2name(rid).is_ok());
        }
    }

    #[test]
    fn test_rid2name_out_of_bounds() {
        let vcf = Reader::from_path("test/test_string.vcf").expect("Error opening file");
        let header = vcf.header();
        // Exactly at contig_count should fail (was an off-by-one bug)
        assert!(header.rid2name(header.contig_count()).is_err());
        assert!(header.rid2name(header.contig_count() + 1).is_err());
        assert!(header.rid2name(u32::MAX).is_err());
    }

    #[test]
    fn test_id_to_name_valid_roundtrip() {
        let vcf = Reader::from_path("test/test_string.vcf").expect("Error opening file");
        let header = vcf.header();
        // PASS filter always exists as id 0
        let name = header.id_to_name(Id(0)).unwrap();
        assert!(!name.is_empty());
    }

    #[test]
    fn test_id_to_name_out_of_bounds() {
        let vcf = Reader::from_path("test/test_string.vcf").expect("Error opening file");
        let header = vcf.header();
        assert!(header.id_to_name(Id(u32::MAX)).is_err());
    }

    #[test]
    fn test_id_to_sample_valid() {
        let vcf = Reader::from_path("test/test_string.vcf").expect("Error opening file");
        let header = vcf.header();
        assert!(header.sample_count() > 0);
        for i in 0..header.sample_count() {
            let name = header.id_to_sample(Id(i)).unwrap();
            assert!(!name.is_empty());
        }
    }

    #[test]
    fn test_id_to_sample_out_of_bounds() {
        let vcf = Reader::from_path("test/test_string.vcf").expect("Error opening file");
        let header = vcf.header();
        assert!(header.id_to_sample(Id(header.sample_count())).is_err());
    }
}

#[cfg(test)]
mod bcf_hdr_id2int_tests {
    use super::*;
    use crate::bcf::{Format, Read, Reader, Writer};
    use proptest::prelude::*;
    use std::ffi;
    use std::os::raw::c_char;
    use tempfile::NamedTempFile;

    /// Call C bcf_hdr_id2int as the oracle.
    fn id2int_c(hdr: &HeaderView, which: i32, name: &[u8]) -> i32 {
        let c_str = ffi::CString::new(name).unwrap();
        unsafe { htslib::bcf_hdr_id2int(hdr.inner, which, c_str.as_ptr() as *mut c_char) }
    }

    /// Strategy for valid BCF tag/contig/sample names: alphanumeric + underscore, 1-20 chars.
    fn name_strategy() -> impl Strategy<Value = String> {
        "[A-Za-z][A-Za-z0-9_]{0,19}"
    }

    /// Strategy for a set of unique names (1-8 entries).
    fn names_strategy() -> impl Strategy<Value = Vec<String>> {
        proptest::collection::hash_set(name_strategy(), 1..=8).prop_map(|s| s.into_iter().collect())
    }

    /// Build a header with the given contigs, INFO tags, FILTER tags, FORMAT tags, and samples.
    fn build_header(
        contigs: &[String],
        info_tags: &[String],
        filter_tags: &[String],
        format_tags: &[String],
        samples: &[String],
    ) -> (NamedTempFile, Arc<HeaderView>) {
        let tmp = NamedTempFile::new().unwrap();
        let mut header = Header::new();
        for c in contigs {
            header.push_record(format!("##contig=<ID={c},length=1000>").as_bytes());
        }
        for t in info_tags {
            header.push_record(
                format!("##INFO=<ID={t},Number=1,Type=Integer,Description=\"test\">").as_bytes(),
            );
        }
        for t in filter_tags {
            header.push_record(format!("##FILTER=<ID={t},Description=\"test\">").as_bytes());
        }
        for t in format_tags {
            header.push_record(
                format!("##FORMAT=<ID={t},Number=1,Type=Integer,Description=\"test\">").as_bytes(),
            );
        }
        for s in samples {
            header.push_sample(s.as_bytes());
        }
        let writer = Writer::from_path(tmp.path(), &header, true, Format::Bcf).unwrap();
        let hdr = writer.header().clone();
        (tmp, Arc::new(hdr))
    }

    // ---- Differential proptests: Rust HeaderView methods vs C bcf_hdr_id2int ----

    proptest! {
        /// name2rid must match C for all contigs present in the header.
        #[test]
        fn name2rid_matches_c_for_present_contigs(
            contigs in names_strategy(),
        ) {
            let (_tmp, hdr) = build_header(&contigs, &[], &[], &[], &[]);
            for name in &contigs {
                let c_result = id2int_c(&hdr, htslib::BCF_DT_CTG as i32, name.as_bytes());
                let cs = cstr8::CString8::new(name.as_str()).expect("valid CString8");
                let rs_result = hdr.name2rid(&cs);
                match rs_result {
                    Ok(id) => prop_assert_eq!(c_result, id as i32, "contig {}", name),
                    Err(_) => prop_assert_eq!(c_result, -1, "contig {} expected not found", name),
                }
            }
        }

        /// name2rid must match C for absent contigs.
        #[test]
        fn name2rid_matches_c_for_absent_contigs(
            contigs in names_strategy(),
            query in name_strategy(),
        ) {
            let (_tmp, hdr) = build_header(&contigs, &[], &[], &[], &[]);
            if !contigs.contains(&query) {
                let c_result = id2int_c(&hdr, htslib::BCF_DT_CTG as i32, query.as_bytes());
                let cs = cstr8::CString8::new(query.as_str()).expect("valid CString8");
                let rs_result = hdr.name2rid(&cs);
                prop_assert_eq!(c_result, -1);
                prop_assert!(rs_result.is_err());
            }
        }

        /// name_to_id (BCF_DT_ID) must match C for INFO/FILTER/FORMAT tags.
        #[test]
        fn name_to_id_matches_c_for_present_tags(
            info_tags in names_strategy(),
            filter_tags in names_strategy(),
            format_tags in names_strategy(),
        ) {
            let (_tmp, hdr) = build_header(&["chr1".into()], &info_tags, &filter_tags, &format_tags, &[]);
            // All INFO, FILTER, FORMAT tags should be findable
            for tags in [&info_tags, &filter_tags, &format_tags] {
                for name in tags {
                    let c_result = id2int_c(&hdr, htslib::BCF_DT_ID as i32, name.as_bytes());
                    let id = cstr8::CString8::new(name.as_str())
                        .expect("valid CString8");
                    let rs_result = hdr.name_to_id(&id);
                    match rs_result {
                        Ok(id) => prop_assert_eq!(c_result, *id as i32, "tag {}", name),
                        Err(_) => prop_assert_eq!(c_result, -1, "tag {} expected not found", name),
                    }
                }
            }
        }

        /// name_to_id must match C for absent tags.
        #[test]
        fn name_to_id_matches_c_for_absent_tags(
            info_tags in names_strategy(),
            query in name_strategy(),
        ) {
            let (_tmp, hdr) = build_header(&["chr1".into()], &info_tags, &[], &[], &[]);
            if !info_tags.contains(&query) {
                let c_result = id2int_c(&hdr, htslib::BCF_DT_ID as i32, query.as_bytes());
                // Only check queries that C also doesn't find (some built-in tags like PASS exist)
                if c_result == -1 {
                    let id = cstr8::CString8::new(query.as_str()).expect("valid CString8");
                    let rs_result = hdr.name_to_id(&id);
                    prop_assert!(rs_result.is_err());
                }
            }
        }

        /// sample_to_id must match C for all samples present in the header.
        #[test]
        fn sample_to_id_matches_c_for_present_samples(
            samples in names_strategy(),
        ) {
            let (_tmp, hdr) = build_header(&["chr1".into()], &[], &[], &[], &samples);
            for name in &samples {
                let c_result = id2int_c(&hdr, htslib::BCF_DT_SAMPLE as i32, name.as_bytes());
                let cs = cstr8::CString8::new(name.as_str()).expect("valid CString8");
                let rs_result = hdr.sample_to_id(&cs);
                match rs_result {
                    Ok(id) => prop_assert_eq!(c_result, *id as i32, "sample {}", name),
                    Err(_) => prop_assert_eq!(c_result, -1, "sample {} expected not found", name),
                }
            }
        }

        /// sample_to_id must match C for absent samples.
        #[test]
        fn sample_to_id_matches_c_for_absent_samples(
            samples in names_strategy(),
            query in name_strategy(),
        ) {
            let (_tmp, hdr) = build_header(&["chr1".into()], &[], &[], &[], &samples);
            if !samples.contains(&query) {
                let c_result = id2int_c(&hdr, htslib::BCF_DT_SAMPLE as i32, query.as_bytes());
                let cs = cstr8::CString8::new(query.as_str()).expect("valid CString8");
                let rs_result = hdr.sample_to_id(&cs);
                prop_assert_eq!(c_result, -1);
                prop_assert!(rs_result.is_err());
            }
        }
    }

    // ---- Deterministic tests with real VCF files ----

    #[test]
    fn id2int_matches_c_on_real_vcf() {
        let vcf = Reader::from_path("test/test_string.vcf").expect("Error opening file");
        let hdr = vcf.header();

        // Test all contigs
        for rid in 0..hdr.contig_count() {
            let name = hdr.rid2name(rid).unwrap();
            let c_result = id2int_c(hdr, htslib::BCF_DT_CTG as i32, name);
            assert_eq!(c_result, rid as i32, "contig roundtrip for rid {rid}");
        }

        // Test all samples
        for i in 0..hdr.sample_count() {
            let name = hdr.id_to_sample(Id(i)).unwrap();
            let c_result = id2int_c(hdr, htslib::BCF_DT_SAMPLE as i32, &name);
            assert_eq!(c_result, i as i32, "sample roundtrip for {i}");
        }

        // Test absent names
        assert_eq!(
            id2int_c(hdr, htslib::BCF_DT_CTG as i32, b"nonexistent_contig"),
            -1
        );
        assert_eq!(
            id2int_c(hdr, htslib::BCF_DT_SAMPLE as i32, b"nonexistent_sample"),
            -1
        );
        assert_eq!(
            id2int_c(hdr, htslib::BCF_DT_ID as i32, b"nonexistent_tag"),
            -1
        );
    }
}
