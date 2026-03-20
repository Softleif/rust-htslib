// Copyright 2018 Manuel Holtgrewe, Berlin Institute of Health.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Module for working with tabix-indexed text files.
//!
//! This module allows to read tabix-indexed text files (such as BED) in a convenient but in a
//! line-based (and thus format-agnostic way). For accessing tabix-inxed VCF files, using the
//! `bcf` module is probably a better choice as this module gives you lines from the text files
//! which you then have to take care of parsing.
//!
//! In general, for reading tabix-indexed files, first to open the file by creating a `tbx::Reader`
//! objects, possibly translate the chromosome name to its numeric ID in the file, fetch the region
//! of interest using `fetch()`, and finally iterate over the records using `records()`.
//!
//! # Examples
//!
//! ```rust,no_run
//! use cstr8::cstr8;
//! use rust_htslib::tbx::{self, Read};
//!
//! // Create a tabix reader for reading a tabix-indexed BED file.
//! let path_bed = "file.bed.gz";
//! let mut tbx_reader = tbx::Reader::from_path(&path_bed)
//!     .expect(&format!("Could not open {}", path_bed));
//!
//! // Resolve chromosome name to numeric ID.
//! let tid = match tbx_reader.tid(cstr8!("chr1")) {
//!     Ok(tid) => tid,
//!     Err(_) => panic!("Could not resolve 'chr1' to contig ID"),
//! };
//!
//! // Set region to fetch.
//! tbx_reader
//!     .fetch(tid, 0, 100_000)
//!     .expect("Could not seek to chr1:1-100,000");
//!
//! // Read through all records in region.
//! for record in tbx_reader.records() {
//!     // ... actually do some work
//! }
//! ```

use std::ffi;
use std::path::{Path, PathBuf};
use std::ptr;
use url::Url;

use crate::htslib;
use cstr8::{CStr8, CompactCStr8};

#[derive(thiserror::Error, Debug, PartialEq)]
pub enum TbxError {
    #[error("file not found: {0:?}")]
    FileNotFound(PathBuf),
    #[error("non-unicode path: {0:?}")]
    NonUnicodePath(PathBuf),
    #[error("path or sequence name contains interior null byte")]
    NullByte,
    #[error("file is not BGZF-compressed: {0:?}")]
    NotBgzf(PathBuf),
    #[error("invalid or missing tabix index for {0:?}")]
    InvalidIndex(PathBuf),
    #[error("failed to build tabix index for {0:?}")]
    BuildIndexFailed(PathBuf),
    #[error("sequence {name} not found in index")]
    SequenceNotFound { name: String },
    #[error("failed to fetch region")]
    Fetch,
    #[error("no active iterator — call fetch() first")]
    NoIter,
    #[error("truncated tabix record")]
    TruncatedRecord,
    #[error("error setting threads for tabix file reading")]
    SetThreads,
}

type Result<T> = std::result::Result<T, TbxError>;

fn path_to_cstring(path: &Path, must_exist: bool) -> Result<ffi::CString> {
    if must_exist && !path.exists() {
        return Err(TbxError::FileNotFound(path.to_owned()));
    }
    let s = path
        .to_str()
        .ok_or_else(|| TbxError::NonUnicodePath(path.to_owned()))?;
    ffi::CString::new(s).map_err(|_| TbxError::NullByte)
}

/// Preset configurations for common file formats.
#[derive(Debug, Clone, Copy)]
pub enum TabixFormat {
    Bed,
    Gff,
    Sam,
    Vcf,
}

impl TabixFormat {
    fn conf_ptr(self) -> *const htslib::tbx_conf_t {
        // SAFETY: tbx_conf_* are static constants in htslib; taking their
        // address is safe (the unsafe block is overly broad but harmless).
        unsafe {
            match self {
                TabixFormat::Bed => &htslib::tbx_conf_bed,
                TabixFormat::Gff => &htslib::tbx_conf_gff,
                TabixFormat::Sam => &htslib::tbx_conf_sam,
                TabixFormat::Vcf => &htslib::tbx_conf_vcf,
            }
        }
    }
}

/// Build a tabix index for a BGZF-compressed file.
///
/// This reads the file at `path`, builds a `.tbi` index, and writes it
/// alongside the original file (e.g. `file.bed.gz.tbi`).
///
/// # Arguments
///
/// * `path` - path to the BGZF-compressed file
/// * `format` - the file format preset (determines which columns are chrom/start/end)
/// * `n_threads` - number of threads for decompression (0 for single-threaded)
pub fn build_index<P: AsRef<Path>>(path: P, format: TabixFormat, n_threads: u32) -> Result<()> {
    let path = path.as_ref();
    let c_path = path_to_cstring(path, true)?;

    // SAFETY: c_path is a valid null-terminated CString; return checked below.
    let ret = unsafe {
        htslib::tbx_index_build3(
            c_path.as_ptr(),
            ptr::null(),
            0,
            n_threads as i32,
            format.conf_ptr(),
        )
    };

    match ret {
        0 => Ok(()),
        -2 => Err(TbxError::NotBgzf(path.to_owned())),
        _ => Err(TbxError::BuildIndexFailed(path.to_owned())),
    }
}

/// A trait for a Tabix reader with a read method.
pub trait Read: Sized {
    /// Read next line into the given `Vec<u8>` (i.e., ASCII string).
    ///
    /// Use this method in combination with a single allocated record to avoid the reallocations
    /// occurring with the iterator.
    ///
    /// # Arguments
    ///
    /// * `record` - the `Vec<u8>` to be filled
    ///
    /// # Returns
    /// Ok(true) if record was read, Ok(false) if no more record in file
    fn read(&mut self, record: &mut Vec<u8>) -> Result<bool>;

    /// Iterator over the lines/records of the seeked region.
    ///
    /// Note that, while being convenient, this is less efficient than pre-allocating a
    /// `Vec<u8>` and reading into it with the `read()` method, since every iteration involves
    /// the allocation of a new `Vec<u8>`.
    fn records(&mut self) -> Records<'_, Self>;

    /// Return the text headers, split by line.
    fn header(&self) -> &Vec<String>;
}

/// A Tabix file reader.
///
/// This struct and its associated functions are meant for reading plain-text tabix indexed
/// by `tabix`.
///
/// Note that the `tabix` command from `htslib` can actually several more things, including
/// building indices and converting BCF to VCF text output.  Both is out of scope here.
#[derive(Debug)]
pub struct Reader {
    /// The header lines (if any).
    header: Vec<String>,

    /// The file to read from.
    hts_file: *mut htslib::htsFile,
    /// The file format information.
    hts_format: htslib::htsExactFormat,
    /// The tbx_t structure to read from.
    tbx: *mut htslib::tbx_t,
    /// The current buffer.
    buf: htslib::kstring_t,
    /// Iterator over the buffer.
    itr: Option<*mut htslib::hts_itr_t>,

    /// Cached sequence names (populated once at construction).
    cached_seqnames: Vec<CompactCStr8>,
    /// Cached name→tid lookup (populated once at construction).
    tid_by_name: std::collections::HashMap<CompactCStr8, u64>,

    /// The currently fetch region's tid.
    tid: i64,
    /// The currently fetch region's 0-based begin pos.
    start: i64,
    /// The currently fetch region's 0-based end pos.
    end: i64,
}

// SAFETY: Reader owns its htsFile and tbx_t exclusively; htslib tabix
// operations are not tied to a particular thread.
unsafe impl Send for Reader {}

/// Redefinition of `KS_SEP_LINE` from `htslib/kseq.h`.
const KS_SEP_LINE: i32 = 2;

impl Reader {
    /// Create a new Reader from path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let c_path = path_to_cstring(path, true)?;
        Self::new_from_cstring(c_path)
    }

    pub fn from_url(url: &Url) -> Result<Self> {
        let c_path = ffi::CString::new(url.as_str()).map_err(|_| TbxError::NullByte)?;
        Self::new_from_cstring(c_path)
    }

    fn new_from_cstring(c_path: ffi::CString) -> Result<Self> {
        let c_mode = ffi::CString::new("r").unwrap(); // safe: literal
                                                      // SAFETY: c_path and c_mode are valid CStrings; result is null-checked.
        let hts_file = unsafe { htslib::hts_open(c_path.as_ptr(), c_mode.as_ptr()) };
        if hts_file.is_null() {
            return Err(TbxError::InvalidIndex(PathBuf::from(
                c_path.to_string_lossy().as_ref(),
            )));
        }
        // SAFETY: hts_file is non-null (checked above); htsFile.format is a public field.
        let hts_format: u32 = unsafe { (*hts_file).format.format };

        // SAFETY: c_path is a valid CString; result is null-checked.
        let tbx = unsafe { htslib::tbx_index_load(c_path.as_ptr()) };
        if tbx.is_null() {
            return Err(TbxError::InvalidIndex(PathBuf::from(
                c_path.to_string_lossy().as_ref(),
            )));
        }
        let mut header = Vec::new();
        let mut buf = htslib::kstring_t {
            l: 0,
            m: 0,
            s: ptr::null_mut(),
        };
        // SAFETY: hts_file and tbx are non-null (checked above). hts_getline
        // fills buf.s; buf.l > 0 implies buf.s is non-null and points to a
        // valid C string. (*tbx).conf.meta_char is safe because tbx is valid.
        // FIXME: .unwrap() on .to_str() will panic if the header line contains
        // invalid UTF-8; should handle this gracefully.
        unsafe {
            while htslib::hts_getline(hts_file, KS_SEP_LINE, &mut buf) >= 0 {
                if buf.l > 0 && i32::from(*buf.s) == (*tbx).conf.meta_char {
                    header.push(String::from(ffi::CStr::from_ptr(buf.s).to_str().unwrap()));
                } else {
                    break;
                }
            }
        }

        // Build sequence name and tid caches from C tbx_seqnames + tbx_name2id.
        let (cached_seqnames, tid_by_name) = unsafe {
            let mut nseq: i32 = 0;
            let seqs = htslib::tbx_seqnames(tbx, &mut nseq);
            let mut names = Vec::with_capacity(nseq.max(0) as usize);
            let mut tids: std::collections::HashMap<CompactCStr8, u64> =
                std::collections::HashMap::with_capacity(nseq.max(0) as usize);
            for i in 0..nseq {
                let ptr = *seqs.offset(i as isize);
                if !ptr.is_null() {
                    if let Ok(key) = CompactCStr8::from_ptr(ptr as *const u8) {
                        let id = htslib::tbx_name2id(tbx, ptr);
                        if id >= 0 {
                            tids.insert(key.clone(), id as u64);
                        }
                        names.push(key);
                    }
                }
            }
            libc::free(seqs as *mut libc::c_void);
            (names, tids)
        };

        Ok(Reader {
            header,
            hts_file,
            hts_format,
            tbx,
            buf,
            itr: None,
            cached_seqnames,
            tid_by_name,
            tid: -1,
            start: -1,
            end: -1,
        })
    }

    /// Get sequence/target ID from sequence name.
    pub fn tid(&self, name: &CStr8) -> Result<u64> {
        self.tid_by_name
            .get(name)
            .copied()
            .ok_or_else(|| TbxError::SequenceNotFound {
                name: name.as_str().to_owned(),
            })
    }

    /// Fetch region given by numeric sequence number and 0-based begin and end position.
    pub fn fetch(&mut self, tid: u64, start: u64, end: u64) -> Result<()> {
        self.tid = tid as i64;
        self.start = start as i64;
        self.end = end as i64;

        if let Some(itr) = self.itr {
            // SAFETY: itr is non-null (checked via Some).
            unsafe {
                htslib::hts_itr_destroy(itr);
            }
        }
        // SAFETY: self.tbx is non-null; (*self.tbx).idx is valid; result null-checked.
        let itr = unsafe {
            htslib::hts_itr_query(
                (*self.tbx).idx,
                tid as i32,
                start as i64,
                end as i64,
                Some(htslib::tbx_readrec),
            )
        };
        if itr.is_null() {
            self.itr = None;
            Err(TbxError::Fetch)
        } else {
            self.itr = Some(itr);
            Ok(())
        }
    }

    /// Return the sequence contig names.
    pub fn seqnames(&self) -> &[CompactCStr8] {
        &self.cached_seqnames
    }

    /// Activate multi-threaded BGZF read support in htslib. This should permit faster
    /// reading of large BGZF files.
    ///
    /// # Arguments
    ///
    /// * `n_threads` - number of extra background reader threads to use
    pub fn set_threads(&mut self, n_threads: usize) -> Result<()> {
        assert!(n_threads > 0, "n_threads must be > 0");

        // SAFETY: self.hts_file is non-null; n_threads validated > 0.
        let r = unsafe { htslib::hts_set_threads(self.hts_file, n_threads as i32) };
        if r != 0 {
            Err(TbxError::SetThreads)
        } else {
            Ok(())
        }
    }

    pub fn hts_format(&self) -> htslib::htsExactFormat {
        self.hts_format
    }
}

/// Return whether the two given genomic intervals overlap.
fn overlap(tid1: i64, begin1: i64, end1: i64, tid2: i64, begin2: i64, end2: i64) -> bool {
    (tid1 == tid2) && (begin1 < end2) && (begin2 < end1)
}

impl Read for Reader {
    fn read(&mut self, record: &mut Vec<u8>) -> Result<bool> {
        match self.itr {
            Some(itr) => {
                loop {
                    // Try to read next line.
                    // SAFETY: all pointers non-null (hts_file from constructor,
                    // itr from Some, buf is &mut, tbx from constructor). Return checked.
                    let ret = unsafe {
                        htslib::hts_itr_next(
                            (*self.hts_file).fp.bgzf,
                            itr,
                            &mut self.buf as *mut htslib::kstring_t as *mut libc::c_void,
                            self.tbx as *mut libc::c_void,
                        )
                    };
                    // Handle errors first.
                    if ret == -1 {
                        return Ok(false);
                    } else if ret == -2 {
                        return Err(TbxError::TruncatedRecord);
                    } else if ret < 0 {
                        panic!("Return value should not be <0 but was: {}", ret);
                    }
                    // Return first overlapping record (loop will stop when `hts_itr_next(...)`
                    // returns `< 0`).
                    // SAFETY: itr is non-null (checked via Some); reading struct fields.
                    let (tid, start, end) =
                        unsafe { ((*itr).curr_tid, (*itr).curr_beg, (*itr).curr_end) };
                    // XXX: Careful with this tid conversion!!!
                    if overlap(self.tid, self.start, self.end, tid as i64, start, end) {
                        // SAFETY: hts_itr_next returned >= 0, so buf.s is non-null
                        // and contains a valid C string.
                        // FIXME: .unwrap() on .to_str() will panic on non-UTF-8 data.
                        *record =
                            unsafe { Vec::from(ffi::CStr::from_ptr(self.buf.s).to_str().unwrap()) };
                        return Ok(true);
                    }
                }
            }
            _ => Err(TbxError::NoIter),
        }
    }

    fn records(&mut self) -> Records<'_, Self> {
        Records { reader: self }
    }

    fn header(&self) -> &Vec<String> {
        &self.header
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        // SAFETY: all pointers were allocated by htslib during construction;
        // itr checked via Some; destroy/close are symmetric with open/load.
        unsafe {
            if let Some(itr) = self.itr {
                htslib::hts_itr_destroy(itr);
            }
            htslib::tbx_destroy(self.tbx);
            htslib::hts_close(self.hts_file);
        }
    }
}

/// Iterator over the lines of a tabix file.
#[derive(Debug)]
pub struct Records<'a, R: Read> {
    reader: &'a mut R,
}

impl<R: Read> Iterator for Records<'_, R> {
    type Item = Result<Vec<u8>>;

    #[allow(clippy::read_zero_byte_vec)]
    fn next(&mut self) -> Option<Result<Vec<u8>>> {
        let mut record = Vec::new();
        match self.reader.read(&mut record) {
            Ok(false) => None,
            Ok(true) => Some(Ok(record)),
            Err(err) => Some(Err(err)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cstr8::cstr8;

    #[test]
    fn bed_basic() {
        let reader =
            Reader::from_path("test/tabix_reader/test_bed3.bed.gz").expect("Error opening file.");

        // Check sequence name vector.
        assert_eq!(reader.seqnames().len(), 2);
        assert_eq!(reader.seqnames()[0], *"chr1");
        assert_eq!(reader.seqnames()[1], *"chr2");

        // Check mapping between name and idx.
        assert_eq!(reader.tid(cstr8!("chr1")).unwrap(), 0);
        assert_eq!(reader.tid(cstr8!("chr2")).unwrap(), 1);
        assert!(reader.tid(cstr8!("chr3")).is_err());
    }

    #[test]
    fn bed_fetch_from_chr1_read_api() {
        let mut reader =
            Reader::from_path("test/tabix_reader/test_bed3.bed.gz").expect("Error opening file.");

        let chr1_id = reader.tid(cstr8!("chr1")).unwrap();
        assert!(reader.fetch(chr1_id, 1000, 1003).is_ok());

        let mut record = Vec::new();
        assert!(reader.read(&mut record).is_ok());
        assert_eq!(record, Vec::from("chr1\t1001\t1002"));
        assert_eq!(reader.read(&mut record), Ok(false)); // EOF
    }

    #[test]
    fn bed_fetch_from_chr1_iterator_api() {
        let mut reader =
            Reader::from_path("test/tabix_reader/test_bed3.bed.gz").expect("Error opening file.");

        let chr1_id = reader.tid(cstr8!("chr1")).unwrap();
        assert!(reader.fetch(chr1_id, 1000, 1003).is_ok());

        let records: Vec<Vec<u8>> = reader.records().map(|r| r.unwrap()).collect();
        assert_eq!(records, vec![Vec::from("chr1\t1001\t1002")]);
    }

    #[test]
    fn test_fails_on_bam() {
        let reader = Reader::from_path("test/test.bam");
        assert!(reader.is_err());
    }

    #[test]
    fn test_fails_on_non_existiant() {
        let reader = Reader::from_path("test/no_such_file");
        assert!(reader.is_err());
    }

    #[test]
    fn test_fails_on_vcf() {
        let reader = Reader::from_path("test/test_left.vcf");
        assert!(reader.is_err());
    }

    #[test]
    fn test_text_header_regions() {
        // This file has chromosome, start, and end positions with a header line.
        Reader::from_path("test/tabix_reader/genomic_regions_header.txt.gz")
            .expect("Error opening file.");
    }

    #[test]
    fn test_text_header_positions() {
        // This file has chromosome and position with a header line, indexed with
        // `tabix -b2 -e2 <file>`.
        Reader::from_path("test/tabix_reader/genomic_positions_header.txt.gz")
            .expect("Error opening file.");
    }

    #[test]
    fn test_text_bad_header() {
        // This is a duplicate of the above file but the index file is nonsense text.
        Reader::from_path("test/tabix_reader/bad_header.txt.gz")
            .expect_err("Invalid index file should fail.");
    }
}

#[cfg(test)]
mod tbx_accessor_tests {
    use super::*;
    use std::convert::TryFrom;

    /// Call C tbx_name2id as oracle.
    fn tid_c(reader: &Reader, name: &str) -> i32 {
        let c_str = ffi::CString::new(name).unwrap();
        unsafe { htslib::tbx_name2id(reader.tbx, c_str.as_ptr()) }
    }

    /// Call C tbx_seqnames as oracle.
    fn seqnames_c(reader: &Reader) -> Vec<String> {
        let mut nseq: i32 = 0;
        let seqs = unsafe { htslib::tbx_seqnames(reader.tbx, &mut nseq) };
        let mut result = Vec::new();
        for i in 0..nseq {
            unsafe {
                let name = ffi::CStr::from_ptr(*seqs.offset(i as isize))
                    .to_str()
                    .unwrap()
                    .to_owned();
                result.push(name);
            }
        }
        unsafe { libc::free(seqs as *mut libc::c_void) };
        result
    }

    #[test]
    fn tid_matches_c_for_all_seqnames() {
        let reader =
            Reader::from_path("test/tabix_reader/test_bed3.bed.gz").expect("Error opening file.");
        let names = reader.seqnames();
        for name in names {
            let c_tid = tid_c(&reader, name.as_str());
            let rs_tid = reader.tid(name);
            assert_eq!(rs_tid.unwrap(), c_tid as u64, "tid mismatch for {name}");
        }
    }

    #[test]
    fn tid_matches_c_for_absent_names() {
        let reader =
            Reader::from_path("test/tabix_reader/test_bed3.bed.gz").expect("Error opening file.");
        for name in ["chrZ", "nonexistent", "MT"] {
            let c_tid = tid_c(&reader, name);
            assert_eq!(c_tid, -1);
            let cname = CompactCStr8::try_from(name).unwrap();
            assert!(reader.tid(&cname).is_err());
        }
    }

    #[test]
    fn seqnames_matches_c() {
        let reader =
            Reader::from_path("test/tabix_reader/test_bed3.bed.gz").expect("Error opening file.");
        let c_names = seqnames_c(&reader);
        let rs_names: Vec<String> = reader
            .seqnames()
            .iter()
            .map(|s| s.as_str().to_owned())
            .collect();
        assert_eq!(c_names, rs_names);
    }

    #[test]
    fn seqnames_and_tid_roundtrip() {
        let reader =
            Reader::from_path("test/tabix_reader/test_bed3.bed.gz").expect("Error opening file.");
        for (i, name) in reader.seqnames().iter().enumerate() {
            assert_eq!(reader.tid(name).unwrap(), i as u64);
        }
    }

    #[test]
    fn hts_get_format_matches_direct_access() {
        let reader =
            Reader::from_path("test/tabix_reader/test_bed3.bed.gz").expect("Error opening file.");
        let c_format = unsafe {
            let fmt_ptr = htslib::hts_get_format(reader.hts_file);
            (*fmt_ptr).format
        };
        let rs_format = unsafe { (*reader.hts_file).format.format };
        assert_eq!(c_format, rs_format);
    }

    #[test]
    fn hts_get_bgzfp_matches_direct_access() {
        let reader =
            Reader::from_path("test/tabix_reader/test_bed3.bed.gz").expect("Error opening file.");
        let c_bgzf = unsafe { htslib::hts_get_bgzfp(reader.hts_file) };
        let rs_bgzf = unsafe { (*reader.hts_file).fp.bgzf };
        assert_eq!(c_bgzf, rs_bgzf);
    }
}
