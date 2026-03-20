// Copyright 2020 Manuel Landesfeind, Evotec International GmbH
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//!
//! Module for working with faidx-indexed FASTA files.
//!

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use cstr8::CString8;

use crate::htslib;

#[derive(thiserror::Error, Debug, PartialEq)]
pub enum FaidxError {
    #[error("FASTA file not found: {0:?}")]
    FileNotFound(PathBuf),
    #[error("non-unicode path: {0:?}")]
    NonUnicodePath(PathBuf),
    #[error("failed to open FASTA/FAI: {path}")]
    Open { path: String },
    #[error("failed to build FASTA index for {path:?}")]
    Build { path: PathBuf },
    #[error("sequence not found: {name}")]
    SequenceNotFound { name: String },
    #[error("position too large for htslib (must fit in i64)")]
    PositionTooLarge,
    #[error("sequence name at index {index} not found or not valid UTF-8")]
    InvalidSequenceName { index: i32 },
    #[error("sequence name contains interior null byte")]
    NullByteName,
    #[error("fetched sequence is not valid UTF-8")]
    InvalidUtf8,
}

/// A Fasta reader.
#[derive(Debug)]
pub struct Reader {
    inner: *mut htslib::faidx_t,
    /// Cached sequence names, indexed by sequence ID.
    cached_names: Vec<String>,
    /// Cached sequence lengths, keyed by sequence name.
    cached_lengths: HashMap<String, u64>,
}

/// Convert a path to a C string suitable for htslib, with precise errors.
fn path_to_cstr8(path: &Path, must_exist: bool) -> Result<CString8, FaidxError> {
    if must_exist && !path.exists() {
        return Err(FaidxError::FileNotFound(path.to_owned()));
    }
    let path_str = path
        .to_str()
        .ok_or_else(|| FaidxError::NonUnicodePath(path.to_owned()))?;
    CString8::new(path_str).map_err(|_| FaidxError::NonUnicodePath(path.to_owned()))
}

/// Build a faidx for input path.
///
///```
/// use rust_htslib::faidx::build;
/// let path = std::path::PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"),"/test/test_cram.fa"));
/// build(&path).expect("Failed to build fasta index");
///```
pub fn build<P: AsRef<Path>>(path: P) -> Result<(), FaidxError> {
    let path = path.as_ref();
    let cpath = path_to_cstr8(path, true)?;
    // SAFETY: cpath is a valid null-terminated CString; return value checked below.
    let rc = unsafe { htslib::fai_build(cpath.as_ptr().cast()) };
    if rc < 0 {
        Err(FaidxError::Build {
            path: path.to_owned(),
        })
    } else {
        Ok(())
    }
}

impl Reader {
    /// Build the sequence name and length caches from the C faidx handle.
    ///
    /// Called once at construction. After this, `n_seqs()`, `seq_name()`, and
    /// `fetch_seq_len()` are served from pure Rust data structures.
    ///
    /// # Safety
    /// `inner` must be a valid, non-null pointer to an initialized `faidx_t`.
    unsafe fn build_caches(inner: *mut htslib::faidx_t) -> (Vec<String>, HashMap<String, u64>) {
        let n = htslib::faidx_nseq(inner).max(0) as usize;
        let mut names = Vec::with_capacity(n);
        let mut lengths = HashMap::with_capacity(n);
        for i in 0..n {
            let ptr = htslib::faidx_iseq(inner, i as i32);
            if ptr.is_null() {
                continue;
            }
            let name = match std::ffi::CStr::from_ptr(ptr).to_str() {
                Ok(s) => s.to_owned(),
                Err(_) => continue,
            };
            let cname = CString8::new(name.as_str()).unwrap();
            let len = htslib::faidx_seq_len64(inner, cname.as_ptr().cast());
            if len >= 0 {
                lengths.insert(name.clone(), len as u64);
            }
            names.push(name);
        }
        (names, lengths)
    }

    /// Create a new Reader from a path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FaidxError> {
        let path = path.as_ref();
        let cpath = path_to_cstr8(path, true)?;
        // SAFETY: cpath is a valid null-terminated CString; result is null-checked.
        let inner = unsafe { htslib::fai_load(cpath.as_ptr().cast()) };
        if inner.is_null() {
            return Err(FaidxError::Open {
                path: path.to_string_lossy().into_owned(),
            });
        }
        let (cached_names, cached_lengths) = unsafe { Self::build_caches(inner) };
        Ok(Self {
            inner,
            cached_names,
            cached_lengths,
        })
    }

    /// Create a new Reader from a URL.
    pub fn from_url(url: &url::Url) -> Result<Self, FaidxError> {
        let url_str = url.as_str();
        let cpath = CString8::new(url_str).map_err(|_| FaidxError::Open {
            path: url_str.to_owned(),
        })?;
        // SAFETY: cpath is a valid null-terminated CString; result is null-checked.
        let inner = unsafe { htslib::fai_load(cpath.as_ptr().cast()) };
        if inner.is_null() {
            return Err(FaidxError::Open {
                path: url_str.to_owned(),
            });
        }
        let (cached_names, cached_lengths) = unsafe { Self::build_caches(inner) };
        Ok(Self {
            inner,
            cached_names,
            cached_lengths,
        })
    }

    /// Fetch the sequence as a byte array.
    ///
    /// Coordinates are 0-based inclusive: `[begin, end]`.
    pub fn fetch_seq<N: AsRef<str>>(
        &self,
        name: N,
        begin: usize,
        end: usize,
    ) -> Result<Vec<u8>, FaidxError> {
        if begin > i64::MAX as usize || end > i64::MAX as usize {
            return Err(FaidxError::PositionTooLarge);
        }
        let cname = CString8::new(name.as_ref()).map_err(|_| FaidxError::NullByteName)?;
        let mut len_out: htslib::hts_pos_t = 0;
        // SAFETY: self.inner is valid (from constructor null-check); cname is a
        // valid CString; begin/end validated to fit in i64 above; result ptr and
        // len_out are checked below.
        let ptr = unsafe {
            htslib::faidx_fetch_seq64(
                self.inner,
                cname.as_ptr().cast(),
                begin as htslib::hts_pos_t,
                end as htslib::hts_pos_t,
                &mut len_out,
            )
        };
        if ptr.is_null() || len_out < 0 {
            return Err(FaidxError::SequenceNotFound {
                name: name.as_ref().to_owned(),
            });
        }
        // Copy out of C-allocated buffer and free with libc::free.
        // We must not use Vec::from_raw_parts because the pointer was allocated
        // by htslib's malloc, not Rust's global allocator. If a custom allocator
        // is active (e.g. mimalloc), dropping the Vec would be undefined behavior.
        let len = len_out as usize;
        // SAFETY: ptr is non-null (checked above), len is non-negative (checked
        // above). Immediately copied to Vec; ptr is then freed via libc::free
        // (matching htslib's malloc). Not using Vec::from_raw_parts because ptr
        // was allocated by htslib, not Rust's global allocator.
        let vec = unsafe { std::slice::from_raw_parts(ptr as *const u8, len) }.to_vec();
        unsafe { libc::free(ptr as *mut libc::c_void) };
        Ok(vec)
    }

    /// Fetches the sequence and returns it as a string.
    ///
    /// Coordinates are 0-based inclusive: `[begin, end]`.
    pub fn fetch_seq_string<N: AsRef<str>>(
        &self,
        name: N,
        begin: usize,
        end: usize,
    ) -> Result<String, FaidxError> {
        let bytes = self.fetch_seq(name, begin, end)?;
        String::from_utf8(bytes).map_err(|_| FaidxError::InvalidUtf8)
    }

    /// Fetches the number of sequences in the fai index.
    pub fn n_seqs(&self) -> u64 {
        self.cached_names.len() as u64
    }

    /// Fetches the i-th sequence name.
    pub fn seq_name(&self, i: i32) -> Result<String, FaidxError> {
        if i < 0 {
            return Err(FaidxError::InvalidSequenceName { index: i });
        }
        self.cached_names
            .get(i as usize)
            .cloned()
            .ok_or(FaidxError::InvalidSequenceName { index: i })
    }

    /// Fetches the length of the given sequence name.
    ///
    /// Returns `None` if the sequence is not found in the index.
    pub fn fetch_seq_len<N: AsRef<str>>(&self, name: N) -> Option<u64> {
        self.cached_lengths.get(name.as_ref()).copied()
    }

    /// Returns all sequence names.
    ///
    ///```
    /// use rust_htslib::faidx;
    /// let path = concat!(env!("CARGO_MANIFEST_DIR"),"/test/test_cram.fa");
    /// faidx::build(&path).expect("Failed to build fasta index");
    /// let reader = faidx::Reader::from_path(path).expect("Failed to open faidx");
    /// assert_eq!(reader.seq_names(), Ok(vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()]));
    ///```
    pub fn seq_names(&self) -> Result<Vec<String>, FaidxError> {
        let num_seq = self.n_seqs();
        let mut ret = Vec::with_capacity(num_seq as usize);
        for seq_id in 0..num_seq {
            ret.push(self.seq_name(seq_id as i32)?);
        }
        Ok(ret)
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        // SAFETY: self.inner was allocated by fai_load; fai_destroy is symmetric.
        unsafe {
            htslib::fai_destroy(self.inner);
        }
    }
}

// SAFETY: Reader owns its faidx_t exclusively; htslib faidx operations are
// not tied to a particular thread.
unsafe impl Send for Reader {}

#[cfg(test)]
mod tests {
    use super::*;

    fn open_reader() -> Reader {
        Reader::from_path(format!("{}/test/test_cram.fa", env!("CARGO_MANIFEST_DIR"))).unwrap()
    }

    #[test]
    fn faidx_open() {
        open_reader();
    }

    #[test]
    fn faidx_read_chr_first_base() {
        let r = open_reader();

        let bseq = r.fetch_seq("chr1", 0, 0).unwrap();
        assert_eq!(bseq.len(), 1);
        assert_eq!(bseq, b"G");

        let seq = r.fetch_seq_string("chr1", 0, 0).unwrap();
        assert_eq!(seq.len(), 1);
        assert_eq!(seq, "G");
    }

    #[test]
    fn faidx_read_chr_start() {
        let r = open_reader();

        let bseq = r.fetch_seq("chr1", 0, 9).unwrap();
        assert_eq!(bseq.len(), 10);
        assert_eq!(bseq, b"GGGCACAGCC");

        let seq = r.fetch_seq_string("chr1", 0, 9).unwrap();
        assert_eq!(seq.len(), 10);
        assert_eq!(seq, "GGGCACAGCC");
    }

    #[test]
    fn faidx_read_chr_between() {
        let r = open_reader();

        let bseq = r.fetch_seq("chr1", 4, 14).unwrap();
        assert_eq!(bseq.len(), 11);
        assert_eq!(bseq, b"ACAGCCTCACC");

        let seq = r.fetch_seq_string("chr1", 4, 14).unwrap();
        assert_eq!(seq.len(), 11);
        assert_eq!(seq, "ACAGCCTCACC");
    }

    #[test]
    fn faidx_read_chr_end() {
        let r = open_reader();

        let bseq = r.fetch_seq("chr1", 110, 120).unwrap();
        assert_eq!(bseq.len(), 10);
        assert_eq!(bseq, b"CCCCTCCGTG");

        let seq = r.fetch_seq_string("chr1", 110, 120).unwrap();
        assert_eq!(seq.len(), 10);
        assert_eq!(seq, "CCCCTCCGTG");
    }

    #[test]
    fn faidx_read_twice_string() {
        let r = open_reader();
        let seq = r.fetch_seq_string("chr1", 110, 120).unwrap();
        assert_eq!(seq.len(), 10);
        assert_eq!(seq, "CCCCTCCGTG");

        let seq = r.fetch_seq_string("chr1", 5, 9).unwrap();
        assert_eq!(seq.len(), 5);
        assert_eq!(seq, "CAGCC");
    }

    #[test]
    fn faidx_read_twice_bytes() {
        let r = open_reader();
        let seq = r.fetch_seq("chr1", 110, 120).unwrap();
        assert_eq!(seq.len(), 10);
        assert_eq!(seq, b"CCCCTCCGTG");

        let seq = r.fetch_seq("chr1", 5, 9).unwrap();
        assert_eq!(seq.len(), 5);
        assert_eq!(seq, b"CAGCC");
    }

    #[test]
    fn faidx_position_too_large() {
        let r = open_reader();
        let position_too_large = i64::MAX as usize;
        let res = r.fetch_seq("chr1", position_too_large, position_too_large + 1);
        assert!(matches!(res, Err(FaidxError::PositionTooLarge)));
    }

    #[test]
    fn faidx_n_seqs() {
        let r = open_reader();
        assert_eq!(r.n_seqs(), 3);
    }

    #[test]
    fn faidx_seq_name() {
        let r = open_reader();
        let n = r.seq_name(1).unwrap();
        assert_eq!(n, "chr2");
    }

    #[test]
    fn faidx_seq_name_all_valid_indices() {
        let r = open_reader();
        assert_eq!(r.seq_name(0).unwrap(), "chr1");
        assert_eq!(r.seq_name(1).unwrap(), "chr2");
        assert_eq!(r.seq_name(2).unwrap(), "chr3");
    }

    #[test]
    fn faidx_seq_name_boundary() {
        // n_seqs() == 3, so index 2 is the last valid and 3 is the first invalid.
        // faidx_iseq does no bounds checking in C, so calling it with an
        // out-of-bounds index would segfault without our Rust-side guard.
        let r = open_reader();
        assert!(r.seq_name(2).is_ok());
        assert!(matches!(
            r.seq_name(3),
            Err(FaidxError::InvalidSequenceName { index: 3 })
        ));
    }

    #[test]
    fn faidx_seq_name_i32_extremes() {
        let r = open_reader();
        assert!(matches!(
            r.seq_name(i32::MAX),
            Err(FaidxError::InvalidSequenceName { index: i32::MAX })
        ));
        assert!(matches!(
            r.seq_name(i32::MIN),
            Err(FaidxError::InvalidSequenceName { index: i32::MIN })
        ));
    }

    #[test]
    fn faidx_get_seq_len() {
        let r = open_reader();
        assert_eq!(r.fetch_seq_len("chr1"), Some(120));
        assert_eq!(r.fetch_seq_len("chr2"), Some(120));
    }

    #[test]
    fn open_many_readers() {
        for _ in 0..500_000 {
            let reader = open_reader();
            drop(reader);
        }
    }

    #[test]
    fn faidx_open_nonexistent_errors() {
        let result = Reader::from_path("/does/not/exist.fa");
        assert!(matches!(result, Err(FaidxError::FileNotFound(_))));
    }

    #[test]
    fn faidx_fetch_nonexistent_seq_errors() {
        let r = open_reader();
        let result = r.fetch_seq("nonexistent_chromosome", 0, 10);
        assert!(matches!(result, Err(FaidxError::SequenceNotFound { .. })));
    }

    #[test]
    fn faidx_fetch_nonexistent_seq_string_errors() {
        let r = open_reader();
        let result = r.fetch_seq_string("nonexistent_chromosome", 0, 10);
        assert!(matches!(result, Err(FaidxError::SequenceNotFound { .. })));
    }

    #[test]
    fn faidx_seq_name_out_of_bounds_errors() {
        let r = open_reader();
        assert!(matches!(
            r.seq_name(3),
            Err(FaidxError::InvalidSequenceName { index: 3 })
        ));
        assert!(matches!(
            r.seq_name(999),
            Err(FaidxError::InvalidSequenceName { index: 999 })
        ));
        assert!(matches!(
            r.seq_name(-1),
            Err(FaidxError::InvalidSequenceName { index: -1 })
        ));
    }

    #[test]
    fn faidx_fetch_seq_begin_after_end() {
        let r = open_reader();
        // htslib docs say behavior is undefined when begin > end,
        // but it should not segfault
        let result = r.fetch_seq("chr1", 10, 5);
        drop(result);
    }

    #[test]
    fn faidx_fetch_seq_past_chromosome_end() {
        let r = open_reader();
        // chr1 is 120bp; fetch beyond that — htslib clamps
        let result = r.fetch_seq("chr1", 0, 999);
        assert!(result.is_ok());
    }

    #[test]
    fn faidx_fetch_empty_name() {
        let r = open_reader();
        let result = r.fetch_seq("", 0, 10);
        assert!(result.is_err());
    }

    #[test]
    fn faidx_seq_len_nonexistent() {
        let r = open_reader();
        assert_eq!(r.fetch_seq_len("nonexistent"), None);
    }

    #[test]
    fn faidx_fetch_many_times_no_leak() {
        let r = open_reader();
        for _ in 0..100_000 {
            let seq = r.fetch_seq("chr1", 0, 119).unwrap();
            assert_eq!(seq.len(), 120);
        }
    }
}

#[cfg(test)]
mod faidx_cache_tests {
    use super::*;

    fn open_reader() -> Reader {
        Reader::from_path("test/test_cram.fa").expect("Error opening faidx")
    }

    /// Compare cached n_seqs against C faidx_nseq.
    #[test]
    fn n_seqs_matches_c() {
        let r = open_reader();
        let c_result = unsafe { htslib::faidx_nseq(r.inner) }.max(0) as u64;
        assert_eq!(r.n_seqs(), c_result);
    }

    /// Compare every cached seq_name against C faidx_iseq.
    #[test]
    fn seq_name_matches_c_for_all_indices() {
        let r = open_reader();
        let n = unsafe { htslib::faidx_nseq(r.inner) }.max(0);
        for i in 0..n {
            let c_ptr = unsafe { htslib::faidx_iseq(r.inner, i) };
            assert!(!c_ptr.is_null());
            let c_name = unsafe { std::ffi::CStr::from_ptr(c_ptr) }
                .to_str()
                .unwrap()
                .to_owned();
            let rs_name = r.seq_name(i).unwrap();
            assert_eq!(c_name, rs_name, "seq_name mismatch at index {i}");
        }
    }

    /// Compare cached fetch_seq_len against C faidx_seq_len64 for all sequences.
    #[test]
    fn fetch_seq_len_matches_c_for_all_sequences() {
        let r = open_reader();
        let n = r.n_seqs();
        for i in 0..n {
            let name = r.seq_name(i as i32).unwrap();
            let cname = CString8::new(name.as_str()).unwrap();
            let c_len = unsafe { htslib::faidx_seq_len64(r.inner, cname.as_ptr().cast()) };
            let rs_len = r.fetch_seq_len(&name);
            if c_len < 0 {
                assert_eq!(rs_len, None, "expected None for {name}");
            } else {
                assert_eq!(rs_len, Some(c_len as u64), "length mismatch for {name}");
            }
        }
    }

    /// Absent sequences must return None (matching C returning -1).
    #[test]
    fn fetch_seq_len_absent_matches_c() {
        let r = open_reader();
        let cname = CString8::new("nonexistent_chr").unwrap();
        let c_len = unsafe { htslib::faidx_seq_len64(r.inner, cname.as_ptr().cast()) };
        assert_eq!(c_len, -1);
        assert_eq!(r.fetch_seq_len("nonexistent_chr"), None);
    }

    /// Out-of-bounds seq_name must return error (matching C behavior).
    #[test]
    fn seq_name_out_of_bounds() {
        let r = open_reader();
        let n = r.n_seqs();
        assert!(r.seq_name(n as i32).is_err());
        assert!(r.seq_name(-1).is_err());
        assert!(r.seq_name(i32::MAX).is_err());
    }
}
