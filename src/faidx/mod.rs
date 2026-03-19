// Copyright 2020 Manuel Landesfeind, Evotec International GmbH
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//!
//! Module for working with faidx-indexed FASTA files.
//!

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
    /// Create a new Reader from a path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FaidxError> {
        let path = path.as_ref();
        let cpath = path_to_cstr8(path, true)?;
        let inner = unsafe { htslib::fai_load(cpath.as_ptr().cast()) };
        if inner.is_null() {
            // path_to_cstr8 succeeded, so to_string_lossy is lossless here
            return Err(FaidxError::Open {
                path: path.to_string_lossy().into_owned(),
            });
        }
        Ok(Self { inner })
    }

    /// Create a new Reader from a URL.
    pub fn from_url(url: &url::Url) -> Result<Self, FaidxError> {
        let url_str = url.as_str();
        let cpath = CString8::new(url_str).map_err(|_| FaidxError::Open {
            path: url_str.to_owned(),
        })?;
        let inner = unsafe { htslib::fai_load(cpath.as_ptr().cast()) };
        if inner.is_null() {
            return Err(FaidxError::Open {
                path: url_str.to_owned(),
            });
        }
        Ok(Self { inner })
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
        let n = unsafe { htslib::faidx_nseq(self.inner) };
        n.max(0) as u64
    }

    /// Fetches the i-th sequence name.
    pub fn seq_name(&self, i: i32) -> Result<String, FaidxError> {
        let ptr = unsafe { htslib::faidx_iseq(self.inner, i) };
        if ptr.is_null() {
            return Err(FaidxError::InvalidSequenceName { index: i });
        }
        let cname = unsafe { std::ffi::CStr::from_ptr(ptr) };
        cname
            .to_str()
            .map(|s| s.to_owned())
            .map_err(|_| FaidxError::InvalidSequenceName { index: i })
    }

    /// Fetches the length of the given sequence name.
    ///
    /// Returns `None` if the sequence is not found in the index.
    pub fn fetch_seq_len<N: AsRef<str>>(&self, name: N) -> Option<u64> {
        let cname = CString8::new(name.as_ref()).ok()?;
        let seq_len = unsafe { htslib::faidx_seq_len64(self.inner, cname.as_ptr().cast()) };
        if seq_len < 0 {
            None
        } else {
            Some(seq_len as u64)
        }
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
        unsafe {
            htslib::fai_destroy(self.inner);
        }
    }
}

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
