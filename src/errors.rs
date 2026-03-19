use std::path::PathBuf;
use thiserror::Error;

/// Generic result type for functions in this crate with
/// a global error class.
pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Error, Debug, PartialEq)]
pub enum Error {
    // General errors (used by bgzf and other modules not yet refactored)
    #[error("file not found: {path}")]
    FileNotFound { path: PathBuf },
    #[error("file could not be opened: {path}")]
    FileOpen { path: String },
    #[error("invalid (non-unicode) characters in path")]
    NonUnicodePath,
    #[error("failed to create htslib thread pool")]
    ThreadPool,

    #[error("invalid compression level {level}")]
    BgzfInvalidCompressionLevel { level: i8 },

    // Module-specific errors (transparent wrappers)
    #[error(transparent)]
    Faidx(#[from] crate::faidx::FaidxError),
    #[error(transparent)]
    Tbx(#[from] crate::tbx::TbxError),
    #[error(transparent)]
    Bcf(#[from] crate::bcf::BcfError),
    #[error(transparent)]
    Bam(#[from] crate::bam::BamError),
}
