//! Error types for the bwa-mem3-rs crate.

use std::path::PathBuf;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("failed to load index at {path}: {msg}")]
    IndexLoad { path: PathBuf, msg: String },

    #[error("invalid options: {0}")]
    InvalidOpts(String),

    #[error("invalid input: {0}")]
    InvalidInput(String),

    #[error("bwa-mem3 shim: {0}")]
    Shim(String),

    #[error("I/O: {0}")]
    Io(#[from] std::io::Error),

    /// An option field held a value this crate does not recognize — almost
    /// always because the vendored bwa-mem3 gained a new enum variant that the
    /// Rust mirror in `opts.rs` has not been taught yet.
    #[error(
        "unrecognized {kind} value {value} from bwa-mem3 (vendored upstream added a variant?)"
    )]
    UnrecognizedEnum { kind: &'static str, value: i32 },
}

impl From<std::ffi::NulError> for Error {
    fn from(e: std::ffi::NulError) -> Self {
        Error::InvalidInput(format!("path contains NUL byte: {e}"))
    }
}

/// Fetch the shim's thread-local last-error string and wrap it into [`Error::Shim`].
pub(crate) fn shim_err(context: &str) -> Error {
    let msg = unsafe {
        let ptr = bwa_mem3_sys::bwa_shim_last_error();
        if ptr.is_null() {
            "unknown shim error".to_string()
        } else {
            std::ffi::CStr::from_ptr(ptr).to_string_lossy().into_owned()
        }
    };
    Error::Shim(format!("{context}: {msg}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_invalid_input() {
        let e = Error::InvalidInput("empty batch".into());
        assert_eq!(e.to_string(), "invalid input: empty batch");
    }
}
