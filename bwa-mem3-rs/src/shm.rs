//! Shared-memory index lifecycle.
//!
//! bwa-mem3 can stage a packed index in POSIX shared memory under
//! `/bwaidx-<basename>` (with a `/bwactl` registry segment), letting many
//! processes attach a single ~6 GB index instead of each loading from
//! disk.
//!
//! Lifecycle: [`stage`] publishes a prefix, [`destroy`] drops every
//! staged segment, [`list`] prints `<name>\t<bytes>` for each entry to
//! stdout (matching `bwa shm -l`), [`is_staged`] probes whether a given
//! prefix is present.
//!
//! [`crate::BwaIndex::load`] transparently attaches when a segment is
//! present and falls back to disk when not. Callers that want to require
//! the segment should probe with [`is_staged`] first — there is no
//! atomic "attach or fail" primitive (`stage` and `destroy` are racy
//! against an in-flight check).

use std::ffi::CString;
use std::path::Path;

use crate::error::{shim_err, Error, Result};

fn cprefix(prefix: impl AsRef<Path>) -> Result<CString> {
    let p = prefix.as_ref();
    let s = p
        .to_str()
        .ok_or_else(|| Error::InvalidInput("prefix must be valid UTF-8".into()))?;
    Ok(CString::new(s)?)
}

/// Returns true if an index segment keyed by `prefix`'s basename is staged.
pub fn is_staged(prefix: impl AsRef<Path>) -> Result<bool> {
    let c = cprefix(&prefix)?;
    let rc = unsafe { bwa_mem3_sys::bwa_shim_shm_test(c.as_ptr()) };
    match rc {
        1 => Ok(true),
        0 => Ok(false),
        _ => Err(shim_err("shm test")),
    }
}

/// Loads the index at `prefix` from disk, packs it, and stages it.
///
/// No-op if the prefix is already staged.
pub fn stage(prefix: impl AsRef<Path>) -> Result<()> {
    let c = cprefix(&prefix)?;
    let rc = unsafe { bwa_mem3_sys::bwa_shim_shm_stage(c.as_ptr()) };
    if rc < 0 {
        return Err(shim_err("shm stage"));
    }
    Ok(())
}

/// Drops every staged index segment and the control segment. Idempotent.
///
/// Safe against in-flight aligners on the same host: `shm_unlink` removes
/// the registry entry without invalidating existing mappings, so any
/// already-attached [`crate::BwaIndex`] handle keeps working until it is
/// dropped. Only *new* attaches see the segment as gone.
pub fn destroy() -> Result<()> {
    let rc = unsafe { bwa_mem3_sys::bwa_shim_shm_destroy() };
    if rc < 0 {
        return Err(shim_err("shm destroy"));
    }
    Ok(())
}

/// Prints staged segments to stdout, one `<basename>\t<bytes>` per line.
///
/// Matches `bwa shm -l`. The output is written via C's `printf` (upstream
/// implementation), which bypasses Rust's stdout buffering — `cargo test`
/// stdout capture will not see it. Mirroring the CLI semantics avoids
/// duplicating segment-format parsing in the shim; a structured iterator
/// is a follow-up.
pub fn list() -> Result<()> {
    let rc = unsafe { bwa_mem3_sys::bwa_shim_shm_list() };
    if rc < 0 {
        return Err(shim_err("shm list"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn null_prefix_is_invalid() {
        let p = "\0invalid";
        assert!(matches!(is_staged(p), Err(Error::InvalidInput(_))));
    }

    #[test]
    fn missing_prefix_is_not_staged() {
        // A prefix that almost certainly doesn't exist on any host.
        let p = "/nonexistent/bwa-mem3-rs-shm-test-prefix";
        assert!(!is_staged(p).expect("registry probe ok"));
    }
}
