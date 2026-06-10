//! `BwaIndex::load` error paths (issue #19).
//!
//! bwa-mem3's loader historically called `exit(1)` when an on-disk index file
//! was missing, aborting the host process before the safe wrapper could
//! return. These tests pin the contract that a missing or incomplete index
//! surfaces a recoverable [`Err`] instead. They need no prebuilt index and
//! therefore always run.

use bwa_mem3_rs::BwaIndex;
use std::fs;

#[test]
fn load_missing_index_returns_err() {
    let dir = tempfile::tempdir().expect("tempdir");
    let prefix = dir.path().join("absent");

    let msg = match BwaIndex::load(&prefix) {
        Ok(_) => panic!("a missing index must return Err, not abort"),
        Err(e) => e.to_string(),
    };
    assert!(
        msg.contains("absent"),
        "error should name the failing prefix, got: {msg}"
    );
}

#[test]
fn load_partial_index_names_the_missing_file() {
    // Every required companion present except `.ann`: the loader would
    // otherwise read the FM-index, then exit() partway through.
    let dir = tempfile::tempdir().expect("tempdir");
    let prefix = dir.path().join("partial");
    for suffix in [".bwt.2bit.64", ".amb", ".pac"] {
        fs::write(format!("{}{suffix}", prefix.display()), b"stub").expect("write stub");
    }

    let msg = match BwaIndex::load(&prefix) {
        Ok(_) => panic!("an incomplete index must return Err"),
        Err(e) => e.to_string(),
    };
    assert!(
        msg.contains(".ann"),
        "error should name the missing companion file, got: {msg}"
    );
}

#[test]
fn load_rejects_non_regular_companion() {
    // A directory sitting where a regular index file is expected passes a
    // bare `exists()` check but would still trip the vendored loader's
    // `exit(1)`. The pre-flight must treat non-regular paths as unusable.
    let dir = tempfile::tempdir().expect("tempdir");
    let prefix = dir.path().join("nonreg");
    for suffix in [".amb", ".ann", ".pac"] {
        fs::write(format!("{}{suffix}", prefix.display()), b"stub").expect("write stub");
    }
    fs::create_dir(format!("{}.bwt.2bit.64", prefix.display())).expect("mkdir companion");

    let msg = match BwaIndex::load(&prefix) {
        Ok(_) => panic!("a non-regular index file must return Err"),
        Err(e) => e.to_string(),
    };
    assert!(
        msg.contains(".bwt.2bit.64"),
        "error should name the non-regular file, got: {msg}"
    );
}
