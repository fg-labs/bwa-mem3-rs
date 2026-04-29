//! Round-trip exercise of the shared-memory lifecycle.
//!
//! Opt-in — destroys global shm state and stages a multi-GB segment, so
//! only runs when both `BWA_MEM3_RS_TEST_REF` (the index prefix) and
//! `BWA_MEM3_RS_TEST_RUN_SHM=1` are set. CI does not enable it by
//! default; run locally with:
//!
//! ```text
//! BWA_MEM3_RS_TEST_REF=/path/to/ref.fa \
//! BWA_MEM3_RS_TEST_RUN_SHM=1 \
//!     cargo test -p bwa-mem3-rs --test shm_lifecycle
//! ```
//!
//! The test stages the index, loads from the staged segment (transparent
//! attach), runs a tiny `align_batch` against it, and then drops
//! everything. It does NOT preserve any pre-existing staged state.

use std::path::{Path, PathBuf};

use bwa_mem3_rs::{align_batch, shm, BwaIndex, MemOpts, ReadPair};

const INDEX_EXTS: &[&str] = &[".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"];

fn ref_prefix() -> Option<String> {
    let p = std::env::var("BWA_MEM3_RS_TEST_REF").ok()?;
    if std::path::Path::new(&format!("{p}.bwt.2bit.64")).exists() {
        Some(p)
    } else {
        eprintln!("skip: no bwa-mem3 index at BWA_MEM3_RS_TEST_REF={p}");
        None
    }
}

fn shm_enabled() -> bool {
    std::env::var("BWA_MEM3_RS_TEST_RUN_SHM").as_deref() == Ok("1")
}

/// macOS's POSIX shm caps segment names at 31 chars (`PSHMNAMLEN`).
/// The shm registry keys segments as `/bwaidx-<basename>`, so any prefix
/// whose basename is longer than ~22 chars hits the limit. Mirror the
/// index files into a tempdir with a short basename so the round-trip
/// works on macOS too.
fn short_alias_prefix(real_prefix: &str) -> (PathBuf, tempfile::TempDir) {
    let dir = tempfile::tempdir().expect("tempdir");
    let alias = dir.path().join("ref");
    for ext in INDEX_EXTS {
        let src = format!("{real_prefix}{ext}");
        let dst = format!("{}{}", alias.display(), ext);
        let src_path = Path::new(&src);
        if !src_path.exists() {
            continue;
        }
        // Absolutize: a relative `BWA_MEM3_RS_TEST_REF` would otherwise yield
        // a symlink target resolved relative to the tempdir (where the link
        // lives), not the test cwd, and dereferencing it would fail.
        let src_abs = src_path
            .canonicalize()
            .unwrap_or_else(|e| panic!("canonicalize {src}: {e}"));
        std::os::unix::fs::symlink(&src_abs, &dst)
            .unwrap_or_else(|e| panic!("symlink {} -> {dst}: {e}", src_abs.display()));
    }
    (alias, dir)
}

/// Calls `shm::destroy` on Drop so a panic mid-test doesn't leave the
/// host's shm registry populated for the next run.
struct ShmCleanupGuard;
impl Drop for ShmCleanupGuard {
    fn drop(&mut self) {
        let _ = shm::destroy();
    }
}

#[test]
fn stage_load_drop_roundtrip() {
    let Some(real_prefix) = ref_prefix() else {
        return;
    };
    if !shm_enabled() {
        eprintln!("skip: BWA_MEM3_RS_TEST_RUN_SHM not set");
        return;
    }

    // Symlink the index files under a short basename so the segment name
    // fits inside macOS's PSHMNAMLEN limit. Linux is fine with the real
    // basename but the alias works there too.
    let (alias, _td) = short_alias_prefix(&real_prefix);
    let prefix = alias.to_string_lossy().into_owned();

    let _guard = ShmCleanupGuard;

    // Clean slate.
    shm::destroy().expect("initial destroy");
    assert!(!shm::is_staged(&prefix).expect("test"));

    shm::stage(&prefix).expect("stage");
    assert!(shm::is_staged(&prefix).expect("test"));

    // BwaIndex::load attaches transparently when the segment is staged.
    let idx = BwaIndex::load(&prefix).expect("load");
    assert!(idx.n_contigs() > 0);

    let mut opts = MemOpts::new().expect("opts");
    opts.set_pe(true);
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let qual = vec![b'I'; seq.len()];
    let pairs = vec![ReadPair {
        name_r1: b"r",
        seq_r1: seq,
        qual_r1: Some(&qual),
        name_r2: b"r",
        seq_r2: seq,
        qual_r2: Some(&qual),
    }];
    let (aln, _pes) = align_batch(&idx, &opts, &pairs, None).expect("align");
    assert!(aln.len() >= 2);

    // Drop the index handle before destroying the segment so the shm-aware
    // free path runs in the expected order.
    drop(idx);
    shm::destroy().expect("destroy");
    assert!(!shm::is_staged(&prefix).expect("test"));
}
