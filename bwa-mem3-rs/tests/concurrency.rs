//! Verify concurrent `align_batch` calls on the same `BwaIndex` produce
//! byte-identical records whether run serially or in parallel.
//!
//! This is the regression guard for the pac-fetch reference path: extension
//! reconstructs each reference window from the resident `.pac` into a per-thread
//! `thread_local` scratch buffer that the next fetch on that thread overwrites
//! (see `bns_get_seq_v2` in the vendored `bntseq.cpp`). A data race there could
//! corrupt POS/CIGAR/SEQ while leaving the record *count* intact, so this test
//! asserts record IDENTITY between the serial and parallel runs, not just count.
//!
//! It builds a small PhiX index at test time (the pattern used by
//! `three_phase.rs`), so it runs in CI wherever `bwa-mem3` is on PATH instead of
//! being gated on a full index CI never provides. Reads are drawn from the PhiX
//! sequence itself so they actually map and exercise extension — random reads
//! would be all-unmapped and never touch the pac-fetch scratch. A missing
//! `bwa-mem3` skips the test, unless `BWA_MEM3_RS_REQUIRE_TOOLS` is set, in
//! which case the absence is a hard failure (the same skip-vs-panic convention
//! the sys/cli crates use).

use bwa_mem3_rs::{align_batch, AlignmentBatch, BwaIndex, MemOpts, ReadPair};
use rayon::prelude::*;
use std::io::Write;
use std::path::Path;
use std::process::Command;
use std::sync::Arc;

#[path = "../../bwa-mem3-rs-cli/tests/phix_seq.rs"]
mod phix_seq;

/// Locate `bwa-mem3` via `BWA_MEM3_BIN` or `PATH`. Returns `None` (with a skip
/// message) when it is absent, unless `BWA_MEM3_RS_REQUIRE_TOOLS` is set, which
/// turns the absence into a hard failure.
fn find_bwa_mem3() -> Option<String> {
    if let Ok(p) = std::env::var("BWA_MEM3_BIN") {
        if Path::new(&p).exists() {
            return Some(p);
        }
    }
    let out = Command::new("which").arg("bwa-mem3").output().ok();
    let found = out.and_then(|o| {
        o.status
            .success()
            .then(|| String::from_utf8_lossy(&o.stdout).trim().to_string())
            .filter(|p| !p.is_empty())
    });
    if found.is_none() {
        assert!(
            std::env::var_os("BWA_MEM3_RS_REQUIRE_TOOLS").is_none(),
            "BWA_MEM3_RS_REQUIRE_TOOLS is set but bwa-mem3 was not found; \
             set BWA_MEM3_BIN or install it on PATH"
        );
        eprintln!("skip: bwa-mem3 not on PATH (set BWA_MEM3_BIN)");
    }
    found
}

/// Build a PhiX index in a temp dir and load it. The `TempDir` is returned so
/// the caller keeps the on-disk index files alive for the index's lifetime.
fn build_phix_index() -> Option<(tempfile::TempDir, Arc<BwaIndex>)> {
    let bwa = find_bwa_mem3()?;
    let dir = tempfile::tempdir().expect("tempdir");
    let fa = dir.path().join("phix.fa");
    let mut f = std::fs::File::create(&fa).unwrap();
    writeln!(f, ">phix").unwrap();
    for chunk in phix_seq::PHIX_SEQ.as_bytes().chunks(72) {
        f.write_all(chunk).unwrap();
        writeln!(f).unwrap();
    }
    drop(f);
    let status = Command::new(&bwa)
        .arg("index")
        .arg(&fa)
        .status()
        .expect("run bwa-mem3 index");
    assert!(status.success(), "bwa-mem3 index failed");
    let idx = Arc::new(BwaIndex::load(&fa).expect("load PhiX index"));
    Some((dir, idx))
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => *other,
        })
        .collect()
}

/// One simulated read pair: (name, R1 bases, R2 bases).
type Pair = (String, Vec<u8>, Vec<u8>);

/// A record is mapped when the BAM FLAG's `0x4` (unmapped) bit is clear. The
/// record is a packed BAM entry prefixed with a 4-byte block size, so the FLAG
/// (`u16` at core offset 14) sits at byte offset 18. See `align_smoke.rs` for
/// the full field layout.
fn is_mapped(rec: &[u8]) -> bool {
    rec.len() >= 20 && (u16::from_le_bytes([rec[18], rec[19]]) & 0x4) == 0
}

fn canonicalize(batch: &AlignmentBatch) -> Vec<(usize, Vec<u8>)> {
    let mut v: Vec<_> = batch
        .iter()
        .map(|r| (r.pair_idx, r.bytes.to_vec()))
        .collect();
    v.sort();
    v
}

#[test]
fn par_iter_matches_serial() {
    let Some((_dir, idx)) = build_phix_index() else {
        return;
    };
    let mut opts = MemOpts::new().expect("opts");
    opts.set_pe(true);

    // Reads drawn from PhiX itself so they map and trigger extension (which is
    // what fetches reference windows into the per-thread pac-fetch scratch). Each
    // pair is a proper FR pair: R1 forward at `off`, R2 the reverse complement of
    // a window `insert` bases downstream. Offsets are spread across the genome so
    // many distinct windows are fetched, maximising contention on the scratch.
    let genome = phix_seq::PHIX_SEQ.as_bytes();
    const READ_LEN: usize = 100;
    const INSERT: usize = 300;
    const N_BATCHES: usize = 32;
    const PAIRS_PER_BATCH: usize = 16;
    let qual = vec![b'I'; READ_LEN];
    let max_off = genome.len() - INSERT;
    let n_pairs = N_BATCHES * PAIRS_PER_BATCH;
    let step = (max_off / n_pairs).max(1);

    // Owned read data: (name, r1_seq, r2_seq), grouped into batches.
    let batches: Vec<Vec<Pair>> = (0..N_BATCHES)
        .map(|batch_i| {
            (0..PAIRS_PER_BATCH)
                .map(|pair_i| {
                    let idx = batch_i * PAIRS_PER_BATCH + pair_i;
                    let off = (idx * step) % max_off;
                    let r1 = genome[off..off + READ_LEN].to_vec();
                    let r2 = revcomp(&genome[off + INSERT - READ_LEN..off + INSERT]);
                    (format!("b{batch_i}p{pair_i}"), r1, r2)
                })
                .collect()
        })
        .collect();

    let do_align = |batch: &Vec<Pair>| -> Vec<(usize, Vec<u8>)> {
        let pairs: Vec<_> = batch
            .iter()
            .map(|(name, r1, r2)| ReadPair {
                name_r1: name.as_bytes(),
                seq_r1: r1,
                qual_r1: Some(&qual),
                name_r2: name.as_bytes(),
                seq_r2: r2,
                qual_r2: Some(&qual),
            })
            .collect();
        let (aln, _) = align_batch(&idx, &opts, &pairs, None).expect("align");
        canonicalize(&aln)
    };

    let serial: Vec<_> = batches.iter().map(do_align).collect();
    let parallel: Vec<_> = batches.par_iter().map(do_align).collect();

    // align_batch is a pure function of (index, opts, pairs): identical inputs
    // must yield identical records regardless of the thread they run on. A
    // difference here means the parallel run corrupted shared per-thread state —
    // exactly the pac-fetch scratch race this test exists to catch. Asserting
    // full byte identity (not just count) is what makes a silent POS/CIGAR/SEQ
    // corruption visible.
    assert_eq!(serial.len(), parallel.len());
    let mut total = 0usize;
    let mut mapped = 0usize;
    for (i, (s, p)) in serial.iter().zip(parallel.iter()).enumerate() {
        assert_eq!(s, p, "batch {i} differs between serial and parallel runs");
        total += s.len();
        mapped += s.iter().filter(|(_, rec)| is_mapped(rec)).count();
    }
    // The identity assertion is only meaningful if the reads actually mapped and
    // extended — that is what fetches reference windows into the pac-fetch
    // scratch. Unmapped reads still emit records, so assert on the MAPPED count
    // (FLAG & 0x4 clear), not the total, or the guard is vacuous.
    assert!(
        mapped >= n_pairs,
        "expected reads to map and exercise extension, but only {mapped}/{total} records mapped"
    );
    eprintln!(
        "serial == parallel across {} batches, {mapped}/{total} records mapped",
        serial.len()
    );
}
