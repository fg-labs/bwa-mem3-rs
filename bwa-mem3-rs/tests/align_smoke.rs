//! End-to-end smoke test: load a real bwa-mem3 index and run `align_batch`.
//!
//! Skipped unless `BWA_MEM3_RS_TEST_REF` points at a prebuilt bwa-mem3 index
//! (i.e. the prefix such that `<prefix>.bwt.2bit.64` exists).
//!
//! NOTE: the current shim emits SAM lines (not packed BAM) in
//! `Record::bytes`; this test verifies structure, not byte format.

use bwa_mem3_rs::{align_batch, BwaIndex, MemOpts, ReadPair};

fn ref_prefix() -> Option<String> {
    let p = std::env::var("BWA_MEM3_RS_TEST_REF").ok()?;
    if std::path::Path::new(&format!("{p}.bwt.2bit.64")).exists() {
        Some(p)
    } else {
        eprintln!("skip: no bwa-mem3 index at BWA_MEM3_RS_TEST_REF={p}");
        None
    }
}

#[test]
fn align_batch_returns_records() {
    let Some(prefix) = ref_prefix() else {
        return;
    };
    let idx = BwaIndex::load(&prefix).expect("load index");
    assert!(idx.n_contigs() > 0);

    let mut opts = MemOpts::new().expect("opts");
    opts.set_pe(true);

    // Two tiny pairs of random-looking synthetic reads. Expected behavior:
    // align_batch returns something (possibly all-unmapped); no crash.
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let qual = vec![b'I'; seq.len()];
    let pairs = vec![ReadPair {
        name_r1: b"read1",
        seq_r1: seq,
        qual_r1: Some(&qual),
        name_r2: b"read1",
        seq_r2: seq,
        qual_r2: Some(&qual),
    }];

    let (aln, _pes) = align_batch(&idx, &opts, &pairs, None).expect("align");
    eprintln!("aligned {} records over {} pairs", aln.len(), pairs.len());
    assert!(
        aln.len() >= 2,
        "expected at least 2 records, got {}",
        aln.len()
    );
    for r in aln.iter() {
        // Packed BAM record: [u32 le block_size][block_size bytes of record data].
        assert!(
            r.bytes.len() >= 4 + 32,
            "record too short: {} bytes",
            r.bytes.len()
        );
        let block_size = u32::from_le_bytes([r.bytes[0], r.bytes[1], r.bytes[2], r.bytes[3]]);
        assert_eq!(
            block_size as usize + 4,
            r.bytes.len(),
            "block_size does not match buffer length"
        );

        // Parse core BAM fields directly from packed bytes (per BAM spec 4.2).
        let core = &r.bytes[4..];
        let ref_id = i32::from_le_bytes([core[0], core[1], core[2], core[3]]);
        let pos = i32::from_le_bytes([core[4], core[5], core[6], core[7]]);
        let l_read_name = core[8] as usize;
        let mapq = core[9];
        let n_cigar = u16::from_le_bytes([core[12], core[13]]) as usize;
        let flag = u16::from_le_bytes([core[14], core[15]]);
        let l_seq = i32::from_le_bytes([core[16], core[17], core[18], core[19]]) as usize;

        // read_name starts at offset 32, NUL-terminated.
        let name_bytes = &core[32..32 + l_read_name - 1];
        let name = std::str::from_utf8(name_bytes).expect("utf8 name");

        eprintln!(
            "  pair={} len={} ref_id={} pos={} mapq={} n_cigar={} flag=0x{:x} l_seq={} name={}",
            r.pair_idx,
            r.bytes.len(),
            ref_id,
            pos,
            mapq,
            n_cigar,
            flag,
            l_seq,
            name,
        );

        assert_eq!(name, "read1");
        // l_seq is original seq length (60) for primary, less for hard-clipped
        // supplementary. Either way it must be positive.
        assert!(l_seq > 0 && l_seq <= 60);
        assert_eq!(flag & 0x1, 0x1); // paired
    }
}

/// The footgun guard: enabling `--meth` on opts but aligning against a plain
/// (non-meth) index must fail fast rather than silently emit garbage.
#[test]
fn meth_opts_against_non_meth_index_errors() {
    let Some(prefix) = ref_prefix() else {
        return;
    };
    let idx = BwaIndex::load(&prefix).expect("load index");
    assert!(!idx.is_meth(), "test ref must be a plain (non-meth) index");

    let mut opts = MemOpts::new().expect("opts");
    opts.set_pe(true).set_meth(true);

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

    // (AlignmentBatch isn't Debug, so avoid expect_err's Debug bound.)
    let Err(err) = align_batch(&idx, &opts, &pairs, None) else {
        panic!("meth opts + non-meth index must error");
    };
    // Matched on the variant, not just a substring: `Error` has six variants,
    // and several carry free-form strings that could incidentally contain
    // "meth" (an `IndexLoad` whose path does, for instance), so a bare
    // `contains` would not distinguish the guard firing from an unrelated
    // failure on the same call.
    let bwa_mem3_rs::Error::InvalidInput(msg) = &err else {
        panic!("expected Error::InvalidInput from the meth-mismatch guard, got: {err:?}");
    };
    assert!(
        msg.contains("--meth enabled but the index is not a meth dual index"),
        "expected the opts-meth/plain-index message, got: {msg}"
    );
}
