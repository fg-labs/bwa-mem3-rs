//! Non-meth CLI-parity over a MULTI-MAPPING reference — pins the shim's
//! `XA:Z` folding (`pair_and_emit` mirrors upstream `mem_reg2sam`).
//!
//! The reference embeds one 160 bp motif twice, so reads drawn from it map to
//! both copies. Upstream folds the second hit into the primary's `XA:Z` tag
//! rather than emitting a second record; the shim must do the same. Earlier the
//! shim emitted a record per hit, which only the (count-based) `cli_parity`
//! library test caught — nothing verified the `XA:Z` *content*. This asserts
//! every record (and its `XA:Z`) matches `bwa-mem3 mem` exactly, and that at
//! least one record actually carries a non-empty `XA:Z` (so the fold is
//! genuinely exercised).
//!
//! Requires `bwa-mem3` + `samtools`; set `BWA_MEM3_BIN` if not on PATH. Skips
//! gracefully otherwise.

mod common;

use common::{random_dna, record_key_fields, samtools_view, Rng};
use std::process::Command;

#[test]
fn cli_parity_xa_folding() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    if !common::require_samtools() {
        return;
    }

    // Reference: bg | MOTIF | bg | MOTIF | bg  (the 160 bp motif appears twice).
    let mut rng = Rng(0x5eed_1234);
    let motif = random_dna(&mut rng, 160);
    let seq = format!(
        "{}{motif}{}{motif}{}",
        random_dna(&mut rng, 1200),
        random_dna(&mut rng, 1200),
        random_dna(&mut rng, 1200),
    );
    let motif1_start = 1200usize;

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let ref_fa = common::setup_phix_index(dir, &bwa, &seq);

    // Pairs whose R1 sits inside motif copy 1 (so R1 multi-maps → XA), with R2
    // the revcomp of a window `insert` downstream. read_len <= motif so R1 is
    // fully inside the repeat.
    let (read_len, insert) = (100usize, 260usize);
    let mut r1_reads = Vec::new();
    let mut r2_reads = Vec::new();
    for i in 0..60 {
        let start = motif1_start + (i % (160 - read_len + 1));
        let r1 = seq.as_bytes()[start..start + read_len].to_vec();
        let r2 = common::revcomp(&seq.as_bytes()[start + insert - read_len..start + insert]);
        r1_reads.push((format!("x{i}"), r1));
        r2_reads.push((format!("x{i}"), r2));
    }

    // One pair drawn from `rng` but never sliced out of `seq` — i.e. it isn't
    // a substring of the reference the index was built from, unlike every
    // other read above. This pins the shim's unmapped path (`fill_unmapped`
    // in `bwa_shim_align.cpp`): the CLI never emits `HN:i` on an unmapped
    // record, but the shim once did (HN:i:0, since fill_unmapped's calloc'd
    // struct leaves `HN` at 0 rather than upstream's -1 sentinel) until that
    // was fixed. Both mates unmapped, each independently, so this exercises
    // `fill_unmapped` on both sides of the pair. Whether it actually fails to
    // map is verified below, not assumed.
    r1_reads.push((
        "unmapped0".to_string(),
        random_dna(&mut rng, read_len).into_bytes(),
    ));
    r2_reads.push((
        "unmapped0".to_string(),
        random_dna(&mut rng, read_len).into_bytes(),
    ));

    let r1_fq = dir.join("r1.fq");
    let r2_fq = dir.join("r2.fq");
    common::write_fastq(&r1_fq, &r1_reads);
    common::write_fastq(&r2_fq, &r2_reads);

    // bwa-rs mem -> BAM.
    let rs_bam = dir.join("rs.bam");
    let status = Command::new(common::cli_bin())
        .args(["mem"])
        .arg(&ref_fa)
        .arg(&r1_fq)
        .arg(&r2_fq)
        .arg("-o")
        .arg(&rs_bam)
        .status()
        .expect("run bwa-rs mem");
    assert!(status.success(), "bwa-rs mem failed");

    // bwa-mem3 mem -> BAM (reference).
    let cli_out = Command::new(&bwa)
        .args(["mem", "-t", "1"])
        .arg(&ref_fa)
        .arg(&r1_fq)
        .arg(&r2_fq)
        .output()
        .expect("run bwa-mem3 mem");
    assert!(cli_out.status.success(), "bwa-mem3 mem failed");
    let cli_bam = dir.join("cli.bam");
    std::fs::write(&cli_bam, &cli_out.stdout).unwrap();

    let rs_lines = samtools_view(&rs_bam);
    let cli_lines = samtools_view(&cli_bam);

    // The fold must actually be exercised: at least one record carries XA:Z.
    let xa_count = cli_lines.iter().filter(|l| l.contains("\tXA:Z:")).count();
    assert!(
        xa_count > 0,
        "reference produced no XA:Z tags — test would not exercise XA folding"
    );

    // The unmapped-path pin must actually be exercised: verify (not assume)
    // that the reference aligner emits both mates of "unmapped0" and reports
    // each as unmapped (SAM flag bit 0x4). Assert the records *exist* first —
    // counting only the mapped ones would pass vacuously if the pair never
    // appeared at all (a renamed or dropped fixture), which is the same
    // failure mode this block exists to rule out. A short random read can
    // also align to a small reference by chance, which would silently unpin
    // this case.
    let unmapped0_flags: Vec<u32> = cli_lines
        .iter()
        .filter(|l| l.split('\t').next() == Some("unmapped0"))
        .map(|l| l.split('\t').nth(1).unwrap().parse::<u32>().unwrap())
        .collect();
    assert_eq!(
        unmapped0_flags.len(),
        2,
        "expected exactly two \"unmapped0\" records (one per mate) from the reference \
         aligner, got flags {unmapped0_flags:?} — the unmapped-record pin below would \
         be vacuous"
    );
    assert!(
        unmapped0_flags.iter().all(|f| f & 0x4 != 0),
        "expected \"unmapped0\" to be unmapped by the reference aligner on both mates \
         (flags {unmapped0_flags:?}) — the HN:i-on-unmapped-record regression would go \
         unpinned"
    );
    assert_eq!(
        (
            unmapped0_flags.iter().filter(|&&f| f & 0x40 != 0).count(),
            unmapped0_flags.iter().filter(|&&f| f & 0x80 != 0).count(),
        ),
        (1, 1),
        "expected one first-in-pair (0x40) and one second-in-pair (0x80) \"unmapped0\" \
         record, got flags {unmapped0_flags:?} — both mates must exercise fill_unmapped"
    );

    // Sorted vectors (not sets): preserve record multiplicity so a duplicated
    // or missing record diverges even when the distinct key set matches.
    let mut cli_recs: Vec<String> = cli_lines.iter().map(|l| record_key_fields(l)).collect();
    let mut rs_recs: Vec<String> = rs_lines.iter().map(|l| record_key_fields(l)).collect();
    cli_recs.sort();
    rs_recs.sort();
    assert_eq!(
        cli_recs, rs_recs,
        "bwa-rs records (incl. XA:Z folding) diverge from the bwa-mem3 CLI"
    );
}
