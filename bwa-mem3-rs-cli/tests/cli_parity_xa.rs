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
    let Some(bwa) = common::find_bwa_mem3() else {
        eprintln!("skip: bwa-mem3 not on PATH (set BWA_MEM3_BIN)");
        return;
    };
    if !common::have_samtools() {
        eprintln!("skip: samtools not on PATH");
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
