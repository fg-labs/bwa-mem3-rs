//! Header parity: `bwa-rs mem`'s `@HD` and `@SQ` must match `bwa-mem3 mem`'s.
//!
//! bwa-mem3 v0.9.0 routes `@HD` through `compat_target_t` (`emit_hd` /
//! `hd_line`) after three hand-written literals in its own writers drifted
//! apart (fg-labs/bwa-mem3#288). This crate writes its own BAM header and had a
//! fourth literal — `@HD VN:1.6 SO:unknown` against the reference's
//! `@HD VN:1.5 SO:unsorted GO:query`. It now builds `@HD` and `@SQ` with
//! upstream's own `bwa_build_sam_hdr` (`BwaIndex::sam_header`), so this test
//! would fail again if the header writer ever went back to spelling them out.
//! (`cli_parity_alt.rs` pins the `AH:*` an ALT contig's `@SQ` carries.)
//!
//! `@PG` is deliberately not compared: ours legitimately names `bwa-rs`.
//! `read_sidecar`, the one compat switch this crate does not honor, is out of
//! scope here — the PhiX fixture has no `.hdr`/`.dict` sidecar, so both sides
//! generate `@SQ` from the index either way (see `bam_header_text`'s docs).
//!
//! Requires `bwa-mem3` + `samtools`; set `BWA_MEM3_BIN` if not on PATH.

mod common;
mod phix_seq;

use std::process::Command;

fn lines_with_prefix<'a>(lines: &'a [String], prefix: &str) -> Vec<&'a str> {
    lines
        .iter()
        .map(String::as_str)
        .filter(|l| l.starts_with(prefix))
        .collect()
}

/// `bwa-rs mem` and `bwa-mem3 mem` (each with `args`) on PhiX pairs against
/// `ref_fa`, as their header lines `(bwa-rs, CLI)`.
fn headers(
    bwa: &str,
    dir: &std::path::Path,
    ref_fa: &std::path::Path,
    args: &[&str],
) -> (Vec<String>, Vec<String>) {
    let pairs = common::simulate_pairs(phix_seq::PHIX_SEQ.as_bytes(), 20, 100, 250, 7);
    let r1_fq = dir.join("r1.fq");
    let r2_fq = dir.join("r2.fq");
    let r1: Vec<(String, Vec<u8>)> = pairs
        .iter()
        .map(|(n, a, _)| (n.clone(), a.clone()))
        .collect();
    let r2: Vec<(String, Vec<u8>)> = pairs
        .iter()
        .map(|(n, _, b)| (n.clone(), b.clone()))
        .collect();
    common::write_fastq(&r1_fq, &r1);
    common::write_fastq(&r2_fq, &r2);

    // bwa-rs mem -> BAM.
    let rs_bam = dir.join("rs.bam");
    let status = Command::new(common::cli_bin())
        .arg("mem")
        .args(args)
        .arg(ref_fa)
        .arg(&r1_fq)
        .arg(&r2_fq)
        .arg("-o")
        .arg(&rs_bam)
        .status()
        .expect("run bwa-rs mem");
    assert!(status.success(), "bwa-rs mem failed");

    // bwa-mem3 mem (reference); samtools reads its SAM or BAM alike.
    let cli_out = Command::new(bwa)
        .args(["mem", "-t", "1"])
        .args(args)
        .arg(ref_fa)
        .arg(&r1_fq)
        .arg(&r2_fq)
        .output()
        .expect("run bwa-mem3 mem");
    assert!(cli_out.status.success(), "bwa-mem3 mem failed");
    let cli_bam = dir.join("cli.bam");
    std::fs::write(&cli_bam, &cli_out.stdout).unwrap();

    (
        common::stored_bam_header_lines(&rs_bam),
        common::samtools_header_lines(&cli_bam),
    )
}

fn assert_hd_and_sq_match(rs: &[String], cli: &[String]) {
    let cli_hd = lines_with_prefix(cli, "@HD");
    assert_eq!(
        cli_hd.len(),
        1,
        "reference aligner emitted {} @HD lines; the comparison below assumes one",
        cli_hd.len()
    );
    assert_eq!(
        lines_with_prefix(rs, "@HD"),
        cli_hd,
        "@HD differs from the reference aligner — is bam_header_text still \
         building the header with BwaIndex::sam_header?"
    );
    // The SAM spec puts @HD first when present; the filtered comparisons
    // above would not notice it moving.
    assert_eq!(
        rs.first().map(String::as_str),
        Some(cli_hd[0]),
        "@HD is not the first stored header line"
    );
    let cli_sq = lines_with_prefix(cli, "@SQ");
    assert!(!cli_sq.is_empty(), "reference aligner emitted no @SQ");
    assert_eq!(
        lines_with_prefix(rs, "@SQ"),
        cli_sq,
        "@SQ differs from the reference aligner"
    );
}

#[test]
fn hd_and_sq_match_the_reference_aligner() {
    let Some(bwa) = common::require_vendored_bwa_mem3() else {
        return;
    };
    if !common::require_samtools() {
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let ref_fa = common::setup_phix_index(tmp.path(), &bwa, phix_seq::PHIX_SEQ);
    let (rs, cli) = headers(&bwa, tmp.path(), &ref_fa, &[]);
    assert_hd_and_sq_match(&rs, &cli);
}

/// Under `--meth` the header lists the ORIGINAL reference's contigs, not the
/// converted seed index's (`BwaIndex::sam_header` on a dual index).
#[test]
fn meth_hd_and_sq_match_the_reference_aligner() {
    let Some(bwa) = common::require_vendored_bwa_mem3() else {
        return;
    };
    if !common::require_samtools() {
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let ref_fa = common::setup_phix_meth_index(tmp.path(), &bwa, phix_seq::PHIX_SEQ);
    let (rs, cli) = headers(&bwa, tmp.path(), &ref_fa, &["--meth"]);
    assert_hd_and_sq_match(&rs, &cli);
}
