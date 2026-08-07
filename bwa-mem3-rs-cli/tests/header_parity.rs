//! Header parity: `bwa-rs mem`'s `@HD` and `@SQ` must match `bwa-mem3 mem`'s.
//!
//! bwa-mem3 v0.9.0 routes `@HD` through `compat_target_t` (`emit_hd` /
//! `hd_line`) after three hand-written literals in its own writers drifted
//! apart (fg-labs/bwa-mem3#288). This crate writes its own BAM header and had a
//! fourth literal — `@HD VN:1.6 SO:unknown` against the reference's
//! `@HD VN:1.5 SO:unsorted GO:query`. It now reads the line off the compat
//! target via `MemOpts::compat_hd_line`, so this test would fail again if the
//! header writer ever went back to spelling it out.
//!
//! `@PG` is deliberately not compared: ours legitimately names `bwa-rs`.
//! `read_sidecar`, the one compat switch this crate does not honor, is out of
//! scope here — the PhiX fixture has no `.hdr`/`.dict` sidecar, so both sides
//! generate `@SQ` from the index either way (see `write_bam_header`'s docs).
//!
//! Requires `bwa-mem3` + `samtools`; set `BWA_MEM3_BIN` if not on PATH.

mod common;
mod phix_seq;

use std::process::Command;

/// The `@`-prefixed header lines of a BAM file, via `samtools view -H`.
fn bam_header_lines(bam: &std::path::Path) -> Vec<String> {
    let out = Command::new("samtools")
        .args(["view", "-H"])
        .arg(bam)
        .output()
        .expect("run samtools view -H");
    assert!(out.status.success(), "samtools view -H failed");
    String::from_utf8_lossy(&out.stdout)
        .lines()
        .map(str::to_owned)
        .collect()
}

fn lines_with_prefix<'a>(lines: &'a [String], prefix: &str) -> Vec<&'a str> {
    lines
        .iter()
        .map(String::as_str)
        .filter(|l| l.starts_with(prefix))
        .collect()
}

#[test]
fn hd_and_sq_match_the_reference_aligner() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    if !common::require_samtools() {
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let ref_fa = common::setup_phix_index(dir, &bwa, phix_seq::PHIX_SEQ);

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

    let rs = bam_header_lines(&rs_bam);
    let cli = bam_header_lines(&cli_bam);

    let cli_hd = lines_with_prefix(&cli, "@HD");
    assert_eq!(
        cli_hd.len(),
        1,
        "reference aligner emitted {} @HD lines; the comparison below assumes one",
        cli_hd.len()
    );
    assert_eq!(
        lines_with_prefix(&rs, "@HD"),
        cli_hd,
        "@HD differs from the reference aligner — is write_bam_header still \
         reading the line off the compat target?"
    );

    assert_eq!(
        lines_with_prefix(&rs, "@SQ"),
        lines_with_prefix(&cli, "@SQ"),
        "@SQ differs from the reference aligner"
    );
}
