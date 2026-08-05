//! End-to-end regression test: run the `bwa-rs` CLI against a PhiX reference
//! with simulated reads, verify the output parses as valid BAM with the
//! expected record count.
//!
//! Requires `bwa-mem3` (for `index`) and `samtools` on PATH. Set
//! `BWA_MEM3_BIN` to point at the CLI if it's not on PATH. The test skips
//! gracefully when either is missing; CI installs them via conda.

mod common;
mod phix_seq;

use std::process::{Command, Stdio};

#[test]
fn phix_e2e_roundtrip() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    if !common::require_samtools() {
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();

    let ref_fa = common::setup_phix_index(dir, &bwa, phix_seq::PHIX_SEQ);

    // Simulate 1000 paired reads.
    let n_pairs = 1000;
    let read_len = 100;
    let insert = 250;
    let pairs =
        common::simulate_pairs(phix_seq::PHIX_SEQ.as_bytes(), n_pairs, read_len, insert, 42);
    let r1_fq = dir.join("r1.fq");
    let r2_fq = dir.join("r2.fq");
    let r1_reads: Vec<_> = pairs
        .iter()
        .map(|(n, s, _)| (n.clone(), s.clone()))
        .collect();
    let r2_reads: Vec<_> = pairs
        .iter()
        .map(|(n, _, s)| (n.clone(), s.clone()))
        .collect();
    common::write_fastq(&r1_fq, &r1_reads);
    common::write_fastq(&r2_fq, &r2_reads);

    let out_bam = dir.join("out.bam");
    let status = Command::new(common::cli_bin())
        .args(["mem"])
        .arg(ref_fa.as_os_str())
        .arg(&r1_fq)
        .arg(&r2_fq)
        .arg("-o")
        .arg(&out_bam)
        .status()
        .expect("run bwa-rs");
    assert!(status.success(), "bwa-rs mem failed");
    assert!(out_bam.exists(), "bwa-rs did not produce output BAM");

    let qc = Command::new("samtools")
        .args(["quickcheck"])
        .arg(&out_bam)
        .output()
        .unwrap();
    assert!(
        qc.status.success(),
        "samtools quickcheck failed: {}",
        String::from_utf8_lossy(&qc.stderr)
    );

    let view = Command::new("samtools")
        .args(["view", "-c"])
        .arg(&out_bam)
        .stdout(Stdio::piped())
        .output()
        .unwrap();
    assert!(view.status.success(), "samtools view -c failed");
    let n_recs: usize = String::from_utf8_lossy(&view.stdout)
        .trim()
        .parse()
        .expect("parse record count");
    let expected_min = 2 * n_pairs;
    assert!(
        n_recs >= expected_min,
        "too few records: got {n_recs}, expected at least {expected_min}"
    );
    eprintln!("e2e: aligned {n_pairs} pairs → {n_recs} records");
}
