//! End-to-end bisulfite (`--meth`) regression + CLI-parity test.
//!
//! Builds a `bwa-mem3 index --meth` dual index over the embedded PhiX
//! reference, simulates deterministic paired reads, aligns them with both
//! `bwa-rs mem --meth` and `bwa-mem3 mem --meth`, and asserts:
//!   - the bwa-rs BAM is valid and every record carries the Bismark `XR:Z`
//!     tag, with `XG:Z` + `XM:Z` on every mapped record;
//!   - **every** record matches the CLI exactly on RNAME/POS/CIGAR and the
//!     `NM`/`MD`/`XG`/`XR`/`XM`/`XA` tags — secondaries and all.
//!
//! Requires `bwa-mem3` (for `index --meth` + the reference aligner) and
//! `samtools` on PATH; set `BWA_MEM3_BIN` if the CLI is not on PATH. Skips
//! gracefully otherwise.

mod common;
mod phix_seq;

use common::{record_key_fields, samtools_view};
use std::process::Command;

#[test]
fn meth_e2e_cli_parity() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    if !common::require_samtools() {
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();

    // Build the meth DUAL index (<ref>.* + <ref>.meth.*).
    let ref_fa = common::setup_phix_meth_index(dir, &bwa, phix_seq::PHIX_SEQ);

    // Simulate deterministic pairs and write FASTQs.
    let pairs = common::simulate_pairs(phix_seq::PHIX_SEQ.as_bytes(), 200, 100, 300, 7);
    let r1: Vec<_> = pairs
        .iter()
        .map(|(n, a, _)| (n.clone(), a.clone()))
        .collect();
    let r2: Vec<_> = pairs
        .iter()
        .map(|(n, _, b)| (n.clone(), b.clone()))
        .collect();
    let r1_fq = dir.join("r1.fq");
    let r2_fq = dir.join("r2.fq");
    common::write_fastq(&r1_fq, &r1);
    common::write_fastq(&r2_fq, &r2);

    // bwa-rs mem --meth -> BAM.
    let rs_bam = dir.join("rs.bam");
    let status = Command::new(common::cli_bin())
        .args(["mem", "--meth"])
        .arg(&ref_fa)
        .arg(&r1_fq)
        .arg(&r2_fq)
        .arg("-o")
        .arg(&rs_bam)
        .status()
        .expect("run bwa-rs mem --meth");
    assert!(status.success(), "bwa-rs mem --meth failed");

    // Validate structure and tag completeness on the bwa-rs output.
    let qc = Command::new("samtools")
        .args(["quickcheck"])
        .arg(&rs_bam)
        .status()
        .expect("run samtools quickcheck");
    assert!(
        qc.success(),
        "samtools quickcheck failed on bwa-rs meth BAM"
    );

    let rs_lines = samtools_view(&rs_bam);
    assert!(!rs_lines.is_empty(), "bwa-rs meth BAM has no records");
    for line in &rs_lines {
        let f: Vec<&str> = line.split('\t').collect();
        let flag: u32 = f[1].parse().unwrap();
        let has = |t: &str| f[11..].iter().any(|x| x.starts_with(t));
        assert!(has("XR:Z:"), "record missing XR:Z: {}", f[0]);
        if flag & 0x4 == 0 {
            assert!(has("XG:Z:"), "mapped record missing XG:Z: {}", f[0]);
            assert!(has("XM:Z:"), "mapped record missing XM:Z: {}", f[0]);
        }
    }

    // bwa-mem3 mem --meth -> BAM (reference implementation).
    let cli_bam = dir.join("cli.bam");
    let cli_out = Command::new(&bwa)
        .args(["mem", "--meth", "-t", "1"])
        .arg(&ref_fa)
        .arg(&r1_fq)
        .arg(&r2_fq)
        .output()
        .expect("run bwa-mem3 mem --meth");
    assert!(cli_out.status.success(), "bwa-mem3 mem --meth failed");
    std::fs::write(&cli_bam, &cli_out.stdout).unwrap();

    // Compare ALL records exactly (coords, MAPQ, CIGAR, NM, MD, XG, XR, XM, XA)
    // as sorted vectors so record multiplicity — not just the distinct key set —
    // must match the CLI, secondaries and XA:Z folding included. (This PhiX
    // fixture maps uniquely, so it does not surface XA:Z; the dedicated
    // repetitive-reference meth XA folding parity lives in meth_cli_parity_xa.)
    let mut cli_recs: Vec<String> = samtools_view(&cli_bam)
        .iter()
        .map(|l| record_key_fields(l))
        .collect();
    let mut rs_recs: Vec<String> = rs_lines.iter().map(|l| record_key_fields(l)).collect();
    cli_recs.sort();
    rs_recs.sort();

    assert!(!cli_recs.is_empty(), "CLI produced no records");
    assert_eq!(
        cli_recs, rs_recs,
        "bwa-rs meth records diverge from the bwa-mem3 CLI"
    );
}
