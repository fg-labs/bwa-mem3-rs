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

use std::collections::BTreeSet;
use std::process::Command;

/// One record reduced to the fields that must match the CLI, keyed for sorting.
fn record_key_fields(sam_line: &str) -> Option<String> {
    let f: Vec<&str> = sam_line.split('\t').collect();
    if f.len() < 11 {
        return None;
    }
    let (mut nm, mut md, mut xg, mut xr, mut xm, mut xa) = ("", "", "", "", "", "");
    for tag in &f[11..] {
        if let Some(v) = tag.strip_prefix("NM:") {
            nm = v;
        } else if let Some(v) = tag.strip_prefix("MD:") {
            md = v;
        } else if let Some(v) = tag.strip_prefix("XG:") {
            xg = v;
        } else if let Some(v) = tag.strip_prefix("XR:") {
            xr = v;
        } else if let Some(v) = tag.strip_prefix("XM:") {
            xm = v;
        } else if let Some(v) = tag.strip_prefix("XA:") {
            xa = v;
        }
    }
    // qname, flag, rname, pos, mapq, cigar + meth/edit/alt tags.
    Some(format!(
        "{}\t{}\t{}\t{}\t{}\t{}\tNM:{nm}\tMD:{md}\tXG:{xg}\tXR:{xr}\tXM:{xm}\tXA:{xa}",
        f[0], f[1], f[2], f[3], f[4], f[5]
    ))
}

/// `samtools view` a BAM file and return its records as SAM text lines.
fn samtools_view(bam: &std::path::Path) -> Vec<String> {
    let out = Command::new("samtools")
        .args(["view"])
        .arg(bam)
        .output()
        .expect("run samtools view");
    assert!(out.status.success(), "samtools view failed");
    String::from_utf8_lossy(&out.stdout)
        .lines()
        .map(str::to_owned)
        .collect()
}

#[test]
fn meth_e2e_cli_parity() {
    let Some(bwa) = common::find_bwa_mem3() else {
        eprintln!("skip: bwa-mem3 not on PATH (set BWA_MEM3_BIN)");
        return;
    };
    if !common::have_samtools() {
        eprintln!("skip: samtools not on PATH");
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();

    // Write the PhiX FASTA and build the meth DUAL index (<ref>.* + <ref>.meth.*).
    let ref_fa = dir.join("phix.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&ref_fa).unwrap();
        writeln!(f, ">phix").unwrap();
        for chunk in phix_seq::PHIX_SEQ.as_bytes().chunks(72) {
            f.write_all(chunk).unwrap();
            writeln!(f).unwrap();
        }
    }
    let status = Command::new(&bwa)
        .args(["index", "--meth"])
        .arg(&ref_fa)
        .status()
        .expect("run bwa-mem3 index --meth");
    assert!(status.success(), "bwa-mem3 index --meth failed");
    assert!(
        ref_fa.with_extension("fa.meth.bwt.2bit.64").exists(),
        "index --meth did not produce the converted seed index"
    );

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
        .args(["-o".as_ref(), rs_bam.as_os_str()])
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

    // Compare ALL records exactly (coords, MAPQ, CIGAR, NM, MD, XG, XR, XM, XA).
    let cli_recs: BTreeSet<String> = samtools_view(&cli_bam)
        .iter()
        .filter_map(|l| record_key_fields(l))
        .collect();
    let rs_recs: BTreeSet<String> = rs_lines
        .iter()
        .filter_map(|l| record_key_fields(l))
        .collect();

    assert!(!cli_recs.is_empty(), "CLI produced no records");
    assert_eq!(
        cli_recs, rs_recs,
        "bwa-rs meth records diverge from the bwa-mem3 CLI"
    );
}
