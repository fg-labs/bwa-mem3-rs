//! Non-meth CLI-parity over pairs with a zero-length mate, through the legacy
//! `align_batch` path that `bwa-rs mem` drives.
//!
//! The bwa-mem3 CLI aligns an empty read as unmapped: it yields no seeds, so it
//! takes the same branch as an unmappable read, and its record carries SEQ/QUAL
//! `*`. When its mate maps, the half-mapped rewrite (`mem_aln2sam`) places the
//! empty record at the mate's coordinates. `three_phase_parity` pins the same
//! inputs, plus an empty single, on the three-phase API.
//!
//! Requires `bwa-mem3` + `samtools`; set `BWA_MEM3_BIN` if not on PATH. Skips
//! gracefully otherwise.

mod common;

use common::{random_dna, record_key_fields, samtools_view, Rng};
use std::process::Command;

/// FLAG, RNAME, POS, RNEXT, PNEXT, TLEN and SEQ of every record named `qname`,
/// sorted by FLAG. The mate columns are outside `record_key_fields`, and the
/// half-mapped rewrite writes exactly those.
fn placement(lines: &[String], qname: &str) -> Vec<[String; 7]> {
    let mut out: Vec<[String; 7]> = lines
        .iter()
        .filter(|l| l.split('\t').next() == Some(qname))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            [f[1], f[2], f[3], f[6], f[7], f[8], f[9]].map(str::to_string)
        })
        .collect();
    out.sort_by_key(|r| r[0].parse::<u32>().unwrap());
    out
}

#[test]
fn cli_parity_empty_mates() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    if !common::require_samtools() {
        return;
    }

    let mut rng = Rng(0xE397_0002);
    let seq = random_dna(&mut rng, 4000);

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let ref_fa = common::setup_phix_index(dir, &bwa, &seq);

    let (read_len, insert) = (100usize, 300usize);
    let mut r1_reads = Vec::new();
    let mut r2_reads = Vec::new();

    // Ordinary pairs, so the run has an insert-size distribution to estimate.
    for i in 0..40 {
        let start = 200 + i * 40;
        let r1 = seq.as_bytes()[start..start + read_len].to_vec();
        let r2 = common::revcomp(&seq.as_bytes()[start + insert - read_len..start + insert]);
        r1_reads.push((format!("ok{i}"), r1));
        r2_reads.push((format!("ok{i}"), r2));
    }

    // The fixtures under test: an empty R2, an empty R1 (mate on the reverse
    // strand, so the strand copy is exercised too), and both mates empty.
    r1_reads.push((
        "e_r2".to_string(),
        seq.as_bytes()[1500..1500 + read_len].to_vec(),
    ));
    r2_reads.push(("e_r2".to_string(), Vec::new()));
    r1_reads.push(("e_r1".to_string(), Vec::new()));
    r2_reads.push((
        "e_r1".to_string(),
        common::revcomp(&seq.as_bytes()[2500..2500 + read_len]),
    ));
    r1_reads.push(("e_both".to_string(), Vec::new()));
    r2_reads.push(("e_both".to_string(), Vec::new()));

    let r1_fq = dir.join("r1.fq");
    let r2_fq = dir.join("r2.fq");
    common::write_fastq(&r1_fq, &r1_reads);
    common::write_fastq(&r2_fq, &r2_reads);

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
    assert!(status.success(), "bwa-rs mem failed on empty mates");

    let cli_out = Command::new(&bwa)
        .args(["mem", "-t", "1"])
        .arg(&ref_fa)
        .arg(&r1_fq)
        .arg(&r2_fq)
        .output()
        .expect("run bwa-mem3 mem");
    assert!(
        cli_out.status.success(),
        "bwa-mem3 mem failed: {}",
        String::from_utf8_lossy(&cli_out.stderr)
    );
    let cli_bam = dir.join("cli.bam");
    std::fs::write(&cli_bam, &cli_out.stdout).unwrap();

    let rs_lines = samtools_view(&rs_bam);
    let cli_lines = samtools_view(&cli_bam);

    // Fixture validity, asserted against the reference aligner: every empty
    // mate is an unmapped record with SEQ `*`.
    for qname in ["e_r2", "e_r1", "e_both"] {
        let recs = placement(&cli_lines, qname);
        assert_eq!(recs.len(), 2, "expected two {qname} records from the CLI");
        let empties: Vec<_> = recs.iter().filter(|r| r[6] == "*").collect();
        assert!(
            !empties.is_empty()
                && empties
                    .iter()
                    .all(|r| r[0].parse::<u32>().unwrap() & 0x4 != 0),
            "{qname}: the CLI must emit each empty mate as unmapped with SEQ `*`: {recs:?}"
        );
        assert_eq!(
            placement(&rs_lines, qname),
            recs,
            "{qname}: FLAG/placement/mate columns/SEQ diverge from the bwa-mem3 CLI"
        );
    }

    // Sorted vectors (not sets): preserve record multiplicity so a duplicated
    // or missing record diverges even when the distinct key set matches.
    let mut cli_recs: Vec<String> = cli_lines.iter().map(|l| record_key_fields(l)).collect();
    let mut rs_recs: Vec<String> = rs_lines.iter().map(|l| record_key_fields(l)).collect();
    cli_recs.sort();
    rs_recs.sort();
    assert_eq!(
        cli_recs, rs_recs,
        "bwa-rs records diverge from the bwa-mem3 CLI on pairs with an empty mate"
    );
}
