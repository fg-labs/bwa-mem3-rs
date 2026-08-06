//! Non-meth CLI-parity over a pair where exactly ONE mate maps.
//!
//! Upstream `mem_aln2sam` (`bwamem.cpp:2457-2462`) rewrites a half-mapped pair
//! before writing it: the unmapped record takes its mate's `rid`/`pos`/`is_rev`
//! and drops its CIGAR, and the mapped record's mate fields take its own. Both
//! copies happen on LOCAL copies (`ptmp`/`mtmp`), so formatting one record
//! never perturbs the other. The `0x10`/`0x20` strand bits are then derived
//! from the possibly-copied values, which is why an unmapped record placed at
//! its mate's coordinates inherits that mate's strand.
//!
//! Carrying the mate's coordinates is what keeps the pair adjacent under
//! coordinate sort, so this is not only a byte-parity concern — most
//! downstream tools depend on it.
//!
//! The existing fixtures cannot catch this: every simulated read is an exact
//! substring of the reference it was indexed from, and `unmapped0` in
//! `cli_parity_xa` has *both* mates unmapped, which takes neither branch.
//!
//! Requires `bwa-mem3` + `samtools`; set `BWA_MEM3_BIN` if not on PATH. Skips
//! gracefully otherwise.

mod common;

use common::{random_dna, record_key_fields, samtools_view, Rng};
use std::process::Command;

/// Fields of a SAM line needed to assert placement of a half-mapped pair.
///
/// Includes the mate columns deliberately: `record_key_fields` excludes
/// RNEXT/PNEXT/TLEN by design, so the whole-record comparison below cannot see
/// them -- and the mate-coordinate rewrite this test exists to pin writes
/// exactly those fields.
#[derive(PartialEq, Eq, Debug)]
struct Rec {
    flag: u32,
    rname: String,
    pos: String,
    rnext: String,
    pnext: String,
    tlen: String,
}

fn recs_named(lines: &[String], qname: &str) -> Vec<Rec> {
    let mut out: Vec<Rec> = lines
        .iter()
        .filter(|l| l.split('\t').next() == Some(qname))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            Rec {
                flag: f[1].parse().unwrap(),
                rname: f[2].to_string(),
                pos: f[3].to_string(),
                rnext: f[6].to_string(),
                pnext: f[7].to_string(),
                tlen: f[8].to_string(),
            }
        })
        .collect();
    out.sort_by_key(|r| r.flag);
    out
}

#[test]
fn cli_parity_half_mapped_pair() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    if !common::require_samtools() {
        return;
    }

    let mut rng = Rng(0x5eed_9a17);
    let seq = random_dna(&mut rng, 4000);

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let ref_fa = common::setup_phix_index(dir, &bwa, &seq);

    let (read_len, insert) = (100usize, 300usize);
    let mut r1_reads = Vec::new();
    let mut r2_reads = Vec::new();

    // A few ordinary pairs so the run has a normal insert-size distribution to
    // estimate from -- pairing behaviour on the half-mapped pair depends on it.
    for i in 0..40 {
        let start = 200 + i * 40;
        let r1 = seq.as_bytes()[start..start + read_len].to_vec();
        let r2 = common::revcomp(&seq.as_bytes()[start + insert - read_len..start + insert]);
        r1_reads.push((format!("ok{i}"), r1));
        r2_reads.push((format!("ok{i}"), r2));
    }

    // The fixture under test: R1 is drawn from the reference (so it maps), R2
    // is random and is not a substring of it (so it does not). This is the
    // branch `unmapped0` cannot reach, since there both mates fail to map.
    let half_start = 1500usize;
    r1_reads.push((
        "half0".to_string(),
        seq.as_bytes()[half_start..half_start + read_len].to_vec(),
    ));
    r2_reads.push((
        "half0".to_string(),
        random_dna(&mut rng, read_len).into_bytes(),
    ));

    // Second half-mapped pair, with the MAPPED mate on the reverse strand.
    // half0 alone cannot exercise the strand half of the rewrite: its mapped
    // mate is forward, so the copied `is_rev` is 0 and the 0x10/0x20 bits
    // agree with the unfixed code by coincidence. Here the mapped mate is
    // reverse, so the unmapped record must inherit 0x10 and the mapped record
    // must gain 0x20 from the copy.
    let rev_start = 2500usize;
    r1_reads.push((
        "half1".to_string(),
        random_dna(&mut rng, read_len).into_bytes(),
    ));
    r2_reads.push((
        "half1".to_string(),
        common::revcomp(&seq.as_bytes()[rev_start..rev_start + read_len]),
    ));

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
    assert!(status.success(), "bwa-rs mem failed");

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

    // Verify the fixture actually produced a HALF-mapped pair, rather than
    // both mates mapping (which would make every assertion below vacuous) or
    // both failing (which is the case `unmapped0` already covers).
    let cli_half = recs_named(&cli_lines, "half0");
    assert_eq!(
        cli_half.len(),
        2,
        "expected exactly two \"half0\" records from the reference aligner, got {}",
        cli_half.len()
    );
    let mapped = cli_half.iter().filter(|r| r.flag & 0x4 == 0).count();
    let unmapped = cli_half.iter().filter(|r| r.flag & 0x4 != 0).count();
    assert_eq!(
        (mapped, unmapped),
        (1, 1),
        "fixture must be HALF-mapped (one mate mapped, one not); got {mapped} mapped / \
         {unmapped} unmapped — the copy-mate-to-alignment branch would go unexercised"
    );

    // The placement rule itself: upstream gives the unmapped record its mate's
    // RNAME and POS rather than `*`/0.
    let cli_un = cli_half.iter().find(|r| r.flag & 0x4 != 0).unwrap();
    let cli_map = cli_half.iter().find(|r| r.flag & 0x4 == 0).unwrap();
    assert_eq!(
        (cli_un.rname.as_str(), cli_un.pos.as_str()),
        (cli_map.rname.as_str(), cli_map.pos.as_str()),
        "reference aligner should place the unmapped mate at the mapped mate's \
         coordinates — if this ever changes upstream, the shim's mirror must follow"
    );

    // half1: the mapped mate is reverse-strand, so the rewrite must propagate
    // that strand to BOTH records -- 0x10 on the unmapped one (its own strand
    // bit, taken from the copied value) and 0x20 on the mapped one (its mate's
    // strand, taken from its own). Asserted against the reference aligner
    // rather than hardcoded, so this tracks upstream if it ever changes.
    let cli_half1 = recs_named(&cli_lines, "half1");
    let rs_half1 = recs_named(&rs_lines, "half1");
    assert_eq!(
        cli_half1.len(),
        2,
        "expected two \"half1\" records from the CLI"
    );
    let cli_h1_un = cli_half1.iter().find(|r| r.flag & 0x4 != 0).unwrap();
    let cli_h1_map = cli_half1.iter().find(|r| r.flag & 0x4 == 0).unwrap();
    assert!(
        cli_h1_map.flag & 0x10 != 0,
        "fixture must map half1's mate on the REVERSE strand to exercise the \
         strand copy; got flag {}",
        cli_h1_map.flag
    );
    assert_eq!(
        cli_h1_un.flag & 0x10,
        0x10,
        "reference aligner gives the unmapped record its mate's strand (0x10)"
    );
    let mut cli_h1_flags: Vec<u32> = cli_half1.iter().map(|r| r.flag).collect();
    let mut rs_h1_flags: Vec<u32> = rs_half1.iter().map(|r| r.flag).collect();
    cli_h1_flags.sort();
    rs_h1_flags.sort();
    assert_eq!(
        cli_h1_flags, rs_h1_flags,
        "half1 flags diverge: the 0x10/0x20 bits must be derived after the \
         mate-coordinate copy, not from the raw values"
    );

    // Mate columns, for both half-mapped pairs. These are outside what
    // `record_key_fields` compares, yet the rewrite under test sets them, so
    // without this the RNEXT/PNEXT half of the change would be unverified.
    for qname in ["half0", "half1"] {
        assert_eq!(
            recs_named(&cli_lines, qname),
            recs_named(&rs_lines, qname),
            "{qname}: mate columns (RNEXT/PNEXT/TLEN) diverge from the bwa-mem3 CLI"
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
        "bwa-rs records diverge from the bwa-mem3 CLI on a half-mapped pair"
    );
}
