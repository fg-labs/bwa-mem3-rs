//! Non-meth CLI-parity over a split (chimeric) read, which yields a
//! supplementary record.
//!
//! Upstream hard-clips supplementary alignments by default. `mem_reg2aln`
//! only ever writes SOFT clips (`bwamem.cpp:2735,2743` build `clip<<4 | 3`,
//! and `3` is `S` in bwa's `MIDSH` opcode table — it never emits `H`), so the
//! rewrite happens at write time in `add_cigar` (`bwamem.cpp:2410-2421`):
//! when `which != 0` (a 2nd+ emitted record), `!(opt->flag&MEM_F_SOFTCLIP)`
//! and `!p->is_alt`, each `S` becomes `H`. `mem_aln2sam` then trims SEQ/QUAL
//! to match (`:2547-2565`), and rewrites the mate CIGAR in `MC:Z` through the
//! same helper (`:2583`).
//!
//! The shim did none of that: it wrote `p->cigar` verbatim and trimmed SEQ
//! only when the opcode was already `H` — a condition `mem_reg2aln` never
//! produces, so that branch was dead and every supplementary came out
//! soft-clipped with a full-length SEQ.
//!
//! No existing fixture reaches this: every simulated read is a contiguous
//! substring of the reference, so nothing splits. `flag_parity` only asserts
//! records are *primary*, and `align_smoke` explicitly tolerates either.
//!
//! Requires `bwa-mem3` + `samtools`; set `BWA_MEM3_BIN` if not on PATH. Skips
//! gracefully otherwise.

mod common;

use common::{random_dna, record_key_fields, samtools_view, Rng};
use std::process::Command;

/// The SAM fields needed to pin a supplementary record's clipping.
///
/// SEQ is included deliberately: `record_key_fields` excludes it by design,
/// and the trim this test exists to pin writes exactly that column.
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug)]
struct Rec {
    flag: u32,
    cigar: String,
    seq_len: usize,
    qual_len: usize,
    mc: Option<String>,
}

fn recs_named(lines: &[String], qname: &str) -> Vec<Rec> {
    let mut out: Vec<Rec> = lines
        .iter()
        .filter(|l| l.split('\t').next() == Some(qname))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            Rec {
                flag: f[1].parse().unwrap(),
                cigar: f[5].to_string(),
                seq_len: f[9].len(),
                qual_len: f[10].len(),
                mc: f
                    .get(11..)
                    .unwrap_or_default()
                    .iter()
                    .find(|t| t.starts_with("MC:Z:"))
                    .map(|t| (*t).to_string()),
            }
        })
        .collect();
    // Total order (derived, so FLAG then CIGAR then the rest): sorting on FLAG
    // alone would leave records that share one comparing in each aligner's
    // emission order.
    out.sort();
    out
}

#[test]
fn cli_parity_supplementary_is_hard_clipped() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    if !common::require_samtools() {
        return;
    }

    let mut rng = Rng(0x5b1d_c0de);
    let seq = random_dna(&mut rng, 4000);

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let ref_fa = common::setup_phix_index(dir, &bwa, &seq);

    let (read_len, insert) = (100usize, 300usize);
    let mut r1_reads = Vec::new();
    let mut r2_reads = Vec::new();

    // Ordinary pairs so the run has a normal insert-size distribution.
    for i in 0..40 {
        let start = 200 + i * 40;
        let r1 = seq.as_bytes()[start..start + read_len].to_vec();
        let r2 = common::revcomp(&seq.as_bytes()[start + insert - read_len..start + insert]);
        r1_reads.push((format!("ok{i}"), r1));
        r2_reads.push((format!("ok{i}"), r2));
    }

    // The fixture under test: R1 is chimeric — two 60bp arms drawn from loci
    // ~2.5kb apart, so neither arm can be explained as an indel and bwa emits
    // a primary plus a supplementary. Each arm scores 60, comfortably above
    // the default `-T 30`, so both survive as non-secondary alnregs.
    let (arm_a, arm_b) = (500usize, 3000usize);
    let arm = 60usize;
    let mut chimera = seq.as_bytes()[arm_a..arm_a + arm].to_vec();
    chimera.extend_from_slice(&seq.as_bytes()[arm_b..arm_b + arm]);
    r1_reads.push(("split0".to_string(), chimera));
    // The mate is anchored near the first arm, but carries a 15bp random
    // prefix so it aligns with a SOFT clip of its own. That is what exercises
    // the `MC:Z` half of the rule: upstream builds the mate CIGAR through the
    // same `add_cigar` helper and passes it the *record's* `which`
    // (`bwamem.cpp:2583`), so on the supplementary record the mate's clip is
    // rewritten to `H` too — even though the mate itself is a primary. Without
    // a clipped mate, `MC:Z` is `100M` on every record and the rewrite goes
    // unverified.
    let mut clipped_mate = random_dna(&mut rng, 15).into_bytes();
    clipped_mate.extend_from_slice(&common::revcomp(
        &seq.as_bytes()[arm_a + insert - 85..arm_a + insert],
    ));
    r2_reads.push(("split0".to_string(), clipped_mate));

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

    let cli_split = recs_named(&cli_lines, "split0");
    let rs_split = recs_named(&rs_lines, "split0");

    // Fixture validity: without an actual supplementary every assertion below
    // is vacuous. This is the guard the rest of the suite lacks.
    let cli_supp: Vec<&Rec> = cli_split.iter().filter(|r| r.flag & 0x800 != 0).collect();
    assert_eq!(
        cli_supp.len(),
        1,
        "fixture must produce exactly one supplementary record from the \
         reference aligner; got {} (records: {cli_split:?})",
        cli_supp.len()
    );

    // The rule itself, asserted against the reference aligner rather than
    // hardcoded, so it tracks upstream if the default ever changes.
    let supp = cli_supp[0];
    assert!(
        supp.cigar.contains('H'),
        "reference aligner should HARD-clip a supplementary by default; got \
         CIGAR {} — if upstream changed this, the shim's mirror must follow",
        supp.cigar
    );
    assert!(
        supp.seq_len < 120,
        "reference aligner should trim SEQ to the hard-clipped span; got a \
         full-length SEQ of {} bases",
        supp.seq_len
    );
    assert_eq!(
        supp.seq_len, supp.qual_len,
        "SEQ and QUAL must be trimmed together"
    );

    // Fixture validity for the MC:Z half: the mate must itself be soft-clipped,
    // and upstream must then hard-clip it in the supplementary's MC:Z. Without
    // both, the MC rewrite is untested.
    let cli_primary = cli_split
        .iter()
        .find(|r| r.flag & 0x900 == 0)
        .expect("a primary record for split0");
    assert!(
        cli_primary.mc.as_deref().is_some_and(|m| m.contains('S')),
        "fixture must give the mate a SOFT clip so the MC:Z rewrite is \
         exercised; primary carries {:?}",
        cli_primary.mc
    );
    assert!(
        supp.mc.as_deref().is_some_and(|m| m.contains('H')),
        "reference aligner rewrites the mate CIGAR through the same rule on a \
         supplementary record, so MC:Z should be hard-clipped there; got {:?}",
        supp.mc
    );

    // The parity assertion: CIGAR, SEQ/QUAL length, and MC:Z on every record
    // of the split read.
    assert_eq!(
        cli_split, rs_split,
        "split0: supplementary clipping diverges from the bwa-mem3 CLI — \
         supplementaries must be hard-clipped, SEQ/QUAL trimmed to match, and \
         MC:Z rewritten through the same rule"
    );

    // Sorted vectors (not sets): preserve record multiplicity so a duplicated
    // or missing record diverges even when the distinct key set matches.
    let mut cli_recs: Vec<String> = cli_lines.iter().map(|l| record_key_fields(l)).collect();
    let mut rs_recs: Vec<String> = rs_lines.iter().map(|l| record_key_fields(l)).collect();
    cli_recs.sort();
    rs_recs.sort();
    assert_eq!(
        cli_recs, rs_recs,
        "bwa-rs records diverge from the bwa-mem3 CLI on a split read"
    );
}
