//! Bisulfite (`--meth`) CLI-parity over a split (chimeric) read.
//!
//! The plain [`cli_parity_supplementary`] test pins the hard-clip rewrite on
//! the non-meth path. It cannot reach the `--meth` side of that change:
//! `append_bam_record` hands `meth_build_xm` its own copy of the CIGAR, and
//! those opcodes have to agree with the emitted (hard-clipped, therefore
//! trimmed) `seq_text` the tag is computed over. Only a 2nd+ emitted record
//! is clipped at all, so only there can the two disagree.
//!
//! The clip must be at the read's 5' end for this to be observable, which is
//! why the fixture uses asymmetric arms. `meth_build_xm` walks the CIGAR under
//! `r < l_emit` (`meth_xm.cpp:152`), consuming read bases for `S` (op 4) and
//! nothing for `H` (op 5). With a TRAILING clip the `M` has already advanced
//! the cursor to `l_emit`, the loop exits, and the clip opcode is never
//! examined -- S and H are indistinguishable. With a LEADING clip the walk
//! spends every emitted slot writing `.` for an unrewritten `S` and never
//! reaches the `M`, so the tag degenerates to all dots.
//!
//! [`meth_e2e`] does not cover any of this: its reads are contiguous
//! substrings of PhiX, so nothing splits and no clipped record is produced.
//!
//! Requires `bwa-mem3` (for `index --meth` + the reference aligner) and
//! `samtools`; set `BWA_MEM3_BIN` if not on PATH. Skips gracefully otherwise.

mod common;

use common::{random_dna, record_key_fields, samtools_view, Rng};
use std::process::Command;

/// A record's FLAG, clipping, and Bismark tags -- all compared for parity.
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug)]
struct MethRec {
    flag: u32,
    cigar: String,
    seq_len: usize,
    xm: Option<String>,
    xg: Option<String>,
    xr: Option<String>,
}

fn tag(fields: &[&str], key: &str) -> Option<String> {
    fields
        .get(11..)
        .unwrap_or_default()
        .iter()
        .find(|t| t.starts_with(key))
        .map(|t| (*t).to_string())
}

fn recs_named(lines: &[String], qname: &str) -> Vec<MethRec> {
    let mut out: Vec<MethRec> = lines
        .iter()
        .filter(|l| l.split('\t').next() == Some(qname))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            MethRec {
                flag: f[1].parse().unwrap(),
                cigar: f[5].to_string(),
                seq_len: f[9].len(),
                xm: tag(&f, "XM:Z:"),
                xg: tag(&f, "XG:Z:"),
                xr: tag(&f, "XR:Z:"),
            }
        })
        .collect();
    // Total order (derived, so FLAG then CIGAR then the rest), not FLAG alone:
    // records sharing a FLAG would otherwise compare in each aligner's
    // emission order.
    out.sort();
    out
}

#[test]
fn meth_cli_parity_supplementary_is_hard_clipped() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    if !common::require_samtools() {
        return;
    }

    let mut rng = Rng(0x5b1d_dead);
    let seq = random_dna(&mut rng, 4000);

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    // Build the meth DUAL index (<ref>.* + <ref>.meth.*).
    let ref_fa = common::setup_phix_meth_index(dir, &bwa, &seq);

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

    // The chimeric read: two arms from loci ~2.5kb apart, so bwa emits a
    // primary plus a non-primary second record.
    //
    // The arms are deliberately ASYMMETRIC (75 then 55). Equal arms give a
    // second record clipped at its 3' end (`60M60H`), and that cannot
    // distinguish the clip rewrite at all: `meth_build_xm` walks the CIGAR
    // under `r < l_emit` (meth_xm.cpp:152), so a trailing clip is never
    // reached -- `M` has already advanced the read cursor to the end -- and S
    // (op 4, consumes read) behaves identically to H (op 5, consumes nothing).
    //
    // A longer first arm makes it the primary and leaves the second record
    // clipped at its 5' end (`75H55M`). There the distinction is total: with
    // the clip left as S, the walk spends all 55 emitted slots writing '.' for
    // the clip and never reaches the M, yielding an all-dots XM:Z.
    let (arm_a, arm_b) = (500usize, 3000usize);
    let (len_a, len_b) = (75usize, 55usize);
    let mut chimera = seq.as_bytes()[arm_a..arm_a + len_a].to_vec();
    chimera.extend_from_slice(&seq.as_bytes()[arm_b..arm_b + len_b]);
    r1_reads.push(("split0".to_string(), chimera));
    r2_reads.push((
        "split0".to_string(),
        common::revcomp(&seq.as_bytes()[arm_a + insert - read_len..arm_a + insert]),
    ));

    let r1_fq = dir.join("r1.fq");
    let r2_fq = dir.join("r2.fq");
    common::write_fastq(&r1_fq, &r1_reads);
    common::write_fastq(&r2_fq, &r2_reads);

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

    let cli_out = Command::new(&bwa)
        .args(["mem", "--meth", "-t", "1"])
        .arg(&ref_fa)
        .arg(&r1_fq)
        .arg(&r2_fq)
        .output()
        .expect("run bwa-mem3 mem --meth");
    assert!(cli_out.status.success(), "bwa-mem3 mem --meth failed");
    let cli_bam = dir.join("cli.bam");
    std::fs::write(&cli_bam, &cli_out.stdout).unwrap();

    let rs_lines = samtools_view(&rs_bam);
    let cli_lines = samtools_view(&cli_bam);

    let cli_split = recs_named(&cli_lines, "split0");
    let rs_split = recs_named(&rs_lines, "split0");

    // Fixture validity: a 2nd+ emitted record for the split read must exist
    // and must be hard-clipped, else the XM:Z-vs-CIGAR agreement this test
    // exists to pin is never exercised.
    //
    // Selected on 0x900, not 0x800: under --meth the CLI marks the split arm
    // SECONDARY (0x100) rather than supplementary. That is the MEM_F_NO_MULTI
    // path -- `pair_and_emit` tags the region with upstream's internal 0x10000
    // and `append_bam_record` remaps it to 0x100 on write (see the comment at
    // bwa_shim_align.cpp's FLAG packing, and CLAUDE.md gotcha #12). The record
    // still goes through `which != 0`, so it is hard-clipped and trimmed
    // exactly like a 0x800 supplementary; which bit upstream chooses is beside
    // the point for the rewrite under test, so match either.
    let cli_supp: Vec<&MethRec> = cli_split.iter().filter(|r| r.flag & 0x900 != 0).collect();
    assert_eq!(
        cli_supp.len(),
        1,
        "fixture must produce exactly one non-primary record for the split \
         read under --meth; got {} (records: {cli_split:?})",
        cli_supp.len()
    );
    let supp = cli_supp[0];
    assert!(
        supp.cigar.contains('H'),
        "reference aligner should hard-clip a non-primary record under --meth too; \
         got CIGAR {}",
        supp.cigar
    );
    // The clip must be LEADING, else the XM:Z assertion below cannot tell a
    // rewritten CIGAR from an unrewritten one -- see the fixture comment on
    // arm asymmetry. Mutation-verified: reverting only the XM CIGAR rewrite
    // fails this test with a leading clip and passes with a trailing one.
    let first_op = supp
        .cigar
        .trim_start_matches(|c: char| c.is_ascii_digit())
        .chars()
        .next();
    assert_eq!(
        first_op,
        Some('H'),
        "fixture must produce a 5'-clipped non-primary record (e.g. 75H55M) so \
         the XM:Z walk actually depends on the clip opcode; got {}",
        supp.cigar
    );

    // The point of this test: XM:Z must be present on the supplementary and
    // must describe exactly the emitted (trimmed) SEQ, not the full read.
    let xm = supp
        .xm
        .as_deref()
        .expect("supplementary must carry XM:Z under --meth");
    let xm_value = xm.strip_prefix("XM:Z:").unwrap();
    assert_eq!(
        xm_value.len(),
        supp.seq_len,
        "XM:Z must have one call per emitted base — {} calls for a {}-base \
         SEQ means the tag was built over the untrimmed read",
        xm_value.len(),
        supp.seq_len
    );

    // FLAG is compared too: `bwa-rs mem --meth` applies upstream's bwameth
    // default bundle (`fastmap.cpp:1513-1527`), whose `MEM_F_NO_MULTI` decides
    // whether the split arm is marked secondary (0x100) or supplementary
    // (0x800). Records are already sorted by flag in `recs_named`.
    assert_eq!(
        cli_split, rs_split,
        "split0 under --meth: FLAG, clipping, or Bismark tags diverge from the \
         bwa-mem3 CLI"
    );

    // Whole-record parity over every read, `split0` included: the FLAG
    // divergence that once forced an exclusion here is what this PR fixes.
    //
    // Sorted vectors (not sets): preserve record multiplicity so a duplicated
    // or missing record diverges even when the distinct key set matches. The
    // non-empty guard mirrors `meth_e2e` -- without it a run where nothing
    // aligned would compare two empty vectors and pass.
    let mut cli_recs: Vec<String> = cli_lines.iter().map(|l| record_key_fields(l)).collect();
    let mut rs_recs: Vec<String> = rs_lines.iter().map(|l| record_key_fields(l)).collect();
    cli_recs.sort();
    rs_recs.sort();
    assert!(
        !cli_recs.is_empty(),
        "the reference aligner produced no records at all"
    );
    assert_eq!(
        cli_recs, rs_recs,
        "bwa-rs meth records diverge from the bwa-mem3 CLI on a split read"
    );
}
