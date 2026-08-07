//! Non-meth CLI-parity on TLEN for pairs that are not canonical FR.
//!
//! Upstream computes TLEN from strand-aware endpoints (`bwamem.cpp:2532-2539`):
//!
//! ```text
//! p0 = p->pos + (p->is_rev ? rlen(p) - 1 : 0)
//! p1 = m->pos + (m->is_rev ? rlen(m) - 1 : 0)
//! tlen = -(p0 - p1 + (p0 > p1 ? 1 : p0 < p1 ? -1 : 0))
//! ```
//!
//! and emits `0` when either mate has no CIGAR. The shim instead computed the
//! span `max(end) - min(begin)` and took its sign from `a_begin <= b_begin`.
//!
//! Those agree on a canonical FR pair, which is all any existing fixture
//! produces — so the divergence was invisible. They disagree on:
//!
//! - **a fully-overlapping pair**, where both mates start at the same base.
//!   `a_begin <= b_begin` is then true for *both* records, so both got a
//!   POSITIVE TLEN and the pair no longer summed to zero — a SAM spec
//!   violation, not merely a parity difference. Amplicon, cfDNA, and
//!   duplex-UMI panels hit this routinely.
//! - **a contained pair**, where one mate lies inside the other's span; the
//!   span formula reports the outer length rather than the insert.
//!
//! `record_key_fields` excludes TLEN by design, so the whole-record comparison
//! cannot see any of this; the assertions below name the column explicitly.
//!
//! Requires `bwa-mem3` + `samtools`; set `BWA_MEM3_BIN` if not on PATH. Skips
//! gracefully otherwise.

mod common;

use common::{random_dna, record_key_fields, samtools_view, Rng};
use std::process::Command;

/// Reference bases a CIGAR consumes: the M/D/N/=/X opcodes, matching
/// upstream's `get_rlen` (`bwamem.cpp:2768-2777`, which sums M and D) plus the
/// spec opcodes `samtools view` can also emit. Used only to check the
/// contained-pair fixture, not to compute anything asserted.
fn ref_span(cigar: &str) -> i64 {
    let mut span = 0;
    let mut len = 0i64;
    for c in cigar.chars() {
        if let Some(d) = c.to_digit(10) {
            len = len * 10 + i64::from(d);
        } else {
            if matches!(c, 'M' | 'D' | 'N' | '=' | 'X') {
                span += len;
            }
            len = 0;
        }
    }
    span
}

/// (flag, tlen) for one record, which is all this test compares.
fn tlens_named(lines: &[String], qname: &str) -> Vec<(u32, i64)> {
    let mut out: Vec<(u32, i64)> = lines
        .iter()
        .filter(|l| l.split('\t').next() == Some(qname))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[1].parse().unwrap(), f[8].parse().unwrap())
        })
        .collect();
    // Sort the whole tuple, not just the FLAG: two records for one QNAME
    // sharing a FLAG would otherwise leave the comparison dependent on each
    // aligner's emission order. The current fixtures cannot collide (the mates
    // differ on 0x40/0x80), so this removes a latent dependency rather than a
    // live bug.
    out.sort();
    out
}

#[test]
fn cli_parity_tlen_on_overlapping_pairs() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    if !common::require_samtools() {
        return;
    }

    let mut rng = Rng(0x71e0_5a17);
    let seq = random_dna(&mut rng, 4000);

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let ref_fa = common::setup_phix_index(dir, &bwa, &seq);

    let (read_len, insert) = (100usize, 300usize);
    let mut r1_reads = Vec::new();
    let mut r2_reads = Vec::new();

    // Ordinary FR pairs: both to give the run an insert-size distribution and
    // to keep a canonical case in the comparison, since the fix must not
    // perturb the orientation the old formula got right.
    for i in 0..40 {
        let start = 200 + i * 40;
        let r1 = seq.as_bytes()[start..start + read_len].to_vec();
        let r2 = common::revcomp(&seq.as_bytes()[start + insert - read_len..start + insert]);
        r1_reads.push((format!("ok{i}"), r1));
        r2_reads.push((format!("ok{i}"), r2));
    }

    // Fully overlapping: both mates cover exactly the same interval, so both
    // records start at the same base. This is the sign-flip case.
    let ovl_start = 1500usize;
    r1_reads.push((
        "ovl0".to_string(),
        seq.as_bytes()[ovl_start..ovl_start + read_len].to_vec(),
    ));
    r2_reads.push((
        "ovl0".to_string(),
        common::revcomp(&seq.as_bytes()[ovl_start..ovl_start + read_len]),
    ));

    // Contained: the reverse mate lies strictly inside the forward mate's
    // span, so the outer span overstates the insert.
    let cont_start = 2500usize;
    r1_reads.push((
        "cont0".to_string(),
        seq.as_bytes()[cont_start..cont_start + 150].to_vec(),
    ));
    r2_reads.push((
        "cont0".to_string(),
        common::revcomp(&seq.as_bytes()[cont_start + 20..cont_start + 80]),
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

    // Fixture validity: two mapped records, both placed at the same base.
    // Without that, the sign-flip branch is never taken and the assertion
    // below passes for the wrong reason.
    let cli_ovl = tlens_named(&cli_lines, "ovl0");
    assert_eq!(
        cli_ovl.len(),
        2,
        "expected two \"ovl0\" records from the CLI"
    );
    let ovl_pos: Vec<&str> = cli_lines
        .iter()
        .filter(|l| l.split('\t').next() == Some("ovl0"))
        .map(|l| l.split('\t').nth(3).unwrap())
        .collect();
    assert_eq!(
        ovl_pos[0], ovl_pos[1],
        "fixture must place both mates at the SAME position to exercise the \
         sign-flip; got {ovl_pos:?}"
    );

    // The spec property the old formula broke: a pair's TLENs sum to zero.
    // Asserted on the reference aligner first, so this tracks upstream.
    let cli_sum: i64 = cli_ovl.iter().map(|(_, t)| *t).sum();
    assert_eq!(
        cli_sum, 0,
        "reference aligner's TLENs must sum to zero across a pair; got {cli_ovl:?}"
    );
    assert_ne!(
        cli_ovl[0].1, 0,
        "fixture must produce a non-zero TLEN, else the sign assertion is vacuous"
    );

    // Fixture validity for the contained case, for the same reason as ovl0
    // above: if the reverse mate ever stops landing strictly inside the
    // forward mate's aligned interval, this stays green while testing nothing.
    let cont: Vec<(u32, i64, i64)> = cli_lines
        .iter()
        .filter(|l| l.split('\t').next() == Some("cont0"))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let pos: i64 = f[3].parse().unwrap();
            (f[1].parse().unwrap(), pos, pos + ref_span(f[5]) - 1)
        })
        .collect();
    assert_eq!(cont.len(), 2, "expected two \"cont0\" records from the CLI");
    let fwd = cont
        .iter()
        .find(|(flag, ..)| flag & 0x10 == 0)
        .expect("a forward-strand cont0 record");
    let rev = cont
        .iter()
        .find(|(flag, ..)| flag & 0x10 != 0)
        .expect("a reverse-strand cont0 record");
    assert!(
        rev.1 >= fwd.1 && rev.2 <= fwd.2,
        "fixture must CONTAIN the reverse mate inside the forward mate: \
         reverse [{}, {}] must lie within forward [{}, {}]",
        rev.1,
        rev.2,
        fwd.1,
        fwd.2
    );

    for qname in ["ovl0", "cont0"] {
        assert_eq!(
            tlens_named(&cli_lines, qname),
            tlens_named(&rs_lines, qname),
            "{qname}: TLEN diverges from the bwa-mem3 CLI — it must come from \
             upstream's strand-aware endpoint formula, not a span"
        );
    }

    // The ordinary FR pairs must be unaffected by the change.
    for i in 0..40 {
        let qname = format!("ok{i}");
        assert_eq!(
            tlens_named(&cli_lines, &qname),
            tlens_named(&rs_lines, &qname),
            "{qname}: TLEN on a canonical FR pair must not regress"
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
        "bwa-rs records diverge from the bwa-mem3 CLI on overlapping pairs"
    );
}
