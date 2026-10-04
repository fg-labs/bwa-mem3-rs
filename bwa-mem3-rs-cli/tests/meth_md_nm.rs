//! `--meth` NM/MD contract on an explicit bisulfite-converted pair.
//!
//! Since bwa-mem3 v0.14.0 (upstream `f0cb34a2`), `--meth` follows the BISCUIT
//! convention: `MD:Z` is literal — it lists every reference base the read
//! differs from, bisulfite conversions included, so CIGAR + SEQ + MD still
//! rebuild the reference — while `NM:i` excludes conversions. The CLI-parity
//! suites (`meth_e2e`, `meth_cli_parity_*`) only prove bwa-rs agrees with the
//! CLI; this test pins the contract itself, against values computed here.
//!
//! The fixture is one directional pair from a PhiX window: the fragment's top
//! strand has every non-CpG C converted to T (OT), R1 reads its start and R2
//! the reverse complement of its end (so R2 carries the OB-projected G->A).
//! R1 also carries one real A->G substitution, which NM must count.
//!
//! Requires `bwa-mem3` (for `index --meth`) and `samtools`; skips otherwise.

mod common;
mod phix_seq;

use common::samtools_view;
use std::process::Command;

const FRAG_START: usize = 1000;
const FRAG_LEN: usize = 300;
const READ_LEN: usize = 100;
/// The real (non-bisulfite) substitution goes on the first reference A at or
/// after this R1 offset, mid-read so it cannot be soft-clipped away.
const SNV_FROM: usize = 50;

/// Bisulfite-convert a top-strand sequence: every C not followed by G becomes T.
fn convert_non_cpg(seq: &[u8]) -> Vec<u8> {
    (0..seq.len())
        .map(|i| {
            if seq[i] == b'C' && seq.get(i + 1) != Some(&b'G') {
                b'T'
            } else {
                seq[i]
            }
        })
        .collect()
}

/// The MD string of a gapless, unclipped alignment of `read` to `reference`,
/// with every literal base difference listed.
fn literal_md(reference: &[u8], read: &[u8]) -> String {
    let mut md = String::new();
    let mut run = 0;
    for (r, q) in reference.iter().zip(read) {
        if r == q {
            run += 1;
        } else {
            md.push_str(&run.to_string());
            md.push(*r as char);
            run = 0;
        }
    }
    md.push_str(&run.to_string());
    md
}

fn tag<'a>(fields: &'a [&'a str], prefix: &str) -> &'a str {
    fields[11..]
        .iter()
        .find_map(|f| f.strip_prefix(prefix))
        .unwrap_or_else(|| panic!("record {} missing {prefix}", fields[0]))
}

#[test]
fn meth_md_is_literal_and_nm_excludes_conversions() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    if !common::require_samtools() {
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let ref_fa = common::setup_phix_meth_index(dir, &bwa, phix_seq::PHIX_SEQ);

    let reference = phix_seq::PHIX_SEQ.as_bytes();
    let frag_ref = &reference[FRAG_START..FRAG_START + FRAG_LEN];
    let frag = convert_non_cpg(frag_ref);

    let r1_ref = &frag_ref[..READ_LEN];
    let mut r1 = frag[..READ_LEN].to_vec();
    let snv = SNV_FROM + r1_ref[SNV_FROM..].iter().position(|&b| b == b'A').unwrap();
    r1[snv] = b'G';

    let r2_ref = &frag_ref[FRAG_LEN - READ_LEN..];
    let r2_fwd = &frag[FRAG_LEN - READ_LEN..];
    let r2 = common::revcomp(r2_fwd);

    let r1_fq = dir.join("r1.fq");
    let r2_fq = dir.join("r2.fq");
    common::write_fastq(&r1_fq, &[("frag".to_string(), r1.clone())]);
    common::write_fastq(&r2_fq, &[("frag".to_string(), r2)]);

    let bam = dir.join("out.bam");
    let status = Command::new(common::cli_bin())
        .args(["mem", "--meth"])
        .arg(&ref_fa)
        .arg(&r1_fq)
        .arg(&r2_fq)
        .arg("-o")
        .arg(&bam)
        .status()
        .expect("run bwa-rs mem --meth");
    assert!(status.success(), "bwa-rs mem --meth failed");

    let lines = samtools_view(&bam);
    let primaries: Vec<Vec<&str>> = lines
        .iter()
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|f| f[1].parse::<u32>().unwrap() & 0x900 == 0)
        .collect();
    assert_eq!(
        primaries.len(),
        2,
        "expected one primary per mate: {lines:?}"
    );

    // (first-in-pair, expected POS, expected reverse, reference window, forward read, expected NM)
    let expected = [
        (true, FRAG_START + 1, false, r1_ref, r1.as_slice(), 1),
        (
            false,
            FRAG_START + FRAG_LEN - READ_LEN + 1,
            true,
            r2_ref,
            r2_fwd,
            0,
        ),
    ];
    for (first, pos, reverse, window, read_fwd, nm) in expected {
        let rec = primaries
            .iter()
            .find(|f| (f[1].parse::<u32>().unwrap() & 0x40 != 0) == first)
            .expect("mate record");
        let flag: u32 = rec[1].parse().unwrap();
        assert_eq!(rec[3].parse::<usize>().unwrap(), pos, "{} POS", rec[0]);
        assert_eq!(flag & 0x10 != 0, reverse, "{} strand", rec[0]);
        assert_eq!(rec[5], format!("{READ_LEN}M"), "{} CIGAR", rec[0]);

        let md = literal_md(window, read_fwd);
        let conversions = window
            .iter()
            .zip(read_fwd)
            .filter(|(r, q)| **r == b'C' && **q == b'T')
            .count();
        assert!(conversions > 0, "fixture: {} has no conversions", rec[0]);
        assert_eq!(
            tag(rec, "MD:Z:"),
            md,
            "{} MD must list every conversion",
            rec[0]
        );
        assert_eq!(
            tag(rec, "NM:i:").parse::<usize>().unwrap(),
            nm,
            "{} NM must exclude the {conversions} conversions",
            rec[0]
        );
    }
}
