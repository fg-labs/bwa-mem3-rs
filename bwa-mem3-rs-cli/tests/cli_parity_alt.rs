//! CLI parity over an ALT-contig fixture.
//!
//! A read that maps equally well to a primary contig and to an `is_alt` contig
//! gets an `alt_sc > 0`, which upstream turns into a `pa:f` tag and (in the
//! paired branch) an ALT supplementary record; the `is_alt` flag also changes
//! how supplementary MAPQ is lowered. None of that is reachable by the
//! single-contig, no-`.alt` fixtures the rest of the suite uses, so the shim's
//! `pa:f` emission and its `is_alt`-guarded paired/supp handling went
//! unexercised — exactly the code this diff added. This fixture builds a
//! two-contig reference whose second contig is a copy of a region of the first
//! and is marked ALT via a `<prefix>.alt` file, then asserts byte-for-byte
//! parity with `bwa-mem3 mem` over the full record (tags included).
//!
//! Requires `bwa-mem3` + `samtools`; set `BWA_MEM3_BIN` if not on PATH. Skips
//! gracefully otherwise.

mod common;

use common::{random_dna, require_bwa_mem3, require_samtools, revcomp, samtools_view, Rng};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

/// A full alignment record reduced to comparable text: the 11 mandatory
/// columns plus every optional tag, tags sorted so emission order can't matter.
fn norm_records(lines: &[String]) -> Vec<String> {
    let mut out: Vec<String> = lines
        .iter()
        .filter(|l| !l.starts_with('@'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let mut cols: Vec<String> = f[..11].iter().map(|s| (*s).to_string()).collect();
            let mut tags: Vec<&str> = f[11..].to_vec();
            tags.sort_unstable();
            cols.push(tags.join("\t"));
            cols.join("\t")
        })
        .collect();
    out.sort();
    out
}

/// Write a two-contig FASTA (`chr` + `chrALT`), index it with `bwa-mem3 index`,
/// and write a `<prefix>.alt` marking `chrALT` as an alternate locus (bwa keys
/// `is_alt` off the first field of each non-`@` line — bntseq.cpp:236-243).
fn setup_alt_index(dir: &Path, bwa: &str, chr: &[u8], chr_alt: &[u8]) -> PathBuf {
    let ref_fa = dir.join("ref.fa");
    {
        let mut f = std::fs::File::create(&ref_fa).unwrap();
        writeln!(f, ">chr").unwrap();
        for c in chr.chunks(72) {
            f.write_all(c).unwrap();
            writeln!(f).unwrap();
        }
        writeln!(f, ">chrALT").unwrap();
        for c in chr_alt.chunks(72) {
            f.write_all(c).unwrap();
            writeln!(f).unwrap();
        }
    }
    let status = Command::new(bwa)
        .arg("index")
        .arg(&ref_fa)
        .status()
        .expect("run bwa-mem3 index");
    assert!(status.success(), "bwa-mem3 index failed");
    // A single line naming the ALT contig is enough: the parser only reads the
    // first tab-delimited field and marks that contig `is_alt`.
    let alt = dir.join("ref.fa.alt");
    std::fs::write(&alt, "chrALT\t0\tchr\t1\t60\t*\t*\t0\t0\t*\t*\n").unwrap();
    ref_fa
}

#[test]
fn cli_parity_alt_contig_pa_and_supplementary() {
    let Some(bwa) = require_bwa_mem3() else {
        return;
    };
    if !require_samtools() {
        return;
    }

    let mut rng = Rng(0x0A17_0C0D_E000_0001);
    let chr = random_dna(&mut rng, 6000).into_bytes();
    // chrALT is a *mutated* copy of a 600 bp window of chr. Reads are drawn from
    // chrALT, so each matches chrALT exactly but chr with a few mismatches ->
    // the chrALT hit scores higher. Because chrALT is marked ALT, upstream still
    // emits the (lower-scoring) chr hit as the non-ALT primary but records the
    // ALT hit's score as `alt_sc` (bwamem.cpp:4337-4339), which is what drives
    // `pa:f`. An exact copy would tie, and the non-ALT hit would win outright
    // with no alt_sc -- so the divergence has to come from a mutated copy.
    let alt_lo = 2000usize;
    let alt_len = 600usize;
    let mut chr_alt = chr[alt_lo..alt_lo + alt_len].to_vec();
    // ~1 substitution per 120 bp: enough that a 150 bp read spans one or two
    // (so chr scores strictly below chrALT), few enough that the chr hit stays
    // well above -T 30.
    for k in (60..alt_len).step_by(120) {
        chr_alt[k] = match chr_alt[k] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let ref_fa = setup_alt_index(dir, &bwa, &chr, &chr_alt);

    let (read_len, insert) = (150usize, 400usize);
    let mut r1: Vec<(String, Vec<u8>)> = Vec::new();
    let mut r2: Vec<(String, Vec<u8>)> = Vec::new();

    // Concordant pairs across the unique part of chr, to seed a normal
    // insert-size model.
    for i in 0..200 {
        let start = 100 + i * 25;
        if start + insert > alt_lo {
            break;
        }
        r1.push((format!("ok{i}"), chr[start..start + read_len].to_vec()));
        r2.push((
            format!("ok{i}"),
            revcomp(&chr[start + insert - read_len..start + insert]),
        ));
    }

    // Pairs drawn from chrALT (insert 400 fits the 600 bp copy): each read
    // matches chrALT exactly and chr with the seeded mismatches, so chrALT
    // scores higher and the emitted chr primary carries alt_sc > 0 -> pa:f, with
    // the ALT hit surfacing as an is_alt supplementary in the paired branch.
    for (i, j) in [20usize, 60, 100, 140, 180].into_iter().enumerate() {
        r1.push((format!("alt{i}"), chr_alt[j..j + read_len].to_vec()));
        r2.push((
            format!("alt{i}"),
            revcomp(&chr_alt[j + insert - read_len..j + insert]),
        ));
    }

    let r1_fq = dir.join("r1.fq");
    let r2_fq = dir.join("r2.fq");
    common::write_fastq(&r1_fq, &r1);
    common::write_fastq(&r2_fq, &r2);

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

    let cli_lines = samtools_view(&cli_bam);
    let rs_lines = samtools_view(&rs_bam);

    // Vacuity guard: the fixture must actually exercise the ALT path, or the
    // comparison says nothing about the code this test exists to pin.
    let cli_pa = cli_lines
        .iter()
        .filter(|l| !l.starts_with('@'))
        .filter(|l| l.split('\t').skip(11).any(|t| t.starts_with("pa:f:")))
        .count();
    assert!(
        cli_pa > 0,
        "fixture produced no pa:f records, so ALT/alt_sc parity is untested \
         (the ALT contig may not be competing as intended)"
    );
    // `pa:f` only proves alt_sc > 0; it does not prove the ALT hit surfaced as a
    // supplementary record. The paired branch emits that hit as an is_alt
    // supplementary (FLAG 0x800) on chrALT, which is the emission this test pins,
    // so require at least one such record or the is_alt-supp path is untested.
    let cli_alt_supp = cli_lines
        .iter()
        .filter(|l| !l.starts_with('@'))
        .filter(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let flag: u32 = f[1].parse().unwrap_or(0);
            f[0].starts_with("alt") && flag & 0x800 != 0 && f[2] == "chrALT"
        })
        .count();
    assert!(
        cli_alt_supp > 0,
        "fixture produced no is_alt supplementary (a 0x800 record on chrALT), so \
         the ALT-supplementary emission path is untested"
    );

    let cli = norm_records(&cli_lines);
    let rs = norm_records(&rs_lines);
    assert_eq!(
        rs.len(),
        cli.len(),
        "record count mismatch: bwa-rs={} CLI={}",
        rs.len(),
        cli.len()
    );
    let diffs: Vec<String> = rs
        .iter()
        .zip(&cli)
        .filter(|(a, b)| a != b)
        .map(|(a, b)| format!("  bwa-rs: {a}\n  CLI:    {b}"))
        .collect();
    assert!(
        diffs.is_empty(),
        "{} of {} records diverge from bwa-mem3 mem on the ALT fixture:\n{}",
        diffs.len(),
        cli.len(),
        diffs.iter().take(8).cloned().collect::<Vec<_>>().join("\n")
    );
}
