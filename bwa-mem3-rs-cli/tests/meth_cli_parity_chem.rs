//! Bisulfite (`--meth`) CLI parity where the CLI's `--meth` defaults change
//! output: SPEC30 seed pruning and the TAPS chemistry.
//!
//! `bwa-mem3 mem --meth` turns on `--meth-seed-prune=spec30`, which drops
//! short multi-hit SMEMs before SA resolution and so changes XS, XA, MAPQ and
//! sometimes placement. On PhiX every seed is unique, nothing is pruned, and the
//! other `--meth` parity suites pass whether or not the crate applies the
//! default; this crate in fact did not, and differed from the CLI on about half
//! the records of an hg38 chr22 EM-seq simulation. This fixture is a reference
//! of diverged repeat families (like Alu copies), where pruning does change the
//! CLI's output, and each test first proves that before comparing.
//!
//! `--meth=taps` flips which base `XM:Z` calls methylated and defaults the
//! scoring to NEUTRAL; it is checked the same way.
//!
//! The comparison adds `XS:i` and `HN:i` to [`record_key_fields`], since most of
//! the records pruning changes differ only in `XS`.
//!
//! Requires `bwa-mem3` (for `index --meth` + the reference aligner) and
//! `samtools`; set `BWA_MEM3_BIN` if not on PATH. Skips gracefully otherwise.

mod common;

use common::{random_dna, record_key_fields, revcomp, samtools_view, Rng};
use std::path::PathBuf;
use std::process::Command;

const N_PAIRS: usize = 2000;
const READ_LEN: usize = 100;
const FRAG_LEN: usize = 300;

/// A ~60 kb reference: 60 copies of three 300 bp repeat families, each copy
/// 8% diverged, separated by 700 bp of unique sequence.
fn repeat_family_reference(rng: &mut Rng) -> String {
    let families: Vec<String> = (0..3).map(|_| random_dna(rng, 300)).collect();
    let mut seq = String::new();
    for i in 0..60 {
        seq.push_str(&random_dna(rng, 700));
        for base in families[i % 3].chars() {
            if rng.next() % 100 < 8 {
                let others: Vec<char> = "ACGT".chars().filter(|&b| b != base).collect();
                seq.push(others[(rng.next() % 3) as usize]);
            } else {
                seq.push(base);
            }
        }
    }
    seq.push_str(&random_dna(rng, 700));
    seq
}

#[derive(Clone, Copy)]
enum Chem {
    EmSeq,
    Taps,
}

impl Chem {
    fn meth_arg(self) -> &'static str {
        match self {
            Chem::EmSeq => "--meth",
            Chem::Taps => "--meth=taps",
        }
    }
}

/// Convert one strand as the chemistry would: EM-seq converts every non-CpG C
/// (CpGs read as methylated); TAPS converts 80% of CpG Cs (the methylated ones).
fn convert(strand: &[u8], chem: Chem, rng: &mut Rng) -> Vec<u8> {
    (0..strand.len())
        .map(|i| {
            let cpg = strand.get(i + 1) == Some(&b'G');
            let converts = match chem {
                Chem::EmSeq => !cpg,
                Chem::Taps => cpg && rng.next() % 100 < 80,
            };
            if strand[i] == b'C' && converts {
                b'T'
            } else {
                strand[i]
            }
        })
        .collect()
}

/// Directional pairs: each fragment comes from the top or bottom strand, is
/// converted, and R1 reads its 5' end while R2 reads the reverse complement of
/// its 3' end.
fn simulate(reference: &[u8], chem: Chem, rng: &mut Rng) -> Vec<(String, Vec<u8>, Vec<u8>)> {
    (0..N_PAIRS)
        .map(|i| {
            let start = (rng.next() as usize) % (reference.len() - FRAG_LEN);
            let top = &reference[start..start + FRAG_LEN];
            let strand = if rng.next() % 2 == 0 {
                top.to_vec()
            } else {
                revcomp(top)
            };
            let frag = convert(&strand, chem, rng);
            let r1 = frag[..READ_LEN].to_vec();
            let r2 = revcomp(&frag[FRAG_LEN - READ_LEN..]);
            (format!("p{i}"), r1, r2)
        })
        .collect()
}

struct Fixture {
    _tmp: tempfile::TempDir,
    dir: PathBuf,
    ref_fa: PathBuf,
    r1: PathBuf,
    r2: PathBuf,
}

fn fixture(bwa: &str, chem: Chem, seed: u64) -> Fixture {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().to_path_buf();
    let mut rng = Rng(seed);
    let reference = repeat_family_reference(&mut rng);
    let ref_fa = common::setup_phix_meth_index(&dir, bwa, &reference);
    let pairs = simulate(reference.as_bytes(), chem, &mut rng);
    let (r1, r2) = (dir.join("r1.fq"), dir.join("r2.fq"));
    let (r1_reads, r2_reads): (Vec<_>, Vec<_>) = pairs
        .into_iter()
        .map(|(name, a, b)| ((name.clone(), a), (name, b)))
        .unzip();
    common::write_fastq(&r1, &r1_reads);
    common::write_fastq(&r2, &r2_reads);
    Fixture {
        _tmp: tmp,
        dir,
        ref_fa,
        r1,
        r2,
    }
}

/// The seed-pruning override both `bwa-mem3 mem` and the vendored library
/// read; cleared so an inherited value cannot mask the SPEC30 default.
const SEED_PRUNE_ENV: &str = "BWAMEM3_METH_SEED_PRUNE";

/// `bwa-mem3 mem` with `args`, plus extra environment, as SAM lines.
fn run_cli(bwa: &str, fx: &Fixture, args: &[&str], env: &[(&str, &str)], tag: &str) -> Vec<String> {
    let out = Command::new(bwa)
        .arg("mem")
        .args(args)
        .args(["-t", "1"])
        .env_remove(SEED_PRUNE_ENV)
        .envs(env.iter().copied())
        .arg(&fx.ref_fa)
        .arg(&fx.r1)
        .arg(&fx.r2)
        .output()
        .expect("run bwa-mem3 mem");
    assert!(out.status.success(), "bwa-mem3 mem {args:?} failed");
    let path = fx.dir.join(format!("cli_{tag}.bam"));
    std::fs::write(&path, &out.stdout).unwrap();
    samtools_view(&path)
}

/// `bwa-rs mem` with `args`, as SAM lines. One batch covers every pair, so
/// both tools fit one insert-size model over the same reads (gotcha #14).
fn run_rs(fx: &Fixture, args: &[&str]) -> Vec<String> {
    let bam = fx.dir.join("rs.bam");
    let status = Command::new(common::cli_bin())
        .arg("mem")
        .args(args)
        .args(["--batch-size", &N_PAIRS.to_string()])
        .env_remove(SEED_PRUNE_ENV)
        .arg(&fx.ref_fa)
        .arg(&fx.r1)
        .arg(&fx.r2)
        .arg("-o")
        .arg(&bam)
        .status()
        .expect("run bwa-rs mem");
    assert!(status.success(), "bwa-rs mem {args:?} failed");
    samtools_view(&bam)
}

/// [`record_key_fields`] plus the `XS:i` and `HN:i` values.
fn key(sam_line: &str) -> String {
    let tag = |name: &str| {
        sam_line
            .split('\t')
            .skip(11)
            .find_map(|f| f.strip_prefix(name))
            .unwrap_or("")
            .to_string()
    };
    format!(
        "{}\tXS:{}\tHN:{}",
        record_key_fields(sam_line),
        tag("XS:i:"),
        tag("HN:i:")
    )
}

/// Sorted keys, so record multiplicity counts and order does not.
fn keys(lines: &[String]) -> Vec<String> {
    let mut k: Vec<String> = lines.iter().map(|l| key(l)).collect();
    k.sort();
    k
}

/// Positions whose sorted keys differ, plus any difference in record count.
fn differing(a: &[String], b: &[String]) -> usize {
    let (ka, kb) = (keys(a), keys(b));
    ka.iter().zip(kb.iter()).filter(|(x, y)| x != y).count() + ka.len().abs_diff(kb.len())
}

fn skip_without_tools() -> Option<String> {
    let bwa = common::require_bwa_mem3()?;
    common::require_samtools().then_some(bwa)
}

fn assert_parity(rs: &[String], cli: &[String], what: &str) {
    assert!(!cli.is_empty(), "CLI produced no records");
    assert_eq!(
        keys(rs),
        keys(cli),
        "bwa-rs diverges from the bwa-mem3 CLI: {what}"
    );
}

#[test]
fn meth_seed_prune_default_matches_cli() {
    let Some(bwa) = skip_without_tools() else {
        return;
    };
    let fx = fixture(&bwa, Chem::EmSeq, 0x5eed_0001);
    let cli = run_cli(&bwa, &fx, &["--meth"], &[], "default");

    // Non-vacuous: on this reference the CLI's SPEC30 default changes output.
    let unpruned = run_cli(&bwa, &fx, &["--meth"], &[(SEED_PRUNE_ENV, "off")], "off");
    let changed = differing(&cli, &unpruned);
    assert!(
        changed > 0,
        "SPEC30 seed pruning changed no CLI record; the fixture no longer tests it"
    );

    assert_parity(
        &run_rs(&fx, &["--meth"]),
        &cli,
        "--meth (SPEC30 seed pruning)",
    );
}

#[test]
fn taps_chemistry_matches_cli() {
    let Some(bwa) = skip_without_tools() else {
        return;
    };
    let fx = fixture(&bwa, Chem::Taps, 0x5eed_0002);
    let cli = run_cli(&bwa, &fx, &[Chem::Taps.meth_arg()], &[], "taps");

    // Non-vacuous: the chemistry changes the CLI's output (at least XM:Z).
    let as_emseq = run_cli(&bwa, &fx, &[Chem::EmSeq.meth_arg()], &[], "emseq");
    assert!(
        differing(&cli, &as_emseq) > 0,
        "--meth=taps produced the same records as --meth; the chemistry is not exercised"
    );

    assert_parity(&run_rs(&fx, &[Chem::Taps.meth_arg()]), &cli, "--meth=taps");
}

/// The CLI accepts upstream's other spellings of the default chemistry.
#[test]
fn emseq_aliases_match_bare_meth() {
    let Some(bwa) = skip_without_tools() else {
        return;
    };
    let fx = fixture(&bwa, Chem::EmSeq, 0x5eed_0003);
    let bare = run_rs(&fx, &["--meth"]);
    for spelling in ["--meth=emseq", "--meth=em-seq", "--meth=bisulfite"] {
        assert_eq!(keys(&run_rs(&fx, &[spelling])), keys(&bare), "{spelling}");
    }
}
