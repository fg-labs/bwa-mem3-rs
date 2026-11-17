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
//! scoring to NEUTRAL; it is checked the same way, as are its `5base` /
//! `5-base` spellings. The em-seq scoring default (GENOMIC as of bwa-mem3
//! v0.17.0, COLLAPSED before) is pinned too: the 8%-diverged repeat copies give
//! reads real mismatches, which COLLAPSED's `-B 2` and GENOMIC's `-B 4` score
//! differently.
//!
//! `--md short|none` (v0.17.0) and `--meth-context cpg` (v0.16.0) are
//! checked the same way, the latter on a TAPS fixture with off-target
//! conversion so the context changes output.
//!
//! The comparison adds `AS:i`, `XS:i` and `HN:i` to [`record_key_fields`], since
//! most of the records pruning changes differ only in `XS`, and a scoring change
//! shows first in `AS`.
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
    /// TAPS plus off-target conversion (damage) at 5% of non-CpG Cs, which
    /// `--meth-context cpg` scores as mismatches and `all` does not.
    TapsDamaged,
}

impl Chem {
    fn meth_arg(self) -> &'static str {
        match self {
            Chem::EmSeq => "--meth",
            Chem::Taps | Chem::TapsDamaged => "--meth=taps",
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
                Chem::TapsDamaged => {
                    if cpg {
                        rng.next() % 100 < 80
                    } else {
                        rng.next() % 100 < 5
                    }
                }
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

/// [`record_key_fields`] plus the `AS:i`, `XS:i` and `HN:i` values.
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
        "{}\tAS:{}\tXS:{}\tHN:{}",
        record_key_fields(sam_line),
        tag("AS:i:"),
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
    // These tests pass flags and expect defaults added up to the vendored
    // version (`--md`, the em-seq scoring default), so an older oracle is
    // reported up front rather than as a generic `run_cli` failure.
    let bwa = common::require_vendored_bwa_mem3()?;
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

/// bwa-mem3 v0.16.0 added `--meth=5base` (Illumina 5-base, which is TAPS
/// chemistry as far as the aligner is concerned) with a `5-base` alias. Each
/// spelling matches both `bwa-rs mem --meth=taps` and the CLI's own `5base`.
#[test]
fn five_base_spellings_match_taps() {
    let Some(bwa) = skip_without_tools() else {
        return;
    };
    let fx = fixture(&bwa, Chem::Taps, 0x5eed_0004);
    let taps = run_rs(&fx, &[Chem::Taps.meth_arg()]);
    // Non-vacuous: TAPS and em-seq differ on this fixture, so a spelling the
    // parser mapped to em-seq could not match TAPS below.
    assert!(
        differing(&taps, &run_rs(&fx, &["--meth"])) > 0,
        "TAPS and em-seq agree on this fixture; it cannot tell the spellings apart"
    );
    for spelling in ["--meth=5base", "--meth=5-base"] {
        assert_eq!(keys(&run_rs(&fx, &[spelling])), keys(&taps), "{spelling}");
    }
    let cli = run_cli(&bwa, &fx, &["--meth=5base"], &[], "5base");
    assert_parity(&taps, &cli, "--meth=5base");
}

/// bwa-mem3 v0.17.0 (#756) made GENOMIC the em-seq `--meth` scoring default in
/// place of COLLAPSED. `bwa-rs mem --meth` inherits it from `mem_opt_init`; this
/// pins that both tools agree on a fixture where the two modes differ.
#[test]
fn emseq_scoring_default_matches_cli() {
    let Some(bwa) = skip_without_tools() else {
        return;
    };
    let fx = fixture(&bwa, Chem::EmSeq, 0x5eed_0005);
    let cli = run_cli(&bwa, &fx, &["--meth"], &[], "default");

    // Non-vacuous: on this fixture COLLAPSED scoring changes the CLI's output,
    // so a crate still defaulting to it could not pass the parity check below.
    let collapsed = run_cli(
        &bwa,
        &fx,
        &["--meth", "--meth-scoring", "collapsed"],
        &[],
        "collapsed",
    );
    assert!(
        differing(&cli, &collapsed) > 0,
        "--meth-scoring collapsed changed no CLI record; the fixture no longer tests the default"
    );

    assert_parity(&run_rs(&fx, &["--meth"]), &cli, "--meth (scoring default)");
}

/// `--md short` (only the mismatches NM counts, which under `--meth` leaves
/// out the conversions) and `--md none` (no MD tag), bwa-mem3 v0.17.0 (#754).
#[test]
fn md_modes_match_cli() {
    let Some(bwa) = skip_without_tools() else {
        return;
    };
    let fx = fixture(&bwa, Chem::EmSeq, 0x5eed_0006);
    let full = run_cli(&bwa, &fx, &["--meth"], &[], "md_full");
    for mode in ["short", "none"] {
        let cli = run_cli(&bwa, &fx, &["--meth", "--md", mode], &[], mode);
        assert!(
            differing(&cli, &full) > 0,
            "--md {mode} changed no CLI record; the fixture no longer tests it"
        );
        assert_parity(&run_rs(&fx, &["--meth", "--md", mode]), &cli, mode);
    }
}

/// `--meth-context cpg` (bwa-mem3 v0.16.0, #663) frees a conversion only at a
/// CpG. The fixture adds off-target (non-CpG) conversion so it changes output.
#[test]
fn meth_context_cpg_matches_cli() {
    let Some(bwa) = skip_without_tools() else {
        return;
    };
    let fx = fixture(&bwa, Chem::TapsDamaged, 0x5eed_0007);
    let args = ["--meth=taps", "--meth-context", "cpg"];
    let cli = run_cli(&bwa, &fx, &args, &[], "cpg");
    let all = run_cli(&bwa, &fx, &["--meth=taps"], &[], "all");
    assert!(
        differing(&cli, &all) > 0,
        "--meth-context cpg changed no CLI record; the fixture no longer tests it"
    );
    assert_parity(&run_rs(&fx, &args), &cli, "--meth=taps --meth-context cpg");
}

/// `bwa-mem3 mem` refuses `--meth-context cpg` with em-seq chemistry, and so
/// does `bwa-rs mem` -- with an error, not upstream's exit(1) from inside the
/// library.
#[test]
fn meth_context_cpg_with_emseq_is_refused() {
    let Some(bwa) = skip_without_tools() else {
        return;
    };
    let fx = fixture(&bwa, Chem::EmSeq, 0x5eed_0008);
    let cli = Command::new(&bwa)
        .args(["mem", "--meth", "--meth-context", "cpg"])
        .arg(&fx.ref_fa)
        .arg(&fx.r1)
        .arg(&fx.r2)
        .output()
        .expect("run bwa-mem3 mem");
    assert!(!cli.status.success(), "the CLI accepted cpg with em-seq");
    let cli_err = String::from_utf8_lossy(&cli.stderr);
    assert!(
        cli_err.contains("--meth-context cpg requires --meth=taps"),
        "the CLI failed for another reason: {cli_err}"
    );
    let rs = Command::new(common::cli_bin())
        .args(["mem", "--meth", "--meth-context", "cpg"])
        .arg(&fx.ref_fa)
        .arg(&fx.r1)
        .arg(&fx.r2)
        .arg("-o")
        .arg(fx.dir.join("refused.bam"))
        .output()
        .expect("run bwa-rs mem");
    assert!(!rs.status.success(), "bwa-rs accepted cpg with em-seq");
    let stderr = String::from_utf8_lossy(&rs.stderr);
    // The library's error (anyhow's `Error: invalid options: ...`), not
    // upstream's `ERROR: ...` printed on the way to exit(1).
    assert!(
        stderr.contains("invalid options: --meth-context cpg requires --meth=taps"),
        "unexpected error: {stderr}"
    );
    assert!(!stderr.contains("ERROR: --meth-context"), "{stderr}");
}

/// The options are validated before the index is loaded, as `bwa-mem3 mem`
/// does: an invalid combination fails at once rather than after a multi-GB
/// load. A prefix with no index on disk tells the two orders apart (loading
/// first would fail on the load instead), and needs neither bwa-mem3 nor
/// samtools.
#[test]
fn invalid_options_fail_before_the_index_is_loaded() {
    let tmp = tempfile::tempdir().unwrap();
    let missing = tmp.path().join("no-such-index.fa");
    let rs = Command::new(common::cli_bin())
        .args(["mem", "--meth", "--meth-context", "cpg"])
        .arg(&missing)
        .arg(tmp.path().join("r1.fq"))
        .arg(tmp.path().join("r2.fq"))
        .arg("-o")
        .arg(tmp.path().join("out.bam"))
        .output()
        .expect("run bwa-rs mem");
    assert!(!rs.status.success(), "bwa-rs accepted cpg with em-seq");
    let stderr = String::from_utf8_lossy(&rs.stderr);
    assert!(
        stderr.contains("invalid options: --meth-context cpg requires --meth=taps"),
        "expected the options error, got: {stderr}"
    );
    assert!(
        !stderr.contains("loading"),
        "the index was loaded first: {stderr}"
    );
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
