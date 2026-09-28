//! Byte parity of the resident memo against `bwa-mem3 mem -t 1 -K <K> -p
//! --cohort-slices 0 --dedup-reads on|off`: the memo path (PairMemo +
//! seed_extend_with_reps + resolve_memo) must emit exactly the CLI's SAM, and
//! the CLI's own memo must have fired (its duplicate count equals ours).
//! `--cohort-slices 0` is load-bearing: the CLI slices its first batch by
//! default and the slice path disarms the memo, which would make "on == off"
//! vacuous. Requires `bwa-mem3` + `samtools`; skips otherwise.

mod common;
mod phix_seq;

use std::path::Path;
use std::process::Command;

use bwa_mem3_rs::{
    AlignScratch, BwaIndex, MemOpts, MemoStats, PairMemo, ReadPair, ResidentCohort, SingleRead,
};
use common::{classify, cut_cohorts, ids_for, Collector, Read};
use rstest::rstest;

/// A position-dependent QUAL string, so a record carrying another read's QUAL
/// shows.
fn qual(seed: usize, n: usize) -> Vec<u8> {
    (0..n)
        .map(|j| b'#' + ((seed * 11 + j * 7) % 40) as u8)
        .collect()
}

/// Push one pair (R1 then R2 under one name) onto `out`.
fn push_pair(out: &mut Vec<Read>, name: String, r1: &[u8], r2: &[u8], q: usize) {
    out.push(Read {
        name: name.clone(),
        seq: r1.to_vec(),
        qual: qual(q, r1.len()),
    });
    out.push(Read {
        name,
        seq: r2.to_vec(),
        qual: qual(q + 1, r2.len()),
    });
}

/// PhiX pairs with duplicates (new name, different QUAL), triplicates, the
/// near-misses that must not merge (swapped mates, one mate changed, lowercase,
/// an `N`), and a single after every fifth pair. Duplicates land both in the
/// same cohort as their original and in later cohorts.
fn dup_reads(seed: u64) -> Vec<Read> {
    let base = common::simulate_pairs(phix_seq::PHIX_SEQ.as_bytes(), 200, 150, 400, seed);
    let mut out = Vec::new();
    let mut rng = common::Rng(seed + 3);
    for (i, (name, r1, r2)) in base.iter().enumerate() {
        push_pair(&mut out, name.clone(), r1, r2, i);
        if i % 3 == 0 {
            let j = i - (i % 6); // duplicate this or an earlier pair
            push_pair(
                &mut out,
                format!("dup{i}"),
                &base[j].1,
                &base[j].2,
                i + 1000,
            );
        }
        if i % 9 == 0 {
            push_pair(&mut out, format!("tri{i}"), r1, r2, i + 2000);
        }
        if i % 5 == 4 {
            out.push(Read {
                name: format!("s{i}"),
                seq: r1[..120].to_vec(),
                qual: qual(i, 120),
            });
        }
        if i % 40 == 7 {
            push_pair(&mut out, format!("swap{i}"), r2, r1, i + 3000);
            let mut changed = r2.clone();
            changed[20] = if changed[20] == b'A' { b'C' } else { b'A' };
            push_pair(&mut out, format!("one{i}"), r1, &changed, i + 4000);
            push_pair(
                &mut out,
                format!("low{i}"),
                &r1.to_ascii_lowercase(),
                r2,
                i + 5000,
            );
            let mut with_n = r1.clone();
            with_n[30] = b'N';
            push_pair(&mut out, format!("n{i}"), &with_n, r2, i + 6000);
        }
    }
    out.push(Read {
        name: "junk".into(),
        seq: common::random_dna(&mut rng, 120).into_bytes(),
        qual: qual(9, 120),
    });
    out
}

/// The CLI's SAM lines and its reported duplicate-pair count (0 when it
/// printed no stats line).
fn cli_run(
    bwa: &str,
    ref_fa: &Path,
    fq: &Path,
    k: u64,
    dedup: &str,
    dir: &Path,
) -> (Vec<String>, u64) {
    let out = Command::new(bwa)
        .args(["mem", "-t", "1", "-K", &k.to_string(), "-p"])
        .args(["--cohort-slices", "0", "--dedup-reads", dedup])
        .arg(ref_fa)
        .arg(fq)
        .env("BWAMEM3_DEDUP_READS_STATS", "1")
        .output()
        .expect("run bwa-mem3 mem");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "bwa-mem3 mem failed: {stderr}");
    let dups = stderr
        .lines()
        .find_map(|l| l.strip_prefix("[dedup-reads-stats] pairs="))
        .and_then(|rest| rest.split(" dup=").nth(1))
        .and_then(|rest| rest.split_whitespace().next())
        .map_or(0, |d| d.parse().unwrap());
    let bam = dir.join(format!("cli-{dedup}.bam"));
    std::fs::write(&bam, &out.stdout).unwrap();
    (common::samtools_view(&bam), dups)
}

/// Every cohort through a `ResidentCohort` with the memo on, `sub` pairs per
/// range. Returns record bodies in input order and the memo stats summed over
/// the cohorts.
fn run_resident_memo(
    idx: &BwaIndex,
    opts: &MemOpts,
    reads: &[Read],
    k: u64,
    sub: usize,
) -> (Vec<Vec<u8>>, MemoStats) {
    let mut all: Vec<(usize, Vec<u8>)> = Vec::new();
    let mut cohort_base = 0u64;
    let mut total = MemoStats::default();
    let mut scratch = AlignScratch::new().unwrap();
    for range in cut_cohorts(reads, k) {
        let (singles, pairs) = classify(reads, range.clone());
        let cohort = ResidentCohort::new(false).unwrap();
        let mut memo = PairMemo::new();
        let pv: Vec<ReadPair<'_>> = pairs
            .iter()
            .map(|&(a, b)| ReadPair {
                name_r1: reads[a].name.as_bytes(),
                seq_r1: &reads[a].seq,
                qual_r1: Some(&reads[a].qual),
                name_r2: reads[b].name.as_bytes(),
                seq_r2: &reads[b].seq,
                qual_r2: Some(&reads[b].qual),
            })
            .collect();
        let mut ranges = Vec::new();
        for (c, chunk) in pv.chunks(sub.max(1)).enumerate() {
            let mut r = cohort.reserve_pairs(chunk.len()).unwrap();
            cohort.write_pairs(&mut r, chunk).unwrap();
            let reps = memo.mark_range(&r, chunk).unwrap();
            cohort
                .seed_extend_with_reps(idx, opts, &mut scratch, &mut r, &reps)
                .unwrap();
            ranges.push((r, c * sub.max(1), chunk.len()));
        }
        let sv: Vec<SingleRead<'_>> = singles
            .iter()
            .map(|&a| SingleRead {
                name: reads[a].name.as_bytes(),
                seq: &reads[a].seq,
                qual: Some(&reads[a].qual),
            })
            .collect();
        let mut sr = cohort.reserve_singles(sv.len()).unwrap();
        cohort.write_singles(&mut sr, &sv).unwrap();
        cohort
            .seed_extend(idx, opts, &mut scratch, &mut sr)
            .unwrap();
        let stats = cohort.resolve_memo(idx, opts, &mut scratch).unwrap();
        total.dup_pairs += stats.dup_pairs;
        total.copied += stats.copied;
        total.fallback_aligned += stats.fallback_aligned;
        let pestat = cohort.infer_cohort(idx, opts).unwrap();
        let n_se = singles.len() as u64;
        for (mut r, p0, n) in ranges {
            let mut sink = Collector {
                out: &mut all,
                single_idx: &[],
                pair_idx: &pairs[p0..p0 + n],
            };
            let ids = ids_for(cohort_base, n_se, 0, p0 as u64);
            cohort
                .pair_emit(
                    idx,
                    opts,
                    &mut scratch,
                    &mut r,
                    Some(&pestat),
                    ids,
                    0,
                    &mut sink,
                )
                .unwrap();
        }
        let mut sink = Collector {
            out: &mut all,
            single_idx: &singles,
            pair_idx: &[],
        };
        let ids = ids_for(cohort_base, n_se, 0, 0);
        cohort
            .pair_emit(idx, opts, &mut scratch, &mut sr, None, ids, 0, &mut sink)
            .unwrap();
        cohort_base += (range.end - range.start) as u64;
    }
    all.sort_by_key(|(k, _)| *k);
    (all.into_iter().map(|(_, b)| b).collect(), total)
}

#[rstest]
#[case::one_cohort(1_000_000_000, 64)]
#[case::several_cohorts(20_000, 16)]
#[case::tiny_ranges(20_000, 1)]
fn memo_matches_cli_dedup_on_and_off(#[case] k: u64, #[case] sub: usize) {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    if !common::require_samtools() {
        return;
    }
    let reads = dup_reads(11);
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let ref_fa = common::setup_ref_index(dir, &bwa, "phix", phix_seq::PHIX_SEQ.as_bytes());
    let fq = dir.join("in.fq");
    common::write_fastq_reads(&fq, &reads);
    let (cli_on, cli_dups) = cli_run(&bwa, &ref_fa, &fq, k, "on", dir);
    let (cli_off, _) = cli_run(&bwa, &ref_fa, &fq, k, "off", dir);
    assert!(
        cli_dups > 0,
        "the CLI memo never fired (K={k}); is --cohort-slices 0 honored?"
    );
    assert_eq!(cli_on, cli_off, "CLI --dedup-reads on != off (K={k})");

    let idx = BwaIndex::load(&ref_fa).unwrap();
    let opts = MemOpts::new().unwrap();
    let (bodies, stats) = run_resident_memo(&idx, &opts, &reads, k, sub);
    let bam = dir.join("rs.bam");
    common::write_bam(&bam, &idx, &opts, &bodies);
    let rs = common::samtools_view(&bam);
    assert_eq!(
        stats,
        MemoStats {
            dup_pairs: cli_dups,
            copied: cli_dups,
            fallback_aligned: 0
        },
        "memo stats vs the CLI's duplicate count (K={k})"
    );
    assert_eq!(rs.len(), cli_on.len(), "record count (K={k}, sub={sub})");
    for (i, (a, b)) in rs.iter().zip(&cli_on).enumerate() {
        assert_eq!(
            a, b,
            "record {i} differs (K={k}, sub={sub})\n rs: {a}\ncli: {b}"
        );
    }
}
