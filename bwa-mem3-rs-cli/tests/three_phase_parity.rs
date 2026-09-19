//! Byte parity of the three-phase API against `bwa-mem3 mem -t 1 -K <K> -p`
//! across several `-K` cohorts, arbitrary sub-batch sizes, several threads,
//! and single-end / paired / mixed inputs.
//!
//! The test reproduces, in plain Rust, the two rules the CLI applies outside
//! the aligner: the cohort cut (`fast_reader_bseq.c:129-137`: after each read,
//! cut when `size >= K && n_reads % 2 == 0`) and the `-p` classification + id
//! bases (`bwa.cpp:258-274`, `fastmap.cpp:924-944`, `bwamem.cpp:2795-2884`).
//! Anything that consumes this crate (fgumi) ports these two functions.
//!
//! Two dedicated single-end branches are pinned here that no other fixture in
//! this crate reaches: a chimeric SE read that produces a `0x800`
//! supplementary record, and a genuinely-unmappable SE read that produces a
//! single `0x4` unmapped record (see `single_end_supplementary_and_unmapped_branches`).
//!
//! Requires `bwa-mem3` + `samtools`; skips otherwise.

mod common;
mod phix_seq;

use std::path::Path;
use std::process::Command;

use bwa_mem3_rs::{
    pair_emit, seed_extend, AlignScratch, AlnRegs, BwaIndex, IdBases, MemOpts, MemPeStat,
    ReadBatch, ReadPair, RecordOrigin, RecordSink, SingleRead,
};
use rayon::prelude::*;
use rstest::{fixture, rstest};

/// One FASTQ record in input order.
#[derive(Clone)]
struct Read {
    name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

/// Reference cohort cut: indices `[start, end)` into the read list.
fn cut_cohorts(reads: &[Read], k: u64) -> Vec<std::ops::Range<usize>> {
    let mut out = Vec::new();
    let (mut start, mut size, mut n) = (0usize, 0u64, 0usize);
    for (i, r) in reads.iter().enumerate() {
        size += r.seq.len() as u64;
        n += 1;
        if size >= k && n % 2 == 0 {
            out.push(start..i + 1);
            start = i + 1;
            size = 0;
            n = 0;
        }
    }
    if start < reads.len() {
        out.push(start..reads.len());
    }
    out
}

/// `bseq_classify`: consecutive same-name reads pair; the rest are singles.
/// Returns (singles as read indices, pairs as (r1 idx, r2 idx)), in order.
fn classify(reads: &[Read], range: std::ops::Range<usize>) -> (Vec<usize>, Vec<(usize, usize)>) {
    let (mut singles, mut pairs) = (Vec::new(), Vec::new());
    let mut i = range.start;
    while i < range.end {
        if i + 1 < range.end && reads[i].name == reads[i + 1].name {
            pairs.push((i, i + 1));
            i += 2;
        } else {
            singles.push(i);
            i += 1;
        }
    }
    (singles, pairs)
}

/// Id bases for a sub-batch: `cohort_base` = global reads before the cohort,
/// `n_se` = singles in the cohort, offsets = singles/pairs before this batch.
fn ids_for(cohort_base: u64, n_se: u64, se_offset: u64, pe_offset: u64) -> IdBases {
    IdBases {
        first_single_id: cohort_base + se_offset,
        first_pair_id: ((cohort_base + n_se) >> 1) + pe_offset,
    }
}

struct Collector<'a> {
    out: &'a mut Vec<(usize, Vec<u8>)>, // (global read index of origin's first read, body)
    single_idx: &'a [usize],
    pair_idx: &'a [(usize, usize)],
}
impl RecordSink for Collector<'_> {
    fn emit(&mut self, origin: RecordOrigin, body: &[u8]) {
        let key = match origin {
            RecordOrigin::Pair(i) => self.pair_idx[i].0,
            RecordOrigin::Single(i) => self.single_idx[i],
        };
        self.out.push((key, body.to_vec()));
    }
}

/// One template of a `-p` batch: a pair (as `(r1_idx, r2_idx)`) or a single
/// (as its read index) into the input read list. Exactly one side is `Some`.
type Template = (Option<(usize, usize)>, Option<usize>);

/// A Phase 3 emit job: one batch's seeded regs, its id bases, and the
/// batch-local single/pair indices used to key emitted records to input reads.
type Job = (AlnRegs, IdBases, Vec<usize>, Vec<(usize, usize)>);

/// Run every cohort through the three phases with `sub` templates per batch on
/// `threads` threads. Records come back keyed by input position so they can
/// be sorted into input order regardless of dispatch.
fn run_three_phase(
    idx: &BwaIndex,
    opts: &MemOpts,
    reads: &[Read],
    k: u64,
    sub: usize,
    threads: usize,
    // When true, force every cohort's id base to 0 instead of the running global
    // read count -- the classic "per-cohort id reset" bug. Correct output uses
    // false; the tandem-repeat test flips this on to prove the cohort-exact
    // `IdBases` is load-bearing (a wrong base changes the id-seeded tie-break in
    // every cohort after the first, so byte parity against the CLI breaks).
    zero_cohort_base: bool,
) -> Vec<Vec<u8>> {
    // A dedicated rayon pool so `threads` bounds the fan-out exactly as the old
    // hand-rolled scoped-thread split did (threads=1 => fully sequential); both
    // parallel phases run inside it. `map_init` gives each worker its own
    // reused `AlignScratch`, and `collect` on the indexed parallel iterators
    // preserves batch order, so no manual index/sort bookkeeping is needed.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .unwrap();
    let mut all: Vec<(usize, Vec<u8>)> = Vec::new();
    let mut cohort_base = 0u64;
    for range in cut_cohorts(reads, k) {
        let (singles, pairs) = classify(reads, range.clone());
        // Sub-batches: consecutive templates (a template = one pair or one single), in input order.
        // A batch mixes kinds only where the input does.
        let mut templates: Vec<Template> = Vec::new();
        let (mut si, mut pi) = (0, 0);
        let mut pos = range.start;
        while pos < range.end {
            if pi < pairs.len() && pairs[pi].0 == pos {
                templates.push((Some(pairs[pi]), None));
                pos += 2;
                pi += 1;
            } else {
                templates.push((None, Some(singles[si])));
                pos += 1;
                si += 1;
            }
        }
        let batches: Vec<&[Template]> = templates.chunks(sub.max(1)).collect();
        // Phase 1: seed+extend every batch in parallel; ordered `collect`.
        let regs: Vec<AlnRegs> = pool.install(|| {
            batches
                .par_iter()
                .map_init(
                    || AlignScratch::new().unwrap(),
                    |sc, b| {
                        let pv: Vec<ReadPair<'_>> = b
                            .iter()
                            .filter_map(|(p, _)| *p)
                            .map(|(a, c)| ReadPair {
                                name_r1: reads[a].name.as_bytes(),
                                seq_r1: &reads[a].seq,
                                qual_r1: Some(&reads[a].qual),
                                name_r2: reads[c].name.as_bytes(),
                                seq_r2: &reads[c].seq,
                                qual_r2: Some(&reads[c].qual),
                            })
                            .collect();
                        let sv: Vec<SingleRead<'_>> = b
                            .iter()
                            .filter_map(|(_, s)| *s)
                            .map(|a| SingleRead {
                                name: reads[a].name.as_bytes(),
                                seq: &reads[a].seq,
                                qual: Some(&reads[a].qual),
                            })
                            .collect();
                        seed_extend(
                            idx,
                            opts,
                            sc,
                            &ReadBatch {
                                pairs: &pv,
                                singles: &sv,
                            },
                        )
                        .unwrap()
                    },
                )
                .collect()
        });
        // Cohort pestat over every batch, in order.
        let pestat = MemPeStat::infer_cohort(idx, opts, &regs).unwrap();
        // Per-batch id bases + kind indices, in batch order.
        let n_se = singles.len() as u64;
        let mut se_off = 0u64;
        let mut pe_off = 0u64;
        let mut jobs: Vec<Job> = Vec::new();
        let base = if zero_cohort_base { 0 } else { cohort_base };
        for (b, r) in batches.iter().zip(regs) {
            let ids = ids_for(base, n_se, se_off, pe_off);
            let bs: Vec<usize> = b.iter().filter_map(|(_, s)| *s).collect();
            let bp: Vec<(usize, usize)> = b.iter().filter_map(|(p, _)| *p).collect();
            se_off += bs.len() as u64;
            pe_off += bp.len() as u64;
            jobs.push((r, ids, bs, bp));
        }
        // Phase 3: emit every batch in parallel; ids per batch from the cohort
        // layout. Ordered `collect` keeps batches in input order.
        let pestat_ref = &pestat;
        let emitted: Vec<Vec<(usize, Vec<u8>)>> = pool.install(|| {
            jobs.into_par_iter()
                .map_init(
                    || AlignScratch::new().unwrap(),
                    |sc, (r, ids, bs, bp)| {
                        let mut out = Vec::new();
                        let has_pairs = r.n_pairs() > 0;
                        let mut sink = Collector {
                            out: &mut out,
                            single_idx: &bs,
                            pair_idx: &bp,
                        };
                        pair_emit(
                            idx,
                            opts,
                            sc,
                            r,
                            has_pairs.then_some(pestat_ref),
                            ids,
                            &mut sink,
                        )
                        .unwrap();
                        out
                    },
                )
                .collect()
        });
        all.extend(emitted.into_iter().flatten());
        cohort_base += (range.end - range.start) as u64;
    }
    // Stable sort by origin read index restores input order; within a template
    // the emitter's order (R1 then R2, primary then supplementary) is kept.
    all.sort_by_key(|(k, _)| *k);
    all.into_iter().map(|(_, b)| b).collect()
}

fn cli_records(bwa: &str, ref_fa: &Path, fq: &Path, k: u64, dir: &Path) -> Vec<String> {
    let out = Command::new(bwa)
        .args(["mem", "-t", "1", "-K", &k.to_string(), "-p"])
        .arg(ref_fa)
        .arg(fq)
        .output()
        .expect("run bwa-mem3 mem");
    assert!(
        out.status.success(),
        "bwa-mem3 mem failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let bam = dir.join("cli.bam");
    std::fs::write(&bam, &out.stdout).unwrap();
    common::samtools_view(&bam)
}

fn fixture_paired(n: usize, seed: u64) -> Vec<Read> {
    let mut v = Vec::new();
    for (name, r1, r2) in common::simulate_pairs(phix_seq::PHIX_SEQ.as_bytes(), n, 150, 400, seed) {
        v.push(Read {
            name: name.clone(),
            seq: r1,
            qual: vec![b'I'; 150],
        });
        v.push(Read {
            name,
            seq: r2,
            qual: vec![b'I'; 150],
        });
    }
    v
}

fn fixture_single(n: usize, seed: u64) -> Vec<Read> {
    let reference = phix_seq::PHIX_SEQ.as_bytes();
    let mut rng = common::Rng(seed);
    (0..n)
        .map(|i| {
            let start = (rng.next() as usize) % (reference.len() - 120);
            let mut seq = reference[start..start + 120].to_vec();
            if rng.next() % 2 == 0 {
                seq = common::revcomp(&seq);
            }
            Read {
                name: format!("s{i}"),
                seq,
                qual: vec![b'I'; 120],
            }
        })
        .collect()
}

/// Pairs with a single inserted after every 4th pair, plus one junk single
/// that maps nowhere, plus an odd read count so an even-parity cut can land
/// inside a pair (the CLI then classifies that pair as two singles).
fn fixture_mixed(seed: u64) -> Vec<Read> {
    let pairs = fixture_paired(600, seed);
    let singles = fixture_single(150, seed + 1);
    let mut out = Vec::new();
    let mut si = singles.into_iter();
    for (i, chunk) in pairs.chunks(2).enumerate() {
        out.extend_from_slice(chunk);
        if i % 4 == 3 {
            if let Some(s) = si.next() {
                out.push(s);
            }
        }
    }
    let mut rng = common::Rng(seed + 7);
    out.push(Read {
        name: "junk".into(),
        seq: common::random_dna(&mut rng, 120).into_bytes(),
        qual: vec![b'I'; 120],
    });
    out
}

/// Single-end reads exercising the two SE `mem_reg2sam` branches that no
/// other fixture in this crate reaches (see the Task 5 review carry-over):
/// `chimera0` is two 60bp arms drawn from PhiX loci 3kb apart, so neither arm
/// is explainable as one gapped alignment and bwa-mem3 emits a primary plus a
/// `0x800` supplementary; `junk0` is 120bp of random DNA with no seedable
/// match anywhere in PhiX's 5.5kb genome, so it fails `-T 30` and emits a
/// single `0x4` unmapped record. A handful of ordinary reads pad the cohort
/// so the run has realistic depth.
fn fixture_se_edge_cases(seed: u64) -> Vec<Read> {
    let reference = phix_seq::PHIX_SEQ.as_bytes();
    let mut rng = common::Rng(seed);
    let mut out: Vec<Read> = (0..30)
        .map(|i| {
            let start = (rng.next() as usize) % (reference.len() - 120);
            let mut seq = reference[start..start + 120].to_vec();
            if rng.next() % 2 == 0 {
                seq = common::revcomp(&seq);
            }
            Read {
                name: format!("ok{i}"),
                seq,
                qual: vec![b'I'; 120],
            }
        })
        .collect();

    let (arm_a, arm_b, arm) = (200usize, 3200usize, 60usize);
    let mut chimera = reference[arm_a..arm_a + arm].to_vec();
    chimera.extend_from_slice(&reference[arm_b..arm_b + arm]);
    out.push(Read {
        name: "chimera0".into(),
        seq: chimera,
        qual: vec![b'I'; 2 * arm],
    });

    out.push(Read {
        name: "junk0".into(),
        seq: common::random_dna(&mut rng, 120).into_bytes(),
        qual: vec![b'I'; 120],
    });
    out
}

/// Runs the three-phase path and the CLI over `reads`, asserts every record's
/// raw SAM line is byte-identical between the two (full-line `assert_eq!`,
/// not just `common::record_key_fields`), and returns the SAM lines of both
/// so a caller can additionally inspect specific records by name. Returns
/// `None` (asserting nothing) when the required tools are missing.
/// The reference index, options, and CLI SAM lines for a parity run: build an
/// index from `ref_seq` (contig `contig`), write the interleaved cohort FASTQ,
/// and align it with `bwa-mem3 mem -p`. The `TempDir` is returned so callers
/// keep the on-disk index + a scratch dir alive. `None` when the tools are
/// missing (or a hard failure under `BWA_MEM3_RS_REQUIRE_TOOLS`).
///
/// Index, options, and CLI SAM lines for a parity run, plus the `TempDir` that
/// keeps the on-disk index and scratch dir alive.
type ParitySetup = (tempfile::TempDir, BwaIndex, MemOpts, Vec<String>);

fn parity_setup(contig: &str, ref_seq: &[u8], reads: &[Read], k: u64) -> Option<ParitySetup> {
    let bwa = common::require_bwa_mem3()?;
    if !common::require_samtools() {
        return None;
    }
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let ref_fa = common::setup_ref_index(dir, &bwa, contig, ref_seq);
    let fq = dir.join("in.fq");
    // `reads` is already in interleaved order (each fixture emits a pair's R1
    // then R2 consecutively), so a plain in-order FASTQ write is exactly what
    // `bwa-mem3 mem -p` consumes -- no separate interleaving step is needed.
    common::write_fastq(
        &fq,
        &reads
            .iter()
            .map(|r| (r.name.clone(), r.seq.clone()))
            .collect::<Vec<_>>(),
    );
    let cli = cli_records(&bwa, &ref_fa, &fq, k, dir);
    let idx = BwaIndex::load(&ref_fa).unwrap();
    let opts = MemOpts::new().unwrap();
    Some((tmp, idx, opts, cli))
}

/// Render three-phase record bodies as `samtools view` SAM lines via a temp BAM.
fn bodies_to_sam(bam: &Path, idx: &BwaIndex, opts: &MemOpts, bodies: &[Vec<u8>]) -> Vec<String> {
    common::write_bam(bam, idx, opts, bodies);
    common::samtools_view(bam)
}

fn check_parity(
    label: &str,
    reads: &[Read],
    k: u64,
    sub: usize,
    threads: usize,
) -> Option<(Vec<String>, Vec<String>)> {
    let (tmp, idx, opts, cli) = parity_setup("phix", phix_seq::PHIX_SEQ.as_bytes(), reads, k)?;
    let bodies = run_three_phase(&idx, &opts, reads, k, sub, threads, false);
    let rs = bodies_to_sam(&tmp.path().join("rs.bam"), &idx, &opts, &bodies);

    let n_cohorts = cut_cohorts(reads, k).len();
    assert_eq!(
        rs.len(),
        cli.len(),
        "{label}: record count (K={k}, sub={sub}, threads={threads}, cohorts={n_cohorts})"
    );
    for (i, (a, b)) in rs.iter().zip(&cli).enumerate() {
        assert_eq!(
            a, b,
            "{label}: record {i} differs (K={k}, sub={sub}, threads={threads}, cohorts={n_cohorts})\n rs: {a}\ncli: {b}"
        );
    }
    Some((rs, cli))
}

/// SAM lines whose QNAME (first column) equals `qname`.
fn named<'a>(lines: &'a [String], qname: &str) -> Vec<&'a String> {
    lines
        .iter()
        .filter(|l| l.split('\t').next() == Some(qname))
        .collect()
}

/// A SAM line's FLAG column.
fn flag(line: &str) -> u32 {
    line.split('\t').nth(1).unwrap().parse().unwrap()
}

/// A SAM line's MAPQ column.
fn mapq(line: &str) -> u32 {
    line.split('\t').nth(4).unwrap().parse().unwrap()
}

#[fixture]
fn paired_reads() -> Vec<Read> {
    fixture_paired(2000, 42) // 600 kb of bases
}

#[rstest]
fn paired_multi_cohort_matrix(
    paired_reads: Vec<Read>,
    #[values(200_000u64, 5_000_000)] k: u64,
    #[values(1usize, 7, 64, 256, 4096)] sub: usize,
    #[values(1usize, 4)] threads: usize,
) {
    check_parity("paired", &paired_reads, k, sub, threads);
}

#[fixture]
fn single_reads() -> Vec<Read> {
    fixture_single(1500, 77)
}

#[rstest]
fn single_end_multi_cohort_matrix(
    single_reads: Vec<Read>,
    #[values(50_000u64, 5_000_000)] k: u64,
    #[values(1usize, 13, 512)] sub: usize,
) {
    check_parity("single", &single_reads, k, sub, 4);
}

/// `fixture_mixed` plus the fixture-validity guard: `k = 30_000` must produce
/// at least one even-parity cohort cut that lands inside a pair (the CLI then
/// classifies that split pair as two singles), or the mid-pair classification
/// path goes unexercised.
#[fixture]
fn mixed_reads() -> Vec<Read> {
    let reads = fixture_mixed(99);
    let cohorts = cut_cohorts(&reads, 30_000);
    let split_pair = cohorts
        .iter()
        .any(|c| c.end < reads.len() && reads[c.end - 1].name == reads[c.end].name);
    assert!(
        split_pair,
        "fixture must produce a cohort cut inside a pair"
    );
    reads
}

#[rstest]
fn mixed_multi_cohort_matrix(
    mixed_reads: Vec<Read>,
    #[values(1usize, 5, 64, 2048)] sub: usize,
    #[values(1usize, 3)] threads: usize,
) {
    // K chosen so at least one even-parity cut falls between the two reads of
    // a pair (asserted in the `mixed_reads` fixture above).
    check_parity("mixed", &mixed_reads, 30_000, sub, threads);
}

/// Carry-over from the Task 5 review: two SE `mem_reg2sam` branches were
/// never pinned byte-for-byte because no other fixture in this crate produces
/// them (PhiX unique mappers only take the "one aligned record" path). This
/// pins both: a chimeric SE read's `0x800` supplementary, and an unmappable SE
/// read's single `0x4` record, across a small slice of the K/sub/threads
/// space (the `paired`/`single`/`mixed` matrices above already cover the
/// dimension exhaustively; this fixture only needs to prove the two branches
/// specifically survive it).
#[rstest]
fn single_end_supplementary_and_unmapped_branches(
    #[values(50_000u64, 5_000_000)] k: u64,
    #[values(1usize, 7)] sub: usize,
    #[values(1usize, 4)] threads: usize,
) {
    let reads = fixture_se_edge_cases(0xC417_0001);
    let Some((rs, cli)) = check_parity("se_edge_cases", &reads, k, sub, threads) else {
        return;
    };

    // Fixture validity, asserted against the reference aligner rather than
    // assumed: without a genuine supplementary and a genuine unmapped record,
    // the branch-specific assertions below would be vacuous.
    let cli_chimera = named(&cli, "chimera0");
    let cli_supp: Vec<&&String> = cli_chimera
        .iter()
        .filter(|l| flag(l) & 0x800 != 0)
        .collect();
    assert_eq!(
        cli_supp.len(),
        1,
        "fixture must produce exactly one SE supplementary record from the \
         reference aligner (K={k}, sub={sub}, threads={threads}); got {} of {} \
         chimera0 records",
        cli_supp.len(),
        cli_chimera.len()
    );

    let cli_junk = named(&cli, "junk0");
    assert_eq!(
        cli_junk.len(),
        1,
        "fixture must produce exactly one junk0 record from the reference \
         aligner (K={k}, sub={sub}, threads={threads})"
    );
    assert_eq!(
        flag(cli_junk[0]) & 0x4,
        0x4,
        "fixture's junk0 read must be genuinely unmapped by the reference \
         aligner (K={k}, sub={sub}, threads={threads})"
    );

    // The branches themselves, on the three-phase path.
    let rs_chimera = named(&rs, "chimera0");
    let rs_supp: Vec<&&String> = rs_chimera.iter().filter(|l| flag(l) & 0x800 != 0).collect();
    assert_eq!(
        rs_supp.len(),
        1,
        "three-phase path must reproduce the SE supplementary branch \
         (K={k}, sub={sub}, threads={threads}); rs chimera0 records: {rs_chimera:?}"
    );

    let rs_junk = named(&rs, "junk0");
    assert_eq!(
        rs_junk.len(),
        1,
        "three-phase path must emit exactly one junk0 record \
         (K={k}, sub={sub}, threads={threads})"
    );
    assert_eq!(
        flag(rs_junk[0]) & 0x4,
        0x4,
        "three-phase path must reproduce the SE unmapped branch \
         (K={k}, sub={sub}, threads={threads})"
    );
}

// ---------------------------------------------------------------------------
// Tandem-repeat parity: the id-seeded tie-break, exercised where it MATTERS.
//
// Every other fixture in this crate draws reads from repeat-free PhiX, where
// each read maps to exactly one locus. There the global pair/single id feeds
// the tie-break hash but changes no output, so a wrong `IdBases` -- even one
// that ignored the cohort base entirely -- would still pass byte parity. These
// tests draw reads from within a single unit of a TANDEM-REPEAT reference, so
// every read maps to `TANDEM_COPIES` equal-score loci and the id-seeded hash
// (`hash_64(id<<1|i)`, `mem_mark_primary_se`/`mem_pair`) alone decides which
// locus is primary and the XA order of the rest. Reads span >= 2 cohorts so a
// NON-ZERO cohort base is actually applied, and both the paired and the
// SE/mixed emit paths are covered.
// ---------------------------------------------------------------------------

/// The 300 bp PhiX window used as the tandem repeat unit (same window the sys
/// crate's `first_pair_id_reaches_the_tie_break_hash` uses).
const TANDEM_UNIT: std::ops::Range<usize> = 500..800;
/// Number of times the unit is repeated in the reference.
const TANDEM_COPIES: usize = 5;

fn tandem_unit() -> &'static [u8] {
    &phix_seq::PHIX_SEQ.as_bytes()[TANDEM_UNIT]
}

fn tandem_reference() -> Vec<u8> {
    tandem_unit().repeat(TANDEM_COPIES)
}

/// FR pairs drawn from WITHIN one repeat unit, so each pair maps to all
/// `TANDEM_COPIES` copies at equal score.
fn fixture_tandem_paired(n: usize, seed: u64) -> Vec<Read> {
    let unit = tandem_unit();
    let (read_len, insert) = (60usize, 120usize);
    let mut rng = common::Rng(seed);
    let mut v = Vec::with_capacity(2 * n);
    for i in 0..n {
        let start = (rng.next() as usize) % (unit.len() - insert + 1);
        let r1 = unit[start..start + read_len].to_vec();
        let r2 = common::revcomp(&unit[start + insert - read_len..start + insert]);
        v.push(Read {
            name: format!("t{i}"),
            seq: r1,
            qual: vec![b'I'; read_len],
        });
        v.push(Read {
            name: format!("t{i}"),
            seq: r2,
            qual: vec![b'I'; read_len],
        });
    }
    v
}

/// Tandem pairs with a single (also drawn from within one unit, so it too
/// multi-maps) inserted after every third pair, so the SE emit path runs
/// alongside the paired one in the same multi-cohort layout.
fn fixture_tandem_mixed(n_pairs: usize, seed: u64) -> Vec<Read> {
    let pairs = fixture_tandem_paired(n_pairs, seed);
    let unit = tandem_unit();
    let mut rng = common::Rng(seed + 1);
    let mut out = Vec::new();
    for (i, chunk) in pairs.chunks(2).enumerate() {
        out.extend_from_slice(chunk);
        if i % 3 == 2 {
            let start = (rng.next() as usize) % (unit.len() - 60 + 1);
            let mut seq = unit[start..start + 60].to_vec();
            if rng.next() % 2 == 0 {
                seq = common::revcomp(&seq);
            }
            out.push(Read {
                name: format!("ts{i}"),
                seq,
                qual: vec![b'I'; 60],
            });
        }
    }
    out
}

/// Parity on the tandem-repeat reference, asserting the id is genuinely
/// load-bearing so the parity check is not vacuous:
///   1. the fixture spans >= 2 cohorts (so a non-zero cohort base is applied)
///      and the CLI genuinely multi-maps it (MAPQ-0 mapped primaries exist);
///   2. the correct cohort-exact `IdBases` reproduce the CLI byte-for-byte; and
///   3. forcing every cohort id base to 0 (the classic per-cohort-reset bug)
///      DIVERGES from the CLI -- direct proof that the id decides the output on
///      this fixture, i.e. (2) would fail if `IdBases` carried the wrong base.
fn check_tie_break_parity(label: &str, reads: &[Read], k: u64, sub: usize, threads: usize) {
    let refseq = tandem_reference();
    let Some((tmp, idx, opts, cli)) = parity_setup("tandem", &refseq, reads, k) else {
        return;
    };
    let dir = tmp.path();

    // (1) fixture validity, asserted rather than assumed.
    let cohorts = cut_cohorts(reads, k).len();
    assert!(
        cohorts >= 2,
        "{label}: need >= 2 cohorts so a non-zero cohort base is exercised \
         (K={k}); got {cohorts}"
    );
    let mapped_multi = cli
        .iter()
        .filter(|l| flag(l) & 0x4 == 0 && mapq(l) == 0)
        .count();
    assert!(
        mapped_multi > 0,
        "{label}: fixture must genuinely multi-map (the CLI produced no mapped \
         MAPQ-0 record, so the tie-break id would be inert)"
    );

    // (2) correct cohort-exact ids must reproduce the CLI byte-for-byte.
    let good = run_three_phase(&idx, &opts, reads, k, sub, threads, false);
    let good_sam = bodies_to_sam(&dir.join("good.bam"), &idx, &opts, &good);
    assert_eq!(
        good_sam.len(),
        cli.len(),
        "{label}: record count (K={k}, sub={sub}, threads={threads}, cohorts={cohorts})"
    );
    for (i, (a, b)) in good_sam.iter().zip(&cli).enumerate() {
        assert_eq!(
            a, b,
            "{label}: record {i} differs (K={k}, sub={sub}, threads={threads}, cohorts={cohorts})\n rs: {a}\ncli: {b}"
        );
    }

    // (3) the broken per-cohort-reset ids must DIVERGE, proving (2) is not vacuous.
    let broken = run_three_phase(&idx, &opts, reads, k, sub, threads, true);
    let broken_sam = bodies_to_sam(&dir.join("broken.bam"), &idx, &opts, &broken);
    assert_ne!(
        broken_sam, cli,
        "{label}: forcing every cohort id base to 0 still matched the CLI, so \
         the cohort-exact IdBases is not observed by the tie-break and the \
         parity check above is vacuous (K={k}, sub={sub}, threads={threads}, cohorts={cohorts})"
    );
}

/// Paired path: the id-seeded tie-break decides each pair's primary locus.
#[rstest]
fn paired_tandem_repeat_tie_break_matches_cli(
    #[values(1usize, 64)] sub: usize,
    #[values(1usize, 4)] threads: usize,
) {
    // 300 pairs of 60 bp = 36 kb; K = 12 kb cuts ~3 cohorts (asserted >= 2).
    let reads = fixture_tandem_paired(300, 0xB165_EED5);
    check_tie_break_parity("tandem-paired", &reads, 12_000, sub, threads);
}

/// SE/mixed path: singles share the multi-cohort layout with the pairs, so the
/// single id's tie-break (`mem_mark_primary_se`) is exercised too.
#[rstest]
fn mixed_tandem_repeat_tie_break_matches_cli(
    #[values(1usize, 64)] sub: usize,
    #[values(1usize, 4)] threads: usize,
) {
    let reads = fixture_tandem_mixed(300, 0x7A4D_E110);
    check_tie_break_parity("tandem-mixed", &reads, 12_000, sub, threads);
}

#[test]
fn cut_rule_matches_fast_reader_on_odd_lengths() {
    // Reads of length 3 with K=7: size hits 9 after the 3rd read (n odd → no
    // cut), 12 after the 4th (n even → cut). Then repeats.
    let reads: Vec<Read> = (0..10)
        .map(|i| Read {
            name: format!("x{i}"),
            seq: b"ACG".to_vec(),
            qual: b"III".to_vec(),
        })
        .collect();
    assert_eq!(cut_cohorts(&reads, 7), vec![0..4, 4..8, 8..10]);
}

#[rstest]
#[case::cohort_start(1000, 3, 0, 0, IdBases { first_single_id: 1000, first_pair_id: 501 })]
#[case::second_batch_with_offsets(1000, 3, 2, 10, IdBases { first_single_id: 1002, first_pair_id: 511 })]
fn id_bases_follow_the_cli_formulas(
    #[case] cohort_base: u64,
    #[case] n_se: u64,
    #[case] se_offset: u64,
    #[case] pe_offset: u64,
    #[case] expected: IdBases,
) {
    assert_eq!(ids_for(cohort_base, n_se, se_offset, pe_offset), expected);
}
