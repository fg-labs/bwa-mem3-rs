//! CLI parity: align a simulated paired-end fixture with both `bwa-mem3 mem`
//! and `align_batch`, then compare the resulting records **field by field**.
//!
//! Skipped unless `BWA_MEM3_BIN` (path to the CLI), `BWA_MEM3_RS_TEST_REF`
//! (index prefix, which must also be the indexed FASTA and have a `.fai`), and
//! `samtools` on PATH are all available.
//!
//! # Why the fixture is simulated from the reference
//!
//! This test used to align four copies of a 60-mer `ACGT` repeat and assert only
//! that the two record *counts* were equal. That passed without exercising
//! anything: the fixture logged `# candidate unique pairs for (FF, FR, RF, RR):
//! (0, 0, 0, 0)`, i.e. no pair contributed to the insert-size model, so the
//! pairing path this crate exists to reproduce was never entered. Reads are now
//! drawn from the reference itself in FR orientation with substitution noise, so
//! they map uniquely, pair properly, and — being real genomic sequence — some
//! land in repeats and acquire `XA:Z`.
//!
//! # Why both sides must see one insert-size cohort
//!
//! `mem_pestat` runs per batch and its model feeds mate rescue, so two runs that
//! split the same reads into different cohorts legitimately disagree on
//! borderline pairs (see gotcha #14 in `CLAUDE.md`). This test aligns every pair
//! in a single `align_batch` call and runs the CLI at `-t 1` on an input small
//! enough to be one chunk, so both sides estimate one model over the same reads.
//! With that held, parity is exact — every field of every record.

use std::collections::BTreeMap;
use std::io::Write;
use std::process::Command;

use bwa_mem3_rs::{align_batch, BwaIndex, MemOpts, ReadPair};

const N_PAIRS: usize = 2000;
const READ_LEN: usize = 150;
const INSERT: usize = 400;
/// Substitution probability per base. Enough to make MAPQ/NM/MD non-trivial
/// without pushing reads off their true locus.
const SUB_RATE: f64 = 0.01;
const SEED: u64 = 42;

/// Tags compared by value. The remaining tags are compared by key presence
/// only, via the `TAGS:` component of [`Record::compare_key`].
const COMPARED_TAGS: &[&str] = &["NM", "MD", "AS", "XS", "XA", "SA", "MC", "MQ", "HN"];

// ---------------------------------------------------------------------------
// Environment gating
// ---------------------------------------------------------------------------

fn skip() -> Option<(String, String)> {
    let bwa = std::env::var("BWA_MEM3_BIN").ok()?;
    let prefix = std::env::var("BWA_MEM3_RS_TEST_REF").ok()?;
    if !std::path::Path::new(&bwa).exists() {
        eprintln!("skip: BWA_MEM3_BIN={bwa} not found");
        return None;
    }
    if !std::path::Path::new(&format!("{prefix}.bwt.2bit.64")).exists() {
        eprintln!("skip: no bwa-mem3 index at {prefix}");
        return None;
    }
    if !std::path::Path::new(&format!("{prefix}.fai")).exists() {
        eprintln!("skip: no FASTA index at {prefix}.fai (needed to draw reads)");
        return None;
    }
    if Command::new("samtools").arg("--version").output().is_err() {
        eprintln!("skip: samtools not on PATH");
        return None;
    }
    Some((bwa, prefix))
}

// ---------------------------------------------------------------------------
// Deterministic fixture
// ---------------------------------------------------------------------------

/// SplitMix64 — a tiny, fully deterministic PRNG so the fixture is identical on
/// every platform and every run. Test data is generated, never committed.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    fn unit(&mut self) -> f64 {
        // Top 53 bits -> [0, 1); exactly representable in f64.
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Fetch a region of the reference via `samtools faidx`, preserving case.
fn fetch_region(fasta: &str, region: &str) -> Vec<u8> {
    let out = Command::new("samtools")
        .args(["faidx", fasta, region])
        .output()
        .expect("run samtools faidx");
    assert!(
        out.status.success(),
        "samtools faidx {region} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout)
        .lines()
        .skip(1)
        .flat_map(|l| l.bytes())
        .collect()
}

/// Pick a region to draw reads from, by *measuring* candidate windows rather
/// than guessing at one.
///
/// The obvious heuristic — the middle of the longest contig — is wrong on a
/// human reference: the middle of chr1 is the centromere, and a fixture drawn
/// from alpha satellite gives every read a MAPQ-0 multi-hit. `mem_pestat` then
/// finds no unique pairs at all, fails every orientation, and **nothing** is
/// flagged properly-paired, so the pairing path goes untested. (Measured:
/// 0/4000 properly paired, 1037 records with `XA:Z`.)
///
/// So scan candidate offsets across the contig and score each window on two
/// signals available without aligning anything: `N` runs (assembly gaps, which
/// yield unmappable reads) and lower-case bases (RepeatMasker soft-masking,
/// which marks exactly the repeat content that destroys unique mapping).
/// References that are not soft-masked simply score 0 everywhere and fall back
/// to the first gap-free window, which is the old behavior minus the centromere.
fn pick_region(idx: &BwaIndex, fasta: &str) -> Option<(String, Vec<u8>)> {
    const WINDOW: usize = 400_000;
    /// Reject a window with any assembly gap: reads drawn from one are
    /// unmappable and would skew the fixture.
    const MAX_N_FRAC: f64 = 0.0;

    let (name, len) = idx
        .contigs()
        .max_by_key(|&(_, len)| len)
        .map(|(n, l)| (n.to_owned(), l as usize))?;
    if len < 4 * INSERT {
        eprintln!("skip: longest contig {name} is only {len} bp");
        return None;
    }
    let window = WINDOW.min(len / 4).max(4 * INSERT);

    let mut best: Option<(f64, String, Vec<u8>)> = None;
    for pct in [10, 20, 30, 40, 60, 70, 80, 90] {
        let start = len * pct / 100;
        if start + window > len {
            continue;
        }
        // samtools faidx regions are 1-based inclusive.
        let region = format!("{name}:{}-{}", start + 1, start + window);
        let bases = fetch_region(fasta, &region);
        if bases.len() < window / 2 {
            continue;
        }
        let n_frac = bases
            .iter()
            .filter(|b| b.eq_ignore_ascii_case(&b'N'))
            .count() as f64
            / bases.len() as f64;
        if n_frac > MAX_N_FRAC {
            continue;
        }
        let repeat_frac =
            bases.iter().filter(|b| b.is_ascii_lowercase()).count() as f64 / bases.len() as f64;
        // `map_or(true, ..)`, not `is_none_or`: the latter is stable since 1.82
        // and this workspace's MSRV is 1.75.
        #[allow(clippy::unnecessary_map_or)]
        if best.as_ref().map_or(true, |(b, _, _)| repeat_frac < *b) {
            best = Some((repeat_frac, region, bases));
        }
        if repeat_frac < 0.2 {
            break; // good enough; stop paying for faidx calls
        }
    }

    let (repeat_frac, region, bases) = best?;
    eprintln!(
        "fixture region {region} ({:.1}% repeat-masked)",
        repeat_frac * 100.0
    );
    Some((
        region,
        bases.iter().map(|b| b.to_ascii_uppercase()).collect(),
    ))
}

/// Simulate FR-oriented pairs: R1 forward from the fragment start, R2
/// reverse-complement from the fragment end. Fragments containing `N` are
/// skipped so every read has a true locus.
fn simulate_pairs(ref_bases: &[u8]) -> Vec<(String, Vec<u8>, Vec<u8>)> {
    let mut rng = Rng(SEED);
    let mut out = Vec::with_capacity(N_PAIRS);
    let mut attempts = 0usize;
    while out.len() < N_PAIRS {
        attempts += 1;
        assert!(
            attempts < N_PAIRS * 100,
            "could not find {N_PAIRS} N-free fragments in the chosen region"
        );
        let start = rng.below(ref_bases.len() - INSERT);
        let frag = &ref_bases[start..start + INSERT];
        if frag.contains(&b'N') {
            continue;
        }
        let r1 = mutate(&frag[..READ_LEN], &mut rng);
        let r2 = mutate(&revcomp(&frag[INSERT - READ_LEN..]), &mut rng);
        out.push((format!("pair{}", out.len()), r1, r2));
    }
    out
}

fn mutate(seq: &[u8], rng: &mut Rng) -> Vec<u8> {
    const BASES: &[u8] = b"ACGT";
    seq.iter()
        .map(|&b| {
            if rng.unit() < SUB_RATE {
                // Pick a *different* base so the substitution always counts.
                let others: Vec<u8> = BASES.iter().copied().filter(|&x| x != b).collect();
                others[rng.below(others.len())]
            } else {
                b
            }
        })
        .collect()
}

fn write_fastq(path: &std::path::Path, reads: impl Iterator<Item = (String, Vec<u8>)>) {
    let mut f = std::fs::File::create(path).unwrap();
    for (name, seq) in reads {
        let qual = vec![b'I'; seq.len()];
        writeln!(f, "@{name}").unwrap();
        f.write_all(&seq).unwrap();
        f.write_all(b"\n+\n").unwrap();
        f.write_all(&qual).unwrap();
        f.write_all(b"\n").unwrap();
    }
}

// ---------------------------------------------------------------------------
// Record model
// ---------------------------------------------------------------------------

/// One alignment record reduced to comparable text, from either source.
#[derive(Debug, PartialEq, Eq)]
struct Record {
    qname: String,
    flag: u16,
    rname: String,
    pos: i32,
    mapq: u8,
    cigar: String,
    rnext: String,
    pnext: i32,
    tlen: i32,
    seq: String,
    qual: String,
    tags: BTreeMap<String, String>,
}

impl Record {
    /// Key identifying the same record on both sides. A read can appear more
    /// than once (supplementaries), and flag distinguishes those.
    fn key(&self) -> (String, u16) {
        (self.qname.clone(), self.flag)
    }

    /// Everything compared, rendered as one string so a failure prints the
    /// whole record rather than the first differing field.
    fn compare_key(&self) -> String {
        let vals: Vec<String> = COMPARED_TAGS
            .iter()
            .map(|t| format!("{t}:{}", self.tags.get(*t).map_or("-", String::as_str)))
            .collect();
        let keys: Vec<&str> = self.tags.keys().map(String::as_str).collect();
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tTAGS:{}",
            self.qname,
            self.flag,
            self.rname,
            self.pos,
            self.mapq,
            self.cigar,
            self.rnext,
            self.pnext,
            self.tlen,
            self.seq,
            self.qual,
            vals.join("\t"),
            keys.join(","),
        )
    }
}

/// Parse one SAM line (no header) as emitted by `bwa-mem3 mem`.
fn record_from_sam(line: &str) -> Record {
    let f: Vec<&str> = line.split('\t').collect();
    assert!(
        f.len() >= 11,
        "malformed SAM line ({} fields, need >= 11): {line}",
        f.len()
    );
    let mut tags = BTreeMap::new();
    for field in &f[11..] {
        let mut parts = field.splitn(3, ':');
        let (Some(key), Some(ty), Some(value)) = (parts.next(), parts.next(), parts.next()) else {
            panic!("malformed SAM aux field (want TAG:TYPE:VALUE): {field:?} in: {line}");
        };
        assert!(
            key.len() == 2 && ty.len() == 1,
            "malformed SAM aux field: {field:?} in: {line}"
        );
        // Store integer types under a single label: SAM text renders every
        // width as a bare number, while BAM picks the narrowest encoding, so
        // comparing the type code would be comparing storage, not content.
        let ty = if "cCsSiI".contains(ty) { "i" } else { ty };
        tags.insert(key.to_owned(), format!("{ty}:{value}"));
    }
    Record {
        qname: f[0].to_owned(),
        flag: f[1].parse().expect("FLAG"),
        rname: f[2].to_owned(),
        pos: f[3].parse().expect("POS"),
        mapq: f[4].parse().expect("MAPQ"),
        cigar: f[5].to_owned(),
        rnext: f[6].to_owned(),
        pnext: f[7].parse().expect("PNEXT"),
        tlen: f[8].parse().expect("TLEN"),
        seq: f[9].to_owned(),
        qual: f[10].to_owned(),
        tags,
    }
}

// ---------------------------------------------------------------------------
// Packed-BAM decoding (BAM spec 4.2)
// ---------------------------------------------------------------------------

const CIGAR_OPS: &[u8] = b"MIDNSHP=X";
const SEQ_NIBBLES: &[u8] = b"=ACMGRSVTWYHKDBN";

struct Cursor<'a> {
    b: &'a [u8],
    i: usize,
}

impl<'a> Cursor<'a> {
    fn u8(&mut self) -> u8 {
        let v = self.b[self.i];
        self.i += 1;
        v
    }
    fn u16(&mut self) -> u16 {
        let v = u16::from_le_bytes(self.b[self.i..self.i + 2].try_into().unwrap());
        self.i += 2;
        v
    }
    fn i32(&mut self) -> i32 {
        let v = i32::from_le_bytes(self.b[self.i..self.i + 4].try_into().unwrap());
        self.i += 4;
        v
    }
    fn take(&mut self, n: usize) -> &'a [u8] {
        let v = &self.b[self.i..self.i + n];
        self.i += n;
        v
    }
    fn cstr(&mut self) -> String {
        let start = self.i;
        while self.b[self.i] != 0 {
            self.i += 1;
        }
        let s = String::from_utf8_lossy(&self.b[start..self.i]).into_owned();
        self.i += 1; // NUL
        s
    }
    fn done(&self) -> bool {
        self.i >= self.b.len()
    }
}

/// Render one BAM aux field the way SAM text would.
fn decode_aux(c: &mut Cursor<'_>) -> (String, String) {
    let tag = String::from_utf8_lossy(c.take(2)).into_owned();
    let ty = c.u8();
    let value = match ty {
        b'A' => format!("A:{}", c.u8() as char),
        b'c' => format!("i:{}", c.u8() as i8),
        b'C' => format!("i:{}", c.u8()),
        b's' => format!("i:{}", c.u16() as i16),
        b'S' => format!("i:{}", c.u16()),
        b'i' => format!("i:{}", c.i32()),
        b'I' => format!("i:{}", c.i32() as u32),
        b'f' => format!("f:{}", f32::from_bits(c.i32() as u32)),
        b'Z' => format!("Z:{}", c.cstr()),
        b'H' => format!("H:{}", c.cstr()),
        b'B' => {
            let sub = c.u8();
            let n = c.i32() as usize;
            let mut parts = vec![(sub as char).to_string()];
            for _ in 0..n {
                parts.push(match sub {
                    b'c' => (c.u8() as i8).to_string(),
                    b'C' => c.u8().to_string(),
                    b's' => (c.u16() as i16).to_string(),
                    b'S' => c.u16().to_string(),
                    b'f' => f32::from_bits(c.i32() as u32).to_string(),
                    _ => c.i32().to_string(),
                });
            }
            format!("B:{}", parts.join(","))
        }
        other => panic!(
            "unknown BAM aux type code {:?} for tag {tag}",
            other as char
        ),
    };
    (tag, value)
}

/// Decode one packed record (`[u32 le block_size][block_size bytes]`) as
/// returned in `Record::bytes`.
fn record_from_packed_bam(bytes: &[u8], idx: &BwaIndex) -> Record {
    assert!(bytes.len() >= 4 + 32, "record too short: {}", bytes.len());
    let block_size = u32::from_le_bytes(bytes[0..4].try_into().unwrap()) as usize;
    assert_eq!(
        block_size + 4,
        bytes.len(),
        "block_size {block_size} does not match buffer length {}",
        bytes.len()
    );

    let mut c = Cursor {
        b: &bytes[4..],
        i: 0,
    };
    let ref_id = c.i32();
    let pos = c.i32();
    let l_read_name = c.u8() as usize;
    let mapq = c.u8();
    let _bin = c.u16();
    let n_cigar_op = c.u16() as usize;
    let flag = c.u16();
    let l_seq = c.i32() as usize;
    let next_ref_id = c.i32();
    let next_pos = c.i32();
    let tlen = c.i32();

    let qname = {
        let raw = c.take(l_read_name);
        String::from_utf8_lossy(&raw[..l_read_name - 1]).into_owned() // strip NUL
    };

    let mut cigar = String::new();
    for _ in 0..n_cigar_op {
        let op = c.i32() as u32;
        cigar.push_str(&(op >> 4).to_string());
        cigar.push(CIGAR_OPS[(op & 0xf) as usize] as char);
    }
    if cigar.is_empty() {
        cigar.push('*');
    }

    let packed_seq = c.take(l_seq.div_ceil(2)).to_vec();
    let mut seq = String::with_capacity(l_seq);
    for i in 0..l_seq {
        let byte = packed_seq[i / 2];
        let nib = if i % 2 == 0 { byte >> 4 } else { byte & 0xf };
        seq.push(SEQ_NIBBLES[nib as usize] as char);
    }
    if seq.is_empty() {
        seq.push('*');
    }

    let raw_qual = c.take(l_seq);
    let qual = if raw_qual.first() == Some(&0xff) {
        "*".to_owned()
    } else {
        raw_qual.iter().map(|&q| (q + 33) as char).collect()
    };

    let mut tags = BTreeMap::new();
    while !c.done() {
        let (k, v) = decode_aux(&mut c);
        tags.insert(k, v);
    }

    let name_of = |r: i32| -> String {
        if r < 0 {
            "*".to_owned()
        } else {
            idx.contig_name(r as usize).to_owned()
        }
    };
    let rnext = if next_ref_id >= 0 && next_ref_id == ref_id {
        "=".to_owned()
    } else {
        name_of(next_ref_id)
    };

    Record {
        qname,
        flag,
        rname: name_of(ref_id),
        pos: pos + 1, // BAM is 0-based, SAM is 1-based
        mapq,
        cigar,
        rnext,
        pnext: next_pos + 1,
        tlen,
        seq,
        qual,
        tags,
    }
}

// ---------------------------------------------------------------------------
// The test
// ---------------------------------------------------------------------------

#[test]
fn cli_and_align_batch_agree_on_every_record() {
    let Some((bwa, prefix)) = skip() else {
        return;
    };

    let idx = BwaIndex::load(&prefix).expect("load index");
    let Some((_region, ref_bases)) = pick_region(&idx, &prefix) else {
        eprintln!("skip: no gap-free window found in the reference");
        return;
    };

    let reads = simulate_pairs(&ref_bases);
    let tmp = tempfile::tempdir().unwrap();
    let r1_path = tmp.path().join("r1.fq");
    let r2_path = tmp.path().join("r2.fq");
    write_fastq(
        &r1_path,
        reads.iter().map(|(n, s, _)| (n.clone(), s.clone())),
    );
    write_fastq(
        &r2_path,
        reads.iter().map(|(n, _, s)| (n.clone(), s.clone())),
    );

    // --- CLI side: one chunk at -t 1. ---
    let out = Command::new(&bwa)
        .args([
            "mem",
            "-t",
            "1",
            &prefix,
            r1_path.to_str().unwrap(),
            r2_path.to_str().unwrap(),
        ])
        .output()
        .expect("run bwa-mem3 mem");
    assert!(
        out.status.success(),
        "bwa-mem3 mem failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let cli_text = String::from_utf8_lossy(&out.stdout);
    let cli: Vec<Record> = cli_text
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(record_from_sam)
        .collect();

    // --- Library side: one align_batch call, hence one insert-size cohort. ---
    let mut opts = MemOpts::new().expect("opts");
    opts.set_pe(true);
    let quals: Vec<Vec<u8>> = reads.iter().map(|(_, s, _)| vec![b'I'; s.len()]).collect();
    let pairs: Vec<ReadPair<'_>> = reads
        .iter()
        .zip(&quals)
        .map(|((name, s1, s2), q)| ReadPair {
            name_r1: name.as_bytes(),
            seq_r1: s1,
            qual_r1: Some(q),
            name_r2: name.as_bytes(),
            seq_r2: s2,
            qual_r2: Some(q),
        })
        .collect();
    let (aln, _pes) = align_batch(&idx, &opts, &pairs, None).expect("align_batch");
    let rs: Vec<Record> = aln
        .iter()
        .map(|r| record_from_packed_bam(r.bytes, &idx))
        .collect();

    // --- Compare. ---
    assert_eq!(
        rs.len(),
        cli.len(),
        "record count mismatch: bwa-rs={} CLI={}",
        rs.len(),
        cli.len()
    );

    // The fixture must actually exercise pairing and multi-mapping, or the
    // comparison below is vacuous. This is the check the old count-only test
    // lacked, and the reason its all-`ACGT` fixture went unnoticed.
    let paired = cli.iter().filter(|r| r.flag & 0x2 != 0).count();
    let with_xa = cli.iter().filter(|r| r.tags.contains_key("XA")).count();
    eprintln!(
        "fixture: {} records, {paired} properly paired, {with_xa} with XA:Z",
        cli.len()
    );
    assert!(
        paired > cli.len() / 2,
        "fixture does not exercise pairing: only {paired}/{} properly paired",
        cli.len()
    );
    assert!(
        with_xa > 0,
        "fixture has no multi-mapping reads, so XA:Z parity is untested"
    );

    let mut cli_by_key: BTreeMap<(String, u16), &Record> =
        cli.iter().map(|r| (r.key(), r)).collect();
    assert_eq!(cli_by_key.len(), cli.len(), "CLI emitted duplicate keys");

    let mut mismatches = Vec::new();
    for r in &rs {
        match cli_by_key.remove(&r.key()) {
            Some(c) if c.compare_key() == r.compare_key() => {}
            Some(c) => mismatches.push(format!(
                "record differs for {:?}:\n  CLI:    {}\n  bwa-rs: {}",
                r.key(),
                c.compare_key(),
                r.compare_key()
            )),
            None => mismatches.push(format!(
                "bwa-rs emitted a record the CLI did not: {:?}",
                r.key()
            )),
        }
    }
    for key in cli_by_key.keys() {
        mismatches.push(format!("CLI emitted a record bwa-rs did not: {key:?}"));
    }

    assert!(
        mismatches.is_empty(),
        "{} of {} records diverge:\n{}",
        mismatches.len(),
        cli.len(),
        mismatches
            .iter()
            .take(10)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}
