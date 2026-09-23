//! Shared test support: reference indexes built with `bwa-mem3 index`, a
//! fixture covering every record shape the structured-fields sink reports, an
//! independent BAM record serializer, and helpers that drive a
//! [`ResidentCohort`] end to end. Included by each integration test binary
//! with `mod support;`; not every binary uses every item.
#![allow(dead_code)]

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Arc, OnceLock};

use bwa_mem3_rs::{
    AlignScratch, AlignedFields, AlignedFieldsSink, BwaIndex, Error, IdBases, Mate, MemOpts,
    MemPeStat, ReadPair, RecordOrigin, RecordSink, RecordVec, ResidentCohort, ResidentRange,
    SingleRead,
};

#[path = "../../../bwa-mem3-rs-cli/tests/phix_seq.rs"]
pub mod phix_seq;

/// A reference index built once per test-binary run, kept alive for the process.
pub struct RefIndex {
    _dir: tempfile::TempDir, // holds the on-disk index files alive
    pub prefix: PathBuf,
    pub idx: Arc<BwaIndex>,
    /// The reference's contigs, as written to the FASTA.
    pub contigs: Vec<(String, Vec<u8>)>,
}

/// Locate `bwa-mem3` via `BWA_MEM3_BIN` or `PATH`. Returns `None` (with a skip
/// message) when it is absent, unless `BWA_MEM3_RS_REQUIRE_TOOLS` is set, which
/// turns the absence into a hard failure.
pub fn find_bwa_mem3() -> Option<String> {
    if let Ok(p) = std::env::var("BWA_MEM3_BIN") {
        if Path::new(&p).exists() {
            return Some(p);
        }
    }
    let out = Command::new("which").arg("bwa-mem3").output().ok();
    let found = out.and_then(|o| {
        o.status
            .success()
            .then(|| String::from_utf8_lossy(&o.stdout).trim().to_string())
            .filter(|p| !p.is_empty())
    });
    if found.is_none() {
        assert!(
            std::env::var_os("BWA_MEM3_RS_REQUIRE_TOOLS").is_none(),
            "BWA_MEM3_RS_REQUIRE_TOOLS is set but bwa-mem3 was not found; \
             set BWA_MEM3_BIN or install it on PATH"
        );
        eprintln!("skip: bwa-mem3 not on PATH (set BWA_MEM3_BIN)");
    }
    found
}

/// Write `contigs` as `ref.fa` in a fresh temp dir, list `alt` contigs in
/// `ref.fa.alt` (which marks them as ALT loci when the index loads), and index
/// it with `bwa-mem3 index <extra_args>`. `None` (skip) without the tools.
pub fn build_reference(
    contigs: &[(&str, &[u8])],
    alt: &[&str],
    extra_args: &[&str],
) -> Option<(tempfile::TempDir, PathBuf)> {
    let bwa = find_bwa_mem3()?;
    let dir = tempfile::tempdir().expect("tempdir");
    let fa = dir.path().join("ref.fa");
    let mut f = std::fs::File::create(&fa).unwrap();
    for (name, seq) in contigs {
        writeln!(f, ">{name}").unwrap();
        for chunk in seq.chunks(72) {
            f.write_all(chunk).unwrap();
            writeln!(f).unwrap();
        }
    }
    drop(f);
    if !alt.is_empty() {
        let mut a = std::fs::File::create(dir.path().join("ref.fa.alt")).unwrap();
        for name in alt {
            writeln!(a, "{name}\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*").unwrap();
        }
    }
    let status = Command::new(&bwa)
        .arg("index")
        .args(extra_args)
        .arg(&fa)
        .status()
        .expect("run bwa-mem3 index");
    assert!(status.success(), "bwa-mem3 index failed");
    Some((dir, fa))
}

/// The PhiX genome as plain bases.
pub fn phix_bases() -> Vec<u8> {
    phix_seq::PHIX_SEQ
        .bytes()
        .filter(u8::is_ascii_alphabetic)
        .collect()
}

/// Build+load a PhiX index once, shared across the whole test binary. `None`
/// (skip) when `bwa-mem3` is unavailable and not required.
pub fn phix() -> Option<&'static RefIndex> {
    static REF: OnceLock<Option<RefIndex>> = OnceLock::new();
    REF.get_or_init(|| {
        let bases = phix_bases();
        let (dir, prefix) = build_reference(&[("phix", &bases)], &[], &[])?;
        let idx = Arc::new(BwaIndex::load(&prefix).expect("load PhiX index"));
        Some(RefIndex {
            _dir: dir,
            prefix,
            idx,
            contigs: vec![("phix".into(), bases)],
        })
    })
    .as_ref()
}

/// A PhiX `--meth` dual index (converted seed index + original reference),
/// built and loaded once per test binary. `None` (skip) without the tools.
pub fn phix_meth() -> Option<&'static RefIndex> {
    static REF: OnceLock<Option<RefIndex>> = OnceLock::new();
    REF.get_or_init(|| {
        let bases = phix_bases();
        let (dir, prefix) = build_reference(&[("phix", &bases)], &[], &["--meth"])?;
        let mut seed = prefix.clone().into_os_string();
        seed.push(".meth");
        let idx = Arc::new(BwaIndex::load_meth(&seed, &prefix).expect("load PhiX meth index"));
        Some(RefIndex {
            _dir: dir,
            prefix,
            idx,
            contigs: vec![("phix".into(), bases)],
        })
    })
    .as_ref()
}

/// A two-contig reference: `primary` (PhiX's first 3 kb) plus `alt`, a copy of
/// primary's bases 500..1500 with a few substitutions, listed in `ref.fa.alt`.
/// Reads from that stretch hit both, so they carry `XA:Z` and, their best hit
/// having an ALT competitor, `pa:f`. `None` (skip) without the tools.
pub fn alt_reference() -> Option<&'static RefIndex> {
    static REF: OnceLock<Option<RefIndex>> = OnceLock::new();
    REF.get_or_init(|| {
        let bases = phix_bases();
        let primary = bases[..3000].to_vec();
        let mut alt = primary[500..1500].to_vec();
        for i in (50..alt.len()).step_by(211) {
            alt[i] = if alt[i] == b'A' { b'C' } else { b'A' };
        }
        let (dir, prefix) =
            build_reference(&[("primary", &primary), ("alt", &alt)], &["alt"], &[])?;
        let idx = Arc::new(BwaIndex::load(&prefix).expect("load ALT index"));
        Some(RefIndex {
            _dir: dir,
            prefix,
            idx,
            contigs: vec![("primary".into(), primary), ("alt".into(), alt)],
        })
    })
    .as_ref()
}

/// The PhiX index prefix (its FASTA path), or `None` (skip) without the tools.
pub fn ref_prefix() -> Option<&'static Path> {
    phix().map(|p| p.prefix.as_path())
}

/// One `BwaIndex` shared by every test that only reads it, loaded at most
/// once per test binary run. See the module doc for why this matters.
pub fn shared_idx() -> Option<Arc<BwaIndex>> {
    phix().map(|p| p.idx.clone())
}

pub fn pairs(n: usize) -> (Vec<String>, Vec<Vec<u8>>, Vec<Vec<u8>>) {
    // Deterministic pseudo-reads: 150-mers over a fixed xorshift stream. Not
    // genomic, so most are unmapped -- enough to exercise ownership + parity.
    let mut x = 0x9E37_79B9u64;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    let mut names = Vec::new();
    let mut r1 = Vec::new();
    let mut r2 = Vec::new();
    for i in 0..n {
        names.push(format!("q{i}"));
        r1.push((0..150).map(|_| b"ACGT"[(next() % 4) as usize]).collect());
        r2.push((0..150).map(|_| b"ACGT"[(next() % 4) as usize]).collect());
    }
    (names, r1, r2)
}

// ---------------------------------------------------------------------------
// Structured-fields sink: `pair_emit_fields` vs packed `pair_emit`.
// ---------------------------------------------------------------------------

/// One input read, owned, as handed to the shim.
pub struct FixtureRead {
    pub name: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

/// Reads built from PhiX so they map, covering every record shape the
/// structured sink reports differently from a plain forward record: FR and RF pairs
/// (reverse strand on either mate), mismatches and indels (NM/MD), chimeric
/// reads (supplementaries, hard clips, SA:Z, clip-mode MC:Z), half-mapped
/// pairs (mate-coordinate copy onto the unmapped read), fully unmapped pairs,
/// missing QUAL (0xFF), and IUPAC/lowercase bases (2-bit N / case folding).
pub struct Fixture {
    pub pairs: Vec<(FixtureRead, FixtureRead)>,
    pub singles: Vec<FixtureRead>,
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// An owned fixture read with position-dependent QUAL (so a missed or doubled
/// reversal shows), or none.
pub fn fixture_read(name: String, seq: Vec<u8>, with_qual: bool) -> FixtureRead {
    let qual = with_qual.then(|| {
        (0..seq.len())
            .map(|i| b'!' + ((i * 7 + seq.len()) % 41) as u8)
            .collect()
    });
    FixtureRead {
        name: name.into_bytes(),
        seq,
        qual,
    }
}

/// [`field_fixture_from`] over PhiX.
pub fn field_fixture() -> Fixture {
    field_fixture_from(&phix_bases())
}

/// Reads drawn from `reference` (at least ~4 kb) covering every record shape
/// the structured sink restructures; see [`Fixture`].
pub fn field_fixture_from(reference: &[u8]) -> Fixture {
    let mut x = 0x2545_F491_4F6C_DD1Du64;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    let len = 150;
    let span = reference.len() - 1000;
    let mut pairs = Vec::new();
    let mut singles = Vec::new();
    for i in 0..160 {
        let p = (next() as usize) % span;
        let insert = 250 + (next() as usize) % 250;
        let fwd = reference[p..p + len].to_vec();
        let rev = revcomp(&reference[p + insert - len..p + insert]);
        let with_qual = i % 11 != 0;
        let (r1, r2) = match i % 8 {
            // FR / RF proper pairs.
            0 | 1 => (fwd, rev),
            2 => (rev, fwd),
            // Substitutions + a 2 bp deletion on R1, a 2 bp insertion on R2.
            3 => {
                let mut r1 = reference[p..p + 60].to_vec();
                r1.extend_from_slice(&reference[p + 62..p + len + 2]);
                r1[20] = if r1[20] == b'A' { b'C' } else { b'A' };
                let mut r2 = rev.clone();
                r2.insert(75, b'G');
                r2.insert(75, b'T');
                r2.truncate(len);
                (r1, r2)
            }
            // Chimeric R1: head from `p`, tail reverse-complemented from far
            // away, so the tail is a reverse-strand supplementary.
            4 => {
                let q = (p + 2000) % span;
                let mut r1 = reference[p..p + 90].to_vec();
                r1.extend(revcomp(&reference[q..q + 60]));
                (r1, rev)
            }
            // Half-mapped: R2 is random, R1 maps on either strand.
            5 => {
                let junk = (0..len).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
                (if i % 16 == 5 { fwd } else { rev }, junk)
            }
            // Both unmapped.
            6 => {
                let a = (0..len).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
                let b = (0..len).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
                (a, b)
            }
            // Lowercase and IUPAC bases on an otherwise proper pair.
            _ => {
                let mut r1 = fwd;
                r1[10] = b'N';
                r1[11] = b'R';
                r1[40..60].make_ascii_lowercase();
                (r1, rev)
            }
        };
        let name = format!("f{i}");
        pairs.push((
            fixture_read(name.clone(), r1, with_qual),
            fixture_read(name, r2, with_qual),
        ));
    }
    for i in 0..60 {
        let p = (next() as usize) % span;
        let seq = match i % 4 {
            0 => reference[p..p + len].to_vec(),
            1 => revcomp(&reference[p..p + len]),
            2 => {
                let q = (p + 2500) % span;
                let mut s = revcomp(&reference[p..p + 80]);
                s.extend_from_slice(&reference[q..q + 70]);
                s
            }
            _ => (0..len).map(|_| b"ACGT"[(next() % 4) as usize]).collect(),
        };
        singles.push(fixture_read(format!("s{i}"), seq, i % 7 != 0));
    }
    Fixture { pairs, singles }
}

/// One record reported to the structured sink, copied out of the callback.
#[derive(Debug)]
pub struct FieldRecord {
    pub origin: RecordOrigin,
    pub mate: Mate,
    pub is_primary: bool,
    /// The record rebuilt from its fields; compared to the packed body.
    pub body: Vec<u8>,
    /// Debug rendering of the fields, for failure messages.
    pub fields: String,
    pub flag: u16,
    pub tid: i32,
    pub query_len: usize,
    pub read_len: usize,
    pub has_sa: bool,
    pub has_mc: bool,
    pub has_xa: bool,
    pub has_pa: bool,
    pub rg: Option<Vec<u8>>,
    /// `(has XG:Z, has XM:Z)` when the record carries `--meth` tags.
    pub meth: Option<(bool, bool)>,
}

/// Smallest-width BAM integer aux value, the rule the packed path uses.
pub fn put_int_tag(out: &mut Vec<u8>, tag: &[u8; 2], v: i32) {
    out.extend_from_slice(tag);
    let v = i64::from(v);
    if (-128..=127).contains(&v) {
        out.push(b'c');
        out.push(v as i8 as u8);
    } else if (0..=255).contains(&v) {
        out.push(b'C');
        out.push(v as u8);
    } else if (-32768..=32767).contains(&v) {
        out.push(b's');
        out.extend_from_slice(&(v as i16).to_le_bytes());
    } else if (0..=65535).contains(&v) {
        out.push(b'S');
        out.extend_from_slice(&(v as u16).to_le_bytes());
    } else {
        out.push(b'i');
        out.extend_from_slice(&(v as i32).to_le_bytes());
    }
}

pub fn put_str_tag(out: &mut Vec<u8>, tag: &[u8; 2], v: &[u8]) {
    out.extend_from_slice(tag);
    out.push(b'Z');
    out.extend_from_slice(v);
    out.push(0);
}

/// BAM 4-bit code for an input base, via bwa's 2-bit alphabet (case-folded;
/// anything but ACGT is N), complemented on the reverse strand.
pub fn bam_base(ascii: u8, complement: bool) -> u8 {
    let two_bit = match ascii.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    };
    let two_bit = if complement && two_bit < 4 {
        3 - two_bit
    } else {
        two_bit
    };
    [1, 2, 4, 8, 15][two_bit]
}

/// BAM 4-bit code of an ASCII base under htslib's `seq_nt16_table`: IUPAC
/// ambiguity codes keep their own code, case-insensitive, anything else is N.
pub fn nt16_base(ascii: u8) -> u8 {
    match ascii.to_ascii_uppercase() {
        b'=' => 0,
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' | b'U' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        _ => 15,
    }
}

/// BAM 4-bit code of one emitted SEQ base from the input read, following the
/// upstream writers: without `--meth` through bwa's 2-bit alphabet; under
/// `--meth` the original base, IUPAC-preserving forward and 2-bit-complemented
/// reverse.
pub fn seq_code(ascii: u8, rev: bool, meth: bool) -> u8 {
    if meth && !rev {
        nt16_base(ascii)
    } else {
        bam_base(ascii, rev)
    }
}

/// Build a packed BAM record body from `f` and the read it belongs to. Pure
/// formatting, written independently of the shim, so agreement with the packed
/// body proves the fields carry everything that body does. `meth` selects the
/// `--meth` SEQ encoding.
pub fn serialize_fields(f: &AlignedFields<'_>, read: &FixtureRead, meth: bool) -> Vec<u8> {
    let query = f.query_start..f.query_end;
    let window: Vec<usize> = if f.is_rev() {
        query.rev().collect()
    } else {
        query.collect()
    };
    let mut out = Vec::new();
    out.extend_from_slice(&f.tid.to_le_bytes());
    out.extend_from_slice(&f.pos.to_le_bytes());
    out.push(u8::try_from(read.name.len() + 1).unwrap());
    out.push(f.mapq);
    out.extend_from_slice(&f.bin.to_le_bytes());
    out.extend_from_slice(&u16::try_from(f.cigar.len()).unwrap().to_le_bytes());
    out.extend_from_slice(&f.flag.to_le_bytes());
    out.extend_from_slice(&i32::try_from(window.len()).unwrap().to_le_bytes());
    out.extend_from_slice(&f.next_tid.to_le_bytes());
    out.extend_from_slice(&f.next_pos.to_le_bytes());
    out.extend_from_slice(&f.tlen.to_le_bytes());
    out.extend_from_slice(&read.name);
    out.push(0);
    for op in f.cigar {
        out.extend_from_slice(&op.to_le_bytes());
    }
    for pair in window.chunks(2) {
        let hi = seq_code(read.seq[pair[0]], f.is_rev(), meth);
        let lo = pair
            .get(1)
            .map_or(0, |&j| seq_code(read.seq[j], f.is_rev(), meth));
        out.push((hi << 4) | lo);
    }
    match &read.qual {
        Some(q) => out.extend(window.iter().map(|&j| q[j] - 33)),
        None => out.extend(window.iter().map(|_| 0xFF)),
    }
    if let Some(v) = f.nm {
        put_int_tag(&mut out, b"NM", v);
    }
    if let Some(v) = f.md {
        put_str_tag(&mut out, b"MD", v);
    }
    if let Some(v) = f.mc {
        put_str_tag(&mut out, b"MC", v);
    }
    if let Some(v) = f.mq {
        put_int_tag(&mut out, b"MQ", v);
    }
    if let Some(v) = f.score {
        put_int_tag(&mut out, b"AS", v);
    }
    if let Some(v) = f.sub {
        put_int_tag(&mut out, b"XS", v);
    }
    if let Some(v) = f.rg {
        put_str_tag(&mut out, b"RG", v);
    }
    if let Some(v) = f.sa {
        put_str_tag(&mut out, b"SA", v);
    }
    if let Some(v) = f.pa {
        out.extend_from_slice(b"paf");
        out.extend_from_slice(&v.to_le_bytes());
    }
    if let Some(v) = f.xa {
        put_str_tag(&mut out, b"XA", v);
    }
    if let Some(v) = f.hn {
        put_int_tag(&mut out, b"HN", v);
    }
    if let Some(m) = f.meth {
        put_str_tag(&mut out, b"XR", m.xr);
        if let Some(v) = m.xg {
            put_str_tag(&mut out, b"XG", v);
        }
        if let Some(v) = m.xm {
            put_str_tag(&mut out, b"XM", v);
        }
    }
    out
}

/// Structured sink that re-serializes every record against the fixture.
pub struct FieldCollector<'a> {
    pub fixture: &'a Fixture,
    /// Whether the run aligns under `--meth` (selects the SEQ encoding).
    pub meth: bool,
    pub records: Vec<FieldRecord>,
}

impl<'a> FieldCollector<'a> {
    pub fn new(fixture: &'a Fixture, meth: bool) -> Self {
        Self {
            fixture,
            meth,
            records: Vec::new(),
        }
    }
}

impl AlignedFieldsSink for FieldCollector<'_> {
    fn emit(&mut self, origin: RecordOrigin, mate: Mate, is_primary: bool, f: AlignedFields<'_>) {
        let read = match origin {
            RecordOrigin::Pair(i) => {
                let (r1, r2) = &self.fixture.pairs[i];
                match mate {
                    Mate::R1 => r1,
                    Mate::R2 => r2,
                }
            }
            RecordOrigin::Single(i) => &self.fixture.singles[i],
        };
        self.records.push(FieldRecord {
            origin,
            mate,
            is_primary,
            body: serialize_fields(&f, read, self.meth),
            fields: format!("{f:?}"),
            flag: f.flag,
            tid: f.tid,
            query_len: f.query_end - f.query_start,
            read_len: read.seq.len(),
            has_sa: f.sa.is_some(),
            has_mc: f.mc.is_some(),
            has_xa: f.xa.is_some(),
            has_pa: f.pa.is_some(),
            rg: f.rg.map(<[u8]>::to_vec),
            meth: f.meth.map(|m| (m.xg.is_some(), m.xm.is_some())),
        });
    }
}

/// Which emission path a cohort run drives.
pub enum EmitMode<'a> {
    Packed(&'a mut dyn RecordSink),
    Fields(&'a mut dyn AlignedFieldsSink),
}

/// Run the whole fixture through a fresh resident cohort in `sub`-pair
/// sub-chunks (singles as one range after the pairs), emitting to `mode`. A
/// fresh cohort per run matters: pair_emit is not idempotent on a range (mate
/// rescue appends regions to the mate's alnreg vector).
pub fn run_cohort(
    idx: &BwaIndex,
    opts: &MemOpts,
    fixture: &Fixture,
    sub: usize,
    mode: EmitMode<'_>,
) {
    let mut scratch = AlignScratch::new().unwrap();
    let cohort = ResidentCohort::new(opts.meth()).unwrap();
    let mut ranges: Vec<(ResidentRange, usize)> = Vec::new();
    for (k, chunk) in fixture.pairs.chunks(sub).enumerate() {
        let mut range = cohort.reserve_pairs(chunk.len()).unwrap();
        for (i, (r1, r2)) in chunk.iter().enumerate() {
            let pair = ReadPair {
                name_r1: &r1.name,
                seq_r1: &r1.seq,
                qual_r1: r1.qual.as_deref(),
                name_r2: &r2.name,
                seq_r2: &r2.seq,
                qual_r2: r2.qual.as_deref(),
            };
            cohort.write_pair(&mut range, i, pair).unwrap();
        }
        cohort
            .seed_extend(idx, opts, &mut scratch, &mut range)
            .unwrap();
        ranges.push((range, k * sub));
    }
    let mut singles = cohort.reserve_singles(fixture.singles.len()).unwrap();
    for (i, r) in fixture.singles.iter().enumerate() {
        let read = SingleRead {
            name: &r.name,
            seq: &r.seq,
            qual: r.qual.as_deref(),
        };
        cohort.write_single(&mut singles, i, read).unwrap();
    }
    cohort
        .seed_extend(idx, opts, &mut scratch, &mut singles)
        .unwrap();
    let pestat = cohort.infer_cohort(idx, opts).unwrap();

    let mut emit =
        |range: &mut ResidentRange, ids: IdBases, base: usize, mode: &mut EmitMode<'_>| {
            match mode {
                EmitMode::Packed(sink) => cohort.pair_emit(
                    idx,
                    opts,
                    &mut scratch,
                    range,
                    Some(&pestat),
                    ids,
                    base,
                    &mut **sink,
                ),
                EmitMode::Fields(sink) => cohort.pair_emit_fields(
                    idx,
                    opts,
                    &mut scratch,
                    range,
                    Some(&pestat),
                    ids,
                    base,
                    &mut **sink,
                ),
            }
            .unwrap();
        };
    let mut mode = mode;
    for (mut range, p0) in ranges {
        let ids = IdBases {
            first_single_id: 0,
            first_pair_id: p0 as u64,
        };
        emit(&mut range, ids, p0, &mut mode);
    }
    let ids = IdBases {
        first_single_id: (2 * fixture.pairs.len()) as u64,
        first_pair_id: 0,
    };
    emit(&mut singles, ids, 0, &mut mode);
}

/// Shapes the fixture must actually produce, so the parity check cannot pass
/// vacuously if a fixture or option change stops reaching a restructured path.
#[derive(Default, Debug)]
pub struct Coverage {
    pub reverse: usize,
    pub supplementary: usize,
    pub split_secondary: usize,
    pub hard_clipped_window: usize,
    pub sa: usize,
    pub mc: usize,
    pub xa: usize,
    pub pa: usize,
    pub rg: usize,
    pub meth_mapped: usize,
    pub placed_unmapped: usize,
    pub unmapped_unplaced: usize,
    pub singles: usize,
}

impl Coverage {
    /// Tally one record's shape.
    pub fn add(&mut self, got: &FieldRecord) {
        let mapped = got.flag & 0x4 == 0;
        self.reverse += usize::from(mapped && got.flag & 0x10 != 0);
        self.supplementary += usize::from(got.flag & 0x800 != 0);
        self.split_secondary += usize::from(got.flag & 0x100 != 0);
        self.hard_clipped_window += usize::from(got.query_len < got.read_len);
        self.sa += usize::from(got.has_sa);
        self.mc += usize::from(got.has_mc);
        self.xa += usize::from(got.has_xa);
        self.pa += usize::from(got.has_pa);
        self.rg += usize::from(got.rg.is_some());
        self.meth_mapped += usize::from(got.meth == Some((true, true)));
        self.placed_unmapped += usize::from(!mapped && got.tid >= 0);
        self.unmapped_unplaced += usize::from(!mapped && got.tid < 0);
        self.singles += usize::from(matches!(got.origin, RecordOrigin::Single(_)));
    }
}

/// Run `fixture` through both emission paths (packed and structured) of fresh
/// cohorts under `opts`, in `sub`-pair sub-chunks, and assert every
/// structured record, rebuilt, is byte-identical to its packed body, with the
/// same origin, and that `mate`/`is_primary` agree with its FLAG. Tallies each
/// record into `coverage`.
pub fn assert_fields_match_packed(
    idx: &BwaIndex,
    opts: &MemOpts,
    fixture: &Fixture,
    sub: usize,
    label: &str,
    coverage: &mut Coverage,
) {
    let mut packed = RecordVec::default();
    run_cohort(idx, opts, fixture, sub, EmitMode::Packed(&mut packed));
    let mut collector = FieldCollector::new(fixture, opts.meth());
    run_cohort(idx, opts, fixture, sub, EmitMode::Fields(&mut collector));

    assert_eq!(
        collector.records.len(),
        packed.records.len(),
        "{label} sub={sub}: record count differs"
    );
    for (n, (got, (origin, body))) in collector.records.iter().zip(&packed.records).enumerate() {
        assert_eq!(got.origin, *origin, "{label} sub={sub} record {n}: origin");
        if got.body != *body {
            let at = got
                .body
                .iter()
                .zip(body)
                .position(|(a, b)| a != b)
                .unwrap_or(got.body.len().min(body.len()));
            panic!(
                "{label} sub={sub} record {n} ({origin:?}, {:?}): rebuilt body differs from \
                 packed at byte {at} (lens {} vs {})\nfields: {}",
                got.mate,
                got.body.len(),
                body.len(),
                got.fields
            );
        }
        let expected_mate = match got.origin {
            RecordOrigin::Pair(_) if got.flag & 0x80 != 0 => Mate::R2,
            _ => Mate::R1,
        };
        assert_eq!(got.mate, expected_mate, "{label} record {n}: mate vs FLAG");
        assert_eq!(
            got.is_primary,
            got.flag & 0x900 == 0,
            "{label} record {n}: is_primary vs FLAG {:#x}",
            got.flag
        );
        coverage.add(got);
    }
}

// ---------------------------------------------------------------------------
// ResidentCohort concurrency + lifecycle.
// ---------------------------------------------------------------------------

/// Owned packed records in emission order.
pub type Records = Vec<(RecordOrigin, Vec<u8>)>;

pub fn as_read_pair(p: &(FixtureRead, FixtureRead)) -> ReadPair<'_> {
    let (r1, r2) = p;
    ReadPair {
        name_r1: &r1.name,
        seq_r1: &r1.seq,
        qual_r1: r1.qual.as_deref(),
        name_r2: &r2.name,
        seq_r2: &r2.seq,
        qual_r2: r2.qual.as_deref(),
    }
}

/// Reserve + write one pair range holding `pairs`.
pub fn written_pair_range(
    cohort: &ResidentCohort,
    pairs: &[(FixtureRead, FixtureRead)],
) -> ResidentRange {
    let mut range = cohort.reserve_pairs(pairs.len()).unwrap();
    for (i, p) in pairs.iter().enumerate() {
        cohort.write_pair(&mut range, i, as_read_pair(p)).unwrap();
    }
    range
}

/// Emit `range` (whose first pair is global pair `p0`) into owned packed records.
pub fn emit_range(
    cohort: &ResidentCohort,
    idx: &BwaIndex,
    opts: &MemOpts,
    scratch: &mut AlignScratch,
    range: &mut ResidentRange,
    pestat: &MemPeStat,
    p0: usize,
) -> Records {
    let ids = IdBases {
        first_single_id: 0,
        first_pair_id: p0 as u64,
    };
    let mut sink = RecordVec::default();
    cohort
        .pair_emit(idx, opts, scratch, range, Some(pestat), ids, p0, &mut sink)
        .unwrap();
    sink.records
}

/// Serial reference: every fixture pair in one range of one cohort. Returns the
/// cohort pestat and the records per global pair index. Emission is
/// independent of how pairs are grouped into ranges (the invariant the
/// sub-batch sweeps above prove), so any partition must reproduce these.
pub fn serial_reference(idx: &BwaIndex, opts: &MemOpts, fixture: &Fixture) -> (MemPeStat, Records) {
    let mut scratch = AlignScratch::new().unwrap();
    let cohort = ResidentCohort::new(false).unwrap();
    let mut range = written_pair_range(&cohort, &fixture.pairs);
    cohort
        .seed_extend(idx, opts, &mut scratch, &mut range)
        .unwrap();
    let pestat = cohort.infer_cohort(idx, opts).unwrap();
    let records = emit_range(&cohort, idx, opts, &mut scratch, &mut range, &pestat, 0);
    (pestat, records)
}

/// Records of `reference` whose origin is a pair in `[p0, p0 + n)`.
pub fn reference_slice(reference: &[(RecordOrigin, Vec<u8>)], p0: usize, n: usize) -> Records {
    reference
        .iter()
        .filter(|(o, _)| matches!(o, RecordOrigin::Pair(i) if (p0..p0 + n).contains(i)))
        .cloned()
        .collect()
}

/// The phrase every resident lifecycle violation's error carries.
pub const LIFECYCLE: &str = "lifecycle violated";

/// Assert `result` is the `InvalidInput` rejection whose message names
/// `needle`, so a test pins *which* misuse was caught, not merely that
/// something failed.
pub fn assert_invalid(result: bwa_mem3_rs::Result<()>, needle: &str, what: &str) {
    match result {
        Err(Error::InvalidInput(msg)) => assert!(
            msg.contains(needle),
            "{what}: expected an error naming {needle:?}, got {msg:?}"
        ),
        other => panic!("{what}: expected InvalidInput naming {needle:?}, got {other:?}"),
    }
}
