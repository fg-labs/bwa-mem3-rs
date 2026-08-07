//! Shared test fixtures: locate `bwa-mem3`, build a PhiX index, simulate
//! deterministic paired reads. Used by the e2e CLI test and the
//! flag-parity library test.

#![allow(dead_code)]

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

pub fn find_bwa_mem3() -> Option<String> {
    if let Ok(p) = std::env::var("BWA_MEM3_BIN") {
        if Path::new(&p).exists() {
            return Some(p);
        }
    }
    let out = Command::new("which").arg("bwa-mem3").output().ok()?;
    if !out.status.success() {
        return None;
    }
    let p = String::from_utf8(out.stdout).ok()?.trim().to_string();
    if p.is_empty() {
        None
    } else {
        Some(p)
    }
}

pub fn have_samtools() -> bool {
    Command::new("samtools").arg("--version").output().is_ok()
}

/// True when the caller has declared that `bwa-mem3` and `samtools` MUST be
/// present — set by the CI job that installs them. Turns a silent skip into a
/// hard failure so a broken CI environment cannot masquerade as a passing run.
///
/// Deliberately NOT `CI`: `cargo ci-test` runs these very targets in a job
/// that has no `bwa-mem3`, so gating on `CI` would fail the build job on every
/// platform. `CI` is also set by `act`, many container images, and direnv.
fn tools_required() -> bool {
    std::env::var_os("BWA_MEM3_RS_REQUIRE_TOOLS").is_some()
}

/// Locate `bwa-mem3`, panicking instead of returning `None` when
/// [`tools_required`] holds. Use this in place of [`find_bwa_mem3`] in tests
/// whose whole purpose is comparing against the reference aligner.
pub fn require_bwa_mem3() -> Option<String> {
    match find_bwa_mem3() {
        Some(p) => Some(p),
        None if tools_required() => panic!(
            "BWA_MEM3_RS_REQUIRE_TOOLS is set but bwa-mem3 was not found; \
             set BWA_MEM3_BIN or install it on PATH"
        ),
        None => {
            eprintln!("skip: bwa-mem3 not on PATH (set BWA_MEM3_BIN)");
            None
        }
    }
}

/// Check for `samtools`, panicking instead of returning `false` when
/// [`tools_required`] holds.
pub fn require_samtools() -> bool {
    if have_samtools() {
        return true;
    }
    assert!(
        !tools_required(),
        "BWA_MEM3_RS_REQUIRE_TOOLS is set but samtools was not found on PATH"
    );
    eprintln!("skip: samtools not on PATH");
    false
}

pub fn cli_bin() -> PathBuf {
    std::env::var("CARGO_BIN_EXE_bwa-rs").map_or_else(
        |_| PathBuf::from(env!("CARGO_BIN_EXE_bwa-rs")),
        PathBuf::from,
    )
}

/// Deterministic xorshift PRNG.
pub struct Rng(pub u64);
impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

pub fn simulate_pairs(
    ref_seq: &[u8],
    n: usize,
    read_len: usize,
    insert: usize,
    seed: u64,
) -> Vec<(String, Vec<u8>, Vec<u8>)> {
    let mut rng = Rng(seed);
    assert!(read_len <= insert, "read_len must be <= insert");
    let max_start = ref_seq
        .len()
        .checked_sub(insert)
        .expect("reference must be at least `insert` bases long");
    (0..n)
        .map(|i| {
            let start = if max_start == 0 {
                0
            } else {
                (rng.next() as usize) % (max_start + 1)
            };
            let r1 = ref_seq[start..start + read_len].to_vec();
            let span = &ref_seq[start + insert - read_len..start + insert];
            let r2 = revcomp(span);
            (format!("r{i}"), r1, r2)
        })
        .collect()
}

/// Deterministic ACGT string of length `n`, drawn from `rng`.
pub fn random_dna(rng: &mut Rng, n: usize) -> String {
    (0..n)
        .map(|_| ['A', 'C', 'G', 'T'][(rng.next() % 4) as usize])
        .collect()
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

/// `samtools view` a BAM file, returning its records as SAM text lines.
pub fn samtools_view(bam: &Path) -> Vec<String> {
    let out = Command::new("samtools")
        .args(["view"])
        .arg(bam)
        .output()
        .expect("run samtools view");
    assert!(out.status.success(), "samtools view failed");
    String::from_utf8_lossy(&out.stdout)
        .lines()
        .map(str::to_owned)
        .collect()
}

/// Aux tag keys that differ between the shim and upstream's SAM writer *by
/// design*, and so must not surface in [`record_key_fields`]' `TAGS:`
/// key-presence component. Add a key here, with a comment justifying it,
/// rather than special-casing it inline.
///
/// **Currently empty, and worth keeping that way.** Every entry is a hole in
/// the parity suite: a key listed here is subtracted from the comparison on
/// both sides, so a genuine divergence involving it cannot fail a test. The
/// list previously held `MQ`, which upstream emits on every record with a mate
/// and the shim did not; the shim now emits it too, so the exemption is gone
/// rather than merely justified.
const DELIBERATELY_ASYMMETRIC_TAG_KEYS: &[&str] = &[];

/// Extract the tag key from one aux field, requiring the SAM `TAG:TYPE:VALUE`
/// shape (a 2-character tag, a 1-character type code, then the value, which
/// may itself be empty for e.g. `XA:Z:`).
///
/// Validating rather than falling back to "the whole field is the key" is
/// load-bearing for the parity comparison, not defensive tidiness: a
/// truncated field reduces to a *plausible* key — and if that key is ever
/// listed in [`DELIBERATELY_ASYMMETRIC_TAG_KEYS`] it would then be dropped, so
/// a corrupt record would compare equal to a clean one and the parity check
/// would pass on it. (That list is empty today, which removes the specific
/// hazard but not the reason to validate.) Fail fast for the same reason the
/// ≥11-column assert in [`record_key_fields`] does.
fn aux_tag_key<'a>(field: &'a str, sam_line: &str) -> &'a str {
    let mut parts = field.splitn(3, ':');
    let (Some(key), Some(ty), Some(_value)) = (parts.next(), parts.next(), parts.next()) else {
        panic!("malformed SAM aux field (want TAG:TYPE:VALUE): {field:?} in: {sam_line}");
    };
    assert!(
        key.len() == 2 && ty.len() == 1,
        "malformed SAM aux field (want a 2-char tag and 1-char type): {field:?} in: {sam_line}"
    );
    key
}

/// Reduce a SAM line to the fields that must match the reference aligner:
/// qname, flag, rname, pos, mapq, cigar, and the `NM`/`MD`/`XG`/`XR`/`XM`/`XA`
/// tags. These are the semantically meaningful fields the shim reproduces
/// byte-for-byte; SEQ/QUAL/RNEXT/PNEXT/TLEN and the shim's own aux ordering
/// (which places the Bismark tags differently from upstream's SAM writer) are
/// excluded on purpose.
///
/// The key also appends a sorted `TAGS:<key,key,...>` component listing every
/// aux tag *key* present on the record, not just the six compared above,
/// except for [`DELIBERATELY_ASYMMETRIC_TAG_KEYS`] (currently empty).
/// Without the `TAGS:` component, a record carrying an extra tag this
/// reduction does not inspect (e.g. a newly-emitted upstream `XB`) would
/// reduce to the same key as a record without it, so a vendor bump that
/// starts emitting a new tag would pass every parity test unnoticed. The
/// exclusion list exists to keep that check from firing on asymmetries that
/// are known and intentional, and is currently empty — the shim emits every
/// tag upstream does. Only the set of keys is compared — values of uncompared tags stay out of scope, and
/// aux ordering stays unpinned, same as before.
///
/// Callers compare the returned keys as *sorted vectors* — preserving record
/// multiplicity so a duplicated or missing record is caught, unlike a set.
/// Panics on a malformed SAM line rather than silently dropping it: `samtools
/// view` (no header) only ever emits ≥11-column records, so a short line is a
/// real bug, not noise to hide. [`aux_tag_key`] applies the same rule to each
/// individual aux field.
pub fn record_key_fields(sam_line: &str) -> String {
    let f: Vec<&str> = sam_line.split('\t').collect();
    assert!(
        f.len() >= 11,
        "malformed SAM line ({} fields, need >= 11): {sam_line}",
        f.len()
    );
    let (mut nm, mut md, mut xg, mut xr, mut xm, mut xa) = ("", "", "", "", "", "");
    for tag in &f[11..] {
        if let Some(v) = tag.strip_prefix("NM:") {
            nm = v;
        } else if let Some(v) = tag.strip_prefix("MD:") {
            md = v;
        } else if let Some(v) = tag.strip_prefix("XG:") {
            xg = v;
        } else if let Some(v) = tag.strip_prefix("XR:") {
            xr = v;
        } else if let Some(v) = tag.strip_prefix("XM:") {
            xm = v;
        } else if let Some(v) = tag.strip_prefix("XA:") {
            xa = v;
        }
    }
    // Sorted set of aux tag KEYS present on the record, minus the keys the
    // shim and upstream disagree on by design (see
    // DELIBERATELY_ASYMMETRIC_TAG_KEYS). The loop above only compares six
    // tags by value; without this set, a record that carries an extra tag we
    // do not inspect (e.g. a newly-emitted upstream XB) reduces to the same
    // key as one that does not, and every parity test passes through the
    // change. Keys only — values of uncompared tags stay deliberately out of
    // scope, and ordering stays unpinned.
    let mut tag_keys: Vec<&str> = f[11..]
        .iter()
        .map(|t| aux_tag_key(t, sam_line))
        .filter(|k| !DELIBERATELY_ASYMMETRIC_TAG_KEYS.contains(k))
        .collect();
    tag_keys.sort_unstable();
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\tNM:{nm}\tMD:{md}\tXG:{xg}\tXR:{xr}\tXM:{xm}\tXA:{xa}\tTAGS:{}",
        f[0],
        f[1],
        f[2],
        f[3],
        f[4],
        f[5],
        tag_keys.join(",")
    )
}

pub fn write_fastq(path: &Path, reads: &[(String, Vec<u8>)]) {
    let mut f = fs::File::create(path).unwrap();
    for (name, seq) in reads {
        let qual = vec![b'I'; seq.len()];
        writeln!(f, "@{name}").unwrap();
        f.write_all(seq).unwrap();
        writeln!(f, "\n+").unwrap();
        f.write_all(&qual).unwrap();
        writeln!(f).unwrap();
    }
}

/// Write the embedded PhiX sequence as a FASTA and build a bwa-mem3 index.
/// Returns the FASTA path (suitable as an index prefix for bwa-mem3-rs).
pub fn setup_phix_index(dir: &Path, bwa_mem3_bin: &str, phix_seq: &str) -> PathBuf {
    setup_phix_index_inner(dir, bwa_mem3_bin, phix_seq, &[], "fa.bwt.2bit.64")
}

/// Like [`setup_phix_index`] but builds the bisulfite dual index
/// (`bwa-mem3 index --meth`), producing `<ref>.*` + `<ref>.meth.*`.
pub fn setup_phix_meth_index(dir: &Path, bwa_mem3_bin: &str, phix_seq: &str) -> PathBuf {
    setup_phix_index_inner(
        dir,
        bwa_mem3_bin,
        phix_seq,
        &["--meth"],
        "fa.meth.bwt.2bit.64",
    )
}

fn setup_phix_index_inner(
    dir: &Path,
    bwa_mem3_bin: &str,
    phix_seq: &str,
    extra_args: &[&str],
    expect_suffix: &str,
) -> PathBuf {
    let ref_fa = dir.join("phix.fa");
    let mut f = fs::File::create(&ref_fa).unwrap();
    writeln!(f, ">phix").unwrap();
    for chunk in phix_seq.as_bytes().chunks(72) {
        f.write_all(chunk).unwrap();
        writeln!(f).unwrap();
    }
    drop(f);

    let status = Command::new(bwa_mem3_bin)
        .arg("index")
        .args(extra_args)
        .arg(&ref_fa)
        .status()
        .expect("run bwa-mem3 index");
    assert!(status.success(), "bwa-mem3 index failed");
    assert!(
        ref_fa.with_extension(expect_suffix).exists(),
        "bwa-mem3 index did not produce {expect_suffix}"
    );
    ref_fa
}
