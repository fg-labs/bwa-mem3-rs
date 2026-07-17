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

/// Reduce a SAM line to the fields that must match the reference aligner:
/// qname, flag, rname, pos, mapq, cigar, and the `NM`/`MD`/`XG`/`XR`/`XM`/`XA`
/// tags. These are the semantically meaningful fields the shim reproduces
/// byte-for-byte; SEQ/QUAL/RNEXT/PNEXT/TLEN and the shim's own aux ordering
/// (which deliberately omits `MQ:i` and places the Bismark tags differently
/// from upstream's SAM writer) are excluded on purpose.
///
/// Callers compare the returned keys as *sorted vectors* — preserving record
/// multiplicity so a duplicated or missing record is caught, unlike a set.
/// Panics on a malformed SAM line rather than silently dropping it: `samtools
/// view` (no header) only ever emits ≥11-column records, so a short line is a
/// real bug, not noise to hide.
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
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\tNM:{nm}\tMD:{md}\tXG:{xg}\tXR:{xr}\tXM:{xm}\tXA:{xa}",
        f[0], f[1], f[2], f[3], f[4], f[5]
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
