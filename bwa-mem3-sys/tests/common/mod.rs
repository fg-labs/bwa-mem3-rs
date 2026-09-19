//! Raw-FFI fixtures for the sys-crate behavior tests. Everything here is
//! `unsafe` by nature (this is the -sys crate); the safe wrapper's tests live
//! in `bwa-mem3-rs`. Test data is generated, never committed.

#![allow(dead_code)]

use std::ffi::{CStr, CString};
use std::io::Write;
use std::path::Path;
use std::process::Command;

use bwa_mem3_sys as sys;

#[path = "../../../bwa-mem3-rs-cli/tests/phix_seq.rs"]
mod phix_seq;
pub use phix_seq::PHIX_SEQ;

fn find_bwa_mem3() -> Option<String> {
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

/// Build a PhiX index in a temp dir. `None` (with a skip message) when
/// `bwa-mem3` is unavailable, unless `BWA_MEM3_RS_REQUIRE_TOOLS` is set.
pub fn phix_index() -> Option<(tempfile::TempDir, CString)> {
    let Some(bwa) = find_bwa_mem3() else {
        assert!(
            std::env::var_os("BWA_MEM3_RS_REQUIRE_TOOLS").is_none(),
            "BWA_MEM3_RS_REQUIRE_TOOLS is set but bwa-mem3 was not found"
        );
        eprintln!("skip: bwa-mem3 not on PATH (set BWA_MEM3_BIN)");
        return None;
    };
    let dir = tempfile::tempdir().expect("tempdir");
    let fa = dir.path().join("phix.fa");
    let mut f = std::fs::File::create(&fa).unwrap();
    writeln!(f, ">phix").unwrap();
    for chunk in PHIX_SEQ.as_bytes().chunks(72) {
        f.write_all(chunk).unwrap();
        writeln!(f).unwrap();
    }
    drop(f);
    let status = Command::new(bwa)
        .arg("index")
        .arg(&fa)
        .status()
        .expect("run bwa-mem3 index");
    assert!(status.success(), "bwa-mem3 index failed");
    let prefix = CString::new(fa.to_str().unwrap()).unwrap();
    Some((dir, prefix))
}

pub fn load_idx(prefix: &CStr) -> *mut sys::BwaIndex {
    let idx = unsafe { sys::bwa_shim_idx_load(prefix.as_ptr()) };
    assert!(!idx.is_null(), "index load failed");
    idx
}

pub fn new_opts() -> *mut sys::mem_opt_t {
    let o = unsafe { sys::bwa_shim_opts_new() };
    assert!(!o.is_null());
    o
}

pub fn new_pestat() -> *mut sys::mem_pestat_t {
    let p = unsafe { sys::bwa_shim_pestat_zero() };
    assert!(!p.is_null());
    p
}

/// Deterministic xorshift PRNG (same as the CLI-crate fixtures).
pub struct Rng(pub u64);
impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Owned read data; `pairs()` borrows it as FFI structs.
pub struct Fixture {
    pub names: Vec<CString>,
    pub r1: Vec<Vec<u8>>,
    pub r2: Vec<Vec<u8>>,
    pub qual: Vec<Vec<u8>>,
}

/// FR pairs drawn from PhiX with ~1% substitutions so MAPQ/NM/XA are non-trivial.
pub fn simulate(n: usize, read_len: usize, insert: usize, seed: u64) -> Fixture {
    let reference = PHIX_SEQ.as_bytes();
    let mut rng = Rng(seed);
    let max_start = reference.len() - insert;
    let mut names = Vec::with_capacity(n);
    let mut r1 = Vec::with_capacity(n);
    let mut r2 = Vec::with_capacity(n);
    let mut qual = Vec::with_capacity(n);
    for i in 0..n {
        let start = (rng.next() as usize) % (max_start + 1);
        let mut a = reference[start..start + read_len].to_vec();
        let mut b = revcomp(&reference[start + insert - read_len..start + insert]);
        for base in a.iter_mut().chain(b.iter_mut()) {
            if rng.next() % 100 == 0 {
                *base = b"ACGT"[(rng.next() % 4) as usize];
            }
        }
        names.push(CString::new(format!("r{i}")).unwrap());
        r1.push(a);
        r2.push(b);
        qual.push(vec![b'I'; read_len]);
    }
    Fixture {
        names,
        r1,
        r2,
        qual,
    }
}

impl Fixture {
    pub fn pairs(&self) -> Vec<sys::BwaReadPair> {
        (0..self.names.len())
            .map(|i| sys::BwaReadPair {
                r1_name: self.names[i].as_ptr(),
                r1_name_len: self.names[i].as_bytes().len(),
                r1_seq: self.r1[i].as_ptr(),
                r1_seq_len: self.r1[i].len(),
                r1_qual: self.qual[i].as_ptr(),
                r2_name: self.names[i].as_ptr(),
                r2_name_len: self.names[i].as_bytes().len(),
                r2_seq: self.r2[i].as_ptr(),
                r2_seq_len: self.r2[i].len(),
                r2_qual: self.qual[i].as_ptr(),
            })
            .collect()
    }
}

/// `(pair_idx, [u32 block_size][body])` per record, in emission order.
pub fn collect_records(b: *const sys::BwaBatch) -> Vec<(usize, Vec<u8>)> {
    let n = unsafe { sys::bwa_shim_batch_n_records(b) };
    (0..n)
        .map(|i| unsafe {
            let ptr = sys::bwa_shim_batch_record_ptr(b, i);
            let len = sys::bwa_shim_batch_record_len(b, i);
            (
                sys::bwa_shim_batch_pair_idx(b, i),
                std::slice::from_raw_parts(ptr, len).to_vec(),
            )
        })
        .collect()
}

/// Run the legacy `bwa_shim_align_batch` and collect its records.
pub fn align_batch_records(
    idx: *const sys::BwaIndex,
    opts: *const sys::mem_opt_t,
    pairs: &[sys::BwaReadPair],
) -> Vec<(usize, Vec<u8>)> {
    let pes = new_pestat();
    let b = unsafe {
        sys::bwa_shim_align_batch(
            idx,
            opts,
            pairs.as_ptr(),
            pairs.len(),
            std::ptr::null(),
            pes,
        )
    };
    assert!(!b.is_null(), "align_batch failed: {}", last_error());
    let recs = collect_records(b);
    unsafe {
        sys::bwa_shim_batch_free(b);
        sys::bwa_shim_pestat_free(pes);
    }
    recs
}

pub fn last_error() -> String {
    unsafe {
        let p = sys::bwa_shim_last_error();
        if p.is_null() {
            "(no error)".into()
        } else {
            CStr::from_ptr(p).to_string_lossy().into_owned()
        }
    }
}
