//! Byte-level CLI parity for the three-phase API's single-end and mixed
//! (`-p` smart-pairing) paths, checked against the real `bwa-mem3` binary.
//!
//! The shim emits packed BAM record bodies; `bwa-mem3 mem` writes SAM text. We
//! reduce both to the same key (qname, flag, rname, pos, mapq, cigar, and every
//! aux tag as `KEY:type:value` with integer subtypes normalized to `i` exactly
//! as SAM does) and compare them as sorted multisets — order-independent,
//! because the CLI reorders records back to input order above the shim while
//! the shim emits pairs-group then singles-group.
//!
//! SEQ/QUAL/RNEXT/PNEXT/TLEN are excluded from the key, matching the
//! established `cli_parity_*` reduction: those are validated for pairs by the
//! all-pairs suite, and singles carry none of them. What this pins that the
//! FFI-level tests cannot: the SE emit *policy* (flags, cigar, NM/MD/XA/AS/XS)
//! is byte-identical to `bwa-mem3`'s single-end group, and a MIXED cohort's
//! pair records (incl. their pestat-dependent MAPQ/flags) still match the CLI's
//! `-p` output when singles share the cohort.
//!
//! Requires `bwa-mem3` on PATH (or `BWA_MEM3_BIN`); skips gracefully otherwise.

mod common;

use bwa_mem3_sys as sys;
use std::io::Write;
use std::path::Path;
use std::process::Command;

/// Locate the `bwa-mem3` binary the same way `common::phix_index` does.
fn bwa_bin() -> Option<String> {
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
    (!p.is_empty()).then_some(p)
}

struct Sink(Vec<(u32, usize, Vec<u8>)>);

/// # Safety
/// `ctx` must be a valid `*mut Sink` that outlives the call, and `body`/`len`
/// must describe a readable buffer of `len` bytes -- both guaranteed by the
/// shim when it invokes this callback from `bwa_shim_pair_emit`.
unsafe extern "C" fn sink_fn(
    ctx: *mut std::ffi::c_void,
    kind: u32,
    idx: usize,
    body: *const u8,
    len: usize,
) {
    // SAFETY: `ctx` is the `&mut Sink` the caller passed as the callback context.
    let sink = &mut *ctx.cast::<Sink>();
    // SAFETY: `body`/`len` describe a valid record buffer for this call (see the
    // `# Safety` section); the bytes are copied out immediately.
    sink.0
        .push((kind, idx, std::slice::from_raw_parts(body, len).to_vec()));
}

/// Contig names of a loaded index, indexed by refID.
fn contig_names(idx: *const sys::BwaIndex) -> Vec<String> {
    // SAFETY: `idx` is a valid non-null index handle alive for the call; the
    // accessor only reads it.
    let n = unsafe { sys::bwa_shim_idx_n_contigs(idx) };
    (0..n)
        // SAFETY: `i < n`, so contig `i` exists; the shim returns a pointer to a
        // NUL-terminated name owned by `idx` and valid while `idx` lives, copied
        // into an owned `String` here.
        .map(|i| unsafe {
            std::ffi::CStr::from_ptr(sys::bwa_shim_idx_contig_name(idx, i))
                .to_string_lossy()
                .into_owned()
        })
        .collect()
}

/// Reduce one aux block (BAM binary) to a sorted list of `KEY:type:value`
/// strings, normalizing integer subtypes (c/C/s/S/i/I) to `i` as SAM does.
fn decode_aux(mut a: &[u8]) -> Vec<String> {
    let mut tags = Vec::new();
    while a.len() >= 3 {
        let key = format!("{}{}", a[0] as char, a[1] as char);
        let ty = a[2] as char;
        a = &a[3..];
        let (rendered, consumed): (String, usize) = match ty {
            'A' => (format!("A:{}", a[0] as char), 1),
            'c' => ((a[0] as i8 as i64).to_string_i(), 1),
            'C' => ((a[0] as i64).to_string_i(), 1),
            's' => ((i16::from_le_bytes([a[0], a[1]]) as i64).to_string_i(), 2),
            'S' => ((u16::from_le_bytes([a[0], a[1]]) as i64).to_string_i(), 2),
            'i' => (
                (i32::from_le_bytes([a[0], a[1], a[2], a[3]]) as i64).to_string_i(),
                4,
            ),
            'I' => (
                (u32::from_le_bytes([a[0], a[1], a[2], a[3]]) as i64).to_string_i(),
                4,
            ),
            'f' => (
                format!("f:{}", f32::from_le_bytes([a[0], a[1], a[2], a[3]])),
                4,
            ),
            'Z' | 'H' => {
                let end = a
                    .iter()
                    .position(|&b| b == 0)
                    .expect("unterminated Z/H aux");
                (
                    format!("{ty}:{}", String::from_utf8_lossy(&a[..end])),
                    end + 1,
                )
            }
            other => panic!("unhandled BAM aux type {other:?}"),
        };
        tags.push(format!("{key}:{rendered}"));
        a = &a[consumed..];
    }
    assert!(a.is_empty(), "trailing bytes in aux block: {a:?}");
    tags.sort();
    tags
}

/// Small helper so integer branches read as `<value>.to_string_i()` → `i:<v>`.
trait ToSamInt {
    fn to_string_i(self) -> String;
}
impl ToSamInt for i64 {
    fn to_string_i(self) -> String {
        format!("i:{self}")
    }
}

/// Reduce one shim BAM record body to the comparison key.
fn bam_key(body: &[u8], contigs: &[String]) -> String {
    let refid = i32::from_le_bytes(body[0..4].try_into().unwrap());
    let pos = i32::from_le_bytes(body[4..8].try_into().unwrap());
    let l_read_name = body[8] as usize;
    let mapq = body[9];
    let n_cigar = u16::from_le_bytes([body[12], body[13]]) as usize;
    let flag = u16::from_le_bytes([body[14], body[15]]);
    let l_seq = i32::from_le_bytes(body[16..20].try_into().unwrap()) as usize;

    let name_off = 32;
    let qname = String::from_utf8_lossy(&body[name_off..name_off + l_read_name - 1]).into_owned();
    let cigar_off = name_off + l_read_name;
    let mut cigar = String::new();
    for i in 0..n_cigar {
        let c = u32::from_le_bytes(
            body[cigar_off + 4 * i..cigar_off + 4 * i + 4]
                .try_into()
                .unwrap(),
        );
        cigar.push_str(&(c >> 4).to_string());
        cigar.push("MIDNSHP=X".as_bytes()[(c & 0xf) as usize] as char);
    }
    if cigar.is_empty() {
        cigar.push('*');
    }
    let seq_off = cigar_off + 4 * n_cigar;
    let qual_off = seq_off + l_seq.div_ceil(2);
    let aux_off = qual_off + l_seq;

    let rname = if refid < 0 {
        "*".to_string()
    } else {
        contigs[refid as usize].clone()
    };
    let pos1 = if pos < 0 { 0 } else { pos + 1 };
    let tags = decode_aux(&body[aux_off..]).join("\t");
    format!("{qname}\t{flag}\t{rname}\t{pos1}\t{mapq}\t{cigar}\t{tags}")
}

/// Reduce one SAM line (from `bwa-mem3 mem`) to the same key.
fn sam_key(line: &str) -> String {
    let f: Vec<&str> = line.split('\t').collect();
    assert!(f.len() >= 11, "short SAM line: {line}");
    let mut tags: Vec<String> = f[11..].iter().map(|s| s.to_string()).collect();
    tags.sort();
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        f[0],
        f[1],
        f[2],
        f[3],
        f[4],
        f[5],
        tags.join("\t")
    )
}

/// SAM record lines (header stripped) from a `bwa-mem3 mem` invocation.
fn cli_sam_keys(args: &[&std::ffi::OsStr]) -> Vec<String> {
    let out = Command::new(args[0])
        .args(&args[1..])
        .output()
        .expect("run bwa-mem3 mem");
    assert!(
        out.status.success(),
        "bwa-mem3 mem failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(sam_key)
        .collect()
}

fn write_fastq(path: &Path, reads: &[(String, Vec<u8>)]) {
    let mut f = std::fs::File::create(path).unwrap();
    for (name, seq) in reads {
        writeln!(f, "@{name}").unwrap();
        f.write_all(seq).unwrap();
        writeln!(f, "\n+").unwrap();
        f.write_all(&vec![b'I'; seq.len()]).unwrap();
        writeln!(f).unwrap();
    }
}

/// Single-end input: `bwa-mem3 mem ref reads.fq` vs the shim's singles group.
#[test]
fn single_end_matches_bwa_mem3_cli() {
    let Some(bwa) = bwa_bin() else {
        eprintln!("skip: bwa-mem3 not found");
        return;
    };
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let contigs = contig_names(idx);

    // Unique single-end reads drawn straight from PhiX (repeat-free → unique
    // mappers, so the tie-break id is irrelevant to the bytes).
    let reference = common::PHIX_SEQ.as_bytes();
    let mut rng = common::Rng(101);
    let read_len = 100usize;
    let mut reads: Vec<(String, Vec<u8>)> = Vec::new();
    for i in 0..60 {
        let start = (rng.next() as usize) % (reference.len() - read_len);
        let mut seq = reference[start..start + read_len].to_vec();
        if rng.next() % 2 == 0 {
            seq = common::revcomp(&seq);
        }
        reads.push((format!("s{i}"), seq));
    }

    // Reference: bwa-mem3 single-end (one fastq file, no -p).
    let dir = tempfile::tempdir().unwrap();
    let fq = dir.path().join("reads.fq");
    write_fastq(&fq, &reads);
    let prefix_os = prefix.to_str().unwrap();
    let mut cli = cli_sam_keys(&[
        bwa.as_ref(),
        "mem".as_ref(),
        "-t".as_ref(),
        "1".as_ref(),
        prefix_os.as_ref(),
        fq.as_os_str(),
    ]);

    // Shim: same reads as a singles-only batch.
    let names: Vec<std::ffi::CString> = reads
        .iter()
        .map(|(n, _)| std::ffi::CString::new(n.as_str()).unwrap())
        .collect();
    let quals: Vec<Vec<u8>> = reads.iter().map(|(_, s)| vec![b'I'; s.len()]).collect();
    let singles: Vec<sys::BwaSingleRead> = (0..reads.len())
        .map(|i| sys::BwaSingleRead {
            name: names[i].as_ptr(),
            name_len: names[i].as_bytes().len(),
            seq: reads[i].1.as_ptr(),
            seq_len: reads[i].1.len(),
            qual: quals[i].as_ptr(),
        })
        .collect();
    let mut shim = run_shim_keys(idx, opts, &[], &singles, 0, 0, &contigs);

    assert!(!cli.is_empty(), "reference produced no records");
    assert_eq!(shim.len(), cli.len(), "record count differs");
    cli.sort();
    shim.sort();
    assert_eq!(shim, cli, "single-end records diverge from bwa-mem3 CLI");

    // SAFETY: `opts`/`idx` are the live owned handles; each freed once and not
    // used afterward.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// Mixed cohort: an interleaved file of consecutive-same-name pairs and unique
/// singletons through `bwa-mem3 mem -p` vs the shim's mixed batch. Pins that
/// the pairs-only cohort pestat matches the CLI's `-p` insert-size model even
/// with singles sharing the cohort.
#[test]
fn mixed_matches_bwa_mem3_cli_smart_pairing() {
    let Some(bwa) = bwa_bin() else {
        eprintln!("skip: bwa-mem3 not found");
        return;
    };
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let contigs = contig_names(idx);

    // FR pairs from PhiX (unique) + unique singletons.
    let fx = common::simulate(20, 100, 300, 55);
    let reference = common::PHIX_SEQ.as_bytes();
    let mut rng = common::Rng(202);
    let read_len = 100usize;
    let n_singles = 15usize;
    let mut single_seqs: Vec<Vec<u8>> = Vec::new();
    for _ in 0..n_singles {
        let start = (rng.next() as usize) % (reference.len() - read_len);
        let mut seq = reference[start..start + read_len].to_vec();
        if rng.next() % 2 == 0 {
            seq = common::revcomp(&seq);
        }
        single_seqs.push(seq);
    }

    // Interleaved reference input: [p_i R1, p_i R2, g_i single] repeated. A
    // single is placed after each pair so bseq_classify (consecutive-same-name)
    // classifies the pair as PE and the singleton (unique name) as SE — the
    // -p SE/PE split the shim mirrors by construction (pairs[] + singles[]).
    let mut inter: Vec<(String, Vec<u8>)> = Vec::new();
    for (i, name_c) in fx.names.iter().enumerate() {
        let name = name_c.to_str().unwrap().to_string();
        inter.push((name.clone(), fx.r1[i].clone()));
        inter.push((name, fx.r2[i].clone()));
        if let Some(s) = single_seqs.get(i) {
            inter.push((format!("g{i}"), s.clone()));
        }
    }
    for (j, s) in single_seqs.iter().enumerate().skip(fx.names.len()) {
        inter.push((format!("g{j}"), s.clone()));
    }

    let dir = tempfile::tempdir().unwrap();
    let fq = dir.path().join("inter.fq");
    write_fastq(&fq, &inter);
    let prefix_os = prefix.to_str().unwrap();
    let mut cli = cli_sam_keys(&[
        bwa.as_ref(),
        "mem".as_ref(),
        "-t".as_ref(),
        "1".as_ref(),
        "-p".as_ref(),
        prefix_os.as_ref(),
        fq.as_os_str(),
    ]);

    // Shim: the same pairs and singletons, split into the two groups. Ids follow
    // the CLI's -p layout (SE ids from n_processed=0; PE ids from n_sep[0]>>1);
    // irrelevant to the bytes here (unique mappers) but faithful to the layout.
    let pairs = fx.pairs();
    let names: Vec<std::ffi::CString> = (0..n_singles)
        .map(|j| std::ffi::CString::new(format!("g{j}")).unwrap())
        .collect();
    let quals: Vec<Vec<u8>> = single_seqs.iter().map(|s| vec![b'I'; s.len()]).collect();
    let singles: Vec<sys::BwaSingleRead> = (0..n_singles)
        .map(|j| sys::BwaSingleRead {
            name: names[j].as_ptr(),
            name_len: names[j].as_bytes().len(),
            seq: single_seqs[j].as_ptr(),
            seq_len: single_seqs[j].len(),
            qual: quals[j].as_ptr(),
        })
        .collect();
    let first_pair_id = (n_singles as u64) >> 1;
    let mut shim = run_shim_keys(idx, opts, &pairs, &singles, 0, first_pair_id, &contigs);

    assert!(!cli.is_empty(), "reference produced no records");
    assert_eq!(shim.len(), cli.len(), "record count differs");
    cli.sort();
    shim.sort();
    assert_eq!(shim, cli, "mixed -p records diverge from bwa-mem3 CLI");

    // SAFETY: `opts`/`idx` are the live owned handles; each freed once and not
    // used afterward.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// Drive the three phases over a pair+single batch and return per-record keys.
fn run_shim_keys(
    idx: *const sys::BwaIndex,
    opts: *const sys::mem_opt_t,
    pairs: &[sys::BwaReadPair],
    singles: &[sys::BwaSingleRead],
    first_single_id: u64,
    first_pair_id: u64,
    contigs: &[String],
) -> Vec<String> {
    // SAFETY: constructor with no arguments; returns an owned scratch handle.
    let sc = unsafe { sys::bwa_shim_scratch_new() };
    let batch = sys::BwaReadBatch {
        pairs: pairs.as_ptr(),
        n_pairs: pairs.len(),
        singles: singles.as_ptr(),
        n_singles: singles.len(),
    };
    // SAFETY: valid `idx`/`opts`/`sc`; `batch` names live `pairs`/`singles`
    // slices with matching `n_pairs`/`n_singles` counts, alive for the call.
    let regs = unsafe { sys::bwa_shim_seed_extend(idx, opts, sc, &batch) };
    assert!(!regs.is_null(), "{}", common::last_error());
    let pes = common::new_pestat();
    let regs_c: [*const sys::BwaRegs; 1] = [regs.cast_const()];
    // SAFETY: valid `idx`/`opts`/`pes`; `regs_c` is a 1-element array of the live
    // `regs` pointer and the count `1` matches it.
    let rc = unsafe { sys::bwa_shim_pestat_cohort(idx, opts, regs_c.as_ptr(), 1, pes) };
    assert_eq!(rc, 0);
    let mut sink = Sink(Vec::new());
    let ids = sys::BwaIdBases {
        first_single_id,
        first_pair_id,
    };
    let pes_arg = if pairs.is_empty() {
        std::ptr::null()
    } else {
        pes.cast_const()
    };
    // SAFETY: valid `idx`/`opts`/`sc`/`regs`; `pes_arg` is the live pestat for a
    // paired batch or null for a pair-free one (the shim's contract); `sink_fn`'s
    // context is `&mut sink`, live for the whole call (see `sink_fn`'s `# Safety`).
    let rc = unsafe {
        sys::bwa_shim_pair_emit(
            idx,
            opts,
            sc,
            regs,
            pes_arg,
            ids,
            Some(sink_fn),
            (&mut sink as *mut Sink).cast(),
        )
    };
    assert_eq!(rc, 0, "{}", common::last_error());
    // SAFETY: `pes`/`sc` are the live owned handles; each freed once, not used
    // afterward.
    unsafe {
        sys::bwa_shim_pestat_free(pes);
        sys::bwa_shim_scratch_free(sc);
    }
    sink.0
        .iter()
        .map(|(_, _, body)| bam_key(body, contigs))
        .collect()
}

/// Regression guard for the batched mate-rescue chunking (finding [1]):
/// `shim_pair_emit`'s `BWAMEM_BATCHED_MATESW` path now processes the batch's
/// pairs in `BATCH_SIZE`-sized chunks, resetting the kswv `pcnt`/`gcnt`/ref-window
/// offset per chunk exactly like the CLI's `kt_for`-dispatched `worker_sam` — so
/// the int32 `SeqPair.idr` seqBuf offset can never run past 2^31 and trip
/// `seqbuf_capacity_fatal()`'s `exit()` across the FFI boundary.
///
/// This drives a SINGLE `pair_emit` call over 1300 pairs (2600 reads). With
/// `BATCH_SIZE` = 512 (x86) that is 6 chunks; with 1024 (arm64) it is 3 chunks —
/// so the internal chunk loop crosses at least two `BATCH_SIZE` boundaries on
/// every target. Mate rescue is FORCED to fire (every 3rd R2 mutated past
/// seedability), which is what enqueues the rescue jobs whose running offset the
/// old code let grow unbounded; a chunk that failed to reset would diverge from
/// the CLI's per-work-item reset. The CLI is the independent oracle here: it
/// chunks in its own `kt_for` code, untouched by this shim change, so a
/// chunk-boundary bug in the shim shows up as a mismatch against it (an
/// `align_batch`-vs-`three_phase` self-comparison would share the buggy code and
/// pass vacuously). `-K` is large enough to keep the CLI in one `-K` cohort, so
/// its pestat matches the shim's single-cohort estimate.
///
/// This is also the batched-vs-scalar mate-rescue divergence guard (finding
/// [19]): CI runs this exact target twice -- once in the default BATCHED build
/// and once with `CXXFLAGS=-DDISABLE_BATCHED_MATESW=1` (the scalar A/B step in
/// `.github/workflows/check.yml`). The batched run asserts `batched == CLI` and
/// the scalar run asserts `scalar == CLI` on this same rescue-firing fixture, so
/// transitively `batched == scalar` where rescue actually fires -- no separate
/// batched-vs-scalar test is needed.
#[test]
fn batched_rescue_spanning_multiple_chunks_matches_bwa_mem3_cli() {
    let Some(bwa) = bwa_bin() else {
        eprintln!("skip: bwa-mem3 not found");
        return;
    };
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let contigs = contig_names(idx);

    // 1300 pairs => 2600 reads: >= 3 chunks at BATCH_SIZE=1024, 6 at 512, so the
    // internal chunk loop crosses >= 2 BATCH_SIZE boundaries on every target.
    let mut fx = common::simulate(1300, 100, 300, 71);
    // Same rescue-forcing recipe as the FFI rescue-heavy test: mutate every 3rd
    // R2 at every 7th base (~14 substitutions in 100 bp) so no 19-mer seed
    // survives and R2 is placeable only by mate rescue from R1. This is what
    // enqueues the batched rescue jobs across the chunk boundaries.
    for i in (0..fx.r2.len()).step_by(3) {
        for j in (0..fx.r2[i].len()).step_by(7) {
            fx.r2[i][j] = match fx.r2[i][j] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
    }

    // CLI reference: interleaved (consecutive-same-name) pairs through `-p`, one
    // -K cohort so the insert-size model matches the shim's single cohort.
    let mut inter: Vec<(String, Vec<u8>)> = Vec::with_capacity(2 * fx.names.len());
    for (i, name_c) in fx.names.iter().enumerate() {
        let name = name_c.to_str().unwrap().to_string();
        inter.push((name.clone(), fx.r1[i].clone()));
        inter.push((name, fx.r2[i].clone()));
    }
    let dir = tempfile::tempdir().unwrap();
    let fq = dir.path().join("inter.fq");
    write_fastq(&fq, &inter);
    let prefix_os = prefix.to_str().unwrap();
    let mut cli = cli_sam_keys(&[
        bwa.as_ref(),
        "mem".as_ref(),
        "-t".as_ref(),
        "1".as_ref(),
        "-K".as_ref(),
        "100000000".as_ref(),
        "-p".as_ref(),
        prefix_os.as_ref(),
        fq.as_os_str(),
    ]);

    // Shim: the same pairs in ONE seed_extend + ONE pair_emit call, so the
    // batched rescue path chunks internally across the BATCH_SIZE boundaries.
    let pairs = fx.pairs();
    let mut shim = run_shim_keys(idx, opts, &pairs, &[], 0, 0, &contigs);

    assert!(!cli.is_empty(), "reference produced no records");
    assert_eq!(
        shim.len(),
        2 * pairs.len(),
        "shim did not emit one record per read"
    );
    assert_eq!(shim.len(), cli.len(), "record count differs");
    cli.sort();
    shim.sort();
    assert_eq!(
        shim, cli,
        "multi-chunk batched rescue diverged from bwa-mem3 CLI"
    );

    // Rescue must actually have fired across the chunks, or the guard is
    // vacuous: count mutated R2 (0x80 set) that are still mapped (0x4 clear) in
    // the CLI SAM. cli is sorted, so recompute keys unsorted here would be
    // fiddly; instead re-run the reduction over the raw SAM.
    let out = Command::new(&bwa)
        .args([
            "mem".as_ref(),
            "-t".as_ref(),
            "1".as_ref(),
            "-K".as_ref(),
            "100000000".as_ref(),
            "-p".as_ref(),
            prefix_os.as_ref(),
            fq.as_os_str(),
        ] as [&std::ffi::OsStr; 8])
        .output()
        .expect("run bwa-mem3 mem");
    let mapped_r2 = String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter(|l| !l.starts_with('@'))
        .filter(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let flag: u16 = f[1].parse().unwrap_or(0);
            // read2 (0x80) and mapped (0x4 clear)
            flag & 0x80 != 0 && flag & 0x4 == 0
        })
        .count();
    assert!(
        mapped_r2 > 100,
        "mate rescue did not fire across chunks ({mapped_r2} mapped R2)"
    );

    // SAFETY: `opts`/`idx` are the live owned handles; each freed once and not
    // used afterward.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}
