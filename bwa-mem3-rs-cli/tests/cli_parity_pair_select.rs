//! Non-meth CLI-parity over a reference that forces the PAIRED-branch
//! primary-selection path (`z[k] != 0`) in `pair_and_emit`.
//!
//! The reference embeds a 160 bp motif twice: copy A verbatim and copy B with a
//! couple of point mutations. R1 is drawn from copy A, so its single-end best
//! hit is copy A (index 0) and copy B is a lower-scoring secondary. R2 is
//! anchored uniquely just downstream of copy B, so `mem_pair` prefers the copy-B
//! placement for R1 — i.e. the paired-selected region `z[0]` is the *secondary*
//! copy-B hit, not the SE-best copy-A hit at index 0.
//!
//! Upstream `mem_sam_pe`'s paired block emits `a[k].a[z[k]]` as the primary and
//! folds the switched-away copy-A hit into its `XA:Z`. The shim's unified
//! `mem_reg2sam`-style emission keys off `secondary` alone, so without the
//! secondary_all-fold suppression it would surface copy A as an extra primary
//! and demote the copy-B pair-primary to a `0x800` supplementary — diverging
//! from `bwa-mem3 mem`. This pins that path to byte-for-byte parity.
//!
//! Requires `bwa-mem3` + `samtools`; set `BWA_MEM3_BIN` if not on PATH. Skips
//! gracefully otherwise.

mod common;

use common::{random_dna, record_key_fields, samtools_view, Rng};
use std::process::Command;

#[test]
fn cli_parity_paired_primary_selection() {
    let Some(bwa) = common::find_bwa_mem3() else {
        eprintln!("skip: bwa-mem3 not on PATH (set BWA_MEM3_BIN)");
        return;
    };
    if !common::have_samtools() {
        eprintln!("skip: samtools not on PATH");
        return;
    }

    // Reference: bg | MOTIF (copy A) | bg | MOTIF' (copy B) | bg.
    // Copy B differs from copy A at two interior positions, so a read from copy
    // A aligns to copy A exactly (higher score, index 0) and to copy B with two
    // mismatches (lower score, secondary).
    let mut rng = Rng(0x9a3f_c105);
    let motif: Vec<u8> = random_dna(&mut rng, 160).into_bytes();
    let mut motif_b = motif.clone();
    // Two point mutations well inside the read window [20, 120): flip to a base
    // that differs from the original at each site.
    for &pos in &[55usize, 95usize] {
        motif_b[pos] = match motif_b[pos] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }
    let bg1 = random_dna(&mut rng, 1200).into_bytes();
    let bg2 = random_dna(&mut rng, 1200).into_bytes();
    let bg3 = random_dna(&mut rng, 1200).into_bytes();
    let mut seq = Vec::new();
    seq.extend_from_slice(&bg1);
    let copy_a_start = seq.len();
    seq.extend_from_slice(&motif);
    seq.extend_from_slice(&bg2);
    let copy_b_start = seq.len();
    seq.extend_from_slice(&motif_b);
    seq.extend_from_slice(&bg3);
    let seq: String = String::from_utf8(seq).unwrap();
    let bytes = seq.as_bytes();

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let ref_fa = common::setup_phix_index(dir, &bwa, &seq);

    let (read_len, insert) = (100usize, 260usize);
    let mut r1_reads = Vec::new();
    let mut r2_reads = Vec::new();

    // Seed pestat with uniquely-mapping FR pairs drawn entirely from the unique
    // bg1 tract, with a spread of insert sizes centered on `insert`. The spread
    // matters: `bwa-mem3`'s proper-pair window is derived from the insert-size
    // std.dev, so a zero-variance seed (all inserts identical) yields a
    // single-point window that the copy-B placement (insert `insert`) can fall
    // just outside of. A ~240..280 spread gives a window that comfortably
    // contains `insert`, so the copy-B pairing is proper. Without any unique
    // pairs there is no insert estimate at all and `mem_pair` never engages.
    for i in 0..120 {
        let start = 100 + (i % 600);
        let ins = 240 + (i % 41); // 240..280, mean ~260
        let r1 = bytes[start..start + read_len].to_vec();
        let r2 = common::revcomp(&bytes[start + ins - read_len..start + ins]);
        r1_reads.push((format!("u{i}"), r1));
        r2_reads.push((format!("u{i}"), r2));
    }

    // R1 sits inside copy A (multi-maps to A and B); R2 is the revcomp of a
    // window `insert` downstream of copy B, landing in the unique bg3 tract so it
    // anchors the pair to copy B. With pes FR ~= insert, the copy-B placement is
    // proper (insert ~260) while copy A is far out of range (~1620), so
    // `mem_pair` selects the copy-B secondary as z[0] != 0. `read_len <= 160`
    // keeps R1 fully in the motif.
    for i in 0..80 {
        let off = 10 + (i % 30); // vary the window but keep both mutations covered
        let a1 = copy_a_start + off;
        let r1 = bytes[a1..a1 + read_len].to_vec();
        // Anchor the mate downstream of copy B: R2 revcomp of the window at the
        // copy-B placement's insert distance.
        let b1 = copy_b_start + off;
        let r2 = common::revcomp(&bytes[b1 + insert - read_len..b1 + insert]);
        r1_reads.push((format!("p{i}"), r1));
        r2_reads.push((format!("p{i}"), r2));
    }
    let r1_fq = dir.join("r1.fq");
    let r2_fq = dir.join("r2.fq");
    common::write_fastq(&r1_fq, &r1_reads);
    common::write_fastq(&r2_fq, &r2_reads);

    // bwa-rs mem -> BAM.
    let rs_bam = dir.join("rs.bam");
    let status = Command::new(common::cli_bin())
        .args(["mem"])
        .arg(&ref_fa)
        .arg(&r1_fq)
        .arg(&r2_fq)
        .arg("-o")
        .arg(&rs_bam)
        .status()
        .expect("run bwa-rs mem");
    assert!(status.success(), "bwa-rs mem failed");

    // bwa-mem3 mem -> BAM (reference).
    let cli_out = Command::new(&bwa)
        .args(["mem", "-t", "1"])
        .arg(&ref_fa)
        .arg(&r1_fq)
        .arg(&r2_fq)
        .output()
        .expect("run bwa-mem3 mem");
    assert!(cli_out.status.success(), "bwa-mem3 mem failed");
    let cli_bam = dir.join("cli.bam");
    std::fs::write(&cli_bam, &cli_out.stdout).unwrap();

    let rs_lines = samtools_view(&rs_bam);
    let cli_lines = samtools_view(&cli_bam);

    // The paired-selection path must actually be exercised: the R1 primary must
    // land on copy B (1-based pos copy_b_start+off+1 range), with copy A folded
    // into its XA:Z. Detect it in the CLI output via a primary R1 record (flag
    // has 0x40 first-in-pair, no 0x100/0x800) whose POS is in the copy-B window
    // while its XA:Z references the copy-A window.
    let copy_a_lo = copy_a_start + 1;
    let copy_a_hi = copy_a_start + 160;
    let copy_b_lo = copy_b_start + 1;
    let copy_b_hi = copy_b_start + 160;
    let mut exercised = false;
    for l in &cli_lines {
        let f: Vec<&str> = l.split('\t').collect();
        let flag: u32 = f[1].parse().unwrap();
        let is_r1_primary = (flag & 0x40) != 0 && (flag & 0x900) == 0;
        if !is_r1_primary {
            continue;
        }
        let pos: usize = f[3].parse().unwrap();
        let in_b = pos >= copy_b_lo && pos <= copy_b_hi;
        let xa_has_a = l.contains("\tXA:Z:")
            && f[11..].iter().any(|t| {
                t.strip_prefix("XA:Z:").is_some_and(|xa| {
                    xa.split(';').filter(|e| !e.is_empty()).any(|e| {
                        // entry: rname,±pos,cigar,NM
                        e.split(',').nth(1).is_some_and(|p| {
                            let ap: usize = p.trim_start_matches(['+', '-']).parse().unwrap_or(0);
                            ap >= copy_a_lo && ap <= copy_a_hi
                        })
                    })
                })
            });
        if in_b && xa_has_a {
            exercised = true;
            break;
        }
    }
    assert!(
        exercised,
        "fixture did not exercise the paired z[k]!=0 path (no R1 primary on \
         copy B with copy A folded into XA:Z) — the test would not cover the fix"
    );

    // Sorted vectors (not sets): preserve record multiplicity so a duplicated or
    // missing record diverges even when the distinct key set matches.
    let mut cli_recs: Vec<String> = cli_lines.iter().map(|l| record_key_fields(l)).collect();
    let mut rs_recs: Vec<String> = rs_lines.iter().map(|l| record_key_fields(l)).collect();
    cli_recs.sort();
    rs_recs.sort();
    assert_eq!(
        cli_recs, rs_recs,
        "bwa-rs records diverge from bwa-mem3 CLI on the paired primary-selection path"
    );
}
