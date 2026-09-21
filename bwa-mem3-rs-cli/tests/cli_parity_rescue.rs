//! CLI parity over a below-`opt->T` mate-rescue fixture.
//!
//! When one mate maps cleanly and the other has only a short (< `opt->T`)
//! alignment near it, upstream's paired branch emits that rescued mate as the
//! paired primary *unconditionally* — the `opt->T` gate applies only to the
//! no-pairing path and to ALT supplementaries (bwamem_pair.cpp:1279-1310). The
//! shim previously ran the general `mem_reg2sam` emit loop for the paired case
//! too, which dropped a below-T rescued primary to an unmapped record. No
//! fixture with full-length substring reads reaches this: every read maps
//! comfortably above `-T 30`. This fixture pins it by giving one mate a 25 bp
//! true match plus a random tail (score ~25 < 30), so it survives only through
//! pairing, and asserts byte parity with `bwa-mem3 mem`.
//!
//! Requires `bwa-mem3` + `samtools`; set `BWA_MEM3_BIN` if not on PATH. Skips
//! gracefully otherwise.

mod common;

use common::{
    random_dna, require_bwa_mem3, require_samtools, revcomp, samtools_view, setup_ref_index, Rng,
};
use std::process::Command;

fn norm_records(lines: &[String]) -> Vec<String> {
    let mut out: Vec<String> = lines
        .iter()
        .filter(|l| !l.starts_with('@'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let mut cols: Vec<String> = f[..11].iter().map(|s| (*s).to_string()).collect();
            let mut tags: Vec<&str> = f[11..].to_vec();
            tags.sort_unstable();
            cols.push(tags.join("\t"));
            cols.join("\t")
        })
        .collect();
    out.sort();
    out
}

#[test]
fn cli_parity_below_t_paired_rescue() {
    let Some(bwa) = require_bwa_mem3() else {
        return;
    };
    if !require_samtools() {
        return;
    }

    let mut rng = Rng(0x0BE1_0F7E_5C0E_0001);
    let chr = random_dna(&mut rng, 6000).into_bytes();

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let ref_fa = setup_ref_index(dir, &bwa, "chr", &chr);

    let (read_len, insert) = (150usize, 400usize);
    let mut r1: Vec<(String, Vec<u8>)> = Vec::new();
    let mut r2: Vec<(String, Vec<u8>)> = Vec::new();

    // Concordant pairs to establish the FR insert-size model that mate rescue
    // reads off.
    for i in 0..200 {
        let start = 100 + i * 25;
        if start + insert > chr.len() {
            break;
        }
        r1.push((format!("ok{i}"), chr[start..start + read_len].to_vec()));
        r2.push((
            format!("ok{i}"),
            revcomp(&chr[start + insert - read_len..start + insert]),
        ));
    }

    // Rescue pairs: R1 maps full length; R2 carries only the fragment's last
    // 25 bp (score ~25, below the default -T 30) followed by a random tail, so
    // R2 aligns only as a short soft-clipped block near R1 and survives solely
    // as the paired primary. 25 bp > the 19 bp min seed, so it seeds.
    let match_len = 25usize;
    for i in 0..4 {
        let start = 500 + i * 700;
        let r1_read = chr[start..start + read_len].to_vec();
        let true_tail = revcomp(&chr[start + insert - match_len..start + insert]);
        let mut r2_read = true_tail;
        r2_read.extend(random_dna(&mut rng, read_len - match_len).into_bytes());
        r1.push((format!("rescue{i}"), r1_read));
        r2.push((format!("rescue{i}"), r2_read));
    }

    let r1_fq = dir.join("r1.fq");
    let r2_fq = dir.join("r2.fq");
    common::write_fastq(&r1_fq, &r1);
    common::write_fastq(&r2_fq, &r2);

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

    let cli_lines = samtools_view(&cli_bam);
    let rs_lines = samtools_view(&rs_bam);

    // Vacuity guard: the reference aligner must actually emit at least one
    // rescue mate as a MAPPED, soft-clipped, below-T record — otherwise the
    // paired below-T path is not being exercised and the parity check is empty.
    let rescued_mapped = cli_lines
        .iter()
        .filter(|l| !l.starts_with('@'))
        .filter(|l| {
            l.split('\t')
                .next()
                .is_some_and(|q| q.starts_with("rescue"))
        })
        .any(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let flag: u32 = f[1].parse().unwrap_or(0);
            let mapped = flag & 0x4 == 0;
            let short_block = f[5].contains('S') && !f[5].starts_with("150M");
            let below_t = f[11..].iter().any(|t| {
                t.strip_prefix("AS:i:")
                    .and_then(|v| v.parse::<i32>().ok())
                    .is_some_and(|s| s < 30)
            });
            mapped && short_block && below_t
        });
    assert!(
        rescued_mapped,
        "fixture did not produce a mapped, soft-clipped, below-T (AS<30) rescue \
         mate from the reference aligner; the paired below-T path is untested"
    );

    let cli = norm_records(&cli_lines);
    let rs = norm_records(&rs_lines);
    assert_eq!(
        rs.len(),
        cli.len(),
        "record count mismatch: bwa-rs={} CLI={}",
        rs.len(),
        cli.len()
    );
    let diffs: Vec<String> = rs
        .iter()
        .zip(&cli)
        .filter(|(a, b)| a != b)
        .map(|(a, b)| format!("  bwa-rs: {a}\n  CLI:    {b}"))
        .collect();
    assert!(
        diffs.is_empty(),
        "{} of {} records diverge from bwa-mem3 mem on the rescue fixture:\n{}",
        diffs.len(),
        cli.len(),
        diffs.iter().take(8).cloned().collect::<Vec<_>>().join("\n")
    );
}
