//! `bwa-rs` CLI — minimal `bwa-mem3 mem`-equivalent wrapper.
//!
//! Reads paired FASTQ(s), aligns via bwa-rs's `align_batch`, writes a
//! BGZF-compressed BAM stream. Output structure is proper BAM: magic +
//! header text + contig table + concatenated packed records + BGZF EOF.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use bwa_mem3_rs::{align_batch, shm, BwaIndex, MemOpts, ReadPair};
use clap::{Parser, Subcommand};
use flate2::read::MultiGzDecoder;
use noodles_bgzf as bgzf;

#[derive(Parser, Debug)]
#[command(
    name = "bwa-rs",
    version,
    about = "Thin CLI over bwa-rs: align paired FASTQ to BAM"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand, Debug)]
enum Cmd {
    /// Align paired-end reads to a reference and emit BAM.
    Mem {
        /// Prefix of a prebuilt bwa-mem3 index (e.g. `ref.fa`).
        prefix: PathBuf,
        /// R1 FASTQ (gz-compressed inputs supported).
        r1: PathBuf,
        /// R2 FASTQ.
        r2: PathBuf,
        /// Output BAM (stdout if omitted).
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
        /// Pairs per alignment batch.
        #[arg(long, default_value_t = 1024)]
        batch_size: usize,
        /// Minimum seed length (`-k`).
        #[arg(short = 'k', long)]
        min_seed_len: Option<i32>,
    },
    /// Manage indexes pinned in POSIX shared memory.
    Shm {
        #[command(subcommand)]
        action: ShmAction,
    },
}

#[derive(Subcommand, Debug)]
enum ShmAction {
    /// Stage an index into shared memory so subsequent loads attach.
    #[command(alias = "stage")]
    Load {
        /// Prefix of a prebuilt bwa-mem3 index (e.g. `ref.fa`).
        prefix: PathBuf,
    },
    /// List staged indexes (one `<basename>\t<bytes>` line per entry).
    List,
    /// Drop every staged segment.
    Drop,
    /// Print whether `<prefix>` is staged. Exits 0 if staged, 1 if not.
    Status {
        /// Prefix of a prebuilt bwa-mem3 index (e.g. `ref.fa`).
        prefix: PathBuf,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Mem {
            prefix,
            r1,
            r2,
            output,
            batch_size,
            min_seed_len,
        } => run_mem(
            &prefix,
            &r1,
            &r2,
            output.as_deref(),
            batch_size,
            min_seed_len,
        ),
        Cmd::Shm { action } => run_shm(action),
    }
}

fn run_shm(action: ShmAction) -> Result<()> {
    match action {
        ShmAction::Load { prefix } => shm::stage(&prefix)
            .with_context(|| format!("staging index {prefix:?} in shared memory")),
        ShmAction::List => shm::list().context("listing staged shared-memory segments"),
        ShmAction::Drop => shm::destroy().context("dropping staged shared-memory segments"),
        ShmAction::Status { prefix } => {
            let staged = shm::is_staged(&prefix)
                .with_context(|| format!("probing shm registry for {prefix:?}"))?;
            println!("{}", if staged { "staged" } else { "not staged" });
            std::process::exit(if staged { 0 } else { 1 });
        }
    }
}

fn run_mem(
    prefix: &std::path::Path,
    r1_path: &std::path::Path,
    r2_path: &std::path::Path,
    output: Option<&std::path::Path>,
    batch_size: usize,
    min_seed_len: Option<i32>,
) -> Result<()> {
    let idx = BwaIndex::load(prefix).with_context(|| format!("loading index {prefix:?}"))?;
    let mut opts = MemOpts::new()?;
    opts.set_pe(true);
    if let Some(k) = min_seed_len {
        opts.set_min_seed_len(k);
    }

    let r1 = open_fastq(r1_path)?;
    let r2 = open_fastq(r2_path)?;
    let mut r1_it = noodles_fastq::io::Reader::new(r1);
    let mut r2_it = noodles_fastq::io::Reader::new(r2);

    let out: Box<dyn Write> = match output {
        Some(p) => Box::new(File::create(p).with_context(|| format!("create {p:?}"))?),
        None => Box::new(io::stdout().lock()),
    };
    let mut bgzf_writer = bgzf::io::Writer::new(out);

    write_bam_header(&mut bgzf_writer, &idx)?;

    // Buffers: we accumulate full Record structures so lifetime-free
    // ownership is easy. Each inner Vec holds name + seq + qual.
    let mut pending: Vec<PendingPair> = Vec::with_capacity(batch_size);
    let mut total_recs: u64 = 0;

    loop {
        pending.clear();
        for _ in 0..batch_size {
            let mut rec1 = noodles_fastq::Record::default();
            let mut rec2 = noodles_fastq::Record::default();
            let got1 = r1_it.read_record(&mut rec1).context("reading R1")?;
            let got2 = r2_it.read_record(&mut rec2).context("reading R2")?;
            match (got1, got2) {
                (0, 0) => break,
                (0, _) | (_, 0) => anyhow::bail!("R1 and R2 have unequal lengths"),
                _ => {}
            }
            pending.push(PendingPair::from_records(&rec1, &rec2));
        }
        if pending.is_empty() {
            break;
        }

        // Build ReadPair views into the pending Vec.
        let pairs: Vec<ReadPair<'_>> = pending
            .iter()
            .map(|p| ReadPair {
                name_r1: &p.name_r1,
                seq_r1: &p.seq_r1,
                qual_r1: Some(&p.qual_r1),
                name_r2: &p.name_r2,
                seq_r2: &p.seq_r2,
                qual_r2: Some(&p.qual_r2),
            })
            .collect();

        let (aln, _) = align_batch(&idx, &opts, &pairs, None).context("align_batch")?;
        for r in aln.iter() {
            bgzf_writer.write_all(r.bytes).context("writing record")?;
            total_recs += 1;
        }
    }

    bgzf_writer.finish().context("finalizing BGZF stream")?;
    eprintln!("bwa-rs: wrote {total_recs} records");
    Ok(())
}

struct PendingPair {
    name_r1: Vec<u8>,
    seq_r1: Vec<u8>,
    qual_r1: Vec<u8>,
    name_r2: Vec<u8>,
    seq_r2: Vec<u8>,
    qual_r2: Vec<u8>,
}

impl PendingPair {
    fn from_records(r1: &noodles_fastq::Record, r2: &noodles_fastq::Record) -> Self {
        Self {
            name_r1: r1.name().to_vec(),
            seq_r1: r1.sequence().to_vec(),
            qual_r1: r1.quality_scores().to_vec(),
            name_r2: r2.name().to_vec(),
            seq_r2: r2.sequence().to_vec(),
            qual_r2: r2.quality_scores().to_vec(),
        }
    }
}

/// Open a FASTQ file, transparently decompressing `.gz` / `.fq.gz`.
fn open_fastq(path: &std::path::Path) -> Result<Box<dyn BufRead>> {
    let f = File::open(path).with_context(|| format!("open {path:?}"))?;
    let is_gz = path
        .extension()
        .and_then(|s| s.to_str())
        .is_some_and(|s| s.eq_ignore_ascii_case("gz"));
    if is_gz {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

/// Emit a minimal BAM header: magic + l_text + @HD + @SQ lines +
/// n_ref + per-ref (l_name + name + NUL + l_ref).
fn write_bam_header<W: Write>(w: &mut W, idx: &BwaIndex) -> Result<()> {
    // Build the SAM header text.
    let mut text = String::new();
    text.push_str("@HD\tVN:1.6\tSO:unknown\n");
    for (name, len) in idx.contigs() {
        text.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
    }
    text.push_str("@PG\tID:bwa-rs\tPN:bwa-rs\tVN:");
    text.push_str(env!("CARGO_PKG_VERSION"));
    text.push('\n');

    let text_bytes = text.as_bytes();
    // magic
    w.write_all(b"BAM\x01")?;
    // l_text
    w.write_all(&u32::try_from(text_bytes.len())?.to_le_bytes())?;
    w.write_all(text_bytes)?;

    // n_ref
    let n_ref = u32::try_from(idx.n_contigs())?;
    w.write_all(&n_ref.to_le_bytes())?;
    for (name, len) in idx.contigs() {
        let name_bytes = name.as_bytes();
        // l_name = strlen + NUL
        let l_name = u32::try_from(name_bytes.len() + 1)?;
        w.write_all(&l_name.to_le_bytes())?;
        w.write_all(name_bytes)?;
        w.write_all(b"\0")?;
        // l_ref
        w.write_all(&u32::try_from(len)?.to_le_bytes())?;
    }

    Ok(())
}

// Silence unused warnings when Read isn't touched directly (noodles handles it).
#[allow(dead_code)]
fn _ensure_read_trait_in_scope<R: Read>(_: R) {}
