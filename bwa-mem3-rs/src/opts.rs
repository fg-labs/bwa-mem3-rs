//! Alignment options and paired-end insert-size model.
//!
//! [`MemOpts`] wraps bwa-mem3's `mem_opt_t` (POD). Fields are directly
//! mutable via typed setters. [`MemPeStat`] wraps a 4-orientation
//! `mem_pestat_t[4]` array used for paired-end insert-size modeling.

use bwa_mem3_sys as sys;

use crate::error::{shim_err, Result};

/// Read-type preset flags (bwa-mem3 CLI `-x`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Pacbio,
    Ont2d,
    Intractg,
}

/// Seed emission order (`--seed-order`), mirroring bwa-mem3's `seed_order_t`.
///
/// [`SeedOrder::Off`] is the default and is byte-identical to the baseline;
/// the other modes reorder seeds before chaining and are opt-in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeedOrder {
    /// Default resolve order — byte-identical to the baseline.
    Off = 0,
    /// Emit seeds longest-first across the whole read.
    GlobalLongest = 1,
    /// Emit seeds longest-first within each local cluster.
    LocalLongest = 2,
    /// Order by how many shorter seeds each seed absorbs.
    AbsorbCount = 3,
    /// Prefer seeds that absorb the most others first.
    MostAbsorb = 4,
}

impl TryFrom<i32> for SeedOrder {
    type Error = crate::Error;

    /// Map the raw `seed_order_t` value to a [`SeedOrder`].
    ///
    /// Returns [`Error::UnrecognizedEnum`](crate::Error::UnrecognizedEnum) for
    /// a value bwa-mem3 does not define. This deliberately does *not* fall
    /// back to [`SeedOrder::Off`]: a silent fallback made a new upstream mode
    /// read back as the default, which is indistinguishable from a
    /// byte-identical baseline run.
    fn try_from(v: i32) -> std::result::Result<Self, Self::Error> {
        match v {
            0 => Ok(SeedOrder::Off),
            1 => Ok(SeedOrder::GlobalLongest),
            2 => Ok(SeedOrder::LocalLongest),
            3 => Ok(SeedOrder::AbsorbCount),
            4 => Ok(SeedOrder::MostAbsorb),
            other => Err(crate::Error::UnrecognizedEnum {
                kind: "seed_emit_order",
                value: other,
            }),
        }
    }
}

/// Bisulfite scoring mode (`--meth-scoring`), mirroring bwa-mem3's
/// `mem_meth_scoring`. Only meaningful under [`MemOpts::set_meth`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MethScoring {
    /// C↔T and G↔A interchangeable (both matrix cells freed) — reproduces
    /// bwameth. bwa-mem3's default through v0.16.0; select it explicitly to
    /// keep that output.
    Collapsed = 0,
    /// Only the bisulfite-conversion cell is freed; the mirror stays a
    /// mismatch, so genuine variants still score as mismatches (variant-aware,
    /// truthful NM/MD). bwa-mem3's default as of v0.17.0.
    Genomic = 1,
    /// Like [`Genomic`](Self::Genomic), but the conversion cell scores `0` —
    /// tolerated rather than rewarded. Added in bwa-mem3 v0.9.0 and made the
    /// automatic default for TAPS chemistry, where conversions are ~20-30×
    /// rarer than under em-seq so a full match over-credits spurious C→T
    /// alignments.
    Neutral = 2,
}

impl TryFrom<i32> for MethScoring {
    type Error = crate::Error;

    /// Map the raw `mem_meth_scoring` value to a [`MethScoring`].
    ///
    /// Returns [`Error::UnrecognizedEnum`](crate::Error::UnrecognizedEnum) for
    /// a value bwa-mem3 does not define. This deliberately does not fall back
    /// to a default: a silent fallback is what made
    /// `MEM_METH_SCORING_NEUTRAL = 2` read back as
    /// [`MethScoring::Collapsed`] before it was mapped here.
    fn try_from(v: i32) -> std::result::Result<Self, Self::Error> {
        match v {
            0 => Ok(MethScoring::Collapsed),
            1 => Ok(MethScoring::Genomic),
            2 => Ok(MethScoring::Neutral),
            other => Err(crate::Error::UnrecognizedEnum {
                kind: "meth_scoring",
                value: other,
            }),
        }
    }
}

/// Methylation chemistry (`--meth=emseq|taps`), mirroring bwa-mem3's
/// `meth_chem_t`. Only meaningful under [`MemOpts::set_meth`].
///
/// Both chemistries turn a C into a T and align the same way, so the
/// chemistry only changes which bases mean "methylated":
///
/// - [`EmSeq`](Self::EmSeq) (bisulfite or EM-seq): *unmethylated* C converts.
///   A retained C is methylated.
/// - [`Taps`](Self::Taps): *methylated* C converts. A retained C is
///   unmethylated.
///
/// That flips the `XM:Z` methylation call. TAPS also changes one default that
/// [`MemOpts::apply_meth_defaults`] applies: [`MethScoring::Neutral`] instead
/// of [`MethScoring::Genomic`] (TAPS conversions are ~20-30× rarer, so a
/// full match over-credits spurious C→T alignments). The SPEC30 seed-pruning
/// default applies to both chemistries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MethChem {
    /// Bisulfite or EM-seq: unmethylated C converts to T. bwa-mem3's default.
    EmSeq = 0,
    /// TET-assisted pyridine borane sequencing: methylated C converts to T.
    Taps = 1,
}

impl TryFrom<i32> for MethChem {
    type Error = crate::Error;

    /// Map the raw `meth_chem_t` value to a [`MethChem`].
    ///
    /// Returns [`Error::UnrecognizedEnum`](crate::Error::UnrecognizedEnum) for
    /// a value bwa-mem3 does not define, rather than falling back to
    /// [`MethChem::EmSeq`], which would flip every `XM:Z` call of a new
    /// chemistry without a trace.
    fn try_from(v: i32) -> std::result::Result<Self, Self::Error> {
        match v {
            0 => Ok(MethChem::EmSeq),
            1 => Ok(MethChem::Taps),
            other => Err(crate::Error::UnrecognizedEnum {
                kind: "meth_chem",
                value: other,
            }),
        }
    }
}

/// Bisulfite seed pruning (`--meth-seed-prune`), mirroring bwa-mem3's
/// `mem_meth_seed_prune`. Only meaningful under [`MemOpts::set_meth`].
///
/// The 3-letter alphabet over-seeds with short, high-multiplicity SMEMs.
/// Pruning them before SA resolution is much faster and **changes output**
/// (XS, XA, MAPQ, and occasionally placement), so it must match the CLI's
/// setting for byte parity. [`MemOpts::apply_meth_defaults`] selects
/// [`Spec30`](Self::Spec30), as `bwa-mem3 mem --meth` does.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MethSeedPrune {
    /// Keep every SMEM.
    Off = 0,
    /// Per read, if the longest SMEM is at least 30 bp, keep an SMEM only if
    /// it is at least 25 bp, or unique and at least 22 bp. The `--meth`
    /// default.
    Spec30 = 1,
    /// Drop an SMEM shorter than 25 bp with more than one hit.
    Baseline = 2,
}

impl TryFrom<i32> for MethSeedPrune {
    type Error = crate::Error;

    /// Map the raw `mem_meth_seed_prune` value to a [`MethSeedPrune`].
    ///
    /// Returns [`Error::UnrecognizedEnum`](crate::Error::UnrecognizedEnum) for
    /// a value bwa-mem3 does not define.
    fn try_from(v: i32) -> std::result::Result<Self, Self::Error> {
        match v {
            0 => Ok(MethSeedPrune::Off),
            1 => Ok(MethSeedPrune::Spec30),
            2 => Ok(MethSeedPrune::Baseline),
            other => Err(crate::Error::UnrecognizedEnum {
                kind: "meth_seed_prune",
                value: other,
            }),
        }
    }
}

/// Where a methylation conversion is free (`--meth-context`), mirroring
/// bwa-mem3's `meth_ctx_t` (v0.16.0+). Only meaningful under
/// [`MemOpts::set_meth`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MethContext {
    /// A conversion (C→T on the OT strand, G→A on OB) is free at every C
    /// (G). bwa-mem3's default.
    All = 0,
    /// A conversion (C→T on OT, G→A on OB) is free only at a CpG, so
    /// off-target conversion, damage or a variant at any other C (G) scores as
    /// a mismatch. Requires
    /// [`MethChem::Taps`] (which covers Illumina 5-base) and a scoring mode
    /// other than [`MethScoring::Collapsed`];
    /// [`MemOpts::apply_meth_defaults`] rejects any other combination.
    Cpg = 1,
}

impl TryFrom<i32> for MethContext {
    type Error = crate::Error;

    /// Map the raw `meth_ctx_t` value to a [`MethContext`].
    ///
    /// Returns [`Error::UnrecognizedEnum`](crate::Error::UnrecognizedEnum) for
    /// a value bwa-mem3 does not define.
    fn try_from(v: i32) -> std::result::Result<Self, Self::Error> {
        match v {
            0 => Ok(MethContext::All),
            1 => Ok(MethContext::Cpg),
            other => Err(crate::Error::UnrecognizedEnum {
                kind: "meth_context",
                value: other,
            }),
        }
    }
}

/// How the `MD:Z` tag is written (`--md`), mirroring bwa-mem3's `mem_md_mode`
/// (v0.17.0+).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MdMode {
    /// No `MD:Z` tag.
    None = 0,
    /// The literal `MD`: every reference base the read differs from. The
    /// default.
    Full = 1,
    /// Only the mismatches the scoring matrix penalises, the ones `NM`
    /// counts. Under `--meth` that leaves out the differences the matrix
    /// frees, not necessarily every conversion: with [`MethContext::Cpg`] an
    /// off-CpG C→T is still penalised, so it counts toward `NM` and appears
    /// here. Without `--meth` every mismatch is penalised and this equals
    /// [`Full`](Self::Full).
    Short = 2,
}

impl TryFrom<i32> for MdMode {
    type Error = crate::Error;

    /// Map the raw `mem_md_mode` value to an [`MdMode`].
    ///
    /// Returns [`Error::UnrecognizedEnum`](crate::Error::UnrecognizedEnum) for
    /// a value bwa-mem3 does not define.
    fn try_from(v: i32) -> std::result::Result<Self, Self::Error> {
        match v {
            0 => Ok(MdMode::None),
            1 => Ok(MdMode::Full),
            2 => Ok(MdMode::Short),
            other => Err(crate::Error::UnrecognizedEnum {
                kind: "md_mode",
                value: other,
            }),
        }
    }
}

/// A [`MemOpts::phase_snapshot`]: the alignment-bearing option bytes a
/// split-phase handle was seeded with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PhaseOpts(Box<[u8]>);

impl PhaseOpts {
    /// `Ok` if `self` (the options a call was given) equals `seeded`, the
    /// snapshot the regions it consumes were seeded with.
    pub(crate) fn check_against(&self, seeded: &PhaseOpts, what: &str) -> Result<()> {
        if self == seeded {
            Ok(())
        } else {
            Err(crate::Error::InvalidOpts(format!(
                "{what}: the options differ from those the reads were seeded with; every \
                 phase of a batch or cohort must use the same alignment options (only the \
                 MD mode and the paired-end flag may change)"
            )))
        }
    }
}

/// bwa-mem3 alignment options (`mem_opt_t`).
pub struct MemOpts {
    pub(crate) handle: *mut sys::mem_opt_t,
    /// The bisulfite knobs the caller set explicitly, as
    /// `sys::BWA_SHIM_METH_SET_*` bits. [`MemOpts::apply_meth_defaults`]
    /// leaves these alone, as `bwa-mem3 mem --meth` leaves a flag the user
    /// passed.
    meth_explicit: u32,
}

impl MemOpts {
    /// Create options with bwa-mem3 defaults (via `mem_opt_init`).
    pub fn new() -> Result<Self> {
        let handle = unsafe { sys::bwa_shim_opts_new() };
        if handle.is_null() {
            return Err(shim_err("opts_new"));
        }
        Ok(MemOpts {
            handle,
            meth_explicit: 0,
        })
    }

    /// The `@HD` header line the active output-compatibility target calls for
    /// (no trailing newline), or `None` when that target emits no `@HD`.
    ///
    /// A caller writing its own SAM/BAM header should prefer
    /// [`BwaIndex::sam_header`](crate::BwaIndex::sam_header), which returns
    /// the whole `@HD` + `@SQ` block from bwa-mem3's own builder (including
    /// `AH:*` on ALT contigs). This accessor is for a caller that needs only
    /// `@HD`, in place of a literal: bwa-mem3's `compat_target_t` is the
    /// single place that decides output shaping, and upstream folded three
    /// hand-written `@HD` literals into it after they drifted apart
    /// (fg-labs/bwa-mem3#288), so a literal here would reintroduce exactly
    /// that.
    ///
    /// The default target (`off`, bwa-mem3's native output) emits one; that is
    /// the only target reachable today, since this crate exposes no `--compat`
    /// selector.
    pub fn compat_hd_line(&self) -> Option<&str> {
        let p = unsafe { sys::bwa_shim_compat_hd_line(self.handle) };
        if p.is_null() {
            return None;
        }
        // Static storage owned by bwa-mem3's compat table; valid for the
        // process lifetime, so borrowing it for `&self` is sound.
        unsafe { std::ffi::CStr::from_ptr(p) }.to_str().ok()
    }

    /// Apply what `bwa-mem3 mem --meth` sets for bisulfite alignment, so this
    /// crate's `--meth` output matches the CLI's:
    ///
    /// - bwameth's compatibility bundle: `-L 10,10`, `-U 100`, `-T 40`, `-M`,
    ///   and (COLLAPSED scoring only) `-B 2`;
    /// - [`MethScoring::Neutral`] under [`MethChem::Taps`], unless
    ///   [`set_meth_scoring`](Self::set_meth_scoring) was called;
    /// - [`MethSeedPrune::Spec30`], unless
    ///   [`set_meth_seed_prune`](Self::set_meth_seed_prune) was called;
    /// - the certified adaptive band off (byte-identical either way).
    ///
    /// Like `bwa-mem3 mem`, it applies nothing unless
    /// [`set_meth(true)`](Self::set_meth) came first; only the scoring
    /// matrices are refilled.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidOpts`](crate::Error::InvalidOpts), leaving the
    /// options unchanged, for a [`MethContext::Cpg`] that `bwa-mem3 mem`
    /// refuses: without `--meth`, with [`MethChem::EmSeq`], or with
    /// [`MethScoring::Collapsed`] scoring. Under `--meth` it also refuses a
    /// [match score](Self::set_match_score) outside `1..=127`, the range
    /// `bwa-mem3 mem` accepts for `-A`: the bundle multiplies it by up to 100.
    ///
    /// Each bundle constant is **scaled by the match score** `a`, because
    /// bwameth quotes them at bwa's default `a == 1` and bwa scales every other
    /// score-derived default the same way. At the default `a == 1` the scaling
    /// is a no-op.
    ///
    /// The bundle is bwa-mem3's own `mem_opt_apply_meth_defaults`, so it cannot
    /// drift from upstream, and the scoring matrices are refilled afterwards
    /// (except while the resulting mismatch penalty is outside `0..=128`, as
    /// for the [`set_mismatch_penalty`](Self::set_mismatch_penalty) setter).
    ///
    /// # Ordering
    ///
    /// The knobs involved are not symmetric, so the order matters in both
    /// directions:
    ///
    /// **Set before** — these are *inputs*:
    /// - [`set_meth(true)`](Self::set_meth): as in `bwa-mem3 mem`, every
    ///   default here is skipped when `--meth` is off;
    /// - [`set_match_score`](Self::set_match_score), since every constant is
    ///   expressed in units of it;
    /// - [`set_meth_chem`](Self::set_meth_chem), since TAPS picks the scoring
    ///   default;
    /// - [`set_meth_scoring`](Self::set_meth_scoring), since the `-B 2` branch
    ///   keys off the resolved mode. An explicit mode always wins over the
    ///   TAPS default.
    ///
    /// **Either order** — [`set_meth_seed_prune`](Self::set_meth_seed_prune):
    /// an explicit choice is remembered and survives this call.
    ///
    /// **Set after** — these are *overwritten* unconditionally:
    /// [`set_minimum_score`](Self::set_minimum_score),
    /// [`set_clip_penalty`](Self::set_clip_penalty),
    /// [`set_unpaired_penalty`](Self::set_unpaired_penalty).
    ///
    /// **Set after, but only under [`MethScoring::Collapsed`]** —
    /// [`set_mismatch_penalty`](Self::set_mismatch_penalty). The bundle's
    /// `-B 2` applies to COLLAPSED alone; [`Genomic`](MethScoring::Genomic) and
    /// [`Neutral`](MethScoring::Neutral) are variant-aware and keep bwa's
    /// default, so under those two a `-B` set *before* the call survives it.
    pub fn apply_meth_defaults(&mut self) -> Result<&mut Self> {
        let rc = unsafe {
            sys::bwa_shim_opts_apply_meth_defaults_masked(self.handle, self.meth_explicit)
        };
        if rc != 0 {
            return Err(crate::Error::InvalidOpts(crate::error::shim_last_message()));
        }
        Ok(self)
    }

    /// Apply a bwa-mem3 `-x` preset on top of current values.
    #[must_use]
    pub fn with_mode(mut self, mode: Mode) -> Self {
        let o = unsafe { &mut *self.handle };
        match mode {
            Mode::Pacbio => {
                o.min_seed_len = 17;
                o.o_del = 1;
                o.o_ins = 1;
                o.e_del = 1;
                o.e_ins = 1;
                o.b = 1;
                o.split_factor = 10.0;
                o.T = 0;
            }
            Mode::Ont2d => {
                o.min_seed_len = 14;
                o.o_del = 1;
                o.o_ins = 1;
                o.e_del = 1;
                o.e_ins = 1;
                o.b = 1;
                o.split_factor = 20.0;
                o.T = 20;
            }
            Mode::Intractg => {
                o.b = 9;
                o.o_del = 16;
                o.o_ins = 16;
                o.pen_clip5 = 5;
                o.pen_clip3 = 5;
            }
        }
        // Every preset changes `b`, so the matrices derived from it are stale
        // until rebuilt -- the same coupling the `a`/`b` setters maintain.
        // Upstream does this too, re-running bwa_fill_scmat after parsing -x
        // (fastmap.cpp:2726, in main_mem).
        self.refill_matrices();
        self
    }

    /// Rebuild the scoring matrices from `a`/`b` (`bwa_fill_scmat`, plus the
    /// `--meth` copies), as every `a`/`b` change must. The shim skips it while
    /// `b` is outside `0..=128`: `bwa_fill_scmat` negates it as a C `int`,
    /// which overflows for `i32::MIN`, and such options cannot align anyway
    /// ([`check_effective_scoring`](Self::check_effective_scoring)). The next
    /// storable `b` rebuilds them from the current `a`.
    fn refill_matrices(&mut self) {
        // SAFETY: `handle` is a live `mem_opt_t` owned by `self`.
        unsafe { sys::bwa_shim_opts_fill_scmat(self.handle) };
    }

    /// Enable or disable paired-end mode.
    pub fn set_pe(&mut self, is_pe: bool) -> &mut Self {
        let o = unsafe { &mut *self.handle };
        if is_pe {
            o.flag |= sys::MEM_F_PE as i32;
        } else {
            o.flag &= !(sys::MEM_F_PE as i32);
        }
        self
    }

    /// Soft-clip supplementary alignments (`-Y`).
    pub fn set_soft_clip_supplementary(&mut self, v: bool) -> &mut Self {
        let o = unsafe { &mut *self.handle };
        if v {
            o.flag |= sys::MEM_F_SOFTCLIP as i32;
        } else {
            o.flag &= !(sys::MEM_F_SOFTCLIP as i32);
        }
        self
    }

    /// Mark shorter split hits as secondary rather than supplementary (`-M`).
    ///
    /// With this set, the 2nd+ emitted record of a split read carries `0x100`
    /// instead of `0x800` — Picard-compatible output. `bwa-mem3 mem` sets it
    /// under `--meth`, so anything aiming for byte-parity with the CLI in
    /// bisulfite mode needs it too (see [`MemOpts::set_meth`]).
    pub fn set_mark_split_secondary(&mut self, v: bool) -> &mut Self {
        let o = unsafe { &mut *self.handle };
        if v {
            o.flag |= sys::MEM_F_NO_MULTI as i32;
        } else {
            o.flag &= !(sys::MEM_F_NO_MULTI as i32);
        }
        self
    }

    /// Whether shorter split hits are marked secondary (`-M`).
    #[must_use]
    pub fn mark_split_secondary(&self) -> bool {
        unsafe { (*self.handle).flag & sys::MEM_F_NO_MULTI as i32 != 0 }
    }

    // ---------- field accessors ----------

    #[must_use]
    pub fn min_seed_len(&self) -> i32 {
        unsafe { (*self.handle).min_seed_len }
    }
    pub fn set_min_seed_len(&mut self, v: i32) -> &mut Self {
        unsafe {
            (*self.handle).min_seed_len = v;
        }
        self
    }

    #[must_use]
    pub fn band_width(&self) -> i32 {
        unsafe { (*self.handle).w }
    }
    pub fn set_band_width(&mut self, v: i32) -> &mut Self {
        unsafe {
            (*self.handle).w = v;
        }
        self
    }

    #[must_use]
    pub fn match_score(&self) -> i32 {
        unsafe { (*self.handle).a }
    }
    /// Set the match score (`-A`).
    ///
    /// Rebuilds the scoring matrices, which several alignment paths read
    /// instead of `a`/`b` directly — see [`Self::set_mismatch_penalty`] —
    /// unless the stored mismatch penalty is outside `0..=128`; the next
    /// storable penalty then rebuilds them from this score.
    pub fn set_match_score(&mut self, v: i32) -> &mut Self {
        unsafe {
            (*self.handle).a = v;
        }
        self.refill_matrices();
        self
    }

    #[must_use]
    pub fn mismatch_penalty(&self) -> i32 {
        unsafe { (*self.handle).b }
    }
    /// Set the mismatch penalty (`-B`).
    ///
    /// Rebuilds the scoring matrices. This is not optional bookkeeping:
    /// `a`/`b` are read directly by some paths (`bwamem.cpp:2200`, `:2296`,
    /// `:3449`, `:4130`, ...) while the Smith-Waterman kernels read only
    /// `opt->mat`, so writing the field alone leaves the two disagreeing —
    /// the score would change for MAPQ and heuristics but not for alignment
    /// itself. Upstream always pairs the write with `bwa_fill_scmat`
    /// (`fastmap.cpp:2726`, in `main_mem`).
    ///
    /// As in `bwa-mem3 mem`, every alignment entry point rejects a penalty
    /// outside `0..=128` (the kernels store its negation in 8 bits), and under
    /// [`set_meth(true)`](Self::set_meth) also `0`, with
    /// [`Error::InvalidOpts`](crate::Error::InvalidOpts). Such a penalty is
    /// stored but does not rebuild the matrices (`bwa_fill_scmat` negates it
    /// as a C `int`, which overflows for `i32::MIN`), so until a storable
    /// penalty is set again [`scoring_matrix`](Self::scoring_matrix) shows
    /// the last storable one; the options cannot align in the meantime.
    pub fn set_mismatch_penalty(&mut self, v: i32) -> &mut Self {
        unsafe {
            (*self.handle).b = v;
        }
        self.refill_matrices();
        self
    }

    /// The 5×5 substitution matrix (`opt->mat`) the Smith-Waterman kernels
    /// actually score with, row-major over `A,C,G,T,N`.
    ///
    /// Derived from `a`/`b` by `bwa_fill_scmat`: the diagonal is `+a`, the
    /// off-diagonal `-b`, and every ambiguous-base cell `-1`. Exposed because
    /// several paths read `a`/`b` directly while the kernels read only this,
    /// so it is the only way to observe whether the two agree — see
    /// [`set_mismatch_penalty`](Self::set_mismatch_penalty). While the stored
    /// penalty is outside `0..=128` it is not rebuilt, and still reflects the
    /// last storable one.
    #[must_use]
    pub fn scoring_matrix(&self) -> [i8; 25] {
        unsafe { (*self.handle).mat }
    }

    /// The per-hypothesis bisulfite matrices `(mat_ot, mat_ob)`, as rebuilt by
    /// `mem_opt_fill_meth_mat`. Only read under [`set_meth`](Self::set_meth).
    ///
    /// Each is a copy of [`scoring_matrix`](Self::scoring_matrix) with the
    /// conversion cell freed — `OT` frees ref-C/read-T, `OB` ref-G/read-A —
    /// to `+a` under `Collapsed`/`Genomic` or `0` under `Neutral`.
    /// `Collapsed` additionally frees the mirror cell, making C/T (and G/A)
    /// mutually interchangeable.
    #[must_use]
    pub fn meth_scoring_matrices(&self) -> ([i8; 25], [i8; 25]) {
        unsafe { ((*self.handle).mat_ot, (*self.handle).mat_ob) }
    }

    pub fn set_gap_open(&mut self, del: i32, ins: i32) -> &mut Self {
        unsafe {
            (*self.handle).o_del = del;
            (*self.handle).o_ins = ins;
        }
        self
    }

    pub fn set_gap_extend(&mut self, del: i32, ins: i32) -> &mut Self {
        unsafe {
            (*self.handle).e_del = del;
            (*self.handle).e_ins = ins;
        }
        self
    }

    pub fn set_clip_penalty(&mut self, five: i32, three: i32) -> &mut Self {
        unsafe {
            (*self.handle).pen_clip5 = five;
            (*self.handle).pen_clip3 = three;
        }
        self
    }

    /// The 5'/3' clipping penalties (`-L`), in the order
    /// [`set_clip_penalty`](Self::set_clip_penalty) takes them.
    #[must_use]
    pub fn clip_penalty(&self) -> (i32, i32) {
        unsafe { ((*self.handle).pen_clip5, (*self.handle).pen_clip3) }
    }

    #[must_use]
    pub fn minimum_score(&self) -> i32 {
        unsafe { (*self.handle).T }
    }
    pub fn set_minimum_score(&mut self, v: i32) -> &mut Self {
        unsafe {
            (*self.handle).T = v;
        }
        self
    }

    #[must_use]
    pub fn max_occurrences(&self) -> i32 {
        unsafe { (*self.handle).max_occ }
    }
    pub fn set_max_occurrences(&mut self, v: i32) -> &mut Self {
        unsafe {
            (*self.handle).max_occ = v;
        }
        self
    }

    pub fn set_xa_max_hits(&mut self, primary: i32, alt: i32) -> &mut Self {
        unsafe {
            (*self.handle).max_XA_hits = primary;
            (*self.handle).max_XA_hits_alt = alt;
        }
        self
    }

    pub fn set_xa_drop_ratio(&mut self, v: f32) -> &mut Self {
        unsafe {
            (*self.handle).XA_drop_ratio = v;
        }
        self
    }

    /// Set the `@RG` ID emitted on records as `RG:Z:...`. `None` clears it.
    ///
    /// **Process-wide:** this writes to bwa-mem3's global `bwa_rg_id` buffer.
    /// Aligning with different read groups from multiple threads requires
    /// serializing the `set_read_group_id` call or using separate processes.
    pub fn set_read_group_id(&mut self, id: Option<&str>) -> Result<&mut Self> {
        use std::ffi::CString;
        match id {
            Some(s) => {
                let c = CString::new(s)?;
                unsafe { bwa_mem3_sys::bwa_shim_set_rg_id(c.as_ptr()) };
            }
            None => unsafe { bwa_mem3_sys::bwa_shim_set_rg_id(std::ptr::null()) },
        }
        Ok(self)
    }

    pub fn set_unpaired_penalty(&mut self, v: i32) -> &mut Self {
        unsafe {
            (*self.handle).pen_unpaired = v;
        }
        self
    }

    /// The penalty for not pairing two mates (`-U`).
    #[must_use]
    pub fn unpaired_penalty(&self) -> i32 {
        unsafe { (*self.handle).pen_unpaired }
    }

    // ---------- bwa-mem3 0.6.0 opt-in knobs ----------
    //
    // All default to off/0 (via mem_opt_init) and, unless noted, are
    // byte-identical to the baseline when left at their defaults.

    /// Minimum seed length to extend (`--min-ext-len`): seeds shorter than
    /// this are not banded-SW extended. `0` (the default) disables the filter.
    /// Opt-in speed lever.
    #[must_use]
    pub fn min_ext_len(&self) -> i32 {
        unsafe { (*self.handle).min_ext_len }
    }
    pub fn set_min_ext_len(&mut self, v: i32) -> &mut Self {
        unsafe {
            (*self.handle).min_ext_len = v;
        }
        self
    }

    /// Cap on chains extended per read (`--max-extend-chains`): keep only the
    /// top-N by weight before banded-SW. `0` (the default) disables the cap.
    /// Opt-in speed lever; **not** byte-identical to the default.
    #[must_use]
    pub fn max_extend_chains(&self) -> i32 {
        unsafe { (*self.handle).max_extend_chains }
    }
    pub fn set_max_extend_chains(&mut self, v: i32) -> &mut Self {
        unsafe {
            (*self.handle).max_extend_chains = v;
        }
        self
    }

    /// Mate-concordant chain-retention window (`--extend-mate-concordant`),
    /// applied only when [`set_max_extend_chains`](Self::set_max_extend_chains)
    /// caps a paired-end read. `0` (the default) is off, `-1` auto-derives the
    /// window from the estimated proper-pair insert bound, and a positive value
    /// is a fixed window in bp. **Not** byte-identical to the default.
    #[must_use]
    pub fn mate_concordant_window(&self) -> i32 {
        unsafe { (*self.handle).mate_concordant_window }
    }
    pub fn set_mate_concordant_window(&mut self, v: i32) -> &mut Self {
        unsafe {
            (*self.handle).mate_concordant_window = v;
        }
        self
    }

    /// Deduplicate fully-identical SMEMs before SA expansion (`--smem-dedup`).
    /// `false` (the default) is byte-identical to the baseline.
    #[must_use]
    pub fn smem_dedup(&self) -> bool {
        unsafe { (*self.handle).smem_dedup != 0 }
    }
    pub fn set_smem_dedup(&mut self, v: bool) -> &mut Self {
        unsafe {
            (*self.handle).smem_dedup = i32::from(v);
        }
        self
    }

    /// Skip banded-SW extension of seeds contained (same diagonal) in a longer
    /// in-chain seed (`--skip-contained-ext`). `false` is the default; enabling
    /// it is byte-identical to the baseline.
    #[must_use]
    pub fn skip_contained_ext(&self) -> bool {
        unsafe { (*self.handle).skip_contained_ext != 0 }
    }
    pub fn set_skip_contained_ext(&mut self, v: bool) -> &mut Self {
        unsafe {
            (*self.handle).skip_contained_ext = i32::from(v);
        }
        self
    }

    /// Adaptive chain-geometry banding start band (`--adaptive-band`). `0` (the
    /// default) disables adaptive banding; a positive value is the starting
    /// band width (bwa-mem3's `--adaptive-band` uses 20). Long-read speed
    /// lever; a no-op on the 8-bit short-read tier.
    #[must_use]
    pub fn band_start(&self) -> i32 {
        unsafe { (*self.handle).band_start }
    }
    pub fn set_band_start(&mut self, v: i32) -> &mut Self {
        unsafe {
            (*self.handle).band_start = v;
        }
        self
    }

    /// Seed emission order (`--seed-order`). [`SeedOrder::Off`] (the default)
    /// is byte-identical to the baseline.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnrecognizedEnum`](crate::Error::UnrecognizedEnum) if
    /// the vendored bwa-mem3 reports a mode this crate does not know, rather
    /// than silently reporting the default.
    pub fn seed_order(&self) -> Result<SeedOrder> {
        SeedOrder::try_from(unsafe { (*self.handle).seed_emit_order })
    }
    pub fn set_seed_order(&mut self, order: SeedOrder) -> &mut Self {
        unsafe {
            (*self.handle).seed_emit_order = order as i32;
        }
        self
    }

    // ---------- bisulfite (--meth) knobs ----------

    /// Enable or disable bisulfite (BS-seq) alignment (`--meth`).
    ///
    /// When enabled, align with a bisulfite dual index
    /// ([`BwaIndex::load_meth`](crate::BwaIndex::load_meth)): reads are
    /// projected to match the converted seed index, alignment runs in original
    /// coordinates, and each record carries Bismark `XR:Z` / `XG:Z` / `XM:Z`
    /// tags. Toggling this rebuilds the per-hypothesis scoring matrices,
    /// unless the mismatch penalty is outside `0..=128`; the next storable
    /// penalty then rebuilds them.
    ///
    /// `set_meth(true)` only flips a flag; it **must** be paired with a meth
    /// dual index ([`BwaIndex::load_meth`](crate::BwaIndex::load_meth)).
    /// Aligning with a mismatch (meth opts + plain index, or vice versa) is
    /// rejected up front with [`Error::InvalidInput`](crate::Error), so the two
    /// cannot silently disagree.
    pub fn set_meth(&mut self, v: bool) -> &mut Self {
        unsafe {
            (*self.handle).meth_mode = i32::from(v);
            sys::bwa_shim_opts_fill_meth_mat(self.handle);
        }
        self
    }
    #[must_use]
    pub fn meth(&self) -> bool {
        unsafe { (*self.handle).meth_mode != 0 }
    }

    /// Bisulfite scoring mode (`--meth-scoring`); only meaningful under
    /// [`set_meth(true)`](Self::set_meth). Rebuilds the scoring matrices,
    /// unless the mismatch penalty is outside `0..=128`; the next storable
    /// penalty then rebuilds them under this mode.
    /// Defaults to [`MethScoring::Genomic`] (bwa-mem3 v0.17.0+;
    /// [`MethScoring::Collapsed`] before), so pass `Collapsed` to keep the
    /// pre-v0.17.0 `--meth` output.
    ///
    /// An explicit mode wins over the [`MethChem::Taps`] default that
    /// [`apply_meth_defaults`](Self::apply_meth_defaults) would otherwise pick.
    pub fn set_meth_scoring(&mut self, m: MethScoring) -> &mut Self {
        unsafe {
            (*self.handle).meth_scoring = m as i32;
            sys::bwa_shim_opts_fill_meth_mat(self.handle);
        }
        self.meth_explicit |= sys::BWA_SHIM_METH_SET_SCORING;
        self
    }
    /// Bisulfite scoring mode (`--meth-scoring`); only meaningful under
    /// [`set_meth(true)`](Self::set_meth). Defaults to
    /// [`MethScoring::Genomic`] (bwa-mem3 v0.17.0+).
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnrecognizedEnum`](crate::Error::UnrecognizedEnum) if
    /// the vendored bwa-mem3 reports a mode this crate does not know.
    pub fn meth_scoring(&self) -> Result<MethScoring> {
        MethScoring::try_from(unsafe { (*self.handle).meth_scoring })
    }

    /// Methylation chemistry (`--meth=emseq|taps`); only meaningful under
    /// [`set_meth(true)`](Self::set_meth). Defaults to [`MethChem::EmSeq`].
    ///
    /// Set it **before** [`apply_meth_defaults`](Self::apply_meth_defaults):
    /// TAPS changes the scoring default that call picks.
    pub fn set_meth_chem(&mut self, chem: MethChem) -> &mut Self {
        unsafe {
            (*self.handle).meth_chem = chem as i32;
        }
        self
    }
    /// Methylation chemistry (`--meth=emseq|taps`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnrecognizedEnum`](crate::Error::UnrecognizedEnum) if
    /// the vendored bwa-mem3 reports a chemistry this crate does not know.
    pub fn meth_chem(&self) -> Result<MethChem> {
        MethChem::try_from(unsafe { (*self.handle).meth_chem })
    }

    /// Where a methylation conversion is free (`--meth-context`); only
    /// meaningful under [`set_meth(true)`](Self::set_meth). Only stored here:
    /// a [`MethContext::Cpg`] that `bwa-mem3 mem` refuses (without `--meth`,
    /// with [`MethChem::EmSeq`], or with [`MethScoring::Collapsed`]) is
    /// rejected by [`apply_meth_defaults`](Self::apply_meth_defaults) and again
    /// by every alignment call, as [`Error::InvalidOpts`](crate::Error::InvalidOpts).
    ///
    /// Unlike the CLI, which refuses any explicit `--meth-context` without
    /// `--meth`, an explicit [`MethContext::All`] without `--meth` is accepted:
    /// it is the default, so it cannot be told apart from an unset context,
    /// and it has no effect.
    pub fn set_meth_context(&mut self, ctx: MethContext) -> &mut Self {
        unsafe {
            (*self.handle).meth_context = ctx as i32;
        }
        self
    }
    /// Where a methylation conversion is free (`--meth-context`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnrecognizedEnum`](crate::Error::UnrecognizedEnum) if
    /// the vendored bwa-mem3 reports a context this crate does not know.
    pub fn meth_context(&self) -> Result<MethContext> {
        MethContext::try_from(unsafe { (*self.handle).meth_context })
    }

    /// Every option check `bwa-mem3 mem` makes before aligning, judged on the
    /// options as they stand ([`check_effective_scoring`](Self::check_effective_scoring)
    /// and [`check_meth_options`](Self::check_meth_options)). Every alignment
    /// entry point runs it, through `check_alignment_preconditions`, so a
    /// combination the CLI refuses cannot reach the aligner however the
    /// setters were ordered.
    pub(crate) fn validate_for_alignment(&self) -> Result<()> {
        self.check_effective_scoring()?;
        self.check_meth_options()
    }

    /// Reject `--meth` options that `bwa-mem3 mem` refuses, judged on the
    /// options as they stand: a non-positive effective mismatch penalty under
    /// `--meth` (`fastmap.cpp`'s post-defaults check: NM is derived from the
    /// scoring matrix, so `-B 0` collapses it to zero), and a
    /// [`MethContext::Cpg`] without `--meth`, with [`MethChem::EmSeq`], or with
    /// [`MethScoring::Collapsed`]. The context checks are the shim's
    /// (`bwa_shim_opts_meth_context_error`), the same copy
    /// [`apply_meth_defaults`](Self::apply_meth_defaults) runs against the
    /// scoring it is about to resolve; repeating them here catches a later
    /// setter, or a skipped `apply_meth_defaults`.
    pub(crate) fn check_meth_options(&self) -> Result<()> {
        if self.meth() {
            let b = self.mismatch_penalty();
            if b <= 0 {
                return Err(crate::Error::InvalidOpts(format!(
                    "--meth requires a positive mismatch penalty, but the effective penalty is \
                     {b}; a non-positive penalty makes every substitution free or rewarded, \
                     which collapses NM to zero and hides real variants (check -B, \
                     MemOpts::set_mismatch_penalty)"
                )));
            }
        }
        // SAFETY: `handle` is a live `mem_opt_t` owned by `self`. The result is
        // NULL or a NUL-terminated string valid until the next shim call on
        // this thread, copied out before any other.
        let msg = unsafe { sys::bwa_shim_opts_meth_context_error(self.handle) };
        if msg.is_null() {
            return Ok(());
        }
        let msg = unsafe { std::ffi::CStr::from_ptr(msg) };
        Err(crate::Error::InvalidOpts(
            msg.to_string_lossy().into_owned(),
        ))
    }

    /// Reject options `bwa-mem3 mem` refuses on every path, judged on the
    /// options as they stand. First the ranges its parser holds for the knobs
    /// this crate can set (`-A` in `1..=127`, `-w` in `0..=i32::MAX >> 2`,
    /// `-O` non-negative, `-E` positive, `-U` non-negative, `-c` positive),
    /// since the setters accept any `i32`; then its `check_effective_scoring`
    /// (`fastmap.cpp`): a mismatch penalty outside `0..=128`, since the kernels
    /// store its negation in 8 bits (256 would wrap to a free mismatch) and a
    /// negative one breaks the 8-bit admission tests and the MAPQ model; and a
    /// gap open plus extension above `2^29 - 1`, which overflows the banded
    /// kernels' gap arithmetic. Every alignment entry point runs this.
    ///
    /// The parse ranges are memory safety, not CLI parity: outside them the
    /// kernels overflow signed `int` arithmetic (`-U`: `score -
    /// pen_unpaired`), divide by zero (`-c` is the SA-resolve step divisor),
    /// or write past the rows a band owns (a negative `-w`). The messages are
    /// upstream's, naming the setter alongside the flag. `-L`/`-T` are absent
    /// because upstream bounds neither.
    pub(crate) fn check_effective_scoring(&self) -> Result<()> {
        const OE_MAX: i64 = (1 << 29) - 1;
        let invalid = |msg: String| Err(crate::Error::InvalidOpts(msg));
        let a = self.match_score();
        if !(1..=127).contains(&a) {
            return invalid(format!(
                "-A match score must be an integer in 1..127 (got {a}): a score of 0 or less \
                 leaves nothing to align, and the alignment kernels store it in 8 bits \
                 (MemOpts::set_match_score)"
            ));
        }
        // SAFETY: `handle` is a live `mem_opt_t` owned by `self`; reads only.
        let (w, o_del, e_del, o_ins, e_ins, pen_unpaired, max_occ) = unsafe {
            let h = &*self.handle;
            (
                h.w,
                h.o_del,
                h.e_del,
                h.o_ins,
                h.e_ins,
                h.pen_unpaired,
                h.max_occ,
            )
        };
        if !(0..=i32::MAX >> 2).contains(&w) {
            return invalid(format!(
                "-w band width must be an integer in 0..{} (got {w}) (MemOpts::set_band_width)",
                i32::MAX >> 2
            ));
        }
        if o_del < 0 || o_ins < 0 {
            return invalid(format!(
                "-O gap-open penalty must be a non-negative integer in 0..{} (got \
                 {o_del},{o_ins}) (MemOpts::set_gap_open)",
                i32::MAX
            ));
        }
        if e_del < 1 || e_ins < 1 {
            return invalid(format!(
                "-E gap-extension penalty must be a positive integer in 1..{} (got \
                 {e_del},{e_ins}) (MemOpts::set_gap_extend)",
                i32::MAX
            ));
        }
        if pen_unpaired < 0 {
            return invalid(format!(
                "-U unpaired penalty must be a non-negative integer in 0..{} (got \
                 {pen_unpaired}) (MemOpts::set_unpaired_penalty)",
                i32::MAX
            ));
        }
        if max_occ < 1 {
            return invalid(format!(
                "-c max occurrences must be a positive integer in 1..{} (got {max_occ}) \
                 (MemOpts::set_max_occurrences)",
                i32::MAX
            ));
        }
        let b = self.mismatch_penalty();
        if !(0..=128).contains(&b) {
            return Err(crate::Error::InvalidOpts(format!(
                "the effective mismatch penalty is {b}, outside 0..128: a negative penalty \
                 scores a mismatch above zero, which the 8-bit alignment kernels and the \
                 mapping-quality model do not allow, and the kernels store its negation in 8 \
                 bits (check -B, MemOpts::set_mismatch_penalty; unlike `bwa-mem3 mem -A`, \
                 MemOpts::set_match_score does not rescale it)"
            )));
        }
        let (del, ins) = (
            i64::from(o_del) + i64::from(e_del),
            i64::from(o_ins) + i64::from(e_ins),
        );
        if del > OE_MAX || ins > OE_MAX {
            return Err(crate::Error::InvalidOpts(format!(
                "the effective gap open plus extension penalty is {del} for deletions and \
                 {ins} for insertions; each must be at most {OE_MAX} so the alignment \
                 kernels' gap scores stay in range (check -O and -E, MemOpts::set_gap_open \
                 and MemOpts::set_gap_extend)"
            )));
        }
        Ok(())
    }

    /// The options as they bear on alignment, for checking that every phase
    /// of a split-phase alignment (seeding, extension, the cohort pestat,
    /// pairing and emission) sees the same ones: the raw `mem_opt_t` bytes
    /// with the fields that cannot change the result blanked. Those are the
    /// emission-only knobs (`md_mode` and the `MEM_F_XB` bit, which only shape
    /// the tags written, and `tie_break_mode`, which only the emit phase reads,
    /// through `mem_tie_break_id`), and the two the shim overrides on its own per-call
    /// copy in every phase (`MEM_F_PE`, set for the pair group and cleared for
    /// singles, and `n_threads`, forced to 1). Every other field is compared,
    /// so a field a future bwa-mem3 adds is covered by default.
    pub(crate) fn phase_snapshot(&self) -> PhaseOpts {
        let n = std::mem::size_of::<sys::mem_opt_t>();
        let mut bytes = vec![0u8; n];
        // The byte range of the field at `p`, from its own type, so a field
        // upstream widens is still blanked whole.
        fn span<T>(base: *const sys::mem_opt_t, p: *const T) -> std::ops::Range<usize> {
            let off = p as usize - base as usize;
            off..off + std::mem::size_of::<T>()
        }
        // SAFETY: `handle` is a live `mem_opt_t` owned by `self`, allocated
        // zeroed by `mem_opt_init` (calloc), so its padding bytes are
        // initialized too; `bytes` holds exactly `n` bytes. The field pointers
        // are taken without creating references.
        let (md, tie, threads, flag) = unsafe {
            std::ptr::copy_nonoverlapping(self.handle.cast::<u8>(), bytes.as_mut_ptr(), n);
            let h = self.handle.cast_const();
            // `flag` is read back as an `i32` to clear two bits; a retyped
            // field must fail to compile here rather than be misread.
            let flag: *const i32 = std::ptr::addr_of!((*h).flag);
            (
                span(h, std::ptr::addr_of!((*h).md_mode)),
                span(h, std::ptr::addr_of!((*h).tie_break_mode)),
                span(h, std::ptr::addr_of!((*h).n_threads)),
                span(h, flag),
            )
        };
        bytes[md].fill(0);
        bytes[tie].fill(0);
        bytes[threads].fill(0);
        let flag_bytes: [u8; 4] = bytes[flag.clone()].try_into().unwrap();
        let flag_val = i32::from_ne_bytes(flag_bytes) & !((sys::MEM_F_XB | sys::MEM_F_PE) as i32);
        bytes[flag].copy_from_slice(&flag_val.to_ne_bytes());
        PhaseOpts(bytes.into_boxed_slice())
    }

    /// Reject options that differ from `seeded`, the snapshot taken when the
    /// regions `what` consumes were seeded: those regions were scored under
    /// the old options, and pairing, rescue and `mem_reg2aln` under the new
    /// ones would mix the two in one record.
    pub(crate) fn check_phase(&self, seeded: &PhaseOpts, what: &str) -> Result<()> {
        self.phase_snapshot().check_against(seeded, what)
    }

    /// How the `MD:Z` tag is written (`--md`). Defaults to [`MdMode::Full`].
    pub fn set_md_mode(&mut self, mode: MdMode) -> &mut Self {
        unsafe {
            (*self.handle).md_mode = mode as i32;
        }
        self
    }
    /// How the `MD:Z` tag is written (`--md`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnrecognizedEnum`](crate::Error::UnrecognizedEnum) if
    /// the vendored bwa-mem3 reports a mode this crate does not know.
    pub fn md_mode(&self) -> Result<MdMode> {
        MdMode::try_from(unsafe { (*self.handle).md_mode })
    }

    /// Bisulfite seed pruning (`--meth-seed-prune`); only meaningful under
    /// [`set_meth(true)`](Self::set_meth). An explicit choice survives
    /// [`apply_meth_defaults`](Self::apply_meth_defaults), whichever is
    /// called first.
    pub fn set_meth_seed_prune(&mut self, prune: MethSeedPrune) -> &mut Self {
        unsafe {
            (*self.handle).meth_seed_prune = prune as i32;
        }
        self.meth_explicit |= sys::BWA_SHIM_METH_SET_SEED_PRUNE;
        self
    }
    /// Bisulfite seed pruning (`--meth-seed-prune`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnrecognizedEnum`](crate::Error::UnrecognizedEnum) if
    /// the vendored bwa-mem3 reports a mode this crate does not know.
    pub fn meth_seed_prune(&self) -> Result<MethSeedPrune> {
        MethSeedPrune::try_from(unsafe { (*self.handle).meth_seed_prune })
    }

    /// bwameth.py-style longest-`M` chimera QC heuristic (`--meth-chimera-qc`);
    /// only meaningful under [`set_meth(true)`](Self::set_meth). Off by default.
    pub fn set_meth_chimera_qc(&mut self, v: bool) -> &mut Self {
        unsafe {
            (*self.handle).meth_chimera_qc = i32::from(v);
        }
        self
    }
    #[must_use]
    pub fn meth_chimera_qc(&self) -> bool {
        unsafe { (*self.handle).meth_chimera_qc != 0 }
    }

    pub(crate) fn as_ptr(&self) -> *const sys::mem_opt_t {
        self.handle
    }
}

impl Drop for MemOpts {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe { sys::bwa_shim_opts_free(self.handle) };
        }
    }
}

// SAFETY: mem_opt_t is POD. MemOpts exclusively owns the heap allocation;
// no aliasing. &MemOpts is passed as `const` across FFI.
unsafe impl Send for MemOpts {}
unsafe impl Sync for MemOpts {}

/// Paired-end orientation key.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PeOrientation {
    Ff = 0,
    Fr = 1,
    Rf = 2,
    Rr = 3,
}

/// Observed insert-size statistics for a single orientation.
#[derive(Debug, Clone, Copy, Default)]
pub struct PeOrient {
    pub low: i32,
    pub high: i32,
    pub failed: bool,
    pub avg: f64,
    pub std: f64,
}

impl From<sys::mem_pestat_t> for PeOrient {
    fn from(p: sys::mem_pestat_t) -> Self {
        Self {
            low: p.low,
            high: p.high,
            failed: p.failed != 0,
            avg: p.avg,
            std: p.std,
        }
    }
}

impl From<PeOrient> for sys::mem_pestat_t {
    fn from(p: PeOrient) -> Self {
        sys::mem_pestat_t {
            low: p.low,
            high: p.high,
            failed: if p.failed { 1 } else { 0 },
            avg: p.avg,
            std: p.std,
        }
    }
}

/// Insert-size model across all four orientations (`mem_pestat_t[4]`).
pub struct MemPeStat {
    pub(crate) handle: *mut sys::mem_pestat_t,
}

impl MemPeStat {
    /// Construct a zeroed insert-size model.
    pub fn zero() -> Result<Self> {
        let handle = unsafe { sys::bwa_shim_pestat_zero() };
        if handle.is_null() {
            return Err(shim_err("pestat_zero"));
        }
        Ok(MemPeStat { handle })
    }

    #[must_use]
    pub fn orientation(&self, o: PeOrientation) -> PeOrient {
        unsafe { (*self.handle.add(o as usize)).into() }
    }

    pub fn set_orientation(&mut self, o: PeOrientation, v: PeOrient) -> &mut Self {
        unsafe {
            *self.handle.add(o as usize) = v.into();
        }
        self
    }

    pub(crate) fn as_ptr(&self) -> *const sys::mem_pestat_t {
        self.handle
    }

    pub(crate) fn as_mut_ptr(&mut self) -> *mut sys::mem_pestat_t {
        self.handle
    }
}

impl Drop for MemPeStat {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe { sys::bwa_shim_pestat_free(self.handle) };
        }
    }
}

unsafe impl Send for MemPeStat {}
unsafe impl Sync for MemPeStat {}

#[cfg(test)]
mod tests {
    use super::*;
    // Only the tests match on the variant; a module-level import would be
    // dead in a non-test build and trip clippy's -D unused-imports.
    use crate::error::Error;

    /// Reading `mat` requires reaching through the raw handle: it is not part
    /// of the safe surface, and these tests are the only consumer.
    fn scoring_matrix(o: &MemOpts) -> [i8; 25] {
        unsafe { (*o.handle).mat }
    }

    #[test]
    fn changing_the_mismatch_penalty_rebuilds_the_scoring_matrix() {
        let mut o = MemOpts::new().unwrap();
        // mem_opt_init defaults: a=1, b=4, so off-diagonal cells are -4.
        assert_eq!(
            scoring_matrix(&o)[1],
            -4,
            "unexpected default mismatch cell"
        );

        o.set_mismatch_penalty(2);

        assert_eq!(o.mismatch_penalty(), 2);
        assert_eq!(
            scoring_matrix(&o)[1],
            -2,
            "set_mismatch_penalty must rebuild opt->mat, not just write opt->b -- \
             the SW kernels read only the matrix"
        );
    }

    #[test]
    fn changing_the_match_score_rebuilds_the_scoring_matrix() {
        let mut o = MemOpts::new().unwrap();
        assert_eq!(scoring_matrix(&o)[0], 1, "unexpected default match cell");

        o.set_match_score(3);

        assert_eq!(o.match_score(), 3);
        assert_eq!(
            scoring_matrix(&o)[0],
            3,
            "set_match_score must rebuild opt->mat"
        );
    }

    #[test]
    fn rebuilding_the_scoring_matrix_refreshes_the_meth_matrices() {
        let mut o = MemOpts::new().unwrap();
        o.set_meth(true);
        o.set_mismatch_penalty(2);

        // mat_ot is a copy of mat with the C->T cell freed to a match, so its
        // ordinary off-diagonal cells must track the new penalty. A stale copy
        // would silently keep -4 and ignore the caller's scoring.
        let mat_ot = unsafe { (*o.handle).mat_ot };
        assert_eq!(
            mat_ot[1], -2,
            "the bisulfite matrices are copies of opt->mat and must be refilled \
             whenever it is rebuilt (mem_opt_fill_meth_mat's own contract)"
        );
    }

    #[test]
    fn mark_split_secondary_round_trips() {
        let mut o = MemOpts::new().unwrap();
        assert!(!o.mark_split_secondary(), "-M is off by default");
        o.set_mark_split_secondary(true);
        assert!(o.mark_split_secondary());
        o.set_mark_split_secondary(false);
        assert!(!o.mark_split_secondary());
    }

    #[test]
    fn construct_defaults() {
        let o = MemOpts::new().unwrap();
        // bwa-mem3 defaults (from mem_opt_init):
        // a=1, b=4, o_del=o_ins=6, e_del=e_ins=1, w=100, T=30, pen_clip5=pen_clip3=5
        assert_eq!(o.match_score(), 1);
        assert_eq!(o.mismatch_penalty(), 4);
        assert_eq!(o.band_width(), 100);
        assert_eq!(o.minimum_score(), 30);
    }

    #[test]
    fn mode_pacbio_changes_defaults() {
        let o = MemOpts::new().unwrap().with_mode(Mode::Pacbio);
        assert_eq!(o.min_seed_len(), 17);
        assert_eq!(o.mismatch_penalty(), 1);
        assert_eq!(o.minimum_score(), 0);
    }

    #[test]
    fn setters_chain() {
        let mut o = MemOpts::new().unwrap();
        o.set_min_seed_len(19).set_band_width(50).set_match_score(2);
        assert_eq!(o.min_seed_len(), 19);
        assert_eq!(o.band_width(), 50);
        assert_eq!(o.match_score(), 2);
    }

    #[test]
    fn new_opt_in_knobs_default_off() {
        let o = MemOpts::new().unwrap();
        assert_eq!(o.min_ext_len(), 0);
        assert_eq!(o.max_extend_chains(), 0);
        assert_eq!(o.mate_concordant_window(), 0);
        assert!(!o.smem_dedup());
        assert_eq!(o.band_start(), 0);
        assert_eq!(o.seed_order().unwrap(), SeedOrder::Off);
    }

    #[test]
    fn skip_contained_ext_defaults_on() {
        // Upstream (fg-labs/bwa-mem3 since the contained-seed extension skip
        // became a byte-identical default) initializes skip_contained_ext = 1:
        // it is a pure speed lever with output identical to the reference path,
        // so it is on by default rather than opt-in. `--keep-contained-ext`
        // (set_skip_contained_ext(false)) opts out.
        let o = MemOpts::new().unwrap();
        assert!(o.skip_contained_ext());
    }

    #[test]
    fn new_opt_in_knobs_roundtrip() {
        let mut o = MemOpts::new().unwrap();
        o.set_min_ext_len(25)
            .set_max_extend_chains(4)
            .set_mate_concordant_window(-1)
            .set_smem_dedup(true)
            // default is on (see skip_contained_ext_defaults_on), so toggle it
            // OFF here to actually exercise the setter (the --keep-contained-ext
            // opt-out direction).
            .set_skip_contained_ext(false)
            .set_band_start(20)
            .set_seed_order(SeedOrder::MostAbsorb);
        assert_eq!(o.min_ext_len(), 25);
        assert_eq!(o.max_extend_chains(), 4);
        assert_eq!(o.mate_concordant_window(), -1);
        assert!(o.smem_dedup());
        assert!(!o.skip_contained_ext());
        assert_eq!(o.band_start(), 20);
        assert_eq!(o.seed_order().unwrap(), SeedOrder::MostAbsorb);
    }

    #[test]
    fn seed_order_all_variants_roundtrip() {
        let mut o = MemOpts::new().unwrap();
        for order in [
            SeedOrder::Off,
            SeedOrder::GlobalLongest,
            SeedOrder::LocalLongest,
            SeedOrder::AbsorbCount,
            SeedOrder::MostAbsorb,
        ] {
            o.set_seed_order(order);
            assert_eq!(o.seed_order().unwrap(), order);
        }
    }

    #[test]
    fn meth_knobs_roundtrip() {
        let mut o = MemOpts::new().unwrap();
        // Defaults: meth off, genomic scoring (bwa-mem3 v0.17.0, #756).
        assert!(!o.meth());
        assert_eq!(o.meth_scoring().unwrap(), MethScoring::Genomic);
        assert!(!o.meth_chimera_qc());
        o.set_meth(true)
            .set_meth_scoring(MethScoring::Collapsed)
            .set_meth_chimera_qc(true);
        assert!(o.meth());
        assert_eq!(o.meth_scoring().unwrap(), MethScoring::Collapsed);
        assert!(o.meth_chimera_qc());
        o.set_meth_scoring(MethScoring::Genomic);
        assert_eq!(o.meth_scoring().unwrap(), MethScoring::Genomic);
        o.set_meth(false);
        assert!(!o.meth());
    }

    #[test]
    fn pe_flag_toggles() {
        let mut o = MemOpts::new().unwrap();
        o.set_pe(true);
        o.set_pe(false);
        o.set_pe(true);
    }

    #[test]
    fn pestat_round_trip() {
        let mut p = MemPeStat::zero().unwrap();
        p.set_orientation(
            PeOrientation::Fr,
            PeOrient {
                low: 100,
                high: 500,
                failed: false,
                avg: 250.0,
                std: 50.0,
            },
        );
        let got = p.orientation(PeOrientation::Fr);
        assert_eq!(got.low, 100);
        assert_eq!(got.high, 500);
        assert!((got.avg - 250.0).abs() < 1e-9);
        assert!(!got.failed);
    }

    #[test]
    fn set_read_group_id_roundtrip() {
        let mut o = MemOpts::new().unwrap();
        o.set_read_group_id(Some("sample42")).unwrap();
        o.set_read_group_id(None).unwrap();
        // No getter yet (the id is a process-wide global in bwa-mem3); this test
        // just verifies neither set nor clear panics.
    }

    #[test]
    fn pestat_orientations_independent() {
        let mut p = MemPeStat::zero().unwrap();
        p.set_orientation(
            PeOrientation::Ff,
            PeOrient {
                low: 1,
                high: 2,
                ..Default::default()
            },
        );
        assert_eq!(p.orientation(PeOrientation::Fr).low, 0);
    }

    #[test]
    fn known_seed_order_values_convert() {
        assert_eq!(SeedOrder::try_from(0).unwrap(), SeedOrder::Off);
        assert_eq!(SeedOrder::try_from(4).unwrap(), SeedOrder::MostAbsorb);
    }

    #[test]
    fn unknown_seed_order_value_is_an_error_not_off() {
        // Upstream adding a 6th mode must not silently read back as Off.
        //
        // Matched on the variant, not on the rendered string: `Error` has six
        // variants, so a regression returning `InvalidOpts("bad value 5")` —
        // or an `IndexLoad` whose path merely contains the character '5' —
        // would satisfy a `to_string().contains('5')` assertion. `kind` is
        // asserted too, so a copy-paste swap with `meth_scoring` fails here.
        let err = SeedOrder::try_from(5).unwrap_err();
        assert!(
            matches!(
                err,
                Error::UnrecognizedEnum {
                    kind: "seed_emit_order",
                    value: 5
                }
            ),
            "expected UnrecognizedEnum{{ kind: \"seed_emit_order\", value: 5 }}, got: {err:?}"
        );
    }

    #[test]
    fn known_meth_scoring_values_convert() {
        assert_eq!(MethScoring::try_from(0).unwrap(), MethScoring::Collapsed);
        assert_eq!(MethScoring::try_from(1).unwrap(), MethScoring::Genomic);
    }

    #[test]
    fn neutral_meth_scoring_maps_to_its_own_variant() {
        // bwa-mem3 v0.9.0 uses MEM_METH_SCORING_NEUTRAL = 2 as the automatic
        // default for TAPS chemistry, so it is a value the vendored library
        // really produces. Before it was mapped it read back as Collapsed,
        // which silently mis-described the scoring in use.
        assert_eq!(MethScoring::try_from(2).unwrap(), MethScoring::Neutral);
    }

    #[test]
    fn unknown_meth_scoring_value_is_an_error_not_a_fallback() {
        // The guarantee is that an unrecognized value ERRORS rather than
        // collapsing to a default -- checked at a value bwa-mem3 does not
        // define, so this keeps testing the guarantee as upstream adds modes.
        let err = MethScoring::try_from(3).unwrap_err();
        assert!(
            matches!(
                err,
                Error::UnrecognizedEnum {
                    kind: "meth_scoring",
                    value: 3
                }
            ),
            "expected UnrecognizedEnum{{ kind: \"meth_scoring\", value: 3 }}, got: {err:?}"
        );
    }

    #[test]
    fn round_trip_through_the_raw_discriminant() {
        for v in [0, 1, 2, 3, 4] {
            let order = SeedOrder::try_from(v).unwrap();
            assert_eq!(order as i32, v, "discriminant must round-trip");
        }
    }

    #[test]
    fn meth_scoring_round_trips_through_the_raw_discriminant() {
        // The setters depend on `m as i32` round-tripping through
        // TryFrom<i32> back to the same variant — mirrors SeedOrder's
        // round-trip test above, which caught the same class of bug there.
        for v in [0, 1] {
            let scoring = MethScoring::try_from(v).unwrap();
            assert_eq!(scoring as i32, v, "discriminant must round-trip");
        }
    }

    #[test]
    fn meth_chem_and_seed_prune_round_trip_through_the_raw_discriminant() {
        for v in [0, 1] {
            assert_eq!(MethChem::try_from(v).unwrap() as i32, v);
        }
        for v in [0, 1, 2] {
            assert_eq!(MethSeedPrune::try_from(v).unwrap() as i32, v);
        }
    }

    #[test]
    fn unknown_meth_chem_and_seed_prune_values_are_errors_not_fallbacks() {
        assert!(matches!(
            MethChem::try_from(2).unwrap_err(),
            Error::UnrecognizedEnum {
                kind: "meth_chem",
                value: 2
            }
        ));
        assert!(matches!(
            MethSeedPrune::try_from(3).unwrap_err(),
            Error::UnrecognizedEnum {
                kind: "meth_seed_prune",
                value: 3
            }
        ));
    }

    /// The raw `mem_opt_t` bytes, for "left exactly as they were" checks: a
    /// refused call must not change any field, not just the few a test names.
    fn raw_bytes(o: &MemOpts) -> Vec<u8> {
        let n = std::mem::size_of::<sys::mem_opt_t>();
        // SAFETY: `handle` is a live, zero-initialized `mem_opt_t` of `n` bytes.
        unsafe { std::slice::from_raw_parts(o.handle.cast::<u8>(), n).to_vec() }
    }

    /// Every scoring matrix the options carry: the base one, the per-strand
    /// `--meth` copies, and the `--meth-context cpg` 7x7 pair
    /// (`meth_model_fill_mats`), all rebuilt together by one refill.
    fn all_matrices(o: &MemOpts) -> Vec<i8> {
        // SAFETY: `handle` is a live `mem_opt_t` owned by `o`; reads only.
        let h = unsafe { &*o.handle };
        [
            &h.mat[..],
            &h.mat_ot[..],
            &h.mat_ob[..],
            &h.mat7_ot[..],
            &h.mat7_ob[..],
        ]
        .concat()
    }

    /// Assert that `r` is `InvalidOpts` with a message containing `needle`.
    fn assert_invalid<T>(r: Result<T>, needle: &str) {
        match r {
            Err(crate::Error::InvalidOpts(m)) => assert!(m.contains(needle), "{m}"),
            Err(e) => panic!("expected InvalidOpts({needle}), got {e}"),
            Ok(_) => panic!("expected InvalidOpts({needle}), got Ok"),
        }
    }

    /// Fresh options with `--meth` on and the chemistry set, as a caller sets
    /// them up before `apply_meth_defaults`.
    fn meth_opts(chem: MethChem) -> MemOpts {
        let mut o = MemOpts::new().unwrap();
        o.set_meth(true).set_meth_chem(chem);
        o
    }

    #[test]
    fn meth_defaults_follow_the_cli_for_emseq() {
        let mut o = meth_opts(MethChem::EmSeq);
        let default_b = o.mismatch_penalty();
        assert_eq!(
            o.meth_seed_prune().unwrap(),
            MethSeedPrune::Off,
            "mem_opt_init default"
        );
        o.apply_meth_defaults().unwrap();
        assert_eq!(o.meth_chem().unwrap(), MethChem::EmSeq);
        assert_eq!(
            o.meth_scoring().unwrap(),
            MethScoring::Genomic,
            "em-seq scores GENOMIC by default as of bwa-mem3 v0.17.0"
        );
        assert_eq!(o.meth_seed_prune().unwrap(), MethSeedPrune::Spec30);
        assert_eq!(
            o.mismatch_penalty(),
            default_b,
            "GENOMIC is variant-aware: no -B 2"
        );
    }

    #[test]
    fn explicit_collapsed_scoring_takes_bwameths_mismatch_penalty() {
        let mut o = meth_opts(MethChem::EmSeq);
        let default_b = o.mismatch_penalty();
        o.set_meth_scoring(MethScoring::Collapsed);
        o.apply_meth_defaults().unwrap();
        assert_eq!(o.meth_scoring().unwrap(), MethScoring::Collapsed);
        assert_eq!(o.mismatch_penalty(), 2, "COLLAPSED takes bwameth's -B 2");
        assert_ne!(default_b, 2, "the -B check above must not hold vacuously");
    }

    #[test]
    fn taps_defaults_to_neutral_scoring_and_keeps_bwas_mismatch_penalty() {
        let mut o = meth_opts(MethChem::Taps);
        let default_b = o.mismatch_penalty();
        o.apply_meth_defaults().unwrap();
        assert_eq!(o.meth_scoring().unwrap(), MethScoring::Neutral);
        assert_eq!(o.meth_seed_prune().unwrap(), MethSeedPrune::Spec30);
        assert_eq!(
            o.mismatch_penalty(),
            default_b,
            "NEUTRAL is variant-aware: no -B 2"
        );
    }

    #[test]
    fn explicit_meth_scoring_wins_over_the_taps_default_in_either_order() {
        let mut before = MemOpts::new().unwrap();
        before
            .set_meth(true)
            .set_meth_scoring(MethScoring::Genomic)
            .set_meth_chem(MethChem::Taps)
            .apply_meth_defaults()
            .unwrap();
        assert_eq!(before.meth_scoring().unwrap(), MethScoring::Genomic);

        let mut after = meth_opts(MethChem::Taps);
        after
            .set_meth_scoring(MethScoring::Collapsed)
            .apply_meth_defaults()
            .unwrap();
        assert_eq!(
            after.meth_scoring().unwrap(),
            MethScoring::Collapsed,
            "an explicit COLLAPSED is not the same as unset"
        );
    }

    #[test]
    fn md_mode_roundtrips_and_defaults_to_full() {
        let mut o = MemOpts::new().unwrap();
        assert_eq!(o.md_mode().unwrap(), MdMode::Full, "mem_opt_init default");
        for mode in [MdMode::None, MdMode::Short, MdMode::Full] {
            o.set_md_mode(mode);
            assert_eq!(o.md_mode().unwrap(), mode);
        }
    }

    #[test]
    fn meth_context_roundtrips_and_defaults_to_all() {
        let mut o = MemOpts::new().unwrap();
        assert_eq!(o.meth_context().unwrap(), MethContext::All);
        o.set_meth_context(MethContext::Cpg);
        assert_eq!(o.meth_context().unwrap(), MethContext::Cpg);
    }

    #[test]
    fn meth_context_and_md_mode_map_raw_values() {
        assert_eq!(MethContext::try_from(0).unwrap(), MethContext::All);
        assert_eq!(MethContext::try_from(1).unwrap(), MethContext::Cpg);
        assert!(matches!(
            MethContext::try_from(2),
            Err(crate::Error::UnrecognizedEnum {
                kind: "meth_context",
                value: 2
            })
        ));
        assert_eq!(MdMode::try_from(0).unwrap(), MdMode::None);
        assert_eq!(MdMode::try_from(1).unwrap(), MdMode::Full);
        assert_eq!(MdMode::try_from(2).unwrap(), MdMode::Short);
        assert!(matches!(
            MdMode::try_from(3),
            Err(crate::Error::UnrecognizedEnum {
                kind: "md_mode",
                value: 3
            })
        ));
    }

    /// The COLLAPSED check judges the scoring as `meth_opts_finalize` would
    /// resolve it: under TAPS an unset scoring becomes NEUTRAL, so a COLLAPSED
    /// value the caller did not choose is not a conflict, while an explicit one
    /// is.
    #[test]
    fn cpg_collapsed_check_uses_the_resolved_taps_scoring() {
        let mut implicit = meth_opts(MethChem::Taps);
        implicit.set_meth_context(MethContext::Cpg);
        // A COLLAPSED value that did not come through set_meth_scoring.
        unsafe { (*implicit.handle).meth_scoring = MethScoring::Collapsed as i32 };
        implicit.apply_meth_defaults().unwrap();
        assert_eq!(implicit.meth_scoring().unwrap(), MethScoring::Neutral);

        let mut explicit = meth_opts(MethChem::Taps);
        explicit
            .set_meth_context(MethContext::Cpg)
            .set_meth_scoring(MethScoring::Collapsed);
        assert_invalid(
            explicit.apply_meth_defaults(),
            "incompatible with --meth-scoring collapsed",
        );
    }

    /// The bundle multiplies the match score by up to 100, so under `--meth` a
    /// score outside `-A`'s `1..=127` is refused before it is applied (a large
    /// one would overflow `int`), with every option left as it was. The message
    /// is `check_effective_scoring`'s, word for word: the shim carries its own
    /// copy for raw C callers. Without `--meth` nothing is multiplied and the
    /// entry points judge the score.
    #[test]
    fn apply_meth_defaults_refuses_an_out_of_range_match_score() {
        for a in [0, 128, 21_474_837, i32::MAX] {
            let mut o = meth_opts(MethChem::Taps);
            o.set_match_score(a);
            let before = raw_bytes(&o);
            let Err(crate::Error::InvalidOpts(scoring_msg)) = o.check_effective_scoring() else {
                panic!("a = {a}: check_effective_scoring accepted it");
            };
            match o.apply_meth_defaults() {
                Err(crate::Error::InvalidOpts(msg)) => assert_eq!(msg, scoring_msg, "a = {a}"),
                Err(e) => panic!("a = {a}: expected InvalidOpts, got {e}"),
                Ok(_) => panic!("a = {a}: accepted"),
            }
            assert!(raw_bytes(&o) == before, "a = {a}: options changed");
        }
        let mut ok = meth_opts(MethChem::Taps);
        ok.set_match_score(127);
        ok.apply_meth_defaults().unwrap();
        assert_eq!(ok.unpaired_penalty(), 12_700);

        let mut off = MemOpts::new().unwrap();
        off.set_match_score(i32::MAX);
        off.apply_meth_defaults().unwrap();
    }

    /// The apply-time check can be skipped, or undone by a later setter, so the
    /// alignment entry points re-check the options as they stand.
    #[test]
    fn check_meth_options_judges_the_options_as_they_stand() {
        let mut o = meth_opts(MethChem::Taps);
        o.set_meth_context(MethContext::Cpg);
        o.apply_meth_defaults().unwrap();
        o.check_meth_options().unwrap();

        o.set_meth_scoring(MethScoring::Collapsed);
        assert_invalid(
            o.check_meth_options(),
            "incompatible with --meth-scoring collapsed",
        );
        o.set_meth_scoring(MethScoring::Neutral)
            .set_meth_chem(MethChem::EmSeq);
        assert_invalid(
            o.check_meth_options(),
            "requires --meth=taps or --meth=5base",
        );

        let mut skipped = meth_opts(MethChem::EmSeq);
        skipped.set_meth_context(MethContext::Cpg);
        assert_invalid(
            skipped.check_meth_options(),
            "requires --meth=taps or --meth=5base",
        );

        let mut off = MemOpts::new().unwrap();
        off.set_meth_context(MethContext::Cpg);
        assert_invalid(off.check_meth_options(), "--meth-context requires --meth");
        off.set_meth_context(MethContext::All);
        off.check_meth_options().unwrap();
    }

    /// `bwa-mem3 mem --meth` refuses a non-positive effective mismatch
    /// penalty (NM is derived from the scoring matrix, so `-B 0` collapses it
    /// to zero); without `--meth` a zero penalty is allowed.
    #[test]
    fn check_meth_options_rejects_a_non_positive_penalty_under_meth() {
        for b in [0, -1] {
            let mut o = meth_opts(MethChem::EmSeq);
            o.apply_meth_defaults().unwrap();
            o.set_mismatch_penalty(b);
            let Err(err) = o.check_meth_options() else {
                panic!("expected an error for -B {b}");
            };
            assert!(
                matches!(&err, crate::Error::InvalidOpts(m)
                    if m.contains("--meth requires a positive mismatch penalty")
                        && m.contains(&format!("effective penalty is {b}"))),
                "{err}"
            );
        }

        let mut off = MemOpts::new().unwrap();
        off.set_mismatch_penalty(0);
        off.check_meth_options().unwrap();
    }

    /// `bwa-mem3 mem`'s `check_effective_scoring`, on every path (not only
    /// `--meth`): the kernels store `-b` in 8 bits, so `b` must lie in
    /// `0..=128` (256 would wrap to a free mismatch), and each `o + e` must
    /// leave the banded kernels' gap arithmetic in range.
    #[test]
    fn check_effective_scoring_mirrors_the_cli() {
        let err_has = |o: &MemOpts, needle: &str| match o.check_effective_scoring() {
            Err(crate::Error::InvalidOpts(m)) => assert!(m.contains(needle), "{m}"),
            other => panic!("expected InvalidOpts({needle}), got {other:?}"),
        };
        for meth in [false, true] {
            for b in [-1, 129, 256] {
                let mut o = MemOpts::new().unwrap();
                o.set_meth(meth).set_mismatch_penalty(b);
                err_has(
                    &o,
                    &format!("effective mismatch penalty is {b}, outside 0..128"),
                );
            }
            for b in [0, 1, 128] {
                let mut o = MemOpts::new().unwrap();
                o.set_meth(meth).set_mismatch_penalty(b);
                o.check_effective_scoring().unwrap();
            }
        }
        // Each of deletions and insertions is bounded on its own.
        let oe_max = (1 << 29) - 1;
        let mut o = MemOpts::new().unwrap();
        o.set_gap_open(oe_max - 1, oe_max - 1).set_gap_extend(1, 1);
        o.check_effective_scoring().unwrap();
        o.set_gap_open(oe_max, 6);
        err_has(&o, &format!("penalty is {} for deletions", oe_max + 1));
        o.set_gap_open(6, oe_max);
        err_has(&o, &format!("{} for insertions", oe_max + 1));

        // The signs `bwa-mem3 mem` holds at parse, checked first: a negative
        // open or a non-positive extension, whose sum `i32::MIN + -1` would
        // pass the upper bound and overflow the kernels' `int` arithmetic.
        let mut o = MemOpts::new().unwrap();
        o.set_gap_open(0, 0).set_gap_extend(1, 1);
        o.check_effective_scoring().unwrap();
        o.set_gap_open(i32::MIN, 6);
        err_has(&o, "-O gap-open penalty must be a non-negative integer");
        o.set_gap_open(6, -1);
        err_has(&o, "-O gap-open penalty must be a non-negative integer");
        o.set_gap_open(6, 6).set_gap_extend(-1, 1);
        err_has(&o, "-E gap-extension penalty must be a positive integer");
        o.set_gap_extend(1, 0);
        err_has(&o, "-E gap-extension penalty must be a positive integer");

        // `-A`: 1..=127, as `bwa-mem3 mem` parses it (the matrix stores it in
        // 8 bits).
        for a in [0, -1, 128] {
            let mut o = MemOpts::new().unwrap();
            o.set_match_score(a);
            err_has(
                &o,
                &format!("-A match score must be an integer in 1..127 (got {a})"),
            );
        }
        let mut o = MemOpts::new().unwrap();
        o.set_match_score(127);
        o.check_effective_scoring().unwrap();

        // The other ranges the parser holds, each guarding the kernels rather
        // than parity: `-w` sizes the band rows, `-U` is subtracted from a
        // score in `int`, and `-c` divides the SA-resolve step.
        let max_w = i32::MAX >> 2;
        for (set, ok, bad, needle) in [
            (
                MemOpts::set_band_width as fn(&mut MemOpts, i32) -> &mut MemOpts,
                [0, max_w],
                [-1, max_w + 1],
                "-w band width must be an integer in 0..",
            ),
            (
                MemOpts::set_unpaired_penalty,
                [0, i32::MAX],
                [-1, i32::MIN],
                "-U unpaired penalty must be a non-negative integer",
            ),
            (
                MemOpts::set_max_occurrences,
                [1, i32::MAX],
                [0, i32::MIN],
                "-c max occurrences must be a positive integer",
            ),
        ] {
            for v in ok {
                let mut o = MemOpts::new().unwrap();
                set(&mut o, v);
                o.check_effective_scoring().unwrap();
            }
            for v in bad {
                let mut o = MemOpts::new().unwrap();
                set(&mut o, v);
                err_has(&o, needle);
                err_has(&o, &format!("(got {v})"));
            }
        }
    }

    /// `set_mismatch_penalty` refills the matrices only for a penalty the
    /// kernels can store, so `i32::MIN` never reaches `bwa_fill_scmat`'s `-b`
    /// (signed overflow); the matrix then still reflects the last storable
    /// penalty, and alignment refuses the options.
    #[test]
    fn an_unstorable_mismatch_penalty_leaves_the_matrix_alone() {
        let mut o = MemOpts::new().unwrap();
        o.set_mismatch_penalty(128);
        assert_eq!(
            o.scoring_matrix()[1],
            -128,
            "128 is storable: off-diagonal -128"
        );
        let at_128 = all_matrices(&o);
        for b in [i32::MIN, -1, 129] {
            o.set_mismatch_penalty(b);
            assert_eq!(o.mismatch_penalty(), b);
            assert_eq!(
                all_matrices(&o),
                at_128,
                "-B {b} must not refill the matrix"
            );
            assert_invalid(
                o.check_effective_scoring(),
                &format!("effective mismatch penalty is {b}"),
            );
        }

        // Nor may a match-score change refill it while the stored penalty is
        // unstorable; the next storable penalty rebuilds from the current -A.
        o.set_mismatch_penalty(i32::MIN).set_match_score(2);
        assert_eq!(
            all_matrices(&o),
            at_128,
            "-A must not refill with -B i32::MIN"
        );
        o.set_mismatch_penalty(4);
        let m = o.scoring_matrix();
        assert_eq!((m[0], m[1]), (2, -4), "rebuilt from -A 2 and -B 4");
    }

    /// `apply_meth_defaults` refills the matrices in the shim, after the
    /// COLLAPSED bundle may have replaced `b`; it must skip that refill on the
    /// same terms as `refill_matrices` when the resulting `b` is unstorable.
    #[test]
    fn apply_meth_defaults_leaves_the_matrix_alone_for_an_unstorable_penalty() {
        for meth in [false, true] {
            for b in [i32::MIN, -1, 129] {
                let mut o = MemOpts::new().unwrap();
                o.set_meth(meth);
                if meth {
                    // GENOMIC keeps the caller's -B through the bundle.
                    o.set_meth_scoring(MethScoring::Genomic);
                }
                let before = all_matrices(&o);
                o.set_mismatch_penalty(b);
                o.apply_meth_defaults().unwrap();
                assert_eq!(o.mismatch_penalty(), b, "meth={meth}: -B {b} kept");
                assert_eq!(
                    all_matrices(&o),
                    before,
                    "meth={meth}: -B {b} must not refill the matrices"
                );
            }
        }

        // TAPS with no explicit scoring resolves to NEUTRAL inside the call,
        // which would rewrite the conversion cells of the meth matrices; that
        // refill is skipped too while -B is unstorable.
        for b in [i32::MIN, -1, 129] {
            let mut o = meth_opts(MethChem::Taps);
            o.set_mismatch_penalty(b);
            let before = all_matrices(&o);
            o.apply_meth_defaults().unwrap();
            assert_eq!(o.meth_scoring().unwrap(), MethScoring::Neutral);
            assert_eq!(
                all_matrices(&o),
                before,
                "implicit TAPS NEUTRAL, -B {b}: must not refill the matrices"
            );
        }

        // COLLAPSED replaces an unstorable -B with 2a, which then refills.
        let mut o = MemOpts::new().unwrap();
        o.set_meth(true).set_meth_scoring(MethScoring::Collapsed);
        o.set_mismatch_penalty(i32::MIN);
        o.apply_meth_defaults().unwrap();
        assert_eq!(o.mismatch_penalty(), 2);
        assert_eq!(o.scoring_matrix()[1], -2, "rebuilt from the bundle's -B 2");
    }

    /// `set_meth` and `set_meth_scoring` rebuild only the meth matrices, which
    /// are copies of the base matrix; while `-B` is unstorable they leave all of
    /// them alone like every other refill, and the next storable `-B` rebuilds
    /// both from the scoring mode chosen in between.
    #[test]
    fn meth_setters_leave_the_matrices_alone_for_an_unstorable_penalty() {
        let mut o = MemOpts::new().unwrap();
        o.set_mismatch_penalty(129);
        let before = all_matrices(&o);
        o.set_meth(true).set_meth_scoring(MethScoring::Neutral);
        assert_eq!(
            all_matrices(&o),
            before,
            "set_meth / set_meth_scoring must not refill with -B 129"
        );

        o.set_mismatch_penalty(4);
        let (ot, ob) = o.meth_scoring_matrices();
        assert_eq!(o.scoring_matrix()[1], -4, "rebuilt from -B 4");
        assert_eq!(ot[5 + 3], 0, "OT C->T cell rebuilt under NEUTRAL");
        assert_eq!(ob[2 * 5], 0, "OB G->A cell rebuilt under NEUTRAL");
    }

    /// The phase snapshot ignores only the emission-only knobs (`md_mode`,
    /// `MEM_F_XB`, `tie_break_mode`) and the two fields the shim overrides;
    /// any alignment-bearing change makes it differ.
    #[test]
    fn phase_snapshot_ignores_only_emission_knobs() {
        let base = MemOpts::new().unwrap();
        let seeded = base.phase_snapshot();
        base.check_phase(&seeded, "t").unwrap();

        let mut emit_only = MemOpts::new().unwrap();
        emit_only.set_md_mode(MdMode::None);
        // SAFETY: a live, uniquely owned mem_opt_t; sets one bit of a POD field.
        unsafe { (*emit_only.handle).flag |= sys::MEM_F_XB as i32 };
        emit_only.check_phase(&seeded, "t").unwrap();
        // --tie-break: read only at emit, by mem_tie_break_id.
        let tie = MemOpts::new().unwrap();
        // SAFETY: a live, uniquely owned mem_opt_t; writes one POD field.
        unsafe { (*tie.handle).tie_break_mode = 1 };
        tie.check_phase(&seeded, "t").unwrap();

        // The two fields the shim overrides on its own per-call copy are
        // blanked too, each on its own.
        let mut pe = MemOpts::new().unwrap();
        pe.set_pe(true);
        pe.check_phase(&seeded, "t").unwrap();
        let threads = MemOpts::new().unwrap();
        // SAFETY: a live, uniquely owned mem_opt_t; writes one POD field.
        unsafe { (*threads.handle).n_threads = 7 };
        threads.check_phase(&seeded, "t").unwrap();

        let changed: [fn(&mut MemOpts); 5] = [
            |o| {
                // A flag bit next to the blanked ones still counts.
                o.set_soft_clip_supplementary(true);
            },
            |o| {
                o.set_mismatch_penalty(5);
            },
            |o| {
                o.set_meth_context(MethContext::Cpg);
            },
            |o| {
                o.set_meth_chem(MethChem::Taps);
            },
            |o| {
                o.set_minimum_score(31);
            },
        ];
        for (i, change) in changed.iter().enumerate() {
            let mut o = MemOpts::new().unwrap();
            change(&mut o);
            let Err(err) = o.check_phase(&seeded, "pair_emit") else {
                panic!("change {i} went unnoticed");
            };
            assert!(
                matches!(&err, crate::Error::InvalidOpts(m)
                    if m.starts_with("pair_emit: the options differ")),
                "{err}"
            );
        }
    }

    #[test]
    fn cpg_context_with_taps_takes_the_neutral_default() {
        let mut o = meth_opts(MethChem::Taps);
        o.set_meth_context(MethContext::Cpg);
        o.apply_meth_defaults().unwrap();
        assert_eq!(o.meth_scoring().unwrap(), MethScoring::Neutral);
        assert_eq!(o.meth_seed_prune().unwrap(), MethSeedPrune::Spec30);
    }

    /// The three `--meth-context cpg` combinations `bwa-mem3 mem` refuses with
    /// `exit(1)` come back as errors, and leave the options as they were (the
    /// defaults are not half-applied).
    #[test]
    fn invalid_cpg_context_is_an_error_and_changes_nothing() {
        let mut no_meth = MemOpts::new().unwrap();
        no_meth.set_meth_context(MethContext::Cpg);
        let mut emseq = meth_opts(MethChem::EmSeq);
        emseq.set_meth_context(MethContext::Cpg);
        let mut collapsed = meth_opts(MethChem::Taps);
        collapsed
            .set_meth_context(MethContext::Cpg)
            .set_meth_scoring(MethScoring::Collapsed);
        for (mut o, needle) in [
            (no_meth, "--meth-context requires --meth"),
            (emseq, "requires --meth=taps or --meth=5base"),
            (collapsed, "incompatible with --meth-scoring collapsed"),
        ] {
            let before = raw_bytes(&o);
            assert_invalid(o.apply_meth_defaults(), needle);
            assert!(raw_bytes(&o) == before, "{needle}: options changed");
        }
    }

    /// As in `bwa-mem3 mem`, the `--meth` defaults apply only under `--meth`.
    #[test]
    fn meth_defaults_are_skipped_without_meth() {
        let mut o = MemOpts::new().unwrap();
        o.set_meth_chem(MethChem::Taps);
        let before = raw_bytes(&o);
        o.apply_meth_defaults().unwrap();
        // Not one field: no TAPS scoring default, no SPEC30, no -T/-L/-U/-M,
        // and the matrix refill rebuilds the same matrices.
        assert!(
            raw_bytes(&o) == before,
            "the defaults applied without --meth"
        );
    }

    #[test]
    fn explicit_meth_seed_prune_survives_the_defaults() {
        for prune in [MethSeedPrune::Off, MethSeedPrune::Baseline] {
            let mut o = meth_opts(MethChem::EmSeq);
            o.set_meth_seed_prune(prune).apply_meth_defaults().unwrap();
            assert_eq!(o.meth_seed_prune().unwrap(), prune);
        }
    }
}
