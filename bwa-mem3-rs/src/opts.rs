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
    GlobalLongest = 1,
    LocalLongest = 2,
    AbsorbCount = 3,
    MostAbsorb = 4,
}

impl SeedOrder {
    /// Map the raw `seed_order_t` value to a [`SeedOrder`], falling back to
    /// [`SeedOrder::Off`] for any value bwa-mem3 does not define.
    fn from_raw(v: i32) -> Self {
        match v {
            1 => SeedOrder::GlobalLongest,
            2 => SeedOrder::LocalLongest,
            3 => SeedOrder::AbsorbCount,
            4 => SeedOrder::MostAbsorb,
            _ => SeedOrder::Off,
        }
    }
}

/// Bisulfite scoring mode (`--meth-scoring`), mirroring bwa-mem3's
/// `mem_meth_scoring`. Only meaningful under [`MemOpts::set_meth`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MethScoring {
    /// C↔T and G↔A interchangeable (both matrix cells freed) — reproduces
    /// bwameth. This is bwa-mem3's default.
    Collapsed = 0,
    /// Only the bisulfite-conversion cell is freed; the mirror stays a
    /// mismatch, so genuine variants still score as mismatches (variant-aware,
    /// truthful NM/MD).
    Genomic = 1,
}

impl MethScoring {
    /// Map the raw `mem_meth_scoring` value to a [`MethScoring`], falling back
    /// to [`MethScoring::Collapsed`] (the default) for undefined values.
    fn from_raw(v: i32) -> Self {
        match v {
            1 => MethScoring::Genomic,
            _ => MethScoring::Collapsed,
        }
    }
}

/// bwa-mem3 alignment options (`mem_opt_t`).
pub struct MemOpts {
    pub(crate) handle: *mut sys::mem_opt_t,
}

impl MemOpts {
    /// Create options with bwa-mem3 defaults (via `mem_opt_init`).
    pub fn new() -> Result<Self> {
        let handle = unsafe { sys::bwa_shim_opts_new() };
        if handle.is_null() {
            return Err(shim_err("opts_new"));
        }
        Ok(MemOpts { handle })
    }

    /// Apply a bwa-mem3 `-x` preset on top of current values.
    #[must_use]
    pub fn with_mode(self, mode: Mode) -> Self {
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
        self
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
    pub fn set_match_score(&mut self, v: i32) -> &mut Self {
        unsafe {
            (*self.handle).a = v;
        }
        self
    }

    #[must_use]
    pub fn mismatch_penalty(&self) -> i32 {
        unsafe { (*self.handle).b }
    }
    pub fn set_mismatch_penalty(&mut self, v: i32) -> &mut Self {
        unsafe {
            (*self.handle).b = v;
        }
        self
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
    #[must_use]
    pub fn seed_order(&self) -> SeedOrder {
        SeedOrder::from_raw(unsafe { (*self.handle).seed_emit_order })
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
    /// tags. Toggling this rebuilds the per-hypothesis scoring matrices.
    ///
    /// Aligning with `meth` enabled but a **non-meth** index (loaded via
    /// [`BwaIndex::load`](crate::BwaIndex::load)) falls back to normal
    /// (non-bisulfite) behavior, so always pair `set_meth(true)` with
    /// `load_meth`.
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
    /// [`set_meth(true)`](Self::set_meth). Rebuilds the scoring matrices.
    pub fn set_meth_scoring(&mut self, m: MethScoring) -> &mut Self {
        unsafe {
            (*self.handle).meth_scoring = m as i32;
            sys::bwa_shim_opts_fill_meth_mat(self.handle);
        }
        self
    }
    #[must_use]
    pub fn meth_scoring(&self) -> MethScoring {
        MethScoring::from_raw(unsafe { (*self.handle).meth_scoring })
    }

    /// bwameth.py-style longest-`M` chimera QC heuristic (`--meth-chimera-qc`);
    /// only meaningful under [`set_meth(true)`](Self::set_meth). Off by default.
    pub fn set_meth_chimera_qc(&mut self, v: bool) -> &mut Self {
        unsafe {
            (*self.handle).meth_chimera_qc = i32::from(v);
        }
        self
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
        assert!(!o.skip_contained_ext());
        assert_eq!(o.band_start(), 0);
        assert_eq!(o.seed_order(), SeedOrder::Off);
    }

    #[test]
    fn new_opt_in_knobs_roundtrip() {
        let mut o = MemOpts::new().unwrap();
        o.set_min_ext_len(25)
            .set_max_extend_chains(4)
            .set_mate_concordant_window(-1)
            .set_smem_dedup(true)
            .set_skip_contained_ext(true)
            .set_band_start(20)
            .set_seed_order(SeedOrder::MostAbsorb);
        assert_eq!(o.min_ext_len(), 25);
        assert_eq!(o.max_extend_chains(), 4);
        assert_eq!(o.mate_concordant_window(), -1);
        assert!(o.smem_dedup());
        assert!(o.skip_contained_ext());
        assert_eq!(o.band_start(), 20);
        assert_eq!(o.seed_order(), SeedOrder::MostAbsorb);
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
            assert_eq!(o.seed_order(), order);
        }
    }

    #[test]
    fn meth_knobs_roundtrip() {
        let mut o = MemOpts::new().unwrap();
        // Defaults: meth off, collapsed scoring.
        assert!(!o.meth());
        assert_eq!(o.meth_scoring(), MethScoring::Collapsed);
        o.set_meth(true)
            .set_meth_scoring(MethScoring::Genomic)
            .set_meth_chimera_qc(true);
        assert!(o.meth());
        assert_eq!(o.meth_scoring(), MethScoring::Genomic);
        o.set_meth_scoring(MethScoring::Collapsed);
        assert_eq!(o.meth_scoring(), MethScoring::Collapsed);
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
}
