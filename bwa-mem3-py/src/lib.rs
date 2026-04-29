//! Python bindings for `bwa-mem3-rs`.
//!
//! Exposes the safe Rust API as a CPython extension module
//! (`bwa_mem3._bwa_mem3`). The thin facade in `python/bwa_mem3/__init__.py`
//! re-exports everything for the canonical `import bwa_mem3` namespace.

// pyo3 0.22's #[pymethods] macro injects an Into::into on the return
// type that newer clippy flags as `useless_conversion`. Allow at crate
// level rather than dotting #[allow] across every method.
#![allow(clippy::useless_conversion)]

use std::path::PathBuf;
use std::sync::Arc;

use bwa_mem3_rs as bwa;
use pyo3::exceptions::{PyRuntimeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyByteArray, PyBytes, PyList, PyTuple};

fn map_err(e: bwa::Error) -> PyErr {
    match e {
        bwa::Error::InvalidInput(m) => PyValueError::new_err(m),
        other => PyRuntimeError::new_err(other.to_string()),
    }
}

// ---------- BwaIndex ----------

/// Reference index handle.
///
/// `BwaIndex(prefix)` loads from disk; if a shared-memory segment for
/// `prefix` is staged, bwa-mem3 transparently attaches to it. Callers
/// that want fail-fast behavior should probe with `shm.is_staged`
/// first.
#[pyclass(name = "BwaIndex", module = "bwa_mem3._bwa_mem3", frozen)]
struct PyBwaIndex {
    inner: Arc<bwa::BwaIndex>,
}

#[pymethods]
impl PyBwaIndex {
    #[new]
    fn new(prefix: PathBuf) -> PyResult<Self> {
        let inner = bwa::BwaIndex::load(&prefix).map_err(map_err)?;
        Ok(Self {
            inner: Arc::new(inner),
        })
    }

    fn n_contigs(&self) -> usize {
        self.inner.n_contigs()
    }

    /// Returns `[(name, length), ...]` for every contig in the index.
    fn contigs<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let list = PyList::empty_bound(py);
        for (name, len) in self.inner.contigs() {
            let t = PyTuple::new_bound(py, [name.into_py(py), len.into_py(py)]);
            list.append(t)?;
        }
        Ok(list)
    }

    fn __repr__(&self) -> String {
        format!("BwaIndex(n_contigs={})", self.inner.n_contigs())
    }
}

// ---------- MemOpts ----------

/// Bwa-mem3 alignment options.
///
/// Constructed with bwa-mem3 defaults; mutate via setters. Setters return
/// `None` (in-place mutation) to follow Python conventions; chain via the
/// constructor + setattr pattern if needed.
#[pyclass(name = "MemOpts", module = "bwa_mem3._bwa_mem3")]
struct PyMemOpts {
    inner: bwa::MemOpts,
}

#[pymethods]
impl PyMemOpts {
    #[new]
    fn new() -> PyResult<Self> {
        let inner = bwa::MemOpts::new().map_err(map_err)?;
        Ok(Self { inner })
    }

    /// Apply a `-x` preset on top of current values.
    /// `mode` is one of `"pacbio"`, `"ont2d"`, `"intractg"`.
    fn apply_mode(&mut self, mode: &str) -> PyResult<()> {
        let m = match mode {
            "pacbio" => bwa::Mode::Pacbio,
            "ont2d" => bwa::Mode::Ont2d,
            "intractg" => bwa::Mode::Intractg,
            other => {
                return Err(PyValueError::new_err(format!(
                    "unknown mode {other:?}; expected one of 'pacbio', 'ont2d', 'intractg'"
                )))
            }
        };
        // Drop-and-replace: with_mode consumes self in Rust.
        let inner = std::mem::replace(&mut self.inner, bwa::MemOpts::new().map_err(map_err)?);
        self.inner = inner.with_mode(m);
        Ok(())
    }

    fn set_pe(&mut self, is_pe: bool) {
        self.inner.set_pe(is_pe);
    }

    fn set_soft_clip_supplementary(&mut self, v: bool) {
        self.inner.set_soft_clip_supplementary(v);
    }

    #[getter]
    fn min_seed_len(&self) -> i32 {
        self.inner.min_seed_len()
    }
    #[setter]
    fn set_min_seed_len(&mut self, v: i32) {
        self.inner.set_min_seed_len(v);
    }

    #[getter]
    fn band_width(&self) -> i32 {
        self.inner.band_width()
    }
    #[setter]
    fn set_band_width(&mut self, v: i32) {
        self.inner.set_band_width(v);
    }

    #[getter]
    fn match_score(&self) -> i32 {
        self.inner.match_score()
    }
    #[setter]
    fn set_match_score(&mut self, v: i32) {
        self.inner.set_match_score(v);
    }

    #[getter]
    fn mismatch_penalty(&self) -> i32 {
        self.inner.mismatch_penalty()
    }
    #[setter]
    fn set_mismatch_penalty(&mut self, v: i32) {
        self.inner.set_mismatch_penalty(v);
    }

    fn set_gap_open(&mut self, del: i32, ins: i32) {
        self.inner.set_gap_open(del, ins);
    }
    fn set_gap_extend(&mut self, del: i32, ins: i32) {
        self.inner.set_gap_extend(del, ins);
    }
    fn set_clip_penalty(&mut self, five: i32, three: i32) {
        self.inner.set_clip_penalty(five, three);
    }

    #[getter]
    fn minimum_score(&self) -> i32 {
        self.inner.minimum_score()
    }
    #[setter]
    fn set_minimum_score(&mut self, v: i32) {
        self.inner.set_minimum_score(v);
    }

    #[getter]
    fn max_occurrences(&self) -> i32 {
        self.inner.max_occurrences()
    }
    #[setter]
    fn set_max_occurrences(&mut self, v: i32) {
        self.inner.set_max_occurrences(v);
    }

    fn set_xa_max_hits(&mut self, primary: i32, alt: i32) {
        self.inner.set_xa_max_hits(primary, alt);
    }
    fn set_xa_drop_ratio(&mut self, v: f32) {
        self.inner.set_xa_drop_ratio(v);
    }
    fn set_unpaired_penalty(&mut self, v: i32) {
        self.inner.set_unpaired_penalty(v);
    }

    /// Set the `@RG` ID emitted as `RG:Z:` on every record. `None` clears.
    /// Note: writes to a process-wide global in bwa-mem3.
    #[pyo3(signature = (id=None))]
    fn set_read_group_id(&mut self, id: Option<&str>) -> PyResult<()> {
        self.inner.set_read_group_id(id).map_err(map_err)?;
        Ok(())
    }
}

// ---------- PeOrient / MemPeStat ----------

/// Insert-size statistics for a single orientation.
#[pyclass(name = "PeOrient", module = "bwa_mem3._bwa_mem3", get_all, set_all)]
#[derive(Clone)]
struct PyPeOrient {
    low: i32,
    high: i32,
    failed: bool,
    avg: f64,
    std: f64,
}

#[pymethods]
impl PyPeOrient {
    #[new]
    #[pyo3(signature = (low=0, high=0, failed=false, avg=0.0, std=0.0))]
    fn new(low: i32, high: i32, failed: bool, avg: f64, std: f64) -> Self {
        Self {
            low,
            high,
            failed,
            avg,
            std,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "PeOrient(low={}, high={}, failed={}, avg={}, std={})",
            self.low, self.high, self.failed, self.avg, self.std
        )
    }
}

impl From<bwa::PeOrient> for PyPeOrient {
    fn from(p: bwa::PeOrient) -> Self {
        Self {
            low: p.low,
            high: p.high,
            failed: p.failed,
            avg: p.avg,
            std: p.std,
        }
    }
}

impl From<PyPeOrient> for bwa::PeOrient {
    fn from(p: PyPeOrient) -> Self {
        Self {
            low: p.low,
            high: p.high,
            failed: p.failed,
            avg: p.avg,
            std: p.std,
        }
    }
}

fn orientation_from_str(s: &str) -> PyResult<bwa::PeOrientation> {
    match s.to_ascii_uppercase().as_str() {
        "FF" => Ok(bwa::PeOrientation::Ff),
        "FR" => Ok(bwa::PeOrientation::Fr),
        "RF" => Ok(bwa::PeOrientation::Rf),
        "RR" => Ok(bwa::PeOrientation::Rr),
        other => Err(PyValueError::new_err(format!(
            "unknown PE orientation {other:?}; expected one of 'FF', 'FR', 'RF', 'RR'"
        ))),
    }
}

/// 4-orientation paired-end insert-size model (`mem_pestat_t[4]`).
#[pyclass(name = "MemPeStat", module = "bwa_mem3._bwa_mem3")]
struct PyMemPeStat {
    inner: bwa::MemPeStat,
}

#[pymethods]
impl PyMemPeStat {
    #[new]
    fn zero() -> PyResult<Self> {
        Ok(Self {
            inner: bwa::MemPeStat::zero().map_err(map_err)?,
        })
    }

    /// `o` is one of `"FF"`, `"FR"`, `"RF"`, `"RR"`.
    fn orientation(&self, o: &str) -> PyResult<PyPeOrient> {
        let p = self.inner.orientation(orientation_from_str(o)?);
        Ok(p.into())
    }

    fn set_orientation(&mut self, o: &str, v: PyPeOrient) -> PyResult<()> {
        self.inner
            .set_orientation(orientation_from_str(o)?, v.into());
        Ok(())
    }
}

// ---------- ReadPair / Record ----------

/// Copy any bytes-like Python object (`bytes`, `bytearray`, `memoryview` over a
/// contiguous u8 buffer, array.array('B'), numpy uint8 arrays, …) into an
/// owned `Box<[u8]>`.
///
/// The copy decouples the resulting buffer's lifetime from any Python object,
/// which is what lets `align_batch` & friends release the GIL during the bwa
/// call — Python is free to mutate or drop the original object in the
/// meantime.
///
/// `bytes` and `bytearray` are handled directly; anything else is routed
/// through `memoryview(obj).tobytes()`. `memoryview` rejects anything
/// that doesn't expose the buffer protocol (including `int` and iterables
/// of ints, which `bytes(obj)` would silently accept and do the wrong
/// thing with).
fn to_owned_bytes(obj: &Bound<'_, PyAny>) -> PyResult<Box<[u8]>> {
    if let Ok(b) = obj.downcast::<PyBytes>() {
        return Ok(b.as_bytes().to_vec().into_boxed_slice());
    }
    if let Ok(b) = obj.downcast::<PyByteArray>() {
        // SAFETY: PyByteArray contents may be mutated under the GIL, but
        // we copy out immediately before any other Python code can run,
        // so the aliased view is stable for the duration of the copy.
        let slice = unsafe { b.as_bytes() };
        return Ok(slice.to_vec().into_boxed_slice());
    }
    let py = obj.py();
    let mv = py
        .import_bound("builtins")
        .and_then(|b| b.getattr("memoryview"))
        .and_then(|mv_ty| mv_ty.call1((obj,)))
        .map_err(|_| {
            PyTypeError::new_err(format!(
                "expected a bytes-like object (bytes, bytearray, memoryview, …); got {}",
                obj.get_type()
                    .name()
                    .map(|s| s.to_string())
                    .unwrap_or_default()
            ))
        })?;
    let coerced = mv.call_method0("tobytes")?;
    let coerced = coerced.downcast_into::<PyBytes>().map_err(|_| {
        PyTypeError::new_err("memoryview.tobytes() did not return bytes (this should not happen)")
    })?;
    Ok(coerced.as_bytes().to_vec().into_boxed_slice())
}

fn to_owned_bytes_opt(obj: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Box<[u8]>>> {
    obj.map(to_owned_bytes).transpose()
}

/// One paired-end read. Owns its name/seq/qual bytes outright (copied from
/// the input bytes-like objects at construction). The copy lets alignment
/// entry points release the GIL during the native call without any
/// borrow-from-Python lifetime concerns.
#[pyclass(name = "ReadPair", module = "bwa_mem3._bwa_mem3", frozen)]
struct PyReadPair {
    name_r1: Box<[u8]>,
    seq_r1: Box<[u8]>,
    qual_r1: Option<Box<[u8]>>,
    name_r2: Box<[u8]>,
    seq_r2: Box<[u8]>,
    qual_r2: Option<Box<[u8]>>,
}

#[pymethods]
impl PyReadPair {
    /// Construct a paired-end read from any bytes-like objects.
    ///
    /// Each field accepts `bytes`, `bytearray`, `memoryview` over a contiguous
    /// `u8` buffer, or any object implementing the buffer protocol with
    /// itemsize 1. Contents are copied into Rust-owned storage; callers may
    /// reuse / mutate the input objects after construction.
    #[new]
    #[pyo3(signature = (name_r1, seq_r1, name_r2, seq_r2, qual_r1=None, qual_r2=None))]
    fn new(
        name_r1: &Bound<'_, PyAny>,
        seq_r1: &Bound<'_, PyAny>,
        name_r2: &Bound<'_, PyAny>,
        seq_r2: &Bound<'_, PyAny>,
        qual_r1: Option<&Bound<'_, PyAny>>,
        qual_r2: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Ok(Self {
            name_r1: to_owned_bytes(name_r1)?,
            seq_r1: to_owned_bytes(seq_r1)?,
            qual_r1: to_owned_bytes_opt(qual_r1)?,
            name_r2: to_owned_bytes(name_r2)?,
            seq_r2: to_owned_bytes(seq_r2)?,
            qual_r2: to_owned_bytes_opt(qual_r2)?,
        })
    }
}

/// One packed BAM record produced by alignment.
#[pyclass(name = "Record", module = "bwa_mem3._bwa_mem3", frozen, get_all)]
struct PyRecord {
    pair_idx: usize,
    bytes: Py<PyBytes>,
}

#[pymethods]
impl PyRecord {
    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "Record(pair_idx={}, bytes=<{} bytes>)",
            self.pair_idx,
            self.bytes.bind(py).len().unwrap_or(0)
        )
    }
}

// ---------- Alignment functions ----------

/// Borrow byte slices from a list of `PyReadPair` into a `Vec<bwa::ReadPair>`
/// suitable for passing to the native align entry points.
///
/// The returned slices reference Rust-owned `Box<[u8]>` storage on each
/// `PyReadPair`. They remain valid as long as every `PyRef` in `pairs` is
/// alive — which keeps the underlying Python objects (and thus their owned
/// byte buffers) rooted — and crucially do **not** depend on the GIL,
/// because the buffer storage is Rust-side. Callers may keep the returned
/// `Vec` across a `py.allow_threads(...)` block, provided `pairs` itself
/// lives across the block (so the `PyRef` guards aren't dropped without
/// the GIL).
fn borrow_pairs<'a>(pairs: &'a [PyRef<'_, PyReadPair>]) -> Vec<bwa::ReadPair<'a>> {
    pairs
        .iter()
        .map(|p| bwa::ReadPair {
            name_r1: &p.name_r1,
            seq_r1: &p.seq_r1,
            qual_r1: p.qual_r1.as_deref(),
            name_r2: &p.name_r2,
            seq_r2: &p.seq_r2,
            qual_r2: p.qual_r2.as_deref(),
        })
        .collect()
}

fn batch_to_records<'py>(
    py: Python<'py>,
    batch: &bwa::AlignmentBatch,
) -> PyResult<Bound<'py, PyList>> {
    let list = PyList::empty_bound(py);
    for rec in batch.iter() {
        let bytes = PyBytes::new_bound(py, rec.bytes).unbind();
        let py_rec = PyRecord {
            pair_idx: rec.pair_idx,
            bytes,
        };
        list.append(Bound::new(py, py_rec)?)?;
    }
    Ok(list)
}

/// Align a batch of read pairs. Returns `(records, pestat_out)`.
///
/// Releases the GIL for the duration of the native alignment call (input
/// bytes are owned by each `ReadPair`, see its constructor docstring), so
/// running `align_batch` concurrently from a
/// `concurrent.futures.ThreadPoolExecutor` against the same `BwaIndex`
/// yields true parallelism — at the cost of one copy of each input
/// name/seq/qual at `ReadPair` construction time. For sharing the index
/// across *process* workers as well, see [`shm.stage`].
#[pyfunction]
#[pyo3(signature = (index, opts, pairs, pestat_in=None))]
fn align_batch<'py>(
    py: Python<'py>,
    index: &PyBwaIndex,
    opts: &PyMemOpts,
    pairs: Vec<PyRef<'py, PyReadPair>>,
    pestat_in: Option<&PyMemPeStat>,
) -> PyResult<(Bound<'py, PyList>, PyMemPeStat)> {
    let borrowed = borrow_pairs(&pairs);
    let pestat_ref = pestat_in.map(|p| &p.inner);
    let index_inner = &index.inner;
    let opts_inner = &opts.inner;
    let (batch, pestat_out) = py
        .allow_threads(|| bwa::align_batch(index_inner, opts_inner, &borrowed, pestat_ref))
        .map_err(map_err)?;
    let records = batch_to_records(py, &batch)?;
    Ok((records, PyMemPeStat { inner: pestat_out }))
}

/// Phase 1: seed only. Returns an opaque handle consumed by `extend_batch`.
///
/// Releases the GIL during the native seeding call (see `align_batch` for
/// the threading rationale).
#[pyfunction]
fn seed_batch<'py>(
    py: Python<'py>,
    index: &PyBwaIndex,
    opts: &PyMemOpts,
    pairs: Vec<PyRef<'py, PyReadPair>>,
) -> PyResult<PySeeds> {
    let borrowed = borrow_pairs(&pairs);
    let index_inner = &index.inner;
    let opts_inner = &opts.inner;
    let seeds = py
        .allow_threads(|| bwa::seed_batch(index_inner, opts_inner, &borrowed))
        .map_err(map_err)?;
    Ok(PySeeds { inner: Some(seeds) })
}

/// Phase 2: extend pre-computed seeds to alignments. Consumes `seeds`.
///
/// Releases the GIL during the native extend call (see `align_batch` for
/// the threading rationale).
#[pyfunction]
#[pyo3(signature = (index, opts, seeds, pairs, pestat_in=None))]
fn extend_batch<'py>(
    py: Python<'py>,
    index: &PyBwaIndex,
    opts: &PyMemOpts,
    seeds: &mut PySeeds,
    pairs: Vec<PyRef<'py, PyReadPair>>,
    pestat_in: Option<&PyMemPeStat>,
) -> PyResult<(Bound<'py, PyList>, PyMemPeStat)> {
    let inner = seeds
        .inner
        .take()
        .ok_or_else(|| PyValueError::new_err("seeds already consumed"))?;
    let borrowed = borrow_pairs(&pairs);
    let pestat_ref = pestat_in.map(|p| &p.inner);
    let index_inner = &index.inner;
    let opts_inner = &opts.inner;
    let (batch, pestat_out) = py
        .allow_threads(|| bwa::extend_batch(index_inner, opts_inner, inner, &borrowed, pestat_ref))
        .map_err(map_err)?;
    let records = batch_to_records(py, &batch)?;
    Ok((records, PyMemPeStat { inner: pestat_out }))
}

/// Run seeding + SE-extension on a pilot batch and return the inferred
/// 4-orientation insert-size model. Discards the alignments.
///
/// Releases the GIL during the native call.
#[pyfunction]
fn estimate_pestat<'py>(
    py: Python<'py>,
    index: &PyBwaIndex,
    opts: &PyMemOpts,
    pairs: Vec<PyRef<'py, PyReadPair>>,
) -> PyResult<PyMemPeStat> {
    let borrowed = borrow_pairs(&pairs);
    let index_inner = &index.inner;
    let opts_inner = &opts.inner;
    let pestat = py
        .allow_threads(|| bwa::estimate_pestat(index_inner, opts_inner, &borrowed))
        .map_err(map_err)?;
    Ok(PyMemPeStat { inner: pestat })
}

/// Opaque seeds handle from `seed_batch`. Single-shot — consumed by `extend_batch`.
#[pyclass(name = "Seeds", module = "bwa_mem3._bwa_mem3")]
struct PySeeds {
    inner: Option<bwa::Seeds>,
}

#[pymethods]
impl PySeeds {
    fn __repr__(&self) -> &'static str {
        if self.inner.is_some() {
            "Seeds(<live>)"
        } else {
            "Seeds(<consumed>)"
        }
    }
}

// ---------- shm submodule ----------

#[pyfunction(name = "is_staged")]
fn shm_is_staged(prefix: PathBuf) -> PyResult<bool> {
    bwa::shm::is_staged(&prefix).map_err(map_err)
}

#[pyfunction(name = "stage")]
fn shm_stage(prefix: PathBuf) -> PyResult<()> {
    bwa::shm::stage(&prefix).map_err(map_err)
}

#[pyfunction(name = "destroy")]
fn shm_destroy() -> PyResult<()> {
    bwa::shm::destroy().map_err(map_err)
}

#[pyfunction(name = "list")]
fn shm_list() -> PyResult<()> {
    bwa::shm::list().map_err(map_err)
}

// ---------- module init ----------

#[pymodule]
fn _bwa_mem3(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBwaIndex>()?;
    m.add_class::<PyMemOpts>()?;
    m.add_class::<PyMemPeStat>()?;
    m.add_class::<PyPeOrient>()?;
    m.add_class::<PyReadPair>()?;
    m.add_class::<PyRecord>()?;
    m.add_class::<PySeeds>()?;
    m.add_function(wrap_pyfunction!(align_batch, m)?)?;
    m.add_function(wrap_pyfunction!(seed_batch, m)?)?;
    m.add_function(wrap_pyfunction!(extend_batch, m)?)?;
    m.add_function(wrap_pyfunction!(estimate_pestat, m)?)?;

    let shm = PyModule::new_bound(m.py(), "shm")?;
    shm.add_function(wrap_pyfunction!(shm_is_staged, &shm)?)?;
    shm.add_function(wrap_pyfunction!(shm_stage, &shm)?)?;
    shm.add_function(wrap_pyfunction!(shm_destroy, &shm)?)?;
    shm.add_function(wrap_pyfunction!(shm_list, &shm)?)?;
    m.add_submodule(&shm)?;
    Ok(())
}
