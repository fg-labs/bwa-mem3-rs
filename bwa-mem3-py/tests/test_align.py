"""
Integration tests that require a real bwa-mem3 index.

Index source priority (see `conftest.py::aligned_index`):
1. `BWA_MEM3_RS_TEST_REF` — a pre-built large index (hg38, etc.).
2. PhiX174 built on the fly via `bwa-mem3 index` (CI default).

Skips if neither path is available.
"""

from __future__ import annotations

import concurrent.futures
import os
import sys
from pathlib import Path

import pytest
from bwa_mem3 import BwaIndex
from bwa_mem3 import MemOpts
from bwa_mem3 import ReadPair
from bwa_mem3 import align_batch
from bwa_mem3 import estimate_pestat
from bwa_mem3 import extend_batch
from bwa_mem3 import seed_batch
from bwa_mem3 import shm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bam import decode as bam_decode  # noqa: E402
from conftest import simulate_pairs  # noqa: E402
from phix_seq import PHIX_SEQ  # noqa: E402


def _synthetic_pair(name: bytes) -> ReadPair:
    seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    qual = b"I" * len(seq)
    return ReadPair(name_r1=name, seq_r1=seq, qual_r1=qual, name_r2=name, seq_r2=seq, qual_r2=qual)


def test_load_index_and_align(aligned_index: str) -> None:
    idx = BwaIndex(aligned_index)
    assert idx.n_contigs() > 0
    contigs = idx.contigs()
    assert all(isinstance(c, tuple) and len(c) == 2 for c in contigs)
    assert all(isinstance(c[0], str) and isinstance(c[1], int) for c in contigs)

    opts = MemOpts()
    opts.set_pe(True)

    pairs = [_synthetic_pair(b"r0"), _synthetic_pair(b"r1")]
    records, _pestat = align_batch(idx, opts, pairs)
    assert len(records) >= 2 * len(pairs)
    for rec in records:
        assert isinstance(rec.bytes, bytes)
        assert rec.pair_idx in (0, 1)


def test_phase_split_matches_align_batch(aligned_index: str) -> None:
    idx = BwaIndex(aligned_index)
    opts = MemOpts()
    opts.set_pe(True)
    pairs = [_synthetic_pair(b"x")]

    seeds = seed_batch(idx, opts, pairs)
    rec_phased, _ = extend_batch(idx, opts, seeds, pairs)

    rec_one_shot, _ = align_batch(idx, opts, pairs)

    # Same record count and same record bytes per pair (orderings match
    # because both paths run a single thread serially).
    assert len(rec_phased) == len(rec_one_shot)
    for a, b in zip(rec_phased, rec_one_shot):
        assert a.pair_idx == b.pair_idx
        assert a.bytes == b.bytes


def test_seeds_consumed_only_once(aligned_index: str) -> None:
    """extend_batch consumes the Seeds; reusing it raises ValueError."""
    idx = BwaIndex(aligned_index)
    opts = MemOpts()
    opts.set_pe(True)
    pairs = [_synthetic_pair(b"x")]

    seeds = seed_batch(idx, opts, pairs)
    extend_batch(idx, opts, seeds, pairs)
    with pytest.raises(ValueError):
        extend_batch(idx, opts, seeds, pairs)


def test_estimate_pestat_then_align_with_pestat_in(aligned_index: str) -> None:
    """estimate_pestat -> align_batch with pestat_in round-trip."""
    idx = BwaIndex(aligned_index)
    opts = MemOpts()
    opts.set_pe(True)
    pairs = [_synthetic_pair(f"r{i}".encode()) for i in range(4)]

    pestat = estimate_pestat(idx, opts, pairs)
    records, _out_pestat = align_batch(idx, opts, pairs, pestat_in=pestat)
    assert len(records) >= 2 * len(pairs)


def test_align_batch_releases_gil_under_thread_pool(aligned_index: str) -> None:
    """
    align_batch releases the GIL under concurrent use.

    Concurrent calls against the same BwaIndex from a ThreadPoolExecutor
    produce identical results to a serial run (and don't deadlock).
    """
    idx = BwaIndex(aligned_index)
    opts = MemOpts()
    opts.set_pe(True)

    def make_batch(seed: int) -> list[ReadPair]:
        return [_synthetic_pair(f"r{seed}-{i}".encode()) for i in range(4)]

    batches = [make_batch(i) for i in range(4)]

    serial = [align_batch(idx, opts, batch) for batch in batches]

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        parallel = list(ex.map(lambda b: align_batch(idx, opts, b), batches))

    assert len(parallel) == len(serial)
    for (serial_recs, _), (parallel_recs, _) in zip(serial, parallel):
        assert len(parallel_recs) == len(serial_recs)
        for a, b in zip(serial_recs, parallel_recs):
            assert a.pair_idx == b.pair_idx
            assert a.bytes == b.bytes


# ---- additional coverage ----


def test_align_batch_empty_pairs(aligned_index: str) -> None:
    """Aligning an empty batch returns no records and doesn't crash."""
    idx = BwaIndex(aligned_index)
    opts = MemOpts()
    opts.set_pe(True)
    records, _pestat = align_batch(idx, opts, [])
    assert list(records) == []


def test_align_batch_single_end(aligned_index: str) -> None:
    """With `set_pe(False)`, alignment still produces records."""
    idx = BwaIndex(aligned_index)
    opts = MemOpts()
    opts.set_pe(False)
    pairs = [_synthetic_pair(b"r0"), _synthetic_pair(b"r1")]
    records, _pestat = align_batch(idx, opts, pairs)
    # Two reads per pair regardless of PE flag; orientation/pairing
    # bookkeeping just changes the emitted flags.
    assert len(records) >= 2 * len(pairs)


def test_read_group_id_emitted_in_records(aligned_index: str) -> None:
    """`set_read_group_id` writes an RG:Z aux tag on every emitted record."""
    rg_id = "bwa-mem3-py-test-rg"
    idx = BwaIndex(aligned_index)
    opts = MemOpts()
    opts.set_pe(True)
    try:
        opts.set_read_group_id(rg_id)
        pairs = [_synthetic_pair(b"r0"), _synthetic_pair(b"r1")]
        records, _ = align_batch(idx, opts, pairs)
        assert len(records) > 0
        for rec in records:
            decoded = bam_decode(rec.bytes)
            assert decoded.aux_z(b"RG") == rg_id, (
                f"expected RG:Z:{rg_id} on every record; got {decoded.aux!r}"
            )
    finally:
        # `set_read_group_id` writes a bwa-mem3 process-wide global; clear it
        # so it doesn't leak into other tests.
        opts.set_read_group_id(None)


def test_mismatched_seq_qual_length_surfaces_error(aligned_index: str) -> None:
    """
    A read pair where len(seq) != len(qual) is rejected at align time.

    `ReadPair.__init__` only copies bytes; the length-mismatch check lives
    in the Rust `ReadPair::validate`, which `align_batch` (and friends) run
    before touching the native aligner. The mismatch surfaces as a Python
    `ValueError` (mapped from `bwa_mem3_rs::Error::InvalidInput`).
    """
    idx = BwaIndex(aligned_index)
    opts = MemOpts()
    opts.set_pe(True)
    bad = ReadPair(
        name_r1=b"r0",
        seq_r1=b"ACGTACGT",
        qual_r1=b"III",  # length mismatch
        name_r2=b"r0",
        seq_r2=b"ACGTACGT",
        qual_r2=b"IIIIIIII",
    )
    with pytest.raises(ValueError):
        align_batch(idx, opts, [bad])


def test_aligns_simulated_phix_reads(aligned_index: str) -> None:
    """
    Most simulated PhiX reads should align (FLAG&UNMAPPED == 0).

    Only meaningful when the index is the bundled PhiX (since the simulated
    reads are sampled from the PhiX sequence). With a different
    `BWA_MEM3_RS_TEST_REF`, the read source mismatches the index and the
    test is skipped.
    """
    if os.environ.get("BWA_MEM3_RS_TEST_REF"):
        pytest.skip("simulated reads are PhiX-specific; large-index path skipped")
    pairs_data = simulate_pairs(PHIX_SEQ.encode(), n=50)
    pairs = [
        ReadPair(
            name_r1=name,
            seq_r1=r1,
            qual_r1=b"I" * len(r1),
            name_r2=name,
            seq_r2=r2,
            qual_r2=b"I" * len(r2),
        )
        for name, r1, r2 in pairs_data
    ]
    idx = BwaIndex(aligned_index)
    opts = MemOpts()
    opts.set_pe(True)
    records, _ = align_batch(idx, opts, pairs)
    decoded = [bam_decode(r.bytes) for r in records]
    n_mapped = sum(1 for d in decoded if not d.is_unmapped)
    # At least 80% of records (across both mates) should align — read_len=100
    # against a 5,386 bp single-contig reference is trivially uniquely placed.
    assert n_mapped >= int(0.8 * len(decoded)), (
        f"only {n_mapped}/{len(decoded)} records mapped; expected ≥ 80%"
    )


@pytest.mark.skipif(sys.platform != "linux", reason="bwa-mem3 shm uses Linux /dev/shm semantics")
def test_shm_stage_destroy_lifecycle(aligned_index: str) -> None:
    """
    Stage the index, see is_staged flip True, destroy, see it flip back.

    `shm.destroy()` is process-global (clears every bwa-mem3 segment owned by
    this user), so the test must own the staged state and clean up on every
    exit path.
    """
    assert shm.is_staged(aligned_index) is False
    try:
        shm.stage(aligned_index)
        assert shm.is_staged(aligned_index) is True
        # Loading the index transparently attaches to the staged segment.
        idx = BwaIndex(aligned_index)
        assert idx.n_contigs() > 0
        shm.list()  # smoke; just verify it doesn't raise
    finally:
        shm.destroy()
    assert shm.is_staged(aligned_index) is False
