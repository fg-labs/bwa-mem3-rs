"""
Smoke tests that don't require a real index.

The integration test in test_align.py loads a real index and aligns a
synthetic batch; that one is gated on BWA_MEM3_RS_TEST_REF.
"""

from __future__ import annotations

import gc
from pathlib import Path

import bwa_mem3
import pytest
from bwa_mem3 import MemOpts
from bwa_mem3 import MemPeStat
from bwa_mem3 import PeOrient
from bwa_mem3 import ReadPair
from bwa_mem3 import shm


def test_module_surface() -> None:
    assert hasattr(bwa_mem3, "BwaIndex")
    assert hasattr(bwa_mem3, "MemOpts")
    assert hasattr(bwa_mem3, "MemPeStat")
    assert hasattr(bwa_mem3, "ReadPair")
    assert hasattr(bwa_mem3, "Record")
    assert hasattr(bwa_mem3, "Seeds")
    assert callable(bwa_mem3.align_batch)
    assert callable(bwa_mem3.seed_batch)
    assert callable(bwa_mem3.extend_batch)
    assert callable(bwa_mem3.estimate_pestat)
    assert hasattr(bwa_mem3.shm, "is_staged")
    assert hasattr(bwa_mem3.shm, "stage")
    assert hasattr(bwa_mem3.shm, "destroy")
    assert hasattr(bwa_mem3.shm, "list")


def test_memopts_defaults() -> None:
    o = MemOpts()
    # bwa-mem3's mem_opt_init defaults.
    assert o.match_score == 1
    assert o.mismatch_penalty == 4
    assert o.band_width == 100
    assert o.minimum_score == 30
    assert o.min_seed_len == 19


def test_memopts_setters() -> None:
    o = MemOpts()
    o.min_seed_len = 21
    assert o.min_seed_len == 21
    o.set_pe(True)
    o.set_pe(False)
    o.set_gap_open(7, 7)
    o.set_xa_max_hits(3, 100)
    o.set_read_group_id("sample42")
    o.set_read_group_id(None)


def test_memopts_apply_mode_pacbio() -> None:
    o = MemOpts()
    o.apply_mode("pacbio")
    assert o.min_seed_len == 17
    assert o.mismatch_penalty == 1


def test_memopts_apply_mode_ont2d() -> None:
    o = MemOpts()
    o.apply_mode("ont2d")
    assert o.min_seed_len == 14
    assert o.mismatch_penalty == 1
    assert o.minimum_score == 20


def test_memopts_apply_mode_intractg() -> None:
    o = MemOpts()
    # `intractg` only touches b, o_del, o_ins, pen_clip5, pen_clip3; the seed
    # length is left at the default.
    default_seed = MemOpts().min_seed_len
    o.apply_mode("intractg")
    assert o.mismatch_penalty == 9
    assert o.min_seed_len == default_seed


def test_memopts_apply_mode_unknown() -> None:
    o = MemOpts()
    with pytest.raises(ValueError):
        o.apply_mode("nonsense")


def test_pestat_orientations() -> None:
    p = MemPeStat()
    p.set_orientation("FR", PeOrient(low=100, high=500, failed=False, avg=250.0, std=50.0))
    fr = p.orientation("FR")
    assert fr.low == 100
    assert fr.high == 500
    assert fr.avg == 250.0
    assert fr.failed is False

    # Untouched orientation is zero.
    rf = p.orientation("RF")
    assert rf.low == 0

    with pytest.raises(ValueError):
        p.orientation("XX")


def test_pestat_all_orientations_round_trip() -> None:
    """Round-trip every orientation independently; overwrites preserve only the new value."""
    p = MemPeStat()
    fixtures = {
        "FF": PeOrient(low=10, high=20, avg=15.0, std=2.5),
        "FR": PeOrient(low=100, high=200, avg=150.0, std=20.0),
        "RF": PeOrient(low=30, high=40, failed=True, avg=35.0, std=5.0),
        "RR": PeOrient(low=50, high=60, avg=55.0, std=5.0),
    }
    for name, value in fixtures.items():
        p.set_orientation(name, value)
    for name, expected in fixtures.items():
        got = p.orientation(name)
        assert got.low == expected.low
        assert got.high == expected.high
        assert got.failed is expected.failed
        assert got.avg == expected.avg
        assert got.std == expected.std

    # Overwrite FR; siblings unchanged.
    p.set_orientation("FR", PeOrient(low=999, high=1000))
    assert p.orientation("FR").low == 999
    assert p.orientation("FF").low == 10  # unaffected
    assert p.orientation("RR").low == 50


def test_readpair_construction() -> None:
    p = ReadPair(
        name_r1=b"r1",
        seq_r1=b"ACGT",
        qual_r1=b"IIII",
        name_r2=b"r2",
        seq_r2=b"ACGT",
        qual_r2=b"IIII",
    )
    # ReadPair is opaque; just verify it constructs.
    assert p is not None


def test_readpair_qual_optional() -> None:
    p = ReadPair(name_r1=b"r1", seq_r1=b"ACGT", name_r2=b"r2", seq_r2=b"ACGT")
    assert p is not None


def test_readpair_owns_input_bytes() -> None:
    """ReadPair retains its bytes after the original references go out of scope."""
    name = bytes(b"r0")
    seq = bytes(b"ACGTACGT")
    qual = bytes(b"IIIIIIII")
    p = ReadPair(name_r1=name, seq_r1=seq, qual_r1=qual, name_r2=name, seq_r2=seq, qual_r2=qual)
    # Drop our local refs; ReadPair copies into its own owned buffers and
    # must retain them under GC pressure.
    del name, seq, qual
    gc.collect()
    assert p is not None


def test_readpair_accepts_bytearray() -> None:
    """
    ReadPair accepts a mutable bytearray and copies it.

    The input may be mutated after construction without affecting the ReadPair.
    """
    seq = bytearray(b"ACGTACGT")
    p = ReadPair(
        name_r1=bytearray(b"r0"),
        seq_r1=seq,
        name_r2=bytearray(b"r0"),
        seq_r2=bytearray(b"ACGTACGT"),
        qual_r1=bytearray(b"IIIIIIII"),
        qual_r2=bytearray(b"IIIIIIII"),
    )
    # Mutate the input after construction; ReadPair owns its copy.
    seq[0] = ord("N")
    assert p is not None


def test_readpair_accepts_memoryview() -> None:
    """ReadPair accepts memoryview over a contiguous u8 buffer."""
    backing = bytearray(b"ACGTACGTACGT")
    mv = memoryview(backing)
    p = ReadPair(
        name_r1=memoryview(bytearray(b"r0")),
        seq_r1=mv,
        name_r2=memoryview(bytearray(b"r1")),
        seq_r2=mv,
    )
    assert p is not None


def test_readpair_rejects_non_bytes_like() -> None:
    """Non bytes-like inputs raise TypeError."""
    with pytest.raises(TypeError):
        ReadPair(name_r1="r0", seq_r1=b"ACGT", name_r2=b"r1", seq_r2=b"ACGT")  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        ReadPair(name_r1=b"r0", seq_r1=12345, name_r2=b"r1", seq_r2=b"ACGT")  # type: ignore[arg-type]


def test_shm_probe_missing() -> None:
    # An almost-certainly-not-staged prefix should report False, not raise.
    assert shm.is_staged("/nonexistent/bwa-mem3-py-test-prefix") is False


def test_load_missing_index_raises() -> None:
    # Historically bwa-mem3 called `exit(1)` when `.bwt.2bit.64` was missing,
    # killing the interpreter (issue #19). The pre-flight in `BwaIndex::load`
    # now surfaces a recoverable error, so a missing prefix must raise.
    with pytest.raises(RuntimeError, match="absent"):
        bwa_mem3.BwaIndex("/nonexistent/absent")


def test_load_partial_index_raises(tmp_path: Path) -> None:
    # Every required companion present except `.ann`: the loader would
    # otherwise exit() partway through. The error must name the missing file.
    for suffix in (".bwt.2bit.64", ".amb", ".pac"):
        (tmp_path / f"partial{suffix}").write_bytes(b"stub")
    with pytest.raises(RuntimeError, match=r"\.ann"):
        bwa_mem3.BwaIndex(str(tmp_path / "partial"))
