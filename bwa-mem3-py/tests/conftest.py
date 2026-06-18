"""
Shared fixtures for the bwa-mem3 Python test suite.

Mirrors the pattern used by `bwa-mem3-rs-cli/tests/common/mod.rs`: locate
`bwa-mem3` on PATH (or via the `BWA_MEM3_BIN` env var), write the embedded
PhiX174 FASTA to a tempdir, run `bwa-mem3 index`, and hand the prefix back
to dependent tests.

Two index sources are supported by `aligned_index`:

* `BWA_MEM3_RS_TEST_REF` — points at a pre-built large index (hg38, etc.);
  preferred for thorough local runs.
* The PhiX174 fixture — small enough to build in CI; gracefully skips when
  `bwa-mem3` isn't reachable.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

# pytest adds the conftest.py's directory to sys.path automatically; sibling
# fixtures import as top-level modules.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from phix_seq import PHIX_SEQ  # noqa: E402


def _find_bwa_mem3() -> str | None:
    env = os.environ.get("BWA_MEM3_BIN")
    if env:
        p = Path(env)
        if p.is_file() and os.access(p, os.X_OK):
            return str(p)
    return shutil.which("bwa-mem3")


def _existing_test_ref() -> str | None:
    p = os.environ.get("BWA_MEM3_RS_TEST_REF")
    if p and Path(f"{p}.bwt.2bit.64").exists():
        return p
    return None


@pytest.fixture(scope="session")
def bwa_mem3_bin() -> str:
    binp = _find_bwa_mem3()
    if binp is None:
        pytest.skip("bwa-mem3 not on PATH and BWA_MEM3_BIN unset")
    return binp


@pytest.fixture(scope="session")
def phix_index(tmp_path_factory: pytest.TempPathFactory, bwa_mem3_bin: str) -> str:
    """
    Build a bwa-mem3 index from the embedded PhiX174 sequence.

    Returns the prefix (the path such that `<prefix>.bwt.2bit.64` exists).
    """
    work = tmp_path_factory.mktemp("phix")
    fasta = work / "phix.fa"
    with fasta.open("w") as f:
        f.write(">phix\n")
        for i in range(0, len(PHIX_SEQ), 72):
            f.write(PHIX_SEQ[i : i + 72])
            f.write("\n")
    subprocess.run(
        [bwa_mem3_bin, "index", str(fasta)],
        check=True,
        timeout=300,
    )
    if not (work / "phix.fa.bwt.2bit.64").exists():
        pytest.skip("bwa-mem3 index did not produce a .bwt.2bit.64 file")
    return str(fasta)


@pytest.fixture(scope="session")
def aligned_index(request: pytest.FixtureRequest) -> str:
    """
    Index prefix used by integration tests.

    Resolves in priority order:
    1. `BWA_MEM3_RS_TEST_REF` (a pre-built large index); preferred locally.
    2. The session-scoped PhiX index; used in CI.

    Skips the dependent test if neither is reachable.
    """
    pre = _existing_test_ref()
    if pre is not None:
        return pre
    # Lazily request phix_index only when BWA_MEM3_RS_TEST_REF isn't set,
    # so users with a large local index don't pay the bwa-mem3 build cost.
    return request.getfixturevalue("phix_index")


def _revcomp(seq: bytes) -> bytes:
    table = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")
    return seq.translate(table)[::-1]


class _Rng:
    """Deterministic xorshift PRNG mirroring `bwa-mem3-rs-cli/tests/common::Rng`."""

    __slots__ = ("state",)

    def __init__(self, seed: int) -> None:
        self.state = seed & ((1 << 64) - 1)

    def next(self) -> int:
        s = self.state
        s ^= (s << 13) & ((1 << 64) - 1)
        s ^= s >> 7
        s ^= (s << 17) & ((1 << 64) - 1)
        self.state = s
        return s


def simulate_pairs(
    ref_seq: bytes,
    n: int,
    read_len: int = 100,
    insert: int = 250,
    seed: int = 42,
) -> list[tuple[bytes, bytes, bytes]]:
    """Simulate `n` paired reads from `ref_seq`: (name, r1, r2-revcomp)."""
    if read_len > insert:
        raise ValueError("read_len must be <= insert")
    max_start = len(ref_seq) - insert
    if max_start < 0:
        raise ValueError("reference shorter than insert size")
    rng = _Rng(seed)
    pairs: list[tuple[bytes, bytes, bytes]] = []
    for i in range(n):
        start = rng.next() % (max_start + 1) if max_start > 0 else 0
        r1 = ref_seq[start : start + read_len]
        r2 = _revcomp(ref_seq[start + insert - read_len : start + insert])
        pairs.append((f"r{i}".encode(), bytes(r1), bytes(r2)))
    return pairs
