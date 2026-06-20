"""
Type stub for the ``bwa_mem3._bwa_mem3`` compiled extension.

Hand-written to mirror the PyO3 bindings in ``src/lib.rs``. Keep the two in
sync: when a class, method, or function signature changes there, update it
here. Shipped in the wheel (alongside ``py.typed``) so downstream callers and
IDEs get type information for the native API.
"""

import os
from collections.abc import Sequence
from typing import Final

# Any object copied into Rust-owned storage at construction. The native code
# also accepts arbitrary buffer-protocol objects (e.g. ``array.array('B')``,
# numpy ``uint8`` arrays); the common cases are spelled out here.
BytesLike = bytes | bytearray | memoryview
StrPath = str | os.PathLike[str]

class BwaIndex:
    """Reference index handle. ``BwaIndex(prefix)`` loads from disk."""

    def __init__(self, prefix: StrPath) -> None: ...
    def n_contigs(self) -> int: ...
    def contigs(self) -> list[tuple[str, int]]: ...
    def __repr__(self) -> str: ...

class MemOpts:
    """Bwa-mem3 alignment options, constructed with bwa-mem3 defaults."""

    def __init__(self) -> None: ...
    def apply_mode(self, mode: str) -> None: ...
    def set_pe(self, is_pe: bool) -> None: ...
    def set_soft_clip_supplementary(self, v: bool) -> None: ...
    min_seed_len: int
    band_width: int
    match_score: int
    mismatch_penalty: int
    minimum_score: int
    max_occurrences: int
    # `del`/`ins` are positional-only (`del` is not a valid Python keyword arg).
    def set_gap_open(self, del_: int, ins: int, /) -> None: ...
    def set_gap_extend(self, del_: int, ins: int, /) -> None: ...
    def set_clip_penalty(self, five: int, three: int) -> None: ...
    def set_xa_max_hits(self, primary: int, alt: int) -> None: ...
    def set_xa_drop_ratio(self, v: float) -> None: ...
    def set_unpaired_penalty(self, v: int) -> None: ...
    def set_read_group_id(self, id: str | None = None) -> None: ...

class PeOrient:
    """Insert-size statistics for a single orientation."""

    low: int
    high: int
    failed: bool
    avg: float
    std: float
    def __init__(
        self,
        low: int = 0,
        high: int = 0,
        failed: bool = False,
        avg: float = 0.0,
        std: float = 0.0,
    ) -> None: ...
    def __repr__(self) -> str: ...

class MemPeStat:
    """4-orientation paired-end insert-size model (`mem_pestat_t[4]`)."""

    def __init__(self) -> None: ...
    def orientation(self, o: str) -> PeOrient: ...
    def set_orientation(self, o: str, v: PeOrient) -> None: ...

class ReadPair:
    """One paired-end read; owns copies of its name/seq/qual bytes."""

    def __init__(
        self,
        name_r1: BytesLike,
        seq_r1: BytesLike,
        name_r2: BytesLike,
        seq_r2: BytesLike,
        qual_r1: BytesLike | None = None,
        qual_r2: BytesLike | None = None,
    ) -> None: ...

class Record:
    """One packed BAM record produced by alignment."""

    pair_idx: int
    bytes: bytes
    def __repr__(self) -> str: ...

class Seeds:
    """Opaque seeds handle from `seed_batch`; consumed once by `extend_batch`."""

    def __repr__(self) -> str: ...

def align_batch(
    index: BwaIndex,
    opts: MemOpts,
    pairs: Sequence[ReadPair],
    pestat_in: MemPeStat | None = None,
) -> tuple[list[Record], MemPeStat]: ...
def seed_batch(index: BwaIndex, opts: MemOpts, pairs: Sequence[ReadPair]) -> Seeds: ...
def extend_batch(
    index: BwaIndex,
    opts: MemOpts,
    seeds: Seeds,
    pairs: Sequence[ReadPair],
    pestat_in: MemPeStat | None = None,
) -> tuple[list[Record], MemPeStat]: ...
def estimate_pestat(index: BwaIndex, opts: MemOpts, pairs: Sequence[ReadPair]) -> MemPeStat: ...

class _Shm:
    """The native `shm` submodule: shared-memory index staging helpers."""

    @staticmethod
    def is_staged(prefix: StrPath) -> bool: ...
    @staticmethod
    def stage(prefix: StrPath) -> None: ...
    @staticmethod
    def destroy() -> None: ...
    @staticmethod
    def list() -> None: ...

shm: Final[_Shm]
