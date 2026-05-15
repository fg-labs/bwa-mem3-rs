"""Python bindings for bwa-mem3 alignment.

Re-exports the CPython extension module ``bwa_mem3._bwa_mem3``. The
canonical entry points::

    from bwa_mem3 import BwaIndex, MemOpts, MemPeStat, ReadPair, align_batch
    from bwa_mem3 import shm

See the project README for a quickstart.
"""

from ._bwa_mem3 import (
    BwaIndex,
    MemOpts,
    MemPeStat,
    PeOrient,
    ReadPair,
    Record,
    Seeds,
    align_batch,
    estimate_pestat,
    extend_batch,
    seed_batch,
    shm,
)

__all__ = [
    "BwaIndex",
    "MemOpts",
    "MemPeStat",
    "PeOrient",
    "ReadPair",
    "Record",
    "Seeds",
    "align_batch",
    "estimate_pestat",
    "extend_batch",
    "seed_batch",
    "shm",
]
