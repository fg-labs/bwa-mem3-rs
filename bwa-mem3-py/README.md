# bwa-mem3 (Python)

Python bindings for [bwa-mem3-rs], built with [pyo3] + [maturin].

```python
from bwa_mem3 import BwaIndex, MemOpts, ReadPair, align_batch

idx = BwaIndex("/path/to/ref.fa")
opts = MemOpts()
opts.set_pe(True)

pairs = [
    ReadPair(name_r1=b"r0", seq_r1=b"ACGT...", qual_r1=b"IIII...",
             name_r2=b"r0", seq_r2=b"TGCA...", qual_r2=b"IIII..."),
]

records, pestat = align_batch(idx, opts, pairs)
for rec in records:
    # rec.bytes is one packed BAM record (BAM block payload, no BGZF).
    write_to_disk(rec.bytes)
```

## Shared memory

`bwa-mem3` can stage the index in POSIX shared memory so multiple
processes attach a single ~10 GB resident segment instead of each
loading from disk:

```python
from bwa_mem3 import BwaIndex, shm

shm.stage("/path/to/ref.fa")
# BwaIndex(...) attaches transparently when a segment is staged.
idx = BwaIndex("/path/to/ref.fa")
# ... align across many processes ...
shm.destroy()
```

For fail-fast behavior when no segment is staged, probe first:

```python
if not shm.is_staged("/path/to/ref.fa"):
    raise RuntimeError("expected staged shm segment")
idx = BwaIndex("/path/to/ref.fa")
```

## Building from source

```bash
cd bwa-mem3-py
pixi install
pixi run develop      # build & install into the env
pixi run pytest tests
```

[bwa-mem3-rs]: https://github.com/fg-labs/bwa-mem3-rs
[pyo3]: https://pyo3.rs
[maturin]: https://www.maturin.rs
