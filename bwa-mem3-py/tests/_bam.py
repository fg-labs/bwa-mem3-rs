"""
Minimal BAM record decoder for test assertions.

`Record.bytes` is laid out as `[u32 le block_size][body]` (the BAM block
format; see `bwa-mem3-sys/shim/bwa_shim_align.cpp:39`). This module decodes
just enough of the body to assert on FLAG, MAPQ, CIGAR length, and aux tags
without pulling in pysam.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass

# BAM FLAG bits (subset).
PAIRED = 0x1
PROPER_PAIR = 0x2
UNMAPPED = 0x4
MATE_UNMAPPED = 0x8
REVERSE = 0x10
MATE_REVERSE = 0x20
READ1 = 0x40
READ2 = 0x80
SECONDARY = 0x100
SUPPLEMENTARY = 0x800


# BAM aux tag value sizes for fixed-width types.
_AUX_FIXED_SIZE: dict[bytes, int] = {
    b"A": 1,
    b"c": 1,
    b"C": 1,
    b"s": 2,
    b"S": 2,
    b"i": 4,
    b"I": 4,
    b"f": 4,
}


@dataclass(frozen=True)
class BamRecord:
    flag: int
    mapq: int
    n_cigar_op: int
    l_seq: int
    read_name: bytes
    aux: dict[bytes, tuple[bytes, bytes]]
    """Aux tags as {tag (2 bytes): (type (1 byte), raw value bytes)}."""

    def aux_z(self, tag: bytes) -> str | None:
        """Return the value of a `Z`-typed aux tag as a UTF-8 string, or None."""
        v = self.aux.get(tag)
        if v is None or v[0] != b"Z":
            return None
        # Trim the terminating NUL.
        return v[1].rstrip(b"\x00").decode("utf-8", errors="replace")

    @property
    def is_unmapped(self) -> bool:
        return bool(self.flag & UNMAPPED)


def decode(record_bytes: bytes) -> BamRecord:
    """Decode one packed BAM record (with the leading u32 block_size)."""
    if len(record_bytes) < 4:
        raise ValueError("record shorter than block_size prefix")
    (block_size,) = struct.unpack_from("<I", record_bytes, 0)
    if 4 + block_size != len(record_bytes):
        raise ValueError(
            f"block_size mismatch: header says {block_size}, body has {len(record_bytes) - 4} bytes"
        )
    body = record_bytes[4:]
    if len(body) < 32:
        raise ValueError("body shorter than fixed BAM header (32 bytes)")
    l_read_name = body[8]
    mapq = body[9]
    n_cigar_op = struct.unpack_from("<H", body, 12)[0]
    flag = struct.unpack_from("<H", body, 14)[0]
    l_seq = struct.unpack_from("<i", body, 16)[0]

    off = 32
    read_name = body[off : off + l_read_name].rstrip(b"\x00")
    off += l_read_name
    off += 4 * n_cigar_op  # CIGAR
    off += (l_seq + 1) // 2  # seq (4-bit packed)
    off += l_seq  # qual

    aux: dict[bytes, tuple[bytes, bytes]] = {}
    while off < len(body):
        if off + 3 > len(body):
            raise ValueError("truncated aux tag header")
        tag = body[off : off + 2]
        typ = body[off + 2 : off + 3]
        off += 3
        if typ in _AUX_FIXED_SIZE:
            n = _AUX_FIXED_SIZE[typ]
            value = body[off : off + n]
            off += n
        elif typ == b"Z" or typ == b"H":
            end = body.index(b"\x00", off) + 1
            value = body[off:end]
            off = end
        elif typ == b"B":
            sub = body[off : off + 1]
            count = struct.unpack_from("<I", body, off + 1)[0]
            sub_sz = _AUX_FIXED_SIZE.get(sub, 0)
            if sub_sz == 0:
                raise ValueError(f"unsupported B-array subtype {sub!r}")
            end = off + 5 + sub_sz * count
            value = body[off:end]
            off = end
        else:
            raise ValueError(f"unsupported aux type {typ!r}")
        aux[tag] = (typ, value)

    return BamRecord(
        flag=flag,
        mapq=mapq,
        n_cigar_op=n_cigar_op,
        l_seq=l_seq,
        read_name=read_name,
        aux=aux,
    )
