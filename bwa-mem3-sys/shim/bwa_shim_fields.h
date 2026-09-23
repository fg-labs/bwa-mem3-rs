/*
 * Structured-fields record sink for the resident-cohort pair-emit calls.
 *
 * Plain C POD, stdint-only, so it can be included both by the public header
 * (bwa_shim.h, for bindgen) and by bwa_shim_align.cpp, which cannot include
 * bwa_shim.h (its POD mem_opt_t collides with upstream's real one). One
 * definition means the two translation units cannot drift apart.
 */
#ifndef BWA_SHIM_FIELDS_H
#define BWA_SHIM_FIELDS_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Everything the packed-BAM record path serializes for one emitted record,
 * as values instead of bytes. The packed path is literally "compute these,
 * then serialize them", so a record built from these fields is byte-identical
 * to the one the packed sink would have received.
 *
 * Every pointer is borrowed from shim scratch (or from bwa's own per-record
 * state) and is valid only for the duration of the sink call -- the same
 * contract as the packed sink's `body`. A NULL string / zero `has_*` means the
 * tag is absent from the packed record. Strings are NUL-terminated. */
typedef struct BwaAlignedFields {
    /* Fixed-width BAM fields, exactly as serialized: tid/pos are the EFFECTIVE
     * placement (a half-mapped pair's unmapped read carries its mate's), and
     * `flag` is the final 16-bit FLAG, whose 0x10 bit is the emitted strand. */
    int32_t  tid;
    int32_t  pos;
    int32_t  next_tid;
    int32_t  next_pos;
    int32_t  tlen;
    uint16_t flag;
    uint16_t bin;
    uint8_t  mapq;

    /* Emitted SEQ/QUAL as a window of the read: bases [query_start, query_end)
     * of the read in as-sequenced orientation, reverse-complemented (QUAL
     * reversed) when flag & 0x10. The window is shorter than the read on a
     * hard-clipped supplementary and empty on a raw-0x100 (-a) secondary.
     * Without --meth the packed SEQ passes through bwa's 2-bit alphabet, so
     * lowercase is folded to uppercase and every base other than A/C/G/T
     * (IUPAC codes included) becomes N. Under --meth the packed path emits
     * the ORIGINAL (unprojected) input bases: forward bases upper-cased with
     * IUPAC codes kept, reverse bases complemented with every base other
     * than A/C/G/T becoming N. */
    int32_t  query_start;
    int32_t  query_end;

    /* CIGAR in BAM opcodes (MIDNSHP=X), with the hard/soft clip rewrite
     * already applied. */
    const uint32_t *cigar;
    uint32_t        n_cigar;

    /* Aux tags, listed in the packed record's write order. */
    uint8_t     has_nm;  int32_t nm;
    const char *md;
    const char *mc;
    uint8_t     has_mq;  int32_t mq;
    uint8_t     has_as;  int32_t score;   /* AS:i */
    uint8_t     has_xs;  int32_t sub;     /* XS:i */
    const char *rg;
    const char *sa;
    uint8_t     has_pa;  float pa;        /* pa:f */
    const char *xa;
    uint8_t     has_hn;  int32_t hn;
    const char *xr;                       /* --meth only */
    const char *xg;                       /* --meth only */
    const char *xm;                       /* --meth only */
} BwaAlignedFields;

/* Structured sink: called once per emitted record, in the same order the
 * packed sink would be. `mate` is 0 for R1 or a single and 1 for R2 (which of
 * a pair's two reads this record belongs to); `is_primary` is 1 for the first
 * record emitted for that read and 0 for each supplementary (or, under -a,
 * secondary) after it. `f` is valid only for the duration of the call. */
typedef void (*BwaFieldSinkFn)(void *ctx, uint32_t origin_kind, size_t origin_idx,
                               uint8_t mate, int is_primary, const BwaAlignedFields *f);

#ifdef __cplusplus
}
#endif

#endif /* BWA_SHIM_FIELDS_H */
