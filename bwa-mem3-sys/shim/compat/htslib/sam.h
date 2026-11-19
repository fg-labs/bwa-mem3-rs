/* Stand-in for htslib's <htslib/sam.h>, carrying only the BAM CIGAR encoding.
 *
 * bwa-mem3 v0.18.0's alt_lift.cpp (the --alt-project ALT -> primary liftover)
 * includes htslib/sam.h for the CIGAR op constants, the bam_cigar_* macros and
 * bam_cigar_table, and nothing else. This crate does not vendor htslib (the
 * refresh prunes ext/htslib), yet alt_lift.cpp cannot be left out of the build
 * the way the htslib-based writers are: mem_kernel2_core and the pairing code
 * call into it. build.rs puts shim/compat on the include path so this header
 * answers the include instead.
 *
 * It carries only the CIGAR encoding: the op codes, the bam_cigar_* macros and
 * the table alt_lift.cpp uses (plus the op codes it does not name, which the
 * table needs). If upstream starts using more of htslib there, the build fails
 * on the missing name rather than linking something half-defined; extend this
 * header only with more of the CIGAR encoding, and stop to rethink if the need
 * is anything beyond it. Delete the header (and the include path in build.rs)
 * once upstream's alt_lift.cpp stops including htslib.
 *
 * The values are the BAM specification's CIGAR encoding (SAMv1 section 4.2),
 * matching htslib's sam.h / sam.c (MIT/Expat, Genome Research Ltd.).
 */
#ifndef BWA_SHIM_COMPAT_HTSLIB_SAM_H
#define BWA_SHIM_COMPAT_HTSLIB_SAM_H

#include <stddef.h>
#include <stdint.h>

#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6
#define BAM_CEQUAL      7
#define BAM_CDIFF       8
#define BAM_CBACK       9

#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  0xf
#define BAM_CIGAR_TYPE  0x3C1A7

#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
#define bam_cigar_gen(l, o) ((l)<<BAM_CIGAR_SHIFT|(o))
/* bit 1: consumes query; bit 2: consumes reference */
#define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3)

/* htslib's bam_cigar_table: CIGAR op character -> BAM op code, -1 for anything
 * else. htslib defines it as an extern array in sam.c; here it is a constexpr
 * table in the header, so the one TU that includes this needs no extra object
 * and no link symbol. Indexing it as bam_cigar_table[(unsigned char) ch] works
 * the same either way. */
struct bwa_shim_cigar_table_t {
    int8_t v[256];
    constexpr bwa_shim_cigar_table_t() : v{} {
        for (int i = 0; i < 256; ++i) v[i] = -1;
        v['M'] = BAM_CMATCH;
        v['I'] = BAM_CINS;
        v['D'] = BAM_CDEL;
        v['N'] = BAM_CREF_SKIP;
        v['S'] = BAM_CSOFT_CLIP;
        v['H'] = BAM_CHARD_CLIP;
        v['P'] = BAM_CPAD;
        v['='] = BAM_CEQUAL;
        v['X'] = BAM_CDIFF;
        v['B'] = BAM_CBACK;
    }
    constexpr int8_t operator[](size_t i) const { return v[i]; }
};
static constexpr bwa_shim_cigar_table_t bam_cigar_table{};

static_assert(bam_cigar_table['M'] == 0 && bam_cigar_table['B'] == 9 &&
              bam_cigar_table['m'] == -1 && bam_cigar_table[0] == -1,
              "compat bam_cigar_table: SAMv1 op codes (self-consistency, not a check against htslib)");
static_assert(bam_cigar_type(BAM_CMATCH) == 3 && bam_cigar_type(BAM_CINS) == 1 &&
              bam_cigar_type(BAM_CDEL) == 2 && bam_cigar_type(BAM_CSOFT_CLIP) == 1 &&
              bam_cigar_type(BAM_CHARD_CLIP) == 0 && bam_cigar_type(BAM_CDIFF) == 3,
              "compat bam_cigar_type: SAMv1 consumes-query/reference bits per op");

#endif
