/* bwa-mem3-sys/shim/bwa_shim_types.h
 *
 * POD struct definitions copied from bwa-mem3's `bwamem.h` (mem_opt_t,
 * mem_pestat_t) so bindgen can emit full Rust struct layouts without
 * having to parse the C++ header graph rooted at bwamem.h.
 *
 * Keep the layouts byte-identical to bwa-mem3's. When refreshing the
 * vendored snapshot, diff against `vendor/bwa-mem3/src/bwamem.h` around
 * lines 95-156 (mem_opt_t) and 239-243 (mem_pestat_t); if either changed,
 * update here. The shim allocates the real struct via upstream's
 * `mem_opt_init()`, so any field-order/size drift silently corrupts the
 * offsets Rust reads/writes through the bindgen view of this header.
 */
#ifndef BWA_SHIM_TYPES_H
#define BWA_SHIM_TYPES_H

#include <stddef.h>
#include <stdint.h>

/* MEM_F_* flag bits (from bwamem.h). */
#define MEM_F_PE             0x2
#define MEM_F_NOPAIRING      0x4
#define MEM_F_ALL            0x8
#define MEM_F_NO_MULTI       0x10
#define MEM_F_NO_RESCUE      0x20
#define MEM_F_REF_HDR        0x100
#define MEM_F_SOFTCLIP       0x200
#define MEM_F_SMARTPE        0x400
#define MEM_F_PRIMARY5       0x800
#define MEM_F_KEEP_SUPP_MAPQ 0x1000
#define MEM_F_XB             0x2000

/* Mirror of bwamem.h:95-156. Layout must match exactly. */
typedef struct mem_opt_t {
    int a, b;
    int o_del, e_del;
    int o_ins, e_ins;
    int pen_unpaired;
    int pen_clip5, pen_clip3;
    int w;
    int zdrop;

    uint64_t max_mem_intv;

    int T;
    int flag;
    int min_seed_len;
    int min_ext_len;
    int max_extend_chains;
    int mate_concordant_window;
    int est_insert_high;
    /* upstream type is `seed_order_t`, a plain (int-sized) enum in bwamem.h;
     * mirrored as int for layout. The Rust API exposes it as the `SeedOrder`
     * enum (bwa-mem3-rs `MemOpts::seed_order`). */
    int seed_emit_order;
    int min_chain_weight;
    int max_chain_extend;
    float split_factor;
    int split_width;
    int max_occ;
    int max_chain_gap;
    int n_threads;
    int64_t chunk_size;
    float mask_level;
    float drop_ratio;
    float XA_drop_ratio;
    float mask_level_redun;
    float mapQ_coef_len;
    int mapQ_coef_fac;
    int max_ins;
    int max_matesw;
    /* v0.9.0 */
    int rescue_kmer;
    int rescue_band;
    int rescue_skip;
    int max_XA_hits, max_XA_hits_alt;
    int8_t mat[25];
    int8_t mat_ot[25];
    int8_t mat_ob[25];
    int    bam_mode;
    int    bam_level;
    int    meth_mode;
    int    meth_scoring;
    /* v0.9.0 */
    int    meth_chem;
    int    meth_tags;
    char   meth_set_as_failed;
    int    meth_chimera_qc;
    /* v0.9.0 */
    int    proper_pair_from_emitted;
    int    supp_rep_hard_cap;
    int    smem_dedup;
    /* v0.9.0 */
    int    alnreg_sort_fast;
    int    skip_contained_ext;
    int    band_start;
    /* v0.9.0. Upstream's type is `const compat_target_t *`; mirrored as an
     * opaque `const void *` because the POD only has to reproduce the LAYOUT
     * (pointer size and alignment), and pulling in compat_target_t would drag
     * the header graph bwa_shim.h exists to avoid. Nothing on the Rust side
     * reads it; bwa_shim_layout_assert.cpp checks the offset against the real
     * struct, which is what keeps this honest. */
    const void *compat;
} mem_opt_t;

/* Mirror of bwamem.h:239-243. */
typedef struct mem_pestat_t {
    int low, high;
    int failed;
    double avg, std;
} mem_pestat_t;

#endif /* BWA_SHIM_TYPES_H */
