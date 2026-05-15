/* bwa-mem3-sys/shim/bwa_shim.cpp
 *
 * Shim implementation for the Rust FFI crate.
 *
 * STATUS: index lifecycle + error TLS + verbosity are wired to real
 * bwa-mem3 public API. The alignment entry points (seed_batch,
 * extend_batch, align_batch, estimate_pestat) are stubs that return
 * NULL/non-zero with a "not yet implemented" error. Filling them in
 * requires replicating the init sequence from fastmap.cpp's
 * `memoryAlloc` and `mem_process_seqs` in a reentrant, single-threaded
 * form, and then emitting packed BAM bytes from `mem_aln_t` values
 * directly.
 *
 * See CLAUDE.md for the intended shim semantics + known gotchas.
 */

#include "bwa_shim.h"

/* Include upstream's real bwamem.h to static-assert our copied mem_opt_t
 * layout still matches. `bwa_shim_types.h` defines mem_opt_t first (via
 * bwa_shim.h), and bwamem.h redefines it — but C++ forbids redefinition of
 * a struct with different body, so we re-namespace upstream's into a nested
 * struct via preprocessor trickery before inclusion.
 *
 * Simpler: extract the field offset checks into a separate TU. For now,
 * rely on build-time review when refreshing vendor/ and document in the
 * types header. (Layout verification is a follow-up task.)
 */
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Forward declarations from vendored bwa-mem3 public headers. We avoid
 * including the full C++ headers here because bwa_shim.cpp is compiled as
 * C++; including "bwa.h" directly pulls in the whole graph. Instead we
 * declare only the symbols we actually call. */

/* bwa-mem3 internals we reference from this TU. All wrapped functions live
 * in bwa_shim_align.cpp which includes upstream's real headers; here we only
 * need the shim_align_* bridge declarations. */
mem_opt_t *mem_opt_init(void);  /* C++ linkage, matches bwamem.h */

extern "C" {
    struct ShimReadPair {
        const char    *r1_name;  size_t r1_name_len;
        const uint8_t *r1_seq;   size_t r1_seq_len;
        const uint8_t *r1_qual;
        const char    *r2_name;  size_t r2_name_len;
        const uint8_t *r2_seq;   size_t r2_seq_len;
        const uint8_t *r2_qual;
    };
    struct ShimAlignOutput;

    void *shim_align_idx_load(const char *prefix);
    void  shim_align_idx_free(void *fmi);
    size_t shim_align_idx_n_contigs(void *fmi);
    const char *shim_align_idx_contig_name(void *fmi, size_t i);
    int64_t shim_align_idx_contig_len(void *fmi, size_t i);

    struct ShimSeeds;

    ShimSeeds       *shim_seed_batch(
        void *fmi, mem_opt_t *opts,
        const ShimReadPair *pairs, size_t n_pairs);
    void             shim_seeds_free(ShimSeeds *seeds);

    ShimAlignOutput *shim_extend_batch(
        void *fmi, ShimSeeds *seeds,
        const mem_pestat_t *pestat_in);

    ShimAlignOutput *shim_align_batch(
        void *fmi, mem_opt_t *opts,
        const ShimReadPair *pairs, size_t n_pairs,
        const mem_pestat_t *pestat_in);

    int shim_estimate_pestat(
        void *fmi, mem_opt_t *opts,
        const ShimReadPair *pairs, size_t n_pairs,
        mem_pestat_t *pestat_out);

    size_t         shim_align_out_n_recs(ShimAlignOutput *out);
    size_t         shim_align_out_pair_idx(ShimAlignOutput *out, size_t i);
    const uint8_t *shim_align_out_rec_ptr(ShimAlignOutput *out, size_t i);
    size_t         shim_align_out_rec_len(ShimAlignOutput *out, size_t i);
    void           shim_align_out_free(ShimAlignOutput *out);
    void           shim_align_get_pestat(ShimAlignOutput *out, mem_pestat_t *dst);
}

/* From bwa-mem3's bntseq.h, exposed via bwaidx_t */
struct bntann1_t { int64_t offset; int32_t len; int32_t n_ambs; uint32_t gi;
                   int32_t is_alt; char *name; char *anno; };
struct bntseq_t { int64_t l_pac; int32_t n_seqs; uint32_t seed;
                  struct bntann1_t *anns; /* ... */ };

struct bwaidx_s { struct bntseq_t *bns; /* ... */ };

/* Global verbosity knob in bwa-mem3 */
extern "C" int bwa_verbose;
/* 256-byte buffer; bwa-mem3 reads `bwa_rg_id[0]` truthiness before emitting RG:Z. */
extern "C" char bwa_rg_id[256];

/* -------- thread-local last-error --------------------------------- */

#if defined(_MSC_VER)
#  define BWA_SHIM_TLS __declspec(thread)
#else
#  define BWA_SHIM_TLS __thread
#endif

static BWA_SHIM_TLS char g_last_err[512];

static void shim_set_err(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(g_last_err, sizeof(g_last_err), fmt, ap);
    va_end(ap);
}
static void shim_clear_err(void) { g_last_err[0] = '\0'; }

extern "C" const char *bwa_shim_last_error(void) {
    return g_last_err[0] ? g_last_err : NULL;
}

/* -------- verbosity ----------------------------------------------- */

extern "C" void bwa_shim_set_verbosity(int level) { bwa_verbose = level; }

extern "C" void bwa_shim_set_rg_id(const char *id) {
    if (!id) { bwa_rg_id[0] = '\0'; return; }
    size_t n = strnlen(id, sizeof(bwa_rg_id) - 1);
    memcpy(bwa_rg_id, id, n);
    bwa_rg_id[n] = '\0';
}

/* -------- shared memory ------------------------------------------- */
/* Forward-declare upstream's bwa_shm_* entry points so we don't have to
 * pull bwa_shm.h into this TU. They're plain C linkage in upstream. */
extern "C" int bwa_shm_test(const char *prefix);
extern "C" int bwa_shm_stage(const char *prefix);
extern "C" int bwa_shm_destroy(void);
extern "C" int bwa_shm_list(void);

extern "C" int bwa_shim_shm_test(const char *prefix) {
    shim_clear_err();
    if (!prefix) {
        shim_set_err("null prefix");
        return -1;
    }
    int rc = bwa_shm_test(prefix);
    if (rc < 0) shim_set_err("bwa_shm_test failed (registry inaccessible?)");
    return rc;
}

extern "C" int bwa_shim_shm_stage(const char *prefix) {
    shim_clear_err();
    if (!prefix) {
        shim_set_err("null prefix");
        return -1;
    }
    int rc = bwa_shm_stage(prefix);
    if (rc < 0) shim_set_err("bwa_shm_stage failed for '%s'", prefix);
    return rc;
}

extern "C" int bwa_shim_shm_destroy(void) {
    shim_clear_err();
    int rc = bwa_shm_destroy();
    if (rc < 0) shim_set_err("bwa_shm_destroy failed");
    return rc;
}

extern "C" int bwa_shim_shm_list(void) {
    shim_clear_err();
    int rc = bwa_shm_list();
    if (rc < 0) shim_set_err("bwa_shm_list failed");
    return rc;
}

/* -------- options + pestat lifecycle ------------------------------ */
/*
 * mem_opt_t and mem_pestat_t are kept opaque to Rust for now. The shim
 * constructs them via `mem_opt_init` (bwa-mem3 public API) and frees them.
 * A future revision will add a second bindgen pass on bwamem.h to expose
 * the full struct so Rust can mutate fields directly.
 */

extern "C" mem_opt_t *bwa_shim_opts_new(void) {
    shim_clear_err();
    mem_opt_t *o = mem_opt_init();
    if (!o) shim_set_err("mem_opt_init returned NULL");
    return o;
}

extern "C" void bwa_shim_opts_free(mem_opt_t *opts) {
    if (opts) free(opts);
}

extern "C" int bwa_shim_opts_set_int(mem_opt_t *opts, const char *key, int value) {
    (void)opts; (void)key; (void)value;
    shim_set_err("bwa_shim_opts_set_int: deprecated; set fields directly via Rust");
    return -1;
}

extern "C" mem_pestat_t *bwa_shim_pestat_zero(void) {
    shim_clear_err();
    /* 4-orientation insert-size model. */
    mem_pestat_t *p = (mem_pestat_t *)calloc(4, sizeof(mem_pestat_t));
    if (!p) shim_set_err("pestat_zero alloc failed");
    return p;
}

extern "C" void bwa_shim_pestat_free(mem_pestat_t *p) {
    if (p) free(p);
}

/* -------- index --------------------------------------------------- */

struct BwaIndex {
    void *fmi;  /* FMI_search* (upstream class); opaque to this TU */
};

extern "C" BwaIndex *bwa_shim_idx_load(const char *prefix) {
    shim_clear_err();
    if (!prefix) {
        shim_set_err("null prefix");
        return NULL;
    }
    BwaIndex *h = (BwaIndex *) calloc(1, sizeof(BwaIndex));
    if (!h) {
        shim_set_err("calloc failed");
        return NULL;
    }
    h->fmi = shim_align_idx_load(prefix);
    if (!h->fmi) {
        shim_set_err("FMI_search load failed for '%s'", prefix);
        free(h);
        return NULL;
    }
    return h;
}

extern "C" void bwa_shim_idx_free(BwaIndex *h) {
    if (!h) return;
    if (h->fmi) {
        shim_align_idx_free(h->fmi);
        h->fmi = NULL;
    }
    free(h);
}

extern "C" size_t bwa_shim_idx_n_contigs(const BwaIndex *h) {
    return h ? shim_align_idx_n_contigs(h->fmi) : 0;
}

extern "C" const char *bwa_shim_idx_contig_name(const BwaIndex *h, size_t i) {
    return shim_align_idx_contig_name(h->fmi, i);
}

extern "C" int64_t bwa_shim_idx_contig_len(const BwaIndex *h, size_t i) {
    return shim_align_idx_contig_len(h->fmi, i);
}

/* -------- seeds / batches ----------------------------------------- */

/* BwaBatch wraps the ShimAlignOutput produced by shim_align_batch.
 * The `bytes` stored inside are SAM lines (temporary; see STATUS doc);
 * full packed-BAM emission is a follow-up. */
struct BwaBatch {
    struct ShimAlignOutput *inner;
};

/* BwaSeeds wraps the ShimSeeds from the align bridge. */
struct BwaSeeds {
    struct ShimSeeds *inner;
};

extern "C" BwaSeeds *bwa_shim_seed_batch(
    const BwaIndex *idx, const mem_opt_t *opts,
    const BwaReadPair *pairs, size_t n_pairs)
{
    shim_clear_err();
    if (!idx || !opts || (n_pairs > 0 && !pairs)) {
        shim_set_err("null arg");
        return NULL;
    }
    ShimSeeds *inner = shim_seed_batch(
        idx->fmi, const_cast<mem_opt_t *>(opts),
        reinterpret_cast<const ShimReadPair *>(pairs), n_pairs);
    if (!inner) {
        if (!bwa_shim_last_error()) shim_set_err("seed_batch failed");
        return NULL;
    }
    BwaSeeds *s = (BwaSeeds *) calloc(1, sizeof(BwaSeeds));
    s->inner = inner;
    return s;
}

extern "C" void bwa_shim_seeds_free(BwaSeeds *seeds) {
    if (!seeds) return;
    if (seeds->inner) shim_seeds_free(seeds->inner);
    free(seeds);
}

extern "C" BwaBatch *bwa_shim_extend_batch(
    const BwaIndex *idx, const mem_opt_t *opts,
    BwaSeeds *seeds,
    const BwaReadPair *pairs, size_t n_pairs,
    const mem_pestat_t *pestat_in,
    mem_pestat_t *pestat_out)
{
    shim_clear_err();
    (void)opts; (void)pairs; (void)n_pairs;
    if (!idx || !seeds || !pestat_out) {
        if (seeds) bwa_shim_seeds_free(seeds);
        shim_set_err("null arg");
        return NULL;
    }
    /* Transfer the ShimSeeds ownership into the shim bridge. */
    ShimSeeds *inner_seeds = seeds->inner;
    seeds->inner = nullptr;
    free(seeds);

    ShimAlignOutput *inner = shim_extend_batch(idx->fmi, inner_seeds, pestat_in);
    if (!inner) {
        if (!bwa_shim_last_error()) shim_set_err("extend_batch failed");
        return NULL;
    }
    shim_align_get_pestat(inner, pestat_out);
    BwaBatch *b = (BwaBatch *) calloc(1, sizeof(BwaBatch));
    b->inner = inner;
    return b;
}

extern "C" BwaBatch *bwa_shim_align_batch(
    const BwaIndex *idx, const mem_opt_t *opts,
    const BwaReadPair *pairs, size_t n_pairs,
    const mem_pestat_t *pestat_in,
    mem_pestat_t *pestat_out)
{
    shim_clear_err();
    if (!idx || !opts || (n_pairs > 0 && !pairs) || !pestat_out) {
        shim_set_err("null arg");
        return NULL;
    }
    /* Cast our BwaReadPair (from bwa_shim.h) to the bridge's ShimReadPair.
     * Both have identical layout (same fields in same order). */
    ShimAlignOutput *inner = shim_align_batch(
        idx->fmi, const_cast<mem_opt_t *>(opts),
        reinterpret_cast<const ShimReadPair *>(pairs), n_pairs,
        pestat_in);
    if (!inner) {
        if (!bwa_shim_last_error()) shim_set_err("align_batch failed");
        return NULL;
    }
    shim_align_get_pestat(inner, pestat_out);
    BwaBatch *b = (BwaBatch *) calloc(1, sizeof(BwaBatch));
    b->inner = inner;
    return b;
}

extern "C" int bwa_shim_estimate_pestat(
    const BwaIndex *idx, const mem_opt_t *opts,
    const BwaReadPair *pairs, size_t n_pairs,
    mem_pestat_t *pestat_out)
{
    shim_clear_err();
    if (!idx || !opts || !pestat_out) {
        shim_set_err("null arg");
        return -1;
    }
    return shim_estimate_pestat(
        idx->fmi, const_cast<mem_opt_t *>(opts),
        reinterpret_cast<const ShimReadPair *>(pairs), n_pairs, pestat_out);
}

extern "C" size_t bwa_shim_batch_n_records(const BwaBatch *b) {
    return b ? shim_align_out_n_recs(b->inner) : 0;
}
extern "C" size_t bwa_shim_batch_pair_idx(const BwaBatch *b, size_t rec) {
    return shim_align_out_pair_idx(b->inner, rec);
}
extern "C" const uint8_t *bwa_shim_batch_record_ptr(const BwaBatch *b, size_t rec) {
    return shim_align_out_rec_ptr(b->inner, rec);
}
extern "C" size_t bwa_shim_batch_record_len(const BwaBatch *b, size_t rec) {
    return shim_align_out_rec_len(b->inner, rec);
}
extern "C" void bwa_shim_batch_free(BwaBatch *b) {
    if (!b) return;
    shim_align_out_free(b->inner);
    free(b);
}
