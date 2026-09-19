/* bwa-mem3-sys/shim/bwa_shim.h
 *
 * Public C API for the bwa-mem3 Rust FFI crate. Consumers of this header
 * should only see:
 *   - opaque handle types (BwaIndex, BwaSeeds, BwaBatch)
 *   - POD input types (BwaReadPair)
 *   - function prototypes
 *
 * All mem_opt_t / mem_pestat_t references are forward-declared as opaque
 * so bindgen does not have to parse the full bwa-mem3 C++ header graph
 * for the shim header alone. The Rust side gets mem_opt_t / mem_pestat_t
 * bindings from a separate bindgen pass on bwamem.h.
 */
#ifndef BWA_SHIM_H
#define BWA_SHIM_H

#include <stddef.h>
#include <stdint.h>

#include "bwa_shim_types.h"  /* mem_opt_t, mem_pestat_t POD layouts */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct BwaIndex BwaIndex;
typedef struct BwaSeeds BwaSeeds;
typedef struct BwaBatch BwaBatch;
typedef struct BwaScratch BwaScratch;
typedef struct BwaRegs    BwaRegs;

typedef struct {
    const char    *r1_name;  size_t r1_name_len;
    const uint8_t *r1_seq;   size_t r1_seq_len;
    const uint8_t *r1_qual;
    const char    *r2_name;  size_t r2_name_len;
    const uint8_t *r2_seq;   size_t r2_seq_len;
    const uint8_t *r2_qual;
} BwaReadPair;

/* A single-end read. Layout mirrors the bridge's ShimSingleRead. */
typedef struct {
    const char    *name;  size_t name_len;
    const uint8_t *seq;   size_t seq_len;
    const uint8_t *qual;
} BwaSingleRead;

/* One work item: pairs (r1/r2 interleaved into seqs[0..2*n_pairs)) followed by
 * singles (seqs[2*n_pairs..)). Either slice may be empty. */
typedef struct {
    const BwaReadPair   *pairs;   size_t n_pairs;
    const BwaSingleRead *singles; size_t n_singles;
} BwaReadBatch;

/* Options lifecycle. `opts_new` returns a heap-allocated mem_opt_t populated
 * with bwa-mem3 defaults (mem_opt_init). `opts_free` releases it. Field-level
 * getters/setters will be added once the full mem_opt_t is bound to Rust. */
mem_opt_t *bwa_shim_opts_new(void);
void              bwa_shim_opts_free(mem_opt_t *opts);

/* D3 (--meth): (re)build the per-hypothesis bisulfite scoring matrices from the
 * base matrix per opts->meth_scoring. Call after changing meth_scoring. */
void              bwa_shim_opts_fill_meth_mat(mem_opt_t *opts);

/* Rebuild the 5x5 scoring matrix from opts->a / opts->b, then refresh the
 * bisulfite matrices derived from it. Call after changing either score.
 *
 * Writing `a`/`b` alone leaves the crate incoherent: several code paths read
 * them directly (bwamem.cpp:2200, 2296, 3449, 4130, ...) while the SW kernels
 * read only opts->mat, which still encodes the previous values. Upstream
 * always pairs the two (bwa_fill_scmat then mem_opt_fill_meth_mat --
 * fastmap.cpp:2726 and :2730, in main_mem), and mem_opt_fill_meth_mat's own
 * contract comment requires it after EVERY rebuild of opts->mat. */
void              bwa_shim_opts_fill_scmat(mem_opt_t *opts);

/* Set common single-integer fields; one function per semantically-distinct knob.
 * Returns 0 on success, non-zero if the key is unknown. */
int bwa_shim_opts_set_int(mem_opt_t *opts, const char *key, int value);

/* The @HD header line the active output-compatibility target calls for, without
 * a trailing newline, or NULL when the target emits no @HD.
 *
 * A caller writing its own SAM/BAM header must use this rather than a literal
 * of its own: `compat_target_t` is the one place that decides output shaping,
 * and upstream collapsed three drifted @HD literals into it for exactly that
 * reason (fg-labs/bwa-mem3#288). The string is static storage owned by the
 * compat table; do not free it. */
const char *bwa_shim_compat_hd_line(const mem_opt_t *opts);

/* Apply bwa-mem3's bwameth-compatibility defaults for `--meth`, then refill the
 * scoring matrices (the bundle can change `b`, so they would otherwise be
 * stale).
 *
 * Wraps upstream's `mem_opt_apply_meth_defaults`, so the constants scale with
 * the match score `a` exactly as upstream scales them.
 *
 * Ordering is not symmetric across the knobs involved:
 *   - `a` (`-A`) is an INPUT -- the constants are expressed in units of it --
 *     so set it BEFORE calling. So is `meth_scoring`: the `-B` branch keys off
 *     the resolved mode.
 *   - `T`/`pen_clip5`/`pen_clip3`/`pen_unpaired` (`-T`/`-L`/`-U`) are
 *     OVERWRITTEN unconditionally, because upstream's "user set this" mask is
 *     passed empty, so set any of those AFTER or they are silently clobbered.
 *   - `b` (`-B`) is overwritten ONLY under COLLAPSED scoring. GENOMIC and
 *     NEUTRAL keep bwa's default (bwamem.cpp:511-515), so a caller-set `-B`
 *     survives the bundle under those two modes. */
void bwa_shim_opts_apply_meth_defaults(mem_opt_t *opts);

/* PE-stats lifecycle. `pestat_zero` returns a zeroed 4-orientation array. */
mem_pestat_t *bwa_shim_pestat_zero(void);
void                 bwa_shim_pestat_free(mem_pestat_t *pestat);

BwaIndex *bwa_shim_idx_load(const char *prefix);
/* D3 (--meth): load a dual index — `seed_prefix` = converted `<ref>.meth`,
 * `orig_prefix` = un-converted `<ref>`. Use with meth_mode set on the opts. */
BwaIndex *bwa_shim_idx_load_meth(const char *seed_prefix, const char *orig_prefix);
/* Non-zero iff `idx` is a --meth dual index (loaded via bwa_shim_idx_load_meth). */
int       bwa_shim_idx_is_meth(const BwaIndex *idx);
void      bwa_shim_idx_free(BwaIndex *idx);
size_t    bwa_shim_idx_n_contigs(const BwaIndex *idx);
const char *bwa_shim_idx_contig_name(const BwaIndex *idx, size_t i);
int64_t   bwa_shim_idx_contig_len(const BwaIndex *idx, size_t i);

BwaSeeds *bwa_shim_seed_batch(
    const BwaIndex *idx, const mem_opt_t *opts,
    const BwaReadPair *pairs, size_t n_pairs);
void bwa_shim_seeds_free(BwaSeeds *seeds);

BwaBatch *bwa_shim_extend_batch(
    const BwaIndex *idx, const mem_opt_t *opts,
    BwaSeeds *seeds,
    const BwaReadPair *pairs, size_t n_pairs,
    const mem_pestat_t *pestat_in,
    mem_pestat_t *pestat_out);

BwaBatch *bwa_shim_align_batch(
    const BwaIndex *idx, const mem_opt_t *opts,
    const BwaReadPair *pairs, size_t n_pairs,
    const mem_pestat_t *pestat_in,
    mem_pestat_t *pestat_out);

int bwa_shim_estimate_pestat(
    const BwaIndex *idx, const mem_opt_t *opts,
    const BwaReadPair *pairs, size_t n_pairs,
    mem_pestat_t *pestat_out);

/* ---- Three-phase API (caller-owned parallelism, cohort-exact output) ---- */

/* Per-thread reusable scratch: banded-SW buffers, SMEM buffers, chain/seed
 * windows, record-building buffers. ~24 MB after first use. Send, not Sync. */
BwaScratch *bwa_shim_scratch_new(void);
void        bwa_shim_scratch_free(BwaScratch *sc);

/* Phase 1: seed + single-end-extend every read of `batch` (the fused
 * `worker_bwt_aln` work). Per-read independent: safe to split a cohort into
 * any number of batches on any number of threads. Returns NULL + last_error
 * on failure. `opts` is never written. */
BwaRegs *bwa_shim_seed_extend(const BwaIndex *idx, const mem_opt_t *opts,
                              BwaScratch *sc, const BwaReadBatch *batch);
void     bwa_shim_regs_free(BwaRegs *r);
size_t   bwa_shim_regs_n_pairs(const BwaRegs *r);
size_t   bwa_shim_regs_n_singles(const BwaRegs *r);
/* Bytes held on the C heap by `r`: copied names/seqs/quals + alnreg arrays. */
size_t   bwa_shim_regs_heap_bytes(const BwaRegs *r);

/* Compatibility: wrap phase-1 output in a BwaSeeds so the legacy
 * bwa_shim_extend_batch can consume it. Takes ownership of `r`. */
BwaSeeds *bwa_shim_seeds_from_regs(BwaRegs *r);

/* Origin kinds for the record sink: whether a record came from the batch's
 * pairs (interleaved R1/R2) or its singles. */
#define BWA_ORIGIN_PAIR   0u
#define BWA_ORIGIN_SINGLE 1u

/* Global read ordinals, reproducing bwa-mem3's worker_sam id formulas
 * (bwamem.cpp:2795-2884): pair i of this batch gets id first_pair_id + i
 * (== (n_processed >> 1) + pos in the CLI); single i gets first_single_id + i
 * (== n_processed + i). The caller derives both from the cohort's global
 * read offset and the SE/PE group layout (fastmap.cpp:924-944). */
typedef struct { uint64_t first_single_id; uint64_t first_pair_id; } BwaIdBases;

/* Called once per emitted record with the packed BAM BODY (no u32 block_size
 * prefix). `origin_idx` indexes the batch's pairs or singles. The pointer is
 * valid only for the duration of the call. */
typedef void (*BwaRecordSinkFn)(void *ctx, uint32_t origin_kind, size_t origin_idx,
                                const uint8_t *body, size_t body_len);

/* Cohort insert-size model over the PE reads of several phase-1 batches
 * (mem_pestat once, over the concatenated per-read alnreg headers, in the
 * order given). Singles are ignored. `out` = mem_pestat_t[4]. Returns 0. */
int bwa_shim_pestat_cohort(const BwaIndex *idx, const mem_opt_t *opts,
                           const BwaRegs *const *regs, size_t n_regs, mem_pestat_t *out);

/* Phase 3: pairing + mate rescue + primary marking + emission (worker_sam).
 * Consumes `regs` (freed on every return path). `pestat` may be NULL only
 * when the batch has no pairs. Records are emitted in input order: pairs
 * (R1 side then R2 side, primary then supplementary), then singles. */
int bwa_shim_pair_emit(const BwaIndex *idx, const mem_opt_t *opts, BwaScratch *sc,
                       BwaRegs *regs, const mem_pestat_t *pestat, BwaIdBases ids,
                       BwaRecordSinkFn sink, void *ctx);

size_t         bwa_shim_batch_n_records (const BwaBatch *b);
size_t         bwa_shim_batch_pair_idx  (const BwaBatch *b, size_t rec);
const uint8_t *bwa_shim_batch_record_ptr(const BwaBatch *b, size_t rec);
size_t         bwa_shim_batch_record_len(const BwaBatch *b, size_t rec);
void           bwa_shim_batch_free      (BwaBatch *b);

const char *bwa_shim_last_error(void);
void        bwa_shim_set_verbosity(int level);

/* Vendored bwa-mem3 version (PACKAGE_VERSION, e.g. "0.9.0") — for `@PG VN:`. */
const char *bwa_shim_version(void);
/* Human-readable build description: "bwa-mem3 <version>; compiler: <line>". */
const char *bwa_shim_build_info(void);

/* Shared-memory index lifecycle. Thin wrappers over bwa-mem3's bwa_shm.h
 * (POSIX shm_open + a control segment named "/bwactl"). The shim's
 * `bwa_shim_idx_load` already attaches transparently when a segment named
 * after the prefix is staged; these entry points expose stage / drop /
 * list / probe to the Rust caller. */

/* Returns 1 if an index keyed by `prefix`'s basename is currently staged,
 * 0 if not, -1 on registry-access error. */
int bwa_shim_shm_test(const char *prefix);

/* Loads the index at `prefix` from disk, packs it, and stages it under
 * `/bwaidx-<basename>`. Returns 0 on success or if the prefix was already
 * staged, -1 on error. */
int bwa_shim_shm_stage(const char *prefix);

/* Drops every staged index segment plus the control segment. Returns 0
 * on success, -1 on error. Idempotent. */
int bwa_shim_shm_destroy(void);

/* Prints `<basename>\t<bytes>\n` for every staged segment to stdout (matches
 * `bwa shm -l`). Returns 0 on success, -1 on registry-access error. */
int bwa_shim_shm_list(void);

/* Set the @RG ID emitted as `RG:Z:` on all records. `id` may be NULL to
 * clear. Internally sets bwa-mem3's `bwa_rg_id[256]` global, which is
 * process-wide; callers aligning with different read groups from
 * multiple threads must serialize or use distinct processes. */
void bwa_shim_set_rg_id(const char *id);

#ifdef __cplusplus
}
#endif

#endif /* BWA_SHIM_H */
