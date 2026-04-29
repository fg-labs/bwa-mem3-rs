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

typedef struct {
    const char    *r1_name;  size_t r1_name_len;
    const uint8_t *r1_seq;   size_t r1_seq_len;
    const uint8_t *r1_qual;
    const char    *r2_name;  size_t r2_name_len;
    const uint8_t *r2_seq;   size_t r2_seq_len;
    const uint8_t *r2_qual;
} BwaReadPair;

/* Options lifecycle. `opts_new` returns a heap-allocated mem_opt_t populated
 * with bwa-mem3 defaults (mem_opt_init). `opts_free` releases it. Field-level
 * getters/setters will be added once the full mem_opt_t is bound to Rust. */
mem_opt_t *bwa_shim_opts_new(void);
void              bwa_shim_opts_free(mem_opt_t *opts);

/* Set common single-integer fields; one function per semantically-distinct knob.
 * Returns 0 on success, non-zero if the key is unknown. */
int bwa_shim_opts_set_int(mem_opt_t *opts, const char *key, int value);

/* PE-stats lifecycle. `pestat_zero` returns a zeroed 4-orientation array. */
mem_pestat_t *bwa_shim_pestat_zero(void);
void                 bwa_shim_pestat_free(mem_pestat_t *pestat);

BwaIndex *bwa_shim_idx_load(const char *prefix);
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

size_t         bwa_shim_batch_n_records (const BwaBatch *b);
size_t         bwa_shim_batch_pair_idx  (const BwaBatch *b, size_t rec);
const uint8_t *bwa_shim_batch_record_ptr(const BwaBatch *b, size_t rec);
size_t         bwa_shim_batch_record_len(const BwaBatch *b, size_t rec);
void           bwa_shim_batch_free      (BwaBatch *b);

const char *bwa_shim_last_error(void);
void        bwa_shim_set_verbosity(int level);

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
