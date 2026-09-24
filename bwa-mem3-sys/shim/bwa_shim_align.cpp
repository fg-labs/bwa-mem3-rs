/* bwa-mem3-sys/shim/bwa_shim_align.cpp
 *
 * Implements the shim's alignment functions by calling into bwa-mem3's
 * public API. Includes upstream's bwamem.h / FMI_search.h directly;
 * exposes a C-linkage bridge interface that bwa_shim.cpp consumes via
 * opaque pointers.
 *
 * Phase split:
 *   seed_batch  -> ShimSeeds  (owns worker_t + chains + copied bseq1_t[])
 *   extend_batch(ShimSeeds)   -> ShimAlignOutput  (packed BAM records)
 *   align_batch = seed + extend  (convenience)
 *   estimate_pestat = seed + SE-extend + mem_pestat  (no pairing, no emission)
 *
 * BAM emission is direct from mem_aln_t (no SAM intermediate).
 */

#include <mutex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "bwa.h"
#include "bwamem.h"
#include "compat_target.h"  /* compat_target_t::emit_hn — see the HN gate below */
#include "utils.h"          /* xassert — survives NDEBUG, unlike assert */
#include "FMI_search.h"
#include "meth_xm.h"   /* meth_build_xm — D3 (--meth) XM:Z tag */
#include "kswv.h"      /* kswr_t / SIMD_WIDTH8 — batched mate-rescue path below */
#include "bwa_shim_fields.h"  /* BwaAlignedFields — the structured-fields sink */

/* sort_classify has external linkage in bwamem.cpp but no header declaration
 * (the CLI's worker_sam calls it from the same TU). The batched mate-rescue
 * path in shim_pair_emit needs it, so forward-declare it here. */
extern int64_t sort_classify(mem_cache *mmc, int64_t pcnt, int tid);

/* Kernel thread slots a scratch may use: every `tid` indexes the scratch's own
 * mem_cache arrays (MAX_THREADS-sized) and the process-global tprof[][LIM_C]
 * profiling counters, so it must stay below both. stride_slot's 32 lines x 8
 * slots layout assumes exactly this many. */
#define SHIM_TID_SLOTS 256
static_assert(SHIM_TID_SLOTS <= MAX_THREADS && SHIM_TID_SLOTS <= LIM_C,
              "a shim tid slot must index both mem_cache and tprof");
static_assert(SHIM_TID_SLOTS == 32 * 8, "stride_slot lays out 32 lines of 8 slots");

/* Kernel thread slots handed to live scratches. The kernels bump
 * tprof[row][tid] counters on hot paths (per seed extension, per SA lookup);
 * tprof rows are tid-contiguous uint64_t, so 8 consecutive tids share a 64-byte
 * line. A new scratch takes the least-shared slot, first in stride-8 order (0,
 * 8, 16, ..., 248, 1, 9, ...), so on 64-byte-line hardware (Graviton, x86) up
 * to 32 live scratches -- one per pool worker on a 32-core box -- sit on 32
 * distinct lines, and up to 256 on distinct slots. On 128-byte-line hardware
 * (Apple Silicon) a 256-slot row spans only 16 lines, so 32 scratches pair up
 * two to a line whatever the layout; this one is still the best available.
 * Beyond 256 live scratches slots are shared, which shares nothing but racy
 * statistics; the per-slot count keeps a shared slot taken until its last
 * holder frees it. */
static std::mutex g_tid_mutex;
static unsigned g_tid_holders[SHIM_TID_SLOTS];

static int stride_slot(unsigned k) {
    return (int) ((k % 32) * 8 + (k / 32) % 8);
}

static int acquire_tid_slot(void) {
    std::lock_guard<std::mutex> lock(g_tid_mutex);
    int best = stride_slot(0);
    for (unsigned k = 0; k < SHIM_TID_SLOTS; ++k) {
        int slot = stride_slot(k);
        if (g_tid_holders[slot] == 0) { best = slot; break; }
        if (g_tid_holders[slot] < g_tid_holders[best]) best = slot;
    }
    ++g_tid_holders[best];
    return best;
}

static void release_tid_slot(int slot) {
    std::lock_guard<std::mutex> lock(g_tid_mutex);
    if (g_tid_holders[slot] > 0) --g_tid_holders[slot];
}

/* worker_alloc / worker_free — per-worker scratch (chaining arrays, BSW
 * buffers, lazy SMEM buffers) sized for one BATCH_SIZE chunk. Upstream defines
 * these in `fastmap.cpp`, but as of bwa-mem3 0.6.0 that TU is transitively
 * coupled to the new `fast_reader` FASTQ path (libdeflate + zlib-ng) and to
 * numa/htslib, none of which this crate compiles (the Rust CLI does its own
 * I/O). We therefore do not build `fastmap.cpp` (see `bwa-mem3-sys/build.rs`)
 * and carry file-local copies of just these two allocators here.
 *
 * The allocation set below is verbatim from upstream's `worker_alloc` /
 * `worker_free`, with ONE deliberate divergence: upstream's per-call
 * "Memory pre-allocation" stderr diagnostics are dropped. Upstream calls
 * worker_alloc once per thread at process start; this crate calls it once per
 * `ShimScratch` (shim_scratch_new), which the caller reuses across batches, so
 * those prints would spam stderr on every scratch.
 *
 * KEEP IN SYNC with `vendor/bwa-mem3/src/fastmap.cpp` on every vendor refresh:
 * the buffer set must match exactly what `mem_kernel1_core` /
 * `mem_kernel2_core` expect (the `worker_t` / `mem_cache` field layout comes
 * from the vendored `bwamem.h`). There is no compile-time guard for this the
 * way there is for the POD layout (bwa_shim_layout_assert.cpp) — the backstop
 * is the real-index integration suite (align_smoke / concurrency /
 * phase_split), so run it with `BWA_MEM3_RS_TEST_REF` set after any refresh.
 * A cleaner long-term fix is to split these into a dependency-light TU in the
 * fg-labs/bwa-mem3 fork and compile that instead of copying. */
/* Allocate the per-thread scratch a fused seed+extend needs, in kernel thread
 * slot `tid` (one ShimScratch per caller thread, so nthreads is 1). Upstream's
 * worker_alloc loops over every thread's slot; the loops below run over this
 * one slot but keep upstream's shape, so a vendor refresh diffs cleanly. Chain/seed
 * windows are BATCH_SIZE-sized because seeding and extension now run
 * chunk-by-chunk inside one call (upstream's v0.9.0 fusion, bwamem.cpp:2782):
 * chains never outlive the chunk that produced them. `regs` are NOT here --
 * they belong to the ShimRegs that outlives the call. */
static void worker_alloc(worker_t &w, int tid)
{
    xassert(tid >= 0 && tid < SHIM_TID_SLOTS, "worker_alloc: tid slot out of range");
    w.nthreads = 1;
    w.regs = NULL;
    /* Vestigial worker_t members that nothing allocates (every real
     * `auxSeedBuf` upstream is a local in test_and_merge, bwamem.cpp:894); zero
     * them so worker_free and any future upstream code cannot act on an
     * indeterminate pointer. */
    w.auxSeedBuf = NULL; w.auxSeedBufSize = 0;
    /* BATCH_SIZE-sized, chunk-local windows. Seeding and extension are fused
     * per BATCH_SIZE chunk in seed_extend_reads, so a chain is dead the moment
     * its chunk's mem_kernel2_core returns -- exactly upstream's v0.9.0 worker
     * layout, which no longer needs the pre-0.9.0 nreads-sized, seq_id-indexed
     * arrays this crate used to carry across the seed/extend phase barrier. */
    w.chain_scratch = (mem_chain_v*) malloc ((size_t)BATCH_SIZE * sizeof(mem_chain_v));
    w.seed_scratch  = (mem_seed_t *) calloc(sizeof(mem_seed_t), (size_t)BATCH_SIZE * AVG_SEEDS_PER_READ);
    /* xassert, not assert: these are ALLOCATION failures, and plain assert
     * compiles out under NDEBUG -- turning an OOM into a null deref in exactly
     * the build most likely to be memory-pressured (fastmap.cpp:306, :341). */
    xassert(w.seed_scratch  != NULL, "out of memory: w.seed_scratch");
    xassert(w.chain_scratch != NULL, "out of memory: w.chain_scratch");

    w.seed_scratch_size = BATCH_SIZE * AVG_SEEDS_PER_READ;

    /* SWA mem allocation. NOTE the basis changed in v0.9.0: upstream now sizes
     * this from AVG_SEEDS_PER_READ, not SEEDS_PER_READ (fastmap.cpp:380). */
    int64_t wsize = BATCH_SIZE * AVG_SEEDS_PER_READ;
    for (int l = tid; l <= tid; l++)
    {
        w.mmc.seqBufLeftRef[l*CACHE_LINE]  = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w.mmc.seqBufLeftQer[l*CACHE_LINE]  = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w.mmc.seqBufRightRef[l*CACHE_LINE] = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w.mmc.seqBufRightQer[l*CACHE_LINE] = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN, 64);

        w.mmc.wsize_buf_ref[l*CACHE_LINE] = wsize * MAX_SEQ_LEN_REF;
        w.mmc.wsize_buf_qer[l*CACHE_LINE] = wsize * MAX_SEQ_LEN_QER;

        xassert(w.mmc.seqBufLeftRef[l*CACHE_LINE]  != NULL, "out of memory: seqBufLeftRef");
        xassert(w.mmc.seqBufLeftQer[l*CACHE_LINE]  != NULL, "out of memory: seqBufLeftQer");
        xassert(w.mmc.seqBufRightRef[l*CACHE_LINE] != NULL, "out of memory: seqBufRightRef");
        xassert(w.mmc.seqBufRightQer[l*CACHE_LINE] != NULL, "out of memory: seqBufRightQer");
    }

    for (int l = tid; l <= tid; l++) {
        w.mmc.seqPairArrayAux[l]      = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.seqPairArrayLeft128[l]  = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.seqPairArrayRight128[l] = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.wsize[l] = wsize;

        xassert(w.mmc.seqPairArrayAux[l] != NULL, "out of memory: seqPairArrayAux");
        xassert(w.mmc.seqPairArrayLeft128[l] != NULL, "out of memory: seqPairArrayLeft128");
        xassert(w.mmc.seqPairArrayRight128[l] != NULL, "out of memory: seqPairArrayRight128");
    }

    // SMEM buffers (matchArray / min_intv_ar / query_pos_ar / enc_qdb / rid)
    // and the lockstep-batch slot buffers are sized from the observed max
    // read length on each batch in mem_collect_smem; they're NULL here and
    // grow on first use. `lim` is still a fixed BATCH_SIZE+32 allocation
    // because its size does not depend on read length.
    for (int l = tid; l <= tid; l++)
    {
        w.mmc.wsize_mem[l]     = 0;
        w.mmc.wsize_mem_s[l]   = 0;
        w.mmc.wsize_mem_r[l]   = 0;
        w.mmc.wsize_qdb[l]     = 0;
        w.mmc.matchArray[l]    = NULL;
        w.mmc.min_intv_ar[l]   = NULL;
        w.mmc.query_pos_ar[l]  = NULL;
        w.mmc.enc_qdb[l]       = NULL;
        w.mmc.rid[l]           = NULL;
        w.mmc.lim[l]           = (int32_t *) _mm_malloc((BATCH_SIZE + 32) * sizeof(int32_t), 64);

        w.mmc.lockstep_prev[l]      = NULL;
        w.mmc.lockstep_match_buf[l] = NULL;
        w.mmc.lockstep_buf_cap[l]   = 0;

        w.mmc.smem_sort_scratch[l].cnt    = NULL;
        w.mmc.smem_sort_scratch[l].cntCap = 0;
        w.mmc.smem_sort_scratch[l].tmp    = NULL;
        w.mmc.smem_sort_scratch[l].tmpCap = 0;
    }
}

static void worker_free(worker_t &w, int tid)
{
    // Catch mismatched alloc/free pairs before they drive out-of-bounds frees.
    assert(w.nthreads == 1);

    free(w.chain_scratch);
    /* w.regs is NOT freed here: during seed_extend_reads it aliases the
     * caller-owned ShimRegs::regs, which shim_regs_free owns. worker_alloc
     * leaves it NULL and nothing in the scratch ever allocates it. */
    free(w.seed_scratch);

    for (int l = tid; l <= tid; l++) {
        _mm_free(w.mmc.seqBufLeftRef[l*CACHE_LINE]);
        _mm_free(w.mmc.seqBufRightRef[l*CACHE_LINE]);
        _mm_free(w.mmc.seqBufLeftQer[l*CACHE_LINE]);
        _mm_free(w.mmc.seqBufRightQer[l*CACHE_LINE]);
    }

    for (int l = tid; l <= tid; l++) {
        free(w.mmc.seqPairArrayAux[l]);
        free(w.mmc.seqPairArrayLeft128[l]);
        free(w.mmc.seqPairArrayRight128[l]);
    }

    // NULL-safe: SMEM buffers are now allocated lazily on first batch;
    // workers that never ran a batch leave them as NULL. _mm_free / free
    // are both well-defined on NULL.
    for (int l = tid; l <= tid; l++) {
        _mm_free(w.mmc.matchArray[l]);
        free(w.mmc.min_intv_ar[l]);
        free(w.mmc.query_pos_ar[l]);
        free(w.mmc.enc_qdb[l]);
        free(w.mmc.rid[l]);
        _mm_free(w.mmc.lim[l]);

        _mm_free(w.mmc.lockstep_prev[l]);
        _mm_free(w.mmc.lockstep_match_buf[l]);

        _mm_free(w.mmc.smem_sort_scratch[l].cnt);
        _mm_free(w.mmc.smem_sort_scratch[l].tmp);
    }
}

extern "C" {

struct ShimReadPair {
    const char    *r1_name;  size_t r1_name_len;
    const uint8_t *r1_seq;   size_t r1_seq_len;
    const uint8_t *r1_qual;
    const char    *r2_name;  size_t r2_name_len;
    const uint8_t *r2_seq;   size_t r2_seq_len;
    const uint8_t *r2_qual;
};

/* Output: concatenated packed-BAM records + index table.
 * Each record: [u32 le block_size][block_size bytes of BAM record data]. */
struct ShimAlignOutput {
    uint8_t *buf;       size_t buf_len;
    size_t  *rec_off;   size_t  *rec_len;
    size_t  *pair_idx;
    size_t   n_recs;    size_t cap;
    size_t   buf_cap;
    mem_pestat_t pes[4];
};

/* Per-thread reusable scratch (public: BwaScratch). */
struct ShimScratch {
    worker_t w;              /* mmc + BATCH_SIZE chain/seed windows; nthreads = 1 */
    int tid;                 /* kernel thread slot: w.mmc's live entry and tprof column */
    uint8_t *rec_buf; size_t rec_cap;   /* one packed record under construction */
    uint8_t *aux_buf; size_t aux_cap;   /* one aux block under construction */
    /* One record's computed fields (compute_record_fields): the emitted CIGAR
     * and the MC:Z / SA:Z texts. Grown in place and reused, like aux_buf. */
    uint32_t *cig_buf; size_t cig_cap;  /* cap in ops */
    uint8_t  *txt_buf; size_t txt_cap;
};

/* Phase-1 output (public: BwaRegs). Owns the read copies and the per-read
 * alnreg arrays; borrows nothing from the scratch. */
struct ShimRegs {
    bseq1_t      *seqs;      int n_seqs;      /* [0,2*n_pairs) pairs, then singles */
    size_t        n_pairs;   size_t n_singles;
    mem_alnreg_v *regs;                       /* n_seqs entries */
    int           meth_mode;
    size_t        heap_bytes;                 /* names+seqs+quals (+meth_orig) */
};

struct ShimSingleRead {
    const char    *name;  size_t name_len;
    const uint8_t *seq;   size_t seq_len;
    const uint8_t *qual;
};
struct ShimReadBatch {
    const ShimReadPair   *pairs;   size_t n_pairs;
    const ShimSingleRead *singles; size_t n_singles;
};

typedef void (*ShimRecordSinkFn)(void *ctx, uint32_t origin_kind, size_t origin_idx,
                                 const uint8_t *body, size_t body_len);
typedef BwaFieldSinkFn ShimFieldSinkFn;
struct ShimIdBases { uint64_t first_single_id; uint64_t first_pair_id; };

/* Where a record goes once built. Replaces the ShimAlignOutput* that
 * append_bam_record used to write into; the legacy output type is now just
 * one particular sink (legacy_out_sink below). Exactly one of `sink` (packed
 * BAM body) and `field_sink` (structured fields) is set; `field_sink` is last
 * so the existing `{ sc, sink, ctx, kind }` initializers leave it NULL. */
struct ShimEmit {
    ShimScratch      *sc;
    ShimRecordSinkFn  sink;
    void             *ctx;
    uint32_t          origin_kind;
    ShimFieldSinkFn   field_sink;
};

/* Legacy sink: append `[u32 block_size][body]` to a ShimAlignOutput and index it. */
static void legacy_out_sink(void *ctx, uint32_t /*kind*/, size_t origin_idx,
                            const uint8_t *body, size_t body_len)
{
    ShimAlignOutput *out = (ShimAlignOutput *) ctx;
    size_t rec_size = 4 + body_len;
    if (out->buf_len + rec_size > out->buf_cap) {
        while (out->buf_len + rec_size > out->buf_cap) out->buf_cap *= 2;
        uint8_t *tmp = (uint8_t *) realloc(out->buf, out->buf_cap);
        xassert(tmp != NULL, "out of memory: legacy_out_sink buf");
        out->buf = tmp;
    }
    if (out->n_recs == out->cap) {
        out->cap *= 2;
        size_t *tmp_off = (size_t *) realloc(out->rec_off,  out->cap * sizeof(size_t));
        xassert(tmp_off != NULL, "out of memory: legacy_out_sink rec_off");
        out->rec_off = tmp_off;
        size_t *tmp_len = (size_t *) realloc(out->rec_len,  out->cap * sizeof(size_t));
        xassert(tmp_len != NULL, "out of memory: legacy_out_sink rec_len");
        out->rec_len = tmp_len;
        size_t *tmp_pidx = (size_t *) realloc(out->pair_idx, out->cap * sizeof(size_t));
        xassert(tmp_pidx != NULL, "out of memory: legacy_out_sink pair_idx");
        out->pair_idx = tmp_pidx;
    }
    uint32_t bs32 = (uint32_t) body_len;
    memcpy(out->buf + out->buf_len, &bs32, 4);
    memcpy(out->buf + out->buf_len + 4, body, body_len);
    out->rec_off[out->n_recs]  = out->buf_len;
    out->rec_len[out->n_recs]  = rec_size;
    out->pair_idx[out->n_recs] = origin_idx;
    out->n_recs++;
    out->buf_len += rec_size;
}

/* Legacy phase-1 handle: now just an owned ShimRegs. */
struct ShimSeeds { ShimRegs *regs; };

/* Index handle: owns the loaded FMI_search (plus the original bns/pac in
 * --meth mode). The reference is never unpacked: extension reconstructs each
 * window from the resident `.pac` on demand (bns_get_seq_v2's
 * `ref_string == NULL` path — the non-meth pac via `fmi`, the meth ORIGINAL pac
 * via `meth_orig_pac`), matching upstream bwa-mem3 and saving the ~3s startup
 * unpack + ~6GB array. */
struct BwaShimIndex {
    FMI_search *fmi;
    /* D3 (--meth) dual-coordinate handles. In meth mode `fmi` is the
     * f/r-doubled CONVERTED seed index (`<ref>.meth.*`) used only for candidate
     * generation, while chaining/extension/output run in ORIGINAL coordinates
     * loaded here from the un-converted `<ref>.*` prefix. Both NULL for a normal
     * (non-meth) index; the original reference is pac-fetched from
     * `meth_orig_pac` on demand, never unpacked. */
    bntseq_t *meth_orig_bns;
    uint8_t  *meth_orig_pac;
};

/* ------------------ Index ------------------ */

void shim_align_idx_free(void *opaque);  /* defined below; used by _load_meth */

/* D3 (--meth): load the ORIGINAL reference's bns + pac from `prefix` as
 * resident handles (distinct from the converted seed FM-index). Ported from
 * upstream fastmap.cpp's meth_orig_ref_load_handles: bns_restore then slurp the
 * whole .pac into memory and close the file. Returns 0 on success (writes
 * *bns_out / *pac_out), -1 on failure (out-params left NULL). */
static int shim_meth_orig_ref_load(const char *prefix,
                                    bntseq_t **bns_out, uint8_t **pac_out) {
    *bns_out = nullptr;
    *pac_out = nullptr;
    bntseq_t *bns = bns_restore(prefix);
    if (bns == nullptr) return -1;
    int64_t pac_bytes = bns->l_pac / 4 + 1;
    uint8_t *pac = (uint8_t *) calloc(pac_bytes, 1);
    if (pac == nullptr) {
        bns_destroy(bns);
        return -1;
    }
    /* bns_restore left .pac open in bns->fp_pac; slurp it whole, then close. */
    err_fread_noeof(pac, 1, pac_bytes, bns->fp_pac);
    err_fclose(bns->fp_pac);
    bns->fp_pac = nullptr;
    *bns_out = bns;
    *pac_out = pac;
    return 0;
}

/* As shim_align_idx_load, loading the FM-index with `n_threads` (>= 1)
 * threads (the CLI's `-t` behavior for index load, fastmap.cpp:2868).
 * `n_threads < 1` is clamped to 1. */
void *shim_align_idx_load_threads(const char *prefix, int n_threads) {
    if (n_threads < 1) n_threads = 1;
    FMI_search *fmi = new FMI_search(prefix);
    fmi->load_index(/*load_pac=*/true, n_threads);
    BwaShimIndex *idx = (BwaShimIndex *) calloc(1, sizeof(BwaShimIndex));
    if (!idx) {
        delete fmi;
        return nullptr;
    }
    idx->fmi = fmi;
    /* pac-fetch: do NOT materialize the unpacked 2*l_pac reference. With no
     * materialized reference, extension reconstructs each window from the
     * resident `.pac` on demand (bns_get_seq_v2's `ref_string == NULL` path,
     * per-thread scratch) — byte-identical to the unpacked array it replaces,
     * and exactly what upstream bwa-mem3's `mem` does (upstream removed the
     * unpacked `.0123` entirely). This drops the ~3s single-threaded startup
     * unpack and the ~6 GB resident array; the `pac` stays owned by `fmi`
     * (loaded with load_pac=true above). The meth path pac-fetches its
     * ORIGINAL reference the same way (see shim_align_idx_load_meth). */
    return static_cast<void *>(idx);
}

void *shim_align_idx_load(const char *prefix) {
    return shim_align_idx_load_threads(prefix, /*n_threads=*/1);
}

/* D3 (--meth): load a dual index. `seed_prefix` is the converted seed index
 * (`<ref>.meth`), `orig_prefix` the un-converted original reference (`<ref>`).
 * The seed FM-index plus the original bns/pac stay resident; the original
 * reference is NOT unpacked — extension pac-fetches it from `meth_orig_pac` on
 * demand (bns_get_seq_v2's `ref_string == NULL` path via mem_kernel2_core's
 * meth aln_pac routing), matching upstream bwa-mem3. */
void *shim_align_idx_load_meth(const char *seed_prefix, const char *orig_prefix) {
    void *opaque = shim_align_idx_load(seed_prefix);
    if (!opaque) return nullptr;
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(opaque);
    if (shim_meth_orig_ref_load(orig_prefix, &idx->meth_orig_bns, &idx->meth_orig_pac) != 0) {
        shim_align_idx_free(opaque);
        return nullptr;
    }
    return opaque;
}

void shim_align_idx_free(void *opaque) {
    if (!opaque) return;
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(opaque);
    if (idx->meth_orig_pac) free(idx->meth_orig_pac);
    if (idx->meth_orig_bns) bns_destroy(idx->meth_orig_bns);
    delete idx->fmi;
    free(idx);
}

/* Contigs reported to callers (used to build the BAM header) come from the
 * ORIGINAL reference in --meth mode — emitted records use original rids/coords,
 * not the f/r-doubled converted seed index. Outside --meth this is the seed
 * index's bns. */
static const bntseq_t *shim_header_bns(void *opaque) {
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(opaque);
    return idx->meth_orig_bns ? idx->meth_orig_bns : idx->fmi->idx->bns;
}

/* Non-zero iff `opaque` was loaded as a --meth dual index (shim_align_idx_load_meth). */
int shim_align_idx_is_meth(void *opaque) {
    return static_cast<BwaShimIndex *>(opaque)->meth_orig_bns != nullptr;
}

size_t shim_align_idx_n_contigs(void *opaque) {
    return (size_t) shim_header_bns(opaque)->n_seqs;
}

const char *shim_align_idx_contig_name(void *opaque, size_t i) {
    return shim_header_bns(opaque)->anns[i].name;
}

int64_t shim_align_idx_contig_len(void *opaque, size_t i) {
    return shim_header_bns(opaque)->anns[i].len;
}

/* The @HD line the selected compat target wants, or NULL when it suppresses
 * @HD entirely (bwa-mem2 emits none; see compat_target.cpp).
 *
 * This must be read from the target rather than spelled out by the caller.
 * Upstream consolidated three hardcoded @HD literals into the single
 * BWAMEM3_DEFAULT_HD_LINE macro precisely because they had drifted apart
 * (fg-labs/bwa-mem3#288: the SAM-text writer said "VN:1.5 SO:unsorted
 * GO:query" while the BAM writers said "VN:1.6 SO:unsorted"), so a fourth
 * literal living in this crate's Rust header writer would recreate the bug
 * upstream just finished fixing -- and it did: bwa-rs emitted
 * "@HD VN:1.6 SO:unknown" against the CLI's "@HD VN:1.5 SO:unsorted GO:query".
 *
 * `compat` is set by mem_opt_init, but it is a pointer on a struct callers can
 * memset, so it is NULL-checked like the emit_mq / emit_hn gates. */
const char *shim_compat_hd_line(const mem_opt_t *opt) {
    if (opt == NULL || opt->compat == NULL || !opt->compat->emit_hd) return NULL;
    return opt->compat->hd_line;
}

/* Apply upstream's bwameth-compatibility defaults for --meth.
 *
 * Delegates to mem_opt_apply_meth_defaults (bwamem.cpp:504) rather than
 * replicating the bundle. The replicated version this replaces applied
 * bwameth's constants FLAT, but they are quoted at bwameth's match score
 * (a == 1) and upstream scales each by opt->a -- so `--meth` combined with a
 * non-default -A silently discarded it, leaving T at 40 while the alignment
 * scores it gates had doubled. Gotcha #13 exists for this class of drift:
 * where upstream factors a policy out, call it.
 *
 * The `opt0` sentinel is upstream's "did the user set this explicitly" mask
 * (non-zero field == user supplied). A zeroed one means "nothing was set", so
 * every default applies. Taking no mask is the deliberate simplification:
 * expressing one would mean tracking per-field "was set" state on MemOpts.
 *
 * That makes the ordering contract asymmetric, in three ways:
 *
 *   - `a` and `meth_scoring` are INPUTS -- the constants are expressed in units
 *     of `a`, and the -B branch keys off the resolved scoring mode -- so both
 *     must be set BEFORE this call.
 *   - `T`, `pen_clip5`, `pen_clip3` and `pen_unpaired` are written
 *     unconditionally (the empty mask says nobody set them), so a caller
 *     wanting its own must set them AFTER.
 *   - `b` is written only under COLLAPSED. GENOMIC and NEUTRAL are
 *     variant-aware and keep bwa's default, because their mirror cell must stay
 *     a real mismatch (bwamem.cpp:511-515), so a caller-set -B survives the
 *     bundle under those two modes and is clobbered under COLLAPSED.
 *
 * Refilling the matrices is part of the operation, not the caller's job: the
 * COLLAPSED branch can change opt->b, and a stale mat/mat_ot/mat_ob would
 * score every subsequent alignment with the pre-default penalty. Upstream
 * refills at the same point (fastmap.cpp:2726, :2730).
 *
 * NOT replicated here: upstream's TAPS => NEUTRAL scoring default, which it
 * applies just BEFORE this call (fastmap.cpp:2695-2696). It is unreachable for
 * this crate -- mem_opt_init defaults meth_chem to METH_CHEM_EMSEQ and no API
 * exposes meth_chem -- but it must be added here if one ever does, because the
 * COLLAPSED -B branch below keys off the resolved scoring mode. */
void shim_opts_apply_meth_defaults(mem_opt_t *opt) {
    if (opt == NULL) return;
    mem_opt_t opt0;
    memset(&opt0, 0, sizeof(opt0));
    mem_opt_apply_meth_defaults(opt, &opt0);
    bwa_fill_scmat(opt->a, opt->b, opt->mat);
    mem_opt_fill_meth_mat(opt);
}

/* ------------------ Helpers ------------------ */

/* Forward-declared so the copy helpers can unwind a partially-built array on an
 * allocation failure (free_seqs is NULL-safe on the un-built calloc'd tail). */
static void free_seqs(bseq1_t *seqs, int nseqs);

/* The --meth read setup shared by every path that copies reads in, applied
 * after the original bases were saved to s->meth_orig_seq: record the read's
 * chemistry and project s->seq in place. `role` 0 = R1 or a single (OT, C->T),
 * 1 = R2 (OB, G->A). */
static void meth_project_in_place(bseq1_t *s, int role) {
    /* v0.9.0: read-number chemistry, consumed by the seed filter in
     * meth_seed_to_orig (bwamem.cpp:1878-1881), which drops any seed
     * whose genomic strand does not match the read's. R1 = OT = 1,
     * R2 = OB = 0, exactly as upstream's ingest sets it
     * (fastmap.cpp:820).
     *
     * Leaving it unset is NOT inert: the field is zero from the
     * calloc, and 0 is the VALID value for OB, so the filter runs
     * and silently discards every R1 seed. The filter only
     * self-disables at < 0. That presented as exactly half the
     * reads unmapped -- every R1, no R2. */
    s->meth_base_ot = (role == 0) ? 1 : 0;
    char from = (role == 0) ? 'C' : 'G';
    char to   = (role == 0) ? 'T' : 'A';
    for (int j = 0; j < s->l_seq; ++j)
        if (s->seq[j] == from || s->seq[j] == (char)(from + 32))
            s->seq[j] = to;
}

/* Decode one read into `s` from borrowed bytes: copy name/seq/qual into owned
 * buffers and, under --meth, keep the original bases and apply the bisulfite
 * projection. `role` 0 = R1 or a single (OT, C->T), 1 = R2 (OB, G->A).
 * Returns 0, or -1 on OOM; every buffer it did allocate is recorded in `s`, so
 * the caller's free_seqs / segment teardown (NULL-safe on the fields left
 * unbuilt) releases a partial decode. Used by the legacy batch path and the
 * per-slot resident writes; decode_read_arena is its arena-backed twin, and the
 * two share meth_project_in_place so they cannot project differently. */
static int decode_read(bseq1_t *s, const char *name, size_t name_len,
                       const uint8_t *seq, size_t seq_len, const uint8_t *qual,
                       int meth_mode, int role, size_t *heap_bytes) {
    s->l_seq = (int) seq_len;
    s->name = (char *) malloc(name_len + 1);
    if (!s->name) return -1;
    memcpy(s->name, name, name_len); s->name[name_len] = '\0';
    s->seq = (char *) malloc(seq_len + 1);
    if (!s->seq) return -1;
    memcpy(s->seq, seq, seq_len); s->seq[seq_len] = '\0';
    if (qual) {
        s->qual = (char *) malloc(seq_len + 1);
        if (!s->qual) return -1;
        memcpy(s->qual, qual, seq_len); s->qual[seq_len] = '\0';
    } else {
        s->qual = nullptr;
    }
    *heap_bytes += (name_len + 1) + (seq_len + 1) + (qual ? seq_len + 1 : 0);
    s->sam = nullptr;
    /* D3 (--meth): retain the ORIGINAL (unconverted) read bases before
     * projecting seq to match the converted seed index, then apply the
     * per-strand bisulfite projection in place. R1 is the OT/CT read (C→T), R2
     * the OB/GA read (G→A) — matching upstream fastmap.cpp's YC assignment; a
     * single is treated as R1, matching upstream's single-file ingest.
     * mem_kernel1_core 2-bit-encodes the projected seq; mem_reg2aln /
     * meth_build_xm use meth_orig_seq for CIGAR/NM/MD and the XM call string,
     * and serialize_record emits it as SEQ. Same orientation as seq. */
    if (meth_mode) {
        /* Copy exactly l_seq bytes, not strdup: a base byte of 0 would make
         * strdup stop short and every later reader of the l_seq original bases
         * run off the end -- and diverge from decode_read_arena's copy. */
        s->meth_orig_seq = (char *) malloc(seq_len + 1);
        if (!s->meth_orig_seq) return -1;
        memcpy(s->meth_orig_seq, s->seq, seq_len); s->meth_orig_seq[seq_len] = '\0';
        *heap_bytes += (size_t) s->l_seq + 1;
        meth_project_in_place(s, role);
    }
    return 0;
}

/* Bump cursor over one resident segment's read arena. */
struct ArenaCursor { uint8_t *p; };

/* Place `len` bytes plus a NUL at the cursor and advance it. */
static char *arena_put(ArenaCursor *cur, const void *src, size_t len) {
    char *dst = (char *) cur->p;
    memcpy(dst, src, len);
    dst[len] = '\0';
    cur->p += len + 1;
    return dst;
}

/* Bytes decode_read_arena places for one read: the same total decode_read
 * adds to its heap_bytes, so the two paths account identically. */
static size_t read_arena_bytes(size_t name_len, size_t seq_len, const uint8_t *qual,
                               int meth_mode) {
    return (name_len + 1) + (seq_len + 1) + (qual ? seq_len + 1 : 0)
         + (meth_mode ? seq_len + 1 : 0);
}

/* decode_read, placing every string at `cur` instead of in its own malloc:
 * the same bytes and terminators, and under --meth the same original-base copy
 * (taken before the projection) and the same projection. */
static void decode_read_arena(bseq1_t *s, ArenaCursor *cur, const char *name, size_t name_len,
                              const uint8_t *seq, size_t seq_len, const uint8_t *qual,
                              int meth_mode, int role) {
    s->l_seq = (int) seq_len;
    s->name = arena_put(cur, name, name_len);
    s->seq = arena_put(cur, seq, seq_len);
    s->qual = qual ? arena_put(cur, qual, seq_len) : nullptr;
    s->sam = nullptr;
    if (meth_mode) {
        s->meth_orig_seq = arena_put(cur, s->seq, seq_len);
        meth_project_in_place(s, role);
    }
}

static bseq1_t *copy_pairs_to_seqs(const ShimReadPair *pairs, size_t n_pairs,
                                   int meth_mode, size_t *heap_bytes) {
    *heap_bytes = 0;
    int nseqs = (int)(2 * n_pairs);
    /* NULL is an ALLOCATION FAILURE only when nseqs > 0. calloc(0, n) may
     * legally return NULL, so on an empty batch a NULL return here means "no
     * reads", not "out of memory" -- the caller distinguishes the two on
     * n_pairs. Downstream is NULL-safe either way: every loop over the array is
     * bounded by n_seqs, and free_seqs() returns early on NULL. */
    bseq1_t *seqs = (bseq1_t *) calloc(nseqs, sizeof(bseq1_t));
    if (!seqs) return nullptr;
    for (size_t i = 0; i < n_pairs; ++i) {
        const ShimReadPair *p = &pairs[i];
        /* An OOM unwinds the partially-built array and returns NULL (the
         * caller's `!pairs_only` path); free_seqs() is NULL-safe on the entries
         * not yet built (calloc leaves them NULL). */
        if (decode_read(&seqs[2*i], p->r1_name, p->r1_name_len, p->r1_seq, p->r1_seq_len,
                        p->r1_qual, meth_mode, 0, heap_bytes) != 0
            || decode_read(&seqs[2*i + 1], p->r2_name, p->r2_name_len, p->r2_seq,
                           p->r2_seq_len, p->r2_qual, meth_mode, 1, heap_bytes) != 0) {
            free_seqs(seqs, nseqs);
            return nullptr;
        }
    }
    return seqs;
}

/* Append `n` single reads at seqs[base..base+n). Same ownership rules as
 * copy_pairs_to_seqs. Under --meth a single is treated as an R1/OT read
 * (C→T, meth_base_ot = 1), matching upstream's single-file ingest.
 *
 * Returns 0 on success, -1 on allocation failure. On failure the caller owns
 * `seqs` (it is the caller's r->seqs) and frees it via free_seqs, which is
 * NULL-safe on the entries this function had not yet built. */
static int copy_singles_to_seqs(bseq1_t *seqs, size_t base, const ShimSingleRead *singles,
                                size_t n, int meth_mode, size_t *heap_bytes)
{
    for (size_t i = 0; i < n; ++i) {
        const ShimSingleRead *p = &singles[i];
        if (decode_read(&seqs[base + i], p->name, p->name_len, p->seq, p->seq_len, p->qual,
                        meth_mode, 0, heap_bytes) != 0)
            return -1;
    }
    return 0;
}

static void free_seqs(bseq1_t *seqs, int nseqs) {
    if (!seqs) return;
    for (int i = 0; i < nseqs; ++i) {
        free(seqs[i].name);
        free(seqs[i].seq);
        free(seqs[i].qual);
        free(seqs[i].sam);
        free(seqs[i].meth_orig_seq);  /* NULL outside --meth; free() is NULL-safe */
    }
    free(seqs);
}

/* ------------------ BAM emission (direct from mem_aln_t) ------------------ */

static inline uint16_t reg2bin(int beg, int end) {
    end--;
    if (beg >> 14 == end >> 14) return ((1 << 15) - 1) / 7 + (beg >> 14);
    if (beg >> 17 == end >> 17) return ((1 << 12) - 1) / 7 + (beg >> 17);
    if (beg >> 20 == end >> 20) return ((1 <<  9) - 1) / 7 + (beg >> 20);
    if (beg >> 23 == end >> 23) return ((1 <<  6) - 1) / 7 + (beg >> 23);
    if (beg >> 26 == end >> 26) return ((1 <<  3) - 1) / 7 + (beg >> 26);
    return 0;
}

/* BAM 4-bit base encoding for bwa-mem3's 2-bit-encoded bytes.
 *
 * `mem_kernel1_core` rewrites `s->seq` in place via nst_nt4_table: the bytes
 * become 0=A, 1=C, 2=G, 3=T, 4=N. We emit them as BAM 4-bit nibbles:
 *   A=1, C=2, G=4, T=8, N=15.
 *
 * For reverse-strand emission we also need the complement map (in bwa's
 * 2-bit space: 0<->3, 1<->2, N stays N).
 */
static const uint8_t bwa2_to_bam4[8] = {
    /* 0 A */ 1,
    /* 1 C */ 2,
    /* 2 G */ 4,
    /* 3 T */ 8,
    /* 4 N */ 15,
    /* 5-7 */ 15, 15, 15,
};
static inline uint8_t bwa2_complement(uint8_t b) {
    return (b < 4) ? (uint8_t)(3 - b) : (uint8_t)4;
}

/* bwa-mem3 uses a 5-char CIGAR opcode table "MIDSH" (S=3, H=4), incompatible
 * with the BAM spec's "MIDNSHP=X" (S=4, H=5). Remap when emitting to BAM. */
static inline uint32_t bwa_cigar_to_bam(uint32_t op) {
    static const uint8_t BWA_TO_BAM[5] = {0, 1, 2, 4, 5}; /* M I D S H */
    uint32_t len = op >> 4;
    uint32_t o = op & 0xf;
    uint32_t bam_o = (o < 5) ? BWA_TO_BAM[o] : o; /* defensive: passthrough unknown */
    return (len << 4) | bam_o;
}

/* Mirror of upstream's `add_cigar` clip rewrite (bwamem.cpp:2410-2421).
 * bwa-mem3 HARD-clips supplementary records by default: on a 2nd+ emitted
 * record (`which != 0`) every clip opcode becomes H, and S otherwise --
 * unless -Y (MEM_F_SOFTCLIP) or the alignment is on an ALT contig, in which
 * case the soft clip stands. `mem_reg2aln` only ever writes S (bwamem.cpp:
 * 2735,2743 build `clip<<4 | 3`), so it is this rewrite, not the aligner,
 * that produces every H the CLI emits.
 *
 * Opcodes here are bwa's MIDSH table (S=3, H=4), NOT the BAM spec's; feed the
 * result through bwa_cigar_to_bam before writing it. */
static inline uint32_t bwa_apply_clip_mode(uint32_t cig, const mem_opt_t *opt,
                                           int is_alt, int which) {
    uint32_t o = cig & 0xf;
    if (!(opt->flag & MEM_F_SOFTCLIP) && !is_alt && (o == 3 || o == 4))
        o = which ? 4 : 3;
    return (cig & ~(uint32_t) 0xf) | o;
}

/* BAM 4-bit code of an ASCII base, as htslib's seq_nt16_table assigns it
 * (IUPAC ambiguity codes keep their own code; case-insensitive; anything else
 * is N). */
static inline uint8_t ascii_to_bam4(unsigned char c) {
    switch (c & ~0x20) {   /* fold to upper case */
        case 'A': return 1;  case 'C': return 2;  case 'M': return 3;  case 'G': return 4;
        case 'R': return 5;  case 'S': return 6;  case 'V': return 7;  case 'T': return 8;
        case 'U': return 8;  case 'W': return 9;  case 'Y': return 10; case 'H': return 11;
        case 'K': return 12; case 'D': return 13; case 'B': return 14;
        default:  return (c == '=') ? 0 : 15;
    }
}

/* One emitted --meth SEQ base, as upstream's meth writer (meth_bam.cpp)
 * builds it from the ORIGINAL (unprojected) read so MethylDackel sees real
 * C/Ts: the forward strand keeps the base (upper-cased, so an IUPAC code
 * survives), the reverse strand complements through nst_nt4_table and maps
 * anything but A/C/G/T to N. */
static inline uint8_t meth_seq_bam4(unsigned char c, int rev) {
    if (!rev) return ascii_to_bam4(c);
    static const uint8_t RC[5] = { 8 /* A->T */, 4 /* C->G */, 2 /* G->C */, 1 /* T->A */, 15 };
    int bi = nst_nt4_table[c];
    return RC[bi < 4 ? bi : 4];
}

/* ASCII base complement (A<->T, C<->G, case-insensitive); anything else -> N.
 * Used to build the D3 (--meth) XM:Z seq_text in emitted-SEQ orientation. */
static inline char ascii_complement(char b) {
    switch (b) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default:            return 'N';
    }
}

/* CIGAR op length consuming reference. Operates on bwa-mem3's 5-op encoding
 * (M=0, I=1, D=2, S=3, H=4) where reference-consumers are M and D only. */
static int cigar_ref_len(int n_cigar, const uint32_t *cigar) {
    int ref = 0;
    for (int i = 0; i < n_cigar; ++i) {
        int op = cigar[i] & 0xf;
        int len = (int)(cigar[i] >> 4);
        if (op == 0 || op == 2) ref += len; /* M or D */
    }
    return ref;
}

/* Append a byte block to a dynamically-growing buffer. */
static void buf_append(uint8_t **buf, size_t *len, size_t *cap, const void *src, size_t n) {
    if (*len + n > *cap) {
        while (*len + n > *cap) *cap = *cap ? *cap * 2 : 1024;
        uint8_t *tmp = (uint8_t *) realloc(*buf, *cap);
        xassert(tmp != NULL, "out of memory: buf_append");
        *buf = tmp;
    }
    memcpy(*buf + *len, src, n);
    *len += n;
}

/* Emit the smallest-width BAM-encoded integer aux value. */
static void aux_put_i(uint8_t **buf, size_t *len, size_t *cap,
                      const char tag[2], int64_t v) {
    uint8_t tmp[3 + 4];
    tmp[0] = (uint8_t)tag[0]; tmp[1] = (uint8_t)tag[1];
    if (v >= INT8_MIN && v <= INT8_MAX) {
        tmp[2] = 'c'; int8_t x = (int8_t)v; memcpy(tmp + 3, &x, 1);
        buf_append(buf, len, cap, tmp, 4);
    } else if (v >= 0 && v <= UINT8_MAX) {
        tmp[2] = 'C'; uint8_t x = (uint8_t)v; memcpy(tmp + 3, &x, 1);
        buf_append(buf, len, cap, tmp, 4);
    } else if (v >= INT16_MIN && v <= INT16_MAX) {
        tmp[2] = 's'; int16_t x = (int16_t)v; memcpy(tmp + 3, &x, 2);
        buf_append(buf, len, cap, tmp, 5);
    } else if (v >= 0 && v <= UINT16_MAX) {
        tmp[2] = 'S'; uint16_t x = (uint16_t)v; memcpy(tmp + 3, &x, 2);
        buf_append(buf, len, cap, tmp, 5);
    } else {
        tmp[2] = 'i'; int32_t x = (int32_t)v; memcpy(tmp + 3, &x, 4);
        buf_append(buf, len, cap, tmp, 7);
    }
}

static void aux_put_Z(uint8_t **buf, size_t *len, size_t *cap,
                      const char tag[2], const char *s) {
    uint8_t header[3] = {(uint8_t)tag[0], (uint8_t)tag[1], 'Z'};
    buf_append(buf, len, cap, header, 3);
    size_t n = strlen(s);
    buf_append(buf, len, cap, s, n);
    uint8_t zero = 0;
    buf_append(buf, len, cap, &zero, 1);
}

/* BAM float ('f') aux field: tag[2] + 'f' + 4-byte little-endian float. */
static void aux_put_f(uint8_t **buf, size_t *len, size_t *cap,
                      const char tag[2], float v) {
    uint8_t tmp[3 + 4];
    tmp[0] = (uint8_t)tag[0]; tmp[1] = (uint8_t)tag[1]; tmp[2] = 'f';
    memcpy(tmp + 3, &v, 4);
    buf_append(buf, len, cap, tmp, 7);
}

/* Append the SA:Z text (NUL-terminated, no tag header) for a primary alignment
 * that has supplementary hits. */
static void build_sa_text(uint8_t **buf, size_t *len, size_t *cap,
                          const bntseq_t *bns,
                          const mem_aln_t *list, int n, int which) {
    // Format matches mem_aln2sam:
    //   chrom,pos+1,[+-],CIGAR,mapq,NM;chrom,...
    // The contig name is appended directly: its length is unbounded, so it
    // must not pass through the fixed `tmp`, whose snprintf results would then
    // be copied past its end. Every snprintf into `tmp` below has a bounded
    // numeric length.
    char tmp[64];
    for (int i = 0; i < n; ++i) {
        if (i == which || (list[i].flag & 0x100)) continue;
        const mem_aln_t *r = &list[i];
        const char *chrom = bns->anns[r->rid].name;
        buf_append(buf, len, cap, chrom, strlen(chrom));
        int m = snprintf(tmp, sizeof(tmp), ",%lld,%c,",
                         (long long)(r->pos + 1), r->is_rev ? '-' : '+');
        buf_append(buf, len, cap, tmp, m);
        for (int k = 0; k < r->n_cigar; ++k) {
            m = snprintf(tmp, sizeof(tmp), "%u%c",
                         r->cigar[k] >> 4, "MIDSH"[r->cigar[k] & 0xf]);
            buf_append(buf, len, cap, tmp, m);
        }
        m = snprintf(tmp, sizeof(tmp), ",%u,%u;", r->mapq, r->NM);
        buf_append(buf, len, cap, tmp, m);
    }
    uint8_t zero = 0;
    buf_append(buf, len, cap, &zero, 1);
}

/* Append the MC:Z text (mate CIGAR; NUL-terminated, no tag header) from a mate
 * mem_aln_t.
 *
 * `which` is the EMITTING record's index, not the mate's: upstream builds
 * MC:Z through the same `add_cigar` helper and hands it the record's own
 * `which` (bwamem.cpp:2583), so on a supplementary record the mate's clips
 * are rendered hard even though the mate itself is a primary. Mirrored here
 * deliberately -- it looks like an upstream quirk, and reproducing it is the
 * point. `m->is_alt` (not the record's) gates it, matching `add_cigar`'s use
 * of the mem_aln_t it was passed. */
static void build_mc_text(uint8_t **buf, size_t *len, size_t *cap,
                          const mem_opt_t *opt, const mem_aln_t *m, int which) {
    char tmp[32];
    for (int i = 0; i < m->n_cigar; ++i) {
        uint32_t cig = bwa_apply_clip_mode(m->cigar[i], opt, m->is_alt, which);
        int n = snprintf(tmp, sizeof(tmp), "%u%c",
                         cig >> 4, "MIDSH"[cig & 0xf]);
        buf_append(buf, len, cap, tmp, n);
    }
    uint8_t zero = 0;
    buf_append(buf, len, cap, &zero, 1);
}

/* Compute TLEN between two aligned mem_aln_t, mirroring upstream's SAM writer
 * (bwamem.cpp:2532-2539).
 *
 * Upstream takes each mate's OUTERMOST coordinate in its own direction -- the
 * leftmost base on the forward strand, the rightmost on the reverse -- and
 * signs the difference by which comes first, with a unit nudge that makes the
 * two mates' values exact negations of each other:
 *
 *     p0   = pos + (is_rev ? rlen - 1 : 0)
 *     tlen = -(p0 - p1 + sgn(p0 - p1))
 *
 * This is NOT the span between the pair's outer edges, which is what this
 * function used to return. The two agree on a canonical FR pair -- all any
 * fixture produced -- and diverge everywhere else. On a fully-overlapping
 * pair the span's `a_begin <= b_begin` sign test is true for BOTH mates, so
 * both got a positive TLEN and the pair violated the SAM requirement that
 * they sum to zero; on a contained pair the span reports the outer length
 * rather than the insert.
 *
 * `a_ref_len`/`b_ref_len` come from cigar_ref_len(), which sums the same M/D
 * opcodes as upstream's get_rlen() (bwamem.cpp:2768-2777).
 *
 * Returns 0 when either mate is unplaced, the two are on different contigs,
 * or either lacks a CIGAR. That last guard is upstream's `m->n_cigar == 0 ||
 * p->n_cigar == 0` and is what zeroes TLEN on a half-mapped pair: upstream
 * reaches it because its mate-coordinate rewrite sets `n_cigar = 0` on the
 * copied side, while this runs on the pre-rewrite mem_aln_t, where the same
 * pair is caught one line earlier by `rid < 0`. Both guards are kept so the
 * result matches on either path. */
static int32_t compute_tlen(const mem_aln_t *a, int a_ref_len,
                            const mem_aln_t *b, int b_ref_len) {
    if (!a || !b || a->rid < 0 || b->rid < 0 || a->rid != b->rid) return 0;
    if (a->n_cigar == 0 || b->n_cigar == 0) return 0;
    int64_t p0 = a->pos + (a->is_rev ? (int64_t) a_ref_len - 1 : 0);
    int64_t p1 = b->pos + (b->is_rev ? (int64_t) b_ref_len - 1 : 0);
    int64_t nudge = (p0 > p1) ? 1 : (p0 < p1 ? -1 : 0);
    return (int32_t) -(p0 - p1 + nudge);
}

/* Compute every field the packed-BAM record for one mem_aln_t carries, into
 * `f`. This is the single source of truth for record CONTENT: the packed path
 * (serialize_record) only turns `f` into bytes, and the structured-fields path
 * hands `f` to the caller as-is, so the two cannot disagree. `f`'s pointers
 * borrow `sc`'s cig_buf/txt_buf, the mem_aln_t's own buffers (MD, XA), the
 * global RG id and meth_build_xm's thread-local buffer; all stay valid until
 * the next record is computed on this scratch/thread.
 *
 * `opt`/`pac`/`is_r2` are used only for D3 (--meth) tag emission (is_r2: 0 =
 * R1/OT read, 1 = R2/OB read). */
static void compute_record_fields(ShimScratch *sc, BwaAlignedFields *f,
                                  const mem_opt_t *opt, const bntseq_t *bns,
                                  const uint8_t *pac, const bseq1_t *s,
                                  const mem_aln_t *p, int n_list,
                                  const mem_aln_t *list, int which,
                                  const mem_aln_t *m, int is_r2)
{
    memset(f, 0, sizeof(*f));
    int l_seq = s->l_seq;
    int n_cigar = p->n_cigar;
    int ref_len = cigar_ref_len(n_cigar, p->cigar);

    /* BAM field-width validation at the shim boundary. `l_read_name` is written
     * into a uint8 (max 255, NUL included) but the name is copied at full width,
     * and `n_cigar` into a uint16 (max 65535). A read name >= 255 bytes or a
     * CIGAR with >= 65536 ops is pathological but must fail loudly rather than
     * silently emit a truncated/corrupt record. This TU is compiled without
     * access to bwa_shim.cpp's shim_set_err and these emit helpers return void,
     * so mirror worker_alloc: xassert aborts with a message (survives NDEBUG),
     * matching the surrounding error handling here. Never fires on valid input,
     * so output stays byte-identical. Checked here, not in the serializer, so
     * the structured-fields path enforces the same limits. */
    xassert((int) strlen(s->name) + 1 <= 255, "read name too long for BAM (>= 255 bytes)");
    xassert(n_cigar <= 0xffff, "too many CIGAR ops for BAM (> 65535)");

    /* Upstream's mem_aln2sam prologue (bwamem.cpp:2457-2462), reproduced here
     * on LOCAL copies exactly as it does with ptmp/mtmp: for a half-mapped
     * pair the unmapped record is rewritten to carry its mate's placement, and
     * the mapped record's mate fields are rewritten to carry its own.
     *
     * Operating on copies is load-bearing, not stylistic. Mutating lists[k]
     * would leave the rewrite visible when the OTHER record of the pair is
     * emitted, and the second copy branch would then see a mate that already
     * looks mapped -- reporting both records as placed.
     *
     * This is more than byte-parity: giving the unmapped record its mate's
     * coordinates is what keeps the pair adjacent under coordinate sort, which
     * most downstream tools rely on.
     *
     * 0x4/0x8 are deliberately NOT recomputed here. The pair resolve already
     * set them from the ORIGINAL rids, which matches upstream's ordering: it
     * derives those two before the copy and 0x10/0x20 after it. */
    int32_t eff_rid    = p->rid;
    int64_t eff_pos    = p->pos;
    int     eff_is_rev = p->is_rev;
    int32_t mate_rid    = m ? m->rid : -1;
    int64_t mate_pos    = m ? m->pos : 0;
    int     mate_is_rev = m ? m->is_rev : 0;
    if (eff_rid < 0 && m && mate_rid >= 0) {          /* copy mate to alignment */
        eff_rid = mate_rid; eff_pos = mate_pos; eff_is_rev = mate_is_rev;
    }
    if (m && mate_rid < 0 && eff_rid >= 0) {          /* copy alignment to mate */
        mate_rid = eff_rid; mate_pos = eff_pos; mate_is_rev = eff_is_rev;
    }

    /* Emitted-SEQ range: hard-clipped ends are dropped, so the emitted
     * sequence is a subset of s->seq. Computed up front because the --meth
     * XM:Z tag below needs it.
     *
     * The predicate is upstream's (bwamem.cpp:2548-2551, 2562-2565): trim on a
     * supplementary (`which != 0`) unless -Y or an ALT alignment, keyed off the
     * clip opcode being S *or* H. This previously tested for H alone, which
     * `mem_reg2aln` never emits (it writes `clip<<4 | 3`, always S), so the
     * branch was dead and every supplementary carried a full-length SEQ.
     *
     * No is_rev swap here, despite upstream swapping which end qb/qe take
     * (bwamem.cpp:2563-2565): these are offsets into s->seq on the forward
     * path and offsets from the read's far end on the reverse path, because
     * the reverse emitter indexes `s->seq[l_seq - 1 - (seq_start + i)]`. That
     * already maps seq_start onto upstream's `qe` and seq_end onto its `qb`,
     * so adding a swap would double-correct. */
    int clip_hard = (n_cigar > 0 && which
                     && !(opt->flag & MEM_F_SOFTCLIP) && !p->is_alt);
    int seq_start = 0;
    int seq_end = l_seq;
    if (clip_hard) {
        uint32_t first_op = p->cigar[0] & 0xf;
        uint32_t last_op  = p->cigar[n_cigar - 1] & 0xf;
        if (first_op == 3 || first_op == 4) seq_start = (int)(p->cigar[0] >> 4);
        if (last_op  == 3 || last_op  == 4) seq_end  -= (int)(p->cigar[n_cigar - 1] >> 4);
    }
    int emit_len = seq_end - seq_start;
    if (emit_len < 0) emit_len = 0;
    /* A secondary record (raw 0x100: an -a secondary, not a MEM_F_NO_MULTI
     * split, which carries the internal 0x10000 until the flag remap below)
     * gets no SEQ/QUAL, exactly as both upstream writers decide it:
     * `emit_seq = !(p.flag & 0x100)` (bam_writer.cpp:346, meth_bam.cpp). */
    if (p->flag & 0x100) emit_len = 0;
    /* Report the window in as-sequenced read coordinates: the forward emitter
     * reads s->seq[seq_start + i] and the reverse one s->seq[l_seq - 1 -
     * (seq_start + i)], i.e. [l_seq - seq_start - emit_len, l_seq - seq_start)
     * read backwards. serialize_record walks exactly this window. */
    if (eff_is_rev) {
        f->query_start = l_seq - seq_start - emit_len;
        f->query_end   = l_seq - seq_start;
    } else {
        f->query_start = seq_start;
        f->query_end   = seq_start + emit_len;
    }

    /* Fixed-width fields. */
    f->tid  = eff_rid;
    f->pos  = (int32_t) eff_pos;
    f->mapq = (uint8_t) p->mapq;

    /* The packed FLAG is 16 bits, but the pair resolve marks MEM_F_NO_MULTI
     * split hits with upstream's internal 0x10000 (bit 16). Remap it to the BAM
     * secondary bit 0x100 here, exactly as mem_aln2sam does at write time.
     * Keeping the marker at 0x10000 internally (not 0x100) is deliberate: it
     * lets build_sa_text still list NO_MULTI splits in SA:Z (its skip test is
     * `flag & 0x100`), matching upstream. */
    /* 0x10/0x20 are recomputed from the POST-copy strands (upstream derives
     * them after the rewrite above), so an unmapped record placed at its
     * mate's coordinates inherits that mate's strand. Clearing first matters:
     * the pair resolve set them from the raw values. */
    uint32_t flag_raw = (p->flag & ~(uint32_t)(0x10 | 0x20));
    if (eff_is_rev)       flag_raw |= 0x10;
    if (m && mate_is_rev) flag_raw |= 0x20;
    f->flag = (uint16_t)((flag_raw & 0xffff) | ((flag_raw & 0x10000) ? 0x100 : 0));

    /* BAM bin over [pos, pos + rlen), computed as htslib's bam_set1 -- which
     * both upstream writers build records with -- computes it: rlen is the
     * CIGAR's reference span for a mapped record and 1 otherwise. So an
     * unmapped read placed at its mate's coordinates gets reg2bin(pos, pos+1),
     * not 4680; only an unplaced read (pos == -1) lands on 4680, which
     * reg2bin(-1, 0) yields through its arithmetic shift. */
    int bin_rlen = ((f->flag & 0x4) || ref_len == 0) ? 1 : ref_len;
    f->bin = reg2bin((int)eff_pos, (int)eff_pos + bin_rlen);

    f->next_tid = m ? mate_rid : -1;
    f->next_pos = m ? (int32_t) mate_pos : -1;
    f->tlen     = m ? compute_tlen(p, ref_len, m, cigar_ref_len(m->n_cigar, m->cigar)) : 0;

    /* Emitted CIGAR: clip rewrite, then bwa's MIDSH -> BAM opcode remap. */
    if ((size_t) n_cigar > sc->cig_cap) {
        size_t cap = sc->cig_cap ? sc->cig_cap : 64;
        while ((size_t) n_cigar > cap) cap *= 2;
        uint32_t *tmp = (uint32_t *) realloc(sc->cig_buf, cap * sizeof(uint32_t));
        xassert(tmp != NULL, "out of memory: cig_buf");
        sc->cig_buf = tmp; sc->cig_cap = cap;
    }
    for (int i = 0; i < n_cigar; ++i)
        sc->cig_buf[i] = bwa_cigar_to_bam(
            bwa_apply_clip_mode(p->cigar[i], opt, p->is_alt, which));
    f->cigar   = sc->cig_buf;
    f->n_cigar = (uint32_t) n_cigar;

    /* Aux, in the order serialize_record writes it. Every string not owned by
     * this call's bwa state (MC:Z, SA:Z, RG:Z, XM:Z) is copied into the
     * scratch's txt_buf, so each borrowed pointer handed to a field sink points
     * into memory only this exclusively-borrowed scratch owns: nothing the sink
     * can do (set the global read group, run meth_build_xm on its thread) can
     * rewrite a field it is still reading. txt_buf may realloc as it grows, so
     * record offsets and turn them into pointers only after the last append, at
     * the end of this function. */
    uint8_t *txt = sc->txt_buf;
    size_t txt_len = 0, txt_cap = sc->txt_cap;
    size_t mc_off = SIZE_MAX, sa_off = SIZE_MAX, rg_off = SIZE_MAX, xm_off = SIZE_MAX;
    if (p->n_cigar) {
        f->has_nm = 1; f->nm = p->NM;
        /* MD string is stored right after the CIGAR array in p->cigar. */
        f->md = (const char *)(p->cigar + p->n_cigar);
    }
    if (m && m->n_cigar) {
        mc_off = txt_len;
        build_mc_text(&txt, &txt_len, &txt_cap, opt, m, which);
    }
    /* MQ: the mate's MAPQ. Gated on `m` alone -- NOT `m->n_cigar` like MC:Z
     * above -- so it is emitted even when the mate is unmapped, matching
     * upstream exactly (bwamem.cpp:3484, bam_writer.cpp:395, meth_bam.cpp:585,
     * which all use `m && opt->compat->emit_mq`). Emitted between MC and AS,
     * because this crate's aux order mirrors upstream's write order.
     *
     * Not a bwa-mem3 invention: bwa emits MQ (bwamem.c:935, lh3/bwa#330) and
     * bwa-mem2 does not, having forked at 0.7.17 before that landed -- hence
     * the compat switch, whose default target (COMPAT_TARGET_OFF) sets
     * emit_mq = 1. */
    if (m && opt->compat != NULL && opt->compat->emit_mq) { f->has_mq = 1; f->mq = m->mapq; }
    if (p->score >= 0) { f->has_as = 1; f->score = p->score; }
    if (p->sub >= 0)   { f->has_xs = 1; f->sub = p->sub; }
    if (bwa_rg_id[0]) {
        rg_off = txt_len;
        buf_append(&txt, &txt_len, &txt_cap, bwa_rg_id, strlen(bwa_rg_id) + 1);
    }
    /* SA: if this is a primary (flag 0x100 not set) and other primary hits exist */
    if (!(p->flag & 0x100) && n_list > 1) {
        int i;
        for (i = 0; i < n_list; ++i)
            if (i != which && !(list[i].flag & 0x100)) break;
        if (i < n_list) {
            sa_off = txt_len;
            build_sa_text(&txt, &txt_len, &txt_cap, bns, list, n_list, which);
        }
    }
    /* pa:f — mirrors the vendored BAM writer (bam_writer.cpp:468-474) and
     * mem_aln2sam: emitted for a non-secondary record whose read has an ALT
     * competitor (alt_sc > 0), after SA:Z and before XA:Z. `append_bam_record`
     * previously omitted this entirely, so any ALT-adjacent read diverged from
     * `bwa-mem3 mem` -- invisible to the PhiX/tandem-repeat fixtures, which have
     * no `.alt` and so never set alt_sc. Value is the shared bwa_pa_tag_value,
     * so the --bam and SAM-text writers agree to the last float bit. */
    if (!(p->flag & 0x100) && p->alt_sc > 0) {
        f->has_pa = 1; f->pa = bwa_pa_tag_value(p->score, p->alt_sc);
    }
    if (p->XA) f->xa = p->XA;
    /* HN: total # of hits clustered with this primary under XA_drop_ratio
     * (set by mem_gen_alt above). -1 is upstream's "not computed" sentinel
     * (bwamem.h / bwamem.cpp's mem_reg2aln), so guard >= 0.
     *
     * v0.9.0 replaced the old `!meth_mode` condition. That condition existed
     * because upstream's --meth writer (meth_bam.cpp) never emitted HN, so
     * suppressing it was how we matched the CLI. As of v0.9.0 BOTH writers
     * emit it under a compat-target switch -- meth_bam.cpp:663 and
     * bam_writer.cpp:458 use the identical `p.HN >= 0 && opt->compat->emit_hn`
     * -- and COMPAT_TARGET_OFF (the default, "bwa-mem3's native output") sets
     * emit_hn = 1. Keeping the old gate left HN off every --meth record while
     * the CLI emitted it.
     *
     * `compat` is set by mem_opt_init, but it is a pointer on a struct callers
     * can memset, so it is NULL-checked rather than trusted. */
    if (p->HN >= 0 && opt->compat != NULL && opt->compat->emit_hn) { f->has_hn = 1; f->hn = p->HN; }

    /* D3 (--meth) Bismark tags. XR:Z (read conversion) on every record; XG:Z
     * (genome strand) and XM:Z (per-base methylation call) on mapped records.
     * XG / is_top_strand come from the winning hypothesis (p->meth_hypothesis:
     * OT/1 => top/CT, OB/0 or -1 => bottom/GA), NOT the 0x10 RC flag — matching
     * upstream meth_bam.cpp. XM is built from the ORIGINAL (unconverted) read
     * bases in emitted-SEQ orientation. */
    if (opt->meth_mode) {
        f->xr = is_r2 ? "GA" : "CT";
        if (p->rid >= 0) {
            int is_top = (p->meth_hypothesis >= 0 && (p->meth_hypothesis & 1)) ? 1 : 0;
            f->xg = is_top ? "CT" : "GA";
            if (s->meth_orig_seq && emit_len > 0 && n_cigar > 0) {
                char *seq_text = (char *) malloc((size_t)emit_len + 1);
                if (seq_text) {
                    /* `p->is_rev` (not `eff_is_rev`, as the primary SEQ RC
                     * uses) is correct here: this block is `p->rid >= 0`-gated,
                     * and the half-mapped is_rev copy only rewrites an *unmapped*
                     * record (rid < 0), so eff_is_rev == p->is_rev on every path
                     * that reaches the meth XM builder. */
                    for (int i = 0; i < emit_len; ++i) {
                        seq_text[i] = p->is_rev
                            ? ascii_complement(s->meth_orig_seq[l_seq - 1 - (seq_start + i)])
                            : s->meth_orig_seq[seq_start + i];
                    }
                    seq_text[emit_len] = '\0';
                    /* The XM:Z builder walks the emitted CIGAR (clip rewrite
                     * applied) alongside `seq_text`, which is the hard-clipped
                     * span, so the two must agree on whether the clip consumes
                     * query bases.
                     *
                     * v0.9.0 added the chemistry selector. It flips only the
                     * methylated/unmethylated polarity of the call, never the
                     * CpG/CHG/CHH context classification, and comes off
                     * mem_opt_t exactly as upstream's own writer sources it
                     * (meth_bam.cpp:526). mem_opt_init defaults it to
                     * METH_CHEM_EMSEQ, so this is unchanged behavior unless a
                     * caller sets it. The result aliases meth_build_xm's
                     * thread-local buffer, which the next call on this thread
                     * overwrites, so it is copied into txt_buf. */
                    const char *xm = meth_build_xm(bns, pac, p->rid, (int64_t)p->pos,
                                                   is_top, f->cigar, n_cigar, seq_text,
                                                   emit_len, (meth_chem_t) opt->meth_chem);
                    if (xm) {
                        xm_off = txt_len;
                        buf_append(&txt, &txt_len, &txt_cap, xm, strlen(xm) + 1);
                    }
                }
                free(seq_text);
            }
        }
    }

    /* The last txt_buf append is done: keep the (possibly grown) buffer for the
     * next record and turn the recorded offsets into pointers. */
    sc->txt_buf = txt; sc->txt_cap = txt_cap;
    if (mc_off != SIZE_MAX) f->mc = (const char *)(txt + mc_off);
    if (sa_off != SIZE_MAX) f->sa = (const char *)(txt + sa_off);
    if (rg_off != SIZE_MAX) f->rg = (const char *)(txt + rg_off);
    if (xm_off != SIZE_MAX) f->xm = (const char *)(txt + xm_off);
}

/* Serialize `f` into the packed BAM record BODY (no u32 block_size prefix) for
 * read `s`, in the scratch's reusable record buffer, and hand it to the packed
 * sink. Pure formatting: every decision was made by compute_record_fields. */
static void serialize_record(ShimEmit *e, size_t origin_idx, const bseq1_t *s,
                             const BwaAlignedFields *f)
{
    int l_read_name = (int) strlen(s->name) + 1;
    int n_cigar = (int) f->n_cigar;
    int emit_len = f->query_end - f->query_start;
    int is_rev = (f->flag & 0x10) != 0;

    /* Build aux first so we know its size. Reuse the scratch's aux buffer: it
     * is grown in place by buf_append and written back at the end of the record
     * so the next record reuses the larger allocation. */
    uint8_t *aux = e->sc->aux_buf;
    size_t aux_len = 0, aux_cap = e->sc->aux_cap;
    if (f->has_nm) aux_put_i(&aux, &aux_len, &aux_cap, "NM", f->nm);
    if (f->md)     aux_put_Z(&aux, &aux_len, &aux_cap, "MD", f->md);
    if (f->mc)     aux_put_Z(&aux, &aux_len, &aux_cap, "MC", f->mc);
    if (f->has_mq) aux_put_i(&aux, &aux_len, &aux_cap, "MQ", f->mq);
    if (f->has_as) aux_put_i(&aux, &aux_len, &aux_cap, "AS", f->score);
    if (f->has_xs) aux_put_i(&aux, &aux_len, &aux_cap, "XS", f->sub);
    if (f->rg)     aux_put_Z(&aux, &aux_len, &aux_cap, "RG", f->rg);
    if (f->sa)     aux_put_Z(&aux, &aux_len, &aux_cap, "SA", f->sa);
    if (f->has_pa) aux_put_f(&aux, &aux_len, &aux_cap, "pa", f->pa);
    if (f->xa)     aux_put_Z(&aux, &aux_len, &aux_cap, "XA", f->xa);
    if (f->has_hn) aux_put_i(&aux, &aux_len, &aux_cap, "HN", f->hn);
    if (f->xr)     aux_put_Z(&aux, &aux_len, &aux_cap, "XR", f->xr);
    if (f->xg)     aux_put_Z(&aux, &aux_len, &aux_cap, "XG", f->xg);
    if (f->xm)     aux_put_Z(&aux, &aux_len, &aux_cap, "XM", f->xm);

    size_t seq_packed = (size_t)(emit_len + 1) / 2;
    size_t block_size = (size_t)32 + l_read_name + 4 * (size_t)n_cigar
                      + seq_packed + (size_t)emit_len + aux_len;

    if (block_size > e->sc->rec_cap) {
        while (block_size > e->sc->rec_cap) e->sc->rec_cap *= 2;
        uint8_t *tmp = (uint8_t *) realloc(e->sc->rec_buf, e->sc->rec_cap);
        xassert(tmp != NULL, "out of memory: rec_buf");
        e->sc->rec_buf = tmp;
    }
    uint8_t *w = e->sc->rec_buf;

    memcpy(w, &f->tid, 4); w += 4;
    memcpy(w, &f->pos, 4); w += 4;
    *w++ = (uint8_t) l_read_name;
    *w++ = f->mapq;
    memcpy(w, &f->bin, 2); w += 2;
    uint16_t nc = (uint16_t) n_cigar;
    memcpy(w, &nc, 2); w += 2;
    memcpy(w, &f->flag, 2); w += 2;
    int32_t ls = emit_len;
    memcpy(w, &ls, 4); w += 4;
    memcpy(w, &f->next_tid, 4); w += 4;
    memcpy(w, &f->next_pos, 4); w += 4;
    memcpy(w, &f->tlen, 4); w += 4;

    memcpy(w, s->name, l_read_name - 1); w += l_read_name - 1;
    *w++ = 0;

    if (n_cigar) { memcpy(w, f->cigar, 4 * (size_t) n_cigar); w += 4 * (size_t) n_cigar; }

    /* SEQ/QUAL over the window [query_start, query_end) of the read, read
     * backwards and complemented on the reverse strand. Keys off the emitted
     * strand (0x10, i.e. `eff_is_rev`), not `p->is_rev`: for a half-mapped pair
     * upstream overwrites the unmapped read's `p->is_rev = m->is_rev`
     * (bwamem.cpp:4575) before mem_aln2sam writes SEQ, so the unmapped mate's
     * bases come out reverse-complemented to the mapped mate's strand. Keying
     * off the original `p->is_rev` (0 for an unmapped record) left the unmapped
     * mate forward while its 0x10 flag said reverse -- a SEQ/QUAL vs FLAG
     * contradiction, invisible to fixtures whose unmapped mates happen to be
     * forward. For a mapped read eff_is_rev == p->is_rev, so this is a no-op.
     *
     * s->seq holds bwa's 2-bit encoding (0-4); remap to BAM 4-bit nibbles
     * (1/2/4/8/15). Under --meth s->seq is the bisulfite-PROJECTED read the
     * seeds matched, so SEQ comes from the original bases instead, as the
     * upstream meth writer emits it (meth_seq_bam4). */
    if (s->meth_orig_seq) {
        const unsigned char *orig = (const unsigned char *) s->meth_orig_seq;
        for (int i = 0; i < emit_len; i += 2) {
            uint8_t hi = meth_seq_bam4(is_rev ? orig[f->query_end - 1 - i]
                                              : orig[f->query_start + i], is_rev);
            uint8_t lo = 0;
            if (i + 1 < emit_len)
                lo = meth_seq_bam4(is_rev ? orig[f->query_end - 2 - i]
                                          : orig[f->query_start + i + 1], is_rev);
            *w++ = (uint8_t)((hi << 4) | lo);
        }
    } else for (int i = 0; i < emit_len; i += 2) {
        uint8_t b0 = is_rev ? bwa2_complement((uint8_t) s->seq[f->query_end - 1 - i])
                            : (uint8_t) s->seq[f->query_start + i];
        uint8_t hi = bwa2_to_bam4[b0 & 7];
        uint8_t lo = 0;
        if (i + 1 < emit_len) {
            uint8_t b1 = is_rev ? bwa2_complement((uint8_t) s->seq[f->query_end - 2 - i])
                                : (uint8_t) s->seq[f->query_start + i + 1];
            lo = bwa2_to_bam4[b1 & 7];
        }
        *w++ = (uint8_t)((hi << 4) | lo);
    }

    /* qual: ASCII - 33, or 0xFF if missing. */
    if (s->qual) {
        for (int i = 0; i < emit_len; ++i) {
            char q = is_rev ? s->qual[f->query_end - 1 - i] : s->qual[f->query_start + i];
            *w++ = (uint8_t)(q - 33);
        }
    } else if (emit_len > 0) {
        memset(w, 0xFF, emit_len); w += emit_len;
    }

    if (aux_len) { memcpy(w, aux, aux_len); w += aux_len; }
    /* Retain the (possibly grown) aux buffer for the next record instead of
     * freeing it -- buf_append reallocs in place. */
    e->sc->aux_buf = aux; e->sc->aux_cap = aux_cap;

    assert((size_t)(w - e->sc->rec_buf) == block_size);
    e->sink(e->ctx, e->origin_kind, origin_idx, e->sc->rec_buf, block_size);
}

/* Emit one record for mem_aln_t `p` (the `which`-th of `list`, with mate
 * anchor `m`) of read `s` to `e`'s sink: its structured fields to a field
 * sink, or its packed BAM body to a packed sink. `is_r2` doubles as the field
 * sink's `mate` discriminant, and `which == 0` marks the read's first
 * (primary) record -- every caller emits a read's records in list order. */
static void append_bam_record(ShimEmit *e, size_t origin_idx,
                              const mem_opt_t *opt, const bntseq_t *bns,
                              const uint8_t *pac, const bseq1_t *s,
                              const mem_aln_t *p, int n_list,
                              const mem_aln_t *list, int which,
                              const mem_aln_t *m, int is_r2)
{
    BwaAlignedFields f;
    compute_record_fields(e->sc, &f, opt, bns, pac, s, p, n_list, list, which, m, is_r2);
    if (e->field_sink)
        e->field_sink(e->ctx, e->origin_kind, origin_idx, (uint8_t) is_r2, which == 0, &f);
    else
        serialize_record(e, origin_idx, s, &f);
}

/* ------------------ Phase 1: scratch + fused seed+extend ------------------ */

/* Reads per bwa-mem3 kernel batch (macro.h's BATCH_SIZE: 1024 on aarch64, 512
 * elsewhere): the chunk seed_extend_reads and pair_emit_pairs_chunked run each
 * kernel call on, and the work item the CLI's kt_for hands a worker. */
size_t shim_kernel_batch_size(void) {
    return (size_t) BATCH_SIZE;
}

/* The kernel thread slot `sc` runs the kernels in, or -1 for NULL. */
int shim_scratch_tid(const ShimScratch *sc) {
    return sc ? sc->tid : -1;
}

ShimScratch *shim_scratch_new(void) {
    ShimScratch *sc = (ShimScratch *) calloc(1, sizeof(ShimScratch));
    if (!sc) return nullptr;
    sc->tid = acquire_tid_slot();
    worker_alloc(sc->w, sc->tid);
    sc->rec_cap = 4096; sc->rec_buf = (uint8_t *) malloc(sc->rec_cap);
    sc->aux_cap = 1024; sc->aux_buf = (uint8_t *) malloc(sc->aux_cap);
    /* The emission path grows rec_buf/aux_buf from these caps and dereferences
     * them without re-checking, so a failed malloc here must not yield a
     * scratch with a non-zero cap but a NULL buffer. Fail the whole alloc. */
    if (!sc->rec_buf || !sc->aux_buf) {
        worker_free(sc->w, sc->tid);
        release_tid_slot(sc->tid);
        free(sc->rec_buf); free(sc->aux_buf);
        free(sc);
        return nullptr;
    }
    return sc;
}

void shim_scratch_free(ShimScratch *sc) {
    if (!sc) return;
    worker_free(sc->w, sc->tid);
    release_tid_slot(sc->tid);
    free(sc->rec_buf); free(sc->aux_buf);
    free(sc->cig_buf); free(sc->txt_buf);   /* lazily grown; NULL until first use */
    free(sc);
}

size_t shim_regs_n_pairs(const ShimRegs *r)   { return r ? r->n_pairs : 0; }
size_t shim_regs_n_singles(const ShimRegs *r) { return r ? r->n_singles : 0; }
size_t shim_regs_heap_bytes(const ShimRegs *r) {
    if (!r) return 0;
    size_t b = r->heap_bytes + (size_t)r->n_seqs * (sizeof(bseq1_t) + sizeof(mem_alnreg_v));
    for (int i = 0; i < r->n_seqs; ++i) b += (size_t)r->regs[i].m * sizeof(mem_alnreg_t);
    return b;
}

void shim_regs_free(ShimRegs *r) {
    if (!r) return;
    for (int i = 0; i < r->n_seqs; ++i) free(r->regs[i].a);
    free(r->regs);
    free_seqs(r->seqs, r->n_seqs);
    free(r);
}

/* Fused seed + SE-extend over reads `seqs[0, n)` / `regs[0, n)` in BATCH_SIZE
 * chunks, the shape of upstream's worker_bwt_aln (bwamem.cpp:2782-2786). `opt`
 * is the per-group copy (MEM_F_PE set for pairs, cleared for singles).
 * Chain/seed windows are chunk-local so the BATCH_SIZE scratch suffices; the
 * pre-0.9.0 tail `seedBufSz` shrink is dropped as upstream did ("output-dead
 * either way", worker_bwt). Shared by the legacy ShimRegs path and the resident
 * segments. */
static void seed_extend_reads(ShimScratch *sc, BwaShimIndex *idx, const mem_opt_t *opt,
                              bseq1_t *seqs, mem_alnreg_v *regs, int n)
{
    worker_t &w = sc->w;
    FMI_search *fmi = idx->fmi;
    w.opt = opt; w.fmi = fmi; w.seqs = seqs; w.regs = regs;
    w.n_processed = 0; w.pes = nullptr; w.nreads = n;
    w.meth_orig_bns = idx->meth_orig_bns;
    w.meth_orig_pac = idx->meth_orig_pac;
    /* No materialized reference on either path: a NULL ref_string selects
     * bns_get_seq_v2's pac-fetch, which mem_kernel2_core routes to the non-meth
     * pac (via `fmi`) or the meth ORIGINAL pac (`meth_orig_pac`) as appropriate. */
    w.meth_orig_ref_string = nullptr;
    w.ref_string = nullptr;
    for (int seq_id = 0; seq_id < n; seq_id += BATCH_SIZE) {
        int bs = n - seq_id;
        if (bs > BATCH_SIZE) bs = BATCH_SIZE;
        mem_kernel1_core(fmi, opt, seqs + seq_id, bs,
                         w.chain_scratch, w.seed_scratch, w.seed_scratch_size,
                         &w.mmc, sc->tid, idx->meth_orig_bns, idx->meth_orig_pac);
        mem_kernel2_core(fmi, opt, seqs + seq_id, regs + seq_id, bs,
                         w.chain_scratch, &w.mmc, w.ref_string, sc->tid,
                         idx->meth_orig_bns, idx->meth_orig_pac);
    }
}

ShimRegs *shim_seed_extend(void *idx_opaque, const mem_opt_t *opts, ShimScratch *sc,
                           const ShimReadBatch *batch)
{
    if (!idx_opaque || !opts || !sc || !batch) return nullptr;
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(idx_opaque);

    ShimRegs *r = (ShimRegs *) calloc(1, sizeof(ShimRegs));
    if (!r) return nullptr;
    r->n_pairs = batch->n_pairs;
    r->n_singles = batch->n_singles;
    r->n_seqs = (int)(2 * batch->n_pairs + batch->n_singles);
    r->meth_mode = opts->meth_mode;

    /* One contiguous seqs array: [0, 2*n_pairs) hold the pairs (R1/R2
     * interleaved), then [2*n_pairs, 2*n_pairs+n_singles) hold the singles.
     * copy_pairs_to_seqs owns its own array, so build the pairs into a temp and
     * shallow-copy the bseq1_t structs (their name/seq/qual heap now belongs to
     * r->seqs); copy_singles_to_seqs writes in place after them. */
    r->seqs = r->n_seqs > 0 ? (bseq1_t *) calloc((size_t)r->n_seqs, sizeof(bseq1_t)) : nullptr;
    if (!r->seqs && r->n_seqs > 0) { free(r); return nullptr; }
    if (batch->n_pairs > 0) {
        bseq1_t *pairs_only = copy_pairs_to_seqs(batch->pairs, batch->n_pairs, opts->meth_mode, &r->heap_bytes);
        if (!pairs_only) { free(r->seqs); free(r); return nullptr; }
        memcpy(r->seqs, pairs_only, 2 * batch->n_pairs * sizeof(bseq1_t));
        free(pairs_only);   /* shallow: the strings now belong to r->seqs */
    }
    if (copy_singles_to_seqs(r->seqs, 2 * batch->n_pairs, batch->singles, batch->n_singles,
                             opts->meth_mode, &r->heap_bytes) != 0) {
        free_seqs(r->seqs, r->n_seqs);   /* frees the built pairs + partial singles */
        free(r);
        return nullptr;
    }

    r->regs = r->n_seqs > 0 ? (mem_alnreg_v *) calloc((size_t)r->n_seqs, sizeof(mem_alnreg_v)) : nullptr;
    if (!r->regs && r->n_seqs > 0) { free_seqs(r->seqs, r->n_seqs); free(r); return nullptr; }

    /* Seed/extend each group under its own MEM_F_PE setting on a per-call copy
     * (fastmap.cpp:912-921): pairs with the flag SET, singles with it CLEARED. */
    mem_opt_t opt_pe = *opts; opt_pe.n_threads = 1; opt_pe.flag |=  MEM_F_PE;
    mem_opt_t opt_se = *opts; opt_se.n_threads = 1; opt_se.flag &= ~MEM_F_PE;
    if (r->n_pairs > 0)   seed_extend_reads(sc, idx, &opt_pe, r->seqs, r->regs, (int)(2 * r->n_pairs));
    if (r->n_singles > 0) seed_extend_reads(sc, idx, &opt_se, r->seqs + 2 * r->n_pairs,
                                            r->regs + 2 * r->n_pairs, (int)r->n_singles);
    return r;
}

/* ---- legacy phase-1 on top of the fused primitive ---- */

ShimSeeds *shim_seed_batch(void *idx_opaque, const mem_opt_t *opts,
                           const ShimReadPair *pairs, size_t n_pairs)
{
    ShimReadBatch b = { pairs, n_pairs, nullptr, 0 };
    ShimScratch *sc = shim_scratch_new();
    if (!sc) return nullptr;
    ShimRegs *r = shim_seed_extend(idx_opaque, opts, sc, &b);
    shim_scratch_free(sc);
    if (!r) return nullptr;
    ShimSeeds *s = (ShimSeeds *) calloc(1, sizeof(ShimSeeds));
    if (!s) { shim_regs_free(r); return nullptr; }
    s->regs = r;
    return s;
}

void shim_seeds_free(ShimSeeds *s) {
    if (!s) return;
    shim_regs_free(s->regs);
    free(s);
}

/* ------------------ Phase 3: cohort pestat + pair/emit ------------------ */

/* Defined after emit_resolved_pair; forward-declared so shim_pair_emit can
 * drive it after the pair loop. */
static void single_and_emit(ShimEmit *e, size_t origin_idx, uint64_t id,
                            const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                            bseq1_t *s, mem_alnreg_v *a);

/* Defined below; forward-declared so the batched mate-rescue path in
 * shim_pair_emit can emit each resolved pair. */
static void emit_resolved_pair(ShimEmit *e, size_t origin_idx,
                               const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                               bseq1_t *s, mem_alnreg_v *a, const mem_pestat_t pes[4],
                               const int n_pri[2], const int z[2], const int q_se[2],
                               int extra_flag, int paired);

/* mem_pestat over the PE reads of several batches, exactly as the CLI runs it
 * once per -K cohort (mem_process_seqs, bwamem.cpp:3030-3050). Only the
 * 24-byte mem_alnreg_v headers are gathered; the alnreg payloads stay put.
 * mem_pestat sorts internally, so gather order only needs to keep each pair's
 * two reads adjacent.
 *
 * MIXED cohorts (pairs + singles): only the PAIR regs (2*n_pairs per batch) are
 * gathered — singles are excluded. This matches the CLI's `-p` path exactly:
 * bseq_classify splits the cohort into an SE and a PE group, and only the PE
 * group's mem_process_seqs call estimates insert size (over the PE reads' regs,
 * bwamem.cpp:3030-3042); the SE group runs with MEM_F_PE cleared and never
 * touches pestat (fastmap.cpp:919-946). Singles have no mate and so contribute
 * nothing to the insert-size distribution either way. Verified against
 * `bwa-mem3 mem -p` on a mixed interleaved input (three_phase_cli_parity.rs). */
/* mem_pestat over `total` gathered PE-read alnreg headers (R1/R2 interleaved),
 * with the per-call PE opt copy. Shared by the legacy cohort gather below and
 * the resident cohort's. */
static void pestat_over_headers(BwaShimIndex *idx, const mem_opt_t *opts,
                                const mem_alnreg_v *cat, size_t total, mem_pestat_t *out)
{
    const bntseq_t *bns = idx->meth_orig_bns ? idx->meth_orig_bns : idx->fmi->idx->bns;
    mem_opt_t opt_pe = *opts; opt_pe.n_threads = 1; opt_pe.flag |= MEM_F_PE;
    mem_pestat(&opt_pe, bns->l_pac, (int)total, cat, out);
}

int shim_pestat_cohort(void *idx_opaque, const mem_opt_t *opts,
                       const ShimRegs *const *regs, size_t n_regs, mem_pestat_t *out)
{
    if (!idx_opaque || !opts || (n_regs > 0 && !regs) || !out) return -1;
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(idx_opaque);
    size_t total = 0;
    for (size_t k = 0; k < n_regs; ++k) total += 2 * regs[k]->n_pairs;
    mem_alnreg_v *cat = total ? (mem_alnreg_v *) malloc(total * sizeof(mem_alnreg_v)) : nullptr;
    if (total && !cat) return -1;
    size_t at = 0;
    for (size_t k = 0; k < n_regs; ++k) {
        size_t n = 2 * regs[k]->n_pairs;
        if (n) memcpy(cat + at, regs[k]->regs, n * sizeof(mem_alnreg_v));
        at += n;
    }
    pestat_over_headers(idx, opts, cat, total, out);
    free(cat);
    return 0;
}

/* The CLI's look-ahead hint for the CIGAR-regeneration reference window of an
 * alignment starting at `rb`: a verbatim copy of mem_prefetch_cigar_ref, which
 * is static in bwamem.cpp (see its comment there for the strand split and the
 * locality choice). A pure hint -> byte-identical output. */
static inline void shim_prefetch_cigar_ref(const bntseq_t *bns, const uint8_t *pac, int64_t rb)
{
    if (rb < 0) return;
    int64_t l_pac = bns->l_pac;
    int64_t pf = (rb < l_pac) ? (rb >> 2) : (((l_pac << 1) - 1 - rb) >> 2);
    __builtin_prefetch(&pac[pf], 0, 0);
}

/* Free one read's alignment regions and leave an empty vector behind, so
 * every later free (the legacy shim_regs_free, a resident segment's teardown)
 * is a no-op on it. */
static void release_regs(mem_alnreg_v *v) {
    free(v->a);
    v->a = nullptr;
    v->n = v->m = 0;
}

/* Pair, mate-rescue and emit `n_pairs` pairs (reads `seqs[0, 2*n_pairs)`, R1/R2
 * interleaved) — the CLI's worker_sam path on AVX2/AVX-512/NEON
 * (bwamem.cpp:2826-2879): gather every pair's rescue jobs, run them through the
 * SIMD kswv kernel once, then resolve + emit per pair. Pair `i` gets read
 * ordinal `first_pair_id + i` and sink origin `origin_base + i`. `opt_pe` is the
 * per-call PE copy. Kernels run in the scratch's own tid slot. Shared by the legacy
 * and resident paths.
 *
 * Batched mate rescue is the only mate-rescue path upstream ships: the scalar
 * mem_sam_pe / mem_pair_resolve pairing path (and the BWAMEM_BATCHED_MATESW /
 * DISABLE_BATCHED_MATESW gate) were removed in fg-labs/bwa-mem3#513, which
 * requires AVX2+ on x86 (NEON on arm).
 *
 * Chunk the batched mate-rescue exactly as the CLI's kt_for dispatches
 * worker_sam: BATCH_SIZE reads (= BATCH_SIZE/2 pairs) per work item
 * (kthread.cpp:109-118, bwamem.cpp worker_sam ~L2788). The pre-loop's running
 * pcnt/gcnt/maxRefLen/maxQerLen and, with them, the kswv seqBuf ref-window
 * offset (SeqPair.idr, an int32 derived from the monotonic pcnt) reset at every
 * chunk boundary, so the offset can never exceed int32 no matter how large the
 * -K cohort or how long the reads (seqbuf_grow_capacity's
 * SEQBUF_CAPACITY_OVERFLOW → seqbuf_capacity_fatal exit(), which would abort
 * across the FFI boundary). Without the chunk loop a single running pcnt
 * spanned all the pairs, exactly the unbounded-offset the CLI avoids by
 * resetting per BATCH_SIZE.
 *
 * Chunking is purely internal kswv/seqBuf batch bookkeeping. The per-pair
 * rescue result is independent of how pairs are grouped, and the GLOBAL
 * read-ordinal id and the emit order (index i, ascending) are unchanged — only
 * pcnt/gcnt/maxRefLen/maxQerLen/aln/myaln reset per chunk. So output is
 * byte-identical to a single batch and to `bwa-mem3 mem -t 1 -p`, however the
 * caller splits the pairs. For n_pairs <= BATCH_SIZE/2 this is a single
 * iteration. */
static void pair_emit_pairs_chunked(ShimEmit *e, const mem_opt_t *opt_pe,
                                    const bntseq_t *bns, const uint8_t *pac,
                                    const mem_pestat_t *pestat, uint64_t first_pair_id,
                                    bseq1_t *seqs, mem_alnreg_v *regs, size_t n_pairs,
                                    size_t origin_base)
{
    worker_t &w = e->sc->w;
    const int tid = e->sc->tid;
    const size_t pairs_per_chunk = (size_t)BATCH_SIZE / 2;
    for (size_t chunk_start = 0; chunk_start < n_pairs; chunk_start += pairs_per_chunk) {
        size_t chunk_end = chunk_start + pairs_per_chunk;
        if (chunk_end > n_pairs) chunk_end = n_pairs;

        int32_t maxRefLen = 0, maxQerLen = 0, gcnt = 0;
        int64_t pcnt = 0;
        for (size_t i = chunk_start; i < chunk_end; ++i)
            mem_sam_pe_batch_pre(opt_pe, bns, pac, pestat, first_pair_id + (uint64_t)i,
                                 seqs + 2*i, regs + 2*i, &w.mmc, pcnt, gcnt,
                                 maxRefLen, maxQerLen, tid);
        int64_t pcnt8 = sort_classify(&w.mmc, pcnt, tid);
        kswr_t *aln = (kswr_t *) _mm_malloc((pcnt + SIMD_WIDTH8) * sizeof(kswr_t), 64);
        xassert(aln != NULL, "out of memory: aln");
        mem_sam_pe_batch(opt_pe, &w.mmc, pcnt, pcnt8, aln, maxRefLen, maxQerLen, tid);
        gcnt = 0;
        kswr_t *myaln = aln;
        for (size_t i = chunk_start; i < chunk_end; ++i) {
            /* Look-ahead, as worker_sam does (bwamem.cpp:3945-3958): prefetch
             * the next pair's two emit windows while this pair resolves. */
            if (i + 1 < chunk_end) {
                if (regs[2*i + 2].n > 0) shim_prefetch_cigar_ref(bns, pac, regs[2*i + 2].a[0].rb);
                if (regs[2*i + 3].n > 0) shim_prefetch_cigar_ref(bns, pac, regs[2*i + 3].a[0].rb);
            }
            int n_pri[2], z[2], q_se[2], extra_flag, paired;
            mem_pair_resolve_batch_post(opt_pe, bns, pac, pestat, first_pair_id + (uint64_t)i,
                                        seqs + 2*i, regs + 2*i, &myaln, &w.mmc, gcnt, tid,
                                        n_pri, z, q_se, &extra_flag, &paired);
            emit_resolved_pair(e, origin_base + i, opt_pe, bns, pac, seqs + 2*i, regs + 2*i,
                               pestat, n_pri, z, q_se, extra_flag, paired);
            /* Nothing reads a pair's regions once it is emitted: release them
             * now, while they are still hot, as worker_sam does right after
             * mem_sam_pe_batch_post (bwamem.cpp:3970-3971). */
            release_regs(&regs[2*i]);
            release_regs(&regs[2*i + 1]);
        }
        _mm_free(aln);
    }
}

/* SE resolve + emit of `n` singles: single `i` gets read ordinal
 * `first_single_id + i` and sink origin `origin_base + i`, under the per-call
 * SE opt copy `opt_se`. Shared by the legacy and resident paths. */
static void emit_singles(ShimEmit *e, const mem_opt_t *opt_se, const bntseq_t *bns,
                         const uint8_t *pac, uint64_t first_single_id, bseq1_t *seqs,
                         mem_alnreg_v *regs, size_t n, size_t origin_base)
{
    for (size_t i = 0; i < n; ++i) {
        /* Look-ahead, as worker_sam's SE loop does (bwamem.cpp:3984-3985). */
        if (i + 1 < n && regs[i + 1].n > 0)
            shim_prefetch_cigar_ref(bns, pac, regs[i + 1].a[0].rb);
        single_and_emit(e, origin_base + i, first_single_id + (uint64_t)i,
                        opt_se, bns, pac, seqs + i, regs + i);
        release_regs(&regs[i]);
    }
}

/* Phase 3: pairing + mate rescue + primary marking + emission for one batch.
 * Consumes `r` (freed on every return path). Returns -1 on null args, -2 when
 * the batch has pairs but no pestat (the caller maps -2 to a "pestat required"
 * error). Records stream to `sink` as packed BAM bodies in input order. */
int shim_pair_emit(void *idx_opaque, const mem_opt_t *opts, ShimScratch *sc, ShimRegs *r,
                   const mem_pestat_t *pestat, ShimIdBases ids,
                   ShimRecordSinkFn sink, void *ctx)
{
    if (!r) return -1;
    if (!idx_opaque || !opts || !sc || !sink) { shim_regs_free(r); return -1; }
    if (r->n_pairs > 0 && !pestat) { shim_regs_free(r); return -2; }   /* caller maps to "pestat required" */
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(idx_opaque);
    const bntseq_t *bns = idx->meth_orig_bns ? idx->meth_orig_bns : idx->fmi->idx->bns;
    const uint8_t  *pac = idx->meth_orig_pac ? idx->meth_orig_pac : idx->fmi->idx->pac;

    mem_opt_t opt_pe = *opts; opt_pe.n_threads = 1; opt_pe.flag |= MEM_F_PE;
    ShimEmit e = { sc, sink, ctx, 0u /* BWA_ORIGIN_PAIR */, nullptr };
    pair_emit_pairs_chunked(&e, &opt_pe, bns, pac, pestat, ids.first_pair_id,
                            r->seqs, r->regs, r->n_pairs, 0);
    /* Singles after pairs: SE group emitted with MEM_F_PE cleared, ids from
     * first_single_id, one origin_idx per single (index into batch->singles).
     * BWA_ORIGIN_SINGLE (1) distinguishes them from pairs at the sink. */
    mem_opt_t opt_se = *opts; opt_se.n_threads = 1; opt_se.flag &= ~MEM_F_PE;
    e.origin_kind = 1u /* BWA_ORIGIN_SINGLE */;
    size_t k = 2 * r->n_pairs;
    emit_singles(&e, &opt_se, bns, pac, ids.first_single_id, r->seqs + k, r->regs + k,
                 r->n_singles, 0);
    shim_regs_free(r);
    return 0;
}

/* ------------------ Phase 2: extend_batch ------------------ */

/* Defined in the output-accessors section below; forward-declared so
 * shim_extend_batch can release a partially-built output on error. */
void shim_align_out_free(ShimAlignOutput *out);

static ShimAlignOutput *alloc_align_output(size_t n_pairs) {
    ShimAlignOutput *out = (ShimAlignOutput *) calloc(1, sizeof(ShimAlignOutput));
    if (!out) return nullptr;
    out->cap = n_pairs > 0 ? n_pairs * 4 : 64;
    out->rec_off  = (size_t *) malloc(out->cap * sizeof(size_t));
    out->rec_len  = (size_t *) malloc(out->cap * sizeof(size_t));
    out->pair_idx = (size_t *) malloc(out->cap * sizeof(size_t));
    out->buf_cap = 4096;
    out->buf = (uint8_t *) malloc(out->buf_cap);
    /* Any partial failure unwinds the whole output (shim_align_out_free is
     * NULL-safe on the fields not yet allocated). The legacy callers check the
     * NULL return. */
    if (!out->rec_off || !out->rec_len || !out->pair_idx || !out->buf) {
        shim_align_out_free(out);
        return nullptr;
    }
    return out;
}

/* Emission half of a resolved pair: everything after the pair resolve
 * (mem_pair_resolve_batch_post) has decided primaries/flags. On the paired
 * branch extra_flag already includes 0x2 (if the paired alignment was
 * preferred); on the no_pairing branch this half ORs it in after running
 * mem_infer_dir itself, matching mem_sam_pe's no_pairing block. */
static void emit_resolved_pair(ShimEmit *e, size_t origin_idx,
                               const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                               bseq1_t *s, mem_alnreg_v *a, const mem_pestat_t pes[4],
                               const int n_pri[2], const int z[2], const int q_se[2],
                               int extra_flag, int paired)
{
    /* PAIRED branch: emit exactly what upstream mem_sam_pe_batch_post's paired
     * block emits (bwamem_pair.cpp:1279-1310) -- the paired-selected primary
     * a[k].a[z[k]] UNCONDITIONALLY (no opt->T gate), plus at most one ALT
     * supplementary a[k].a[n_pri[k]] gated on score >= T && secondary < 0 &&
     * is_alt. The general mem_reg2sam-style emit loop below this block is the
     * NO-PAIRING / SE emission path; applying it to the paired case (as this
     * shim previously did) wrongly (a) dropped a paired primary whose score is
     * below opt->T -- e.g. a rescued partial mate, 28M122S ~= score 28 < T=30 --
     * to an unmapped record, and (b) emitted a different supplementary/XA set
     * than upstream on ALT-hit reads (its per-region filter + mate-of-first-
     * region logic diverges from upstream's z[k]+n_pri[k] structure). Both
     * surface only against a real reference (partial rescues + a `.alt`), which
     * no PhiX/tandem-repeat fixture exercises. Mirror upstream verbatim. */
    if (paired) {
        char **XA[2] = { nullptr, nullptr };
        int   *HN[2] = { nullptr, nullptr };
        if (!(opt->flag & MEM_F_ALL)) {
            for (int k = 0; k < 2; ++k)
                XA[k] = mem_gen_alt(opt, bns, pac, &a[k], s[k].l_seq, s[k].seq,
                                    &HN[k], s[k].meth_orig_seq);
        }
        mem_aln_t h[2]; memset(h, 0, sizeof(h));
        mem_aln_t g[2]; memset(g, 0, sizeof(g));
        mem_aln_t aa[2][2]; memset(aa, 0, sizeof(aa));
        int n_aa[2] = { 0, 0 };
        for (int k = 0; k < 2; ++k) {
            /* Paired primary: a[k].a[z[k]], emitted unconditionally. */
            h[k] = mem_reg2aln(opt, bns, pac, s[k].l_seq, s[k].seq,
                               &a[k].a[z[k]], s[k].meth_orig_seq);
            h[k].mapq = q_se[k];
            h[k].flag |= (0x40 << k) | extra_flag;   /* 0x1 paired (+0x2) + first/last */
            h[k].XA = XA[k] ? XA[k][z[k]] : nullptr;
            h[k].HN = HN[k] ? HN[k][z[k]] : -1;
            aa[k][n_aa[k]++] = h[k];
            if (n_pri[k] < (int) a[k].n) {   /* the read has ALT hits */
                mem_alnreg_t *p = &a[k].a[n_pri[k]];
                if (p->score < opt->T || p->secondary >= 0 || !p->is_alt) continue;
                g[k] = mem_reg2aln(opt, bns, pac, s[k].l_seq, s[k].seq, p,
                                   s[k].meth_orig_seq);
                g[k].flag |= 0x800 | (0x40 << k) | extra_flag;   /* supplementary */
                g[k].XA = XA[k] ? XA[k][n_pri[k]] : nullptr;
                g[k].HN = HN[k] ? HN[k][n_pri[k]] : -1;
                if (opt->supp_rep_hard_cap > 0 && p->chain_n_hits >= opt->supp_rep_hard_cap)
                    g[k].mapq = 0;   /* fg-labs: force repetitive-supp MAPQ to 0 */
                aa[k][n_aa[k]++] = g[k];
            }
        }
        /* Emit each side's records with the OTHER side's paired primary as the
         * mate anchor (upstream: aa[0] with &h[1], aa[1] with &h[0]).
         * append_bam_record derives 0x10/0x20 from the strands and, for a
         * half-mapped pair, the mate-coordinate copy; both mates are mapped on
         * the paired branch, so 0x4/0x8 stay clear. */
        for (int k = 0; k < 2; ++k) {
            mem_aln_t *mate = &h[!k];
            for (int j = 0; j < n_aa[k]; ++j)
                append_bam_record(e, origin_idx, opt, bns, pac, &s[k],
                                  &aa[k][j], n_aa[k], aa[k], j, mate, k);
        }
        /* aa[k][*] are shallow copies sharing h[k]/g[k]'s cigar buffers, so the
         * cigars are freed once via h[k]/g[k] here (never via aa). */
        for (int k = 0; k < 2; ++k) {
            free(h[k].cigar);
            free(g[k].cigar);
            free(HN[k]);
            if (XA[k]) {
                for (int j = 0; j < (int) a[k].n; ++j) free(XA[k][j]);
                free(XA[k]);
            }
        }
        return;
    }

    /* No-pairing proper-pair selection, computed exactly as upstream
     * mem_sam_pe's no_pairing block does (bwamem_pair.cpp:1651-1658): the
     * emitted primary of side k is region 0 if it clears T, else the n_pri[k]
     * region if THAT clears T, else none (-1). This drives mem_proper_pair_extra_flag
     * on the no-pairing branch below; `z` is documented undefined when paired==0,
     * so it must not be used there. Under the default (proper_pair_from_emitted
     * off) mem_proper_pair_extra_flag ignores `which` entirely and keys off
     * region 0, so passing the correct `which` here is byte-identical to the old
     * `z` on every input that does not set --proper-pair-from-emitted; it fixes
     * the latent garbage-index read only that option would expose. `h_rid[k]` is
     * the rid of the selected region (== upstream's h[i].rid), captured from
     * lists[k] below and used to gate the flag exactly as upstream does. */
    int which[2] = { -1, -1 };
    int h_rid[2] = { -1, -1 };
    for (int k = 0; k < 2; ++k) {
        if (a[k].n) {
            if (a[k].a[0].score >= opt->T) which[k] = 0;
            else if (n_pri[k] < (int)a[k].n && a[k].a[n_pri[k]].score >= opt->T)
                which[k] = n_pri[k];
        }
    }
    /* Build a synthetic unmapped mem_aln_t in `dst[0]` (a single-record list).
     * upstream's own unmapped record (mem_reg2aln, bwamem.cpp:2632-2645)
     * comes from `memset(&a, 0, sizeof(mem_aln_t))` followed by an
     * immediate early return for `ar == 0`, EXCEPT for two fields it
     * explicitly overwrites first: `a.HN = -1` and `a.meth_hypothesis = -1`
     * (bwamem.cpp:2640-2641) — both writer-side sentinels, since 0 is a
     * valid non-sentinel value for either. Every other field, including
     * `score` and `sub`, is left at the memset default of 0 (not -1). Two
     * consequences that must be mirrored exactly, not inferred from "what
     * looks like a sane default":
     *   - HN and meth_hypothesis here MUST be -1, matching upstream's
     *     explicit overwrite. Leaving them at calloc's 0 previously leaked
     *     HN:i:0 onto unmapped records, which the CLI never emits (see
     *     append_bam_record's `p->HN >= 0` guard below); meth_hypothesis
     *     has no emission effect today (XG/XM are gated on `p->rid >= 0`,
     *     never true here) but is set for the same reason so this lambda
     *     doesn't silently drift from upstream again.
     *   - score and sub here MUST be 0, matching upstream's memset (NOT the
     *     -1 sentinel append_bam_record's `p->score >= 0` / `p->sub >= 0`
     *     guards use to suppress AS:i/XS:i elsewhere): the CLI's own
     *     unmapped records carry AS:i:0 and XS:i:0 for exactly this reason,
     *     confirmed against a real `bwa-mem3 mem` run on a read that fails
     *     to map (see the `unmapped0` fixture in `cli_parity_xa.rs`). */
    auto fill_unmapped = [](mem_aln_t *dst) {
        dst[0].rid = -1;
        dst[0].pos = -1;
        dst[0].flag = 4;
        dst[0].n_cigar = 0;
        dst[0].cigar = nullptr;
        dst[0].mapq = 0;
        dst[0].NM = 0;
        dst[0].score = 0;
        dst[0].sub = 0;
        dst[0].HN = -1;
        dst[0].meth_hypothesis = -1;
    };

    /* Convert each side's regions to mem_aln_t arrays for emission and mark
     * which ones actually get emitted, mirroring mem_reg2sam: secondary
     * alignments are folded into the primary's XA:Z tag (via mem_gen_alt, which
     * requires the mem_mark_primary_se that mem_pair_resolve_batch_post ran) rather than
     * emitted as their own records, and sub-threshold regions are dropped.
     * Without this the shim emits every surviving alnreg — over-emitting
     * secondaries on multi-mapping reads (including --meth's collapsed-scoring
     * hits) that the CLI folds into XA. `lists[k]` stays 1:1 with `a[k]` so the
     * pairing indices (z[k]) and mate/SA logic below are unaffected; `emit[k]`
     * gates the final append. */
    mem_aln_t *lists[2] = {nullptr, nullptr};
    int n_lists[2] = {0, 0};
    bool *emit[2] = {nullptr, nullptr};
    char **XA[2] = {nullptr, nullptr};
    int   *HN[2] = {nullptr, nullptr};
    for (int k = 0; k < 2; ++k) {
        if (a[k].n == 0) {
            lists[k] = (mem_aln_t *) calloc(1, sizeof(mem_aln_t));
            xassert(lists[k] != NULL, "out of memory: lists[k]");
            fill_unmapped(lists[k]);
            n_lists[k] = 1;
            emit[k] = (bool *) malloc(sizeof(bool));
            xassert(emit[k] != NULL, "out of memory: emit[k]");
            emit[k][0] = true;
            continue;
        }
        if (!(opt->flag & MEM_F_ALL))
            /* v0.9.0 added the trailing `meth_orig_query`, which mem_gen_alt
             * forwards to mem_reg2aln for every XA sub-entry so they regenerate
             * under the same NM/MD policy as the primary record. It DEFAULTS TO
             * NULL, so omitting it compiles and silently keeps the legacy
             * regen: each XA sub-entry is then scored against the PROJECTED
             * (C->T / G->A) read, and its NM counts every bisulfite conversion.
             * That showed up as `XA:Z:phix,+1201,100M,20;` where the CLI emits
             * `...,0;`. Pass it exactly as upstream does (bwamem.cpp:3237);
             * s[k].meth_orig_seq is NULL outside --meth, which selects the same
             * legacy behavior the non-meth path always had. */
            XA[k] = mem_gen_alt(opt, bns, pac, &a[k], s[k].l_seq, s[k].seq, &HN[k],
                                s[k].meth_orig_seq);
        lists[k] = (mem_aln_t *) calloc(a[k].n, sizeof(mem_aln_t));
        xassert(lists[k] != NULL, "out of memory: lists[k]");
        emit[k] = (bool *) malloc(a[k].n * sizeof(bool));
        xassert(emit[k] != NULL, "out of memory: emit[k]");
        n_lists[k] = (int)a[k].n;
        int n_emit = 0;
        for (int j = 0; j < (int)a[k].n; ++j) {
            const mem_alnreg_t *ar = &a[k].a[j];
            /* mem_reg2sam emit filter. */
            bool e = ar->score >= opt->T;
            if (e && ar->secondary >= 0 && (ar->is_alt || !(opt->flag & MEM_F_ALL)))
                e = false;
            if (e && ar->secondary >= 0 && ar->secondary < INT_MAX
                && ar->score < a[k].a[ar->secondary].score * opt->drop_ratio)
                e = false;
            /* Convert a region only when something reads it: an emitted
             * record, region 0 (the mate-flag loop below reads lists[!k][0]),
             * or which[k] (h_rid and the mate anchor). mem_reg2aln regenerates
             * the CIGAR with a global alignment, and a repetitive read reaching
             * this branch can carry dozens of secondaries that are only ever
             * folded into XA:Z, so converting all of them made this loop the
             * dominant ksw_global2 caller. Upstream mem_reg2sam likewise
             * converts only what it emits. The skipped entries stay zeroed
             * from the calloc above (null cigar, so the frees below are
             * no-ops) and are never emitted or used as an anchor. */
            if (e || j == 0 || j == which[k]) {
                /* D3 (--meth): pass the ORIGINAL read bases so mem_reg2aln
                 * regenerates CIGAR/NM/MD against the original ref (NULL and a
                 * no-op outside --meth). */
                lists[k][j] = mem_reg2aln(opt, bns, pac, s[k].l_seq, s[k].seq,
                                          ar, s[k].meth_orig_seq);
                lists[k][j].XA = XA[k] ? XA[k][j] : nullptr;
                lists[k][j].HN = HN[k] ? HN[k][j] : -1;
            }
            /* (The paired-branch secondary_all fold that once lived here was
             * removed with the rest of the dead paired code: this tail now runs
             * only when paired==0, where `z[]` is undefined and must not be
             * read. The paired case is handled by the early-return block above.) */
            emit[k][j] = e;
            if (e) ++n_emit;
        }
        /* Capture the selected region's rid while lists[k] is still 1:1 with
         * a[k] (the n_emit==0 branch below replaces lists[k] with a synthetic
         * single unmapped, so lists[k][which[k]] would be out of bounds after).
         * lists[k][which[k]] is mem_reg2aln(a[k].a[which[k]]) with the same
         * args upstream uses, so this rid == upstream's h[k].rid exactly. */
        if (which[k] >= 0) h_rid[k] = lists[k][which[k]].rid;
        /* No region cleared the threshold → emit a single unmapped record
         * (mem_reg2sam's `aa.n == 0` branch). */
        if (n_emit == 0) {
            for (int j = 0; j < n_lists[k]; ++j) free(lists[k][j].cigar);
            free(lists[k]);
            free(emit[k]);
            if (XA[k]) {
                for (int j = 0; j < (int)a[k].n; ++j) free(XA[k][j]);
                free(XA[k]);
                XA[k] = nullptr;
            }
            free(HN[k]);
            HN[k] = nullptr;
            lists[k] = (mem_aln_t *) calloc(1, sizeof(mem_aln_t));
            xassert(lists[k] != NULL, "out of memory: lists[k]");
            fill_unmapped(lists[k]);
            n_lists[k] = 1;
            emit[k] = (bool *) malloc(sizeof(bool));
            xassert(emit[k] != NULL, "out of memory: emit[k]");
            emit[k][0] = true;
        }
    }

    /* Set paired-end flags per mem_aln2sam's flag-propagation rules. */
    for (int k = 0; k < 2; ++k) {
        for (int j = 0; j < n_lists[k]; ++j) {
            mem_aln_t *p = &lists[k][j];
            if (opt->flag & MEM_F_PE) {
                p->flag |= 0x1;                          /* paired */
                p->flag |= (k == 0) ? 0x40 : 0x80;       /* first / last in pair */
                if (n_lists[!k] > 0 && lists[!k][0].rid >= 0) {
                    if (lists[!k][0].is_rev) p->flag |= 0x20;  /* mate reverse */
                } else {
                    p->flag |= 0x8;                            /* mate unmapped */
                }
                if (p->rid < 0) p->flag |= 0x4;                /* self unmapped */
                if (p->is_rev)  p->flag |= 0x10;               /* self reverse  */
            }
        }
    }

    /* No-pairing branch 0x2 decision (mirrors mem_sam_pe's post-resolve block):
     * when mem_pair_resolve_batch_post didn't take the paired branch, decide the
     * proper-pair bit from the two sides' regions.
     *
     * v0.9.0 factored this out of mem_sam_pe into mem_proper_pair_extra_flag
     * explicitly so "this block and its verbatim twin cannot drift apart"
     * (bwamem_pair.cpp:352-357). This shim was a THIRD twin, so call the one
     * definition rather than keep replicating it.
     *
     * It also settles which region the bit derives from: v0.9.0 restored the
     * top-scoring a[0] as the default and made the emitted a[which] opt-in via
     * --proper-pair-from-emitted. Our previous hand-rolled copy used a[0] and so
     * already matched the new default -- but only by coincidence, and it could
     * not honour the option at all. */
    /* This tail is the no-pairing path (paired==1 returns early above), so the
     * proper-pair bit is always derived here from which[]/h_rid[] -- the former
     * `!paired &&` guard was always true and has been dropped. */
    if ((opt->flag & MEM_F_PE) && !(opt->flag & MEM_F_NOPAIRING)
        && which[0] >= 0 && which[1] >= 0
        && h_rid[0] == h_rid[1] && h_rid[0] >= 0) {
        extra_flag |= mem_proper_pair_extra_flag(opt, bns->l_pac, a, which, pes);
    }

    /* Apply extra_flag (0x1 paired + optional 0x2 proper-pair) to every
     * emitted record, matching mem_reg2sam's behavior in the no_pairing
     * branch. (The paired branch applies extra_flag itself, above.) */
    for (int k = 0; k < 2; ++k) {
        for (int j = 0; j < n_lists[k]; ++j) {
            lists[k][j].flag |= extra_flag;
        }
    }

    /* Build the compact emitted list per side (mirrors mem_reg2sam's `aa`):
     * only regions with emit[k][j] set, with the supplementary flag and lowered
     * MAPQ applied up front. append_bam_record's SA:Z serialization walks this
     * compact array, so SA entries reference only emitted records — never a
     * secondary (folded into XA) or a sub-threshold/drop-ratio region that was
     * filtered out. Entries are shallow copies sharing each region's cigar
     * buffer with lists[k]; those buffers are freed once via lists[k] below, so
     * aa[k] entries are never cigar-freed. */
    mem_aln_t *aa[2] = {nullptr, nullptr};
    int n_aa[2] = {0, 0};
    for (int k = 0; k < 2; ++k) {
        aa[k] = (mem_aln_t *) calloc((size_t)n_lists[k], sizeof(mem_aln_t));
        xassert(aa[k] != NULL, "out of memory: aa[k]");
        int l = 0;               /* emitted-so-far on this side */
        for (int j = 0; j < n_lists[k]; ++j) {
            if (!emit[k][j]) continue;   /* secondary folded into XA, or dropped */
            mem_aln_t p = lists[k][j];   /* shallow copy (shares cigar ptr) */
            if (l > 0) {
                /* 2nd+ emitted region (all secondary<0 here) is supplementary;
                 * lower its mapq to the primary's unless -5/-q, per mem_reg2sam.
                 * 0x10000 (not 0x100) under MEM_F_NO_MULTI matches upstream's
                 * internal marker; append_bam_record remaps it to 0x100 on write.
                 *
                 * The lowering and the repetitive-supp hard cap key off the
                 * ALNREG fields (a[k].a[j]), exactly as mem_reg2sam
                 * (bwamem.cpp:4489-4499). The `!is_alt` guard is load-bearing:
                 * a supplementary on an ALT/decoy contig (is_alt set via the
                 * .alt file) keeps its own MAPQ and is NOT lowered to the
                 * primary's -- omitting it zeroed those, which surfaced only
                 * against a real `.alt` (both in the supp record's own MAPQ and
                 * in every SA:Z tag referencing it). lists[k] is 1:1 with a[k]
                 * on this branch (the n_emit==0 path never reaches l>0). */
                p.flag |= (opt->flag & MEM_F_NO_MULTI) ? 0x10000 : 0x800;
                if (!(opt->flag & MEM_F_KEEP_SUPP_MAPQ) && !a[k].a[j].is_alt
                    && p.mapq > aa[k][0].mapq)
                    p.mapq = aa[k][0].mapq;
                if (opt->supp_rep_hard_cap > 0 && a[k].a[j].secondary < 0
                    && a[k].a[j].chain_n_hits >= opt->supp_rep_hard_cap)
                    p.mapq = 0;
            }
            aa[k][l++] = p;
        }
        n_aa[k] = l;
    }

    /* Emit records. This is the NO-PAIRING branch (the paired branch returns
     * above), so the mate anchor is the other side's *selected* primary
     * which[!k] -- exactly upstream's `&h[!k]`, where h[i] =
     * mem_reg2aln(a[i].a[which[i]]) drives the no-pairing mem_reg2sam calls
     * (bwamem_pair.cpp:1346-1361). which[!k] is region 0 in the common case,
     * but n_pri[!k] when region 0 falls below opt->T (e.g. a read whose best
     * hit is an ALT/decoy) -- using lists[!k][0] there drove RNEXT/PNEXT/TLEN
     * and the 0x20 mate-strand bit off the wrong region. When which[!k] < 0 the
     * other side is unmapped: n_emit[!k]==0 left lists[!k] as the synthetic
     * unmapped at [0], so index 0 is the correct anchor. */
    for (int k = 0; k < 2; ++k) {
        int mate_idx = (which[!k] >= 0 && which[!k] < n_lists[!k]) ? which[!k] : 0;
        mem_aln_t *mate = (n_lists[!k] > 0) ? &lists[!k][mate_idx] : nullptr;
        for (int j = 0; j < n_aa[k]; ++j) {
            append_bam_record(e, origin_idx, opt, bns, pac, &s[k],
                              &aa[k][j], n_aa[k], aa[k], j, mate, k);
        }
    }

    /* Free the cigar buffers (owned by lists[k]), the compact copies, the emit
     * masks, and the XA/HN arrays. aa[k] entries share cigar with lists[k] and
     * must not be cigar-freed. */
    for (int k = 0; k < 2; ++k) {
        if (lists[k]) {
            for (int j = 0; j < n_lists[k]; ++j) free(lists[k][j].cigar);
            free(lists[k]);
        }
        free(aa[k]);
        free(emit[k]);
        if (XA[k]) {
            for (int j = 0; j < (int)a[k].n; ++j) free(XA[k][j]);
            free(XA[k]);
        }
        free(HN[k]);
    }
}

/* Single-end resolve + emit: the SE branch of worker_sam (bwamem.cpp:2880-2887)
 * followed by mem_reg2sam's emit policy (bwamem.cpp:3227-3287), in BAM. `m` is
 * always NULL for a single — append_bam_record then writes mtid/mpos = -1,
 * tlen = 0 and sets no mate/pair flag bits (0x1/0x40/0x80/0x8/0x20). */
static void single_and_emit(ShimEmit *e, size_t origin_idx, uint64_t id,
                            const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                            bseq1_t *s, mem_alnreg_v *a)
{
    mem_mark_primary_se(opt, a->n, a->a, (int64_t)id);
#if V17
    if (opt->flag & MEM_F_PRIMARY5) mem_reorder_primary5(opt->T, a);
#endif
    char **XA = nullptr; int *HN = nullptr;
    if (!(opt->flag & MEM_F_ALL))
        XA = mem_gen_alt(opt, bns, pac, a, s->l_seq, s->seq, &HN, s->meth_orig_seq);

    /* Compact emitted list `aa`, exactly mem_reg2sam's loop. */
    mem_aln_t *aa = (mem_aln_t *) calloc(a->n ? a->n : 1, sizeof(mem_aln_t));
    xassert(aa != NULL, "out of memory: aa");
    int l = 0;
    for (int k = 0; k < (int)a->n; ++k) {
        const mem_alnreg_t *p = &a->a[k];
        if (p->score < opt->T) continue;
        if (p->secondary >= 0 && (p->is_alt || !(opt->flag & MEM_F_ALL))) continue;
        if (p->secondary >= 0 && p->secondary < INT_MAX
            && p->score < a->a[p->secondary].score * opt->drop_ratio) continue;
        mem_aln_t q = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, p, s->meth_orig_seq);
        q.XA = XA ? XA[k] : nullptr;
        q.HN = HN ? HN[k] : -1;
        if (p->secondary >= 0) q.sub = -1;
        if (l && p->secondary < 0) q.flag |= (opt->flag & MEM_F_NO_MULTI) ? 0x10000 : 0x800;
        if (!(opt->flag & MEM_F_KEEP_SUPP_MAPQ) && l && !p->is_alt && q.mapq > aa[0].mapq)
            q.mapq = aa[0].mapq;
        if (opt->supp_rep_hard_cap > 0 && l && p->secondary < 0
            && p->chain_n_hits >= opt->supp_rep_hard_cap)
            q.mapq = 0;
        aa[l++] = q;
    }
    if (l == 0) {
        mem_aln_t t = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, nullptr, s->meth_orig_seq);
        append_bam_record(e, origin_idx, opt, bns, pac, s, &t, 1, &t, 0, nullptr, 0);
        free(t.cigar);
    } else {
        for (int k = 0; k < l; ++k)
            append_bam_record(e, origin_idx, opt, bns, pac, s, &aa[k], l, aa, k, nullptr, 0);
        for (int k = 0; k < l; ++k) free(aa[k].cigar);
    }
    free(aa);
    if (XA) { for (int k = 0; k < (int)a->n; ++k) free(XA[k]); free(XA); }
    free(HN);
}

/* Legacy emit half: per-batch pestat + pair/emit into a ShimAlignOutput, over a
 * caller-supplied scratch so seed and emit can share ONE ~24 MB ShimScratch on
 * the fused align_batch path (Task [16]). CONSUMES `r` (shim_pair_emit frees it
 * on every path) but NOT `sc` (the caller owns and frees it). Returns NULL on
 * OOM or emit failure. */
static ShimAlignOutput *legacy_emit_with_scratch(void *idx_opaque, const mem_opt_t *opts,
                                                 ShimScratch *sc, ShimRegs *r,
                                                 const mem_pestat_t *pestat_in)
{
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(idx_opaque);
    /* D3 (--meth): after the seed→original remap in the kernels, every
     * downstream consumer (insert-size, pairing, mem_reg2aln, output coords)
     * runs in ORIGINAL-reference coordinates — use the original bns, not the
     * converted seed index's. Outside --meth this is the seed index. */
    const bntseq_t *bns = idx->meth_orig_bns ? idx->meth_orig_bns : idx->fmi->idx->bns;

    mem_opt_t opt_pe = *opts;
    opt_pe.n_threads = 1;
    opt_pe.flag |= MEM_F_PE;

    mem_pestat_t pes[4];
    /* Auto-pestat over the PAIR reads only (2*n_pairs), never r->n_seqs: a mixed
     * cohort's singles are not mates and must not be fed to mem_pestat as fake
     * pairs (they would corrupt the insert-size histogram). Mirrors
     * shim_pestat_cohort. On the legacy path n_singles is always 0, so this is
     * byte-identical to the former r->n_seqs there. (Task [7]) */
    if (pestat_in) memcpy(pes, pestat_in, sizeof(pes));
    else           mem_pestat(&opt_pe, bns->l_pac, 2 * (int)r->n_pairs, r->regs, pes);

    ShimAlignOutput *out = alloc_align_output(r->n_pairs);
    if (!out) { shim_regs_free(r); return nullptr; }
    /* Pass `r` straight through — no throwaway ShimSeeds wrap/unwrap (Task [13]).
     * shim_pair_emit consumes (frees) `r` on every return path. */
    ShimIdBases ids = { 0, 0 };
    int rc = shim_pair_emit(idx_opaque, opts, sc, r, pes, ids, legacy_out_sink, out);
    if (rc != 0) { shim_align_out_free(out); return nullptr; }
    memcpy(out->pes, pes, sizeof(pes));
    return out;
}

ShimAlignOutput *shim_extend_batch(void *idx_opaque, const mem_opt_t *opts, ShimSeeds *s,
                                   const mem_pestat_t *pestat_in)
{
    if (!s) return nullptr;
    ShimRegs *r = s->regs; s->regs = nullptr; free(s);
    /* Standalone extend from pre-seeded regs: the seed phase's scratch is long
     * gone, so allocate one for the emit phase. */
    ShimScratch *sc = shim_scratch_new();
    if (!sc) { shim_regs_free(r); return nullptr; }
    ShimAlignOutput *out = legacy_emit_with_scratch(idx_opaque, opts, sc, r, pestat_in);
    shim_scratch_free(sc);
    return out;
}

/* ------------------ Convenience: align_batch ------------------ */

ShimAlignOutput *shim_align_batch(void *idx_opaque, const mem_opt_t *opts,
                                  const ShimReadPair *pairs, size_t n_pairs,
                                  const mem_pestat_t *pestat_in)
{
    /* Fused seed+emit shares ONE ShimScratch instead of seed_batch's and
     * extend_batch's two separate ~24 MB alloc+free cycles (Task [16]). */
    ShimReadBatch b = { pairs, n_pairs, nullptr, 0 };
    ShimScratch *sc = shim_scratch_new();
    if (!sc) return nullptr;
    ShimRegs *r = shim_seed_extend(idx_opaque, opts, sc, &b);
    if (!r) { shim_scratch_free(sc); return nullptr; }
    ShimAlignOutput *out = legacy_emit_with_scratch(idx_opaque, opts, sc, r, pestat_in);
    shim_scratch_free(sc);
    return out;
}

/* ------------------ estimate_pestat ------------------ */

int shim_estimate_pestat(void *idx_opaque, const mem_opt_t *opts,
                         const ShimReadPair *pairs, size_t n_pairs,
                         mem_pestat_t *pestat_out)
{
    ShimSeeds *s = shim_seed_batch(idx_opaque, opts, pairs, n_pairs);
    if (!s) return -1;
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(idx_opaque);
    const bntseq_t *bns = idx->meth_orig_bns ? idx->meth_orig_bns : idx->fmi->idx->bns;
    mem_opt_t opt_pe = *opts; opt_pe.n_threads = 1; opt_pe.flag |= MEM_F_PE;
    mem_pestat(&opt_pe, bns->l_pac, s->regs->n_seqs, s->regs->regs, pestat_out);
    shim_seeds_free(s);
    return 0;
}

/* ------------------ Output accessors ------------------ */

size_t shim_align_out_n_recs(ShimAlignOutput *out)         { return out ? out->n_recs : 0; }
size_t shim_align_out_pair_idx(ShimAlignOutput *out, size_t i)   { return out->pair_idx[i]; }
const uint8_t *shim_align_out_rec_ptr(ShimAlignOutput *out, size_t i) {
    return out->buf + out->rec_off[i];
}
size_t shim_align_out_rec_len(ShimAlignOutput *out, size_t i)    { return out->rec_len[i]; }
void shim_align_get_pestat(ShimAlignOutput *out, mem_pestat_t *dst) {
    memcpy(dst, out->pes, sizeof(out->pes));
}
void shim_align_out_free(ShimAlignOutput *out) {
    if (!out) return;
    free(out->buf); free(out->rec_off); free(out->rec_len); free(out->pair_idx);
    free(out);
}

/* ==================================================================
 * Resident-cohort aligner
 *
 * Keeps a -K cohort's decoded reads (bseq1_t[]) and their alignment
 * regions (mem_alnreg_v[]) RESIDENT across the seed_extend -> pestat ->
 * pair_emit phases, instead of the self-contained-per-sub-batch ShimRegs
 * that shim_seed_extend/shim_pair_emit allocate and free on every call.
 * Each sub-chunk owns one reserved SEGMENT of the resident cohort; no
 * sub-chunk allocates or frees the cohort's header arrays, and each segment
 * releases its own reads and regions once it is emitted. A segment is a separate
 * allocation, so a later reserve() -- which only grows the table of segment
 * pointers -- never moves it.
 *
 * The per-segment calls run the same helpers as the legacy path
 * (decode_read, seed_extend_reads, pestat_over_headers,
 * pair_emit_pairs_chunked, emit_singles), so the two paths cannot diverge.
 * The safe wrapper is bwa-mem3-rs's `ResidentCohort`.
 * ================================================================== */

/* One reserved sub-chunk of a resident region: its reads, their alnreg
 * headers, and its lifecycle state. Allocated separately (one per reserve) so
 * its address is stable for the cohort's whole lifetime -- a later reserve only
 * grows the region's table of segment POINTERS, never moves a segment. Every
 * per-segment call takes the segment directly and never consults the table,
 * which is what lets those calls run concurrently with a reserve.
 *
 * Only the one caller holding a segment's range touches it (the bwa-mem3-rs
 * wrapper enforces that with a move-only range token); the whole-cohort call
 * (pestat) reads every segment and must exclude per-segment calls. */
struct ShimResidentSegment {
    bseq1_t      *seqs;       /* n reads; pairs R1/R2 interleaved */
    mem_alnreg_v *regs;       /* parallel per-read alnreg headers */
    uint8_t      *arena;      /* every read's strings, when batch-written; else NULL */
    size_t        n;          /* reads */
    size_t        base;       /* inclusive-start read offset within its region */
    int           is_pairs;
    int           meth_mode;
    int           extended;   /* seed_extend has run; the reads are 2-bit encoded */
    int           emitted;    /* pair_emit has run; mate rescue has grown its regs */
};

/* One resident region (pairs or singles) for one bwa-mem3 -p group: a growable
 * table of segment pointers, in reserve order. */
struct ShimResidentRegion {
    ShimResidentSegment **segs;
    size_t                n_segments;
    size_t                cap_segments;
    size_t                len;        /* total reserved reads */
};

struct ShimResidentCohort {
    ShimResidentRegion pairs;    /* 2*n_pairs reads, R1/R2 interleaved */
    ShimResidentRegion singles;
    int meth_mode;
};

/* Release a segment's reads: its arena (or per-read strings) and any regions
 * left, zeroing the headers. Runs at the end of the segment's emit -- safe
 * because `emitted` makes every later call on the segment a lifecycle error
 * and pestat refuses a cohort with an emitted segment, so the cohort then holds
 * only the zeroed headers -- and as the first step of resident_segment_free,
 * where it is a no-op for an already-released segment. */
static void resident_segment_release_reads(ShimResidentSegment *sg) {
    for (size_t i = 0; i < sg->n; ++i) {
        /* Same per-read teardown as free_seqs + shim_regs_free's regs loop:
         * free() is NULL-safe on the fields a partial write left unbuilt. A
         * batch-written segment's strings live in its arena instead. */
        if (!sg->arena) {
            free(sg->seqs[i].name); free(sg->seqs[i].seq); free(sg->seqs[i].qual);
            free(sg->seqs[i].meth_orig_seq);
        }
        free(sg->seqs[i].sam);
        free(sg->regs[i].a);
    }
    free(sg->arena);
    sg->arena = nullptr;
    if (sg->n) {
        memset(sg->seqs, 0, sg->n * sizeof(bseq1_t));
        memset(sg->regs, 0, sg->n * sizeof(mem_alnreg_v));
    }
}

/* Whether `sg` still holds any read string or alignment region -- false once
 * it has been emitted (and released). A test hook. */
int shim_resident_segment_holds_reads(const ShimResidentSegment *sg) {
    if (!sg) return 0;
    if (sg->arena) return 1;
    for (size_t i = 0; i < sg->n; ++i)
        if (sg->seqs[i].name || sg->seqs[i].seq || sg->regs[i].a) return 1;
    return 0;
}

static void resident_segment_free(ShimResidentSegment *sg) {
    resident_segment_release_reads(sg);
    free(sg->seqs); free(sg->regs); free(sg);
}

static void resident_region_free(ShimResidentRegion *rg) {
    for (size_t s = 0; s < rg->n_segments; ++s) resident_segment_free(rg->segs[s]);
    free(rg->segs);
    memset(rg, 0, sizeof(*rg));
}

/* Append a new segment of `n_reads` reads and return it (NULL on allocation
 * failure, leaving the region unchanged). Single-writer: must not run
 * concurrently with another reserve or with pestat on this cohort, but MAY run
 * concurrently with per-segment calls on already-returned segments. */
static ShimResidentSegment *resident_region_reserve(ShimResidentRegion *rg, size_t n_reads,
                                                    int is_pairs, int meth_mode) {
    if (rg->n_segments == rg->cap_segments) {
        size_t newcap = rg->cap_segments ? rg->cap_segments * 2 : 8;
        ShimResidentSegment **ns =
            (ShimResidentSegment **) realloc(rg->segs, newcap * sizeof(*ns));
        if (!ns) return nullptr;
        rg->segs = ns;
        rg->cap_segments = newcap;
    }
    ShimResidentSegment *sg = (ShimResidentSegment *) calloc(1, sizeof(ShimResidentSegment));
    if (!sg) return nullptr;
    /* calloc(0,...) may legally return NULL; guard n_reads==0 so that is not
     * mistaken for OOM (matches copy_pairs_to_seqs's empty-batch handling). */
    sg->seqs = (bseq1_t *) calloc(n_reads ? n_reads : 1, sizeof(bseq1_t));
    sg->regs = (mem_alnreg_v *) calloc(n_reads ? n_reads : 1, sizeof(mem_alnreg_v));
    if (!sg->seqs || !sg->regs) { free(sg->seqs); free(sg->regs); free(sg); return nullptr; }
    sg->n = n_reads;
    sg->base = rg->len;
    sg->is_pairs = is_pairs;
    sg->meth_mode = meth_mode;
    rg->segs[rg->n_segments++] = sg;
    rg->len += n_reads;
    return sg;
}

ShimResidentCohort *shim_resident_cohort_new(int meth_mode) {
    ShimResidentCohort *c = (ShimResidentCohort *) calloc(1, sizeof(ShimResidentCohort));
    if (!c) return nullptr;
    c->meth_mode = meth_mode;   /* pairs/singles are zeroed by calloc */
    return c;
}

void shim_resident_cohort_free(ShimResidentCohort *c) {
    if (!c) return;
    resident_region_free(&c->pairs);
    resident_region_free(&c->singles);
    free(c);
}


/* Heap bytes a reserved read costs before anything is written into it: its
 * bseq1_t and mem_alnreg_v headers. */
size_t shim_resident_read_overhead(void) {
    return sizeof(bseq1_t) + sizeof(mem_alnreg_v);
}

ShimResidentSegment *shim_resident_reserve_pairs(ShimResidentCohort *c, size_t n_reads,
                                                 size_t *first_out) {
    if (!c || !first_out || n_reads % 2 != 0) return nullptr;
    ShimResidentSegment *sg = resident_region_reserve(&c->pairs, n_reads, 1, c->meth_mode);
    if (sg) *first_out = sg->base;
    return sg;
}
ShimResidentSegment *shim_resident_reserve_singles(ShimResidentCohort *c, size_t n_reads,
                                                   size_t *first_out) {
    if (!c || !first_out) return nullptr;
    ShimResidentSegment *sg = resident_region_reserve(&c->singles, n_reads, 0, c->meth_mode);
    if (sg) *first_out = sg->base;
    return sg;
}

/* Write pair `i` of a pair segment (reads 2i, 2i+1), adding the bytes it
 * copies to *added. Returns 0, -1 on a null arg / bad index / wrong region /
 * OOM, or -3 when the segment was already seed-extended or the slot already
 * written (rewriting would leak the old buffers, and an ASCII read under 2-bit
 * regs would feed the kernels garbage). */
int shim_resident_write_pair(ShimResidentSegment *sg, size_t i, const ShimReadPair *pair,
                             size_t *added) {
    if (!sg || !pair || !added || !sg->is_pairs || i >= sg->n / 2) return -1;
    bseq1_t *r1 = &sg->seqs[2 * i], *r2 = &sg->seqs[2 * i + 1];
    if (sg->extended || sg->arena || r1->name || r2->name) return -3;
    if (decode_read(r1, pair->r1_name, pair->r1_name_len, pair->r1_seq, pair->r1_seq_len,
                    pair->r1_qual, sg->meth_mode, 0, added) != 0) return -1;
    if (decode_read(r2, pair->r2_name, pair->r2_name_len, pair->r2_seq, pair->r2_seq_len,
                    pair->r2_qual, sg->meth_mode, 1, added) != 0) return -1;
    return 0;
}
int shim_resident_write_single(ShimResidentSegment *sg, size_t i, const ShimSingleRead *single,
                               size_t *added) {
    if (!sg || !single || !added || sg->is_pairs || i >= sg->n) return -1;
    bseq1_t *s = &sg->seqs[i];
    if (sg->extended || sg->arena || s->name) return -3;
    return decode_read(s, single->name, single->name_len, single->seq, single->seq_len,
                       single->qual, sg->meth_mode, 0, added);
}

/* Over-allocation past a segment arena's last string, so a vectorized read
 * that runs past the final read's bytes stays inside the allocation. */
#define RESIDENT_ARENA_SLACK 64

/* 0 when `sg` can take a batch write (nothing written, not extended), else -3. */
static int resident_batch_writable(const ShimResidentSegment *sg) {
    if (sg->extended || sg->arena) return -3;
    for (size_t i = 0; i < sg->n; ++i)
        if (sg->seqs[i].name) return -3;
    return 0;
}

/* Write every pair of an unwritten pair segment in one call, placing all the
 * reads' strings in ONE allocation (the segment arena), as the CLI's reader
 * does per chunk, instead of three per read. `n` must be the segment's pair
 * count. Adds the bytes copied to *added -- the same total the per-slot
 * writes would add. Returns 0, -1 on a null arg / wrong region / count
 * mismatch / OOM, or -3 when any slot was already written or the segment was
 * extended. */
int shim_resident_write_pairs(ShimResidentSegment *sg, const ShimReadPair *pairs, size_t n,
                              size_t *added) {
    if (!sg || !added || !sg->is_pairs || (n > 0 && !pairs) || n != sg->n / 2) return -1;
    int rc = resident_batch_writable(sg);
    if (rc != 0) return rc;
    size_t total = 0;
    for (size_t i = 0; i < n; ++i) {
        const ShimReadPair *p = &pairs[i];
        total += read_arena_bytes(p->r1_name_len, p->r1_seq_len, p->r1_qual, sg->meth_mode)
               + read_arena_bytes(p->r2_name_len, p->r2_seq_len, p->r2_qual, sg->meth_mode);
    }
    sg->arena = (uint8_t *) malloc(total + RESIDENT_ARENA_SLACK);
    if (!sg->arena) return -1;
    ArenaCursor cur = { sg->arena };
    for (size_t i = 0; i < n; ++i) {
        const ShimReadPair *p = &pairs[i];
        decode_read_arena(&sg->seqs[2 * i], &cur, p->r1_name, p->r1_name_len, p->r1_seq,
                          p->r1_seq_len, p->r1_qual, sg->meth_mode, 0);
        decode_read_arena(&sg->seqs[2 * i + 1], &cur, p->r2_name, p->r2_name_len, p->r2_seq,
                          p->r2_seq_len, p->r2_qual, sg->meth_mode, 1);
    }
    *added += total;
    return 0;
}

/* shim_resident_write_pairs for a single segment: `n` must be its read count. */
int shim_resident_write_singles(ShimResidentSegment *sg, const ShimSingleRead *reads, size_t n,
                                size_t *added) {
    if (!sg || !added || sg->is_pairs || (n > 0 && !reads) || n != sg->n) return -1;
    int rc = resident_batch_writable(sg);
    if (rc != 0) return rc;
    size_t total = 0;
    for (size_t i = 0; i < n; ++i)
        total += read_arena_bytes(reads[i].name_len, reads[i].seq_len, reads[i].qual,
                                  sg->meth_mode);
    sg->arena = (uint8_t *) malloc(total + RESIDENT_ARENA_SLACK);
    if (!sg->arena) return -1;
    ArenaCursor cur = { sg->arena };
    for (size_t i = 0; i < n; ++i)
        decode_read_arena(&sg->seqs[i], &cur, reads[i].name, reads[i].name_len, reads[i].seq,
                          reads[i].seq_len, reads[i].qual, sg->meth_mode, 0);
    *added += total;
    return 0;
}

/* Seed + SE-extend a whole segment (pairs or singles). Returns 0, -1 on a null
 * arg, or -3 when it was already extended or any slot is unwritten (a NULL seq
 * would crash the seeding kernel). */
int shim_resident_seed_extend(void *idx_opaque, const mem_opt_t *opts, ShimScratch *sc,
                              ShimResidentSegment *sg) {
    if (!idx_opaque || !opts || !sc || !sg) return -1;
    if (sg->extended) return -3;
    for (size_t i = 0; i < sg->n; ++i)
        if (!sg->seqs[i].seq) return -3;
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(idx_opaque);
    /* The per-group opt copy shim_seed_extend makes (MEM_F_PE on for pairs, off
     * for singles). */
    mem_opt_t opt = *opts; opt.n_threads = 1;
    if (sg->is_pairs) opt.flag |= MEM_F_PE; else opt.flag &= ~MEM_F_PE;
    seed_extend_reads(sc, idx, &opt, sg->seqs, sg->regs, (int) sg->n);
    sg->extended = 1;
    return 0;
}

/* mem_pestat over the WHOLE resident pair region (every read there is a PE
 * read). Gathers the 24-byte alnreg headers across segments into one contiguous
 * array, like shim_pestat_cohort's gather, but off the resident segments, so a
 * caller no longer has to collect every sub-batch's regs to call it. A
 * whole-cohort call: must not run concurrently with any per-segment call or
 * reserve. Returns 0, -1 on a null arg / OOM, or -3 when some pair segment is
 * not seed-extended yet (the model would silently omit its reads) or has
 * already been emitted (mate rescue has appended to its regs, so the model
 * would no longer be the one the CLI computes before pairing). */
int shim_resident_pestat_cohort(void *idx_opaque, const mem_opt_t *opts,
                                const ShimResidentCohort *c, mem_pestat_t *out) {
    if (!idx_opaque || !opts || !c || !out) return -1;
    for (size_t s = 0; s < c->pairs.n_segments; ++s)
        if (!c->pairs.segs[s]->extended || c->pairs.segs[s]->emitted) return -3;
    size_t total = c->pairs.len;
    mem_alnreg_v *cat = total ? (mem_alnreg_v *) malloc(total * sizeof(mem_alnreg_v)) : nullptr;
    if (total && !cat) return -1;
    size_t at = 0;
    for (size_t s = 0; s < c->pairs.n_segments; ++s) {
        const ShimResidentSegment *sg = c->pairs.segs[s];
        if (sg->n) memcpy(cat + at, sg->regs, sg->n * sizeof(mem_alnreg_v));
        at += sg->n;
    }
    pestat_over_headers(static_cast<BwaShimIndex *>(idx_opaque), opts, cat, total, out);
    free(cat);
    return 0;
}

/* Pair + mate-rescue + emit a pair segment, or SE-emit a single segment, then
 * release its reads and regions (its headers stay until the cohort is freed).
 * `ids.first_pair_id` / `first_single_id` is the GLOBAL
 * read ordinal of the segment's first pair/single; `origin_base` is added to
 * the local index for the sink's origin_idx. Exactly one of `sink` /
 * `field_sink` is set. Returns 0, -1 on a null arg, -2 when a pair segment has
 * pairs and no pestat, or -3 when it was not seed-extended or was already
 * emitted (mate rescue appended to its regs on the first emission, so a second
 * would not reproduce it). */
static int resident_pair_emit(void *idx_opaque, const mem_opt_t *opts, ShimScratch *sc,
                              ShimResidentSegment *sg, const mem_pestat_t *pestat,
                              ShimIdBases ids, size_t origin_base, ShimRecordSinkFn sink,
                              ShimFieldSinkFn field_sink, void *ctx) {
    if (!idx_opaque || !opts || !sc || !sg || (!sink && !field_sink)) return -1;
    if (sg->is_pairs && sg->n > 0 && !pestat) return -2;
    if (!sg->extended || sg->emitted) return -3;
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(idx_opaque);
    const bntseq_t *bns = idx->meth_orig_bns ? idx->meth_orig_bns : idx->fmi->idx->bns;
    const uint8_t  *pac = idx->meth_orig_pac ? idx->meth_orig_pac : idx->fmi->idx->pac;
    sg->emitted = 1;
    if (sg->is_pairs) {
        mem_opt_t opt_pe = *opts; opt_pe.n_threads = 1; opt_pe.flag |= MEM_F_PE;
        ShimEmit e = { sc, sink, ctx, 0u /* BWA_ORIGIN_PAIR */, field_sink };
        pair_emit_pairs_chunked(&e, &opt_pe, bns, pac, pestat, ids.first_pair_id,
                                sg->seqs, sg->regs, sg->n / 2, origin_base);
    } else {
        mem_opt_t opt_se = *opts; opt_se.n_threads = 1; opt_se.flag &= ~MEM_F_PE;
        ShimEmit e = { sc, sink, ctx, 1u /* BWA_ORIGIN_SINGLE */, field_sink };
        emit_singles(&e, &opt_se, bns, pac, ids.first_single_id, sg->seqs, sg->regs, sg->n,
                     origin_base);
    }
    resident_segment_release_reads(sg);
    return 0;
}

int shim_resident_pair_emit(void *idx_opaque, const mem_opt_t *opts, ShimScratch *sc,
                            ShimResidentSegment *sg, const mem_pestat_t *pestat,
                            ShimIdBases ids, size_t origin_base,
                            ShimRecordSinkFn sink, void *ctx) {
    if (!sink) return -1;
    return resident_pair_emit(idx_opaque, opts, sc, sg, pestat, ids, origin_base,
                              sink, nullptr, ctx);
}

/* shim_resident_pair_emit, reporting each record's structured fields to
 * `field_sink` instead of its packed BAM body. Same resolution, same records,
 * same order: both variants share resident_pair_emit and differ only in
 * append_bam_record's last step. */
int shim_resident_pair_emit_fields(void *idx_opaque, const mem_opt_t *opts, ShimScratch *sc,
                                   ShimResidentSegment *sg, const mem_pestat_t *pestat,
                                   ShimIdBases ids, size_t origin_base,
                                   ShimFieldSinkFn field_sink, void *ctx) {
    if (!field_sink) return -1;
    return resident_pair_emit(idx_opaque, opts, sc, sg, pestat, ids, origin_base,
                              nullptr, field_sink, ctx);
}

} /* extern "C" */
