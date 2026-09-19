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

/* sort_classify has external linkage in bwamem.cpp but no header declaration
 * (the CLI's worker_sam calls it from the same TU). The batched mate-rescue
 * path in shim_pair_emit needs it, so forward-declare it here. */
extern int64_t sort_classify(mem_cache *mmc, int64_t pcnt, int tid);

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
/* Allocate the per-thread scratch a fused seed+extend needs. `nthreads` is
 * always 1 here (one ShimScratch per caller thread), kept as a parameter so
 * the loops stay textually close to upstream's worker_alloc. Chain/seed
 * windows are BATCH_SIZE-sized because seeding and extension now run
 * chunk-by-chunk inside one call (upstream's v0.9.0 fusion, bwamem.cpp:2782):
 * chains never outlive the chunk that produced them. `regs` are NOT here --
 * they belong to the ShimRegs that outlives the call. */
static void worker_alloc(worker_t &w, int32_t nthreads)
{
    assert(nthreads > 0);
    if (nthreads < 1) nthreads = 1;
    w.nthreads = nthreads;
    w.regs = NULL;
    /* Vestigial worker_t members that nothing allocates (every real
     * `auxSeedBuf` upstream is a local in test_and_merge, bwamem.cpp:894); zero
     * them so worker_free and any future upstream code cannot act on an
     * indeterminate pointer. */
    w.auxSeedBuf = NULL; w.auxSeedBufSize = 0;
    /* BATCH_SIZE-sized, chunk-local windows. Seeding and extension are fused
     * per BATCH_SIZE chunk in seed_extend_range, so a chain is dead the moment
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
    for(int l=0; l<nthreads; l++)
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

    for(int l=0; l<nthreads; l++) {
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
    for (int l=0; l<nthreads; l++)
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
    }
}

static void worker_free(worker_t &w, int32_t nthreads)
{
    assert(nthreads > 0);
    // Catch mismatched alloc/free pairs before they drive out-of-bounds frees.
    assert(w.nthreads == nthreads);

    free(w.chain_scratch);
    /* w.regs is NOT freed here: during seed_extend_range it aliases the
     * caller-owned ShimRegs::regs, which shim_regs_free owns. worker_alloc
     * leaves it NULL and nothing in the scratch ever allocates it. */
    free(w.seed_scratch);

    for(int l=0; l<nthreads; l++) {
        _mm_free(w.mmc.seqBufLeftRef[l*CACHE_LINE]);
        _mm_free(w.mmc.seqBufRightRef[l*CACHE_LINE]);
        _mm_free(w.mmc.seqBufLeftQer[l*CACHE_LINE]);
        _mm_free(w.mmc.seqBufRightQer[l*CACHE_LINE]);
    }

    for(int l=0; l<nthreads; l++) {
        free(w.mmc.seqPairArrayAux[l]);
        free(w.mmc.seqPairArrayLeft128[l]);
        free(w.mmc.seqPairArrayRight128[l]);
    }

    // NULL-safe: SMEM buffers are now allocated lazily on first batch;
    // workers that never ran a batch leave them as NULL. _mm_free / free
    // are both well-defined on NULL.
    for(int l=0; l<nthreads; l++) {
        _mm_free(w.mmc.matchArray[l]);
        free(w.mmc.min_intv_ar[l]);
        free(w.mmc.query_pos_ar[l]);
        free(w.mmc.enc_qdb[l]);
        free(w.mmc.rid[l]);
        _mm_free(w.mmc.lim[l]);

        _mm_free(w.mmc.lockstep_prev[l]);
        _mm_free(w.mmc.lockstep_match_buf[l]);
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
    uint8_t *rec_buf; size_t rec_cap;   /* one packed record under construction */
    uint8_t *aux_buf; size_t aux_cap;   /* one aux block under construction */
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
struct ShimIdBases { uint64_t first_single_id; uint64_t first_pair_id; };

/* Where a record goes once built. Replaces the ShimAlignOutput* that
 * append_bam_record used to write into; the legacy output type is now just
 * one particular sink (legacy_out_sink below). */
struct ShimEmit {
    ShimScratch      *sc;
    ShimRecordSinkFn  sink;
    void             *ctx;
    uint32_t          origin_kind;
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

/* Index handle: owns the loaded FMI_search plus the one-time-unpacked
 * 2*l_pac reference string used by mem_chain2aln_across_reads_V2. The
 * ref_string is built once at load time so every seed_batch borrows it
 * instead of rebuilding (~6GB alloc + full pass over packed ref per call
 * on hs38). Forward declaration of build_ref_string below. */
static uint8_t *build_ref_string(const bntseq_t *bns, const uint8_t *pac);

struct BwaShimIndex {
    FMI_search *fmi;
    uint8_t    *ref_string;
    /* D3 (--meth) dual-coordinate handles. In meth mode `fmi` is the
     * f/r-doubled CONVERTED seed index (`<ref>.meth.*`) used only for candidate
     * generation, while chaining/extension/output run in ORIGINAL coordinates
     * loaded here from the un-converted `<ref>.*` prefix. All NULL for a normal
     * (non-meth) index. `meth_orig_ref_string` is the original 2*l_pac unpacked
     * reference (analogue of `ref_string` for the original ref). */
    bntseq_t *meth_orig_bns;
    uint8_t  *meth_orig_pac;
    uint8_t  *meth_orig_ref_string;
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
    idx->ref_string = build_ref_string(fmi->idx->bns, fmi->idx->pac);
    if (!idx->ref_string) {
        delete idx->fmi;
        free(idx);
        return nullptr;
    }
    return static_cast<void *>(idx);
}

void *shim_align_idx_load(const char *prefix) {
    return shim_align_idx_load_threads(prefix, /*n_threads=*/1);
}

/* D3 (--meth): load a dual index. `seed_prefix` is the converted seed index
 * (`<ref>.meth`), `orig_prefix` the un-converted original reference (`<ref>`).
 * Both the seed FM-index and the original bns/pac/ref_string stay resident. */
void *shim_align_idx_load_meth(const char *seed_prefix, const char *orig_prefix) {
    void *opaque = shim_align_idx_load(seed_prefix);
    if (!opaque) return nullptr;
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(opaque);
    if (shim_meth_orig_ref_load(orig_prefix, &idx->meth_orig_bns, &idx->meth_orig_pac) != 0) {
        shim_align_idx_free(opaque);
        return nullptr;
    }
    idx->meth_orig_ref_string = build_ref_string(idx->meth_orig_bns, idx->meth_orig_pac);
    if (!idx->meth_orig_ref_string) {
        shim_align_idx_free(opaque);
        return nullptr;
    }
    return opaque;
}

void shim_align_idx_free(void *opaque) {
    if (!opaque) return;
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(opaque);
    if (idx->ref_string) _mm_free(idx->ref_string);
    if (idx->meth_orig_ref_string) _mm_free(idx->meth_orig_ref_string);
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
        const char    *nm[2] = {p->r1_name,   p->r2_name};
        size_t         nl[2] = {p->r1_name_len, p->r2_name_len};
        const uint8_t *sq[2] = {p->r1_seq, p->r2_seq};
        size_t         sl[2] = {p->r1_seq_len, p->r2_seq_len};
        const uint8_t *ql[2] = {p->r1_qual, p->r2_qual};
        for (int k = 0; k < 2; ++k) {
            bseq1_t *s = &seqs[2*i + k];
            s->l_seq = (int)sl[k];
            /* NULL-check every per-read allocation: an OOM here must unwind the
             * partially-built array and return NULL (the caller's `!pairs_only`
             * path), not deref NULL. free_seqs() is NULL-safe on the entries not
             * yet built (calloc leaves them NULL). */
            s->name = (char *) malloc(nl[k] + 1);
            if (!s->name) { free_seqs(seqs, nseqs); return nullptr; }
            memcpy(s->name, nm[k], nl[k]);
            s->name[nl[k]] = '\0';
            s->seq = (char *) malloc(sl[k] + 1);
            if (!s->seq) { free_seqs(seqs, nseqs); return nullptr; }
            memcpy(s->seq, sq[k], sl[k]);
            s->seq[sl[k]] = '\0';
            if (ql[k]) {
                s->qual = (char *) malloc(sl[k] + 1);
                if (!s->qual) { free_seqs(seqs, nseqs); return nullptr; }
                memcpy(s->qual, ql[k], sl[k]);
                s->qual[sl[k]] = '\0';
            } else {
                s->qual = nullptr;
            }
            *heap_bytes += (nl[k] + 1) + (sl[k] + 1) + (ql[k] ? sl[k] + 1 : 0);
            s->sam = nullptr;
            /* D3 (--meth): retain the ORIGINAL (unconverted) read bases before
             * projecting seq to match the converted seed index, then apply the
             * per-strand bisulfite projection in place. R1 (k==0) is the OT/CT
             * read (C→T), R2 (k==1) is the OB/GA read (G→A) — matching upstream
             * fastmap.cpp's YC assignment. mem_kernel1_core 2-bit-encodes the
             * projected seq; mem_reg2aln / meth_build_xm use meth_orig_seq for
             * CIGAR/NM/MD and the XM call string. Same orientation as seq. */
            if (meth_mode) {
                s->meth_orig_seq = strdup(s->seq);
                if (!s->meth_orig_seq) { free_seqs(seqs, nseqs); return nullptr; }
                *heap_bytes += (size_t)s->l_seq + 1;
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
                s->meth_base_ot = (k == 0) ? 1 : 0;
                char from = (k == 0) ? 'C' : 'G';
                char to   = (k == 0) ? 'T' : 'A';
                for (int j = 0; j < s->l_seq; ++j) {
                    if (s->seq[j] == from || s->seq[j] == (char)(from + 32))
                        s->seq[j] = to;
                }
            }
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
        bseq1_t *s = &seqs[base + i];
        s->l_seq = (int)p->seq_len;
        s->name = (char *) malloc(p->name_len + 1);
        if (!s->name) return -1;
        memcpy(s->name, p->name, p->name_len); s->name[p->name_len] = '\0';
        s->seq = (char *) malloc(p->seq_len + 1);
        if (!s->seq) return -1;
        memcpy(s->seq, p->seq, p->seq_len); s->seq[p->seq_len] = '\0';
        if (p->qual) {
            s->qual = (char *) malloc(p->seq_len + 1);
            if (!s->qual) return -1;
            memcpy(s->qual, p->qual, p->seq_len); s->qual[p->seq_len] = '\0';
        } else {
            s->qual = nullptr;
        }
        s->sam = nullptr;
        *heap_bytes += (p->name_len + 1) + (p->seq_len + 1) + (p->qual ? p->seq_len + 1 : 0);
        if (meth_mode) {
            s->meth_orig_seq = strdup(s->seq);
            if (!s->meth_orig_seq) return -1;
            s->meth_base_ot = 1;
            for (int j = 0; j < s->l_seq; ++j)
                if (s->seq[j] == 'C' || s->seq[j] == 'c') s->seq[j] = 'T';
            *heap_bytes += (size_t)s->l_seq + 1;
        }
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

/* Unpack 2-bit packed reference into a 1-byte-per-base array with the
 * reverse-complement appended (total = 2 * l_pac). Required by
 * mem_chain2aln_across_reads_V2. */
static uint8_t *build_ref_string(const bntseq_t *bns, const uint8_t *pac) {
    int64_t ref_len = bns->l_pac * 2;
    uint8_t *s = (uint8_t *) _mm_malloc(ref_len, 64);
    /* ~2*l_pac bytes (~6 GB on hs38): the largest single alloc in the crate and
     * the one most likely to fail under memory pressure. Return NULL (both
     * callers check) rather than deref it below. */
    if (!s) return nullptr;
    for (int64_t i = 0; i < bns->l_pac; ++i) {
        uint8_t b = (pac[i >> 2] >> ((~i & 3) << 1)) & 3;
        s[i] = b;
        s[ref_len - 1 - i] = (uint8_t)(3 - b);
    }
    return s;
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

/* Emit an SA:Z tag for primary alignments that have supplementary hits. */
static void emit_sa_tag(uint8_t **buf, size_t *len, size_t *cap,
                        const bntseq_t *bns,
                        const mem_aln_t *list, int n, int which) {
    // Format matches mem_aln2sam:
    //   chrom,pos+1,[+-],CIGAR,mapq,NM;chrom,...
    uint8_t header[3] = {'S', 'A', 'Z'};
    buf_append(buf, len, cap, header, 3);
    char tmp[4096];
    for (int i = 0; i < n; ++i) {
        if (i == which || (list[i].flag & 0x100)) continue;
        const mem_aln_t *r = &list[i];
        int m = snprintf(tmp, sizeof(tmp), "%s,%lld,%c,",
                         bns->anns[r->rid].name, (long long)(r->pos + 1),
                         r->is_rev ? '-' : '+');
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

/* Emit MC:Z aux (mate CIGAR) from a mate mem_aln_t.
 *
 * `which` is the EMITTING record's index, not the mate's: upstream builds
 * MC:Z through the same `add_cigar` helper and hands it the record's own
 * `which` (bwamem.cpp:2583), so on a supplementary record the mate's clips
 * are rendered hard even though the mate itself is a primary. Mirrored here
 * deliberately -- it looks like an upstream quirk, and reproducing it is the
 * point. `m->is_alt` (not the record's) gates it, matching `add_cigar`'s use
 * of the mem_aln_t it was passed. */
static void emit_mc_tag(uint8_t **buf, size_t *len, size_t *cap,
                        const mem_opt_t *opt, const mem_aln_t *m, int which) {
    uint8_t header[3] = {'M', 'C', 'Z'};
    buf_append(buf, len, cap, header, 3);
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

/* Append one packed BAM record to `out`'s buffer. `opt`/`pac`/`is_r2` are used
 * only for D3 (--meth) tag emission (is_r2: 0 = R1/OT read, 1 = R2/OB read). */
static void append_bam_record(ShimEmit *e, size_t origin_idx,
                              const mem_opt_t *opt, const bntseq_t *bns,
                              const uint8_t *pac, const bseq1_t *s,
                              const mem_aln_t *p, int n_list,
                              const mem_aln_t *list, int which,
                              const mem_aln_t *m, int is_r2)
{
    int l_read_name = (int) strlen(s->name) + 1;
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
     * so output stays byte-identical. */
    xassert(l_read_name <= 255, "read name too long for BAM (>= 255 bytes)");
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
     * 0x4/0x8 are deliberately NOT recomputed here. pair_and_emit already set
     * them from the ORIGINAL rids, which matches upstream's ordering: it
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

    /* Build aux first so we know its size. Reuse the scratch's aux buffer: it
     * is grown in place by buf_append and written back at the end of the record
     * so the next record reuses the larger allocation. */
    uint8_t *aux = e->sc->aux_buf;
    size_t aux_len = 0, aux_cap = e->sc->aux_cap;
    if (p->n_cigar) {
        aux_put_i(&aux, &aux_len, &aux_cap, "NM", p->NM);
        /* MD string is stored right after the CIGAR array in p->cigar. */
        const char *md = (const char *)(p->cigar + p->n_cigar);
        aux_put_Z(&aux, &aux_len, &aux_cap, "MD", md);
    }
    if (m && m->n_cigar) emit_mc_tag(&aux, &aux_len, &aux_cap, opt, m, which);
    /* MQ: the mate's MAPQ. Gated on `m` alone -- NOT `m->n_cigar` like MC:Z
     * above -- so it is emitted even when the mate is unmapped, matching
     * upstream exactly (bwamem.cpp:3484, bam_writer.cpp:395, meth_bam.cpp:585,
     * which all use `m && opt->compat->emit_mq`). Emitted here, between MC and
     * AS, because this crate's aux order mirrors upstream's write order.
     *
     * Not a bwa-mem3 invention: bwa emits MQ (bwamem.c:935, lh3/bwa#330) and
     * bwa-mem2 does not, having forked at 0.7.17 before that landed -- hence
     * the compat switch, whose default target (COMPAT_TARGET_OFF) sets
     * emit_mq = 1. */
    if (m && opt->compat != NULL && opt->compat->emit_mq)
        aux_put_i(&aux, &aux_len, &aux_cap, "MQ", m->mapq);
    if (p->score >= 0) aux_put_i(&aux, &aux_len, &aux_cap, "AS", p->score);
    if (p->sub >= 0)   aux_put_i(&aux, &aux_len, &aux_cap, "XS", p->sub);
    if (bwa_rg_id[0])  aux_put_Z(&aux, &aux_len, &aux_cap, "RG", bwa_rg_id);
    /* SA: if this is a primary (flag 0x100 not set) and other primary hits exist */
    if (!(p->flag & 0x100) && n_list > 1) {
        int i;
        for (i = 0; i < n_list; ++i)
            if (i != which && !(list[i].flag & 0x100)) break;
        if (i < n_list) emit_sa_tag(&aux, &aux_len, &aux_cap, bns, list, n_list, which);
    }
    if (p->XA) aux_put_Z(&aux, &aux_len, &aux_cap, "XA", p->XA);
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
    if (p->HN >= 0 && opt->compat != NULL && opt->compat->emit_hn)
        aux_put_i(&aux, &aux_len, &aux_cap, "HN", p->HN);

    /* D3 (--meth) Bismark tags. XR:Z (read conversion) on every record; XG:Z
     * (genome strand) and XM:Z (per-base methylation call) on mapped records.
     * XG / is_top_strand come from the winning hypothesis (p->meth_hypothesis:
     * OT/1 => top/CT, OB/0 or -1 => bottom/GA), NOT the 0x10 RC flag — matching
     * upstream meth_bam.cpp. XM is built from the ORIGINAL (unconverted) read
     * bases in emitted-SEQ orientation. */
    if (opt->meth_mode) {
        aux_put_Z(&aux, &aux_len, &aux_cap, "XR", is_r2 ? "GA" : "CT");
        if (p->rid >= 0) {
            int is_top = (p->meth_hypothesis >= 0 && (p->meth_hypothesis & 1)) ? 1 : 0;
            aux_put_Z(&aux, &aux_len, &aux_cap, "XG", is_top ? "CT" : "GA");
            if (s->meth_orig_seq && emit_len > 0 && n_cigar > 0) {
                char *seq_text = (char *) malloc((size_t)emit_len + 1);
                uint32_t *bam_cig = (uint32_t *) malloc((size_t)n_cigar * sizeof(uint32_t));
                if (seq_text && bam_cig) {
                    for (int i = 0; i < emit_len; ++i) {
                        seq_text[i] = p->is_rev
                            ? ascii_complement(s->meth_orig_seq[l_seq - 1 - (seq_start + i)])
                            : s->meth_orig_seq[seq_start + i];
                    }
                    seq_text[emit_len] = '\0';
                    /* Same clip rewrite as the emitted CIGAR below: the XM:Z
                     * builder walks this alongside `seq_text`, which is the
                     * hard-clipped span, so the two must agree on whether the
                     * clip consumes query bases. */
                    for (int i = 0; i < n_cigar; ++i)
                        bam_cig[i] = bwa_cigar_to_bam(
                            bwa_apply_clip_mode(p->cigar[i], opt, p->is_alt, which));
                    /* v0.9.0 added the chemistry selector. It flips only the
                     * methylated/unmethylated polarity of the call, never the
                     * CpG/CHG/CHH context classification, and comes off
                     * mem_opt_t exactly as upstream's own writer sources it
                     * (meth_bam.cpp:526). mem_opt_init defaults it to
                     * METH_CHEM_EMSEQ, so this is unchanged behavior unless a
                     * caller sets it. */
                    char *xm = meth_build_xm(bns, pac, p->rid, (int64_t)p->pos,
                                             is_top, bam_cig, n_cigar, seq_text, emit_len,
                                             (meth_chem_t) opt->meth_chem);
                    if (xm) aux_put_Z(&aux, &aux_len, &aux_cap, "XM", xm);
                }
                free(seq_text);
                free(bam_cig);
            }
        }
    }

    /* Reverse-complement the bwa-2bit-encoded query if is_rev; otherwise
     * point at the forward slice. Quality scores mirror the sequence. */
    uint8_t *emit_seq_buf = nullptr;
    char *emit_qual_buf = nullptr;
    const uint8_t *emit_seq = nullptr;
    const char *emit_qual = nullptr;
    if (emit_len > 0) {
        if (p->is_rev) {
            emit_seq_buf = (uint8_t *) malloc(emit_len);
            xassert(emit_seq_buf != NULL, "out of memory: emit_seq_buf");
            for (int i = 0; i < emit_len; ++i) {
                uint8_t b = (uint8_t) s->seq[l_seq - 1 - (seq_start + i)];
                emit_seq_buf[i] = bwa2_complement(b);
            }
            emit_seq = emit_seq_buf;
            if (s->qual) {
                emit_qual_buf = (char *) malloc(emit_len);
                xassert(emit_qual_buf != NULL, "out of memory: emit_qual_buf");
                for (int i = 0; i < emit_len; ++i)
                    emit_qual_buf[i] = s->qual[l_seq - 1 - (seq_start + i)];
                emit_qual = emit_qual_buf;
            }
        } else {
            emit_seq = (const uint8_t *) s->seq + seq_start;
            if (s->qual) emit_qual = s->qual + seq_start;
        }
    }

    size_t seq_packed = (size_t)(emit_len + 1) / 2;
    size_t block_size = (size_t)32 + l_read_name + 4 * (size_t)n_cigar
                      + seq_packed + (size_t)emit_len + aux_len;

    /* Build the record BODY (no u32 block_size prefix) into the scratch's
     * reusable record buffer; the sink prepends the prefix if it wants one. */
    if (block_size > e->sc->rec_cap) {
        while (block_size > e->sc->rec_cap) e->sc->rec_cap *= 2;
        uint8_t *tmp = (uint8_t *) realloc(e->sc->rec_buf, e->sc->rec_cap);
        xassert(tmp != NULL, "out of memory: rec_buf");
        e->sc->rec_buf = tmp;
    }
    uint8_t *w = e->sc->rec_buf;

    int32_t ref_id = eff_rid;
    memcpy(w, &ref_id, 4); w += 4;
    int32_t pos32 = (int32_t) eff_pos;
    memcpy(w, &pos32, 4); w += 4;

    *w++ = (uint8_t) l_read_name;
    *w++ = (uint8_t) p->mapq;

    uint16_t bin = (ref_id < 0 || ref_len == 0)
                       ? 4680
                       : reg2bin((int)eff_pos, (int)eff_pos + ref_len);
    memcpy(w, &bin, 2); w += 2;

    uint16_t nc = (uint16_t) n_cigar;
    memcpy(w, &nc, 2); w += 2;

    /* The packed FLAG is 16 bits, but pair_and_emit marks MEM_F_NO_MULTI split
     * hits with upstream's internal 0x10000 (bit 16). Remap it to the BAM
     * secondary bit 0x100 here, exactly as mem_aln2sam does at write time.
     * Keeping the marker at 0x10000 internally (not 0x100) is deliberate: it
     * lets emit_sa_tag still list NO_MULTI splits in SA:Z (its skip test is
     * `flag & 0x100`), matching upstream. */
    /* 0x10/0x20 are recomputed from the POST-copy strands (upstream derives
     * them after the rewrite above), so an unmapped record placed at its
     * mate's coordinates inherits that mate's strand. Clearing first matters:
     * pair_and_emit set them from the raw values. */
    uint32_t flag_raw = (p->flag & ~(uint32_t)(0x10 | 0x20));
    if (eff_is_rev)       flag_raw |= 0x10;
    if (m && mate_is_rev) flag_raw |= 0x20;
    uint16_t flag16 = (uint16_t)((flag_raw & 0xffff) | ((flag_raw & 0x10000) ? 0x100 : 0));
    memcpy(w, &flag16, 2); w += 2;

    int32_t ls = emit_len;
    memcpy(w, &ls, 4); w += 4;

    int32_t mtid = m ? mate_rid : -1;
    int32_t mpos = m ? (int32_t) mate_pos : -1;
    int32_t tlen = m ? compute_tlen(p, ref_len, m, cigar_ref_len(m->n_cigar, m->cigar)) : 0;
    memcpy(w, &mtid, 4); w += 4;
    memcpy(w, &mpos, 4); w += 4;
    memcpy(w, &tlen, 4); w += 4;

    memcpy(w, s->name, l_read_name - 1); w += l_read_name - 1;
    *w++ = 0;

    for (int i = 0; i < n_cigar; ++i) {
        uint32_t c = bwa_cigar_to_bam(
            bwa_apply_clip_mode(p->cigar[i], opt, p->is_alt, which));
        memcpy(w, &c, 4); w += 4;
    }

    /* 4-bit packed seq. emit_seq bytes are bwa's 2-bit encoding (0-4); remap
     * to BAM 4-bit nibbles (1/2/4/8/15). */
    for (int i = 0; i < emit_len; i += 2) {
        uint8_t hi = emit_seq ? bwa2_to_bam4[emit_seq[i] & 7] : 15;
        uint8_t lo = 0;
        if (i + 1 < emit_len) {
            lo = emit_seq ? bwa2_to_bam4[emit_seq[i + 1] & 7] : 15;
        }
        *w++ = (uint8_t)((hi << 4) | lo);
    }

    /* qual: ASCII - 33, or 0xFF if missing. */
    if (emit_qual) {
        for (int i = 0; i < emit_len; ++i) *w++ = (uint8_t)(emit_qual[i] - 33);
    } else if (emit_len > 0) {
        memset(w, 0xFF, emit_len); w += emit_len;
    }

    if (aux_len) { memcpy(w, aux, aux_len); w += aux_len; }
    /* Retain the (possibly grown) aux buffer for the next record instead of
     * freeing it -- buf_append reallocs in place. */
    e->sc->aux_buf = aux; e->sc->aux_cap = aux_cap;

    free(emit_seq_buf);
    free(emit_qual_buf);
    assert((size_t)(w - e->sc->rec_buf) == block_size);
    e->sink(e->ctx, e->origin_kind, origin_idx, e->sc->rec_buf, block_size);
}

/* ------------------ Phase 1: scratch + fused seed+extend ------------------ */

ShimScratch *shim_scratch_new(void) {
    ShimScratch *sc = (ShimScratch *) calloc(1, sizeof(ShimScratch));
    if (!sc) return nullptr;
    worker_alloc(sc->w, 1);
    sc->rec_cap = 4096; sc->rec_buf = (uint8_t *) malloc(sc->rec_cap);
    sc->aux_cap = 1024; sc->aux_buf = (uint8_t *) malloc(sc->aux_cap);
    /* The emission path grows rec_buf/aux_buf from these caps and dereferences
     * them without re-checking, so a failed malloc here must not yield a
     * scratch with a non-zero cap but a NULL buffer. Fail the whole alloc. */
    if (!sc->rec_buf || !sc->aux_buf) {
        worker_free(sc->w, 1);
        free(sc->rec_buf); free(sc->aux_buf);
        free(sc);
        return nullptr;
    }
    return sc;
}

void shim_scratch_free(ShimScratch *sc) {
    if (!sc) return;
    worker_free(sc->w, 1);
    free(sc->rec_buf); free(sc->aux_buf);
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

/* Fused seed + SE-extend over seqs[first, first+n) in BATCH_SIZE chunks, the
 * shape of upstream's worker_bwt_aln (bwamem.cpp:2782-2786). `opt` is the
 * per-group copy (MEM_F_PE set for pairs). Chain/seed windows are chunk-local
 * so the BATCH_SIZE scratch suffices; the pre-0.9.0 tail `seedBufSz` shrink is
 * dropped as upstream did ("output-dead either way", worker_bwt). */
static void seed_extend_range(ShimScratch *sc, ShimRegs *r, BwaShimIndex *idx,
                              const mem_opt_t *opt, int first, int n)
{
    worker_t &w = sc->w;
    FMI_search *fmi = idx->fmi;
    w.opt = opt; w.fmi = fmi; w.seqs = r->seqs; w.regs = r->regs;
    w.n_processed = 0; w.pes = nullptr; w.nreads = r->n_seqs;
    w.meth_orig_bns = idx->meth_orig_bns;
    w.meth_orig_pac = idx->meth_orig_pac;
    w.meth_orig_ref_string = idx->meth_orig_ref_string;
    w.ref_string = idx->meth_orig_ref_string ? idx->meth_orig_ref_string : idx->ref_string;
    for (int seq_id = first; seq_id < first + n; seq_id += BATCH_SIZE) {
        int bs = first + n - seq_id;
        if (bs > BATCH_SIZE) bs = BATCH_SIZE;
        mem_kernel1_core(fmi, opt, r->seqs + seq_id, bs,
                         w.chain_scratch, w.seed_scratch, w.seed_scratch_size,
                         &w.mmc, 0 /* tid */, idx->meth_orig_bns, idx->meth_orig_pac);
        mem_kernel2_core(fmi, opt, r->seqs + seq_id, r->regs + seq_id, bs,
                         w.chain_scratch, &w.mmc, w.ref_string, 0 /* tid */,
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
    if (r->n_pairs > 0)   seed_extend_range(sc, r, idx, &opt_pe, 0, (int)(2 * r->n_pairs));
    if (r->n_singles > 0) seed_extend_range(sc, r, idx, &opt_se, (int)(2 * r->n_pairs), (int)r->n_singles);
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

/* Defined below (just before shim_extend_batch); forward-declared so
 * shim_pair_emit can drive it. */
static void pair_and_emit(ShimEmit *e, size_t origin_idx, uint64_t id,
                          const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                          bseq1_t *s, mem_alnreg_v *a, const mem_pestat_t pes[4]);

/* Defined after emit_resolved_pair; forward-declared so shim_pair_emit can
 * drive it after the pair loop. */
static void single_and_emit(ShimEmit *e, size_t origin_idx, uint64_t id,
                            const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                            bseq1_t *s, mem_alnreg_v *a);

/* Defined below (with pair_resolve_scalar); forward-declared so the batched
 * mate-rescue path in shim_pair_emit can emit each resolved pair. */
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
int shim_pestat_cohort(void *idx_opaque, const mem_opt_t *opts,
                       const ShimRegs *const *regs, size_t n_regs, mem_pestat_t *out)
{
    if (!idx_opaque || !opts || (n_regs > 0 && !regs) || !out) return -1;
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(idx_opaque);
    const bntseq_t *bns = idx->meth_orig_bns ? idx->meth_orig_bns : idx->fmi->idx->bns;
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
    mem_opt_t opt_pe = *opts; opt_pe.n_threads = 1; opt_pe.flag |= MEM_F_PE;
    mem_pestat(&opt_pe, bns->l_pac, (int)total, cat, out);
    free(cat);
    return 0;
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
    ShimEmit e = { sc, sink, ctx, 0u /* BWA_ORIGIN_PAIR */ };
#if BWAMEM_BATCHED_MATESW
    /* The CLI's worker_sam path on AVX2/AVX-512/NEON (bwamem.cpp:2826-2879):
     * gather every pair's rescue jobs, run them through the SIMD kswv kernel
     * once, then resolve + emit per pair. tid = 0: one ShimScratch per thread. */
    if (r->n_pairs > 0) {
        worker_t &w = sc->w;
        /* Chunk the batched mate-rescue exactly as the CLI's kt_for dispatches
         * worker_sam: BATCH_SIZE reads (= BATCH_SIZE/2 pairs) per work item
         * (kthread.cpp:109-118, bwamem.cpp worker_sam ~L2788). The pre-loop's
         * running pcnt/gcnt/maxRefLen/maxQerLen and, with them, the kswv seqBuf
         * ref-window offset (SeqPair.idr, an int32 derived from the monotonic
         * pcnt) reset at every chunk boundary, so the offset can never exceed
         * int32 no matter how large the -K cohort or how long the reads
         * (seqbuf_grow_capacity's SEQBUF_CAPACITY_OVERFLOW → seqbuf_capacity_fatal
         * exit(), which would abort across the FFI boundary). Without the chunk
         * loop a single running pcnt spanned all r->n_pairs, exactly the
         * unbounded-offset the CLI avoids by resetting per BATCH_SIZE.
         *
         * Chunking is purely internal kswv/seqBuf batch bookkeeping. The
         * per-pair rescue result is independent of how pairs are grouped, and
         * the GLOBAL read-ordinal id (ids.first_pair_id + global pair index) and
         * the emit order (global index i, ascending) are unchanged — only
         * pcnt/gcnt/maxRefLen/maxQerLen/aln/myaln reset per chunk. So output
         * stays byte-identical to the former single-batch path and to
         * `bwa-mem3 mem -t 1 -p`. For n_pairs <= BATCH_SIZE/2 this is a single
         * iteration, identical to before. */
        const size_t pairs_per_chunk = (size_t)BATCH_SIZE / 2;
        for (size_t chunk_start = 0; chunk_start < r->n_pairs; chunk_start += pairs_per_chunk) {
            size_t chunk_end = chunk_start + pairs_per_chunk;
            if (chunk_end > r->n_pairs) chunk_end = r->n_pairs;

            int32_t maxRefLen = 0, maxQerLen = 0, gcnt = 0;
            int64_t pcnt = 0;
            for (size_t i = chunk_start; i < chunk_end; ++i)
                mem_sam_pe_batch_pre(&opt_pe, bns, pac, pestat, ids.first_pair_id + (uint64_t)i,
                                     r->seqs + 2*i, r->regs + 2*i, &w.mmc, pcnt, gcnt,
                                     maxRefLen, maxQerLen, 0);
            int64_t pcnt8 = sort_classify(&w.mmc, pcnt, 0);
            kswr_t *aln = (kswr_t *) _mm_malloc((pcnt + SIMD_WIDTH8) * sizeof(kswr_t), 64);
            xassert(aln != NULL, "out of memory: aln");
            mem_sam_pe_batch(&opt_pe, &w.mmc, pcnt, pcnt8, aln, maxRefLen, maxQerLen, 0);
            gcnt = 0;
            kswr_t *myaln = aln;
            for (size_t i = chunk_start; i < chunk_end; ++i) {
                int n_pri[2], z[2], q_se[2], extra_flag, paired;
                mem_pair_resolve_batch_post(&opt_pe, bns, pac, pestat, ids.first_pair_id + (uint64_t)i,
                                            r->seqs + 2*i, r->regs + 2*i, &myaln, &w.mmc, gcnt, 0,
                                            n_pri, z, q_se, &extra_flag, &paired);
                emit_resolved_pair(&e, i, &opt_pe, bns, pac, r->seqs + 2*i, r->regs + 2*i, pestat,
                                   n_pri, z, q_se, extra_flag, paired);
            }
            _mm_free(aln);
        }
    }
#else
    for (size_t i = 0; i < r->n_pairs; ++i)
        pair_and_emit(&e, i, ids.first_pair_id + (uint64_t)i,
                      &opt_pe, bns, pac, r->seqs + 2*i, r->regs + 2*i, pestat);
#endif
    /* Singles after pairs: SE group emitted with MEM_F_PE cleared, ids from
     * first_single_id, one origin_idx per single (index into batch->singles).
     * BWA_ORIGIN_SINGLE (1) distinguishes them from pairs at the sink. */
    mem_opt_t opt_se = *opts; opt_se.n_threads = 1; opt_se.flag &= ~MEM_F_PE;
    e.origin_kind = 1u /* BWA_ORIGIN_SINGLE */;
    for (size_t i = 0; i < r->n_singles; ++i) {
        size_t k = 2 * r->n_pairs + i;
        single_and_emit(&e, i, ids.first_single_id + (uint64_t)i,
                        &opt_se, bns, pac, r->seqs + k, r->regs + k);
    }
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

/* Resolve half of the former pair_and_emit: run upstream's full pairing
 * decision (mate-rescue SW + mem_mark_primary_se + optional MEM_F_PRIMARY5
 * reorder + mem_pair + is_multi + q_pe/q_se + secondary<->primary switch),
 * writing the resolved indices/flags into the out-params. Split out so Task 10
 * can swap the scalar resolve for the batched one. On the paired branch
 * extra_flag already includes 0x2 (if the paired alignment was preferred); on
 * the no_pairing branch the emit half ORs it in after running mem_infer_dir
 * itself, matching mem_sam_pe's no_pairing block. */
static void pair_resolve_scalar(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                                const mem_pestat_t pes[4], uint64_t id,
                                bseq1_t s[2], mem_alnreg_v a[2],
                                int n_pri[2], int z[2], int q_se[2],
                                int *extra_flag, int *paired)
{
    *extra_flag = 0; *paired = 0;
    mem_pair_resolve(opt, bns, pac, pes, id, s, a, n_pri, z, q_se, extra_flag, paired);
}

/* Emission half of the former pair_and_emit: everything after mem_pair_resolve.
 * Split out so Task 10 can swap the scalar resolve for the batched one. */
static void emit_resolved_pair(ShimEmit *e, size_t origin_idx,
                               const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                               bseq1_t *s, mem_alnreg_v *a, const mem_pestat_t pes[4],
                               const int n_pri[2], const int z[2], const int q_se[2],
                               int extra_flag, int paired)
{
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
     * requires the mem_mark_primary_se that mem_pair_resolve ran) rather than
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
            /* D3 (--meth): pass the ORIGINAL read bases so mem_reg2aln
             * regenerates CIGAR/NM/MD against the original ref (NULL and a
             * no-op outside --meth). */
            lists[k][j] = mem_reg2aln(opt, bns, pac, s[k].l_seq, s[k].seq,
                                      ar, s[k].meth_orig_seq);
            lists[k][j].XA = XA[k] ? XA[k][j] : nullptr;
            lists[k][j].HN = HN[k] ? HN[k][j] : -1;
            /* mem_reg2sam emit filter. */
            bool e = ar->score >= opt->T;
            if (e && ar->secondary >= 0 && (ar->is_alt || !(opt->flag & MEM_F_ALL)))
                e = false;
            if (e && ar->secondary >= 0 && ar->secondary < INT_MAX
                && ar->score < a[k].a[ar->secondary].score * opt->drop_ratio)
                e = false;
            /* Paired branch, z[k] != 0: mem_pair_resolve promoted the
             * paired-selected region a[k].a[z[k]] (secondary set to -2) and ran
             * the secondary_all switch, which reassigns the old SE-primary
             * (region 0) to z[k]'s group — leaving it with secondary < 0 but
             * secondary_all >= 0. mem_gen_alt folds that region into z[k]'s XA:Z,
             * and upstream mem_sam_pe's paired block emits ONLY z[k] as primary,
             * never the switched-away region. Our emit filter keys off
             * `secondary` alone, so without this it would surface the old primary
             * as an extra record and demote z[k] to a 0x800 supplementary. Drop
             * the folded region so array order leaves z[k] first (primary). Never
             * fires when z[k] == 0 (no switch) or on the no_pairing branch (which
             * returns before the switch, leaving secondary_all == secondary). */
            if (e && paired && j != z[k] && ar->secondary < 0
                && ar->secondary_all >= 0)
                e = false;
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
     * when mem_pair_resolve didn't take the paired branch, decide the
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
    if (!paired && (opt->flag & MEM_F_PE) && !(opt->flag & MEM_F_NOPAIRING)
        && which[0] >= 0 && which[1] >= 0
        && h_rid[0] == h_rid[1] && h_rid[0] >= 0) {
        extra_flag |= mem_proper_pair_extra_flag(opt, bns->l_pac, a, which, pes);
    }

    /* Apply extra_flag (0x1 paired + optional 0x2 proper-pair) to every
     * emitted record, matching mem_reg2sam's behavior in the no_pairing
     * branch. On the paired branch this over-applies 0x2 to non-primary
     * records relative to mem_sam_pe's stricter primary-only application;
     * that matches how downstream tools treat 0x2 as a pair-level flag. */
    for (int k = 0; k < 2; ++k) {
        for (int j = 0; j < n_lists[k]; ++j) {
            lists[k][j].flag |= extra_flag;
        }
    }

    /* On the paired branch, apply the q_se mapq to the chosen primary
     * (mirrors `h[i].mapq = q_se[i]` in mem_sam_pe's paired emission). */
    if (paired) {
        for (int k = 0; k < 2; ++k) {
            if ((int)a[k].n > 0 && z[k] < n_lists[k]) {
                lists[k][z[k]].mapq = q_se[k];
            }
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
                 * internal marker; append_bam_record remaps it to 0x100 on write. */
                p.flag |= (opt->flag & MEM_F_NO_MULTI) ? 0x10000 : 0x800;
                if (!(opt->flag & MEM_F_KEEP_SUPP_MAPQ) && p.mapq > aa[k][0].mapq)
                    p.mapq = aa[k][0].mapq;
            }
            aa[k][l++] = p;
        }
        n_aa[k] = l;
    }

    /* Emit records. For each side k, the pair mate is the paired-selected
     * primary of the other side (z[!k] on the paired branch, 0 otherwise —
     * matching bwamem_pair.cpp:545-556's use of &a[!i].a[z[!i]] as the mate
     * anchor). Without this the paired branch would use lists[!k][0] even when
     * mem_pair_resolve picked a non-zero primary, driving RNEXT/PNEXT/TLEN/MC
     * off the wrong mate alignment. */
    for (int k = 0; k < 2; ++k) {
        int mate_idx = 0;
        if (paired && z[!k] >= 0 && z[!k] < n_lists[!k]) mate_idx = z[!k];
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

/* One pair: resolve then emit. `id` is the GLOBAL pair ordinal
 * ((n_processed >> 1) + pos in worker_sam, bwamem.cpp:2819); it seeds the
 * hash_64 tie-breaks in mem_mark_primary_se / mem_pair, so a different id can
 * legitimately pick a different primary among equal-score hits. */
static void pair_and_emit(ShimEmit *e, size_t origin_idx, uint64_t id,
                          const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                          bseq1_t *s, mem_alnreg_v *a, const mem_pestat_t pes[4])
{
    int n_pri[2], z[2], q_se[2], extra_flag, paired;
    pair_resolve_scalar(opt, bns, pac, pes, id, s, a, n_pri, z, q_se, &extra_flag, &paired);
    emit_resolved_pair(e, origin_idx, opt, bns, pac, s, a, pes, n_pri, z, q_se, extra_flag, paired);
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

} /* extern "C" */
