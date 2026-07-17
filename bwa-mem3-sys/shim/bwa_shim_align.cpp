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
#include "FMI_search.h"

/* worker_alloc / worker_free — per-worker scratch (chaining arrays, BSW
 * buffers, lazy SMEM buffers) sized for the batch. Upstream defines these in
 * `fastmap.cpp`, but as of bwa-mem3 0.6.0 that TU is transitively coupled to
 * the new `fast_reader` FASTQ path (libdeflate + zlib-ng) and to numa/htslib,
 * none of which this crate compiles (the Rust CLI does its own I/O). We
 * therefore do not build `fastmap.cpp` (see `bwa-mem3-sys/build.rs`) and carry
 * file-local copies of just these two allocators here.
 *
 * The allocation set below is verbatim from upstream's `worker_alloc` /
 * `worker_free`, with ONE deliberate divergence: upstream's per-call
 * "Memory pre-allocation" stderr diagnostics are dropped. Upstream calls
 * worker_alloc once per thread at process start; this crate calls it once per
 * `shim_seed_batch` (per read chunk), so those prints would spam stderr on
 * every batch.
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
static void worker_alloc(const mem_opt_t *opt, worker_t &w, int32_t nreads, int32_t nthreads)
{
    assert(opt != NULL);
    assert(nreads >= 0);
    assert(nthreads > 0);

    // Record the thread count on the worker so worker_free can validate the
    // paired call and the per-thread loops can never walk past the slots
    // populated here.
    w.nthreads = nthreads;

    int32_t memSize = nreads;

    /* Mem allocation section for core kernels */
    w.regs = NULL; w.chain_ar = NULL; w.seedBuf = NULL;

    w.regs = (mem_alnreg_v *) calloc(memSize, sizeof(mem_alnreg_v));
    w.chain_ar = (mem_chain_v*) malloc (memSize * sizeof(mem_chain_v));
    w.seedBuf = (mem_seed_t *) calloc(sizeof(mem_seed_t),  memSize * AVG_SEEDS_PER_READ);

    assert(w.seedBuf  != NULL);
    assert(w.regs     != NULL);
    assert(w.chain_ar != NULL);

    w.seedBufSize = BATCH_SIZE * AVG_SEEDS_PER_READ;

    /* SWA mem allocation */
    int64_t wsize = BATCH_SIZE * SEEDS_PER_READ;
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

        assert(w.mmc.seqBufLeftRef[l*CACHE_LINE]  != NULL);
        assert(w.mmc.seqBufLeftQer[l*CACHE_LINE]  != NULL);
        assert(w.mmc.seqBufRightRef[l*CACHE_LINE] != NULL);
        assert(w.mmc.seqBufRightQer[l*CACHE_LINE] != NULL);
    }

    for(int l=0; l<nthreads; l++) {
        w.mmc.seqPairArrayAux[l]      = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.seqPairArrayLeft128[l]  = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.seqPairArrayRight128[l] = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.wsize[l] = wsize;

        assert(w.mmc.seqPairArrayAux[l] != NULL);
        assert(w.mmc.seqPairArrayLeft128[l] != NULL);
        assert(w.mmc.seqPairArrayRight128[l] != NULL);
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

    free(w.chain_ar);
    free(w.regs);
    free(w.seedBuf);

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

/* Phase-1 state carried into extend_batch. Owns worker and seqs; borrows
 * ref_string from the index handle. extend_batch (or estimate_pestat) frees
 * the owned pieces. */
struct ShimSeeds {
    worker_t w;
    bseq1_t *seqs;
    int      n_seqs;
    mem_opt_t *opts;
    FMI_search *fmi;
    uint8_t *ref_string;  /* borrowed from BwaShimIndex; not freed here */
};

/* Index handle: owns the loaded FMI_search plus the one-time-unpacked
 * 2*l_pac reference string used by mem_chain2aln_across_reads_V2. The
 * ref_string is built once at load time so every seed_batch borrows it
 * instead of rebuilding (~6GB alloc + full pass over packed ref per call
 * on hs38). Forward declaration of build_ref_string below. */
static uint8_t *build_ref_string(const bntseq_t *bns, const uint8_t *pac);

struct BwaShimIndex {
    FMI_search *fmi;
    uint8_t    *ref_string;
};

/* ------------------ Index ------------------ */

void *shim_align_idx_load(const char *prefix) {
    FMI_search *fmi = new FMI_search(prefix);
    fmi->load_index();
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

void shim_align_idx_free(void *opaque) {
    if (!opaque) return;
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(opaque);
    if (idx->ref_string) _mm_free(idx->ref_string);
    delete idx->fmi;
    free(idx);
}

size_t shim_align_idx_n_contigs(void *opaque) {
    return (size_t) static_cast<BwaShimIndex *>(opaque)->fmi->idx->bns->n_seqs;
}

const char *shim_align_idx_contig_name(void *opaque, size_t i) {
    return static_cast<BwaShimIndex *>(opaque)->fmi->idx->bns->anns[i].name;
}

int64_t shim_align_idx_contig_len(void *opaque, size_t i) {
    return static_cast<BwaShimIndex *>(opaque)->fmi->idx->bns->anns[i].len;
}

/* ------------------ Helpers ------------------ */

static bseq1_t *copy_pairs_to_seqs(const ShimReadPair *pairs, size_t n_pairs) {
    int nseqs = (int)(2 * n_pairs);
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
            s->name = (char *) malloc(nl[k] + 1);
            memcpy(s->name, nm[k], nl[k]);
            s->name[nl[k]] = '\0';
            s->seq = (char *) malloc(sl[k] + 1);
            memcpy(s->seq, sq[k], sl[k]);
            s->seq[sl[k]] = '\0';
            if (ql[k]) {
                s->qual = (char *) malloc(sl[k] + 1);
                memcpy(s->qual, ql[k], sl[k]);
                s->qual[sl[k]] = '\0';
            } else {
                s->qual = nullptr;
            }
            s->sam = nullptr;
        }
    }
    return seqs;
}

static void free_seqs(bseq1_t *seqs, int nseqs) {
    if (!seqs) return;
    for (int i = 0; i < nseqs; ++i) {
        free(seqs[i].name);
        free(seqs[i].seq);
        free(seqs[i].qual);
        free(seqs[i].sam);
    }
    free(seqs);
}

/* Unpack 2-bit packed reference into a 1-byte-per-base array with the
 * reverse-complement appended (total = 2 * l_pac). Required by
 * mem_chain2aln_across_reads_V2. */
static uint8_t *build_ref_string(const bntseq_t *bns, const uint8_t *pac) {
    int64_t ref_len = bns->l_pac * 2;
    uint8_t *s = (uint8_t *) _mm_malloc(ref_len, 64);
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
        *buf = (uint8_t *) realloc(*buf, *cap);
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

/* Emit MC:Z aux (mate CIGAR) from a mate mem_aln_t. */
static void emit_mc_tag(uint8_t **buf, size_t *len, size_t *cap,
                        const mem_aln_t *m) {
    uint8_t header[3] = {'M', 'C', 'Z'};
    buf_append(buf, len, cap, header, 3);
    char tmp[32];
    for (int i = 0; i < m->n_cigar; ++i) {
        int n = snprintf(tmp, sizeof(tmp), "%u%c",
                         m->cigar[i] >> 4, "MIDSH"[m->cigar[i] & 0xf]);
        buf_append(buf, len, cap, tmp, n);
    }
    uint8_t zero = 0;
    buf_append(buf, len, cap, &zero, 1);
}

/* Compute tlen between two aligned mem_aln_t. */
static int32_t compute_tlen(const mem_aln_t *a, int a_ref_len,
                            const mem_aln_t *b, int b_ref_len) {
    if (!a || !b || a->rid < 0 || b->rid < 0 || a->rid != b->rid) return 0;
    int64_t a_begin = a->pos;
    int64_t b_begin = b->pos;
    int64_t a_end = a->pos + a_ref_len;
    int64_t b_end = b->pos + b_ref_len;
    int64_t begin = a_begin < b_begin ? a_begin : b_begin;
    int64_t end   = a_end   > b_end   ? a_end   : b_end;
    int64_t tlen  = end - begin;
    return a_begin <= b_begin ? (int32_t)tlen : -(int32_t)tlen;
}

/* Append one packed BAM record to `out`'s buffer. */
static void append_bam_record(ShimAlignOutput *out, size_t pair_idx,
                              const bntseq_t *bns, const bseq1_t *s,
                              const mem_aln_t *p, int n_list,
                              const mem_aln_t *list, int which,
                              const mem_aln_t *m)
{
    int l_read_name = (int) strlen(s->name) + 1;
    int l_seq = s->l_seq;
    int n_cigar = p->n_cigar;
    int ref_len = cigar_ref_len(n_cigar, p->cigar);

    /* Build aux first so we know its size. */
    uint8_t *aux = nullptr;
    size_t aux_len = 0, aux_cap = 0;
    if (p->n_cigar) {
        aux_put_i(&aux, &aux_len, &aux_cap, "NM", p->NM);
        /* MD string is stored right after the CIGAR array in p->cigar. */
        const char *md = (const char *)(p->cigar + p->n_cigar);
        aux_put_Z(&aux, &aux_len, &aux_cap, "MD", md);
    }
    if (m && m->n_cigar) emit_mc_tag(&aux, &aux_len, &aux_cap, m);
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

    /* Choose seq/qual source. For supplementary (flag 0x800) and secondary
     * (0x100) without soft-clip-supplementary opt, bwa emits hard-clipped:
     * the sequence and qual are empty-or-subset. We honor p->n_cigar's
     * H ops by emitting only the non-H-clipped range from s->seq. */
    int seq_start = 0;
    int seq_end = l_seq;
    if (n_cigar > 0) {
        /* bwa-mem3 H opcode = 4 (not 5 as in BAM spec). */
        if ((p->cigar[0] & 0xf) == 4) {            /* leading H */
            seq_start = (int)(p->cigar[0] >> 4);
        }
        if ((p->cigar[n_cigar - 1] & 0xf) == 4) {  /* trailing H */
            seq_end -= (int)(p->cigar[n_cigar - 1] >> 4);
        }
    }
    int emit_len = seq_end - seq_start;
    if (emit_len < 0) emit_len = 0;

    /* Reverse-complement the bwa-2bit-encoded query if is_rev; otherwise
     * point at the forward slice. Quality scores mirror the sequence. */
    uint8_t *emit_seq_buf = nullptr;
    char *emit_qual_buf = nullptr;
    const uint8_t *emit_seq = nullptr;
    const char *emit_qual = nullptr;
    if (emit_len > 0) {
        if (p->is_rev) {
            emit_seq_buf = (uint8_t *) malloc(emit_len);
            for (int i = 0; i < emit_len; ++i) {
                uint8_t b = (uint8_t) s->seq[l_seq - 1 - (seq_start + i)];
                emit_seq_buf[i] = bwa2_complement(b);
            }
            emit_seq = emit_seq_buf;
            if (s->qual) {
                emit_qual_buf = (char *) malloc(emit_len);
                for (int i = 0; i < emit_len; ++i)
                    emit_qual_buf[i] = s->qual[l_seq - 1 - (seq_start + i)];
                emit_qual = emit_qual_buf;
            }
        } else {
            emit_seq = (const uint8_t *) s->seq + seq_start;
            if (s->qual) emit_qual = s->qual + seq_start;
        }
    }

    int seq_packed = (emit_len + 1) / 2;
    size_t block_size = (size_t)32 + l_read_name + 4 * (size_t)n_cigar
                      + (size_t)seq_packed + (size_t)emit_len + aux_len;
    size_t rec_size = 4 + block_size;

    /* Reserve in out->buf. */
    if (out->buf_len + rec_size > out->buf_cap) {
        while (out->buf_len + rec_size > out->buf_cap) out->buf_cap *= 2;
        out->buf = (uint8_t *) realloc(out->buf, out->buf_cap);
    }
    if (out->n_recs == out->cap) {
        out->cap *= 2;
        out->rec_off  = (size_t *) realloc(out->rec_off,  out->cap * sizeof(size_t));
        out->rec_len  = (size_t *) realloc(out->rec_len,  out->cap * sizeof(size_t));
        out->pair_idx = (size_t *) realloc(out->pair_idx, out->cap * sizeof(size_t));
    }

    uint8_t *w = out->buf + out->buf_len;
    size_t rec_off = out->buf_len;

    uint32_t bs32 = (uint32_t) block_size;
    memcpy(w, &bs32, 4); w += 4;

    int32_t ref_id = p->rid;
    memcpy(w, &ref_id, 4); w += 4;
    int32_t pos32 = (int32_t) p->pos;
    memcpy(w, &pos32, 4); w += 4;

    *w++ = (uint8_t) l_read_name;
    *w++ = (uint8_t) p->mapq;

    uint16_t bin = (ref_id < 0 || ref_len == 0) ? 4680
                                                : reg2bin((int)p->pos, (int)p->pos + ref_len);
    memcpy(w, &bin, 2); w += 2;

    uint16_t nc = (uint16_t) n_cigar;
    memcpy(w, &nc, 2); w += 2;

    uint16_t flag16 = (uint16_t) p->flag;
    memcpy(w, &flag16, 2); w += 2;

    int32_t ls = emit_len;
    memcpy(w, &ls, 4); w += 4;

    int32_t mtid = m ? m->rid  : -1;
    int32_t mpos = m ? (int32_t) m->pos : -1;
    int32_t tlen = m ? compute_tlen(p, ref_len, m, cigar_ref_len(m->n_cigar, m->cigar)) : 0;
    memcpy(w, &mtid, 4); w += 4;
    memcpy(w, &mpos, 4); w += 4;
    memcpy(w, &tlen, 4); w += 4;

    memcpy(w, s->name, l_read_name - 1); w += l_read_name - 1;
    *w++ = 0;

    for (int i = 0; i < n_cigar; ++i) {
        uint32_t c = bwa_cigar_to_bam(p->cigar[i]);
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
    free(aux);

    free(emit_seq_buf);
    free(emit_qual_buf);

    out->rec_off[out->n_recs]  = rec_off;
    out->rec_len[out->n_recs]  = rec_size;
    out->pair_idx[out->n_recs] = pair_idx;
    out->n_recs++;
    out->buf_len = rec_off + rec_size;
}

/* ------------------ Phase 1: seed_batch ------------------ */

ShimSeeds *shim_seed_batch(void *idx_opaque, mem_opt_t *opts,
                           const ShimReadPair *pairs, size_t n_pairs)
{
    if (!idx_opaque || !opts) return nullptr;
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(idx_opaque);
    FMI_search *fmi = idx->fmi;

    ShimSeeds *s = (ShimSeeds *) calloc(1, sizeof(ShimSeeds));
    s->opts = opts;
    s->fmi  = fmi;
    s->n_seqs = (int)(2 * n_pairs);
    s->seqs = copy_pairs_to_seqs(pairs, n_pairs);
    if (!s->seqs) { free(s); return nullptr; }

    /* Force single-thread; ensure PE flag. */
    opts->n_threads = 1;
    opts->flag |= MEM_F_PE;

    /* Allocate per-worker scratch using upstream's public helper so our
     * layout stays in sync with the main bwa-mem3 pipeline. Sets w.nthreads
     * internally and asserts on bad input / OOM; we fill in the non-scratch
     * fields afterward. */
    worker_alloc(opts, s->w, s->n_seqs, 1);
    s->w.opt = opts;
    s->w.nreads = s->n_seqs;
    s->w.fmi = fmi;
    s->w.seqs = s->seqs;
    s->w.n_processed = 0;
    s->w.pes = nullptr;

    /* Borrow the unpacked reference string from the index handle
     * (built once at load time; see BwaShimIndex / shim_align_idx_load). */
    s->ref_string = idx->ref_string;

    /* Seed in BATCH_SIZE chunks by calling mem_kernel1_core directly. */
    int memSize = s->n_seqs;
    for (int seq_id = 0; seq_id < s->n_seqs; seq_id += BATCH_SIZE) {
        int batch_size = s->n_seqs - seq_id;
        if (batch_size > BATCH_SIZE) batch_size = BATCH_SIZE;
        int seedBufSz = s->w.seedBufSize;
        if (batch_size < BATCH_SIZE) {
            seedBufSz = (memSize - seq_id) * AVG_SEEDS_PER_READ;
        }
        mem_kernel1_core(s->w.fmi, s->w.opt,
                         s->w.seqs + seq_id,
                         batch_size,
                         s->w.chain_ar + seq_id,
                         s->w.seedBuf + (size_t)seq_id * AVG_SEEDS_PER_READ,
                         seedBufSz,
                         &s->w.mmc,
                         0 /* tid */);
    }
    return s;
}

void shim_seeds_free(ShimSeeds *s) {
    if (!s) return;
    if (s->w.chain_ar) {
        /* chains allocated by mem_kernel1_core; their inner seed arrays live in
         * w.seedBuf which alloc_worker owns. No extra freeing needed here. */
    }
    free_seqs(s->seqs, s->n_seqs);
    /* s->ref_string is borrowed from BwaShimIndex; do not free. */
    worker_free(s->w, 1);
    free(s);
}

/* ------------------ Phase 2: extend_batch ------------------ */

static ShimAlignOutput *alloc_align_output(size_t n_pairs) {
    ShimAlignOutput *out = (ShimAlignOutput *) calloc(1, sizeof(ShimAlignOutput));
    out->cap = n_pairs > 0 ? n_pairs * 4 : 64;
    out->rec_off  = (size_t *) malloc(out->cap * sizeof(size_t));
    out->rec_len  = (size_t *) malloc(out->cap * sizeof(size_t));
    out->pair_idx = (size_t *) malloc(out->cap * sizeof(size_t));
    out->buf_cap = 4096;
    out->buf = (uint8_t *) malloc(out->buf_cap);
    return out;
}

/* Run SE extension into w.regs by calling mem_kernel2_core in BATCH_SIZE chunks. */
static void run_se_extension(ShimSeeds *s) {
    s->w.ref_string = s->ref_string;
    const bntseq_t *bns = s->fmi->idx->bns;
    const uint8_t *pac  = s->fmi->idx->pac;
    (void)bns; (void)pac;  /* passed via w.fmi->idx */
    for (int seq_id = 0; seq_id < s->n_seqs; seq_id += BATCH_SIZE) {
        int batch_size = s->n_seqs - seq_id;
        if (batch_size > BATCH_SIZE) batch_size = BATCH_SIZE;
        mem_kernel2_core(s->w.fmi, s->w.opt,
                         s->w.seqs + seq_id,
                         s->w.regs + seq_id,
                         batch_size,
                         s->w.chain_ar + seq_id,
                         &s->w.mmc,
                         s->w.ref_string,
                         0 /* tid */);
    }
}

/* Core pairing + BAM emission for one interleaved pair (r1=seqs[2i], r2=seqs[2i+1]). */
static void pair_and_emit(ShimAlignOutput *out, size_t pair_idx,
                          mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                          bseq1_t *s, mem_alnreg_v *a, const mem_pestat_t pes[4])
{
    /* Run upstream's full pairing decision: mate-rescue SW + mem_mark_primary_se
     * + optional MEM_F_PRIMARY5 reorder + mem_pair + is_multi + q_pe/q_se +
     * secondary<->primary switch. All verbatim from the bwa-mem3 pipeline — this
     * replaces the shim's previously hand-rolled composition. On the paired
     * branch extra_flag already includes 0x2 (if the paired alignment was
     * preferred); on the no_pairing branch we OR it in below after running
     * mem_infer_dir ourselves, matching mem_sam_pe's no_pairing block. */
    int n_pri[2], z[2], q_se[2], extra_flag = 0, paired = 0;
    mem_pair_resolve(opt, bns, pac, pes, (uint64_t)pair_idx,
                     s, a, n_pri, z, q_se, &extra_flag, &paired);

    /* Convert each side's regions to mem_aln_t arrays for emission. */
    mem_aln_t *lists[2] = {nullptr, nullptr};
    int n_lists[2] = {0, 0};
    for (int k = 0; k < 2; ++k) {
        if (a[k].n == 0) {
            /* Unmapped read still needs a record. Build a synthetic mem_aln_t. */
            lists[k] = (mem_aln_t *) calloc(1, sizeof(mem_aln_t));
            lists[k][0].rid = -1;
            lists[k][0].pos = -1;
            lists[k][0].flag = 4;
            lists[k][0].n_cigar = 0;
            lists[k][0].cigar = nullptr;
            lists[k][0].mapq = 0;
            lists[k][0].NM = 0;
            lists[k][0].score = -1;
            lists[k][0].sub = -1;
            n_lists[k] = 1;
        } else {
            lists[k] = (mem_aln_t *) calloc(a[k].n, sizeof(mem_aln_t));
            for (int j = 0; j < (int)a[k].n; ++j) {
                lists[k][j] = mem_reg2aln(opt, bns, pac, s[k].l_seq, s[k].seq, &a[k].a[j]);
            }
            n_lists[k] = (int)a[k].n;
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
     * when mem_pair_resolve didn't take the paired branch, check the top
     * region of each side via mem_infer_dir and OR 0x2 into extra_flag if
     * the inferred distance falls within pes[dir].low..high. */
    if (!paired && (opt->flag & MEM_F_PE) && !(opt->flag & MEM_F_NOPAIRING)
        && a[0].n > 0 && a[1].n > 0
        && n_lists[0] > 0 && n_lists[1] > 0
        && lists[0][0].rid >= 0 && lists[0][0].rid == lists[1][0].rid) {
        int64_t dist;
        int dir = mem_infer_dir(bns->l_pac, a[0].a[0].rb, a[1].a[0].rb, &dist);
        if (!pes[dir].failed && dist >= pes[dir].low && dist <= pes[dir].high) {
            extra_flag |= 0x2;
        }
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

    /* Emit records. For each side k, emit all n_lists[k] entries; pair mate
     * is the paired-selected primary of the other side (z[!k] on the paired
     * branch, 0 otherwise — matching bwamem_pair.cpp:545-556's use of
     * &a[!i].a[z[!i]] as the mate anchor). Without this the paired branch
     * would use lists[!k][0] even when mem_pair_resolve picked a non-zero
     * primary, driving RNEXT/PNEXT/TLEN/MC off the wrong mate alignment. */
    for (int k = 0; k < 2; ++k) {
        int mate_idx = 0;
        if (paired && z[!k] >= 0 && z[!k] < n_lists[!k]) mate_idx = z[!k];
        mem_aln_t *mate = (n_lists[!k] > 0) ? &lists[!k][mate_idx] : nullptr;
        for (int j = 0; j < n_lists[k]; ++j) {
            append_bam_record(out, pair_idx, bns, &s[k],
                              &lists[k][j], n_lists[k], lists[k], j, mate);
        }
    }

    /* Free mem_aln_t cigar buffers. */
    for (int k = 0; k < 2; ++k) {
        if (!lists[k]) continue;
        for (int j = 0; j < n_lists[k]; ++j) {
            free(lists[k][j].cigar);
        }
        free(lists[k]);
    }
}

ShimAlignOutput *shim_extend_batch(void *idx_opaque, ShimSeeds *s,
                                   const mem_pestat_t *pestat_in)
{
    if (!s) return nullptr;
    BwaShimIndex *idx = static_cast<BwaShimIndex *>(idx_opaque);
    FMI_search *fmi = idx ? idx->fmi : s->fmi;

    mem_opt_t *opt = s->opts;
    const bntseq_t *bns = fmi->idx->bns;
    const uint8_t  *pac = fmi->idx->pac;

    /* SE extension. */
    run_se_extension(s);

    /* Insert-size. */
    mem_pestat_t pes[4];
    if (pestat_in) memcpy(pes, pestat_in, sizeof(pes));
    else           mem_pestat(opt, bns->l_pac, s->n_seqs, s->w.regs, pes);

    /* Per-pair emission. */
    ShimAlignOutput *out = alloc_align_output((size_t)(s->n_seqs / 2));
    size_t n_pairs = (size_t)(s->n_seqs / 2);
    for (size_t i = 0; i < n_pairs; ++i) {
        mem_alnreg_v ra[2] = { s->w.regs[2*i], s->w.regs[2*i + 1] };
        bseq1_t ss[2]      = { s->seqs[2*i],   s->seqs[2*i + 1]   };
        pair_and_emit(out, i, opt, bns, pac, ss, ra, pes);
        /* Copy back any modified regs (mate-rescue may append). */
        s->w.regs[2*i]     = ra[0];
        s->w.regs[2*i + 1] = ra[1];
    }

    memcpy(out->pes, pes, sizeof(pes));

    /* Free per-read regs. */
    for (int i = 0; i < s->n_seqs; ++i) free(s->w.regs[i].a);
    shim_seeds_free(s);
    return out;
}

/* ------------------ Convenience: align_batch ------------------ */

ShimAlignOutput *shim_align_batch(void *idx_opaque, mem_opt_t *opts,
                                  const ShimReadPair *pairs, size_t n_pairs,
                                  const mem_pestat_t *pestat_in)
{
    ShimSeeds *s = shim_seed_batch(idx_opaque, opts, pairs, n_pairs);
    if (!s) return nullptr;
    return shim_extend_batch(idx_opaque, s, pestat_in);
}

/* ------------------ estimate_pestat ------------------ */

int shim_estimate_pestat(void *idx_opaque, mem_opt_t *opts,
                         const ShimReadPair *pairs, size_t n_pairs,
                         mem_pestat_t *pestat_out)
{
    ShimSeeds *s = shim_seed_batch(idx_opaque, opts, pairs, n_pairs);
    if (!s) return -1;
    FMI_search *fmi = s->fmi;
    run_se_extension(s);
    mem_pestat(opts, fmi->idx->bns->l_pac, s->n_seqs, s->w.regs, pestat_out);
    for (int i = 0; i < s->n_seqs; ++i) free(s->w.regs[i].a);
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
