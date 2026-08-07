/* bwa-mem3-sys/shim/bwa_shim_layout_assert.cpp
 *
 * Compile-time guard that the POD copies of `mem_opt_t` / `mem_pestat_t` in
 * `bwa_shim_types.h` stay byte-identical to upstream's real definitions in
 * `bwamem.h`. This TU exports no symbols — it exists only for its
 * `static_assert`s, so any field-order/size drift on a future vendor refresh
 * fails `cargo build` loudly instead of silently corrupting the offsets Rust
 * reads/writes through the bindgen view of the POD (the shim allocates the
 * real struct via `mem_opt_init()`; see the maintenance docs gotcha #2).
 *
 * The bindgen `bindgen_test_layout_*` tests only check the generated Rust
 * struct against the C compiler's view of the POD — both derived from
 * `bwa_shim_types.h` — so they cannot catch POD-vs-upstream drift. This does.
 *
 * The POD and upstream both spell the struct `mem_opt_t`, which C++ forbids
 * defining twice with different bodies in one TU, so the POD is pulled in
 * first under renamed tags before the real `bwamem.h`.
 */
#include <cstddef>
#include <cstdint>

#define mem_opt_t    pod_mem_opt_t
#define mem_pestat_t pod_mem_pestat_t
#include "bwa_shim_types.h"
#undef mem_opt_t
#undef mem_pestat_t

/* ---------------------------------------------------------------------------
 * MEM_F_* bit mirror — capture step.
 *
 * `bwa_shim_types.h` hardcodes the flag bits and bindgen's
 * `allowlist_var("MEM_F_.*")` reads that POD copy, never upstream — so a
 * renumbered flag would silently corrupt `opt->flag` with no signal on either
 * side (maintenance docs gotcha #2, extended).
 *
 * The struct-tag rename trick above does NOT work for macros: they have one
 * namespace and no rename mechanism, so `#define MEM_F_PE POD_MEM_F_PE`
 * before including the POD merely redefines MEM_F_PE and leaves one value,
 * making any comparison a tautology.
 *
 * This capture MUST happen here, before `#include "bwamem.h"` below. That
 * include defines its own MEM_F_* macros; since bwa_shim_types.h already
 * defined them above, the preprocessor's "last definition wins" redefinition
 * rule means bwamem.h's value silently replaces the POD's from this point
 * forward (with a -Wmacro-redefined warning if they differ, but otherwise no
 * signal). Capturing *after* that include — which is where the struct
 * static_asserts below sit — would therefore compare upstream's value
 * against itself, not against the POD's: a tautology, verified by compiling
 * that exact ordering with a deliberately wrong MEM_F_PE and observing no
 * error. Capturing here, before the clobber, is the only ordering that
 * actually reads the POD's value.
 *
 * Note this cannot detect a NEW upstream flag — there is no assert for a
 * name you do not know about. `scripts/bwa-mem3-drift-report.sh` compares
 * the flag *set* on a vendor refresh to cover that direction.
 *
 * Do not add -Wno-macro-redefined to the warning-silencing list in build.rs:
 * until this guard existed, that warning was the only drift signal at all.
 * ------------------------------------------------------------------------- */
namespace pod_flags {
constexpr int MEM_F_PE_v = MEM_F_PE;
constexpr int MEM_F_NOPAIRING_v = MEM_F_NOPAIRING;
constexpr int MEM_F_ALL_v = MEM_F_ALL;
constexpr int MEM_F_NO_MULTI_v = MEM_F_NO_MULTI;
constexpr int MEM_F_NO_RESCUE_v = MEM_F_NO_RESCUE;
constexpr int MEM_F_REF_HDR_v = MEM_F_REF_HDR;
constexpr int MEM_F_SOFTCLIP_v = MEM_F_SOFTCLIP;
constexpr int MEM_F_SMARTPE_v = MEM_F_SMARTPE;
constexpr int MEM_F_PRIMARY5_v = MEM_F_PRIMARY5;
constexpr int MEM_F_KEEP_SUPP_MAPQ_v = MEM_F_KEEP_SUPP_MAPQ;
constexpr int MEM_F_XB_v = MEM_F_XB;
} // namespace pod_flags

#undef MEM_F_PE
#undef MEM_F_NOPAIRING
#undef MEM_F_ALL
#undef MEM_F_NO_MULTI
#undef MEM_F_NO_RESCUE
#undef MEM_F_REF_HDR
#undef MEM_F_SOFTCLIP
#undef MEM_F_SMARTPE
#undef MEM_F_PRIMARY5
#undef MEM_F_KEEP_SUPP_MAPQ
#undef MEM_F_XB

#include "bwamem.h"

#define BWA_SHIM_CK(real, pod, field)                            \
    static_assert(offsetof(real, field) == offsetof(pod, field), \
                  #real " offset drift: " #field)
#define BWA_SHIM_CK_OPT(f) BWA_SHIM_CK(mem_opt_t, pod_mem_opt_t, f)
#define BWA_SHIM_CK_PES(f) BWA_SHIM_CK(mem_pestat_t, pod_mem_pestat_t, f)

BWA_SHIM_CK_OPT(a);
BWA_SHIM_CK_OPT(b);
BWA_SHIM_CK_OPT(o_del);
BWA_SHIM_CK_OPT(e_del);
BWA_SHIM_CK_OPT(o_ins);
BWA_SHIM_CK_OPT(e_ins);
BWA_SHIM_CK_OPT(pen_unpaired);
BWA_SHIM_CK_OPT(pen_clip5);
BWA_SHIM_CK_OPT(pen_clip3);
BWA_SHIM_CK_OPT(w);
BWA_SHIM_CK_OPT(zdrop);
BWA_SHIM_CK_OPT(max_mem_intv);
BWA_SHIM_CK_OPT(T);
BWA_SHIM_CK_OPT(flag);
BWA_SHIM_CK_OPT(min_seed_len);
BWA_SHIM_CK_OPT(min_ext_len);
BWA_SHIM_CK_OPT(max_extend_chains);
BWA_SHIM_CK_OPT(mate_concordant_window);
BWA_SHIM_CK_OPT(est_insert_high);
BWA_SHIM_CK_OPT(seed_emit_order);
BWA_SHIM_CK_OPT(min_chain_weight);
BWA_SHIM_CK_OPT(max_chain_extend);
BWA_SHIM_CK_OPT(split_factor);
BWA_SHIM_CK_OPT(split_width);
BWA_SHIM_CK_OPT(max_occ);
BWA_SHIM_CK_OPT(max_chain_gap);
BWA_SHIM_CK_OPT(n_threads);
BWA_SHIM_CK_OPT(chunk_size);
BWA_SHIM_CK_OPT(mask_level);
BWA_SHIM_CK_OPT(drop_ratio);
BWA_SHIM_CK_OPT(XA_drop_ratio);
BWA_SHIM_CK_OPT(mask_level_redun);
BWA_SHIM_CK_OPT(mapQ_coef_len);
BWA_SHIM_CK_OPT(mapQ_coef_fac);
BWA_SHIM_CK_OPT(max_ins);
BWA_SHIM_CK_OPT(max_matesw);
BWA_SHIM_CK_OPT(rescue_kmer);
BWA_SHIM_CK_OPT(rescue_band);
BWA_SHIM_CK_OPT(rescue_skip);
BWA_SHIM_CK_OPT(max_XA_hits);
BWA_SHIM_CK_OPT(max_XA_hits_alt);
BWA_SHIM_CK_OPT(mat);
BWA_SHIM_CK_OPT(mat_ot);
BWA_SHIM_CK_OPT(mat_ob);
BWA_SHIM_CK_OPT(bam_mode);
BWA_SHIM_CK_OPT(bam_level);
BWA_SHIM_CK_OPT(meth_mode);
BWA_SHIM_CK_OPT(meth_scoring);
BWA_SHIM_CK_OPT(meth_chem);
BWA_SHIM_CK_OPT(meth_tags);
BWA_SHIM_CK_OPT(meth_set_as_failed);
BWA_SHIM_CK_OPT(meth_chimera_qc);
BWA_SHIM_CK_OPT(proper_pair_from_emitted);
BWA_SHIM_CK_OPT(supp_rep_hard_cap);
BWA_SHIM_CK_OPT(smem_dedup);
BWA_SHIM_CK_OPT(alnreg_sort_fast);
BWA_SHIM_CK_OPT(skip_contained_ext);
BWA_SHIM_CK_OPT(band_start);
BWA_SHIM_CK_OPT(compat);
static_assert(sizeof(mem_opt_t) == sizeof(pod_mem_opt_t), "mem_opt_t size drift");

BWA_SHIM_CK_PES(low);
BWA_SHIM_CK_PES(high);
BWA_SHIM_CK_PES(failed);
BWA_SHIM_CK_PES(avg);
BWA_SHIM_CK_PES(std);
static_assert(sizeof(mem_pestat_t) == sizeof(pod_mem_pestat_t), "mem_pestat_t size drift");

/* MEM_F_* bit mirror — comparison step. `bwamem.h` above has already
 * (re)defined every MEM_F_* macro to its real value (the #undef block
 * before that include made the redefinition clean rather than clobbering),
 * so `f` below reads upstream's value while `pod_flags::f##_v` holds what
 * was captured from the POD before the clobber. See the capture-step
 * comment above `#include "bwamem.h"` for why the capture has to happen
 * there and not here. */
#define BWA_SHIM_CK_FLAG(f)                            \
    static_assert(pod_flags::f##_v == f,               \
                  "MEM_F_* drift: " #f " differs from upstream bwamem.h")

BWA_SHIM_CK_FLAG(MEM_F_PE);
BWA_SHIM_CK_FLAG(MEM_F_NOPAIRING);
BWA_SHIM_CK_FLAG(MEM_F_ALL);
BWA_SHIM_CK_FLAG(MEM_F_NO_MULTI);
BWA_SHIM_CK_FLAG(MEM_F_NO_RESCUE);
BWA_SHIM_CK_FLAG(MEM_F_REF_HDR);
BWA_SHIM_CK_FLAG(MEM_F_SOFTCLIP);
BWA_SHIM_CK_FLAG(MEM_F_SMARTPE);
BWA_SHIM_CK_FLAG(MEM_F_PRIMARY5);
BWA_SHIM_CK_FLAG(MEM_F_KEEP_SUPP_MAPQ);
BWA_SHIM_CK_FLAG(MEM_F_XB);
