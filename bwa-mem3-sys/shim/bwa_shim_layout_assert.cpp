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
BWA_SHIM_CK_OPT(max_XA_hits);
BWA_SHIM_CK_OPT(max_XA_hits_alt);
BWA_SHIM_CK_OPT(mat);
BWA_SHIM_CK_OPT(mat_ot);
BWA_SHIM_CK_OPT(mat_ob);
BWA_SHIM_CK_OPT(bam_mode);
BWA_SHIM_CK_OPT(bam_level);
BWA_SHIM_CK_OPT(meth_mode);
BWA_SHIM_CK_OPT(meth_scoring);
BWA_SHIM_CK_OPT(meth_set_as_failed);
BWA_SHIM_CK_OPT(meth_chimera_qc);
BWA_SHIM_CK_OPT(supp_rep_hard_cap);
BWA_SHIM_CK_OPT(smem_dedup);
BWA_SHIM_CK_OPT(skip_contained_ext);
BWA_SHIM_CK_OPT(band_start);
static_assert(sizeof(mem_opt_t) == sizeof(pod_mem_opt_t), "mem_opt_t size drift");

BWA_SHIM_CK_PES(low);
BWA_SHIM_CK_PES(high);
BWA_SHIM_CK_PES(failed);
BWA_SHIM_CK_PES(avg);
BWA_SHIM_CK_PES(std);
static_assert(sizeof(mem_pestat_t) == sizeof(pod_mem_pestat_t), "mem_pestat_t size drift");
