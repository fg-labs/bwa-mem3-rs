/*
 * `explicit_mask` bits for bwa_shim_opts_apply_meth_defaults_masked.
 *
 * Plain C macros, so they can be included both by the public header
 * (bwa_shim.h, for bindgen) and by bwa_shim_align.cpp, which cannot include
 * bwa_shim.h (its POD mem_opt_t collides with upstream's real one). One
 * definition means the two translation units cannot drift apart.
 *
 * Each bit says the caller set that knob itself, so the `--meth` defaults must
 * leave it alone. This is upstream's `opt0` "the user passed this flag" mask
 * (fastmap.cpp, main_mem) for the two knobs whose `--meth` default depends on
 * it.
 */
#ifndef BWA_SHIM_METH_H
#define BWA_SHIM_METH_H

#define BWA_SHIM_METH_SET_SCORING    0x1u /* meth_scoring (--meth-scoring)       */
#define BWA_SHIM_METH_SET_SEED_PRUNE 0x2u /* meth_seed_prune (--meth-seed-prune) */

#endif /* BWA_SHIM_METH_H */
