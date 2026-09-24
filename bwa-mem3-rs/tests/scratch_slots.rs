//! Kernel thread slots of live scratches. In its own test binary: the slots
//! are process-global, so scratches other tests create concurrently would share
//! the pool this test inspects.

use bwa_mem3_rs::AlignScratch;

/// Scratches that are alive together get distinct kernel thread slots, and 32
/// of them land on distinct 64-byte lines of bwa-mem3's per-thread profiling
/// counters (8 slots per line), so concurrent workers never write the same
/// counter line on 64-byte-line hardware; freed slots are handed out again.
#[test]
fn live_scratches_get_distinct_slots_on_distinct_lines_and_freed_slots_return() {
    let scratches: Vec<AlignScratch> = (0..32).map(|_| AlignScratch::new().unwrap()).collect();
    let slots: Vec<usize> = scratches.iter().map(AlignScratch::tid_slot).collect();
    let distinct = |v: &[usize]| v.iter().collect::<std::collections::BTreeSet<_>>().len();
    assert_eq!(distinct(&slots), 32, "slots: {slots:?}");
    let lines: Vec<usize> = slots.iter().map(|s| s / 8).collect();
    assert_eq!(distinct(&lines), 32, "lines: {lines:?}");
    assert!(
        slots.iter().all(|&s| s < 256),
        "slots stay inside LIM_C/MAX_THREADS: {slots:?}"
    );

    // A freed scratch's slot is reused, so a long-running process that creates
    // and drops scratches never runs out of distinct slots. (One test, run in
    // sequence: the slot pool is process-global.)
    drop(scratches);
    let reused: Vec<AlignScratch> = (0..32).map(|_| AlignScratch::new().unwrap()).collect();
    let mut again: Vec<usize> = reused.iter().map(AlignScratch::tid_slot).collect();
    let mut first = slots.clone();
    first.sort_unstable();
    again.sort_unstable();
    assert_eq!(again, first, "the freed slots are handed out again");
}
