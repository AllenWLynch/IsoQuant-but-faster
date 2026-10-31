//! The stage interface the scheduler drives (design §5.4).
//!
//! The scheduler ([`crate::engine::execute`]) knows nothing about reads, models or files; it
//! drives an implementation of [`Pipeline`] through the two-phase run:
//!
//! ```text
//! Phase A   pass_a(chr)                                 -> blocks of chr      (one task per chr)
//! Phase C   collect_block(chr, b, block)                -> Store              (per work unit, blocks in order)
//! Barrier   barrier_chr(chr, blocks, &mut stores)       -> ReadTotals         (per chr, after all its units)
//!           decide(Σ ReadTotals)                        -> Global             (once, after every chr)
//! Phase M   construct_block(chr, b, block, store, &mut UnitState, &Global) -> BlockOut
//!                                                                             (per unit, blocks in order)
//! Phase W   merge_chr(chr, blocks, outs in block order) -> ChrOut             (per chr, when its units are done)
//!           finish(ChrOut of every chr, in chr index order)                   (once)
//! ```
//!
//! Every method may be called concurrently from rayon workers (for *different* units /
//! chromosomes) and may itself use rayon (`par_iter`, `join`): nested parallelism inside a
//! giant block uses idle workers through work stealing.
//!
//! **Determinism contract.** Outputs must depend only on the method arguments (plus
//! read-only pipeline state), never on call order across units or chromosomes. The scheduler
//! guarantees: blocks of one unit are collected and constructed sequentially in genomic order
//! with one `UnitState` threaded through construct; `merge_chr` receives the block outputs of
//! a chromosome in genomic order; `finish` receives chromosomes in index order.

use iq_types::Block;

/// Sample-wide read totals used by the barrier decision (4.1 `total_assignments`,
/// `polya_assignments` after multimapper resolution, `processed_read_manager.py:60-87`).
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct ReadTotals {
    pub total: u64,
    pub polya: u64,
}

impl std::ops::AddAssign for ReadTotals {
    fn add_assign(&mut self, o: Self) {
        self.total += o.total;
        self.polya += o.polya;
    }
}

/// Default per-read byte constants of the admission estimates (design §5.4; calibrated in
/// `docs/vtg_proto.md`: ~99 B/read with a stub assignment; U21 real runs measured 187-249 B/read
/// resident `len_bytes` with real assignments, PacBio and ONT).
pub const EST_COLLECT_BYTES_PER_READ: u64 = 224;
pub const EST_CONSTRUCT_BYTES_PER_READ: u64 = 48;

pub trait Pipeline: Sync {
    /// Collected per-block data (a `BlockStore` plus whatever the barrier needs).
    type Store: Send;
    /// State threaded through the blocks of one unit in Phase M (`detected_known`).
    type UnitState: Default + Send;
    /// Sample-wide decisions taken at the barrier (`requires_polya_for_construction`).
    type Global: Send + Sync;
    /// Per-block construct + quant output, held until its chromosome merges.
    type BlockOut: Send;
    /// Per-chromosome merged output (Phase W).
    type ChrOut: Send;

    fn n_chromosomes(&self) -> usize;
    fn chr_name(&self, chr: usize) -> &str;
    /// Pass-A cost estimate (≈ records), used to order pass-A tasks.
    fn chr_weight(&self, chr: usize) -> u64;

    /// Phase A: blocks of `chr` in genomic order (all clusters concatenated).
    fn pass_a(&self, chr: usize) -> Result<Vec<Block>, String>;
    /// Estimated transient bytes of a pass-A task.
    fn est_pass_a(&self, chr: usize) -> u64 {
        let _ = chr;
        128 << 20
    }

    /// Phase C: collect one block.
    fn collect_block(&self, chr: usize, block_idx: usize, block: &Block) -> Result<Self::Store, String>;
    /// Measured resident bytes of a collected store (length-based; see vtg_proto finding #9).
    fn store_bytes(&self, store: &Self::Store) -> u64;
    /// Whether [`Pipeline::spill`] is implemented.
    fn can_spill(&self) -> bool {
        false
    }
    /// Move a store out of anonymous memory (e.g. write it to a file and mmap it). Returns the
    /// new resident bytes. Called only when [`Pipeline::can_spill`].
    fn spill(&self, chr: usize, block_idx: usize, store: &mut Self::Store) -> Result<u64, String> {
        let _ = (chr, block_idx);
        Ok(self.store_bytes(store))
    }

    /// Barrier, per chromosome: duplicate read names (U9, chromosome-local in 4.1's default
    /// `NoSecondary` manager) are resolved across `stores`; returns the chromosome's totals.
    fn barrier_chr(&self, chr: usize, blocks: &[Block], stores: &mut [Self::Store]) -> Result<ReadTotals, String>;
    /// Barrier, global: the sample-wide decision from the summed totals.
    fn decide(&self, totals: ReadTotals) -> Self::Global;

    /// Phase M: construct and quantify one block; consumes its store.
    fn construct_block(
        &self,
        chr: usize,
        block_idx: usize,
        block: &Block,
        store: Self::Store,
        state: &mut Self::UnitState,
        global: &Self::Global,
    ) -> Result<Self::BlockOut, String>;

    /// Phase W: merge the block outputs of one chromosome (genomic order), renumbering IDs.
    fn merge_chr(&self, chr: usize, blocks: &[Block], outs: Vec<Self::BlockOut>) -> Result<Self::ChrOut, String>;
    /// Writes the final files. `outs[i]` belongs to chromosome `i`.
    fn finish(&self, outs: Vec<Self::ChrOut>, totals: ReadTotals) -> Result<(), String>;

    /// Bytes outside the scheduler's accounting (annotation, per-thread scratch), subtracted
    /// from the budget.
    fn fixed_bytes(&self) -> u64 {
        0
    }
    /// Estimated transient bytes of collecting `block` (store under construction + batches).
    fn est_collect(&self, block: &Block) -> u64 {
        block.owned_reads * EST_COLLECT_BYTES_PER_READ
    }
    /// Estimated transient bytes of constructing `block` (store already counted as resident).
    fn est_construct(&self, block: &Block) -> u64 {
        block.owned_reads * EST_CONSTRUCT_BYTES_PER_READ
    }
    /// Estimated CPU cost (LPT priority), ≈ reads.
    fn est_cost(&self, block: &Block) -> u64 {
        block.owned_reads.max(1)
    }
}
