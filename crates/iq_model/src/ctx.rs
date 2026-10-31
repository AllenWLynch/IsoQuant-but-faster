//! Read-only inputs shared by every construction stage (4.1 `ModelConstructionContext`,
//! `model_construction/context.py:113-224`).

use iq_types::indexmap::IndexSet;
use iq_types::{Assignment, BlockStore, CombinedProfile, IntronGraph, IsoformIndex, Params, PathStorage, Predictions, TxIdx};

/// `detected_known_isoforms`, scoped per work unit (D-M2), insertion-ordered.
pub type DetectedKnown = IndexSet<TxIdx>;

/// Abstraction over U6 so stages can run against recorded assignments
/// (`ReplayAssigner`, design §3 U13) or the live assigner.
pub trait ReadAssigner: Sync {
    /// Assign a profile against `index` (`quick_mode` as in 4.1 `assign_to_isoform`).
    /// `read_key` identifies the call for replay (read index or provisional model key).
    fn assign(&self, index: &IsoformIndex, profile: &CombinedProfile, quick_mode: bool, read_key: u64) -> Option<Assignment>;
}

/// Per-block construction context.
pub struct Ctx<'a> {
    pub params: &'a Params,
    pub index: &'a IsoformIndex,
    pub store: &'a BlockStore,
    pub graph: &'a IntronGraph,
    pub paths: &'a PathStorage,
    pub predictions: &'a Predictions,
    pub assigner: &'a dyn ReadAssigner,
}
