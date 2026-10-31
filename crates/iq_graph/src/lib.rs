//! # iq_graph — U11 (intron graph) and U12 (path threading)
//!
//! Design §7.2-7.3. Results are the contract types [`iq_types::IntronGraph`] and
//! [`iq_types::PathStorage`]; graph queries other units need (component coverage, ...)
//! live here as functions over `IntronGraph`.
//! Fixture runners: unit `U11` (fn `IntronGraph.__init__`) and `U12` (fns
//! `IntronPathStorage.fill`, `get_known_spliced_isoforms`), see [`register_fixtures`].
//!
//! U11 modules: [`collector`] (`IntronCollector`), [`graph`] (`IntronGraph.__init__`),
//! [`query`] (component coverage). U12 lives in [`paths`].

pub mod collector;
pub mod graph;
pub mod paths;
mod pyds;
pub mod query;

pub use collector::IntronCollector;
pub use graph::{build_graph, build_graph_with, GraphParams};
pub use paths::thread_paths;
pub use query::{connected_component, max_component_coverage, overlapping_component_max_coverage};

/// Registers this crate's fixture runners.
pub fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    graph::register_fixtures(reg);
    paths::register_fixtures(reg);
}
