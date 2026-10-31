//! # iq_quant — U20: transcript model counters, GTF/TSV writers, ID renumbering
//!
//! Design §8, §5.6.
//! - [`counter`] — per-block `AssignedFeatureCounter` (`unique_only`), keyed by `ModelIdx`.
//! - [`ids`] — novel transcript/gene numbering (D-M1, `ExcludingIdDistributor`), exon ids.
//! - [`gtf`] — `GFFPrinter.dump`.
//! - [`output`] — Phase W per chromosome ([`ChrPhaseW`]) and the final files.
//!
//! Production flow (U21): per block, `TxCounter::add` for each U19 `FinalAssign` in read
//! order; per chromosome, `ChrPhaseW::add_block` over the blocks in genomic order, then
//! `finish`; finally `write_outputs` with every chromosome's `ChrOutput`.
//!
//! Fixture runners: unit `U20`, fns `gff_dump`, `counter_add_read_info`, `counter_dump`
//! (see [`fixture`]).

pub mod counter;
pub mod fixture;
pub mod gtf;
pub mod ids;
pub mod output;

pub use counter::{FeatureCounter, ReadFate, Strategy, TxCounter};
pub use gtf::{GtfChrWriter, PrintModel};
pub use ids::{BlockNames, NovelNamer};
pub use output::{fmt_fixed2, write_outputs, BlockOutput, ChrCounts, ChrOutput, ChrPhaseW, OutputSpec};

/// Registers this crate's fixture runners.
pub fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    fixture::register_fixtures(reg);
}
