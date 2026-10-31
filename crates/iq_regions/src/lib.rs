//! # iq_regions — U3: pass A (coverage clusters and blocks) and read ownership
//!
//! Design §5.1-5.2. Entry points:
//! - [`PassA`] / [`pass_a`]: stream a chromosome's records once → [`Cluster`]s of [`Block`]s
//!   (gap bridging by shared genes, coverage bins, valley splits, coverage cap, gene lists,
//!   owned-record counts) + duplicate-name candidates.
//! - [`split_coverage_regions`]: the 4.1 valley splitter, verbatim.
//! - [`owner_block`]: max-overlap owner block (cost-estimate partition; U8 processes a read in every block whose fetch returns it, as 4.1).
//! - [`GeneOracle`]: the gene overlap query; `iq_annot::AnnotationIndex` implements it via a
//!   per-chromosome adapter (`genes_for_region0`, 0-based query as in 4.1).
//!
//! Fixture runners: unit `U3`, fns `split_coverage_regions`, `process_alignments_in_region`,
//! `cluster`, `genes_in_region` ([`register_fixtures`]).

pub mod dups;
pub mod owner;
pub mod pass_a;
pub mod split;

mod fixture;

pub use iq_types::{Block, Cluster};
pub use owner::{overlapping_blocks, owner_block, owner_of};
pub use pass_a::{pass_a, GeneOracle, NoGenes, PassA, PassAConfig, PassARecord, PassAResult, PassAStats};
pub use split::{coverage_cutoff, is_small_chr, skip_read_fraction, split_coverage_regions, CoverageBins};

/// Registers this crate's fixture runners.
pub fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    fixture::register_fixtures(reg);
}
