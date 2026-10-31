//! # iq_annot — U1: annotation parsing and `IsoformIndex` construction
//!
//! Design §4: native GFF3/GTF parsing into a flat [`AnnotationIndex`] (§2.2, §4.1-4.2),
//! gene overlap queries (§4.3, D-R1 in [`AnnotationIndex::genes_for_region0`]), and the
//! `IsoformIndex` builders (`GeneInfo.__init__` / `from_models` / `from_region`).
//!
//! Fixture runners (unit `U1`, fns `dump_annotation`, `get_gene_info_for_region`,
//! `from_models`) live in [`fixture`].

pub mod attrs;
pub mod builder;
pub mod fixture;
pub mod index;
pub mod parse;

pub use builder::{index_from_genes, index_from_models, index_from_models_named, index_from_region, split_exons};
pub use index::{AnnotationIndex, GeneTypes};
pub use parse::{parse_gff3, parse_gff3_reader, parse_gtf, parse_gtf_reader, AnnotError};

/// Registers this crate's fixture runners.
pub fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    fixture::register_fixtures(reg);
}
