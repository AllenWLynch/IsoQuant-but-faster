//! # iq_read — U4 (CIGAR blocks, polyA/T, deletions) and U7 (exon correction, strand,
//! antisense, read groups)
//!
//! Fixture runners: unit `U4` (fns `process_alignment`,
//! `extract_relevant_deletions`) and `U7` (fns `correct_assigned_read`,
//! `get_assignment_strand`, `check_antisense`), see [`register_fixtures`].
#![allow(unused_variables)] // stub: remove once implemented

pub mod cigar;
pub mod correct;
pub mod deletions;
pub mod groups;
pub mod polya;

use iq_types::{AlignedRead, Params, RecordLite};

pub use correct::{assignment_strand, check_antisense, correct_read};

/// `fn process_alignment(&RecordLite) -> Option<AlignedRead>` (U4): CIGAR → blocks,
/// polyA/T detection and fixing, boundary deletions. `None` when the
/// alignment has no exons. Allocates only the output vectors; use
/// [`process_alignment_into`] to reuse them across reads.
pub fn process_alignment(rec: &RecordLite<'_>, params: &Params) -> Option<AlignedRead> {
    let mut ar = AlignedRead {
        exons: Vec::new(),
        read_blocks: Vec::new(),
        cigar_blocks: Vec::new(),
        polya: iq_types::NO_POLYA,
        exons_changed: false,
        deletions: Vec::new(),
    };
    process_alignment_into(rec, &polya::PolyAParams::from_params(params), polya::ExternalPolyA::default(), &mut ar)
        .then_some(ar)
}

/// [`process_alignment`] into a reused `AlignedRead` (every field is overwritten), with
/// explicit polyA parameters and external polyA hint (`Default` = `--polya_trimmed none`).
/// Returns `false` (and leaves `out` with empty exons) when the alignment has no exons.
pub fn process_alignment_into(
    rec: &RecordLite<'_>,
    p: &polya::PolyAParams,
    ext: polya::ExternalPolyA,
    out: &mut AlignedRead,
) -> bool {
    cigar::read_blocks_into(rec.pos0, rec.cigar, &mut out.exons, &mut out.read_blocks, &mut out.cigar_blocks);
    out.deletions.clear();
    out.polya = iq_types::NO_POLYA;
    out.exons_changed = false;
    if out.exons.is_empty() {
        return false;
    }
    let (info, changed) =
        polya::add_polya_info(rec, p, ext, &mut out.exons, &mut out.read_blocks, &mut out.cigar_blocks);
    out.polya = info;
    out.exons_changed = changed;
    // Matches IsoQuant 4.1: deletions are anchored at the post-polyT-trim first exon start
    // (`alignment_processor.py:422,499`), not at `reference_start + 1`.
    let anchor = out.exons[0].s;
    deletions::extract_relevant_deletions_into(
        rec.cigar,
        anchor,
        &out.exons,
        deletions::BOUNDARY_MARGIN,
        &mut out.deletions,
    );
    true
}

/// Registers this crate's fixture runners. Each module registers its own runners.
/// U7: add `correct::register_fixtures(reg);` here.
pub fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    polya::register_fixtures(reg);
    deletions::register_fixtures(reg);
    correct::register_fixtures(reg);
}
