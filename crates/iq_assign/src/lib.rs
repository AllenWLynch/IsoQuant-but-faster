//! # iq_assign — U5 (read profiles) and U6 (read→isoform assignment)
//!
//! One assigner, used against reference and model indices.
//! Fixture runners: unit `U5` (fns `construct_profiles`, `match_genomic_features`) and
//! `U6` (fn `assign_to_isoform`), see [`register_fixtures`].
#![allow(unused_variables)] // stub: remove once implemented

pub mod assigner;
pub mod inconsistent;
pub mod junction;
pub mod polya_verify;
pub mod profile;

pub use assigner::{assign, AssignParams};
pub use profile::{build_profiles, match_genomic_features, Scratch};

/// Registers this crate's fixture runners.
/// Each module registers its own runners, so U5 and U6 owners never edit this file.
pub fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    profile::register_fixtures(reg);
    assigner::register_fixtures(reg);
}
