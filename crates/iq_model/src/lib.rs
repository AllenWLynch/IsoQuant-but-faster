//! # iq_model — U13-U19: model construction, filtering, end processing, splice-site
//! correction, reassignment, gene joining, final per-read assignment
//!
//! Stages mutate a [`iq_types::ModelStore`] given a read-only
//! [`Ctx`]. Nested U6 calls go through the [`ReadAssigner`] trait so stages can be tested
//! with a replaying assigner before U6 passes.
//! Fixture runners: units `U13`-`U19` (see SCHEMA §4), [`register_fixtures`].
#![allow(unused_variables)] // stub: remove once implemented

pub mod assign_based;
pub mod ctx;
pub mod ends;
pub mod filter;
pub mod final_assign;
pub mod fixture_util;
pub mod fl;
pub mod joiner;
pub mod reassign;
pub mod splice_site;

pub use ctx::{Ctx, DetectedKnown, ReadAssigner};

/// Registers this crate's fixture runners.
pub fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    assign_based::register_fixtures(reg);
    fl::register_fixtures(reg);
    filter::register_fixtures(reg);
    ends::register_fixtures(reg);
    splice_site::register_fixtures(reg);
    reassign::register_fixtures(reg);
    final_assign::register_fixtures(reg);
    joiner::register_fixtures(reg);
}
