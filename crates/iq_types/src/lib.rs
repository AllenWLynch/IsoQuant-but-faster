//! # iq_types — the iqrs contract crate (unit U0, milestone M0)
//!
//! Every other crate depends on this one, and units U1–U21 are developed against it in
//! parallel. Contents:
//!
//! | module | contents |
//! |---|---|
//! | [`prim`] | `Pos`, `Iv`, `Iv0`, `Strand`, `PolyAInfo`, index aliases |
//! | [`enums`] | IsoQuant enums with their numeric codes, cost table, event sets |
//! | [`assign`] | `MatchEvent`, `IsoformMatch`, `Assignment`, `GeneKey`, `IsoRef` |
//! | [`read`] | `ReadRecord`, `RecordLite`, `AlignedRead`, `RefSlice`, `Block` |
//! | [`store`] | `BlockStore`, `ChainTable`, `ReadView` |
//! | [`index`] | `IsoformIndex`, `ProfMatrix`, read profiles |
//! | [`model`] | `Vertex`, `IntronGraph`, `PathStorage`, `TranscriptModel`, `ModelStore` |
//! | [`intern`] | `Interner`, `GroupPool` |
//! | [`params`] | `Params` and the two built-in presets |
//! | [`fixture`] | fixture IO, schema DTOs, `UnitFixture`, `Verdict`, diff, registry |
//!
//! Conventions: genomic coordinates are 1-based closed (`Iv`); maps whose order can reach an
//! output are `IndexMap` (= Python dict insertion order); enums carry IsoQuant codes.
//!
//! ## For unit owners: adding a fixture runner
//! ```ignore
//! use iq_types::fixture::{schema::*, FixResult, FixtureCtx, UnitFixture, TxSpace};
//! pub struct U6Assign;
//! impl UnitFixture for U6Assign {
//!     const UNIT: &'static str = "U6";
//!     const FN: &'static str = "assign_to_isoform";          // = in.fn
//!     type In = U6AssignIn;
//!     type Out = U6AssignOut;
//!     fn run(i: &U6AssignIn, ctx: &FixtureCtx) -> FixResult<U6AssignOut> {
//!         let index = ctx.index(&i.index)?;                    // IsoformIndex from index/<sha1>.json
//!         let a = /* crate entry fn */;
//!         Ok(U6AssignOut { assignment: a.map(|a| AssignmentFixture::from_core(&a, &ctx.ids())) })
//!     }
//! }
//! // in your crate:  pub fn register_fixtures(reg: &mut Registry) { reg.add::<U6Assign>(); }
//! ```
//! Then `cargo run --release -p iqrs-fixture -- run U6 fixtures/`. Compare in fixture (string)
//! space; override `UnitFixture::compare` for set semantics or tolerances (`DiffOptions`).

pub mod assign;
pub mod enums;
pub mod fixture;
pub mod index;
pub mod intern;
pub mod model;
pub mod params;
pub mod prim;
pub mod read;
pub mod store;

pub use assign::{derive_gene_atype, Assignment, Events, GeneKey, IsoRef, IsoformMatch, MatchEvent, Matches};
pub use enums::{
    elongation_cost, match_consts, AmbiguityResolving, AssignmentType, EventType, MatchClass, MultimapStrategy,
    PolyATrimmed, PolyAUsage, StrandnessReporting, TranscriptModelType, ABSENT_POSITION,
};
pub use index::{CombinedProfile, GeneRec, IsoRec, IsoformIndex, MappedReadProfile, ProfMatrix, SortedSpans};
pub use intern::{GroupPool, Interner};
pub use model::{
    vertex_code, FinalAssign, IntronGraph, ModelKind, ModelStore, Path, PathStorage, Peak, Predictions,
    TranscriptModel, VId, Vertex,
};
pub use params::{DataType, Params};
pub use prim::{
    junctions_from_blocks, ChrIdx, GeneIdx, Iv, Iv0, ModelIdx, PolyAInfo, Pos, ReadIdx, Strand, TxIdx, NO_POLYA,
};
pub use read::{AlignedRead, Block, CigarOp, Cluster, ReadRecord, RecordLite, RecordOwned, RefSlice, TagValue};
pub use store::{name_hash, BlockStore, ChainTable, MatchRef, ReadView};

/// Re-exports for downstream crates so they share exact versions.
pub use indexmap;
pub use serde;
pub use serde_json;
pub use smallvec;
