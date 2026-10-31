//! Assignment records: `MatchEvent`, `IsoformMatch`, `Assignment`, and
//! the id-free references `GeneKey` / `IsoRef`.
//!
//! These are the **core** (in-memory) types. They hold numeric keys, not IsoQuant id
//! strings; fixture JSON carries strings and is converted through `fixture::IdMaps`.

use crate::enums::{match_consts, AssignmentType, EventType, MatchClass};
use crate::prim::{GeneIdx, ModelIdx, Strand, TxIdx};
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;

/// IsoQuant `MatchEvent` (IA:422-430).
///
/// Regions are `i64` because the sentinel `undefined_position = 1 << 31` does not fit `i32`.
/// Defaults (`MatchEvent::new(ty)`) are those of the Python constructor:
/// `undefined_region` for both regions and `event_info = 0`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct MatchEvent {
    pub ty: EventType,
    pub iso_region: (i64, i64),
    pub read_region: (i64, i64),
    pub info: i64,
}

impl MatchEvent {
    /// `MatchEvent(event_type)` with Python defaults.
    pub const fn new(ty: EventType) -> Self {
        MatchEvent {
            ty,
            iso_region: match_consts::UNDEFINED_REGION,
            read_region: match_consts::UNDEFINED_REGION,
            info: 0,
        }
    }
    pub const fn with_regions(ty: EventType, iso_region: (i64, i64), read_region: (i64, i64)) -> Self {
        MatchEvent { ty, iso_region, read_region, info: 0 }
    }
    pub const fn with_info(ty: EventType, info: i64) -> Self {
        MatchEvent { info, ..MatchEvent::new(ty) }
    }
}

/// Gene reference. `Ref` = annotated gene (index into the annotation, or into the fixture
/// gene interner). `Novel(k)` = block-local novel gene number `k`, renumbered to
/// `novel_gene_<chr>_<N>` at output (Phase W).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum GeneKey {
    Ref(GeneIdx),
    Novel(u32),
}

/// Isoform reference: an annotated transcript or a block-local model (creation order).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum IsoRef {
    Ref(TxIdx),
    Model(ModelIdx),
}

/// Inline capacity of the event list of one match (1-2 events per read is typical).
pub type Events = SmallVec<[MatchEvent; 2]>;
/// Inline capacity of the match list of one assignment.
pub type Matches = SmallVec<[IsoformMatch; 1]>;

/// IsoQuant `IsoformMatch` (IA:455-538).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct IsoformMatch {
    /// `assigned_gene` (`None` in Python → `None`).
    pub gene: Option<GeneKey>,
    /// `assigned_transcript`.
    pub tx: Option<IsoRef>,
    /// `transcript_strand` (default `'.'`).
    pub tx_strand: Strand,
    /// `match_classification`.
    pub class: MatchClass,
    /// `penalty_score` (f64, summed in event-list order).
    pub penalty: f64,
    /// `match_subclassifications`, in order.
    pub events: Events,
}

impl IsoformMatch {
    /// Python constructor with a list of subclassifications: events of type `none` are
    /// filtered out (IA:477-483).
    pub fn new(
        class: MatchClass,
        gene: Option<GeneKey>,
        tx: Option<IsoRef>,
        events: impl IntoIterator<Item = MatchEvent>,
        tx_strand: Strand,
        penalty: f64,
    ) -> Self {
        IsoformMatch {
            gene,
            tx,
            tx_strand,
            class,
            penalty,
            events: events.into_iter().filter(|e| e.ty != EventType::None).collect(),
        }
    }

    /// Python constructor with a single subclassification object (not a list): it is kept
    /// even if its type is `none` (IA:484-485).
    pub fn new_single(
        class: MatchClass,
        gene: Option<GeneKey>,
        tx: Option<IsoRef>,
        event: MatchEvent,
        tx_strand: Strand,
        penalty: f64,
    ) -> Self {
        let mut events = Events::new();
        events.push(event);
        IsoformMatch { gene, tx, tx_strand, class, penalty, events }
    }

    /// `add_subclassification` (IA:526-531): replaces a lone `undefined`/`none` event,
    /// otherwise appends.
    pub fn add_event(&mut self, ev: MatchEvent) {
        if self.events.len() == 1 && matches!(self.events[0].ty, EventType::Undefined | EventType::None) {
            self.events.clear();
        }
        self.events.push(ev);
    }

    /// `monoexon_is_consistent` (IA:536-538).
    pub fn monoexon_is_consistent(&self) -> bool {
        self.events.iter().all(|e| matches!(e.ty, EventType::None | EventType::MonoExonic | EventType::MonoExonMatch))
    }
}

/// The assignment part of IsoQuant `ReadAssignment`: `assignment_type`,
/// `gene_assignment_type`, `isoform_matches`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Assignment {
    pub atype: AssignmentType,
    pub gene_atype: AssignmentType,
    pub matches: Matches,
}

impl Assignment {
    /// Builds an assignment the way `ReadAssignment.__init__` does, deriving
    /// `gene_assignment_type` from the type and the matches' genes (IA:841-861).
    pub fn new(atype: AssignmentType, matches: Matches) -> Self {
        let gene_atype = derive_gene_atype(atype, &matches);
        Assignment { atype, gene_atype, matches }
    }

    /// An assignment with no matches (`ReadAssignment(read_id, atype)`).
    pub fn empty(atype: AssignmentType) -> Self {
        Assignment::new(atype, Matches::new())
    }
}

/// `gene_assignment_type` derivation of `ReadAssignment.__init__` (IA:841-861).
pub fn derive_gene_atype(atype: AssignmentType, matches: &[IsoformMatch]) -> AssignmentType {
    use AssignmentType as A;
    // number of distinct values of `m.assigned_gene`, counting `None` as a value
    let n_distinct = |include_none: bool| {
        let mut seen: SmallVec<[Option<GeneKey>; 4]> = SmallVec::new();
        for m in matches {
            if (include_none || m.gene.is_some()) && !seen.contains(&m.gene) {
                seen.push(m.gene);
            }
        }
        seen.len()
    };
    match atype {
        A::Ambiguous => {
            if n_distinct(true) > 1 {
                A::Ambiguous
            } else {
                A::Unique
            }
        }
        A::InconsistentAmbiguous => {
            if n_distinct(true) > 1 {
                A::InconsistentAmbiguous
            } else {
                A::Inconsistent
            }
        }
        A::Noninformative | A::Intergenic => match n_distinct(false) {
            0 => atype,
            1 => A::InconsistentGenic,
            _ => A::InconsistentMultigenic,
        },
        _ => atype,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(gene: Option<GeneKey>) -> IsoformMatch {
        IsoformMatch::new(MatchClass::Undefined, gene, None, [], Strand::Unknown, 0.0)
    }

    #[test]
    fn gene_atype() {
        use AssignmentType as A;
        let g1 = Some(GeneKey::Ref(1));
        let g2 = Some(GeneKey::Ref(2));
        assert_eq!(derive_gene_atype(A::Ambiguous, &[m(g1), m(g1)]), A::Unique);
        assert_eq!(derive_gene_atype(A::Ambiguous, &[m(g1), m(g2)]), A::Ambiguous);
        assert_eq!(derive_gene_atype(A::Ambiguous, &[m(g1), m(None)]), A::Ambiguous);
        assert_eq!(derive_gene_atype(A::InconsistentAmbiguous, &[m(g1)]), A::Inconsistent);
        assert_eq!(derive_gene_atype(A::Noninformative, &[m(None)]), A::Noninformative);
        assert_eq!(derive_gene_atype(A::Intergenic, &[m(g1), m(None)]), A::InconsistentGenic);
        assert_eq!(derive_gene_atype(A::Noninformative, &[m(g1), m(g2)]), A::InconsistentMultigenic);
        assert_eq!(derive_gene_atype(A::Unique, &[m(g1), m(g2)]), A::Unique);
    }

    #[test]
    fn events_filter_and_add() {
        let mut x = IsoformMatch::new(
            MatchClass::Undefined,
            None,
            None,
            [MatchEvent::new(EventType::None), MatchEvent::new(EventType::Fsm)],
            Strand::Plus,
            0.0,
        );
        assert_eq!(x.events.len(), 1);
        let mut y = IsoformMatch::new_single(
            MatchClass::Undefined,
            None,
            None,
            MatchEvent::new(EventType::Undefined),
            Strand::Plus,
            0.0,
        );
        y.add_event(MatchEvent::new(EventType::AlignedPolyaTail));
        assert_eq!(y.events.len(), 1);
        x.add_event(MatchEvent::new(EventType::AlignedPolyaTail));
        assert_eq!(x.events.len(), 2);
        assert_eq!(MatchEvent::new(EventType::Fsm).iso_region.0, 1i64 << 31);
    }
}
