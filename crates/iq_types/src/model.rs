//! Construct-side types: intron-graph state and paths, transcript models and
//! the per-block `ModelStore`, terminal predictions (U10), final assignments.

use crate::assign::GeneKey;
use crate::enums::{AssignmentType, TranscriptModelType};
use crate::prim::{Iv, ModelIdx, Pos, ReadIdx, Strand, TxIdx};
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;

// ------------------------------------------------------------------------------------------
// Graph
// ------------------------------------------------------------------------------------------

/// Intron-graph vertex (`intron_graph.py:18-31`): an intron `(start, end)` with `start >= 0`,
/// or a terminal vertex `(code, pos)` with `code < 0` (see [`vertex_code`]).
///
/// Kept as a tuple so that `Ord` equals Python tuple ordering (needed for FL path sorting,
/// `fl_graph_constructor.py:94-95`). Serialized as `[a, b]`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Vertex(pub i64, pub i64);

impl Vertex {
    #[inline]
    pub const fn intron(iv: Iv) -> Self {
        Vertex(iv.s as i64, iv.e as i64)
    }
    #[inline]
    pub const fn is_terminal(&self) -> bool {
        self.0 < 0
    }
    /// The intron, if this is an intron vertex.
    pub const fn as_intron(&self) -> Option<Iv> {
        if self.0 >= 0 {
            Some(Iv::new(self.0 as Pos, self.1 as Pos))
        } else {
            None
        }
    }
}

/// Terminal vertex codes (`intron_graph.py:18-31`).
pub mod vertex_code {
    /// polyA (genomic right end).
    pub const POLYA: i64 = -10;
    pub const READ_END: i64 = -11;
    /// TSS on the genomic right (minus-strand TSS).
    pub const TSS_RIGHT: i64 = -12;
    /// polyT (genomic left end).
    pub const POLYT: i64 = -20;
    pub const READ_START: i64 = -21;
    /// TSS on the genomic left (plus-strand TSS).
    pub const TSS_LEFT: i64 = -22;
}

/// Graph-local vertex id, for unit-internal adjacency tables.
pub type VId = u32;

/// An intron path: a sequence of vertices (terminals included where threaded).
pub type Path = Box<[Vertex]>;

/// Intron-graph state right after `IntronGraph.__init__` (U11 output; schema §2.8 `Graph`).
/// Maps are insertion-ordered (= 4.1 dict order). Adjacency targets are sorted.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct IntronGraph {
    pub clustered_introns: IndexMap<Iv, u32>,
    pub intron_correction_map: IndexMap<Iv, Iv>,
    /// Sorted.
    pub discarded_introns: Vec<Iv>,
    pub outgoing: IndexMap<Vertex, Vec<Vertex>>,
    pub incoming: IndexMap<Vertex, Vec<Vertex>>,
    pub starting_known_positions: IndexMap<Iv, Vec<Pos>>,
    pub terminal_known_positions: IndexMap<Iv, Vec<Pos>>,
    pub max_coverage: u32,
}

/// Path data of one block (U12 output).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PathStorage {
    /// path → read count, dict order (= first read threading it).
    pub paths: IndexMap<Path, u32>,
    /// Full-length paths, sorted.
    pub fl_paths: Vec<Path>,
    /// path → block read indices, dict order, read order.
    pub paths_to_reads: IndexMap<Path, Vec<ReadIdx>>,
    /// `known_isoforms_in_graph`: path → known transcript, dict order.
    pub known_isoforms_in_graph: IndexMap<Path, TxIdx>,
}

// ------------------------------------------------------------------------------------------
// Terminal predictions (U10)
// ------------------------------------------------------------------------------------------

/// A detected peak (`terminal_peaks.py` `Peak`): `(position, count, left, right)`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Peak {
    pub position: Pos,
    pub count: u32,
    pub left: Pos,
    pub right: Pos,
}

/// Per-block predictions (`last_gene_predictions`), list order = 4.1 order.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Predictions {
    pub polya: Vec<(TxIdx, Pos)>,
    /// Empty without `--fl_data` (ONT).
    pub tss: Vec<(TxIdx, Pos)>,
}

// ------------------------------------------------------------------------------------------
// Transcript models
// ------------------------------------------------------------------------------------------

/// Model kind (`TranscriptModelType` + the known transcript).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum ModelKind {
    Known(TxIdx),
    Nic,
    Nnic,
}

impl ModelKind {
    pub const fn model_type(self) -> TranscriptModelType {
        match self {
            ModelKind::Known(_) => TranscriptModelType::Known,
            ModelKind::Nic => TranscriptModelType::NovelInCatalog,
            ModelKind::Nnic => TranscriptModelType::NovelNotInCatalog,
        }
    }
}

/// A transcript model (`gene_info.py:35-52` `TranscriptModel`).
///
/// Note: `intron_path` holds the vertex tuples themselves (not graph-local `VId`s), so a
/// model is meaningful without its graph; it is empty for known and assignment-based models
/// (4.1 `()`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TranscriptModel {
    pub key: ModelIdx,
    pub strand: Strand,
    pub gene: GeneKey,
    pub kind: ModelKind,
    pub exons: Vec<Iv>,
    pub intron_path: Path,
    pub polya_confirmed: bool,
    pub tss_confirmed: bool,
    /// For alt-end NICs derived from a known model.
    pub source_known: Option<TxIdx>,
}

impl TranscriptModel {
    pub fn start(&self) -> Pos {
        self.exons[0].s
    }
    pub fn end(&self) -> Pos {
        self.exons[self.exons.len() - 1].e
    }
    pub fn span(&self) -> Iv {
        Iv::new(self.start(), self.end())
    }
}

/// Per-block model storage (`model_construction/context.py:48-110` `ModelStore`).
///
/// Reads are keyed by **block read index** (D-M5), models by `ModelIdx` (= creation order).
/// Insertion-ordered maps mirror the 4.1 `defaultdict`s, including entries 4.1 creates on
/// read access (callers that emulate such accesses must insert explicitly).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ModelStore {
    /// `transcript_model_storage`, storage order.
    pub models: Vec<TranscriptModel>,
    /// `transcript_read_ids`: model → bound reads (append order; duplicates possible).
    pub transcript_read_ids: IndexMap<ModelIdx, Vec<ReadIdx>>,
    /// `internal_counter`.
    pub internal_counter: IndexMap<ModelIdx, i32>,
    /// `read_assignment_counts`, one slot per block read.
    pub read_assignment_counts: Vec<i32>,
    /// `construction_assignment`, one slot per block read.
    pub construction_assignment: Vec<Option<ModelIdx>>,
    /// `nonfl_known_reads`: known transcript → reads used to build it as a non-FL model.
    pub nonfl_known_reads: IndexMap<TxIdx, Vec<ReadIdx>>,
    /// Next fresh model key.
    pub next_key: ModelIdx,
}

impl ModelStore {
    /// Empty store for a block of `n_reads` reads.
    pub fn new(n_reads: usize) -> Self {
        ModelStore {
            read_assignment_counts: vec![0; n_reads],
            construction_assignment: vec![None; n_reads],
            ..Default::default()
        }
    }
    /// Allocates a fresh `ModelIdx` (creation order).
    pub fn new_key(&mut self) -> ModelIdx {
        let k = self.next_key;
        self.next_key += 1;
        k
    }
    pub fn model(&self, key: ModelIdx) -> Option<&TranscriptModel> {
        self.models.iter().find(|m| m.key == key)
    }
    /// `save_assigned_read` (`context.py:76-82`).
    pub fn save_assigned_read(&mut self, read: ReadIdx, model: ModelIdx) {
        self.transcript_read_ids.entry(model).or_default().push(read);
        *self.internal_counter.entry(model).or_default() += 1;
        self.read_assignment_counts[read as usize] += 1;
        self.construction_assignment[read as usize] = Some(model);
    }
    /// `move_read_to_model` (`context.py:84-93`): removes the first occurrence of `read`
    /// from `source` (`list.remove`).
    pub fn move_read_to_model(&mut self, read: ReadIdx, source: ModelIdx, target: ModelIdx) {
        let v = self.transcript_read_ids.entry(source).or_default();
        let pos = v.iter().position(|&r| r == read).expect("list.remove(x): x not in list");
        v.remove(pos);
        *self.internal_counter.entry(source).or_default() -= 1;
        self.read_assignment_counts[read as usize] -= 1;
        self.save_assigned_read(read, target);
    }
    /// `delete_from_storage` (`context.py:95-99`). Does not remove the model from `models`
    /// (4.1 callers do that separately).
    pub fn delete_from_storage(&mut self, model: ModelIdx) {
        let reads = self.transcript_read_ids.entry(model).or_default().clone();
        for r in reads {
            self.read_assignment_counts[r as usize] -= 1;
        }
        self.transcript_read_ids.shift_remove(&model);
        self.internal_counter.shift_remove(&model);
    }
}

/// Compact final per-read assignment against models (U19 → U20).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FinalAssign {
    pub atype: AssignmentType,
    pub gene_atype: AssignmentType,
    pub tx: SmallVec<[ModelIdx; 1]>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vertex_order_is_python_tuple_order() {
        let mut v = vec![Vertex(100, 200), Vertex(-10, 500), Vertex(-21, 3), Vertex(100, 150)];
        v.sort();
        assert_eq!(v, vec![Vertex(-21, 3), Vertex(-10, 500), Vertex(100, 150), Vertex(100, 200)]);
        assert_eq!(serde_json::to_string(&Vertex(-10, 5)).unwrap(), "[-10,5]");
    }

    #[test]
    fn store_bookkeeping() {
        let mut s = ModelStore::new(3);
        let a = s.new_key();
        let b = s.new_key();
        s.save_assigned_read(0, a);
        s.save_assigned_read(1, a);
        s.move_read_to_model(0, a, b);
        assert_eq!(s.transcript_read_ids[&a], vec![1]);
        assert_eq!(s.internal_counter[&a], 1);
        assert_eq!(s.read_assignment_counts, vec![1, 1, 0]);
        s.delete_from_storage(a);
        assert_eq!(s.read_assignment_counts, vec![1, 0, 0]);
        assert!(!s.transcript_read_ids.contains_key(&a));
    }
}
