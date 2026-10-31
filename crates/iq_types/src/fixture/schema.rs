//! Serde DTOs that mirror `harness/iqcap/SCHEMA.md` (v1) **exactly**: field names, field
//! order (= JSON key order), value encodings. They carry IsoQuant id *strings*; convert to
//! and from the core types with [`super::convert`] and an [`super::IdMaps`].
//!
//! Naming: shared object types end in `Fixture` (`ReadRecordFixture`, ...); per-unit
//! payloads are `<Unit><Fn>In` / `<Unit><Fn>Out`. Every `In` starts with `fn_` (JSON `"fn"`).
//! Core types are reused where the JSON shape is identical: `Iv`, `Iv0`, `Strand`,
//! `PolyAInfo` (`[4 ints]`), the enums (integer codes), `Vertex` (`[a,b]`),
//! `MappedReadProfile` (`Prof`) and `CombinedProfile` (`CProf`).
//!
//! Section numbers below refer to SCHEMA.md.

use super::pyjson::PyF64;
use crate::enums::{AmbiguityResolving, AssignmentType, EventType, MatchClass, MultimapStrategy, TranscriptModelType};
use crate::index::{CombinedProfile, MappedReadProfile};
use crate::model::Vertex;
use crate::prim::{Iv, Iv0, PolyAInfo, Pos, Strand};
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};

/// `[key, value]` pair list in Python dict order.
pub type Pairs<K, V> = Vec<(K, V)>;
/// CIGAR as `[[op, len]]`.
pub type CigarFixture = Vec<(u8, u32)>;
/// `(start, end, featuretype)` triples (`other_features`).
pub type OtherFeatures = Vec<(Pos, Pos, String)>;
/// An intron-graph path `[[a,b], ...]`.
pub type PathFixture = Vec<Vertex>;

// =========================================================================================
// §0.2 envelope exceptions
// =========================================================================================

/// One fixture line. `I`/`O` default to raw JSON.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Case<I = serde_json::Value, O = serde_json::Value> {
    pub unit: String,
    pub v: u32,
    pub case: String,
    #[serde(rename = "in")]
    pub input: I,
    pub out: O,
}

/// `out` of a call that raised.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExceptionOut {
    pub exception: String,
}

// =========================================================================================
// §1 manifest (read-only; unknown fields ignored, not round-tripped)
// =========================================================================================

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Manifest {
    pub v: u32,
    #[serde(default)]
    pub sample: String,
    #[serde(default)]
    pub slice: String,
    #[serde(default)]
    pub isoquant_args: Vec<String>,
    #[serde(default)]
    pub units: IndexMap<String, ManifestUnit>,
    #[serde(default)]
    pub blocks: Vec<ManifestBlock>,
    #[serde(default)]
    pub skipped_blocks: Vec<ManifestSkipped>,
    #[serde(default)]
    pub cpu: serde_json::Value,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ManifestUnit {
    #[serde(default)]
    pub files: Vec<String>,
    #[serde(default)]
    pub cases: u64,
    #[serde(default)]
    pub by_fn: IndexMap<String, u64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ManifestBlock {
    pub chr: String,
    pub block_no: u32,
    pub region: Iv,
    pub n_reads: u64,
    pub captured: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ManifestSkipped {
    pub chr: String,
    pub block_no: u32,
    pub reason: String,
    pub n_reads: u64,
}

// =========================================================================================
// §2 shared objects
// =========================================================================================

/// §2.1 `index/<sha1>.json`: a dumped `GeneInfo`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct IndexDump {
    /// `"genes"`, `"models"`, `"model"`, `"region"`.
    pub kind: String,
    pub chr: Option<String>,
    pub start: Pos,
    pub end: Pos,
    pub delta: i32,
    pub genes: Vec<String>,
    pub gene_regions: Pairs<String, Iv>,
    pub isoform_exons: Pairs<String, Vec<Iv>>,
    pub isoform_introns: Pairs<String, Vec<Iv>>,
    pub isoform_strands: Pairs<String, Strand>,
    pub gene_strands: Pairs<String, Strand>,
    pub gene_id_map: Pairs<String, String>,
    pub intron_features: Vec<Iv>,
    pub intron_profiles: Pairs<String, Vec<i8>>,
    pub intron_ranges: Pairs<String, (i32, i32)>,
    pub exon_features: Vec<Iv>,
    pub exon_profiles: Pairs<String, Vec<i8>>,
    pub exon_ranges: Pairs<String, (i32, i32)>,
    pub split_exon_features: Vec<Iv>,
    pub split_exon_profiles: Pairs<String, Vec<i8>>,
    pub split_exon_ranges: Pairs<String, (i32, i32)>,
    pub sources: Pairs<String, String>,
    pub other_features: Pairs<String, OtherFeatures>,
    pub feature_attributes: Pairs<String, String>,
}

/// §2.2 header line of `block/<sample>/<chr>/<block_no>.jsonl.gz`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BlockHeader {
    /// Always `"block_header"`.
    pub kind: String,
    pub v: u32,
    pub sample: String,
    pub slice: String,
    pub chr: String,
    pub block_no: u32,
    pub region: Iv,
    pub index: String,
    pub gene_ids: Vec<String>,
    pub n_reads: u32,
    pub collect_region: Option<Iv>,
    pub skip_read_fraction: u32,
}

/// §2.2 `ReadRecord`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ReadRecordFixture {
    pub idx: u32,
    pub assignment_id: i64,
    pub name: String,
    pub chr: String,
    pub region: Iv,
    pub exons: Vec<Iv>,
    pub corrected_exons: Vec<Iv>,
    pub boundary_deletions: Vec<(Pos, u32)>,
    pub polya: PolyAInfo,
    pub polya_found: bool,
    pub cage_found: bool,
    pub strand: Strand,
    pub mapped_strand: Strand,
    pub mapq: u8,
    pub multimapper: bool,
    pub read_group: Vec<String>,
    pub assignment: AssignmentFixture,
    pub introns_match: bool,
    pub additional_attributes: IndexMap<String, String>,
}

/// §2.3 inline `params` of U5/U6/U7.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AssignParamsFixture {
    pub delta: i32,
    pub minor_exon_extension: i32,
    pub major_exon_extension: i32,
    pub min_abs_exon_overlap: i32,
    pub min_rel_exon_overlap: f64,
    pub minimal_exon_overlap: i32,
    pub minimal_intron_absence_overlap: i32,
    pub micro_intron_length: i32,
    pub max_fake_terminal_exon_len: i32,
    pub max_intron_shift: i32,
    pub max_missed_exon_len: i32,
    pub max_intron_abs_diff: i32,
    pub max_intron_rel_diff: f64,
    pub max_suspicious_intron_abs_len: i32,
    pub max_suspicious_intron_rel_len: f64,
    /// Integer code (agreed with the schema owner, 2026-10-29).
    pub resolve_ambiguous: AmbiguityResolving,
    pub apa_delta: i32,
    pub count_exons: bool,
    pub correct_fuzzy_junctions: bool,
    pub correct_intron_shifts: bool,
    pub correct_skipped_exons: bool,
    pub correct_fake_terminal_exons: bool,
    pub correct_terminal_exons: bool,
    pub correct_microintron_retention: bool,
}

/// §2.4 reference window: `seq = chr[start0 : start0+len(seq)]` (0-based half-open).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct RefWin {
    pub start0: i64,
    pub seq: String,
}

impl RefWin {
    /// Borrowed 1-based view (`RefSlice.start = start0 + 1`).
    pub fn as_slice(&self) -> crate::read::RefSlice<'_> {
        crate::read::RefSlice { start: (self.start0 + 1) as Pos, seq: self.seq.as_bytes() }
    }
}

/// §2.5 `Assignment`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AssignmentFixture {
    pub atype: AssignmentType,
    pub gene_atype: AssignmentType,
    pub matches: Vec<IsoformMatchFixture>,
}

/// §2.5 `IsoformMatch`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct IsoformMatchFixture {
    pub gene: Option<String>,
    pub tx: Option<String>,
    pub strand: Strand,
    pub class: MatchClass,
    pub penalty: f64,
    pub events: Vec<MatchEventFixture>,
}

/// §2.5 `MatchEvent`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct MatchEventFixture {
    pub t: EventType,
    pub iso: (i64, i64),
    pub read: (i64, i64),
    pub info: i64,
}

/// §2.6 `Model` (TranscriptModel).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ModelFixture {
    pub id: String,
    pub gene: String,
    pub chr: String,
    pub strand: Strand,
    #[serde(rename = "type")]
    pub ty: TranscriptModelType,
    pub exons: Vec<Iv>,
    pub intron_path: PathFixture,
    pub polya_confirmed: bool,
    pub tss_confirmed: bool,
    pub source: String,
    pub other_features: OtherFeatures,
    pub additional_info: Pairs<String, String>,
}

/// §2.7 `Store` (ModelStore snapshot).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StoreFixture {
    pub models: Vec<ModelFixture>,
    pub transcript_read_ids: Pairs<String, Vec<u32>>,
    pub internal_counter: Pairs<String, i64>,
    /// Keyed by read **name** (4.1).
    pub read_assignment_counts: Pairs<String, i64>,
    /// read name → model id.
    pub construction_assignment: Pairs<String, String>,
    /// `read_assignment_counts` re-keyed by `read_idx`: per name entry (dict order), one
    /// pair per block read with that name (ascending idx); duplicates share the value.
    /// Optional only for tolerance of pre-change captures.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub read_assignment_counts_idx: Option<Pairs<u32, i64>>,
    /// `construction_assignment` re-keyed the same way.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub construction_assignment_idx: Option<Pairs<u32, String>>,
    pub nonfl_known_reads: Pairs<String, Vec<u32>>,
    pub detected_known_isoforms: Vec<String>,
    pub id_counter: i64,
}

/// §2.8 `Graph` (IntronGraph after `__init__`).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct GraphFixture {
    pub clustered_introns: Pairs<Iv, u32>,
    pub intron_correction_map: Pairs<Iv, Iv>,
    pub discarded_introns: Vec<Iv>,
    pub outgoing: Pairs<Vertex, Vec<Vertex>>,
    pub incoming: Pairs<Vertex, Vec<Vertex>>,
    pub starting_known_positions: Pairs<Iv, Vec<Pos>>,
    pub terminal_known_positions: Pairs<Iv, Vec<Pos>>,
    pub max_coverage: u32,
}

/// §2.9 block reference.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct BlockRef {
    pub sample: String,
    pub slice: String,
    pub chr: String,
    pub block_no: u32,
}

/// §0.3 read identity: collect side `{"qname","aln_no"}`, construct side `{"read_idx"}`,
/// FL/similar `{"id"}`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum ReadRef {
    Aln { qname: String, aln_no: u32 },
    Idx { read_idx: u32 },
    Id { id: String },
}

/// `{"read_idx": int}`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReadIdxRef {
    pub read_idx: u32,
}

/// §4 `Aln`: an alignment as IsoQuant sees it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlnFixture {
    pub qname: String,
    pub flag: u16,
    pub pos0: i64,
    pub ref_end0: i64,
    pub mapq: u8,
    pub is_reverse: bool,
    pub cigar: CigarFixture,
    pub seq: Option<String>,
}

impl AlnFixture {
    /// Raw BAM CIGAR words.
    pub fn cigar_raw(&self) -> Vec<u32> {
        self.cigar.iter().map(|&(op, len)| (len << 4) | op as u32).collect()
    }
    /// Owned record (no tag).
    pub fn to_record(&self) -> crate::read::RecordOwned {
        crate::read::RecordOwned {
            qname: self.qname.clone(),
            flag: self.flag,
            pos0: self.pos0,
            mapq: self.mapq,
            cigar: self.cigar_raw(),
            seq: self.seq.clone().unwrap_or_default().into_bytes(),
            tag: None,
        }
    }
}

/// `Prof` = 4.1 `MappedReadProfile` (same JSON shape as the core type).
pub type ProfFixture = MappedReadProfile;
/// `CProf` = 4.1 `CombinedReadProfiles` (same JSON shape as the core type).
pub type CProfFixture = CombinedProfile;

// =========================================================================================
// §4 per-unit payloads
// =========================================================================================

// ---- U1 ---------------------------------------------------------------------------------

/// U1 (a) `fn = "dump_annotation"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U1DumpAnnotationIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub db: String,
    pub gff3: String,
    pub chr: String,
    pub gene_types: String,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U1DumpAnnotationOut {
    pub chr: String,
    pub genes: Vec<GeneFixture>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GeneFixture {
    pub id: String,
    pub start: Pos,
    pub end: Pos,
    pub strand: Strand,
    pub source: String,
    pub attributes: Pairs<String, Vec<String>>,
    pub attrs: String,
    pub transcripts: Vec<TxFixture>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TxFixture {
    pub id: String,
    #[serde(rename = "type")]
    pub ty: String,
    pub start: Pos,
    pub end: Pos,
    pub strand: Strand,
    pub source: String,
    pub attributes: Pairs<String, Vec<String>>,
    pub attrs: String,
    pub exons: Vec<Iv>,
    pub other_features: OtherFeatures,
    pub children_order_index: u32,
}

/// U1 (b) `fn = "get_gene_info_for_region"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U1GeneInfoForRegionIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub chr: String,
    /// As passed by 4.1 (0-based closed storage region).
    pub region: Iv0,
    pub delta: i32,
    pub needs_reference: bool,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U1GeneInfoForRegionOut {
    pub gene_ids: Vec<String>,
    pub index: String,
    pub empty: bool,
}

/// U1 (c) `fn = "from_models"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U1FromModelsIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub models: Vec<ModelFixture>,
    pub delta: i32,
    pub call_site: String,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct IndexOut {
    pub index: String,
}

// ---- U2 ---------------------------------------------------------------------------------

/// U2 `fn = "fetch"` (`iqcap.dump_bam`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U2FetchIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub bam: String,
    /// 1-based closed samtools region string, e.g. `"8:1001-2000"`.
    pub region: String,
    pub n: u64,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U2FetchOut {
    pub qname: String,
    pub flag: u16,
    pub ref_id: i32,
    pub ref_name: Option<String>,
    pub pos0: i64,
    pub ref_end0: Option<i64>,
    pub mapq: u8,
    pub cigar: CigarFixture,
    pub seq: Option<String>,
    pub qual_len: i64,
    pub tags: IndexMap<String, TagFixture>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TagFixture {
    /// SAM type letter (`"Z"`, `"i"`, `"A"`, `"f"`, ...).
    #[serde(rename = "type")]
    pub ty: String,
    pub value: TagValueFixture,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum TagValueFixture {
    Int(i64),
    Float(f64),
    Str(String),
}

// ---- U3 ---------------------------------------------------------------------------------

/// U3 `fn = "split_coverage_regions"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U3SplitCoverageIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub chr: String,
    pub region: Iv0,
    pub read_count: u64,
    /// `[bin, count]`, sorted by bin.
    pub coverage: Pairs<i64, u64>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U3SplitCoverageOut {
    pub splits: Vec<(Iv0, u64)>,
}

/// U3 `fn = "process_alignments_in_region"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U3ProcessRegionIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub chr: String,
    pub region: Iv0,
    pub max_coverage: u64,
    pub chr_len: u64,
    pub small_chr: bool,
    pub coverage_cutoff: i64,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U3ProcessRegionOut {
    pub block_no: u32,
    pub skip_read_fraction: u32,
    pub gene_ids: Vec<String>,
    pub index: String,
    pub genic: bool,
}

/// U3 `fn = "cluster"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U3ClusterIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub chr: String,
    pub region: Iv0,
    pub read_count: u64,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U3ClusterOut {
    pub n_blocks: u32,
    pub first_block_no: u32,
}

/// U3 `fn = "genes_in_region"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U3GenesInRegionIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub chr: String,
    pub region: Iv0,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GeneIdsOut {
    pub gene_ids: Vec<String>,
}

// ---- U4 ---------------------------------------------------------------------------------

/// U4 `fn = "process_alignment"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U4ProcessAlignmentIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub read: ReadRef,
    pub aln: AlnFixture,
    pub polya_window: i32,
    pub polya_fraction: f64,
    pub max_fake_terminal_exon_len: i32,
    pub external_strand: Option<String>,
    pub external_no_tail: bool,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U4BlocksFixture {
    pub read_exons: Vec<Iv>,
    pub read_blocks: Vec<(i32, i32)>,
    pub cigar_blocks: Vec<(i32, i32)>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U4ProcessAlignmentOut {
    pub init: U4BlocksFixture,
    pub polya: Option<PolyAInfo>,
    pub read_exons: Vec<Iv>,
    pub read_blocks: Vec<(i32, i32)>,
    pub cigar_blocks: Vec<(i32, i32)>,
    pub exons_changed: bool,
    pub polya_called: bool,
}

/// U4 `fn = "extract_relevant_deletions"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U4ExtractDeletionsIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub read: ReadRef,
    pub cigar: CigarFixture,
    pub ref_start: i64,
    pub exons: Vec<Iv>,
    pub margin: i32,
    /// Alignment's 0-based `reference_start` (schema v1 addendum, capture batch 2; absent in older dirs).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pos0: Option<i64>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U4ExtractDeletionsOut {
    pub deletions: Vec<(Pos, u32)>,
}

// ---- U5 ---------------------------------------------------------------------------------

/// U5 `fn = "construct_profiles"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U5ConstructProfilesIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub call_site: String,
    /// `null` for provisional models.
    pub read: Option<ReadRef>,
    pub index: String,
    pub params: AssignParamsFixture,
    pub blocks: Vec<Iv>,
    pub polya: PolyAInfo,
}
/// `CProf` without `polya`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U5ConstructProfilesOut {
    pub intron: ProfFixture,
    pub exon: Option<ProfFixture>,
    pub split_exon: ProfFixture,
}

/// U5 `fn = "match_genomic_features"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U5MatchGenomicFeaturesIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub call_site: String,
    pub read: ReadRef,
    pub known_features: Vec<Iv>,
    pub delta: i32,
    pub read_features: Vec<Iv>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U5MatchGenomicFeaturesOut {
    pub features: Vec<Iv>,
}

// ---- U6 ---------------------------------------------------------------------------------

/// U6 `fn = "assign_to_isoform"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U6AssignIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    /// `collect`, `fl`, `similar`, `reassign`, `final_single`, `final_full`, ...
    pub call_site: String,
    pub read: ReadRef,
    pub read_id: String,
    pub index: String,
    pub params: AssignParamsFixture,
    pub quick_mode: bool,
    pub profile: CProfFixture,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U6AssignOut {
    pub assignment: Option<AssignmentFixture>,
}

// ---- U7 ---------------------------------------------------------------------------------

/// U7 `fn = "correct_assigned_read"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U7CorrectIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub read: ReadRef,
    pub index: String,
    pub params: AssignParamsFixture,
    pub aln: AlnFixture,
    pub read_exons: Vec<Iv>,
    pub read_start: Pos,
    pub read_end: Pos,
    pub intron_profile: ProfFixture,
    pub assignment: AssignmentFixture,
    /// One window per read exon, `[s-32, e+32]` (1-based), merged and clipped; `[]` without
    /// a reference.
    #[serde(rename = "ref")]
    pub ref_: Vec<RefWin>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U7CorrectOut {
    pub corrected_exons: Vec<Iv>,
}

/// U7 `fn = "get_assignment_strand"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U7StrandIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub read: ReadRef,
    pub assignment: AssignmentFixture,
    pub polya: PolyAInfo,
    pub exons: Vec<Iv>,
    pub corrected_introns: Vec<Iv>,
    /// `(intron, donor, acceptor)` dinucleotides, verbatim.
    pub intron_sites: Vec<(Iv, String, String)>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U7StrandOut {
    pub strand: Strand,
}

/// U7 `fn = "check_antisense"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U7AntisenseIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub read: ReadRef,
    pub strand: Strand,
    pub assignment: AssignmentFixture,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AssignmentOut {
    pub assignment: AssignmentFixture,
}

// ---- U8 ---------------------------------------------------------------------------------

/// U8 `fn = "collect_block"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U8CollectBlockIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub chr: String,
    pub block_no: u32,
    pub region: Iv0,
    pub index: String,
    pub genic: bool,
    pub skip_read_fraction: u32,
    pub max_coverage: u64,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U8CollectBlockOut {
    pub reads: Vec<ReadRecordFixture>,
}

// ---- U9 ---------------------------------------------------------------------------------

/// U9 `fn = "resolve"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U9ResolveIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub strategy: MultimapStrategy,
    pub assignments: Vec<BasicAssignmentFixture>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U9ResolveOut {
    pub assignments: Vec<BasicAssignmentFixture>,
}
/// `BasicReadAssignment`. `genes`/`isoforms` are set iteration order: compare as sets.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BasicAssignmentFixture {
    pub assignment_id: i64,
    pub read_id: String,
    pub chr: String,
    pub start: Pos,
    pub end: Pos,
    pub region: Iv,
    pub multimapper: bool,
    pub polya_found: bool,
    pub atype: AssignmentType,
    pub gene_atype: AssignmentType,
    pub penalty: f64,
    pub genes: Vec<String>,
    pub isoforms: Vec<String>,
}

// ---- U10a -------------------------------------------------------------------------------

/// U10a `fn = "detect_peaks_batch"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U10aDetectPeaksIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    /// `"polya"` or `"tss"`.
    pub model: String,
    /// Each histogram as `[[pos, count]]` in dict order (load-bearing), or null.
    pub histograms: Vec<Option<Pairs<Pos, u32>>>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U10aDetectPeaksOut {
    pub records: Vec<Option<PeakRecordFixture>>,
    pub rows: Vec<PeakRowFixture>,
    /// Per histogram: `[[position, count, left, right]]`.
    pub peaks: Vec<Vec<(Pos, u32, Pos, Pos)>>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PeakRecordFixture {
    pub start: Pos,
    pub mode: Pos,
    pub var: f64,
    pub skew: PyF64,
    pub entropy: f64,
    pub mean_height: f64,
    pub peak_count: u32,
    /// Padded histogram.
    pub histogram: Vec<u32>,
    pub raw_peaks: Option<RawPeaksFixture>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RawPeaksFixture {
    pub peak_location: Vec<i64>,
    pub peak_width: Vec<f64>,
    pub peak_left: Vec<i64>,
    pub peak_right: Vec<i64>,
    pub peak_heights: Vec<i64>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PeakRowFixture {
    pub hist: u32,
    /// 8 features in FEATURE_COLUMNS order (f64 before the f32 cast).
    pub features: Vec<PyF64>,
    pub proba: f64,
    pub margin: f64,
    pub label: i32,
}

/// U10a `fn = "xgb_predict"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U10aXgbPredictIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub model: String,
    pub caller: String,
    pub features: Vec<Vec<PyF64>>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U10aXgbPredictOut {
    pub labels: Vec<i32>,
    pub proba: Vec<f64>,
    pub margin: Vec<f64>,
}

// ---- U10b -------------------------------------------------------------------------------

/// U10b `fn = "flush"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U10bFlushIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    /// `"polya"` or `"tss"`.
    pub kind: String,
    pub block: BlockRef,
    pub transcripts: Pairs<String, TerminalAccFixture>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TerminalAccFixture {
    pub chr: String,
    pub gene_id: String,
    /// Read positions in add order.
    pub data: Vec<Pos>,
    pub annotated: i64,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PredictionsOut {
    pub predictions: Pairs<String, Pos>,
}

// ---- U11 / U12 --------------------------------------------------------------------------

/// U11 `fn = "IntronGraph.__init__"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U11GraphIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub block: BlockRef,
    pub index: String,
    pub polya_predictions: Option<Pairs<String, Pos>>,
    pub tss_predictions: Option<Pairs<String, Pos>>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U11GraphOut {
    pub graph: GraphFixture,
}

/// U12 `fn = "IntronPathStorage.fill"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U12FillIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub block: BlockRef,
    /// Case key of the U11 case holding the graph.
    pub graph: String,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U12FillOut {
    pub paths: Pairs<PathFixture, u32>,
    pub fl_paths: Vec<PathFixture>,
    pub paths_to_reads: Pairs<PathFixture, Vec<u32>>,
}

/// U12 `fn = "get_known_spliced_isoforms"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U12KnownSplicedIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub block: BlockRef,
    pub graph: String,
    pub index: String,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U12KnownSplicedOut {
    pub known: Pairs<PathFixture, String>,
}

// ---- U13-U19 ----------------------------------------------------------------------------

/// Common input of the store-snapshot stages (all U13-U19 inputs begin with
/// `fn, block, index, graph, paths`; U18 has no graph/paths) (U13, U14, U15a pre_filter/filter, U15b
/// add_known/drop_duplicate, U16 correct_transcript_splice_sites, U17): `block`, `index`,
/// optional `graph`/`paths` case keys, `store`, and optional `reads`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StageStoreIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub block: BlockRef,
    pub index: String,
    /// U11 case key, or null.
    pub graph: Option<String>,
    /// U12 fill case key, or null.
    pub paths: Option<String>,
    pub store: StoreFixture,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reads: Option<Vec<u32>>,
    /// U13 only (batch 2): `[vertex, donor, acceptor]` per vertex of each FL path.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub intron_sites: Option<Vec<(Iv, String, String)>>,
    /// U16 store-level only (batch 2): merged reference windows around exon boundaries.
    #[serde(default, rename = "ref", skip_serializing_if = "Option::is_none")]
    pub ref_: Option<Vec<RefWin>>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StoreOut {
    pub store: StoreFixture,
}

/// U15b per-model functions (`correct_novel_transcript_ends`, `_mark_terminal_confirmation`,
/// `derive_alternative_end_models`): model + reads (+ `id_counter` for derive), no store.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StageModelIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub block: BlockRef,
    pub index: String,
    /// U11 case key, or null.
    pub graph: Option<String>,
    /// U12 fill case key, or null.
    pub paths: Option<String>,
    pub model: ModelFixture,
    pub reads: Vec<u32>,
    /// Only for `derive_alternative_end_models` (value before the call).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id_counter: Option<i64>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ModelOut {
    pub model: ModelFixture,
}
/// U15b `derive_alternative_end_models` output.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U15bDeriveAltEndOut {
    pub result: Vec<(ModelFixture, Vec<u32>)>,
    pub id_counter: i64,
}

/// U15a `detect_similar_isoforms`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U15aDetectSimilarIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub block: BlockRef,
    pub index: String,
    /// U11 case key, or null.
    pub graph: Option<String>,
    /// U12 fill case key, or null.
    pub paths: Option<String>,
    pub models: Vec<ModelFixture>,
    pub remove_unconfirmed: bool,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U15aDetectSimilarOut {
    pub to_substitute: Pairs<String, String>,
}

/// U16 `_correct_one_model_splice_sites`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U16CorrectOneIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub block: BlockRef,
    pub index: String,
    pub graph: Option<String>,
    pub paths: Option<String>,
    pub model: ModelFixture,
    pub reads: Vec<U16ReadFixture>,
    /// One window `[b-65, b+64)` per exon boundary.
    #[serde(rename = "ref")]
    pub ref_: Vec<RefWin>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct U16ReadFixture {
    pub idx: u32,
    pub start: Pos,
    pub end: Pos,
    pub boundary_deletions: Vec<(Pos, u32)>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U16CorrectOneOut {
    pub exons: Option<Vec<Iv>>,
}

/// U18 `join_transcripts`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U18JoinIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub block: BlockRef,
    pub index: String,
    pub models: Vec<ModelFixture>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ModelsOut {
    pub models: Vec<ModelFixture>,
}

/// U19 `build_model_read_assignments` output (input is [`StageStoreIn`] with `reads`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U19AssignmentsOut {
    pub assignments: Vec<U19ReadAssignment>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U19ReadAssignment {
    pub read_idx: u32,
    pub assignment: AssignmentFixture,
}

// ---- U20 --------------------------------------------------------------------------------

/// U20 `fn = "gff_dump"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U20GffDumpIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub block: BlockRef,
    /// `"transcript_models"` or `"extended_annotation"`.
    pub printer: String,
    pub index: String,
    pub models: Vec<ModelFixture>,
    pub printed_gene_ids: Vec<String>,
    pub exon_id_counter: i64,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U20GffDumpOut {
    pub lines: Vec<String>,
    pub exon_id_counter: i64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CounterFixture {
    pub name: String,
    pub grouped: bool,
    pub group_index: i64,
}

/// U20 `fn = "counter_add_read_info"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U20CounterAddIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub counter: CounterFixture,
    pub read: Option<ReadIdxRef>,
    pub assignment: AssignmentFixture,
    pub read_group: Vec<String>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U20CounterAddOut {
    /// feature → `[[group, before, after]]`.
    pub features: Pairs<String, Vec<(String, f64, f64)>>,
    pub ambiguous: (i64, i64),
    pub not_assigned: (i64, i64),
    pub confirmed_added: Vec<String>,
}

/// U20 `fn = "counter_dump"`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U20CounterDumpIn {
    #[serde(rename = "fn")]
    pub fn_: String,
    pub counter: CounterFixture,
    pub feature_counter: Pairs<String, Vec<(String, f64)>>,
    pub all_features: Vec<String>,
    pub confirmed_features: Vec<String>,
    pub ambiguous: i64,
    pub not_assigned: i64,
    pub not_aligned: i64,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct U20CounterDumpOut {
    /// `[[basename, text]]`.
    pub files: Pairs<String, String>,
}

/// `(unit, fn)` → Rust payload types, for reference and for the generic round-trip test.
/// Units with one hook list it under its `fn` name.
pub const PAYLOAD_TABLE: &[(&str, &str, &str, &str)] = &[
    ("U1", "dump_annotation", "U1DumpAnnotationIn", "U1DumpAnnotationOut"),
    ("U1", "get_gene_info_for_region", "U1GeneInfoForRegionIn", "U1GeneInfoForRegionOut"),
    ("U1", "from_models", "U1FromModelsIn", "IndexOut"),
    ("U2", "fetch", "U2FetchIn", "U2FetchOut"),
    ("U3", "split_coverage_regions", "U3SplitCoverageIn", "U3SplitCoverageOut"),
    ("U3", "process_alignments_in_region", "U3ProcessRegionIn", "U3ProcessRegionOut"),
    ("U3", "cluster", "U3ClusterIn", "U3ClusterOut"),
    ("U3", "genes_in_region", "U3GenesInRegionIn", "GeneIdsOut"),
    ("U4", "process_alignment", "U4ProcessAlignmentIn", "U4ProcessAlignmentOut"),
    ("U4", "extract_relevant_deletions", "U4ExtractDeletionsIn", "U4ExtractDeletionsOut"),
    ("U5", "construct_profiles", "U5ConstructProfilesIn", "U5ConstructProfilesOut"),
    ("U5", "match_genomic_features", "U5MatchGenomicFeaturesIn", "U5MatchGenomicFeaturesOut"),
    ("U6", "assign_to_isoform", "U6AssignIn", "U6AssignOut"),
    ("U7", "correct_assigned_read", "U7CorrectIn", "U7CorrectOut"),
    ("U7", "get_assignment_strand", "U7StrandIn", "U7StrandOut"),
    ("U7", "check_antisense", "U7AntisenseIn", "AssignmentOut"),
    ("U8", "collect_block", "U8CollectBlockIn", "U8CollectBlockOut"),
    ("U9", "resolve", "U9ResolveIn", "U9ResolveOut"),
    ("U10a", "detect_peaks_batch", "U10aDetectPeaksIn", "U10aDetectPeaksOut"),
    ("U10a", "xgb_predict", "U10aXgbPredictIn", "U10aXgbPredictOut"),
    ("U10b", "flush", "U10bFlushIn", "PredictionsOut"),
    ("U11", "IntronGraph.__init__", "U11GraphIn", "U11GraphOut"),
    ("U12", "IntronPathStorage.fill", "U12FillIn", "U12FillOut"),
    ("U12", "get_known_spliced_isoforms", "U12KnownSplicedIn", "U12KnownSplicedOut"),
    ("U13", "construct_fl_isoforms", "StageStoreIn", "StoreOut"),
    ("U14", "construct_assignment_based_isoforms", "StageStoreIn", "StoreOut"),
    ("U15a", "pre_filter_transcripts", "StageStoreIn", "StoreOut"),
    ("U15a", "filter_transcripts", "StageStoreIn", "StoreOut"),
    ("U15a", "detect_similar_isoforms", "U15aDetectSimilarIn", "U15aDetectSimilarOut"),
    ("U15b", "correct_novel_transcript_ends", "StageModelIn", "ModelOut"),
    ("U15b", "_mark_terminal_confirmation", "StageModelIn", "ModelOut"),
    ("U15b", "_add_known_alternative_end_models", "StageStoreIn", "StoreOut"),
    ("U15b", "_drop_duplicate_alt_end_models", "StageStoreIn", "StoreOut"),
    ("U15b", "derive_alternative_end_models", "StageModelIn", "U15bDeriveAltEndOut"),
    ("U16", "correct_transcript_splice_sites", "StageStoreIn", "StoreOut"),
    ("U16", "_correct_one_model_splice_sites", "U16CorrectOneIn", "U16CorrectOneOut"),
    ("U17", "assign_reads_to_models", "StageStoreIn", "StoreOut"),
    ("U17", "resolve_read_ownership_conflicts", "StageStoreIn", "StoreOut"),
    ("U18", "join_transcripts", "U18JoinIn", "ModelsOut"),
    ("U19", "build_model_read_assignments", "StageStoreIn", "U19AssignmentsOut"),
    ("U20", "gff_dump", "U20GffDumpIn", "U20GffDumpOut"),
    ("U20", "counter_add_read_info", "U20CounterAddIn", "U20CounterAddOut"),
    ("U20", "counter_dump", "U20CounterDumpIn", "U20CounterDumpOut"),
];
