//! Conversions between fixture DTOs ([`super::schema`]) and core types.
//!
//! Fields that exist only in the fixture form (4.1 bookkeeping such as `assignment_id`,
//! `cage_found`, GFF `source`) are returned as `*Extras` so a core value can be turned back
//! into a byte-identical DTO.

use super::ids::{IdMaps, TxSpace};
use super::schema::*;
use crate::assign::{Assignment, IsoformMatch, MatchEvent};
use crate::enums::TranscriptModelType;
use crate::index::{GeneRec, IsoRec, IsoformIndex, ProfMatrix};
use crate::model::{IntronGraph, ModelKind, ModelStore, PathStorage, TranscriptModel};
use crate::prim::{junctions_from_blocks, Iv, Strand};
use crate::read::ReadRecord;
use indexmap::IndexMap;
use rustc_hash::FxHashMap;

// ------------------------------------------------------------------------------------------
// MatchEvent / IsoformMatch / Assignment
// ------------------------------------------------------------------------------------------

impl From<MatchEventFixture> for MatchEvent {
    fn from(e: MatchEventFixture) -> Self {
        MatchEvent { ty: e.t, iso_region: e.iso, read_region: e.read, info: e.info }
    }
}
impl From<MatchEvent> for MatchEventFixture {
    fn from(e: MatchEvent) -> Self {
        MatchEventFixture { t: e.ty, iso: e.iso_region, read: e.read_region, info: e.info }
    }
}

impl IsoformMatchFixture {
    pub fn to_core(&self, ids: &mut IdMaps, space: TxSpace) -> IsoformMatch {
        IsoformMatch {
            gene: self.gene.as_deref().map(|g| ids.gene_key(g)),
            tx: self.tx.as_deref().map(|t| ids.iso_ref(t, space)),
            tx_strand: self.strand,
            class: self.class,
            penalty: self.penalty,
            events: self.events.iter().map(|&e| e.into()).collect(),
        }
    }
    pub fn from_core(m: &IsoformMatch, ids: &IdMaps) -> Self {
        IsoformMatchFixture {
            gene: m.gene.map(|g| ids.gene_name(g).to_owned()),
            tx: m.tx.map(|t| ids.iso_name(t).to_owned()),
            strand: m.tx_strand,
            class: m.class,
            penalty: m.penalty,
            events: m.events.iter().map(|&e| e.into()).collect(),
        }
    }
}

impl AssignmentFixture {
    /// Transcript ids are interpreted in `space` (reference for collect-side payloads,
    /// models for assignments against `from_models` indices).
    pub fn to_core(&self, ids: &mut IdMaps, space: TxSpace) -> Assignment {
        Assignment {
            atype: self.atype,
            gene_atype: self.gene_atype,
            matches: self.matches.iter().map(|m| m.to_core(ids, space)).collect(),
        }
    }
    pub fn from_core(a: &Assignment, ids: &IdMaps) -> Self {
        AssignmentFixture {
            atype: a.atype,
            gene_atype: a.gene_atype,
            matches: a.matches.iter().map(|m| IsoformMatchFixture::from_core(m, ids)).collect(),
        }
    }
}

// ------------------------------------------------------------------------------------------
// ReadRecord
// ------------------------------------------------------------------------------------------

/// Fixture-only fields of a `ReadRecord`.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct ReadRecordExtras {
    pub assignment_id: i64,
    pub chr: String,
    pub cage_found: bool,
    pub introns_match: bool,
    pub additional_attributes: IndexMap<String, String>,
}

impl ReadRecordFixture {
    /// Core record (reference transcript space) plus the fixture-only extras.
    pub fn to_core(&self, ids: &mut IdMaps) -> (ReadRecord, ReadRecordExtras) {
        let rec = ReadRecord {
            idx: self.idx,
            name: self.name.clone(),
            exons: self.exons.clone(),
            corrected_exons: self.corrected_exons.clone(),
            boundary_deletions: self.boundary_deletions.clone(),
            polya: self.polya,
            polya_found: self.polya_found,
            strand: self.strand,
            mapped_strand: self.mapped_strand,
            mapq: self.mapq,
            multimapper: self.multimapper,
            group: ids.group_id(&self.read_group),
            assignment: self.assignment.to_core(ids, TxSpace::Ref),
            genomic_region: self.region,
        };
        let ex = ReadRecordExtras {
            assignment_id: self.assignment_id,
            chr: self.chr.clone(),
            cage_found: self.cage_found,
            introns_match: self.introns_match,
            additional_attributes: self.additional_attributes.clone(),
        };
        (rec, ex)
    }

    pub fn from_core(r: &ReadRecord, ids: &IdMaps, ex: &ReadRecordExtras) -> Self {
        ReadRecordFixture {
            idx: r.idx,
            assignment_id: ex.assignment_id,
            name: r.name.clone(),
            chr: ex.chr.clone(),
            region: r.genomic_region,
            exons: r.exons.clone(),
            corrected_exons: r.corrected_exons.clone(),
            boundary_deletions: r.boundary_deletions.clone(),
            polya: r.polya,
            polya_found: r.polya_found,
            cage_found: ex.cage_found,
            strand: r.strand,
            mapped_strand: r.mapped_strand,
            mapq: r.mapq,
            multimapper: r.multimapper,
            read_group: ids.group_names(r.group),
            assignment: AssignmentFixture::from_core(&r.assignment, ids),
            introns_match: ex.introns_match,
            additional_attributes: ex.additional_attributes.clone(),
        }
    }
}

// ------------------------------------------------------------------------------------------
// IsoformIndex
// ------------------------------------------------------------------------------------------

impl IsoformIndex {
    /// Builds an `IsoformIndex` from a Python `GeneInfo` dump (schema §2.1). Isoform order =
    /// `isoform_exons` dict order. `kind == "genes"|"region"` puts transcripts in the
    /// reference space, `"models"|"model"` in the model space (see [`IdMaps`]).
    pub fn from_dump(d: &IndexDump, ids: &mut IdMaps) -> Result<IsoformIndex, String> {
        let space = match d.kind.as_str() {
            "genes" | "region" => TxSpace::Ref,
            "models" | "model" => TxSpace::Model,
            k => return Err(format!("unknown index kind {k:?}")),
        };
        let chr = d.chr.as_deref().map(|c| ids.chroms.intern(c) as u16).unwrap_or(0);
        let introns_of: FxHashMap<&str, &Vec<Iv>> = d.isoform_introns.iter().map(|(k, v)| (k.as_str(), v)).collect();
        let strand_of: FxHashMap<&str, Strand> = d.isoform_strands.iter().map(|(k, v)| (k.as_str(), *v)).collect();
        let gene_of: FxHashMap<&str, &str> = d.gene_id_map.iter().map(|(k, v)| (k.as_str(), v.as_str())).collect();

        let mut idx = IsoformIndex { chr, start: d.start, end: d.end, delta: d.delta, ..Default::default() };
        let mut names: Vec<&str> = Vec::with_capacity(d.isoform_exons.len());
        for (tid, exons) in &d.isoform_exons {
            let introns: Vec<Iv> = match introns_of.get(tid.as_str()) {
                Some(v) => (*v).clone(),
                None => junctions_from_blocks(exons),
            };
            let gid = gene_of.get(tid.as_str()).ok_or_else(|| format!("transcript {tid} missing from gene_id_map"))?;
            let span = match (exons.first(), exons.last()) {
                (Some(a), Some(b)) => Iv::new(a.s, b.e),
                _ => Iv::default(),
            };
            idx.iso.push(IsoRec {
                key: ids.iso_ref(tid, space),
                gene: ids.gene_key(gid),
                strand: strand_of.get(tid.as_str()).copied().unwrap_or_default(),
                span,
                exon_off: idx.iso_exons.len() as u32,
                n_exon: exons.len() as u32,
                intron_off: idx.iso_introns.len() as u32,
                n_intron: introns.len() as u32,
            });
            idx.iso_exons.extend_from_slice(exons);
            idx.iso_introns.extend_from_slice(&introns);
            names.push(tid);
        }

        idx.introns = d.intron_features.clone();
        idx.exons = d.exon_features.clone();
        idx.split_exons = d.split_exon_features.clone();
        let n = names.len();
        idx.intron_prof = prof_matrix(&names, &d.intron_profiles, d.intron_features.len(), "intron")?;
        idx.exon_prof = prof_matrix(&names, &d.exon_profiles, d.exon_features.len(), "exon")?;
        idx.split_prof = prof_matrix(&names, &d.split_exon_profiles, d.split_exon_features.len(), "split_exon")?;
        let (ri, re, rs) = (ranges(&d.intron_ranges), ranges(&d.exon_ranges), ranges(&d.split_exon_ranges));
        idx.prof_range = names
            .iter()
            .map(|t| {
                let g = |m: &FxHashMap<&str, (u32, u32)>| m.get(t).copied().unwrap_or((0, 0));
                [g(&ri), g(&re), g(&rs)]
            })
            .collect();
        debug_assert_eq!(idx.prof_range.len(), n);

        // genes
        let gstrand: FxHashMap<&str, Strand> = d.gene_strands.iter().map(|(k, v)| (k.as_str(), *v)).collect();
        if d.kind == "genes" {
            let gregion: FxHashMap<&str, Iv> = d.gene_regions.iter().map(|(k, v)| (k.as_str(), *v)).collect();
            for g in &d.genes {
                idx.genes.push(GeneRec {
                    key: ids.gene_key(g),
                    strand: gstrand.get(g.as_str()).copied().unwrap_or_default(),
                    region: gregion.get(g.as_str()).copied().unwrap_or_default(),
                });
            }
        } else {
            for (g, s) in &d.gene_strands {
                let key = ids.gene_key(g);
                let mut region: Option<Iv> = None;
                for r in idx.iso.iter().filter(|r| r.gene == key) {
                    region = Some(match region {
                        None => r.span,
                        Some(x) => Iv::new(x.s.min(r.span.s), x.e.max(r.span.e)),
                    });
                }
                idx.genes.push(GeneRec { key, strand: *s, region: region.unwrap_or_default() });
            }
        }
        idx.finish_derived(&names);
        Ok(idx)
    }
}

fn ranges(r: &Pairs<String, (i32, i32)>) -> FxHashMap<&str, (u32, u32)> {
    r.iter().map(|(k, v)| (k.as_str(), (v.0.max(0) as u32, v.1.max(0) as u32))).collect()
}

fn prof_matrix(
    names: &[&str],
    profiles: &Pairs<String, Vec<i8>>,
    n_feat: usize,
    what: &str,
) -> Result<ProfMatrix, String> {
    let by: FxHashMap<&str, &Vec<i8>> = profiles.iter().map(|(k, v)| (k.as_str(), v)).collect();
    let mut m = ProfMatrix::new(names.len(), n_feat);
    for (i, t) in names.iter().enumerate() {
        if let Some(row) = by.get(t) {
            if row.len() != n_feat {
                return Err(format!("{what} profile of {t}: {} values for {n_feat} features", row.len()));
            }
            m.row_mut(i).copy_from_slice(row);
        }
    }
    Ok(m)
}

// ------------------------------------------------------------------------------------------
// Models and store snapshots
// ------------------------------------------------------------------------------------------

/// Fixture-only fields of a model.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct ModelExtras {
    pub chr: String,
    pub source: String,
    pub other_features: OtherFeatures,
    pub additional_info: Pairs<String, String>,
}

impl ModelFixture {
    /// Core model. `key` = model-space id; a known model's `kind` refers to its reference
    /// transcript (same id string).
    pub fn to_core(&self, ids: &mut IdMaps) -> (TranscriptModel, ModelExtras) {
        let key = ids.models.intern(&self.id);
        let kind = match self.ty {
            TranscriptModelType::Known => ModelKind::Known(ids.txs.intern(&self.id)),
            TranscriptModelType::NovelInCatalog => ModelKind::Nic,
            TranscriptModelType::NovelNotInCatalog => ModelKind::Nnic,
        };
        let m = TranscriptModel {
            key,
            strand: self.strand,
            gene: ids.gene_key(&self.gene),
            kind,
            exons: self.exons.clone(),
            intron_path: self.intron_path.clone().into_boxed_slice(),
            polya_confirmed: self.polya_confirmed,
            tss_confirmed: self.tss_confirmed,
            source_known: None,
        };
        let ex = ModelExtras {
            chr: self.chr.clone(),
            source: self.source.clone(),
            other_features: self.other_features.clone(),
            additional_info: self.additional_info.clone(),
        };
        (m, ex)
    }

    /// DTO from a core model; ids of models created in Rust that the maps do not know yet
    /// must be registered by the caller (`ids.models.intern(..)`) first.
    pub fn from_core(m: &TranscriptModel, ids: &IdMaps, ex: &ModelExtras) -> Self {
        let id = match m.kind {
            ModelKind::Known(t) => ids.txs.resolve(t).to_owned(),
            _ => ids.models.resolve(m.key).to_owned(),
        };
        ModelFixture {
            id,
            gene: ids.gene_name(m.gene).to_owned(),
            chr: ex.chr.clone(),
            strand: m.strand,
            ty: m.kind.model_type(),
            exons: m.exons.clone(),
            intron_path: m.intron_path.to_vec(),
            polya_confirmed: m.polya_confirmed,
            tss_confirmed: m.tss_confirmed,
            source: ex.source.clone(),
            other_features: ex.other_features.clone(),
            additional_info: ex.additional_info.clone(),
        }
    }
}

/// Fixture-only parts of a store snapshot.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct StoreExtras {
    /// Per model, in storage order.
    pub models: Vec<ModelExtras>,
    /// Run-wide `detected_known_isoforms` (sorted ids).
    pub detected_known_isoforms: Vec<String>,
    pub id_counter: i64,
}

impl StoreFixture {
    /// Core store for a block of `read_names.len()` reads. Per-read maps use the `*_idx`
    /// fields when present, otherwise expand each read name to every block read with that
    /// name (`read_names[idx]`).
    pub fn to_core(&self, ids: &mut IdMaps, read_names: &[&str]) -> Result<(ModelStore, StoreExtras), String> {
        let n = read_names.len();
        let mut st = ModelStore::new(n);
        let mut ex = StoreExtras {
            detected_known_isoforms: self.detected_known_isoforms.clone(),
            id_counter: self.id_counter,
            ..Default::default()
        };
        for m in &self.models {
            let (core, mex) = m.to_core(ids);
            st.models.push(core);
            ex.models.push(mex);
        }
        let check =
            |r: u32| if (r as usize) < n { Ok(r) } else { Err(format!("read_idx {r} out of range ({n} reads)")) };
        for (id, reads) in &self.transcript_read_ids {
            let k = ids.models.intern(id);
            st.transcript_read_ids.insert(k, reads.iter().map(|&r| check(r)).collect::<Result<_, _>>()?);
        }
        for (id, c) in &self.internal_counter {
            st.internal_counter.insert(ids.models.intern(id), *c as i32);
        }
        for (id, reads) in &self.nonfl_known_reads {
            st.nonfl_known_reads.insert(ids.txs.intern(id), reads.iter().map(|&r| check(r)).collect::<Result<_, _>>()?);
        }
        let mut by_name: FxHashMap<&str, Vec<usize>> = FxHashMap::default();
        if self.read_assignment_counts_idx.is_none() || self.construction_assignment_idx.is_none() {
            for (i, nm) in read_names.iter().enumerate() {
                by_name.entry(nm).or_default().push(i);
            }
        }
        match &self.read_assignment_counts_idx {
            Some(v) => {
                for &(r, c) in v {
                    st.read_assignment_counts[check(r)? as usize] = c as i32;
                }
            }
            None => {
                for (nm, c) in &self.read_assignment_counts {
                    for &i in by_name.get(nm.as_str()).map(Vec::as_slice).unwrap_or(&[]) {
                        st.read_assignment_counts[i] = *c as i32;
                    }
                }
            }
        }
        match &self.construction_assignment_idx {
            Some(v) => {
                for (r, id) in v {
                    st.construction_assignment[check(*r)? as usize] = Some(ids.models.intern(id));
                }
            }
            None => {
                for (nm, id) in &self.construction_assignment {
                    let k = ids.models.intern(id);
                    for &i in by_name.get(nm.as_str()).map(Vec::as_slice).unwrap_or(&[]) {
                        st.construction_assignment[i] = Some(k);
                    }
                }
            }
        }
        st.next_key = ids.models.len() as u32;
        Ok((st, ex))
    }

    /// DTO from a core store. The name-keyed per-read maps are emitted in ascending read
    /// order with non-zero / set entries only (4.1's dict order and defaultdict zero entries
    /// are not representable); compare them with [`store_diff`], which treats them as sets.
    pub fn from_core(st: &ModelStore, ids: &IdMaps, ex: &StoreExtras, read_names: &[&str]) -> Self {
        let empty = ModelExtras::default();
        let model_name = |k: u32| -> String {
            match st.models.iter().find(|m| m.key == k) {
                Some(TranscriptModel { kind: ModelKind::Known(t), .. }) => ids.txs.resolve(*t).to_owned(),
                _ => ids.models.resolve(k).to_owned(),
            }
        };
        let rac_idx: Pairs<u32, i64> = st
            .read_assignment_counts
            .iter()
            .enumerate()
            .filter(|(_, &c)| c != 0)
            .map(|(i, &c)| (i as u32, c as i64))
            .collect();
        let ca_idx: Pairs<u32, String> = st
            .construction_assignment
            .iter()
            .enumerate()
            .filter_map(|(i, m)| m.map(|k| (i as u32, model_name(k))))
            .collect();
        StoreFixture {
            models: st
                .models
                .iter()
                .enumerate()
                .map(|(i, m)| ModelFixture::from_core(m, ids, ex.models.get(i).unwrap_or(&empty)))
                .collect(),
            transcript_read_ids: st.transcript_read_ids.iter().map(|(k, v)| (model_name(*k), v.clone())).collect(),
            internal_counter: st.internal_counter.iter().map(|(k, v)| (model_name(*k), *v as i64)).collect(),
            read_assignment_counts: rac_idx.iter().map(|(i, c)| (read_names[*i as usize].to_owned(), *c)).collect(),
            construction_assignment: ca_idx
                .iter()
                .map(|(i, m)| (read_names[*i as usize].to_owned(), m.clone()))
                .collect(),
            read_assignment_counts_idx: Some(rac_idx),
            construction_assignment_idx: Some(ca_idx),
            nonfl_known_reads: st
                .nonfl_known_reads
                .iter()
                .map(|(k, v)| (ids.txs.resolve(*k).to_owned(), v.clone()))
                .collect(),
            detected_known_isoforms: ex.detected_known_isoforms.clone(),
            id_counter: ex.id_counter,
        }
    }
}

/// Store comparison for stage runners: `models`, `transcript_read_ids`,
/// `nonfl_known_reads`, `detected_known_isoforms` exactly (order included);
/// `internal_counter` exactly in order of keys but ignoring zero-valued entries; the
/// per-read maps as sets of non-zero (`*_idx`) entries. `id_counter` is ignored (IDs need
/// not match, D-M1).
pub fn store_diff(
    got: &StoreFixture,
    want: &StoreFixture,
    opts: &super::diff::DiffOptions,
) -> Option<super::diff::Diff> {
    #[derive(serde::Serialize)]
    struct View<'a> {
        models: &'a [ModelFixture],
        transcript_read_ids: &'a Pairs<String, Vec<u32>>,
        internal_counter: Vec<&'a (String, i64)>,
        read_assignment_counts_idx: Vec<(u32, i64)>,
        construction_assignment_idx: Vec<(u32, String)>,
        nonfl_known_reads: &'a Pairs<String, Vec<u32>>,
        detected_known_isoforms: &'a [String],
    }
    fn view(s: &StoreFixture) -> View<'_> {
        let mut rac: Vec<(u32, i64)> = s.read_assignment_counts_idx.clone().unwrap_or_default();
        rac.retain(|x| x.1 != 0);
        rac.sort();
        let mut ca = s.construction_assignment_idx.clone().unwrap_or_default();
        ca.sort();
        View {
            models: &s.models,
            transcript_read_ids: &s.transcript_read_ids,
            internal_counter: s.internal_counter.iter().filter(|x| x.1 != 0).collect(),
            read_assignment_counts_idx: rac,
            construction_assignment_idx: ca,
            nonfl_known_reads: &s.nonfl_known_reads,
            detected_known_isoforms: &s.detected_known_isoforms,
        }
    }
    super::diff::diff_values(&view(got), &view(want), opts)
}

// ------------------------------------------------------------------------------------------
// Graph and paths
// ------------------------------------------------------------------------------------------

impl From<&GraphFixture> for IntronGraph {
    fn from(g: &GraphFixture) -> Self {
        IntronGraph {
            clustered_introns: g.clustered_introns.iter().cloned().collect(),
            intron_correction_map: g.intron_correction_map.iter().cloned().collect(),
            discarded_introns: g.discarded_introns.clone(),
            outgoing: g.outgoing.iter().cloned().collect(),
            incoming: g.incoming.iter().cloned().collect(),
            starting_known_positions: g.starting_known_positions.iter().cloned().collect(),
            terminal_known_positions: g.terminal_known_positions.iter().cloned().collect(),
            max_coverage: g.max_coverage,
        }
    }
}
impl From<&IntronGraph> for GraphFixture {
    fn from(g: &IntronGraph) -> Self {
        GraphFixture {
            clustered_introns: g.clustered_introns.iter().map(|(k, v)| (*k, *v)).collect(),
            intron_correction_map: g.intron_correction_map.iter().map(|(k, v)| (*k, *v)).collect(),
            discarded_introns: g.discarded_introns.clone(),
            outgoing: pairs(&g.outgoing),
            incoming: pairs(&g.incoming),
            starting_known_positions: pairs(&g.starting_known_positions),
            terminal_known_positions: pairs(&g.terminal_known_positions),
            max_coverage: g.max_coverage,
        }
    }
}

fn pairs<K: Copy, V: Clone>(m: &IndexMap<K, V>) -> Vec<(K, V)> {
    m.iter().map(|(k, v)| (*k, v.clone())).collect()
}

impl U12FillOut {
    /// Core paths (without `known_isoforms_in_graph`, which comes from the
    /// `get_known_spliced_isoforms` case).
    pub fn to_core(&self) -> PathStorage {
        PathStorage {
            paths: self.paths.iter().map(|(p, c)| (p.clone().into_boxed_slice(), *c)).collect(),
            fl_paths: self.fl_paths.iter().map(|p| p.clone().into_boxed_slice()).collect(),
            paths_to_reads: self
                .paths_to_reads
                .iter()
                .map(|(p, r)| (p.clone().into_boxed_slice(), r.clone()))
                .collect(),
            known_isoforms_in_graph: IndexMap::new(),
        }
    }
    pub fn from_core(p: &PathStorage) -> Self {
        U12FillOut {
            paths: p.paths.iter().map(|(k, c)| (k.to_vec(), *c)).collect(),
            fl_paths: p.fl_paths.iter().map(|k| k.to_vec()).collect(),
            paths_to_reads: p.paths_to_reads.iter().map(|(k, r)| (k.to_vec(), r.clone())).collect(),
        }
    }
}
