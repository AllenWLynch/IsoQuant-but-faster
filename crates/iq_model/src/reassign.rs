//! U17: construction read (re)assignment and read ownership
//! (`ConstructionReadAssigner`, `model_construction/read_assigner.py:33-154`;
//! `reference_unique_known_isoform`, `context.py:36-44`; design §7.8).
//!
//! [`assign_reads_to_models`] runs twice per block (after the pre-filter and after the
//! filter); [`resolve_ownership`] runs once after the second pass.
//!
//! Structure of `assign_reads_to_models` (design §7.8, §5.5 rule 2):
//! 1. `IsoformIndex::from_models(models)` is built once per call.
//! 2. Reads go in chunks of [`CHUNK`]. For every read of a chunk that will reach the
//!    profiling branch (not yet bound, not gated to a reference-unique known isoform), the
//!    profile and the `quick_mode = true` assignment are computed (in parallel with rayon
//!    when there are at least [`PAR_THRESHOLD`] of them).
//! 3. The 4.1 loop then runs serially in read order over the chunk, consuming the
//!    precomputed assignments.
//!
//! Exactness: the loop changes only the current read's `read_assignment_counts` slot (reads
//! are keyed by block index, D-M5), so whether a read reaches the profiling branch is known
//! before the loop, and the assignment is a pure function of the read and the index. A read
//! listed twice in `reads` (never in production: `reads` = all block reads) whose second
//! visit profiles again is recomputed serially.
//!
//! Order: `transcript_read_ids` gains keys in the order of the assignment's matches, which
//! 4.1 sorts by transcript-id string (`id_rank` of the index). Production builds the index
//! with provisional ids (creation order); the fixture runner passes the real ids. This can
//! only change key order, never content.
//!
//! Defaultdict reads that insert keys in 4.1 (`transcript_read_ids[id]`,
//! `internal_counter[id]` in `resolve_read_ownership_conflicts`) are mirrored, since they
//! fix later key order in the store.

use crate::ctx::Ctx;
use iq_assign::profile::{ProfileBuilder, Scratch};
use iq_types::{Assignment, IsoRef, IsoformIndex, ModelIdx, ModelKind, ModelStore, ReadIdx, ReadView, TxIdx};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

/// Minimum number of reads to profile in a chunk before the work goes to rayon.
pub const PAR_THRESHOLD: usize = 256;
/// Reads per chunk (bounds the precomputed assignments held at once).
pub const CHUNK: usize = 1 << 15;

/// `reference_unique_known_isoform` (`context.py:36-44`): the isoform of a unique
/// reference-step assignment (`isoform_matches[0].assigned_transcript`), or `None`.
pub fn reference_unique_known_isoform(r: &ReadView<'_>) -> Option<IsoRef> {
    if r.atype().is_unique() {
        r.first_match().and_then(|m| m.tx())
    } else {
        None
    }
}

/// `assign_reads_to_models(read_assignment_storage)` over all block reads.
pub fn assign_reads_to_models(ms: &mut ModelStore, ctx: &Ctx<'_>) {
    let reads: Vec<ReadIdx> = (0..ctx.store.len() as ReadIdx).collect();
    assign_reads(ms, ctx, &reads);
}

/// `assign_reads_to_models(reads)` (the 4.1 argument, in order). Builds
/// `IsoformIndex::from_models(models)` with provisional ids.
pub fn assign_reads(ms: &mut ModelStore, ctx: &Ctx<'_>, reads: &[ReadIdx]) {
    if ms.models.is_empty() {
        assign_reads_with_index(ms, ctx, reads, &IsoformIndex::default(), PAR_THRESHOLD);
        return;
    }
    let index = iq_annot::index_from_models(ctx.index.chr, &ms.models, ctx.params.matching.delta);
    assign_reads_with_index(ms, ctx, reads, &index, PAR_THRESHOLD);
}

/// What the serial loop does with a read, decided before it.
enum Plan {
    /// `read_assignment_counts > 0` at the start: skipped.
    Bound,
    /// Gated to a reference-unique known isoform: its model, or `None` if not built.
    Gate(Option<ModelIdx>),
    /// Profiled; slot in the chunk's assignment vector.
    Profile(usize),
}

/// [`assign_reads`] against a prebuilt `index` (= `from_models(ms.models)`; unused when
/// there are no models). Parallel above `par_threshold` reads to profile per chunk
/// (`0` = always, `usize::MAX` = never); the result does not depend on it.
pub fn assign_reads_with_index(ms: &mut ModelStore, ctx: &Ctx<'_>, reads: &[ReadIdx], index: &IsoformIndex, par_threshold: usize) {
    if ms.models.is_empty() {
        for &r in reads {
            ms.read_assignment_counts[r as usize] = 0;
        }
        return;
    }
    // model_ids: a known reference isoform keeps its id, so the known models are the
    // surviving known isoforms (dict semantics: the last model with an id wins, as in
    // from_models).
    let known_model: HashMap<TxIdx, ModelIdx> = ms
        .models
        .iter()
        .filter_map(|m| match m.kind {
            ModelKind::Known(t) => Some((t, m.key)),
            _ => None,
        })
        .collect();
    let builder = ProfileBuilder::new(index, &ctx.params.matching);
    let assign_one = |r: ReadIdx, scratch: &mut Scratch| -> Option<Assignment> {
        let v = ctx.store.view(r as usize);
        let exons = v.corrected_exons();
        let prof = builder.build(&exons, v.polya(), scratch);
        ctx.assigner.assign(index, &prof, true, r as u64)
    };

    let mut plan: Vec<Plan> = Vec::new();
    let mut todo: Vec<ReadIdx> = Vec::new();
    let mut scratch = Scratch::new();
    for chunk in reads.chunks(CHUNK) {
        plan.clear();
        todo.clear();
        for &r in chunk {
            let p = if ms.read_assignment_counts[r as usize] > 0 {
                Plan::Bound
            } else {
                match reference_unique_known_isoform(&ctx.store.view(r as usize)) {
                    Some(IsoRef::Ref(t)) => Plan::Gate(known_model.get(&t).copied()),
                    // a block-index assignment never names a model
                    Some(IsoRef::Model(_)) => Plan::Gate(None),
                    None => {
                        todo.push(r);
                        Plan::Profile(todo.len() - 1)
                    }
                }
            };
            plan.push(p);
        }
        let results: Vec<Option<Assignment>> = if todo.len() >= par_threshold {
            todo.par_iter().map_init(Scratch::new, |s, &r| assign_one(r, s)).collect()
        } else {
            todo.iter().map(|&r| assign_one(r, &mut scratch)).collect()
        };
        let mut results: Vec<Option<Option<Assignment>>> = results.into_iter().map(Some).collect();

        for (&r, p) in chunk.iter().zip(&plan) {
            let ri = r as usize;
            // re-checked: a read listed twice may have been bound by its first visit
            if ms.read_assignment_counts[ri] > 0 {
                continue;
            }
            match *p {
                Plan::Bound => {
                    // count was > 0 and only this read's own visits change it
                    unreachable!("bound read lost its count")
                }
                Plan::Gate(Some(k)) => {
                    ms.transcript_read_ids.entry(k).or_default().push(r);
                    *ms.internal_counter.entry(k).or_default() += 1;
                    ms.read_assignment_counts[ri] += 1;
                }
                Plan::Gate(None) => ms.read_assignment_counts[ri] = 0,
                Plan::Profile(slot) => {
                    let a = match results[slot].take() {
                        Some(a) => a,
                        None => assign_one(r, &mut scratch), // duplicate visit
                    };
                    apply_assignment(ms, r, a.as_ref());
                }
            }
        }
    }
}

/// The consistent / inconsistent branch of the 4.1 loop (`read_assigner.py:137-152`).
fn apply_assignment(ms: &mut ModelStore, r: ReadIdx, a: Option<&Assignment>) {
    let ri = r as usize;
    let Some(a) = a.filter(|a| a.atype.is_consistent()) else {
        ms.read_assignment_counts[ri] = 0;
        return;
    };
    let model_of = |m: &iq_types::IsoformMatch| match m.tx {
        Some(IsoRef::Model(k)) => Some(k),
        _ => None,
    };
    if a.matches.len() == 1 {
        if let Some(k) = model_of(&a.matches[0]) {
            *ms.internal_counter.entry(k).or_default() += 1;
        }
    }
    for m in a.matches.iter() {
        ms.read_assignment_counts[ri] += 1;
        // a consistent model assignment always names a model (4.1 would key `None`)
        if let Some(k) = model_of(m) {
            ms.transcript_read_ids.entry(k).or_default().push(r);
        }
    }
}

/// `resolve_read_ownership_conflicts` (`read_assigner.py:33-92`).
pub fn resolve_ownership(ms: &mut ModelStore) {
    if ms.nonfl_known_reads.is_empty() {
        return;
    }
    // reads currently owned by each model
    let mut owned: HashMap<ModelIdx, u32> = HashMap::new();
    for m in ms.construction_assignment.iter().flatten() {
        *owned.entry(*m).or_default() += 1;
    }
    let is_nonfl_known = |kind: ModelKind| matches!(kind, ModelKind::Known(t) if ms.nonfl_known_reads.contains_key(&t));
    // read -> a surviving winner model that should reclaim it (setdefault: first wins)
    let mut read_to_fl: HashMap<ReadIdx, ModelIdx> = HashMap::new();
    for i in 0..ms.models.len() {
        let (key, kind) = (ms.models[i].key, ms.models[i].kind);
        if is_nonfl_known(kind) {
            continue; // assignment-based known: a potential loser, never a winner
        }
        if matches!(kind, ModelKind::Known(_)) && owned.get(&key).copied().unwrap_or(0) > 0 {
            continue; // FL-path known that still owns reads: co-expressed
        }
        for &r in ms.transcript_read_ids.entry(key).or_default().iter() {
            read_to_fl.entry(r).or_insert(key);
        }
    }

    // surviving known models by transcript id
    let surviving: HashMap<TxIdx, ModelIdx> = ms
        .models
        .iter()
        .filter_map(|m| match m.kind {
            ModelKind::Known(t) => Some((t, m.key)),
            _ => None,
        })
        .collect();
    let mut emptied: Vec<(TxIdx, ModelIdx)> = Vec::new();
    let nonfl = std::mem::take(&mut ms.nonfl_known_reads);
    for (&tx, reads) in &nonfl {
        let Some(&known) = surviving.get(&tx) else { continue };
        for &r in reads {
            let ri = r as usize;
            if ms.construction_assignment[ri] != Some(known) {
                continue; // no longer owned by this known (e.g. moved to an alt-end NIC)
            }
            let Some(&fl) = read_to_fl.get(&r) else { continue };
            if fl == known {
                continue;
            }
            ms.construction_assignment[ri] = Some(fl);
            let list = ms.transcript_read_ids.entry(known).or_default();
            if let Some(pos) = list.iter().position(|&x| x == r) {
                list.remove(pos);
                *ms.internal_counter.entry(known).or_default() -= 1;
                ms.read_assignment_counts[ri] -= 1;
            }
        }
        if *ms.internal_counter.entry(known).or_default() <= 0 {
            emptied.push((tx, known));
        }
    }
    ms.nonfl_known_reads = nonfl;

    if !emptied.is_empty() {
        // 4.1 iterates a set here; deleting is order-independent in content
        for &(_, k) in &emptied {
            ms.delete_from_storage(k);
        }
        let drop: HashSet<TxIdx> = emptied.iter().map(|&(t, _)| t).collect();
        ms.models.retain(|m| !matches!(m.kind, ModelKind::Known(t) if drop.contains(&t)));
    }
}

/// Registers the U17 fixture runners.
pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U17Assign>().add::<fixture::U17Resolve>();
}

pub(crate) mod fixture {
    //! U17 runners (plumbing: [`crate::fixture_util`]). Nested assignments run live
    //! ([`crate::fl::LiveAssigner`], U6 exact). The assign runner runs both the serial and
    //! the forced-parallel path and fails the case if they differ.
    use super::*;
    use crate::fixture_util::{compare_stores, StageInput};
    use crate::fl::LiveAssigner;
    use iq_assign::AssignParams;
    use iq_types::fixture::schema::{StageStoreIn, StoreOut};
    use iq_types::fixture::{FixError, FixResult, FixtureCtx, UnitFixture, Verdict};
    use iq_types::{DevId, Predictions};

    pub struct U17Assign;

    impl UnitFixture for U17Assign {
        const UNIT: &'static str = "U17";
        const FN: &'static str = "assign_reads_to_models";
        const DEVIATIONS: &'static [DevId] = &[DevId::M5];
        type In = StageStoreIn;
        type Out = StoreOut;

        fn run(input: &StageStoreIn, ctx: &FixtureCtx) -> FixResult<StoreOut> {
            let si = StageInput::load(input, ctx)?;
            let p = &si.params;
            let assigner = LiveAssigner { params: AssignParams::new(p.matching.clone(), p.splice.clone(), false) };
            let preds = Predictions::default();
            let c = Ctx {
                params: p,
                index: &si.index,
                store: &si.block.store,
                graph: &si.graph,
                paths: &si.paths,
                predictions: &preds,
                assigner: &assigner,
            };
            let reads: Vec<ReadIdx> = match &input.reads {
                Some(r) => r.clone(),
                None => (0..si.block.store.len() as ReadIdx).collect(),
            };
            // from_models with the real transcript ids (4.1 sorts matches by id string)
            let names: Vec<&str> = si.store.models.iter().map(|m| si.ids.models.resolve(m.key)).collect();
            let index = iq_annot::index_from_models_named(si.index.chr, &si.store.models, p.matching.delta, &names);
            let mut serial = si.store.clone();
            assign_reads_with_index(&mut serial, &c, &reads, &index, usize::MAX);
            let mut par = si.store.clone();
            assign_reads_with_index(&mut par, &c, &reads, &index, 0);
            // compare snapshots (IndexMap equality ignores key order)
            let out = si.store_out(&serial, &si.detected_known);
            if out != si.store_out(&par, &si.detected_known) {
                return Err(FixError("U17: serial and parallel reassignment differ".into()));
            }
            Ok(StoreOut { store: out })
        }

        fn compare(got: &StoreOut, want: &StoreOut, ctx: &FixtureCtx) -> Verdict {
            compare_stores(&got.store, &want.store, ctx)
        }
    }

    pub struct U17Resolve;

    impl UnitFixture for U17Resolve {
        const UNIT: &'static str = "U17";
        const FN: &'static str = "resolve_read_ownership_conflicts";
        const DEVIATIONS: &'static [DevId] = &[DevId::M5];
        type In = StageStoreIn;
        type Out = StoreOut;

        fn run(input: &StageStoreIn, ctx: &FixtureCtx) -> FixResult<StoreOut> {
            let si = StageInput::load(input, ctx)?;
            let mut ms = si.store.clone();
            resolve_ownership(&mut ms);
            Ok(StoreOut { store: si.store_out(&ms, &si.detected_known) })
        }

        fn compare(got: &StoreOut, want: &StoreOut, ctx: &FixtureCtx) -> Verdict {
            compare_stores(&got.store, &want.store, ctx)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use iq_types::{GeneKey, Iv, Strand, TranscriptModel};

    fn model(key: ModelIdx, kind: ModelKind) -> TranscriptModel {
        TranscriptModel {
            key,
            strand: Strand::Plus,
            gene: GeneKey::Novel(0),
            kind,
            exons: vec![Iv::new(100, 200), Iv::new(300, 400)],
            intron_path: Default::default(),
            polya_confirmed: false,
            tss_confirmed: false,
            source_known: None,
        }
    }

    fn store(kinds: &[ModelKind], n_reads: usize) -> ModelStore {
        let mut ms = ModelStore::new(n_reads);
        for &k in kinds {
            let key = ms.new_key();
            ms.models.push(model(key, k));
        }
        ms
    }

    fn keys<V>(m: &iq_types::indexmap::IndexMap<ModelIdx, V>) -> Vec<ModelIdx> {
        m.keys().copied().collect()
    }

    #[test]
    fn fl_winners_reclaim_reads() {
        // 0 novel FL, 1 phantom FL known (owns nothing), 2 non-FL known, 3 FL known owning reads
        let mut ms = store(&[ModelKind::Nic, ModelKind::Known(10), ModelKind::Known(20), ModelKind::Known(30)], 6);
        for (r, m) in [(1, 0), (2, 1), (1, 1), (0, 2), (1, 2), (2, 2), (3, 2), (5, 3), (0, 3)] {
            ms.save_assigned_read(r, m);
        }
        ms.nonfl_known_reads.insert(20, vec![0, 1, 2, 3]);
        resolve_ownership(&mut ms);
        // read 0: owned by model 3 now -> untouched; read 1 -> novel 0; read 2 -> phantom 1; read 3 no winner
        assert_eq!(ms.construction_assignment, vec![Some(3), Some(0), Some(1), Some(2), None, Some(3)]);
        assert_eq!(ms.transcript_read_ids[&2], vec![0, 3]);
        assert_eq!(ms.internal_counter[&2], 2);
        assert_eq!(ms.read_assignment_counts, vec![2, 2, 1, 1, 0, 1]);
        // winners' own lists are not touched
        assert_eq!(ms.transcript_read_ids[&0], vec![1]);
        assert_eq!(ms.transcript_read_ids[&1], vec![2, 1]);
        assert_eq!(ms.models.len(), 4);
    }

    #[test]
    fn emptied_known_is_dropped() {
        let mut ms = store(&[ModelKind::Nic, ModelKind::Known(20)], 1);
        ms.save_assigned_read(0, 0);
        ms.save_assigned_read(0, 1);
        ms.nonfl_known_reads.insert(20, vec![0]);
        resolve_ownership(&mut ms);
        assert_eq!(ms.models.iter().map(|m| m.key).collect::<Vec<_>>(), vec![0]);
        assert_eq!(keys(&ms.transcript_read_ids), vec![0]);
        assert_eq!(keys(&ms.internal_counter), vec![0]);
        assert_eq!(ms.read_assignment_counts, vec![1]);
        assert_eq!(ms.construction_assignment, vec![Some(0)]);
        // nonfl_known_reads itself is left as is
        assert_eq!(ms.nonfl_known_reads[&20], vec![0]);
    }

    #[test]
    fn defaultdict_reads_insert_keys() {
        // model 1 is a winner candidate with no entry: 4.1 `transcript_read_ids[id]` inserts []
        let mut ms = store(&[ModelKind::Known(20), ModelKind::Nic], 1);
        ms.save_assigned_read(0, 0);
        ms.nonfl_known_reads.insert(20, vec![0]);
        resolve_ownership(&mut ms);
        assert_eq!(keys(&ms.transcript_read_ids), vec![0, 1]);
        assert!(ms.transcript_read_ids[&1].is_empty());
        assert_eq!(ms.construction_assignment, vec![Some(0)]);
        assert_eq!(ms.models.len(), 2);
    }

    #[test]
    fn no_nonfl_knowns_is_a_no_op() {
        let mut ms = store(&[ModelKind::Nic], 1);
        let before = ms.clone();
        resolve_ownership(&mut ms);
        assert_eq!(keys(&ms.transcript_read_ids), keys(&before.transcript_read_ids));
        assert_eq!(ms, before);
    }

    #[test]
    fn known_not_surviving_is_skipped() {
        // non-FL known 20 has no model any more (filtered): its reads stay where they are
        let mut ms = store(&[ModelKind::Nic], 1);
        ms.save_assigned_read(0, 0);
        ms.nonfl_known_reads.insert(20, vec![0]);
        let before = ms.clone();
        resolve_ownership(&mut ms);
        assert_eq!(ms, before);
    }
}
