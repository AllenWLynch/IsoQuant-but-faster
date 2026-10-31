//! # iq_collect — U8 (block collector, composes U4-U7) and U9 (multimapper resolution)
//!
//! Design §6.2, §5.2, §5.4-5.5. Entry points:
//! - [`collect_block`]: one block's BAM records → [`BlockStore`] (pure function).
//! - [`BlockCollector`]: the same, streaming: feed record batches with
//!   [`BlockCollector::push`], then [`BlockCollector::finish`] for the store and
//!   [`CollectStats`] (read and polyA counts for the §5.4 barrier). Per-read work inside a batch runs
//!   data-parallel on rayon; results are appended in BAM order, so the output does not
//!   depend on the thread count.
//! - U9: [`multimap::resolve`] and [`multimap::resolve_chr_duplicates`].
//!
//! Per-read flow (`alignment_processor.py:377-536`):
//! ```text
//! owner filter (fetch window) → flag filter (secondary/supplementary) → min_mapq →
//!   genic:      counter/subsample → U4 (exons, polyA, deletions; no exons → drop) → U5 → U6 →
//!               inconsistent & mapq < 5 → drop → aligned_polya_tail event → U7 correct →
//!               strand → antisense → ReadRecord
//!   intergenic: no exons → drop; ≤2 raw exons & mapq < 1 → drop; counter/subsample → U4 →
//!               intergenic assignment, corrected = exons → strand → ReadRecord
//! ```
//! The filters and the subsampling counter run sequentially before any per-read work, so
//! which reads are evaluated is known before the parallel part starts.
//!
//! Fixture runners: unit `U8` (fn `collect_block`) and `U9` (fn `resolve`), see
//! [`register_fixtures`].

pub mod fixture;
pub mod multimap;

use iq_assign::profile::ProfileBuilder;
use iq_assign::{assigner::assign_with, AssignParams, Scratch as ProfScratch};
use iq_read::correct::{correct_read_with, introns_of, CorrectParams};
use iq_read::polya::{ExternalPolyA, PolyAParams};
use iq_read::{assignment_strand, check_antisense, process_alignment_into};
use iq_types::read::flags;
use iq_types::{
    AlignedRead, Assignment, AssignmentType, Block, BlockStore, CombinedProfile, EventType, GroupPool,
    IsoformIndex, IsoformMatch, Iv, Iv0, MatchClass, MatchEvent, Matches, Params, ReadRecord, RecordLite, RefSlice,
    Strand, TagValue, NO_POLYA,
};
use rayon::prelude::*;
use std::sync::RwLock;
use std::time::Instant;

// ------------------------------------------------------------------------------------------
// Ownership (4.1 double processing; D-R2 retired)
// ------------------------------------------------------------------------------------------

/// Decides whether a record is processed by this block. Since D-R2 was retired
/// (`validation/deviation_effects.md`), both builds follow 4.1: a record is processed by
/// every block whose fetch returns it, and U9 keeps the best copy.
pub trait Owner: Sync {
    fn owns(&self, rec: &RecordLite<'_>) -> bool;
    /// True when the record is also processed by another block of the chromosome. The
    /// collector records the name hashes of such reads so that U9 can resolve their copies
    /// without scanning every read (a record collected twice is not a pass-A duplicate).
    fn shared(&self, _rec: &RecordLite<'_>) -> bool {
        false
    }
}

/// Every record given to the collector belongs to the block (the caller already did the
/// routing, or the block is the only one of its cluster).
#[derive(Copy, Clone, Debug, Default)]
pub struct AllRecords;

impl Owner for AllRecords {
    fn owns(&self, _: &RecordLite<'_>) -> bool {
        true
    }
}

/// Ownership among the blocks of a chromosome: `regions` are the block regions of the
/// cluster (or of the whole chromosome; sorted, disjoint), `me` the index of this block.
///
/// 4.1 semantics (`alignment_processor.py:160-167, 284-292`): a record belongs to every
/// block whose fetch window `fetch(region.s, region.e + 1)` returns it. A record that also
/// overlaps a neighbouring region is [`Owner::shared`] (regions are sorted and disjoint, so
/// a contiguous span that overlaps any other region overlaps an adjacent one).
#[derive(Copy, Clone, Debug)]
pub struct ClusterOwner<'a> {
    pub regions: &'a [Iv0],
    pub me: usize,
}

/// 0-based closed reference span of a record (`end0 - 1`, at least `pos0`).
#[inline]
pub fn record_span0(rec: &RecordLite<'_>) -> Iv0 {
    let e = rec.end0().max(rec.pos0 + 1) - 1;
    Iv0 { s: rec.pos0, e }
}

/// htslib `fetch(beg, end)` test: `pos < end && endpos > beg` (endpos ≥ pos + 1).
#[inline]
fn fetched_by(span: Iv0, r: Iv0) -> bool {
    span.s <= r.e && span.e >= r.s
}

impl Owner for ClusterOwner<'_> {
    fn owns(&self, rec: &RecordLite<'_>) -> bool {
        fetched_by(record_span0(rec), self.regions[self.me])
    }
    fn shared(&self, rec: &RecordLite<'_>) -> bool {
        let span = record_span0(rec);
        let r = self.regions[self.me];
        (span.s < r.s && self.me > 0 && fetched_by(span, self.regions[self.me - 1]))
            || (span.e > r.e && self.me + 1 < self.regions.len() && fetched_by(span, self.regions[self.me + 1]))
    }
}

// ------------------------------------------------------------------------------------------
// Read groups
// ------------------------------------------------------------------------------------------

/// Maps a record's read-group tag to a group-pool id (`ReadRecord::group`).
pub trait GroupIds: Sync {
    fn group_id(&self, tag: Option<&TagValue>) -> u32;
}

/// No `--read_group`: every read gets [`GroupPool::NONE`].
#[derive(Copy, Clone, Debug, Default)]
pub struct NoGroups;

impl GroupIds for NoGroups {
    fn group_id(&self, _: Option<&TagValue>) -> u32 {
        GroupPool::NONE
    }
}

/// A shared group pool for one tag grouper (`str(value)`, missing → `NA`; D-O2). Ids are
/// assigned in first-seen order, which depends on scheduling when blocks are collected in
/// parallel; the writer must order groups by name (or remap) for deterministic output.
#[derive(Debug, Default)]
pub struct SharedGroups(pub RwLock<GroupPool>);

impl SharedGroups {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn into_pool(self) -> GroupPool {
        self.0.into_inner().unwrap_or_else(|e| e.into_inner())
    }
}

impl GroupIds for SharedGroups {
    fn group_id(&self, tag: Option<&TagValue>) -> u32 {
        let g = iq_read::groups::group_string(tag);
        if let Some(id) = self.0.read().unwrap_or_else(|e| e.into_inner()).get(&g) {
            return id;
        }
        self.0.write().unwrap_or_else(|e| e.into_inner()).intern(&g)
    }
}

/// Any `Fn(Option<&TagValue>) -> u32` (used by the fixture runner).
pub struct FnGroups<F>(pub F);

impl<F: Fn(Option<&TagValue>) -> u32 + Sync> GroupIds for FnGroups<F> {
    fn group_id(&self, tag: Option<&TagValue>) -> u32 {
        (self.0)(tag)
    }
}

// ------------------------------------------------------------------------------------------
// Stats and timing
// ------------------------------------------------------------------------------------------

/// Per-block counters. `stored`/`polya_found` feed the sample-wide polyA fraction (§5.4).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CollectStats {
    /// Records given to the collector.
    pub records: u64,
    /// Not returned by this block's 4.1 fetch window.
    pub not_owned: u64,
    /// Supplementary / secondary.
    pub flag_filtered: u64,
    /// Below `--min_mapq`.
    pub min_mapq: u64,
    /// Skipped by the coverage-cap subsampling counter.
    pub subsampled: u64,
    /// No aligned exons.
    pub no_exons: u64,
    /// Intergenic `≤ 2 exons and mapq < simple_alignments_mapq_cutoff`.
    pub simple_dropped: u64,
    /// Genic, not unique/umd/ambiguous, and mapq < `inconsistent_mapq_cutoff`.
    pub inconsistent_mapq: u64,
    /// Reads stored.
    pub stored: u64,
    /// Stored reads with any polyA/T (`polyA_found`).
    pub polya_found: u64,
}

impl CollectStats {
    pub fn add(&mut self, o: &CollectStats) {
        self.records += o.records;
        self.not_owned += o.not_owned;
        self.flag_filtered += o.flag_filtered;
        self.min_mapq += o.min_mapq;
        self.subsampled += o.subsampled;
        self.no_exons += o.no_exons;
        self.simple_dropped += o.simple_dropped;
        self.inconsistent_mapq += o.inconsistent_mapq;
        self.stored += o.stored;
        self.polya_found += o.polya_found;
    }
}

/// Wall time per stage in nanoseconds (serial collection with timing enabled only).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StageTimes {
    /// Owner/flag/mapq filters, subsampling counter (and intergenic raw-exon count).
    pub select: u64,
    /// U4: CIGAR blocks, polyA/T, deletions.
    pub u4: u64,
    /// U5: profiles.
    pub u5: u64,
    /// U6: assignment.
    pub u6: u64,
    /// U7: correction, strand, antisense.
    pub u7: u64,
    /// ReadRecord fill + BlockStore push.
    pub store: u64,
}

impl StageTimes {
    pub fn total(&self) -> u64 {
        self.select + self.u4 + self.u5 + self.u6 + self.u7 + self.store
    }
    pub fn add(&mut self, o: &StageTimes) {
        self.select += o.select;
        self.u4 += o.u4;
        self.u5 += o.u5;
        self.u6 += o.u6;
        self.u7 += o.u7;
        self.store += o.store;
    }
}

struct Clock<'t> {
    t: Option<(&'t mut StageTimes, Instant)>,
}

impl<'t> Clock<'t> {
    #[inline]
    fn new(t: Option<&'t mut StageTimes>) -> Self {
        Clock { t: t.map(|t| (t, Instant::now())) }
    }
    #[inline]
    fn lap(&mut self, f: fn(&mut StageTimes) -> &mut u64) {
        if let Some((t, last)) = &mut self.t {
            let now = Instant::now();
            *f(t) += (now - *last).as_nanos() as u64;
            *last = now;
        }
    }
}

// ------------------------------------------------------------------------------------------
// Per-read evaluation
// ------------------------------------------------------------------------------------------

/// Per-thread buffers (reused across reads and blocks).
pub struct ReadScratch {
    ar: AlignedRead,
    prof: ProfScratch,
    profile: CombinedProfile,
    corrected: Vec<Iv>,
}

impl Default for ReadScratch {
    fn default() -> Self {
        ReadScratch {
            ar: AlignedRead {
                exons: Vec::new(),
                read_blocks: Vec::new(),
                cigar_blocks: Vec::new(),
                polya: NO_POLYA,
                exons_changed: false,
                deletions: Vec::new(),
            },
            prof: ProfScratch::new(),
            profile: CombinedProfile::default(),
            corrected: Vec::new(),
        }
    }
}

/// Result of evaluating one selected read.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Outcome {
    Kept,
    NoExons,
    InconsistentMapq,
}

fn blank_record() -> ReadRecord {
    ReadRecord {
        idx: 0,
        name: String::new(),
        exons: Vec::new(),
        corrected_exons: Vec::new(),
        boundary_deletions: Vec::new(),
        polya: NO_POLYA,
        polya_found: false,
        strand: Strand::Unknown,
        mapped_strand: Strand::Unknown,
        mapq: 0,
        multimapper: false,
        group: 0,
        assignment: Assignment::empty(AssignmentType::Noninformative),
        genomic_region: Iv::default(),
    }
}

/// Block-invariant state shared by all workers.
struct Ctx<'a> {
    index: &'a IsoformIndex,
    reference: &'a RefSlice<'a>,
    groups: &'a dyn GroupIds,
    genic: bool,
    region: Iv,
    prof: Option<ProfileBuilder<'a>>,
    ap: AssignParams,
    pp: PolyAParams,
    cp: CorrectParams,
    inconsistent_mapq_cutoff: u8,
}

impl Ctx<'_> {
    /// Evaluates one selected record into `out` (every field overwritten when `Kept`).
    fn process(&self, rec: &RecordLite<'_>, sc: &mut ReadScratch, out: &mut ReadRecord, clock: &mut Clock<'_>) -> Outcome {
        if !process_alignment_into(rec, &self.pp, ExternalPolyA::default(), &mut sc.ar) {
            clock.lap(|t| &mut t.u4);
            return Outcome::NoExons;
        }
        clock.lap(|t| &mut t.u4);
        let ar = &sc.ar;
        let mut a;
        if let Some(prof) = &self.prof {
            prof.build_into(&ar.exons, ar.polya, &mut sc.prof, &mut sc.profile);
            clock.lap(|t| &mut t.u5);
            a = assign_with(self.index, &self.ap, &sc.profile, ar.polya, &mut sc.prof.assign)
                .unwrap_or_else(|| Assignment::empty(AssignmentType::Noninformative));
            clock.lap(|t| &mut t.u6);
            if !matches!(
                a.atype,
                AssignmentType::Unique | AssignmentType::UniqueMinorDifference | AssignmentType::Ambiguous
            ) && rec.mapq < self.inconsistent_mapq_cutoff
            {
                return Outcome::InconsistentMapq;
            }
            if ar.exons_changed {
                for m in a.matches.iter_mut() {
                    m.add_event(MatchEvent::new(EventType::AlignedPolyaTail));
                }
            }
            sc.corrected = correct_read_with(&ar.exons, rec, &a, self.index, self.reference, &self.cp);
        } else {
            let mut ms = Matches::new();
            ms.push(IsoformMatch::new(MatchClass::Intergenic, None, None, [], Strand::Unknown, 0.0));
            a = Assignment::new(AssignmentType::Intergenic, ms);
            if ar.exons_changed {
                for m in a.matches.iter_mut() {
                    m.add_event(MatchEvent::new(EventType::AlignedPolyaTail));
                }
            }
            sc.corrected.clear();
            sc.corrected.extend_from_slice(&ar.exons);
        }
        let strand = assignment_strand(&a, ar.exons.len(), introns_of(&sc.corrected), ar.polya, self.reference);
        if self.genic {
            check_antisense(&mut a, strand);
        }
        clock.lap(|t| &mut t.u7);

        out.name.clear();
        match std::str::from_utf8(rec.qname) {
            Ok(s) => out.name.push_str(s),
            Err(_) => out.name.push_str(&String::from_utf8_lossy(rec.qname)),
        }
        out.exons.clear();
        out.exons.extend_from_slice(&ar.exons);
        std::mem::swap(&mut out.corrected_exons, &mut sc.corrected);
        out.boundary_deletions.clear();
        out.boundary_deletions.extend_from_slice(&ar.deletions);
        out.polya = ar.polya;
        out.polya_found = ar.polya.ext_a != -1 || ar.polya.ext_t != -1 || ar.polya.int_a != -1 || ar.polya.int_t != -1;
        out.strand = strand;
        out.mapped_strand = if rec.is_reverse() { Strand::Minus } else { Strand::Plus };
        out.mapq = rec.mapq;
        out.multimapper = rec.is_secondary();
        out.group = self.groups.group_id(rec.tag);
        out.assignment = a;
        out.genomic_region = self.region;
        Outcome::Kept
    }
}

// ------------------------------------------------------------------------------------------
// BlockCollector
// ------------------------------------------------------------------------------------------

/// Selected reads per batch below which the batch is evaluated serially.
pub const PAR_MIN_READS: usize = 64;
/// Reads per rayon task.
pub const PAR_CHUNK: usize = 32;

/// Streaming collector for one block (see the crate docs).
pub struct BlockCollector<'a> {
    ctx: Ctx<'a>,
    owner: &'a dyn Owner,
    params_min_mapq: Option<u8>,
    simple_mapq: u8,
    skip: u32,
    counter: u64,
    parallel: bool,
    store: BlockStore,
    stats: CollectStats,
    times: Option<StageTimes>,
    sel: Vec<u32>,
    /// Name hashes of selected records also processed by another block ([`Owner::shared`]).
    shared: Vec<u64>,
    sc: ReadScratch,
    rec: ReadRecord,
    raw_blocks: (Vec<Iv>, Vec<(i32, i32)>, Vec<(i32, i32)>),
}

impl<'a> BlockCollector<'a> {
    /// `index` is the block's reference index (`get_gene_info_for_region`); an index without
    /// exon features means an intergenic block (`GeneInfo.empty()`). `reference` must
    /// cover every record's span plus the U7 windows (±(12 + delta) bp); the whole
    /// chromosome is simplest.
    pub fn new(
        block: &Block,
        index: &'a IsoformIndex,
        reference: &'a RefSlice<'a>,
        params: &Params,
        owner: &'a dyn Owner,
        groups: &'a dyn GroupIds,
    ) -> Self {
        let genic = !index.exons.is_empty();
        let region = Iv { s: block.region.s as i32, e: block.region.e as i32 };
        let ctx = Ctx {
            index,
            reference,
            groups,
            genic,
            region,
            prof: genic.then(|| ProfileBuilder::new(index, &params.matching)),
            ap: AssignParams::new(params.matching.clone(), params.splice.clone(), false),
            pp: PolyAParams::from_params(params),
            cp: CorrectParams::from_params(params),
            inconsistent_mapq_cutoff: params.filter.inconsistent_mapq_cutoff,
        };
        BlockCollector {
            ctx,
            owner,
            params_min_mapq: params.filter.min_mapq,
            simple_mapq: params.filter.simple_alignments_mapq_cutoff,
            skip: block.skip.max(1),
            counter: 0,
            parallel: true,
            store: BlockStore::new(region),
            stats: CollectStats::default(),
            times: None,
            sel: Vec::new(),
            shared: Vec::new(),
            sc: ReadScratch::default(),
            rec: blank_record(),
            raw_blocks: Default::default(),
        }
    }

    /// Evaluate batches on rayon (default) or serially on the calling thread.
    pub fn parallel(mut self, on: bool) -> Self {
        self.parallel = on;
        self
    }

    /// Record per-stage wall times (forces serial evaluation).
    pub fn timed(mut self) -> Self {
        self.times = Some(StageTimes::default());
        self.parallel = false;
        self
    }

    /// True for a genic block.
    pub fn is_genic(&self) -> bool {
        self.ctx.genic
    }

    /// Sequential pre-pass: owner, flag and mapq filters and the subsampling counter
    /// (`alignment_processor.py:384-402, 457-468`). Fills `self.sel`.
    fn select(&mut self, recs: &[RecordLite<'_>]) {
        self.sel.clear();
        for (i, rec) in recs.iter().enumerate() {
            self.stats.records += 1;
            if !self.owner.owns(rec) {
                self.stats.not_owned += 1;
                continue;
            }
            // 4.1 tests `reference_id == -1`, which a region fetch never returns; placed
            // unmapped records go on and are dropped for having no exons.
            if rec.flag & (flags::SUPPLEMENTARY | flags::SECONDARY) != 0 {
                self.stats.flag_filtered += 1;
                continue;
            }
            if let Some(m) = self.params_min_mapq {
                if m > 0 && rec.mapq < m {
                    self.stats.min_mapq += 1;
                    continue;
                }
            }
            if !self.ctx.genic {
                let (ex, rb, cb) = &mut self.raw_blocks;
                iq_read::cigar::read_blocks_into(rec.pos0, rec.cigar, ex, rb, cb);
                if ex.is_empty() {
                    self.stats.no_exons += 1;
                    continue;
                }
                if ex.len() <= 2 && (rec.is_secondary() || rec.mapq < self.simple_mapq) {
                    self.stats.simple_dropped += 1;
                    continue;
                }
            }
            self.counter += 1;
            if self.skip > 1 && !self.counter.is_multiple_of(self.skip as u64) {
                self.stats.subsampled += 1;
                continue;
            }
            if self.owner.shared(rec) {
                self.shared.push(iq_types::name_hash(rec.qname));
            }
            self.sel.push(i as u32);
        }
    }

    fn push_record(store: &mut BlockStore, stats: &mut CollectStats, rec: &mut ReadRecord) {
        rec.idx = store.len() as u32;
        stats.stored += 1;
        stats.polya_found += rec.polya_found as u64;
        store.push(rec);
    }

    /// Processes one batch of records (BAM order; consecutive batches must continue the
    /// order).
    pub fn push(&mut self, recs: &[RecordLite<'_>]) {
        let t0 = self.times.is_some().then(Instant::now);
        self.select(recs);
        if let (Some(t), Some(t0)) = (&mut self.times, t0) {
            t.select += t0.elapsed().as_nanos() as u64;
        }
        if self.parallel && self.sel.len() >= PAR_MIN_READS && rayon::current_num_threads() > 1 {
            self.push_parallel(recs);
            return;
        }
        let sel = std::mem::take(&mut self.sel);
        for &i in &sel {
            let mut clock = Clock::new(self.times.as_mut());
            match self.ctx.process(&recs[i as usize], &mut self.sc, &mut self.rec, &mut clock) {
                Outcome::Kept => {
                    Self::push_record(&mut self.store, &mut self.stats, &mut self.rec);
                    clock.lap(|t| &mut t.store);
                }
                Outcome::NoExons => self.stats.no_exons += 1,
                Outcome::InconsistentMapq => self.stats.inconsistent_mapq += 1,
            }
        }
        self.sel = sel;
    }

    fn push_parallel(&mut self, recs: &[RecordLite<'_>]) {
        let ctx = &self.ctx;
        // quirk switches are thread-local (quirks41 fixture runs): carry them into workers
        let quirks = iq_types::current_quirks();
        let results: Vec<(Vec<ReadRecord>, u64, u64)> = self
            .sel
            .par_chunks(PAR_CHUNK)
            .map_init(ReadScratch::default, |sc, chunk| iq_types::with_quirks(quirks, || {
                let mut kept = Vec::with_capacity(chunk.len());
                let (mut no_ex, mut mq) = (0u64, 0u64);
                let mut clock = Clock::new(None);
                for &i in chunk {
                    let mut r = blank_record();
                    match ctx.process(&recs[i as usize], sc, &mut r, &mut clock) {
                        Outcome::Kept => kept.push(r),
                        Outcome::NoExons => no_ex += 1,
                        Outcome::InconsistentMapq => mq += 1,
                    }
                }
                (kept, no_ex, mq)
            }))
            .collect();
        for (kept, no_ex, mq) in results {
            self.stats.no_exons += no_ex;
            self.stats.inconsistent_mapq += mq;
            for mut r in kept {
                Self::push_record(&mut self.store, &mut self.stats, &mut r);
            }
        }
    }

    /// Stage times so far (`None` unless [`Self::timed`]).
    pub fn times(&self) -> Option<StageTimes> {
        self.times
    }

    pub fn stats(&self) -> &CollectStats {
        &self.stats
    }

    /// The block's store and counters.
    pub fn finish(self) -> (BlockStore, CollectStats) {
        (self.store, self.stats)
    }

    /// [`Self::finish`] plus the name hashes of selected reads that another block also
    /// processes (unsorted; a superset of the stored ones). Feed their union with the pass-A
    /// duplicate candidates to [`multimap::resolve_chr_duplicates`].
    pub fn finish_shared(self) -> (BlockStore, CollectStats, Vec<u64>) {
        (self.store, self.stats, self.shared)
    }
}

/// Records per internal batch of [`collect_block`].
pub const COLLECT_BATCH: usize = 4096;

/// `fn collect_block(&Block, records, &IsoformIndex, &RefSlice, &Params, owner) -> BlockStore`
/// (U8): `records` are the block's fetched records in BAM order. Read groups are not
/// tagged (use [`BlockCollector`] with a [`GroupIds`] for that).
pub fn collect_block<'r>(
    block: &Block,
    records: impl IntoIterator<Item = RecordLite<'r>>,
    index: &IsoformIndex,
    reference: &RefSlice<'_>,
    params: &Params,
    owner: &dyn Owner,
) -> BlockStore {
    collect_block_with(block, records, index, reference, params, owner, &NoGroups).0
}

/// [`collect_block`] with read-group tagging, returning the stats too.
pub fn collect_block_with<'r>(
    block: &Block,
    records: impl IntoIterator<Item = RecordLite<'r>>,
    index: &IsoformIndex,
    reference: &RefSlice<'_>,
    params: &Params,
    owner: &dyn Owner,
    groups: &dyn GroupIds,
) -> (BlockStore, CollectStats) {
    let reference = *reference;
    let mut c = BlockCollector::new(block, index, &reference, params, owner, groups);
    let mut buf: Vec<RecordLite<'r>> = Vec::with_capacity(COLLECT_BATCH);
    for r in records {
        buf.push(r);
        if buf.len() == COLLECT_BATCH {
            c.push(&buf);
            buf.clear();
        }
    }
    if !buf.is_empty() {
        c.push(&buf);
    }
    c.finish()
}

/// Registers this crate's fixture runners.
pub fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    fixture::register_fixtures(reg);
}
