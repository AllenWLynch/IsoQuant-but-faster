//! U12: `IntronPathStorage.fill`, `IntronPathProcessor.thread_*` (`intron_path.py:17-163`),
//! `get_known_spliced_isoforms` (`context.py:202-224`).
//!
//! Design §7.3. Differences from 4.1 that keep the outputs identical:
//! - [`PathProcessor::new`] precomputes, per intron vertex and side, everything
//!   `thread_ends`/`thread_starts` derive from `get_outgoing`/`get_incoming` on every call
//!   (sorted polyA/TSS position lists, the outermost neighbouring exon boundary, and the last
//!   two entries of the stable-sorted terminal list). 4.1 rebuilds and sorts these lists up
//!   to six times per read.
//! - [`fill_paths`] threads the introns of each distinct interned chain once
//!   (`thread_introns` only depends on the chain); terminals are then looked up per read.
//!   Reads are keyed by `(intron-path id, start vertex, end vertex)` so no path vector is
//!   built per read; the 4.1 dicts are materialized at the end in first-insertion order.

use iq_types::indexmap::IndexMap;
use iq_types::model::vertex_code as vc;
use iq_types::{BlockStore, ChainTable, IntronGraph, IsoRef, IsoformIndex, Iv, Params, Path, PathStorage, Strand, Vertex};
use std::collections::{HashMap, HashSet};

/// Terminal-side summary of one intron vertex (one of outgoing / incoming).
#[derive(Clone, Debug, Default)]
struct Side {
    /// Sorted positions of the confirmed polyA (right) / polyT (left) vertices.
    tail: Vec<i64>,
    /// Sorted positions of the TSS vertices (`tss_right` / `tss_left`).
    tss: Vec<i64>,
    /// Right side: `max(outgoing intron starts)`; left side: `min(incoming intron ends)`.
    intron_bound: Option<i64>,
    /// Right side: last element of the stable-sorted `read_end + tss_right + polya` list;
    /// left side: first element of `read_start + tss_left + polyt`.
    extreme: Option<Vertex>,
    /// Position of the element next to `extreme` in that sorted list (`[-2]` / `[1]`).
    next_pos: Option<i64>,
}

impl Side {
    /// First (in sorted order) of `positions` with `|p - x| <= d`.
    #[inline]
    fn first_within(positions: &[i64], x: i64, d: i64) -> Option<i64> {
        let k = positions.partition_point(|&p| p < x - d);
        positions.get(k).copied().filter(|&p| p <= x + d)
    }
}

/// `IntronPathProcessor` (`intron_path.py:66-163`) over a finished graph.
pub struct PathProcessor<'g> {
    graph: &'g IntronGraph,
    discarded: HashSet<Iv>,
    right: HashMap<Vertex, Side>,
    left: HashMap<Vertex, Side>,
    apa_delta: i64,
    delta: i64,
}

/// Builds the summary of one adjacency list. `right` selects the outgoing (read-end) side.
fn summarize(targets: &[Vertex], right: bool) -> Side {
    let (tail_c, bare_c, tss_c) =
        if right { (vc::POLYA, vc::READ_END, vc::TSS_RIGHT) } else { (vc::POLYT, vc::READ_START, vc::TSS_LEFT) };
    let mut s = Side::default();
    let mut bare: Vec<i64> = Vec::new();
    for v in targets {
        match v.0 {
            c if c >= 0 => {
                let b = if right { v.0 } else { v.1 };
                s.intron_bound = Some(match s.intron_bound {
                    None => b,
                    Some(x) if right => x.max(b),
                    Some(x) => x.min(b),
                });
            }
            c if c == tail_c => s.tail.push(v.1),
            c if c == bare_c => bare.push(v.1),
            c if c == tss_c => s.tss.push(v.1),
            _ => {}
        }
    }
    s.tail.sort_unstable();
    s.tss.sort_unstable();
    bare.sort_unstable();
    // 4.1: sorted(bare + tss + tail, key=pos), a stable sort. Ties at the extreme position are
    // resolved by concatenation order: last of the list on the right, first on the left.
    let mut all: Vec<Vertex> = Vec::with_capacity(bare.len() + s.tss.len() + s.tail.len());
    all.extend(bare.iter().map(|&p| Vertex(bare_c, p)));
    all.extend(s.tss.iter().map(|&p| Vertex(tss_c, p)));
    all.extend(s.tail.iter().map(|&p| Vertex(tail_c, p)));
    all.sort_by_key(|v| v.1);
    if right {
        s.extreme = all.last().copied();
        s.next_pos = all.len().checked_sub(2).map(|k| all[k].1);
    } else {
        s.extreme = all.first().copied();
        s.next_pos = all.get(1).map(|v| v.1);
    }
    s
}

impl<'g> PathProcessor<'g> {
    /// Precomputes per-vertex terminal summaries (once per graph).
    pub fn new(graph: &'g IntronGraph, params: &Params) -> Self {
        let right = graph.outgoing.iter().filter(|(k, _)| !k.is_terminal()).map(|(k, t)| (*k, summarize(t, true))).collect();
        let left = graph.incoming.iter().filter(|(k, _)| !k.is_terminal()).map(|(k, t)| (*k, summarize(t, false))).collect();
        PathProcessor {
            graph,
            discarded: graph.discarded_introns.iter().copied().collect(),
            right,
            left,
            apa_delta: params.matching.apa_delta as i64,
            delta: params.matching.delta as i64,
        }
    }

    /// `thread_introns` (`:87-93`): `None` if an intron is discarded, else the substituted
    /// introns, appended to `out` (cleared first).
    pub fn thread_introns_into(&self, introns: &[Iv], out: &mut Vec<Vertex>) -> bool {
        out.clear();
        for iv in introns {
            if self.discarded.contains(iv) {
                return false;
            }
            let sub = self.graph.intron_correction_map.get(iv).copied().unwrap_or(*iv);
            out.push(Vertex::intron(sub));
        }
        true
    }

    /// `thread_ends` (`:95-136`).
    pub fn thread_ends(&self, intron: Vertex, end: i64, trusted: bool) -> Option<Vertex> {
        let s = self.right.get(&intron)?;
        let d = self.apa_delta;
        if trusted {
            if let Some(p) = Side::first_within(&s.tail, end, d) {
                return Some(Vertex(vc::POLYA, p));
            }
        }
        if let Some(b) = s.intron_bound {
            // rightmost_exon_end = max(outgoing intron starts) - 1
            if !trusted && end <= b - 1 + self.delta {
                return None;
            }
        }
        if let Some(p) = Side::first_within(&s.tss, end, d) {
            return Some(Vertex(vc::TSS_RIGHT, p));
        }
        let ext = s.extreme?;
        if trusted {
            (end >= ext.1 && ext.0 == vc::READ_END).then_some(ext)
        } else {
            (end <= ext.1 + d && s.next_pos.is_none_or(|p| end > p)).then_some(ext)
        }
    }

    /// `thread_starts` (`:138-163`).
    pub fn thread_starts(&self, intron: Vertex, start: i64, trusted: bool) -> Option<Vertex> {
        let s = self.left.get(&intron)?;
        let d = self.apa_delta;
        if trusted {
            if let Some(p) = Side::first_within(&s.tail, start, d) {
                return Some(Vertex(vc::POLYT, p));
            }
        }
        if let Some(b) = s.intron_bound {
            // leftmost_exon_start = min(incoming intron ends) + 1
            if !trusted && start >= b + 1 - self.delta {
                return None;
            }
        }
        if let Some(p) = Side::first_within(&s.tss, start, d) {
            return Some(Vertex(vc::TSS_LEFT, p));
        }
        let ext = s.extreme?;
        if trusted {
            (start <= ext.1 && ext.0 == vc::READ_START).then_some(ext)
        } else {
            (start >= ext.1 && s.next_pos.is_none_or(|p| start < p)).then_some(ext)
        }
    }
}

/// Per-slot accumulator of [`fill_paths`].
struct Slot {
    ipath: u32,
    start: Option<Vertex>,
    end: Option<Vertex>,
    fl: bool,
    reads: Vec<u32>,
}

/// `IntronPathStorage.fill` (`intron_path.py:26-55`): `paths`, `fl_paths`, `paths_to_reads`
/// (`known_isoforms_in_graph` is left empty; see [`known_spliced_isoforms`]).
pub fn fill_paths(store: &BlockStore, proc: &PathProcessor<'_>, params: &Params) -> PathStorage {
    const UNSEEN: u32 = u32::MAX;
    const DISCARDED: u32 = u32::MAX - 1;
    // chain id -> intron-path id (or UNSEEN / DISCARDED; empty paths are DISCARDED too)
    let mut chain_path: Vec<u32> = vec![UNSEEN; store.chains.len()];
    let mut ipaths: Vec<Box<[Vertex]>> = Vec::new();
    let mut ipath_ids: HashMap<Box<[Vertex]>, u32> = HashMap::new();
    let mut slots: Vec<Slot> = Vec::new();
    let mut slot_of: HashMap<(u32, Option<Vertex>, Option<Vertex>), u32> = HashMap::new();
    let mut buf: Vec<Vertex> = Vec::new();

    for r in store.iter() {
        if r.multimapper() {
            continue;
        }
        let c = r.chain_id();
        if c == ChainTable::NONE {
            continue;
        }
        let pid = match chain_path[c as usize] {
            UNSEEN => {
                let ok = match r.corrected_introns_slice() {
                    Some(sl) => proc.thread_introns_into(sl, &mut buf),
                    None => proc.thread_introns_into(&r.corrected_introns(), &mut buf),
                };
                let id = if !ok || buf.is_empty() {
                    DISCARDED
                } else if let Some(&id) = ipath_ids.get(buf.as_slice()) {
                    id
                } else {
                    let id = ipaths.len() as u32;
                    let p: Box<[Vertex]> = buf.as_slice().into();
                    ipaths.push(p.clone());
                    ipath_ids.insert(p, id);
                    id
                };
                chain_path[c as usize] = id;
                id
            }
            id => id,
        };
        if pid == DISCARDED {
            continue;
        }
        let ip = &ipaths[pid as usize];
        let pa = r.polya();
        let strand = r.strand();
        let end_trusted = strand == Strand::Plus && (pa.ext_a != -1 || pa.int_a != -1);
        let term = proc.thread_ends(ip[ip.len() - 1], r.end() as i64, end_trusted);
        let start_trusted = strand == Strand::Minus && (pa.ext_t != -1 || pa.int_t != -1);
        let first = proc.thread_starts(ip[0], r.start() as i64, start_trusted);

        let sid = *slot_of.entry((pid, first, term)).or_insert_with(|| {
            slots.push(Slot { ipath: pid, start: first, end: term, fl: false, reads: Vec::new() });
            (slots.len() - 1) as u32
        });
        let slot = &mut slots[sid as usize];
        if let (Some(t), Some(s)) = (term, first) {
            if !params.requires_polya_for_construction || t.0 == vc::POLYA || s.0 == vc::POLYT {
                slot.fl = true;
            }
        }
        slot.reads.push(r.idx());
    }

    let mut out = PathStorage::default();
    out.paths.reserve(slots.len());
    out.paths_to_reads.reserve(slots.len());
    for s in slots {
        let ip = &ipaths[s.ipath as usize];
        let mut p: Vec<Vertex> = Vec::with_capacity(ip.len() + 2);
        p.extend(s.start);
        p.extend_from_slice(ip);
        p.extend(s.end);
        let p: Path = p.into_boxed_slice();
        if s.fl {
            out.fl_paths.push(p.clone());
        }
        out.paths.insert(p.clone(), s.reads.len() as u32);
        out.paths_to_reads.insert(p, s.reads);
    }
    out.fl_paths.sort_unstable();
    out
}

/// `get_known_spliced_isoforms` (`context.py:202-210`): threaded intron path of every
/// reference isoform (in index order) → isoform. A later isoform with the same path
/// overwrites the value but keeps the first insertion position (Python dict semantics).
pub fn known_spliced_isoforms(index: &IsoformIndex, proc: &PathProcessor<'_>) -> IndexMap<Path, u32> {
    let mut out: IndexMap<Path, u32> = IndexMap::new();
    let mut buf = Vec::new();
    for (i, rec) in index.iso.iter().enumerate() {
        if !proc.thread_introns_into(index.iso_introns(i), &mut buf) || buf.is_empty() {
            continue;
        }
        let tx = match rec.key {
            IsoRef::Ref(t) => t,
            // model indices carry no reference transcript; 4.1 only calls this on GeneInfo
            IsoRef::Model(m) => m,
        };
        out.insert(buf.as_slice().into(), tx);
    }
    out
}

/// `fn thread_paths(&BlockStore, &IntronGraph, &IsoformIndex, &Params) -> PathStorage`
/// (`context.py:212-224` `build_graph` after the graph: fill + known isoforms).
pub fn thread_paths(store: &BlockStore, graph: &IntronGraph, index: &IsoformIndex, params: &Params) -> PathStorage {
    let proc = PathProcessor::new(graph, params);
    let mut ps = fill_paths(store, &proc, params);
    ps.known_isoforms_in_graph = known_spliced_isoforms(index, &proc);
    ps
}

pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U12Fill>();
    reg.add::<fixture::U12Known>();
}

pub mod fixture {
    use super::*;
    use iq_types::fixture::schema::{GraphFixture, U11GraphOut, U12FillIn, U12FillOut, U12KnownSplicedIn, U12KnownSplicedOut};
    use iq_types::fixture::{read_lines, unit_files, Case, FixError, FixResult, FixtureCtx, UnitFixture};
    use iq_types::serde_json;
    use std::path::PathBuf;
    use std::sync::{Arc, Mutex, OnceLock};

    type GraphCache = Mutex<HashMap<PathBuf, Arc<HashMap<String, Arc<IntronGraph>>>>>;

    /// Loads the U11 graph referenced by `key` from the capture dir's U11 files (all U11
    /// cases of a capture dir are loaded once and cached).
    pub fn graph(ctx: &FixtureCtx, key: &str) -> FixResult<Arc<IntronGraph>> {
        static CACHE: OnceLock<GraphCache> = OnceLock::new();
        let cache = CACHE.get_or_init(Default::default);
        let all = {
            let mut g = cache.lock().unwrap();
            match g.get(&ctx.root) {
                Some(m) => m.clone(),
                None => {
                    let mut m = HashMap::new();
                    for f in unit_files(&ctx.root, "U11")? {
                        for case in read_lines::<Case>(&f)? {
                            let case = case?;
                            // cases recording a Python exception carry no graph
                            if let Ok(out) = serde_json::from_value::<U11GraphOut>(case.out) {
                                m.insert(case.case, Arc::new(IntronGraph::from(&out.graph as &GraphFixture)));
                            }
                        }
                    }
                    let m = Arc::new(m);
                    g.insert(ctx.root.clone(), m.clone());
                    m
                }
            }
        };
        all.get(key).cloned().ok_or_else(|| FixError(format!("U11 case {key} not found")))
    }

    pub struct U12Fill;
    impl UnitFixture for U12Fill {
        const UNIT: &'static str = "U12";
        const FN: &'static str = "IntronPathStorage.fill";
        type In = U12FillIn;
        type Out = U12FillOut;
        fn run(input: &U12FillIn, ctx: &FixtureCtx) -> FixResult<U12FillOut> {
            let g = graph(ctx, &input.graph)?;
            let block = ctx.block(&input.block)?;
            let params = ctx.params()?;
            let proc = PathProcessor::new(&g, &params);
            Ok(U12FillOut::from_core(&fill_paths(&block.store, &proc, &params)))
        }
    }

    pub struct U12Known;
    impl UnitFixture for U12Known {
        const UNIT: &'static str = "U12";
        const FN: &'static str = "get_known_spliced_isoforms";
        type In = U12KnownSplicedIn;
        type Out = U12KnownSplicedOut;
        fn run(input: &U12KnownSplicedIn, ctx: &FixtureCtx) -> FixResult<U12KnownSplicedOut> {
            let g = graph(ctx, &input.graph)?;
            let index = ctx.index(&input.index)?;
            let params = ctx.params()?;
            let proc = PathProcessor::new(&g, &params);
            let known = known_spliced_isoforms(&index, &proc);
            let ids = ctx.ids();
            let space = index.iso.first().map(|r| matches!(r.key, IsoRef::Model(_))).unwrap_or(false);
            Ok(U12KnownSplicedOut {
                known: known
                    .iter()
                    .map(|(p, &t)| {
                        let r = if space { IsoRef::Model(t) } else { IsoRef::Ref(t) };
                        (p.to_vec(), ids.iso_name(r).to_owned())
                    })
                    .collect(),
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Literal port of 4.1 `get_outgoing/get_incoming` + `thread_ends/thread_starts`.
    #[allow(clippy::if_same_then_else)] // kept branch-for-branch with the Python
    mod naive {
        use super::*;
        fn get(m: &IndexMap<Vertex, Vec<Vertex>>, v: Vertex, t: Option<i64>) -> Vec<Vertex> {
            let mut r: Vec<Vertex> = match m.get(&v) {
                None => return vec![],
                Some(l) => l.iter().copied().filter(|x| match t { None => x.0 >= 0, Some(t) => x.0 == t }).collect(),
            };
            r.sort();
            r
        }
        pub fn ends(g: &IntronGraph, intron: Vertex, end: i64, trusted: bool, apa: i64, delta: i64) -> Option<Vertex> {
            let polyas = get(&g.outgoing, intron, Some(vc::POLYA));
            if trusted {
                if let Some(v) = polyas.iter().find(|v| (v.1 - end).abs() <= apa) {
                    return Some(*v);
                }
            }
            let oi = get(&g.outgoing, intron, None);
            if !oi.is_empty() {
                let r = oi.iter().map(|v| v.0).max().unwrap() - 1;
                if !trusted && end <= r + delta {
                    return None;
                }
            }
            let tss = get(&g.outgoing, intron, Some(vc::TSS_RIGHT));
            if let Some(v) = tss.iter().find(|v| (v.1 - end).abs() <= apa) {
                return Some(*v);
            }
            let mut all = get(&g.outgoing, intron, Some(vc::READ_END));
            all.extend(tss);
            all.extend(polyas);
            all.sort_by_key(|v| v.1);
            let r = *all.last()?;
            if trusted && end >= r.1 && r.0 == vc::READ_END {
                Some(r)
            } else if !trusted && end <= r.1 + apa && (all.len() <= 1 || end > all[all.len() - 2].1) {
                Some(r)
            } else {
                None
            }
        }
        pub fn starts(g: &IntronGraph, intron: Vertex, start: i64, trusted: bool, apa: i64, delta: i64) -> Option<Vertex> {
            let polyts = get(&g.incoming, intron, Some(vc::POLYT));
            if trusted {
                if let Some(v) = polyts.iter().find(|v| (v.1 - start).abs() <= apa) {
                    return Some(*v);
                }
            }
            let ii = get(&g.incoming, intron, None);
            if !ii.is_empty() {
                let l = ii.iter().map(|v| v.1).min().unwrap() + 1;
                if !trusted && start >= l - delta {
                    return None;
                }
            }
            let tss = get(&g.incoming, intron, Some(vc::TSS_LEFT));
            if let Some(v) = tss.iter().find(|v| (v.1 - start).abs() <= apa) {
                return Some(*v);
            }
            let mut all = get(&g.incoming, intron, Some(vc::READ_START));
            all.extend(tss);
            all.extend(polyts);
            all.sort_by_key(|v| v.1);
            let l = *all.first()?;
            if trusted && start <= l.1 && l.0 == vc::READ_START {
                Some(l)
            } else if !trusted && start >= l.1 && (all.len() <= 1 || start < all[1].1) {
                Some(l)
            } else {
                None
            }
        }
    }

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self, n: u64) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (self.0 >> 33) % n
        }
    }

    fn params() -> Params {
        Params::pacbio()
    }

    /// Random small graphs with dense, tied terminal positions: the precomputed summaries must
    /// reproduce the literal 4.1 functions on every (vertex, position, trusted) query.
    #[test]
    fn thread_terminals_match_literal_port() {
        let p = params();
        let (apa, delta) = (p.matching.apa_delta as i64, p.matching.delta as i64);
        let mut rng = Lcg(42);
        for _ in 0..2000 {
            let mut g = IntronGraph::default();
            let v = Vertex(1000, 1100);
            let mut out = Vec::new();
            let mut inc = Vec::new();
            for _ in 0..rng.next(3) {
                let s = 1100 + rng.next(200) as i64;
                out.push(Vertex(s, s + 50));
            }
            for _ in 0..rng.next(3) {
                let e = 900 - rng.next(200) as i64;
                inc.push(Vertex(e - 50, e));
            }
            for _ in 0..rng.next(6) {
                let c = [vc::POLYA, vc::READ_END, vc::TSS_RIGHT][rng.next(3) as usize];
                out.push(Vertex(c, 1100 + rng.next(150) as i64));
            }
            for _ in 0..rng.next(6) {
                let c = [vc::POLYT, vc::READ_START, vc::TSS_LEFT][rng.next(3) as usize];
                inc.push(Vertex(c, 1000 - rng.next(150) as i64));
            }
            out.sort();
            out.dedup();
            inc.sort();
            inc.dedup();
            g.outgoing.insert(v, out);
            g.incoming.insert(v, inc);
            let proc = PathProcessor::new(&g, &p);
            for x in 1050..1350 {
                for t in [false, true] {
                    assert_eq!(proc.thread_ends(v, x, t), naive::ends(&g, v, x, t, apa, delta), "end {x} {t} {g:?}");
                }
            }
            for x in 750..1050 {
                for t in [false, true] {
                    assert_eq!(proc.thread_starts(v, x, t), naive::starts(&g, v, x, t, apa, delta), "start {x} {t} {g:?}");
                }
            }
            // vertex missing from the graph
            assert_eq!(proc.thread_ends(Vertex(1, 2), 1200, true), None);
            assert_eq!(proc.thread_starts(Vertex(1, 2), 900, false), None);
        }
    }

    #[test]
    fn thread_introns_substitutes_and_discards() {
        let mut g = IntronGraph::default();
        g.intron_correction_map.insert(Iv::new(10, 20), Iv::new(11, 20));
        g.discarded_introns = vec![Iv::new(50, 60)];
        let p = params();
        let proc = PathProcessor::new(&g, &p);
        let mut buf = Vec::new();
        assert!(proc.thread_introns_into(&[Iv::new(10, 20), Iv::new(30, 40)], &mut buf));
        assert_eq!(buf, vec![Vertex(11, 20), Vertex(30, 40)]);
        assert!(!proc.thread_introns_into(&[Iv::new(10, 20), Iv::new(50, 60)], &mut buf));
    }
}
