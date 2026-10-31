//! `IntronGraph.__init__` (`intron_graph.py:148-809`): construction, simplification and
//! terminal-vertex attachment.
//!
//! The 4.1 graph is a pair of `defaultdict(set)`s plus the collector's dicts, and its
//! captured state depends on dict order, including the empty entries that `defaultdict`
//! *reads* insert (`is_isolated`, `signleton_dead_end`, `is_start_internal`,
//! `clustered_introns[v]` on collapsed vertices, ...). Every such read is mirrored with
//! [`PyDict::entry_default`] at the same point in the same order where it affects content
//! (e.g. a 0-count `clustered_introns` entry read through a stale set).
//!
//! Adjacency sets are iterated in sorted order, while 4.1 iterates in CPython hash order. That
//! only changes the dict *key order* of `outgoing_edges`/`incoming_edges` (empty entries
//! inserted by `signleton_dead_end` reads), never content: every removal phase decides
//! before it mutates, and later loops use `sorted(keys)` or order-free folds. The fixture
//! captures that key order, but no output depends on it: the only downstream readers of
//! the dicts in key order are `intron_path.py:64-67` (folds the values into a set) and
//! `assignment_based_constructor.py:108-135` (terminal-exon lists used only in an
//! existence check). The U11 runner therefore compares adjacency as unordered maps.

use crate::collector::{chain_introns, IntronCollector};
use crate::pyds::{find_closest, PyDict, VSet};
use iq_types::indexmap::IndexMap;
use iq_types::{
    vertex_code as vc, BlockStore, ChainTable, GeneKey, IntronGraph, IsoRef, IsoformIndex, Iv, Params, Pos,
    Predictions, Strand, TxIdx, Vertex,
};
use std::collections::{BTreeMap, BinaryHeap, HashMap, HashSet};

/// Parameters `IntronGraph` reads (`self.params.*`).
#[derive(Clone, Copy, Debug)]
pub struct GraphParams {
    pub delta: i64,
    pub apa_delta: i64,
    pub min_novel_intron_count: u32,
    pub graph_clustering_distance: i64,
    pub graph_clustering_ratio: f64,
    pub singleton_adjacent_cov: u32,
    pub min_novel_isolated_intron_abs: u32,
    pub terminal_position_abs: u32,
    pub terminal_position_rel: f64,
    pub terminal_internal_position_rel: f64,
}

impl GraphParams {
    pub fn from_params(p: &Params) -> Self {
        GraphParams {
            delta: p.matching.delta as i64,
            apa_delta: p.matching.apa_delta as i64,
            min_novel_intron_count: p.model.min_novel_intron_count,
            graph_clustering_distance: p.model.graph_clustering_distance as i64,
            graph_clustering_ratio: p.model.graph_clustering_ratio,
            singleton_adjacent_cov: p.model.singleton_adjacent_cov,
            min_novel_isolated_intron_abs: p.model.min_novel_isolated_intron_abs,
            terminal_position_abs: p.model.terminal_position_abs,
            terminal_position_rel: p.model.terminal_position_rel,
            terminal_internal_position_rel: p.model.terminal_internal_position_rel,
        }
    }
}

type Adj = PyDict<Vertex, VSet>;
/// Read terminal positions of one intron: position → count, insertion order.
type PosCounts = PyDict<i64, u32>;

/// Mutable graph state during construction (mirrors the 4.1 object).
struct Builder<'a> {
    p: GraphParams,
    index: &'a IsoformIndex,
    ic: IntronCollector,
    out: Adj,
    inc: Adj,
    starting_known: PyDict<Vertex, Vec<i64>>,
    terminal_known: PyDict<Vertex, Vec<i64>>,
    max_coverage: u32,
}

/// `fn build_graph(&BlockStore, &IsoformIndex, &Predictions, &Params) -> IntronGraph` (U11).
pub fn build_graph(store: &BlockStore, index: &IsoformIndex, preds: &Predictions, params: &Params) -> IntronGraph {
    build_graph_with(store, index, preds, &GraphParams::from_params(params))
}

/// [`build_graph`] with explicit graph parameters.
pub fn build_graph_with(store: &BlockStore, index: &IsoformIndex, preds: &Predictions, p: &GraphParams) -> IntronGraph {
    let mut b = Builder {
        p: *p,
        index,
        ic: IntronCollector::new(&index.introns, p.delta as i32),
        out: Adj::default(),
        inc: Adj::default(),
        starting_known: PyDict::default(),
        terminal_known: PyDict::default(),
        max_coverage: 0,
    };
    // __init__ :174-180
    for iso in &index.iso {
        if iso.n_intron == 0 {
            continue;
        }
        let introns = &index.iso_introns[iso.intron_off as usize..(iso.intron_off + iso.n_intron) as usize];
        let exons = &index.iso_exons[iso.exon_off as usize..(iso.exon_off + iso.n_exon) as usize];
        b.starting_known.entry_default(Vertex::intron(introns[0])).push(exons[0].s as i64);
        b.terminal_known.entry_default(Vertex::intron(introns[introns.len() - 1])).push(exons[exons.len() - 1].e as i64);
    }
    // IntronCollector.process
    let all = IntronCollector::collect_introns(store);
    b.ic.cluster_introns(&all, p.min_novel_intron_count);
    b.construct(store);
    b.simplify();
    b.attach_terminal_positions(store, preds);
    b.finish()
}

impl Builder<'_> {
    // ---------------------------------------------------------------------------- build

    /// `add_edge` + `construct` (`:192-199`, `:277-288`). Edges come from each distinct
    /// chain once, in order of first occurrence: repeated chains never create dict keys.
    fn construct(&mut self, store: &BlockStore) {
        let mut seen = vec![false; store.chains.len()];
        let mut chain: Vec<Vertex> = Vec::new();
        for r in store.iter() {
            let c = r.chain_id();
            if r.multimapper() || c == ChainTable::NONE || seen[c as usize] {
                continue;
            }
            seen[c as usize] = true;
            chain.clear();
            chain.extend(chain_introns(&store.chains, c));
            if chain.iter().any(|v| self.ic.discarded.contains(v)) {
                continue;
            }
            for w in chain.windows(2) {
                let v1 = self.ic.substitute(w[0]);
                let v2 = self.ic.substitute(w[1]);
                self.out.entry_default(v1).insert(v2);
                self.inc.entry_default(v2).insert(v1);
            }
        }
        self.max_coverage = self.ic.clustered.iter().map(|(_, &c)| c).max().unwrap_or(0);
    }

    // ------------------------------------------------------------------------- simplify

    fn simplify(&mut self) {
        self.clean_tips_and_bulges();
        self.remove_singleton_dead_ends();
        self.remove_isolates();
        self.ic.simplify_correction_map();
    }

    /// `collapse_vertex` (`:264-275`), with 4.1's evaluation order of the dict reads.
    fn collapse_vertex(&mut self, x: Vertex, s: Vertex) {
        self.out.entry_default(s);
        let out_x = self.out.entry_default(x).clone();
        self.out.get_mut(&s).expect("created").extend(&out_x);
        for i in out_x.iter() {
            let set = self.inc.entry_default(*i);
            set.remove(&x);
            set.insert(s);
        }
        self.inc.entry_default(s);
        let inc_x = self.inc.entry_default(x).clone();
        self.inc.get_mut(&s).expect("created").extend(&inc_x);
        for i in inc_x.iter() {
            let set = self.out.entry_default(*i);
            set.remove(&x);
            set.insert(s);
        }
        self.ic.add_substitute(x, s);
    }

    /// `collapse_vertex_set` (`:407-434`): returns `substitute_dict` sorted by key.
    ///
    /// Same result as 4.1 in O(n log n + n·w): the approved set is a sorted vector, and for
    /// each vertex only approved introns with `|Δstart| < distance` (a contiguous range) are
    /// scanned. The strict `<` thresholds are kept (design §7.2).
    fn collapse_vertex_set(&mut self, set: &[Vertex]) -> Vec<(Vertex, Vertex)> {
        if set.len() <= 1 {
            return Vec::new();
        }
        let d = self.p.graph_clustering_distance;
        let ratio = self.p.graph_clustering_ratio;
        let mut counts: Vec<(u32, Vertex)> = set.iter().map(|&v| (*self.ic.clustered.entry_default(v), v)).collect();
        counts.sort_unstable_by(|a, b| b.cmp(a));
        let mut approved: Vec<(Vertex, u32)> = Vec::new();
        let mut subs: Vec<(Vertex, Vertex)> = Vec::new();
        for &(count, v) in &counts {
            let lo = approved.partition_point(|a| a.0 .0 <= v.0 - d);
            let hi = approved.partition_point(|a| a.0 .0 < v.0 + d);
            let mut best: Option<(i64, Vertex)> = None;
            for &(a, ca) in &approved[lo..hi] {
                let sd = (a.0 - v.0).abs();
                let ed = (a.1 - v.1).abs();
                if sd < d && ed < d && (count as f64) < (ca as f64) * ratio {
                    let key = (sd + ed, a);
                    if best.is_none_or(|b| key < b) {
                        best = Some(key);
                    }
                }
            }
            match best {
                None => {
                    let at = approved.partition_point(|a| a.0 < v);
                    approved.insert(at, (v, count));
                }
                Some((_, a)) => subs.push((v, a)),
            }
        }
        subs.sort_unstable();
        subs
    }

    /// One direction of `clean_tips_and_bulges` (`:297-328`).
    fn clean_side(&mut self, outgoing: bool) {
        let keys = if outgoing { self.out.keys_sorted() } else { self.inc.keys_sorted() };
        let mut to_remove: HashSet<Vertex> = HashSet::new();
        for cur in keys {
            let adj = if outgoing { &mut self.out } else { &mut self.inc };
            let set: Vec<Vertex> = adj.entry_default(cur).to_vec();
            for (i, s) in self.collapse_vertex_set(&set) {
                if !to_remove.insert(i) {
                    continue;
                }
                self.collapse_vertex(i, s);
            }
        }
        for i in &to_remove {
            self.out.remove(i);
            self.inc.remove(i);
        }
    }

    fn clean_tips_and_bulges(&mut self) {
        self.clean_side(true);
        self.clean_side(false);
    }

    /// `signleton_dead_end` / `signleton_dead_start` (`:211-233`).
    fn singleton_dead_path(&mut self, mut v: Vertex, forward: bool) -> Vec<Vertex> {
        let mut path: Vec<Vertex> = Vec::new();
        while *self.ic.clustered.entry_default(v) == 1 {
            if path.contains(&v) {
                break; // a cycle: 4.1 would not terminate (needs a self loop, never seen)
            }
            path.push(v);
            let adj = if forward { &mut self.out } else { &mut self.inc };
            let next = adj.entry_default(v);
            if next.len() != 1 {
                break;
            }
            v = *next.iter().next().expect("len 1");
        }
        let adj = if forward { &mut self.out } else { &mut self.inc };
        if adj.entry_default(v).is_empty() {
            path
        } else {
            Vec::new()
        }
    }

    /// One direction of `remove_singleton_dead_ends` (`:330-370`).
    fn dead_end_side(&mut self, forward: bool) {
        let keys = if forward { self.out.keys_sorted() } else { self.inc.keys_sorted() };
        let mut to_clean: Vec<(Vertex, Vec<Vertex>)> = Vec::new();
        for cur in keys {
            if *self.ic.clustered.entry_default(cur) < self.p.singleton_adjacent_cov {
                continue;
            }
            let adj = if forward { &self.out } else { &self.inc };
            let nbrs: Vec<Vertex> = adj.get(&cur).map(VSet::to_vec).unwrap_or_default();
            let paths: Vec<Vec<Vertex>> = nbrs.iter().map(|&i| self.singleton_dead_path(i, forward)).collect();
            if paths.iter().any(Vec::is_empty) {
                continue;
            }
            to_clean.push((cur, paths.concat()));
        }
        for (intron, clean) in to_clean {
            let adj = if forward { &mut self.out } else { &mut self.inc };
            adj.insert(intron, VSet::default());
            for i in &clean {
                self.out.remove(i);
                self.inc.remove(i);
            }
        }
    }

    fn remove_singleton_dead_ends(&mut self) {
        self.dead_end_side(true);
        self.dead_end_side(false);
    }

    /// `remove_isolates` (`:372-405`).
    fn remove_isolates(&mut self) {
        let keys: Vec<Vertex> = self.ic.clustered.keys().collect();
        let mut isolated: Vec<Vertex> = Vec::new();
        for v in keys {
            // `not out[v] and not in[v]`: short-circuit, so in[v] is read only if out[v] is empty
            if self.out.entry_default(v).is_empty() && self.inc.entry_default(v).is_empty() {
                isolated.push(v);
            }
        }
        isolated.sort_unstable();
        let subs = self.collapse_vertex_set(&isolated);
        for &(i, s) in &subs {
            self.collapse_vertex(i, s);
        }
        for (i, _) in &subs {
            self.out.remove(i);
            self.inc.remove(i);
        }
        let cutoff = self.p.min_novel_isolated_intron_abs;
        let mut to_remove: Vec<Vertex> = Vec::new();
        for &v in &isolated {
            let Some(&c) = self.ic.clustered.get(&v) else { continue };
            if self.ic.known_introns.contains(&v) {
                continue;
            }
            if c < cutoff {
                self.ic.discard(v);
                to_remove.push(v);
            }
        }
        for i in &to_remove {
            self.out.remove(i);
            self.inc.remove(i);
        }
    }

    // ------------------------------------------------------------------------ terminals

    /// `attach_terminal_positions` (`:436-441`).
    fn attach_terminal_positions(&mut self, store: &BlockStore, preds: &Predictions) {
        let t = self.collect_terminal_positions(store);
        let introns = self.ic.clustered.keys_sorted();
        let gene_of = tx_genes(self.index);
        self.attach_side(&introns, &t.polya_ends, &t.read_ends, true, preds, &gene_of);
        self.attach_side(&introns, &t.polyt_starts, &t.read_starts, false, preds, &gene_of);
    }

    /// `collect_terminal_positions` (`:622-666`). Per-chain facts (discarded intron,
    /// substituted first/last intron, internal-check thresholds) are computed once per
    /// distinct chain; the per-read work is the position counting.
    fn collect_terminal_positions(&mut self, store: &BlockStore) -> TerminalPositions {
        #[derive(Clone, Copy)]
        enum ChainInfo {
            Unknown,
            Skip,
            Ok(ChainEnds),
        }
        #[derive(Clone, Copy)]
        struct ChainEnds {
            first: Vertex,
            last: Vertex,
            /// `min(inc[1] for inc in incoming[first])`, read lazily (creates the entry).
            min_inc_end: Option<Option<i64>>,
            /// `max(out[0] for out in outgoing[last])`, read lazily (creates the entry).
            max_out_start: Option<Option<i64>>,
        }
        let mut info = vec![ChainInfo::Unknown; store.chains.len()];
        let mut t = TerminalPositions::default();
        let delta = self.p.delta;
        for r in store.iter() {
            let c = r.chain_id();
            if r.multimapper() || c == ChainTable::NONE {
                continue;
            }
            if let ChainInfo::Unknown = info[c as usize] {
                let mut it = chain_introns(&store.chains, c);
                info[c as usize] = match it.next() {
                    None => ChainInfo::Skip,
                    Some(f) => {
                        let mut last = f;
                        let mut disc = self.ic.discarded.contains(&f);
                        for v in it {
                            disc |= self.ic.discarded.contains(&v);
                            last = v;
                        }
                        if disc {
                            ChainInfo::Skip
                        } else {
                            ChainInfo::Ok(ChainEnds {
                                first: self.ic.substitute(f),
                                last: self.ic.substitute(last),
                                min_inc_end: None,
                                max_out_start: None,
                            })
                        }
                    }
                };
            }
            let ChainInfo::Ok(ce) = &mut info[c as usize] else { continue };
            let read_start = r.start() as i64;
            if read_start >= ce.first.0 {
                continue; // substituted intron shifted left: skips the end side too
            }
            let pa = r.polya();
            let strand = r.strand();
            if strand == Strand::Minus && (pa.ext_t != -1 || pa.int_t != -1) {
                *t.polyt_starts.entry(ce.first).or_default().entry_default(read_start) += 1;
            } else {
                // is_start_internal: any(inc[1] - delta <= read_start)
                let first = ce.first;
                let m = *ce.min_inc_end.get_or_insert_with(|| self.inc.entry_default(first).iter().map(|v| v.1).min());
                if m.is_none_or(|m| m - delta > read_start) {
                    *t.read_starts.entry(ce.first).or_default().entry_default(read_start) += 1;
                }
            }
            let read_end = r.end() as i64;
            if read_end <= ce.last.1 {
                continue;
            }
            if strand == Strand::Plus && (pa.ext_a != -1 || pa.int_a != -1) {
                *t.polya_ends.entry(ce.last).or_default().entry_default(read_end) += 1;
            } else {
                // is_end_internal: any(out[0] + delta >= read_end)
                let last = ce.last;
                let m = *ce.max_out_start.get_or_insert_with(|| self.out.entry_default(last).iter().map(|v| v.0).max());
                if m.is_none_or(|m| m + delta < read_end) {
                    *t.read_ends.entry(ce.last).or_default().entry_default(read_end) += 1;
                }
            }
        }
        t
    }

    /// `_gene_terminal_introns` (`:443-456`).
    fn gene_terminal_introns(&self, read_end: bool) -> HashMap<GeneKey, HashSet<Vertex>> {
        let mut g: HashMap<GeneKey, HashSet<Vertex>> = HashMap::new();
        for iso in &self.index.iso {
            if iso.n_intron == 0 {
                continue;
            }
            let k = if read_end { iso.intron_off + iso.n_intron - 1 } else { iso.intron_off };
            let ti = Vertex::intron(self.index.iso_introns[k as usize]);
            g.entry(iso.gene).or_default().insert(self.ic.substitute(ti));
        }
        g
    }

    /// `_admissible_predictions` (`:458-503`): intron → sorted admissible positions.
    fn admissible(
        &self,
        introns: &HashSet<Vertex>,
        preds: &[(TxIdx, Pos)],
        read_end: bool,
        gene_of: &HashMap<TxIdx, GeneKey>,
    ) -> HashMap<Vertex, Vec<i64>> {
        let mut adm: HashMap<Vertex, Vec<i64>> = HashMap::new();
        if preds.is_empty() {
            return adm;
        }
        let terms = self.gene_terminal_introns(read_end);
        for &(tx, pos) in preds {
            let Some(g) = gene_of.get(&tx) else { continue };
            let Some(gi) = terms.get(g) else { continue };
            let pos = pos as i64;
            for v in gi {
                if !introns.contains(v) {
                    continue;
                }
                if read_end && pos <= v.1 || !read_end && pos >= v.0 {
                    continue;
                }
                adm.entry(*v).or_default().push(pos);
            }
        }
        for l in adm.values_mut() {
            l.sort_unstable();
            l.dedup();
        }
        adm
    }

    /// `_refine_positions` (`:608-620`).
    fn refine(&self, clustered: PosCounts, predicted: Option<&Vec<i64>>) -> PosCounts {
        let Some(pred) = predicted.filter(|p| !p.is_empty()) else { return clustered };
        if clustered.is_empty() {
            return clustered;
        }
        let mut refined = PosCounts::default();
        for (&pos, &count) in clustered.iter() {
            let new_pos = match find_closest(pos, pred) {
                Some((n, d)) if d <= self.p.apa_delta => n,
                _ => pos,
            };
            *refined.entry_default(new_pos) += count;
        }
        refined
    }

    /// `cluster_polya_positions` (`:668-700`): greedy max-count clustering with ties to the
    /// first-inserted position, as a lazy-deletion heap + ordered position map
    /// (O(P log P) instead of O(P²)).
    fn cluster_polya_positions(&mut self, positions: Option<&PosCounts>, intron: Vertex, read_end: bool) -> PosCounts {
        let mut clustered = PosCounts::default();
        let Some(positions) = positions.filter(|p| !p.is_empty()) else { return clustered };
        let known: Vec<i64> = if read_end {
            self.terminal_known.entry_default(intron).clone()
        } else {
            self.starting_known.entry_default(intron).clone()
        };
        let entries: Vec<(i64, u32)> = positions.iter().map(|(&p, &c)| (p, c)).collect();
        let mut alive = vec![true; entries.len()];
        let mut by_pos: BTreeMap<i64, usize> = entries.iter().enumerate().map(|(i, &(p, _))| (p, i)).collect();
        let mut heap: BinaryHeap<(u32, std::cmp::Reverse<usize>)> =
            entries.iter().enumerate().map(|(i, &(_, c))| (c, std::cmp::Reverse(i))).collect();
        let apa = self.p.apa_delta;
        while let Some((_, std::cmp::Reverse(i))) = heap.pop() {
            if !alive[i] {
                continue;
            }
            let mut top = entries[i].0;
            if let Some((n, d)) = find_closest(top, &known) {
                // `if nearest_position and ...`: position 0 is falsy in 4.1
                if n != 0 && d <= apa {
                    top = n;
                }
            }
            // 4.1 asserts top > intron[1] (read end) / top < intron[0]; it holds because
            // known positions lie in the terminal exon and read positions were filtered.
            let in_range: Vec<i64> = by_pos.range(top - apa..=top + apa).map(|(&p, _)| p).collect();
            let mut total = 0u32;
            for p in in_range {
                let k = by_pos.remove(&p).expect("present");
                alive[k] = false;
                total += entries[k].1;
            }
            clustered.insert(top, total);
        }
        let max_count = clustered.iter().map(|(_, &c)| c).max().unwrap_or(0);
        if max_count == self.p.terminal_position_abs {
            return clustered;
        }
        let cutoff = (max_count as f64 * self.p.terminal_position_rel).max(self.p.terminal_position_abs as f64);
        let mut result = PosCounts::default();
        for (&k, &v) in clustered.iter() {
            if v as f64 >= cutoff {
                result.insert(k, v);
            }
        }
        result
    }

    /// `_attach_side` (`:505-606`).
    fn attach_side(
        &mut self,
        introns: &[Vertex],
        polya_confirmed: &HashMap<Vertex, PosCounts>,
        read_terminal: &HashMap<Vertex, PosCounts>,
        read_end: bool,
        preds: &Predictions,
        gene_of: &HashMap<TxIdx, GeneKey>,
    ) {
        let intron_set: HashSet<Vertex> = introns.iter().copied().collect();
        let polya_adm = self.admissible(&intron_set, &preds.polya, read_end, gene_of);
        let tss_adm = self.admissible(&intron_set, &preds.tss, read_end, gene_of);
        let apa = self.p.apa_delta;

        // Step 1: polyA / polyT confirmed positions.
        let mut polya_pos: Vec<PosCounts> = Vec::with_capacity(introns.len());
        for &v in introns {
            let cl = self.cluster_polya_positions(polya_confirmed.get(&v), v, read_end);
            polya_pos.push(self.refine(cl, polya_adm.get(&v)));
        }

        // Step 1.5: TSS-confirmed positions; the rest stays for the cutoff path.
        let mut tss_pos: Vec<PosCounts> = Vec::with_capacity(introns.len());
        let mut remaining: Vec<Vec<(i64, u32)>> = Vec::with_capacity(introns.len());
        for &v in introns {
            let mut confirmed = PosCounts::default();
            let mut rem: Vec<(i64, u32)> = Vec::new();
            let tss = tss_adm.get(&v).map(Vec::as_slice).unwrap_or(&[]);
            if let Some(pc) = read_terminal.get(&v) {
                for (&pos, &count) in pc.iter() {
                    match find_closest(pos, tss) {
                        Some((n, d)) if d <= apa => *confirmed.entry_default(n) += count,
                        _ => rem.push((pos, count)),
                    }
                }
            }
            tss_pos.push(confirmed);
            remaining.push(rem);
        }

        // Step 2: cutoffs and extra (non-confirmed) positions. Step 3: cluster + refine.
        let term_adm = if read_end { &polya_adm } else { &tss_adm };
        let mut term_pos: Vec<PosCounts> = Vec::with_capacity(introns.len());
        for (k, &v) in introns.iter().enumerate() {
            let mut cutoff = self.p.terminal_position_abs as f64;
            let mut confirmed: HashMap<i64, u32> = polya_pos[k].iter().map(|(&p, &c)| (p, c)).collect();
            for (&p, &c) in tss_pos[k].iter() {
                *confirmed.entry(p).or_insert(0) += c;
            }
            let extra: Vec<(i64, u32)> = if !confirmed.is_empty() {
                let maxc = *confirmed.values().max().expect("non-empty");
                cutoff = cutoff.max(maxc as f64 * self.p.terminal_position_rel);
                let furthest =
                    if read_end { *confirmed.keys().max().expect("ne") } else { *confirmed.keys().min().expect("ne") };
                remaining[k]
                    .iter()
                    .copied()
                    .filter(|&(p, _)| if read_end { p >= furthest + apa } else { p <= furthest - apa })
                    .collect()
            } else {
                std::mem::take(&mut remaining[k])
            };
            let nbrs: Vec<Vertex> = {
                let adj = if read_end { &self.out } else { &self.inc };
                adj.get(&v).map(VSet::to_vec).unwrap_or_default()
            };
            if !nbrs.is_empty() {
                let cov = nbrs.iter().map(|&n| *self.ic.clustered.entry_default(n)).max().expect("ne");
                cutoff = cutoff.max(cov as f64 * self.p.terminal_internal_position_rel);
            }
            // cluster_terminal_positions (`:702-711`)
            let mut cl = PosCounts::default();
            if !extra.is_empty() {
                let total: u32 = extra.iter().map(|x| x.1).sum();
                let pos = if read_end {
                    extra.iter().map(|x| x.0).max().expect("ne")
                } else {
                    extra.iter().map(|x| x.0).min().expect("ne")
                };
                if (total as f64) >= cutoff {
                    cl.insert(pos, total);
                }
            }
            term_pos.push(self.refine(cl, term_adm.get(&v)));
        }

        // Step 4: attach terminal vertices.
        let (pv, rv, tv) = if read_end {
            (vc::POLYA, vc::READ_END, vc::TSS_RIGHT)
        } else {
            (vc::POLYT, vc::READ_START, vc::TSS_LEFT)
        };
        let edges = if read_end { &mut self.out } else { &mut self.inc };
        for (k, &v) in introns.iter().enumerate() {
            for (code, pc) in [(pv, &polya_pos[k]), (tv, &tss_pos[k]), (rv, &term_pos[k])] {
                for (&pos, _) in pc.iter() {
                    edges.entry_default(v).insert(Vertex(code, pos));
                }
            }
        }
    }

    // --------------------------------------------------------------------------- output

    fn finish(self) -> IntronGraph {
        let iv = |v: &Vertex| Iv::new(v.0 as Pos, v.1 as Pos);
        let mut discarded: Vec<Iv> = self.ic.discarded.iter().map(iv).collect();
        discarded.sort_unstable();
        let adj = |a: &Adj| -> IndexMap<Vertex, Vec<Vertex>> {
            a.iter().map(|(k, s)| (*k, s.to_vec())).collect()
        };
        let known = |d: &PyDict<Vertex, Vec<i64>>| -> IndexMap<Iv, Vec<Pos>> {
            d.iter().map(|(k, l)| (iv(k), l.iter().map(|&x| x as Pos).collect())).collect()
        };
        IntronGraph {
            clustered_introns: self.ic.clustered.iter().map(|(k, &c)| (iv(k), c)).collect(),
            intron_correction_map: self.ic.correction.iter().map(|(k, s)| (iv(k), iv(s))).collect(),
            discarded_introns: discarded,
            outgoing: adj(&self.out),
            incoming: adj(&self.inc),
            starting_known_positions: known(&self.starting_known),
            terminal_known_positions: known(&self.terminal_known),
            max_coverage: self.max_coverage,
        }
    }
}

/// Read terminus positions per (substituted) terminal intron.
#[derive(Default)]
struct TerminalPositions {
    polya_ends: HashMap<Vertex, PosCounts>,
    read_ends: HashMap<Vertex, PosCounts>,
    polyt_starts: HashMap<Vertex, PosCounts>,
    read_starts: HashMap<Vertex, PosCounts>,
}

/// `gene_info.gene_id_map` restricted to the index's reference transcripts.
fn tx_genes(index: &IsoformIndex) -> HashMap<TxIdx, GeneKey> {
    index
        .iso
        .iter()
        .filter_map(|i| match i.key {
            IsoRef::Ref(t) => Some((t, i.gene)),
            IsoRef::Model(_) => None,
        })
        .collect()
}

pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U11Graph>();
}

pub(crate) mod fixture {
    use super::*;
    use iq_types::fixture::ids::TxSpace;
    use iq_types::fixture::schema::{GraphFixture, U11GraphIn, U11GraphOut};
    use iq_types::fixture::{FixResult, FixtureCtx, UnitFixture, Verdict};
    use iq_types::fixture::diff::{diff_values, Diff};
    use iq_types::DevId;

    pub struct U11Graph;

    /// Converts the contract graph to the fixture DTO.
    pub fn to_fixture(g: &IntronGraph) -> GraphFixture {
        GraphFixture {
            clustered_introns: g.clustered_introns.iter().map(|(k, v)| (*k, *v)).collect(),
            intron_correction_map: g.intron_correction_map.iter().map(|(k, v)| (*k, *v)).collect(),
            discarded_introns: g.discarded_introns.clone(),
            outgoing: g.outgoing.iter().map(|(k, v)| (*k, v.clone())).collect(),
            incoming: g.incoming.iter().map(|(k, v)| (*k, v.clone())).collect(),
            starting_known_positions: g.starting_known_positions.iter().map(|(k, v)| (*k, v.clone())).collect(),
            terminal_known_positions: g.terminal_known_positions.iter().map(|(k, v)| (*k, v.clone())).collect(),
            max_coverage: g.max_coverage,
        }
    }

    impl UnitFixture for U11Graph {
        const UNIT: &'static str = "U11";
        const FN: &'static str = "IntronGraph.__init__";
        const DEVIATIONS: &'static [DevId] = &[DevId::G1, DevId::G2];
        type In = U11GraphIn;
        type Out = U11GraphOut;

        /// Adjacency (`outgoing`/`incoming`) compares as unordered maps key → sorted
        /// targets: 4.1's key order is a CPython set-iteration artifact that no output
        /// depends on (see module docs). Everything else is exact.
        fn compare(got: &Self::Out, want: &Self::Out, ctx: &FixtureCtx) -> Verdict {
            let norm = |o: &U11GraphOut| {
                let mut g = o.graph.clone();
                g.outgoing.sort_unstable();
                g.incoming.sort_unstable();
                U11GraphOut { graph: g }
            };
            let (g, w) = (norm(got), norm(want));
            if g == w {
                return Verdict::Exact;
            }
            match diff_values(&g, &w, &ctx.diff_opts) {
                Some(d) => Verdict::Mismatch(d),
                None => Verdict::Mismatch(Diff::error("PartialEq differs".to_owned())),
            }
        }

        fn run(input: &Self::In, ctx: &FixtureCtx) -> FixResult<Self::Out> {
            let block = ctx.block(&input.block)?;
            let index = ctx.index(&input.index)?;
            let params = ctx.params()?;
            let conv = |l: &Option<Vec<(String, Pos)>>| -> Vec<(TxIdx, Pos)> {
                let mut ids = ctx.ids();
                l.iter()
                    .flatten()
                    .filter_map(|(t, p)| match ids.iso_ref(t, TxSpace::Ref) {
                        IsoRef::Ref(i) => Some((i, *p)),
                        IsoRef::Model(_) => None,
                    })
                    .collect()
            };
            let preds = Predictions { polya: conv(&input.polya_predictions), tss: conv(&input.tss_predictions) };
            let g = build_graph(&block.store, &index, &preds, &params);
            Ok(U11GraphOut { graph: to_fixture(&g) })
        }
    }
}

#[cfg(test)]
mod fixture_tests {
    use super::fixture::U11Graph;
    use iq_types::fixture::schema::{U11GraphIn, U11GraphOut};
    use iq_types::fixture::{read_cases, unit_files, FixtureCtx, UnitFixture};
    use std::path::PathBuf;

    /// `IQRS_FIXTURES=<root> [IQRS_U11_DUMP=<dir>] cargo test -p iq_graph fixtures_u11`: runs every
    /// U11 case; with DUMP, writes `got`/`want` JSON of mismatching cases for inspection.
    #[test]
    fn fixtures_u11() {
        // parity check: needs the quirks41 build (the default build differs by D-G1)
        let Ok(root) = std::env::var("IQRS_FIXTURES") else { return };
        if !iq_types::QUIRKS41 {
            return;
        }
        let dump = std::env::var("IQRS_U11_DUMP").ok().map(PathBuf::from);
        let mut bad = 0;
        for sample in std::fs::read_dir(&root).unwrap().flatten() {
            for slice in std::fs::read_dir(sample.path()).into_iter().flatten().flatten() {
                let cap = slice.path();
                let Ok(files) = unit_files(&cap, "U11") else { continue };
                let ctx = FixtureCtx::new(&cap);
                for f in files {
                    for case in read_cases::<U11GraphIn, U11GraphOut>(&f).unwrap() {
                        let case = case.unwrap();
                        let got = U11Graph::run(&case.input, &ctx).unwrap();
                        if !U11Graph::compare(&got, &case.out, &ctx).is_ok() {
                            bad += 1;
                            if let Some(d) = &dump {
                                let n = case.case.replace(':', "_");
                                std::fs::write(d.join(format!("{n}.got.json")), iq_types::serde_json::to_string(&got).unwrap()).unwrap();
                                std::fs::write(d.join(format!("{n}.want.json")), iq_types::serde_json::to_string(&case.out).unwrap()).unwrap();
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(bad, 0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> GraphParams {
        GraphParams {
            delta: 6,
            apa_delta: 50,
            min_novel_intron_count: 1,
            graph_clustering_distance: 20,
            graph_clustering_ratio: 0.5,
            singleton_adjacent_cov: 10,
            min_novel_isolated_intron_abs: 3,
            terminal_position_abs: 1,
            terminal_position_rel: 0.01,
            terminal_internal_position_rel: 0.02,
        }
    }

    /// Literal port of `collapse_vertex_set` (O(n·a log a)).
    fn naive(counts: &HashMap<Vertex, u32>, set: &[Vertex], d: i64, ratio: f64) -> Vec<(Vertex, Vertex)> {
        if set.len() <= 1 {
            return vec![];
        }
        let mut vc: Vec<(u32, Vertex)> = set.iter().map(|v| (counts[v], *v)).collect();
        vc.sort_by(|a, b| b.cmp(a));
        let mut approved: Vec<Vertex> = vec![];
        let mut subs = vec![];
        for (c, v) in vc {
            let mut ap = approved.clone();
            ap.sort();
            let mut sim: Vec<(i64, Vertex)> = vec![];
            for i in ap {
                let (sd, ed) = ((i.0 - v.0).abs(), (i.1 - v.1).abs());
                if sd < d && ed < d && (c as f64) < counts[&i] as f64 * ratio {
                    sim.push((sd + ed, i));
                }
            }
            if sim.is_empty() {
                approved.push(v);
            } else {
                sim.sort();
                subs.push((v, sim[0].1));
            }
        }
        subs.sort();
        subs
    }

    #[test]
    fn collapse_vertex_set_matches_naive() {
        let index = IsoformIndex::default();
        let mut seed = 12345u64;
        let mut rnd = |n: u64| {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 33) % n
        };
        for _ in 0..300 {
            let mut b = Builder {
                p: params(),
                index: &index,
                ic: IntronCollector::new(&[], 6),
                out: Adj::default(),
                inc: Adj::default(),
                starting_known: PyDict::default(),
                terminal_known: PyDict::default(),
                max_coverage: 0,
            };
            let n = 1 + rnd(40) as usize;
            let mut counts = HashMap::new();
            let mut set = vec![];
            for _ in 0..n {
                let s = 1000 + rnd(80) as i64;
                let v = Vertex(s, s + 100 + rnd(60) as i64);
                if counts.contains_key(&v) {
                    continue;
                }
                let c = 1 + rnd(30) as u32;
                counts.insert(v, c);
                b.ic.clustered.insert(v, c);
                set.push(v);
            }
            assert_eq!(b.collapse_vertex_set(&set), naive(&counts, &set, 20, 0.5));
        }
    }
}
