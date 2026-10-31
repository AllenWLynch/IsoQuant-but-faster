//! `IntronCollector` (`intron_graph.py:41-146`): intron counting, similarity clustering
//! (4.1 tie and min-count behaviour), the correction map and discarded introns.

use crate::pyds::PyDict;
use iq_types::{BlockStore, ChainTable, Iv, Vertex};
use std::collections::{HashMap, HashSet};

/// Collector state. Keys are intron vertices (`Vertex(s, e)`, `s >= 0`).
#[derive(Clone, Debug, Default)]
pub struct IntronCollector {
    pub(crate) known_introns: HashSet<Vertex>,
    pub(crate) delta: i64,
    /// `clustered_introns` (a `defaultdict(int)`; see [`PyDict::entry_default`]).
    pub(crate) clustered: PyDict<Vertex, u32>,
    pub(crate) correction: PyDict<Vertex, Vertex>,
    pub(crate) discarded: HashSet<Vertex>,
}

/// The real introns (gaps with `s <= e`) of chain `c` (`junctions_from_blocks`).
#[inline]
pub(crate) fn chain_introns(chains: &ChainTable, c: u32) -> impl Iterator<Item = Vertex> + '_ {
    let gaps: &[Iv] = if c == ChainTable::NONE { &[] } else { chains.get(c) };
    gaps.iter().filter(|g| g.s <= g.e).map(|g| Vertex::intron(*g))
}

impl IntronCollector {
    pub fn new(known_introns: &[Iv], delta: i32) -> Self {
        IntronCollector {
            known_introns: known_introns.iter().map(|&iv| Vertex::intron(iv)).collect(),
            delta: delta as i64,
            ..Default::default()
        }
    }

    /// `collect_introns` (`:53-60`): counts per intron over non-multimapped reads with
    /// introns. Counted once per distinct chain × multiplicity, no per-read allocation.
    pub fn collect_introns(store: &BlockStore) -> HashMap<Vertex, u32> {
        let mut mult = vec![0u32; store.chains.len()];
        for r in store.iter() {
            let c = r.chain_id();
            if c != ChainTable::NONE && !r.multimapper() {
                mult[c as usize] += 1;
            }
        }
        let mut all: HashMap<Vertex, u32> = HashMap::new();
        for (c, &m) in mult.iter().enumerate() {
            if m == 0 {
                continue;
            }
            for v in chain_introns(&store.chains, c as u32) {
                *all.entry(v).or_insert(0) += m;
            }
        }
        all
    }

    /// `construct_similar_intron_map` (`:62-73`): for each intron, whether it has any
    /// similar intron (`|Δs| <= delta && |Δe| <= delta`), and the similar ones. `ordered`
    /// must be sorted. Returns adjacency lists indexed like `ordered`.
    fn similar_map(&self, ordered: &[Vertex]) -> Vec<Vec<u32>> {
        let mut sim: Vec<Vec<u32>> = vec![Vec::new(); ordered.len()];
        for i in 0..ordered.len() {
            let a = ordered[i];
            let mut j = i + 1;
            while j < ordered.len() && (ordered[j].0 - a.0).abs() <= self.delta {
                if (ordered[j].1 - a.1).abs() <= self.delta {
                    sim[i].push(j as u32);
                    sim[j].push(i as u32);
                }
                j += 1;
            }
        }
        sim
    }

    /// `cluster_introns` (`:75-101`).
    pub fn cluster_introns(&mut self, all: &HashMap<Vertex, u32>, min_count: u32) {
        let mut ordered: Vec<Vertex> = all.keys().copied().collect();
        ordered.sort_unstable();
        let sim = self.similar_map(&ordered);
        // sorted([(count, intron)], reverse=True)
        let mut by_count: Vec<(u32, u32)> = ordered.iter().enumerate().map(|(i, v)| (all[v], i as u32)).collect();
        by_count.sort_unstable_by(|a, b| b.cmp(a)); // index order == intron order
        for &(count, i) in &by_count {
            let intron = ordered[i as usize];
            if self.known_introns.contains(&intron) {
                self.clustered.insert(intron, count);
            } else if !sim[i as usize].is_empty() {
                // the already-clustered similar introns
                let mut best: Option<(u32, Vertex)> = None;
                for &j in &sim[i as usize] {
                    let s = ordered[j as usize];
                    if self.clustered.get(&s).is_none() {
                        continue;
                    }
                    // Matches IsoQuant 4.1: it sorts (count_of_current, similar), so the largest
                    // intron tuple wins regardless of the similar intron's own count.
                    let key = (0, s);
                    if best.is_none_or(|b| key > b) {
                        best = Some(key);
                    }
                }
                match best {
                    Some((_, sub)) => {
                        *self.clustered.entry_default(sub) += count;
                        self.correction.insert(intron, sub);
                    }
                    None => {
                        // Matches IsoQuant 4.1: the min-count test is skipped here.
                        self.clustered.insert(intron, count);
                    }
                }
            } else if count < min_count {
                self.discarded.insert(intron);
            } else {
                self.clustered.insert(intron, count);
            }
        }
    }

    /// `add_substitute` (`:111-114`). `clustered[sub] += clustered[orig]` evaluates
    /// `clustered[sub]` first (both are defaultdict reads).
    pub(crate) fn add_substitute(&mut self, orig: Vertex, sub: Vertex) {
        self.clustered.entry_default(sub);
        let c = *self.clustered.entry_default(orig);
        *self.clustered.entry_default(sub) += c;
        self.clustered.remove(&orig);
        self.correction.insert(orig, sub);
    }

    /// `discard` (`:116-119`).
    pub(crate) fn discard(&mut self, intron: Vertex) {
        self.discarded.insert(intron);
        self.clustered.remove(&intron);
    }

    /// `substitute` (`:121-124`): one level.
    #[inline]
    pub(crate) fn substitute(&self, v: Vertex) -> Vertex {
        self.correction.get(&v).copied().unwrap_or(v)
    }

    /// `simplify_correction_map` (`:126-145`).
    pub(crate) fn simplify_correction_map(&mut self) {
        let keys = self.correction.keys_sorted();
        let mut to_remove: Vec<Vertex> = Vec::new();
        let limit = self.correction.len() + 1;
        for intron in keys {
            let mut subs = self.correction.get(&intron).copied().expect("key");
            if self.discarded.contains(&subs) {
                to_remove.push(intron);
                continue;
            }
            if !self.correction.contains(&subs) {
                continue;
            }
            // 4.1 loops forever on a cycle; a cycle cannot arise (substitutes are always
            // clustered at the time they are chosen), the bound only guards against it.
            let mut steps = 0;
            while let Some(&n) = self.correction.get(&subs) {
                subs = n;
                steps += 1;
                if steps > limit {
                    break;
                }
            }
            if self.discarded.contains(&subs) {
                to_remove.push(intron);
                continue;
            }
            self.correction.insert(intron, subs);
        }
        for intron in to_remove {
            self.discard(intron);
            self.correction.remove(&intron);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(s: i64, e: i64) -> Vertex {
        Vertex(s, e)
    }

    #[test]
    fn cluster_basic() {
        let mut c = IntronCollector::new(&[], 4);
        let all: HashMap<Vertex, u32> = [(v(100, 200), 10), (v(102, 201), 2), (v(500, 600), 1)].into_iter().collect();
        c.cluster_introns(&all, 1);
        assert_eq!(
            c.clustered.iter().map(|(k, v)| (*k, *v)).collect::<Vec<_>>(),
            vec![(v(100, 200), 12), (v(500, 600), 1)]
        );
        assert_eq!(c.correction.get(&v(102, 201)), Some(&v(100, 200)));
    }

    /// Two clustered candidates: 4.1 takes the largest tuple, not the largest count.
    #[test]
    fn cluster_tie_g1() {
        // A=(100,200) count 10, B=(108,200) count 9 (not similar to A: Δs=8 > 4),
        // C=(104,200) count 1, similar to both.
        let all: HashMap<Vertex, u32> = [(v(100, 200), 10), (v(108, 200), 9), (v(104, 200), 1)].into_iter().collect();
        let mut c = IntronCollector::new(&[], 4);
        c.cluster_introns(&all, 1);
        assert_eq!(c.correction.get(&v(104, 200)).copied(), Some(v(108, 200)));
    }

    /// Similar neighbours that are all unclustered: 4.1 keeps the intron regardless of count.
    #[test]
    fn min_count_g2() {
        // A known? no. A=(100,200) count 1, B=(102,200) count 1, min_count 2. A is processed
        // first (larger tuple is B: (1,B) > (1,A)), so B first: similar A not clustered.
        let all: HashMap<Vertex, u32> = [(v(100, 200), 1), (v(102, 200), 1)].into_iter().collect();
        let mut c = IntronCollector::new(&[], 4);
        c.cluster_introns(&all, 2);
        // B kept (no min-count test), A then substitutes into B.
        assert_eq!((c.clustered.len(), c.discarded.len()), (1, 0));
    }
}
