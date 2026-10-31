//! Graph queries on the final [`IntronGraph`] used by later units (U15a).
//!
//! Note: in 4.1 these run on the live object, whose `defaultdict`s may have gained empty
//! entries from reads after `__init__` (e.g. path threading). Such an entry can add an
//! edgeless vertex to `get_overlapping_component_max_coverage`'s scan; callers that need
//! that must reproduce it. These functions use the graph as captured after `__init__`.

use iq_types::{Iv, IntronGraph, Vertex};
use std::collections::{HashSet, VecDeque};

fn count(g: &IntronGraph, v: Vertex) -> u32 {
    v.as_intron().and_then(|iv| g.clustered_introns.get(&iv).copied()).unwrap_or(0)
}

/// `get_connected_component` (`intron_graph.py:731-750`): intron vertices reachable from
/// the path's introns, edges taken in both directions.
pub fn connected_component(g: &IntronGraph, path: &[Vertex]) -> HashSet<Vertex> {
    let mut seen: HashSet<Vertex> = HashSet::new();
    let mut q: VecDeque<Vertex> = VecDeque::new();
    for &v in path {
        if v.0 < 0 {
            continue;
        }
        if seen.insert(v) {
            q.push_back(v);
        }
    }
    while let Some(v) = q.pop_front() {
        for adj in [&g.outgoing, &g.incoming] {
            if let Some(ns) = adj.get(&v) {
                for &n in ns {
                    if n.0 >= 0 && seen.insert(n) {
                        q.push_back(n);
                    }
                }
            }
        }
    }
    seen
}

/// `IntronGraph.get_max_component_coverage` (`intron_graph.py:752-756`).
pub fn max_component_coverage(g: &IntronGraph, path: &[Vertex]) -> u32 {
    connected_component(g, path).into_iter().map(|v| count(g, v)).max().unwrap_or(0)
}

/// `IntronGraph.get_overlapping_component_max_coverage` (`intron_graph.py:758-781`).
pub fn overlapping_component_max_coverage(g: &IntronGraph, iv: Iv) -> u32 {
    let mut processed: HashSet<Vertex> = HashSet::new();
    let mut all: Vec<Vertex> = g.outgoing.keys().chain(g.incoming.keys()).copied().collect();
    all.sort_unstable();
    all.dedup();
    let right_codes = [iq_types::vertex_code::POLYA, iq_types::vertex_code::READ_END, iq_types::vertex_code::TSS_RIGHT];
    for v in all {
        if processed.contains(&v) {
            continue;
        }
        let mut right = v.1;
        for w in g.outgoing.get(&v).into_iter().flatten() {
            right = right.max(if right_codes.contains(&w.0) { w.1 } else { w.0 });
        }
        let mut left = v.0;
        for w in g.incoming.get(&v).into_iter().flatten() {
            left = left.min(w.1);
        }
        let (a0, a1) = (iv.s as i64, iv.e as i64);
        if !(a1 < left || a0 > right) {
            processed.extend(connected_component(g, &[v]));
        }
    }
    processed.into_iter().map(|v| count(g, v)).max().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn components() {
        let mut g = IntronGraph::default();
        let (a, b, c) = (Vertex(100, 200), Vertex(300, 400), Vertex(1000, 1100));
        g.clustered_introns.insert(Iv::new(100, 200), 5);
        g.clustered_introns.insert(Iv::new(300, 400), 7);
        g.clustered_introns.insert(Iv::new(1000, 1100), 2);
        g.outgoing.insert(a, vec![b]);
        g.outgoing.insert(b, vec![Vertex(-10, 500)]);
        g.incoming.insert(b, vec![a]);
        g.incoming.insert(a, vec![Vertex(-21, 50)]);
        g.outgoing.insert(c, vec![]);
        assert_eq!(max_component_coverage(&g, &[Vertex(-21, 50), a]), 7);
        assert_eq!(max_component_coverage(&g, &[c]), 2);
        assert_eq!(max_component_coverage(&g, &[Vertex(-10, 3)]), 0);
        assert_eq!(overlapping_component_max_coverage(&g, Iv::new(450, 460)), 7); // via polyA end 500
        assert_eq!(overlapping_component_max_coverage(&g, Iv::new(40, 60)), 7); // via read start 50
        assert_eq!(overlapping_component_max_coverage(&g, Iv::new(1050, 1060)), 2);
        assert_eq!(overlapping_component_max_coverage(&g, Iv::new(600, 700)), 0);
    }
}
