# U11 progress (intron graph): done

Status: complete. The quirks41 parity run is 67/67 exact (adjacency key order compared unordered). In `--mode default`, 66 cases are exact and 1 is explained by D-G1. The plain default build
shows 1 mismatch, which equals the D-G1 count.

## Files
- `crates/iq_graph/src/lib.rs`: module wiring and re-exports. It no longer carries `#![allow(unused_variables)]`.
- `collector.rs`: IntronCollector, i.e. collection per distinct chain × multiplicity, the similar-intron window scan, and cluster_introns.
  It also holds the D-G1/D-G2 gates, add_substitute/discard and simplify_correction_map.
- `graph.rs`: `build_graph` / `build_graph_with(GraphParams)`. It covers construct, simplify (tips/bulges, dead ends, isolates)
  and terminal attachment. The fixture runner `U11Graph` lives here, as does the `fixtures_u11` test (IQRS_FIXTURES, plus optional
  IQRS_U11_DUMP for got/want JSON).
- `pyds.rs`: `PyDict` (Python dict order with O(1) del), `VSet` (sorted small set) and `find_closest`.
- `query.rs`: max_component_coverage and overlapping_component_max_coverage (BFS ports, for U15a).

## Findings
- **Adjacency key order (coordinator decision, PySet removed).** 4.1 walks successors in CPython set-hash order, and the
  `signleton_dead_end` defaultdict reads insert `outgoing_edges[v]` keys in that order. With sorted sets, 19/67 cases
  differed **only** in `outgoing`/`incoming` key order. Content was identical (verified on dumps), because every removal
  phase decides before it mutates. The fixture over-specifies this order and no output depends on it: the only readers in
  dict order are `intron_path.py:64-67` (folds the values into a set) and `assignment_based_constructor.py:108-135`
  (existence check). `VSet` is now a sorted vec. The U11 runner's `compare` treats adjacency as unordered maps
  (key → sorted targets); everything else is exact. No DevId is needed.
- A 0-count clustered intron appears in a fixture: `clustered_introns[v]` is read on a collapsed vertex through a stale set
  of an already-collapsed vertex. This is mirrored.
- Performance: collapse_vertex_set scans only the |Δstart| window of a sorted approved vec (a randomized test compares it with
  the literal port). cluster_polya_positions uses a heap and a BTreeMap. Per-read work is limited to position counting, and chain facts
  are cached per chain id.
