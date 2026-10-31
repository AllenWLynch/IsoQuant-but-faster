# U12 progress (path threading)

Status: done. Parity 134/134 exact (fill 67, known 67); default build identical (no deviations involved).

## Design notes
- `paths::PathProcessor::new` precomputes per intron vertex and side: sorted polyA/polyT and TSS positions, the
  outermost neighbouring intron boundary, and the extreme + next position of the stable-sorted terminal list
  (ties broken by 4.1's concatenation order: read_end/start, tss, tail). thread_ends/starts are O(log k).
- `fill_paths` threads each interned chain once (chain id -> intron-path id, deduped by content since different
  chains can substitute to the same path), keys reads by (path id, start, end), materializes dicts at the end.
- `known_spliced_isoforms`: IndexMap insert = Python dict overwrite (keeps first position).
- `thread_paths` = both. Runners call `fill_paths` / `known_spliced_isoforms` (the fill case has no index).
- Fixture graph is loaded from the capture dir's U11 cases by key (cached per capture dir, static in paths::fixture).
- Tests: randomized comparison against a literal port of thread_ends/thread_starts (2000 graphs, tied positions).
