# U14 progress (assignment-based models) — DONE
- crates/iq_model/src/assign_based.rs: construct_assignment_based(_reads); runner U14Construct.
- crates/iq_model/src/fixture_util.rs (shared U13-U19, per coordinator): StageInput::load/store_out,
  compare_stores (placeholder→4.1 id alignment by creation order), graph/fill_paths loaders, next_novel_gene.
- iq_graph/src/paths.rs: `mod fixture` made pub (coordinator-approved).
- quirks41 parity 67/67 exact; default 67/67 exact (both builds). No DevIds.
- Unit tests (6): known mono-exon gates, non-FL gates (graph, ends, requires_polya), novel mono-exon
  (clustering, overlap rejection, internal reads), cluster/is_internal vs naive ports, next_novel_gene.
- Fixture gap: report_novel_unspliced=false everywhere; no known mono-exon emitted in fixtures.
