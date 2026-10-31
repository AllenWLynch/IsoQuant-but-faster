# U15a progress (model filtering)

## State: done (pending final default-build/clippy confirmation)
- `crates/iq_model/src/filter.rs`: pre_filter, filter, detect_similar_isoforms (+ SimilarCache),
  GraphView (live-graph edgeless-vertex emulation on top of iq_graph::query), mapping_quality, 4 unit tests.
- Own `filter::EndProcessor` trait (4 methods filter needs); fixture `ReplayEnds` replays U15b cases by key
  `<sample>:<slice>:<chr>:<block>:<fn>:<model id>.<n>` / `:n.<n>` and checks each call's input (model+reads / store).
- lib.rs: `filter::register_fixtures(reg);`. main added iq_annot dep to iq_model.
- Fixtures quirks41: detect_similar 134/134, filter 67/67, pre_filter 67/67 exact. Default mode: same (no DevIds involved).
- Perf: per-pair assignment cache across both detect calls + skip non-overlapping spans: 2456 vs 6728 assignments (4.1) over fixtures.
- Findings: extra edgeless graph vertices never arise in fixtures (all path introns are graph keys) -> emulation kept, unit-tested.
  chrM hack never exercised by fixtures -> unit test.
