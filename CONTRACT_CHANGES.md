# Pending contract change requests (iq_types / fixture schema v1.x)

Collected by the coordinator; applied in batches between agent waves (contract is frozen while agents run).

| # | From | Request | Kind | Status |
|---|---|---|---|---|
| 1 | U4 | Add `"pos0": int` to U4 `extract_relevant_deletions` `in` (or `deletions` to `process_alignment` `out`) so D-P3 is directly checkable | schema, additive (needs re-capture of U4) | harness applied (capture batch 2, `b2_*` dirs); Rust DTO pending |
| 2 | U4 | `FixtureCtx::sibling(case_key, fn)` helper to look up a sibling case sharing an input | iq_types runner, additive | queued |
| 3 | iq_types | SCHEMA changelog missing U7 `ref` list + U5 null `read` entries | doc | applied (SCHEMA changelog) |
| 4 | iq_types | AssignParams "missing attr → null" vs non-optional DTO fields (no current fixture has null) | schema/doc | queued |
| 5 | U6a | `IsoformIndex` flag for gene-less (model/region) indices: 4.1 `get_gene_regions()` returns {} when gene_db_list is empty, but `from_dump` fills `genes`. Workaround: treat index as gene-less when keys are `IsoRef::Model` | iq_types, additive | queued |
| 6 | U5 | **Harness bug**: U5 `match_genomic_features` records `delta` = constructor kwarg (0 at exon_corrector call site, hooks.py:494) instead of the comparator`s effective delta. Fix: record `self.comparator.keywords["delta"]`; re-capture U5. Runner workaround: `ctx.params()?.matching.delta` (12802/12802 exact) | schema fix + re-capture | harness fixed (capture batch 2): `b2_*` dirs carry the effective delta; no DTO change; keep workaround for M0a dirs |

## Notes for future unit briefs
- U7 owner: add `correct::register_fixtures(reg);` to `crates/iq_read/src/lib.rs` (U4 pilot made lib.rs delegate per module).

## Performance backlog (apply during integration/M6 tuning)
- U6: cache JunctionComparator's running-max-of-known-intron-ends table per IsoformIndex (e.g. in `AssignScratch`, keyed by index identity) instead of rebuilding O(K) per inconsistent read; move per-read `ReadMatches` into scratch. (U6b)
- U6: `Assignment`/match vectors are still allocated per read (U6a).

## Design doc corrections reported by agents
- §6.3 D-P2: unreachable → verdict NONE (applied).
- §6.4.3: "per-intron lookup gives the same boolean" is not exact — with exons shorter than delta a known intron can be consumed by the previous read intron; U6b kept 4.1 behaviour (randomized test vs full sweep). (not yet edited in 07)

## Notes for U8 (collect integration)
- U7 `correct_read` is generic over `RefBases` (impls: `RefSlice` for FASTA mmap, window list for fixtures); `correct_read_with(exons, rec, a, index, ref, &CorrectParams)` takes flags explicitly.
- U7 `assignment_strand(a, n_exons, corrected_introns, polya, ref)`: `n_exons` = len of *uncorrected* exons (4.1 uses `len(exons)`); introns via `introns_of(&corrected_exons)` (no alloc).
- U4 `process_alignment_into` reuses output buffers; U5 `ProfileBuilder::new(&index,&params).build_into(..)` precomputes per block.
- U7 `groups::group_string`: tag → string, missing → "NA".
- U2 `iq_bam`: read-group tag is set on the handle (`BamOptions{tag, threads}`), not per fetch; tag SAM type and ref id are on `RecordBatch`. Use ~4 inflate threads per fetch or parallelise across blocks (consumer-bound beyond 4). `FastaIndex::slice` clips start<1 to 1 (pyfaidx would negative-index) — U8 must not pass negative starts.
- iq_types doc fix (queued, #7): `RecordLite::end0()` returns pos for mapped reads whose CIGAR consumes no reference; pysam `reference_end` gives pos+1. Such reads have no exons; document it.
- Build note: rust-htslib (feature `htslib`) needs `LIBCLANG_PATH=<repo>/env/lib BINDGEN_EXTRA_CLANG_ARGS=-I/usr/lib/gcc/x86_64-redhat-linux/11/include`.
- U3 `iq_regions`: `GeneOracle::genes_for_region0(Iv0,&mut Vec<GeneIdx>)` — adapter at integration calls `AnnotationIndex::genes_for_region0` (D-R1 gate lives in iq_annot, no extra shift in U3). Records with no reference end (unmapped / no CIGAR) are skipped + counted (4.1 crashes) — design §5.1 text saying unmapped-with-position are included needs correcting. Empty coverage → whole region, cov 0 (4.1 crashes). `pass_a(records, genes, PassAConfig) -> PassAResult{clusters, dup_candidates, stats}`; `overlapping_blocks` = 4.1 double-processing for U8 under quirks41.
- Cleanup (later): iq_regions has a test-only BGZF/BAM decoder in `src/fixture.rs` (uses flate2 as normal dep); switch to `iq_bam` and make flate2 a dev-dep.
- Fixture coverage gap: `requires_polya_for_construction=true` never occurs in fixtures (U12 FL check branch untested). Capture a slice/sample where the sample-wide polyA fraction triggers it (dataset_processor.py:283-300) in the next capture batch.
- Decision (coordinator, per user parity policy): Python set/dict iteration order is NOT emulated. U11 adjacency (`outgoing`/`incoming`) key order is compared unordered — it never reaches outputs (only readers: intron_path.py:64-67 → set; assignment_based_constructor.py:108-135 → existence check). Content-affecting defaultdict inserts are still mirrored. Same rule for later units: if hash order would only change listing order, compare unordered; if it changes content, report to main.
- U11 query.rs (component coverage, for U15a): straight port on the captured graph; 4.1's live graph may gain extra empty defaultdict entries after __init__ — U15a must verify against its fixtures.
- Fixture coverage gap (U14): `--report_novel_unspliced` is false in every capture and no fixture emits a known mono-exon model (only 7/67 U14 cases create anything). Next capture batch: add slices run with `--report_novel_unspliced true` (the user's phased pipeline runs use the default true) and a locus with known mono-exon expression.
- Shared module: `crates/iq_model/src/fixture_util.rs` (owner: U14 agent) — U11/U12 fixture loader, next_novel_gene, creation-order id remap. `iq_graph::paths::fixture::graph` made pub for it.

## From vtg memory prototype (docs/vtg_proto.md) — candidate iq_types v1.1 changes
| # | Finding | Proposed change |
|---|---|---|
| 8 | Read names stored inline in BlockStore (20-28 B + 4 B offset/read, ~20% of resident store); design §2.3 assumed cold/spillable | move names to a cold, spillable arena (hash kept hot) |
| 9 | `BlockStore::heap_bytes()` counts capacity → overstates resident by 1.35-2× | **applied**: `len_bytes()` added (coordinator); shrink_to_fit after collect still open |
| 10 | Boundary deletions not rare (PacBio 0.12-0.2/read, ONT 1.1-1.5/read; up to 12 B/read) | fine as is; update design estimate |
| 11 | Pass A of a chromosome is serial (~33 s for vtg chr22 at 4 inflate threads) before its giant unit can collect | scheduler (U21): overlap pass A of next chromosome with collect; consider streaming blocks out of pass A as they close |
| 12 | ONT raw-chain growth after exon correction unmeasured | measure in M3/M5 |
| 13 | U14 | `next_novel_gene` = 1 + max live Novel key → reuses a number if the max-key model is deleted (harmless in prod; fixture compare would misreport). Add a monotonic `next_novel_gene` counter to `ModelStore` if a later unit hits it | iq_types, additive | queued |
| 14 | U13 | U13 fixture `in` lacks reference windows for novel-intron strand (`get_intron_strand`); runner reads `data/ref/genome.fa` (override `$IQRS_REFERENCE`) — load-bearing (11/67 fail without). Add per-FL-intron `intron_sites` to U13 `in` (like U7) | schema, additive + re-capture | harness applied (capture batch 2, `b2_*` dirs); Rust DTO pending |
| 15 | U13 | `select_reference_gene` ties broken on gene-id string; add `IsoformIndex::gene_id_rank` (workaround: `fl::gene_id_ranks` → `FlAux.gene_rank`) | iq_types, additive | queued |
| 16 | U13 | `LiveAssigner` (ReadAssigner over `iq_assign::assign_with`, thread-local scratch) lives in `fl.rs`; move to `iq_model/src/ctx.rs` for U15a/U17/U19 | refactor | queued (later units use `fl::LiveAssigner` for now) |
- U13 entry: `construct_fl(ms, detected_known, ctx, &FlAux{reference: RefSlice, gene_rank})`.
- Design §7.9 correction (U10b): `TerminalCounter._predict_rows` does not call `detect_peaks`; it has its own pandas pipeline whose np.var/stats.skew sum over read positions in read order (vs histogram order in detect_peaks) — last-bit differences in ~24% of random cases, no prediction change observed. U10b keeps a flat (slot,pos) read-order buffer (8 B/accepted read, per block, reused) to match bits. If memory matters on giant blocks (10M reads → 80 MB transient), switching to histogram order is a candidate DevId (no observed output change).
- Fixture artifact: last block of each chromosome has a second `flush` (`.1`) on an emptied accumulator → `[]` (dump() re-calls flush). Not behaviour to port.
- Cleanup: `genome_path`/`read_region` reference helpers are private in `fl.rs` and copied into `splice_site.rs` fixture module; move them to `fixture_util` (with #14 re-capture they become unnecessary).
- Fixture coverage gap (U17): `resolve_read_ownership_conflicts` — only 7/67 cases have non-FL knowns, 1 changes the store, none drops a model (unit-tested). Next capture batch: loci with competing non-FL known models.
- U17 signatures: `assign_reads_to_models(ms, ctx)` (ctx.store = BlockStore), `resolve_ownership(&mut ModelStore)`.
- Fixture coverage gap (U15a): chrM pre-filter hack (>50 models with ≤2 exons) and 4.1 live-graph empty-vertex case never occur; unit-tested only.
- U15a defines `filter::EndProcessor` trait (4 calls used by `filter`); U15b's `Ends` must implement it.
- Perf backlog: U15a detect_similar pair loop serial (rayon over outer models possible); already 2456 vs 4.1's 6728 assignments via pair cache.
- U18 note: on novel+novel gene merges 4.1 keeps the gene with the larger id string; iqrs (no id strings at U18, D-M1) uses creation order in production (`join_genes`), 4.1 string order in the fixture runner (`join_genes_by`). Gene groupings identical in all 67 cases; only the surviving label differs (renamed at output anyway). Theoretical risk: tie between equal-scoring pairs broken differently — watch in M5 E2E. D-M6 changes 3/67 cases, all novel+novel merges (design §7.10 example is novel+annotated; update text).
- Fixture coverage gap (U15b): alt-end derivation thin — 2 derive cases produce models, 1 add-known case changes the store, `novel_apa` never on. Next capture batch: loci with known-model alternative polyA/TSS ends.
- D-E1 added (mono-exon end edit undoes start edit; regression fe84ccdc/0d6e47f2) — iq_types DevId + design §9 row.
| 17 | U19 | `ReadAssigner` trait can't pass a per-call junction delta; single-model path needs `delta = max(delta, graph_clustering_distance)` (`Params::final_single_delta`), so `final_assignments` builds its own `fl::LiveAssigner`. Add optional `MatchingParams` override to `ReadAssigner::assign`, or put assigner delta on `IsoformIndex` | iq_model trait | queued |
- U19 API: `final_assignments(ms, ctx, sink: impl FnMut(u32, FinalAssign))` emits in read order (store from `ctx.store`); `compact(&Assignment) -> FinalAssign`.
- Design §2.2/§8.1 correction (U1): exon_id map keeps the **last** exon row in file order (4.1 FeatureIdStorage overwrites); verified 12,177/12,177 conflicting keys on GRCz11. iq_annot already does this.
- U1 gffutils behaviour: child order = SQLite primary-key (child id string) order; ties under order_by='start' break by child id (`CDS:P_10` < `CDS:P_2`). Attribute values not sorted (sort_attribute_values only affects printing).
- U1 extra D-R1 cases + exon_id check: generators in `crates/iq_annot/tests/gen/`, outputs in `../fixtures/annotation/extra_u1/` (moved from scratchpad).
- Open (U1): model-index id_rank uses creation order (production) vs real id strings (fixtures) — U17/U19 used `_named` in runners; confirm production consumers only use ranks as sets (watch in M5).

- U8/U21 integration: RefSlice per block must cover `[min(cluster.s, index.start)-1000, max(cluster.e, index.end)+1000]` (U7 strand reads index-intron sites on TerminalExonMisalignment); missing bases silently → non-canonical. Read-group ids are first-seen (scheduling-dependent) → writer must order groups by name / remap.
- U9 landed: `resolve_chr_duplicates` removes suspended reads and densely re-indexes stores (== 4.1 loader view); validated on 1701 blocks/22 captures; total/polya == 4.1 log lines.
- Capture batch 2 done (22 `b2_*` dirs, 363 MB): gaps (a)-(f) closed; #1, #3, #6, #14 applied in harness. Coordinator added the 3 optional DTO fields (`U4ExtractDeletionsIn.pos0`, `StageStoreIn.intron_sites`, `StageStoreIn.ref_`) → `iqrs-fixture check` 1,736,452 lines / 5,865 files, 0 non-conforming. **All units on all b2 dirs (quirks41): 0 mismatches, 0 errors** (only D-T1-explained peak ties). Runner workarounds for #1/#6/#14 can now be removed when convenient (old dirs still need them).
- D-T1 revised (user-approved 2026-10-29): `iq_terminal::peaks` now emulates numpy's non-SIMD argsort (`np_argsort`, `Ties::Numpy` default; `Ends::default()` too) in `select_by_peak_distance` and `rank_peaks`; deterministic, no quirk gate. `Ties::PositionAsc/Desc` remain for the D-T1 comparators. `u10a_random_stable` regenerated under `NPY_DISABLE_CPU_FEATURES="AVX512_SPR AVX512_ICL X86_V4 X86_V3"` (gen_u10a.py re-execs itself), bit-exact. Fixtures (quirks41 and default): U10a detect_peaks 215 exact + 66 D-T1 (was 164 + 117), U10b flush 361 + 7 D-T1 (was 364 + 4; SIMD-captured large arrays), U15b all exact. chr8: parity == 4.1 no-SIMD on Iso-Seq 720 and Revio 726; Revio 716 == 4.1. Design §9 D-T1 row updated.

## D-R2 retired (2026-10-29, user-approved; evidence `../validation/deviation_effects.md`)
- `DevId::R2` verdict DROP → **KEEP** (`iq_types::dev`). Both builds now process a valley-straddling read in every block whose 4.1 fetch window returns it, and U9 keeps the best copy (4.1 behaviour). No `quirk(DevId::R2)` call site remains; the variant stays so `IQRS_QUIRKS` specs and labels keep parsing.
- Why: the single-owner rule (max overlap) cost reads their gene. ONT SRR33046254 had 1.22 pp of reads move assigned→no_feature (`seta` lost 289/748), and Revio got read-through models that merged loci (tx/locus +1.4 % on 716, +5.1 % on 726). Default + R2 + M6 passed 15/15 gates on all 5 chr8 samples.
- API changes (additive except one signature):
  - `iq_collect::Owner::shared(&rec) -> bool` (default `false`). `ClusterOwner::owns` is now always the fetch-window test, and `shared` is true when the record also overlaps an adjacent block region.
  - `BlockCollector::finish_shared() -> (BlockStore, CollectStats, Vec<u64>)` returns the name hashes of selected shared reads.
  - **`multimap::resolve_chr_duplicates(stores, candidates: Option<&[u64]>, strategy)`** (was `&[u64]`, with "empty = all"). `Some(c)` must cover the pass-A duplicate candidates ∪ the shared names; `None` scans every read.
  - `iq_types::Block::owned_reads` doc: the cost estimate now adds each boundary bin's coverage to both neighbouring blocks (pass A).
  - `iq_sched::real::RealStore::shared`: the barrier merges it into the U9 candidates.
- Design: §5.2 rewritten, §5.1/§6.2/§9 rows updated (D-R2, D-C1, D-M5).
