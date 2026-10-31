# U1 progress (annotation parser, AnnotationIndex, IsoformIndex builders)

## Status (parity build, quirks41)
- [x] parse_gff3 + AnnotationIndex (crates/iq_annot/src/{parse,index,attrs}.rs)
- [x] dump_annotation 25/25 exact (`run U1 ../fixtures/annotation/U1.annotation.jsonl.gz`; not found by dir scan, no units/)
- [x] get_gene_info_for_region 67/67, from_models 913/913 exact (`run U1 ../fixtures`); default mode: same, 0 D-R1 hits
- [x] D-R1 unit test (index::tests::d_r1_region_shift)
- [x] exon_ids vs gffutils DB (FeatureIdStorage per chr): 320082/320082 equal (examples/annot_bench.rs 3rd arg)
- [x] 300 extra region cases from gffutils+GeneInfo (scratchpad gen_u1.py, 12/chr, incl. gene-boundary and start=0 probes):
  parity 300/300 exact; default 266 exact + 34 explained D-R1; plain default build 34 mismatches (= explained)
- [x] clippy clean both builds; peak RSS trimmed 155→114 MB; contract proposal #5 written below
- [x] final report (b2_* captures under ../fixtures were still being written, partial gzip; not counted)
- Perf: parse 1.33 s warm (2.35 cold), final heap ≈42.5 MB (isoquant) / 46.8 MB (all), peak RSS 155 MB (transient).

## Facts learned
- gffutils `.attributes` are NOT sorted by sort_attribute_values (that only affects str(feature)); values are split on ','
  unless the value contains ", ", then urllib-unquoted.
- **Child order = gffutils primary key order**: SQLite scans relations by its (parent, child, level) PK, so unordered
  `children()` come sorted by child id string, and `order_by='start'` ties break by child id. Ids: ID attr; duplicate
  IDs → `<ID>_<k>` (create_unique, counter per id); no ID → `<featuretype>_<n>` (counter per type, file order, all rows).
  So `CDS:P_10` < `CDS:P_2`, `exon_10` < `exon_9`. This decides other_features order, CDS/exon ties, transcript ties,
  and children_order_index. Implemented as `DbId` in parse.rs.
- Gene order (features_of_type / region()+stable sort): (start, file order) matched all fixtures.
- GFF3 has only chromosomes 1-25 (no scaffold features).
- `from_models` profiles use delta=0 equal_ranges, like from_genes.
- gffutils region(start=0,..) takes a falsy-start branch: condition `start < end_q` only (reproduced under quirk R1).
- All mode on GRCz11: +6438 ncRNA_gene, +314 pseudogene; 59253 transcripts vs 48837 (isoquant).

## Contract proposal for CONTRACT_CHANGES #5 (gene-less indices) — not applied, for `main`
```rust
// iq_types/src/index.rs
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum IndexKind { #[default] Genes, Models, Region }   // 4.1 constructor the index mirrors
pub struct IsoformIndex { /* ... */ pub kind: IndexKind }
impl IsoformIndex {
    /// 4.1 `get_gene_regions()` is non-empty only for `GeneInfo(gene_db_list, db)`.
    pub fn has_gene_regions(&self) -> bool { self.kind == IndexKind::Genes }
}
// fixture/convert.rs from_dump: kind = match d.kind { "genes" => Genes, "models"|"model" => Models, _ => Region }
```
`genes` stays filled for Models (it doubles as `gene_strands`, read by context.py:199, fl_graph_constructor.py:188 and
the joiner); only `region` is meaningless there. Additive (all builders use `..Default::default()`). U1 builders would
set `kind`; U6 replaces its `IsoRef::Model` heuristic (assigner.rs `has_gene_regions`) with the method.

## Current fixture counts (2026-10-29, CARGO_TARGET_DIR=$HOME/.cache/iqrs-target/U1)
| fn | quirks41 parity | default (quirks41 bin, --mode default) | plain default build |
|---|---|---|---|
| dump_annotation | 25/25 exact | 25 exact | 25 exact |
| get_gene_info_for_region | 67/67 exact | 67 exact, 0 D-R1 | 67 exact |
| from_models | 913/913 exact | 913 exact | 913 exact |
Signatures used by iq_model (keep stable): `index_from_models(chr, &[TranscriptModel], delta)`,
`index_from_models_named(chr, &[TranscriptModel], delta, &[&str])`.

## exon_id rule (U20 question)
4.1 keeps the LAST exon row (file order) per (chr,start,end,strand): `FeatureIdStorage` overwrites a dict while
iterating `region(featuretype="exon")` in rowid order. GRCz11: 12,177 keys have differing exon_ids; 4.1 = last in all
12,177 (first in only 1,528). iq_annot already does this; 320,082/320,082 entries equal. Design §2.2/§8.1 ("first in
file order") should be corrected.
