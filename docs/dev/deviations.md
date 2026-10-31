# Behaviour relative to IsoQuant 4.1

iqrs follows IsoQuant 4.1.0. Deviation IDs (`D-..`) appear in code comments and fixture
reports.

## IsoQuant 4.1 quirks that iqrs reproduces

These looked like bugs in 4.1, but iqrs keeps them because they change output. Each site has a
"matches IsoQuant 4.1" comment. Before 0.1.0, a test-only `quirks41` build switched each one
off. That switch is gone.

| ID | 4.1 behaviour | site |
|---|---|---|
| D-R1 | gene queries pass the 0-based region to gffutils (1-bp shift; start 0 = no lower bound) | `iq_annot::AnnotationIndex::genes_for_region0` |
| D-R2 | reads that straddle a coverage valley are processed in every block they overlap; the multimapper resolver keeps the best copy | `iq_collect::ClusterOwner` |
| D-P1 | `find_polya` never tests the last window | `iq_read::polya::find_polya` |
| D-P2 | `count_polya_exons` skips non-qualifying exons (no effect on disjoint exons) | `iq_read::polya` |
| D-P3 | boundary deletions are anchored at the post-polyT-trim first exon | `iq_read::process_alignment_into` |
| D-A3 | no common exon: indexes `split_exons[-1]` | `iq_assign::assigner::categorize_exon_elongation_subtype` |
| D-A4 | `score * 1.5 >= best` with a negative best drops every candidate | `iq_assign::assigner::resolve_by_nucleotide_score` |
| D-C1 | `BasicReadAssignment.penalty_score` is always 0 | `iq_collect::multimap::BasicRead::from_view` |
| D-G1 | intron clustering picks the largest coordinate tuple, not the largest count | `iq_graph::collector` |
| D-G2 | the min-count test is skipped for introns whose similar neighbours are all unclustered | `iq_graph::collector` |
| D-M6 | the gene joiner stops while one pair remains (`len(scores) > 1`) | `iq_model::joiner::join_genes_by` |
| D-E1 | mono-exon end correction undoes the start correction | `iq_model::ends::correct_ends` |
| D-O1 | a novel `exon_id` is printed bare (`N`) when created, `chr.N` when reused | `iq_quant::ids::ExonIds::write_id` |

## Remaining differences

| ID | difference | effect |
|---|---|---|
| D-A1 | `--gene-types all` (opt-in) also loads non-`gene` gene features | none with the default `--gene-types isoquant` |
| D-T1 | iqrs breaks ties between equal peaks exactly as numpy's non-SIMD argsort does. 4.1 on CPUs with x86-simd-sort can order large tied arrays differently | rare differences in polyA/TSS peak choice. Fixture runners accept such cases as `explained: D-T1` |
| D-M1 | no ID is drawn for rejected models, so novel transcript/gene numbers follow emitted models only | novel IDs can differ from 4.1's, but the models are the same |
| D-M2 | `detected_known_isoforms` is scoped per work unit, not per process | none observed |
| D-M5 | model store keyed by block read index, not read name | internal only |
| D-O2 | integer read-group tag values are accepted (4.1 crashes) | iqrs runs where 4.1 fails |
| D-O3 | `--grouped-out PATH` overrides the grouped counts file name | opt-in |
| D-O4 | the GTF header version line reads `4.1.0-iqrs-<version>` | header only |

D-P4, D-A5, D-M3, D-M4, D-T2, D-X1: dead or unreachable 4.1 code paths that were not ported or
have no output effect.
