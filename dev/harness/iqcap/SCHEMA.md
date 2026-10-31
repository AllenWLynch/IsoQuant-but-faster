# iqcap fixture schema (v = 1)

Contract between the Python capture harness `iqcap` (writer) and the Rust crate `iq_types`
(reader). Every fixture is produced from IsoQuant 4.1.0 (commit `cea5337`) run with `--threads 1`.
Source paths `foo.py:N` are relative to `IsoQuant/isoquant_lib/`.

Status: **v1 draft**. Drift rule: every field iqcap emits is documented here. The Rust round-trip test
fails on any undocumented field. That failure is intentional: it is the drift alarm. Breaking changes bump `v` (in every case line and in `fixture_manifest.json`)
and are announced to the maintainers. Additive changes (new optional field, new `fn` value) do not
bump `v`; readers must ignore unknown fields.

---------------------------------------------------------------------------------------------------
## 0. General conventions

### 0.1 Files and directory layout

One capture = one directory `fixtures/<sample>/<slice>/`:

```
fixture_manifest.json            run metadata (§1)
params.json                      effective IsoQuant args (§2.3)
units/<UNIT>.jsonl.gz            one case per line (§0.2); UNIT in U1,U3..U20 (U10a, U10b, U15a, U15b spelled so)
index/<sha1>.json                IsoformIndex dumps, content addressed (§2.1)
block/<sample>/<chr>/<block_no>.jsonl.gz   construct-side read records per block (§2.2)
outputs/                         copies of the IsoQuant final outputs of this run (U20 end-to-end)
slice/                           (optional) the BAM slice used, if produced by iqcap.slice
```

Separate, run-independent dumps:
```
fixtures/annotation/U1.annotation.jsonl.gz     iqcap.dump_annotation (U1 hook a)
fixtures/<sample>/u2/units/U2.jsonl.gz         iqcap.dump_bam (U2), records of one region of that sample's `ordinary` slice BAM (HP:Z synthesized)
```

- All `*.jsonl.gz` are gzip-compressed JSON Lines, UTF-8, one JSON object per line, `\n` terminated.
- With more than one writer process (never with `--threads 1`) a unit file is `units/<UNIT>.<pid>.jsonl.gz`;
  readers must glob `units/<UNIT>*.jsonl.gz`. A single-writer run renames to `units/<UNIT>.jsonl.gz`.
- Within one unit file, cases appear in **call order** of the hooked function in the Python run
  (deterministic with `--threads 1`).

### 0.2 Case envelope

```json
{"unit":"U6","v":1,"case":"SRR33046254:ordinary1:8:12:a37","in":{...},"out":{...}}
```

| field | type | meaning |
|---|---|---|
| `unit` | string | unit id: `U1`..`U20`, `U10a`, `U10b`, `U15a`, `U15b` |
| `v` | int | schema version (this document: 1) |
| `case` | string | deterministic case key (§0.3), unique within the unit file |
| `in` | object | input payload. **Always** contains `"fn"` (string): the hooked Python function, so units with several hooks can be filtered (`"fn":"assign_to_isoform"`, ...) |
| `out` | object | output payload (or the exception, §0.8) |

Key order of every object is the order listed in this document (Python dict insertion order).

**Serialization (pinned).** Every case line, block line and `index/*.json` file is produced by one shared
helper, `iqcap.sink.dumps`:
`json.dumps(obj, separators=(',', ':'), ensure_ascii=True, allow_nan=False, sort_keys=False)`
(so non-ASCII characters appear as `\uXXXX` escapes, there is no whitespace, and keys are in insertion order).
`params.json` uses the same call with `sort_keys=True`, followed by one trailing `\n`.
`fixture_manifest.json` is metadata (written with `indent=1`) and is not part of the round-trip contract.
Case lines and block lines end with `\n`; `index/*.json` has no trailing newline. The SHA-1 is computed over exactly those bytes.

### 0.3 Case keys and identities

`case = "<sample>:<slice>:<chr>:<block_no>:<local>"`

- `sample`, `slice`: as given to `iqcap.run --sample/--slice`.
- `chr`: chromosome name (`"8"`).
- `block_no`: 0-based index of the block (= one `process_alignments_in_region` call = one GeneInfo
  record in the `.save` dump = one `ReadAssignmentLoader.get_next` call) within the chromosome, in
  4.1 processing order. `-` when the hook fires outside a block (e.g. U3 `genes_in_region`).
- `local` = `<fn>:<base>.<k>`, unique within the unit file:
  - `fn`: the `in.fn` value;
  - `base`: `a<aln_no>` (collect-side per-read cases; `aln_no` = 0-based index of the alignment among
    the `AlignmentInfo` objects built in this block, i.e. the read identity at collect time),
    `r<read_idx>` (construct-side per-read cases), `id:<read_id>` (assigner calls on provisional
    models: `fl`/`similar`), a model id (U15b/U16 per-model calls), `polya`/`tss` (U10b), a read name
    (U9), or `n` otherwise;
  - `k`: 0-based counter of (unit, block, fn, base), so repeated calls get `.0`, `.1`, ...
- Exceptions: U9 uses `<sample>:<slice>:<chr>:-:resolve:<read_id>.<k>`. U3 `cluster`/`genes_in_region`
  and U20 `counter_dump` use block `-` (`...:-:counter_dump:<counter file basename>`).
- Example keys: `S:valley:8:0:assign_to_isoform:a37.0`, `S:valley:8:0:assign_to_isoform:r12.1`,
  `S:valley:8:0:construct_fl_isoforms:n.0`.
- Readers must treat keys as opaque strings except for the documented prefix
  `<sample>:<slice>:<chr>:<block_no>`. Cross references (`graph`, `paths`) hold full case keys.

**Read identity.**
- Collect side: `in.read = {"qname": str, "aln_no": int}`.
- Construct side: `read_idx` (int) = position of the ReadAssignment in the list returned by
  `ReadAssignmentLoader.get_next` for that block, i.e. the line number (0-based, after the header)
  in `block/<sample>/<chr>/<block_no>.jsonl.gz`. iqcap tags each object with `_iq_idx` in a
  wrapper around `ReadAssignmentLoader.get_next` (`assignment_loader.py:83`).
- Model identity: the 4.1 `transcript_id` string (e.g. `transcript17.8.nnic`, or a reference id for
  known models). The Rust runner maps ids to `ModelIdx` by creation order.

### 0.4 JSON value encodings

| Python | JSON |
|---|---|
| int | integer (exact; may exceed 32 bits, e.g. `2147483647`, `1073741823`, `2147483648`; use i64) |
| float | number written by Python `json` = `float.__repr__` (shortest round-trip; `0.0`, `0.1`, `1e-05`, `1e+16`, `-0.5`). **Floats are always written with a `.` or an exponent** (a float that happens to be integral is `3.0`, never `3`). Parse as f64. |
| non-finite float | the **strings** `"NaN"`, `"Infinity"`, `"-Infinity"` (only where a field is documented as `f64*`) |
| bool | `true`/`false` |
| None | `null` |
| str | string |
| tuple/list | array |
| interval `(s, e)` | `[s, e]` (`Iv`) |
| dict with string keys | object **only where noted**; otherwise and always for non-string keys: array of `[key, value]` pairs **in Python dict iteration order** (insertion order; this order is load-bearing in 4.1) |
| set | array, **sorted** (Python sort order of the elements), unless noted |
| Enum | its integer `.value` (IsoQuant codes, §3) |
| strand | one-char string `"+"`, `"-"`, `"."` |

Float repr details, for writers that must reproduce bytes (round-trip tests): Python `repr(float)`
uses the shortest digit string that round-trips; positional notation when `-4 <= exp10 < 16`,
otherwise scientific `d[.ddd]e±XX` with at least two exponent digits (`1e-05`, `1.5e+16`); integral
values in positional notation end in `.0`. Type notation below: `f64` = always float; `f64*` = float
or non-finite string.

Penalty scores (`IsoformMatch.penalty_score`) are Python ints (0) or floats in 4.1; iqcap always writes
them as floats (`0.0`).

### 0.5 Coordinates

Unless noted, coordinates are IsoQuant's: **1-based, closed** genomic intervals (`Iv`).
Exceptions (each marked in the field tables):
- `pos0`: BAM 0-based leftmost position (`alignment.reference_start`).
- `ref_end0`: `alignment.reference_end` (0-based, exclusive).
- U3 storage regions / coverage: `AbstractAlignmentStorage.region` is `(reference_start,
  reference_end - 1)` (**0-based closed**); coverage bins are `pos0 // 256`.
  `split_coverage_regions` output regions are mixed (`bin*256+1` clipped by the 0-based input region)
  exactly as 4.1 computes them. `process_alignments_in_region.region` is whichever of the two 4.1
  passed. All recorded verbatim.
- `read_blocks`: 0-based closed positions in the query sequence (`get_read_blocks`).
- `cigar_blocks`: closed index ranges into the CIGAR op list.
- Reference sequence windows (`RefWin`, §2.4): `start0` is a 0-based offset into the chromosome.

### 0.6 Ordering guarantees

- Lists of reads, models, matches, events, exons, profile entries: 4.1 order, never re-sorted.
- `[key, value]` pair lists: Python dict iteration order of the underlying dict at capture time.
- Sets: sorted ascending (Python order; tuples compare lexicographically, strings by code point).
  Exception: where 4.1 iterates a set in its own (hash) order and that order matters, the field is
  documented as "set iteration order".

### 0.7 Sampling

`iqcap.run --sample-rate` takes `R` or `UNIT=R,...,default=R`. Selection uses
`keep(key, R) = int(sha1(key).hexdigest()[:16], 16) / 2**64 < R` (always true for R >= 1):
- **Blocks.** A block (key `<sample>:<slice>:<chr>:<block_no>`) is selected iff `keep(block key, --block-rate)` holds and
  its loader read count is <= `--block-cap` (default 200000). Blocks over the cap are listed in the manifest's
  `skipped_blocks`. For a selected block, all of U8, U10b-U19, U20 `gff_dump` and the `block/` file are written.
  Nothing is sampled within a block, so every cross reference resolves.
- **Per-read / per-call units** (U1 b/c, U3, U4-U7, U9, U10a, U20 `counter_add_read_info`): `keep(k, rate(unit))`.
  `k` is a deterministic function of the case identity: `<block key>:a<aln_no>` for U4 and U7 (so equal rates select
  the same reads), extended with the call site / function for U5/U6. Construct-phase cases (those with call
  sites `fl`, `similar`, `reassign`, `final_*`, U1c, U10a) are written only inside selected blocks.
  For a complete ReplayAssigner record of a block, use `U6=1`. Check `iqcap_args.sample_rate` in the manifest.

Every case is complete on its own (all inputs inline or via the side tables); sampling never leaves
a dangling reference except as stated above.

### 0.8 Exceptions

If a hooked call raises, iqcap writes `out = {"exception": "<type>: <message>"}` and re-raises.

---------------------------------------------------------------------------------------------------
## 1. `fixture_manifest.json`

```json
{
 "v": 1,
 "tool": "iqcap",
 "created": "2026-10-29T12:00:00Z",
 "sample": "SRR33046254", "slice": "ordinary1",
 "isoquant": {"version": "4.1.0", "commit": "cea5337ea4e7...", "path": ".../IsoQuant", "dirty": false},
 "versions": {"python": "3.x.y", "numpy": "...", "scipy": "...", "xgboost": "...", "pandas": "...", "pysam": "...", "gffutils": "..."},
 "cpu": {"model": "...", "avx512f": true, "avx2": true, "flags_hash": "<sha1 of the flags line>"},
 "host": "...",
 "slice_def": {"bam": "...source bam...", "regions": ["8:1000-2000"], "hp_synthesis": "sha1(qname)%3 -> A|B|unphased", "slice_bam": "...", "n_records": 1234},
 "isoquant_args": ["--reference", "...", ...],
 "iqcap_args": {"units": [...], "sample_rate": {...}, "block_rate": 1.0, "block_cap": 200000},
 "units": {"U6": {"files": ["units/U6.jsonl.gz"], "cases": 1234, "by_fn": {"assign_to_isoform": 1234}}},
 "blocks": [{"chr": "8", "block_no": 0, "region": [s, e], "n_reads": 10, "captured": true}],
 "skipped_blocks": [{"chr": "8", "block_no": 3, "reason": "block_cap", "n_reads": 250000}],
 "wall_seconds": 12.3
}
```
`slice_def` may also carry `requested`, `subsample`, `subsample_rule` and (batch 2, additive) `keep_names`:
`{"file", "n", "sha1" (of the sorted names joined by `\n`), "rule"}` when the slice was restricted to a read-name
list (`iqcap.slice --keep-names`, used by `b2_polyarich`).

---------------------------------------------------------------------------------------------------
## 2. Shared side tables and shared object types

### 2.1 `index/<sha1>.json` — IsoformIndex (a dumped `GeneInfo`)

Plain JSON (not gzipped), one object. `<sha1>` = SHA-1 (hex) of the exact file bytes. Written once
per distinct content. Referenced from cases as `"index": "<sha1>"`.

| field | type | source |
|---|---|---|
| `kind` | string | `"genes"` (`GeneInfo.__init__` from gffutils genes), `"models"` (`from_models`), `"model"` (`from_model`), `"region"` (`from_region`, no genes) |
| `chr` | string\|null | `chr_id` (`null` for `from_models([])`) |
| `start`, `end` | int | `start`, `end` (gene-set region, or model span) |
| `delta` | int | `delta` |
| `genes` | [string] | `gene_db_list` ids, in list order (sorted by gene start, gffutils order for ties) ; `[]` for non-`genes` kinds |
| `gene_regions` | [[string, Iv]] | `get_gene_regions()` (gene → (start,end)), `gene_db_list` order; `[]` for non-`genes` kinds |
| `isoform_exons` | [[string, [Iv]]] | `all_isoforms_exons` in dict order (genes by start, transcripts by `children(order_by='start')`; models: storage order) |
| `isoform_introns` | [[string, [Iv]]] | `all_isoforms_introns` in dict order |
| `isoform_strands` | [[string, strand]] | `isoform_strands` in dict order (**note**: for `genes` this dict is filled from `children(featuretype=DEFAULT_RNA_TYPES)` without `order_by`, so its order may differ from `isoform_exons`; and it can contain transcripts that have no exons) |
| `gene_strands` | [[string, strand]] | `gene_strands` in dict order |
| `gene_id_map` | [[string, string]] | `gene_id_map` (transcript → gene) in dict order |
| `intron_features` | [Iv] | `intron_profiles.features` (sorted) |
| `intron_profiles` | [[string, [int]]] | `intron_profiles.profiles` (transcript → list of -2/-1/1, one per feature), dict order |
| `intron_ranges` | [[string, [int,int]]] | `intron_profiles.profile_ranges` (half-open index range) |
| `exon_features`, `exon_profiles`, `exon_ranges` | same | `exon_profiles` |
| `split_exon_features`, `split_exon_profiles`, `split_exon_ranges` | same | `split_exon_profiles` |
| `sources` | [[string, string]] | `sources` (gene and transcript ids → GFF source) in dict order |
| `other_features` | [[string, [[int,int,string]]]] | `other_features` (transcript → `(start, end, featuretype)` CDS/start_codon/stop_codon/UTR, gffutils children order), dict order |
| `feature_attributes` | [[string, string]] | `feature_attributes` (gene id, transcript id, and the quirky `<tx>_<s>_<e>_<strand>` exon keys of `set_gene_attributes`, `gene_info.py:576-596`) → preformatted `k "v"; ` string, dict order. `[]` for non-`genes` kinds. |

Profile values: `1` feature present, `-1` absent but inside the transcript span, `-2` outside the
transcript span (`FeatureProfiles.set_profiles`, `gene_info.py:96-118`).

### 2.2 `block/<sample>/<chr>/<block_no>.jsonl.gz` — construct-side read records

Line 0 is a header, lines 1..N are `ReadRecord`s.

Header:
```json
{"kind":"block_header","v":1,"sample":"...","slice":"...","chr":"8","block_no":12,
 "region":[s,e],"index":"<sha1>","gene_ids":["ENSDARG..."],"n_reads":N,
 "collect_region":[s,e],"skip_read_fraction":1}
```
- `region`: `(gene_info.start, gene_info.end)` of the loader's GeneInfo (gene-set span, or the
  `from_region` region for gene-less blocks).
- `index`: the loader's GeneInfo (`ReadAssignmentLoader.get_next`) — this is the `gene_info` every
  construct-side unit of the block uses.
- `collect_region`: the `process_alignments_in_region` region of the same block (§0.5), or `null`
  if the collect pass was not observed (should not happen).
- Reads: exactly the list returned by `ReadAssignmentLoader.get_next` (post multimapper resolution
  `assignment_loader.py:50-70`, suspended reads already dropped), in that order; `idx` = 0..N-1.

`ReadRecord` (also used inline by U8):

| field | type | source (`ReadAssignment`) |
|---|---|---|
| `idx` | int | block read index (loader order; in U8 `out.reads`: collect order) |
| `assignment_id` | int | `assignment_id` (global counter; multimapper join key for U9) |
| `name` | string | `read_id` (QNAME) |
| `chr` | string | `chr_id` |
| `region` | Iv-like [int,int] | `genomic_region` (the `process_alignments_in_region` region, §0.5) |
| `exons` | [Iv] | `exons` (after polyA/T trimming) |
| `corrected_exons` | [Iv] | `corrected_exons` |
| `boundary_deletions` | [[int,int]] | `boundary_deletions` as (1-based ref pos of first deleted base, length) |
| `polya` | [int,int,int,int] | `polya_info` = `[external_polya_pos, external_polyt_pos, internal_polya_pos, internal_polyt_pos]`, `-1` = absent |
| `polya_found` | bool | `polyA_found` |
| `cage_found` | bool | `cage_found` |
| `strand` | strand | `strand` |
| `mapped_strand` | strand | `mapped_strand` |
| `mapq` | int | `mapping_quality` |
| `multimapper` | bool | `multimapper` |
| `read_group` | [string] | `read_group` (group names, one per grouping strategy; `[]` without `--read_group`) |
| `assignment` | Assignment | §2.5 (`assignment_type`, `gene_assignment_type`, `isoform_matches`) |
| `introns_match` | bool | `introns_match` |
| `additional_attributes` | object | `additional_attributes` (str→str; BAM tags imported with `--bam_tags`; usually `{}`) |

`corrected_introns` is not stored: it is always `junctions_from_blocks(corrected_exons)`.

### 2.3 `params.json`

All attributes of the final IsoQuant `args` namespace (after `set_additional_params`) whose values
are JSON-encodable scalars or lists thereof; Enums by `.value` plus `"<name>__name"` with the enum
member name; other objects as `repr` strings. Keys sorted. Includes in particular `delta`,
`apa_delta`, `minor_exon_extension`, `major_exon_extension`, `min_abs_exon_overlap`,
`min_rel_exon_overlap`, `minimal_exon_overlap`, `minimal_intron_absence_overlap`,
`micro_intron_length`, `max_fake_terminal_exon_len`, `max_intron_shift`, `max_missed_exon_len`,
`max_intron_abs_diff`, `max_intron_rel_diff`, `max_suspicious_intron_abs_len`,
`max_suspicious_intron_rel_len`, `resolve_ambiguous`, `correct_*`, `polya_window`,
`polya_fraction`, `min_novel_intron_count`, `graph_clustering_distance`, `graph_clustering_ratio`,
`terminal_position_abs/rel`, `terminal_internal_position_rel`, `min_known_count`,
`min_novel_count`, `min_novel_count_rel`, `min_mono_count_rel`, `min_mono_exon_coverage`,
`singleton_adjacent_cov`, `require_monointronic_polya`, `require_monoexonic_polya`,
`requires_polya_for_construction`, `report_canonical_strategy`, `report_novel_unspliced`,
`simple_alignments_mapq_cutoff`, `simple_models_mapq_cutoff`, `inconsistent_mapq_cutoff`,
`min_mapq`, `fl_data`, `fl_only`, `data_type`, `model_construction_strategy`.

Units whose function reads `self.params`/`self.args` use these values unless the case carries its
own `params` object (U5/U6/U7 always do, because the assigner is sometimes built with a modified
`delta`, `model_read_counter.py:105-108`).

`AssignParams` (inline `params` of U5/U6/U7; the attributes of the `params` object the
assigner / profile constructor / corrector was built with):
`delta, minor_exon_extension, major_exon_extension, min_abs_exon_overlap, min_rel_exon_overlap,
minimal_exon_overlap, minimal_intron_absence_overlap, micro_intron_length,
max_fake_terminal_exon_len, max_intron_shift, max_missed_exon_len, max_intron_abs_diff,
max_intron_rel_diff, max_suspicious_intron_abs_len, max_suspicious_intron_rel_len,
resolve_ambiguous, apa_delta, count_exons, correct_fuzzy_junctions, correct_intron_shifts,
correct_skipped_exons, correct_fake_terminal_exons, correct_terminal_exons,
correct_microintron_retention` (object, keys in this order). Enum-typed args are written as their
integer `.value` (§0.4): `resolve_ambiguous` is an `AmbiguityResolvingMethod` (§3.5), e.g.
`monoexon_and_fsm` -> `20`. Bools are bools, ints are ints, and a missing attribute is `null`.

### 2.4 Reference windows `RefWin`

`{"start0": int, "seq": string}`: `seq = chr_record[start0 : start0 + len(seq)]` (0-based,
half-open), bases **verbatim from the FASTA** (case preserved; 4.1 compares read bases to
`chr_record[...]` case-sensitively in `AlignmentInfo.get_error_count` and upper-cases in
`get_intron_strand`). Windows are clipped at the chromosome ends.

### 2.5 Assignment, IsoformMatch, MatchEvent

`Assignment` (a `ReadAssignment`'s classification):
```json
{"atype": 1, "gene_atype": 1, "matches": [IsoformMatch, ...]}
```
| field | type | source |
|---|---|---|
| `atype` | int | `assignment_type.value` (§3.1) |
| `gene_atype` | int | `gene_assignment_type.value` |
| `matches` | [IsoformMatch] | `isoform_matches`, in order |

`IsoformMatch`:
```json
{"gene": "ENSDARG...", "tx": "ENSDART...", "strand": "+", "class": 10, "penalty": 0.0, "events": [MatchEvent, ...]}
```
| field | type | source |
|---|---|---|
| `gene` | string\|null | `assigned_gene` |
| `tx` | string\|null | `assigned_transcript` |
| `strand` | strand | `transcript_strand` |
| `class` | int | `match_classification.value` (§3.2) |
| `penalty` | f64 | `penalty_score` |
| `events` | [MatchEvent] | `match_subclassifications`, in order |

`MatchEvent`: `{"t": int, "iso": [int,int], "read": [int,int], "info": int}` =
`event_type.value` (§3.3), `isoform_region`, `read_region`, `event_info`. Regions are index ranges
into the isoform's / read's intron or exon lists, or the special constants (§3.4).

### 2.6 TranscriptModel `Model`

| field | type | source (`TranscriptModel`) |
|---|---|---|
| `id` | string | `transcript_id` |
| `gene` | string | `gene_id` |
| `chr` | string | `chr_id` |
| `strand` | strand | `strand` |
| `type` | int | `transcript_type.value`: known=1, novel_in_catalog=2, novel_not_in_catalog=10 |
| `exons` | [Iv] | `exon_blocks` |
| `intron_path` | [[int,int]] | `intron_path` (graph vertices, §2.8; `[]` for known/assignment-based models) |
| `polya_confirmed`, `tss_confirmed` | bool | |
| `source` | string | `source` (GFF column 2) |
| `other_features` | [[int,int,string]] | `other_features` |
| `additional_info` | [[string,string]] | `additional_info` (OrderedDict) |

### 2.7 ModelStore snapshot `Store`

Snapshot of the per-block `ModelStore` (`context.py:48-110`) plus the shared id state.

| field | type | source |
|---|---|---|
| `models` | [Model] | `transcript_model_storage`, in order |
| `transcript_read_ids` | [[string, [int]]] | `transcript_read_ids` (model id → bound reads), dict order; reads as `read_idx` in list order (duplicates possible) |
| `internal_counter` | [[string, int]] | `internal_counter`, dict order |
| `read_assignment_counts` | [[string, int]] | `read_assignment_counts`, keyed by **read name** as in 4.1, dict order |
| `construction_assignment` | [[string, string]] | read name → model id, dict order |
| `read_assignment_counts_idx` | [[int, int]] | the same as `read_assignment_counts`, re-keyed by `read_idx`: for each name entry (dict order), one pair per block read carrying that name (ascending `read_idx`). Duplicate names in a block therefore share 4.1's single name-keyed value (see D-M5) |
| `construction_assignment_idx` | [[int, string]] | the same re-keying of `construction_assignment` |
| `nonfl_known_reads` | [[string, [int]]] | known id → `read_idx` list, dict order |
| `detected_known_isoforms` | [string] | class-level set `ModelStore.detected_known_isoforms` (shared across all blocks of the run), sorted |
| `id_counter` | int | `ctx.id_distributor.value` (the last transcript/gene number handed out; next id = next value > this not in `forbidden_ids`) |

Note: `transcript_read_ids`, `internal_counter` and `read_assignment_counts` are `defaultdict`s; 4.1
inserts zero entries on read access. Snapshots record the dict as it is (zero entries included).

### 2.8 Intron-graph vertices and paths

A vertex is `[a, b]`: an intron `(start, end)` with `a >= 0`, or a terminal vertex `(code, pos)`
with `code < 0`: polya=-10, read_end=-11, tss_right=-12 (genomic right); polyt=-20,
read_start=-21, tss_left=-22 (genomic left) (`intron_graph.py:18-31`). A path is `[[a,b], ...]`.

`Graph` (state of `IntronGraph` right after `__init__`):

| field | type | source |
|---|---|---|
| `clustered_introns` | [[Iv, int]] | `intron_collector.clustered_introns`, dict order |
| `intron_correction_map` | [[Iv, Iv]] | `intron_collector.intron_correction_map`, dict order |
| `discarded_introns` | [Iv] | `intron_collector.discarded_introns`, sorted |
| `outgoing` | [[vertex, [vertex]]] | `outgoing_edges`, dict order, each target set sorted |
| `incoming` | [[vertex, [vertex]]] | `incoming_edges`, dict order, each source set sorted |
| `starting_known_positions` | [[Iv, [int]]] | dict order, list order |
| `terminal_known_positions` | [[Iv, [int]]] | dict order, list order |
| `max_coverage` | int | `max_coverage` |

Later lazy insertions of empty entries by 4.1 read accessors (`outgoing_edges[v]` on a defaultdict)
are not captured; they never add edges.

### 2.9 Block reference

`"block": {"sample": s, "slice": sl, "chr": c, "block_no": n}` — resolves to
`block/<s>/<c>/<n>.jsonl.gz` of the same fixture directory. Read sets passed to a function are given
as `[read_idx]`.

---------------------------------------------------------------------------------------------------
## 3. Enum codes (IsoQuant's own values)

### 3.1 `ReadAssignmentType` (`isoform_assignment.py:40`)
noninformative=0, unique=1, unique_minor_difference=2, ambiguous=10, intergenic=20,
inconsistent=30, inconsistent_non_intronic=31, inconsistent_ambiguous=32, inconsistent_genic=33,
inconsistent_multigenic=34, discarded=40, suspended=255.

### 3.2 `MatchClassification` (`isoform_assignment.py:84`)
undefined=0, full_splice_match=10, incomplete_splice_match=11, mono_exon_match=12,
novel_in_catalog=20, novel_not_in_catalog=21, genic=30, genic_intron=31, intergenic=32,
antisense=40, fusion=50.

### 3.3 `MatchEventSubtype` (`isoform_assignment.py:129-210`)
none=0, undefined=1, antisense=9, mono_exonic=10, fsm=11, ism_right=13, ism_internal=14, ism_left=15,
mono_exon_match=19, intron_shift=21, exon_misalignment=22, exon_elongation_right=23,
fake_terminal_exon_left=24, exon_elongation_left=25, fake_terminal_exon_right=26,
terminal_exon_misalignment_right=28, terminal_exon_misalignment_left=29, intron_retention=31,
unspliced_intron_retention=32, fake_micro_intron_retention=34, incomplete_intron_retention_left=38,
incomplete_intron_retention_right=39, alt_left_site_known=101, alt_right_site_known=102,
alt_left_site_novel=103, alt_right_site_novel=104, intron_migration=114,
intron_alternation_novel=115, intron_alternation_known=116, mutually_exclusive_exons_novel=121,
mutually_exclusive_exons_known=122, exon_skipping_known=133, exon_skipping_novel=134,
exon_merge_known=135, exon_merge_novel=136, exon_gain_known=145, exon_gain_novel=146,
exon_detach_known=147, exon_detach_novel=148, terminal_exon_shift_known=151,
terminal_exon_shift_novel=152, alternative_structure_novel=181, alternative_structure_known=182,
alternative_polya_site_right=200, alternative_polya_site_left=201, internal_polya_left=202,
internal_polya_right=203, alternative_tss_left=204, alternative_tss_right=205,
correct_polya_site_left=222, correct_polya_site_right=223, terminal_site_match_left=250,
terminal_site_match_right=251, terminal_site_match_left_precise=252,
terminal_site_match_right_precise=253, aligned_polya_tail=300,
major_exon_elongation_right=1003, major_exon_elongation_left=1005, extra_intron_known=1011,
extra_intron_novel=1012, extra_intron_flanking_left=1013, extra_intron_flanking_right=1015.

### 3.4 `SupplementaryMatchConstants` (`isoform_assignment.py:405-412`)
extra_left_mod_position = 1073741823 (2^30-1), extra_right_mod_position = 1073741825 (2^30+1),
undefined_position = 2147483648 (2^31; `undefined_region` = [2^31, 2^31], the MatchEvent default),
absent_position = 2147483647 (2^31-1). Event regions and `info` are i64 in Rust.

### 3.5 Others
- `TranscriptModelType`: known=1, novel_in_catalog=2, novel_not_in_catalog=10.
- `MultimapResolvingStrategy`: ignore_multimapper=0, merge=1, take_best=2.
- `StrandnessReportingLevel`: only_canonical=1, only_stranded=2, all=3, auto=10.
- `AmbiguityResolvingMethod` (`long_read_assigner.py:48`): none=-10, monoexon_only=10, monoexon_and_fsm=20, all=30.
- CIGAR ops (BAM): M=0, I=1, D=2, N=3, S=4, H=5, P=6, `=`=7, X=8.

---------------------------------------------------------------------------------------------------
## 4. Per-unit payloads

Common sub-objects: `Iv` = `[int,int]`; `Prof` (`MappedReadProfile`) =
`{"gene_profile":[int], "read_profile":[int], "read_features":[Iv], "gene_profile_range":[int,int]}`;
`CProf` (`CombinedReadProfiles`) = `{"intron": Prof, "exon": Prof|null, "split_exon": Prof, "polya": [4 ints]}`;
`Aln` (alignment as seen by IsoQuant) =
`{"qname": str, "flag": int, "pos0": int, "ref_end0": int, "mapq": int, "is_reverse": bool, "cigar": [[op,len]], "seq": str|null}`
(`seq` = `query_sequence`, soft clips included, hard clips not; `null` if absent).

### U1 — annotation and IsoformIndex construction

**(a) `fixtures/annotation/U1.annotation.jsonl.gz`** (`iqcap.dump_annotation`, no IsoQuant run).
One case per chromosome, `case = "annotation:<chr>"`.
- `in`: `{"fn":"dump_annotation", "db": path, "gff3": path, "chr": "8", "gene_types": "isoquant"}`
- `out`: `{"chr": "8", "genes": [Gene]}`; genes in `db.features_of_type('gene', order_by='start')`
  order restricted to `chr` (gffutils order for ties = the DB's order, which `get_gene_info_for_region`
  also uses after its stable `sorted(key=start)`).

`Gene`: `{"id", "start", "end", "strand", "source", "attributes": [[k,[v...]]], "attrs": str, "transcripts": [Tx]}`
- `attributes`: raw gffutils attributes, in gffutils order; `attrs`: the string
  `set_gene_attributes` builds for the gene.
`Tx`: `{"id", "type", "start", "end", "strand", "source", "attributes", "attrs", "exons": [Iv], "other_features": [[s,e,type]], "children_order_index": int}`
- transcripts = `db.children(gene, featuretype=DEFAULT_RNA_TYPES, order_by='start')` (the order of
  `set_introns_and_exons`); `children_order_index` = position of this transcript in the *unordered*
  `db.children(gene, featuretype=DEFAULT_RNA_TYPES)` iteration (used by `isoform_strands`,
  `gene_id_map`, `sources`, `other_features` dict order).
- `exons`: after the CDS-as-exon rule (`gene_info.py:600-649`); a transcript with no exon/CDS is
  emitted with `"exons": []` (4.1 skips it in `all_isoforms_exons`, but keeps it in strands/maps).
- `other_features`: `children(t)` features of type CDS/start_codon/stop_codon/UTR in gffutils order.
- `attrs` (transcript): the string `set_gene_attributes` builds for that transcript.

**(b) `fn = "get_gene_info_for_region"`** (`alignment_processor.py:591-605`), `case = <chr>:<block>:n0`
- `in`: `{"fn", "chr", "region": [int,int] (as passed, see §0.5), "delta": int, "needs_reference": bool}`
- `out`: `{"gene_ids": [string] (sorted by start, as passed to GeneInfo), "index": sha1, "empty": bool}`
  (gene query is `genedb.region(seqid=chr, start=region[0], end=region[1], featuretype="gene")`).

**(c) `fn = "from_models"`** (`gene_info.py:283-356`), `case = ...:<block>:from_models:n<k>`
- `in`: `{"fn", "models": [Model], "delta": int, "call_site": string}`
- `out`: `{"index": sha1}`
- Written only in the construct phase, inside selected blocks.
- `call_site`: name of the immediate caller function (`assign_reads_to_models`, `detect_similar_isoforms`,
  `_single_model_assigner`, `build_model_read_assignments`, ...).

### U2 — BAM records (`iqcap.dump_bam`, no IsoQuant run)

`case = "<sample>:<region>:<n>"`, n = 0-based record index in `pysam.fetch(chr, start0, end)` order.
- `in`: `{"fn":"fetch", "bam": path, "region": "8:1001-2000" (1-based closed, as samtools), "n": int}`
- `out`: `{"qname", "flag", "ref_id", "ref_name", "pos0", "ref_end0", "mapq", "cigar": [[op,len]], "seq": str|null, "qual_len": int, "tags": {"HP": {"type": "Z"|"i"|..., "value": str|int}} }`
  (`tags` holds only tags listed by `--tags`, default `HP`; absent tag → key absent).

### U3 — coverage regions and blocks

- **`fn = "split_coverage_regions"`** (`alignment_processor.py:640-682`)
  - `in`: `{"fn", "chr", "region": [int,int] (0-based closed storage region), "read_count": int, "coverage": [[bin, count]] sorted by bin}`
    (coverage copied before the call; 4.1 inserts zero bins while reading).
  - `out`: `{"splits": [[[int,int], int]]}` = list of `((start, end), max_coverage)`.
- **`fn = "process_alignments_in_region"`** (`:337-375`), one per block
  - `in`: `{"fn", "chr", "region": [int,int], "max_coverage": int, "chr_len": int, "small_chr": bool, "coverage_cutoff": int}`
  - `out`: `{"block_no": int, "skip_read_fraction": int, "gene_ids": [string], "index": sha1, "genic": bool}`
- **`fn = "cluster"`** (a `forward_alignments` call = one gap-separated cluster, `:324-335`)
  - `in`: `{"fn", "chr", "region": [int,int] (0-based closed), "read_count": int}`
  - `out`: `{"n_blocks": int, "first_block_no": int}`
- **`fn = "genes_in_region"`** (`get_genes_in_region`, `:583-589`; the gap test of `process`)
  - `in`: `{"fn", "chr", "region": [int,int]}`; `out`: `{"gene_ids": [string] sorted}`.

### U4 — CIGAR→blocks, polyA/T, deletions

- **`fn = "process_alignment"`**: `AlignmentInfo.__init__` + `AlignmentInfo.add_polya_info`
  (written when `add_polya_info` returns; for alignments without exons it is written from `__init__`
  with `polya = null`).
  - `in`: `{"fn", "read": {"qname","aln_no"}, "aln": Aln, "polya_window": int, "polya_fraction": f64, "max_fake_terminal_exon_len": int, "external_strand": str|null, "external_no_tail": bool}`
  - `out`: `{"init": {"read_exons": [Iv], "read_blocks": [[int,int]], "cigar_blocks": [[int,int]]}, "polya": [4 ints]|null, "read_exons": [Iv], "read_blocks": [[int,int]], "cigar_blocks": [[int,int]], "exons_changed": bool, "polya_called": bool}`
    (`polya` = final `polya_info` after detection, `PolyAFixer` correction and shifting;
    `read_*` = after trimming).
- **`fn = "extract_relevant_deletions"`** (`transcript_splice_site_corrector.py:75-100`)
  - `in`: `{"fn", "read": {...}, "cigar": [[op,len]], "ref_start": int, "exons": [Iv], "margin": int, "pos0": int|null}`
    (`ref_start` = the 4.1 anchor actually passed, `read_exons[0][0]`; `pos0` = the alignment's
    `reference_start` (0-based) — the D-P3 anchor. `pos0` is present only in fixtures captured after
    batch 2 (`b2_*` dirs); older fixtures lack it.)
  - `out`: `{"deletions": [[int,int]]}`

### U5 — read profiles

- **`fn = "construct_profiles"`** (`long_read_profiles.py:268-281`)
  - `in`: `{"fn", "call_site": string, "read": {"qname","aln_no"}|{"read_idx": int}|null (provisional models), "index": sha1, "params": AssignParams, "blocks": [Iv], "polya": [4 ints]}`
  - `out`: `CProf` without `polya` → `{"intron": Prof, "exon": Prof|null, "split_exon": Prof}`
- **`fn = "match_genomic_features"`** (`:73-109`), called by the exon corrector
  - `in`: `{"fn", "call_site", "read": {...}, "known_features": [Iv], "delta": int, "read_features": [Iv]}`
  - `out`: `{"features": [Iv]}`
  - `delta`: in fixtures captured after batch 2 (`b2_*` dirs) this is the comparator's **effective** delta,
    `self.comparator.keywords["delta"]` (the `partial(equal_ranges, delta=...)` the matcher uses). Older
    fixtures recorded the constructor kwarg `self.delta`, which is `0` at the `exon_corrector` call site
    (harness bug); for those, use `params.json` `delta` (the runner workaround).

`call_site` (also U6): `collect` (`process_genic`), `fl` (`construct_fl_isoforms`), `similar`
(`detect_similar_isoforms`), `reassign` (`assign_reads_to_models`), `final_single`
(`_model_read_assignment_single`), `final_full` (`build_model_read_assignments`), `sqanti`
(`compare_models_with_known`), `exon_corrector`, or the raw caller function name otherwise.

### U6 — read→isoform assignment

- **`fn = "assign_to_isoform"`** (`long_read_assigner.py:438-507`)
  - `in`: `{"fn", "call_site": string, "read": {...}|{"read_idx": int}|{"id": string}, "read_id": string (the id passed: QNAME, or a provisional transcript id for `fl`/`similar`), "index": sha1, "params": AssignParams, "quick_mode": bool, "profile": CProf}`
  - `out`: `{"assignment": Assignment|null}` (null when 4.1 returns `None`)
  - `read` identity: collect = `{"qname","aln_no"}`; reassign/final = `{"read_idx"}`; fl/similar =
    `{"id": <transcript id>}`.

### U7 — exon correction, strand, antisense

- **`fn = "correct_assigned_read"`** (`exon_corrector.py:29-41`)
  - `in`: `{"fn", "read": {...}, "index": sha1, "params": AssignParams, "aln": Aln, "read_exons": [Iv], "read_start": int, "read_end": int, "intron_profile": Prof, "assignment": Assignment, "ref": [RefWin]}`
    (`ref`: one window per read exon (`read_exons`, i.e. after polyA/T trimming), covering 1-based `[s-32, e+32]`
    (0-based `[s-33, e+32)`), merged when overlapping and clipped at the chromosome ends. `get_error_count` only reads
    positions between a read splice site and the matched reference splice site (within `delta`), so these windows
    cover every base 4.1 reads; `[]` when the run has no reference).
  - `out`: `{"corrected_exons": [Iv]}`
- **`fn = "get_assignment_strand"`** (`alignment_processor.py:538-551`)
  - `in`: `{"fn", "read": {...}, "assignment": Assignment, "polya": [4], "exons": [Iv], "corrected_introns": [Iv], "intron_sites": [[Iv, string, string]]}`
    (`intron_sites`: per corrected intron `(intron, donor, acceptor)` = `chr[s-1:s+1]`, `chr[e-2:e]`
    verbatim, i.e. what `get_intron_strand` reads).
  - `out`: `{"strand": strand}`
- **`fn = "check_antisense"`** (`:684-695`)
  - `in`: `{"fn", "read": {...}, "strand": strand, "assignment": Assignment}`
  - `out`: `{"assignment": Assignment}` (after in-place mutation)

### U8 — collect driver, per block

- **`fn = "collect_block"`**: one case per block (`TmpFileAssignmentPrinter.add_gene_info` + the
  following `add_read_info` calls), `case = ...:<block>:n0`.
  - `in`: `{"fn", "chr", "block_no", "region": [int,int] (process_alignments_in_region region), "index": sha1, "genic": bool, "skip_read_fraction": int, "max_coverage": int}`
  - `out`: `{"reads": [ReadRecord]}` in dump (collect) order, `idx` = 0.. in that order, *before*
    multimapper resolution (compare with `block/` for the post-resolution form).
  The BAM records of the block come from the slice (`U2` dump / the slice BAM): all records
  overlapping `region` (0-based fetch `[region[0], region[1]+1)`).

### U9 — multimapper resolution

- **`fn = "resolve"`** (`multimap_resolver.py:32-48`)
  - `in`: `{"fn", "strategy": int, "assignments": [Basic]}`
  - `out`: `{"assignments": [Basic]}` (the returned list, after mutation)
  - `Basic` (`BasicReadAssignment`): `{"assignment_id", "read_id", "chr", "start", "end", "region": [int,int], "multimapper": bool, "polya_found": bool, "atype": int, "gene_atype": int, "penalty": f64, "genes": [string], "isoforms": [string]}`
    (`genes`/`isoforms` in 4.1 list order = iteration order of a set of pool ints; treat as sets).
  - `case = "<sample>:<slice>:<chr>:-:<read_id>"`.

### U10a — terminal peaks and XGBoost

- **`fn = "detect_peaks_batch"`** (`terminal_peaks.py:84-140`)
  - `in`: `{"fn", "model": "polya"|"tss", "histograms": [ [[pos,count]] | null ]}` (each
    histogram in dict iteration order — load-bearing: it fixes the order of the float sums)
  - `out`: `{"records": [Rec|null], "rows": [Row], "peaks": [[[position,count,left,right]]]}`
    - `Rec` (per histogram): `{"start": int, "mode": int, "var": f64, "skew": f64*, "entropy": f64, "mean_height": f64, "peak_count": int, "histogram": [int] (padded), "raw_peaks": {"peak_location":[int],"peak_width":[f64],"peak_left":[int],"peak_right":[int],"peak_heights":[int]}|null}`
    - `Row` (classifier input row, in batch order): `{"hist": int, "features": [8 x f64*] (FEATURE_COLUMNS order, float64 before the float32 cast), "proba": f64, "margin": f64, "label": int}`
- **`fn = "xgb_predict"`** (every `XGBClassifier.predict` call, incl. `TerminalCounter._predict_rows`)
  - `in`: `{"fn", "model": "polya"|"tss", "caller": string, "features": [[8 x f64*]]}`
  - `out`: `{"labels": [int], "proba": [f64], "margin": [f64]}`
  `proba` = `predict_proba[:,1]`, `margin` = `predict(output_margin=True)`; both side calls, computed on
  the same frame.

### U10b — per-block polyA/TSS predictions

- **`fn = "flush"`** (`terminal_counter.py:230-256`), ungrouped counters only
  - `in`: `{"fn", "kind": "polya"|"tss", "block": BlockRef, "transcripts": [[tx, {"chr": str, "gene_id": str, "data": [int], "annotated": int}]]}`
    (`self.transcripts` accumulator before the call, dict order; `data` = read positions in add order)
  - `out`: `{"predictions": [[tx, int]]}` (`last_gene_predictions`)

### U11 — intron graph

- **`fn = "IntronGraph.__init__"`**, `case = ...:<block>:n0`
  - `in`: `{"fn", "block": BlockRef, "index": sha1, "polya_predictions": [[tx,int]]|null, "tss_predictions": [[tx,int]]|null}`
  - `out`: `{"graph": Graph}`

### U12 — paths

- **`fn = "IntronPathStorage.fill"`**
  - `in`: `{"fn", "block": BlockRef, "graph": "<U11 case key>"}`
  - `out`: `{"paths": [[path, int]] (dict order), "fl_paths": [path] (sorted), "paths_to_reads": [[path, [read_idx]]] (dict order)}`
- **`fn = "get_known_spliced_isoforms"`** (`context.py:202-210`)
  - `in`: `{"fn", "block": BlockRef, "graph": "<U11 case key>", "index": sha1}`
  - `out`: `{"known": [[path, tx]]}` (dict order)

### U13-U19 — construction stages

Except U18, every `in` begins with `{"fn", "block": BlockRef, "index": sha1 (the block GeneInfo), "graph": "<U11 case key>"|null, "paths": "<U12 fill case key>"|null}`
(these are always present, whether or not the stage reads them), followed by the extra fields in the table (the `Store` before the call, when listed); `out`
carries the `Store` after the call and the return value if any. `params.json` holds the args.

| unit | `fn` | extra `in` | `out` |
|---|---|---|---|
| U13 | `construct_fl_isoforms` | `store`, `intron_sites` (b2+) | `store` |
| U14 | `construct_assignment_based_isoforms` | `store`, `reads: [read_idx]` (the argument; normally all) | `store` |
| U15a | `pre_filter_transcripts` | `store` | `store` |
| U15a | `filter_transcripts` | `store` | `store` |
| U15a | `detect_similar_isoforms` | `models: [Model]`, `remove_unconfirmed: bool` | `to_substitute: [[string,string]]` (dict order) |
| U15b | `correct_novel_transcript_ends` | `model: Model`, `reads: [read_idx]` (no store) | `model: Model` (mutated) |
| U15b | `_mark_terminal_confirmation` | `model: Model`, `reads: [read_idx]` (= `store.transcript_read_ids[model.id]`; no store) | `model: Model` |
| U15b | `_add_known_alternative_end_models` | `store` | `store` |
| U15b | `_drop_duplicate_alt_end_models` | `store` | `store` |
| U15b | `derive_alternative_end_models` | `model: Model`, `reads: [read_idx]`, `id_counter: int` (before; no store) | `result: [[Model, [read_idx]]]`, `id_counter: int` |
| U16 | `correct_transcript_splice_sites` | `store`, `ref: [RefWin]` (b2+) | `store` |
| U16 | `_correct_one_model_splice_sites` | `model: Model`, `reads: [{"idx", "start": exons[0][0], "end": exons[-1][1], "boundary_deletions": [[int,int]]}]`, `ref: [RefWin]` (one window `[b-65, b+64)` per exon boundary `b`) | `exons: [Iv]|null` |
| U17 | `assign_reads_to_models` | `store`, `reads: [read_idx]` | `store` |
| U17 | `resolve_read_ownership_conflicts` | `store` | `store` |
| U18 | `join_transcripts` | `in` = `{fn, block, index, models: [Model]}` (joiner input) | `models: [Model]` (returned list) |
| U19 | `build_model_read_assignments` | `store`, `reads: [read_idx]` | `assignments: [{"read_idx": int, "assignment": Assignment}]` (one per input read, in order) |

Batch-2 additions (present only in `b2_*` captures):
- U13 `intron_sites: [[vertex, string, string]]`: for every vertex `[a,b]` of `path[1:-1]` of every FL path
  (`path_storage.fl_paths`), sorted and unique, `(vertex, chr[a-1:a+1], chr[b-2:b])` verbatim, i.e. the two
  slices `get_intron_strand` reads (it upper-cases them). Same layout as U7 `get_assignment_strand` `intron_sites`.
  Terminal vertices that sit inside `path[1:-1]` (code < 0) are included too, sliced with Python semantics
  (negative start → slice from the chromosome end or empty), exactly as 4.1 would. `""` sites when no reference.
- U16 store-level `ref: [RefWin]`: windows `[b-65, b+64)` (0-based half-open, as in `_correct_one_model_splice_sites`)
  around every exon boundary `b` of every non-known model with >= 2 exons in the `store` before the call,
  sorted, merged when overlapping, clipped at the chromosome ends. `[]` when no reference.

Nested calls (U1c `from_models`, U5, U6, U10a, U15b inside U15a, ...) are recorded as their own
cases with the corresponding `call_site`, so stage runners can replay them.

### U20 — counters and writers

- **`fn = "gff_dump"`** (`GFFPrinter.dump`, `transcript_printer.py:56-138`)
  - `in`: `{"fn", "block": BlockRef, "printer": "transcript_models"|"extended_annotation", "index": sha1 (gene_info passed), "models": [Model], "printed_gene_ids": [string] (before, sorted), "exon_id_counter": int (before)}`
  - `out`: `{"lines": [string]}` (exact text written, each with trailing `\n`), `"exon_id_counter": int (after)`
- **`fn = "counter_add_read_info"`** (`AssignedFeatureCounter.add_read_info`) for the transcript-model
  and gene-model counters (`discovered_*`), ungrouped and grouped
  - `in`: `{"fn", "counter": {"name": string (output file basename), "grouped": bool, "group_index": int}, "read": {"read_idx": int}|null, "assignment": Assignment, "read_group": [string]}`
  - `out`: `{"features": [[feature_id, [[group, f64 before, f64 after]]]], "ambiguous": [int,int], "not_assigned": [int,int], "confirmed_added": [string]}`
- **`fn = "counter_dump"`** (`AssignedFeatureCounter.dump`)
  - `in`: `{"fn", "counter": {...}, "feature_counter": [[feature, [[group_name, f64]]]], "all_features": [string], "confirmed_features": [string], "ambiguous": int, "not_assigned": int, "not_aligned": int}`
  - `out`: `{"files": [[basename, text]]}` (contents after `dump()` of the counts, stats and usable files that exist)
  - Written once per non-empty model counter per chromosome. Counters that saw no reads are skipped.
- Final outputs: `outputs/` holds the run's `OUT/*.transcript_models.gtf`,
  `*.discovered_transcript_counts.tsv`, `*.discovered_transcript_grouped_counts.linear.tsv` and
  related grouped files verbatim.

---------------------------------------------------------------------------------------------------
## 5. Capture notes and known quirks of the recorded data

- **Phases.** The collect phase (U3, U4, U5/U6/U7 with `call_site` `collect`/`exon_corrector`, U8) runs for all
  chromosomes first. Multimapper resolution (U9) follows, then the construct phase (U10b-U20, construct-side
  U1c/U5/U6/U10a). With `--threads 1` the run is single-process and deterministic.
- **Valley-straddling reads** (D-R2): 4.1 processes a read overlapping two coverage-split blocks in both
  blocks. These duplicates are what produce U9 cases in slices without secondary alignments.
- **Empty loader blocks**: a block can have `n_reads = 0` in `block/` (all its reads were dropped at collect,
  or suspended by multimapper resolution). Its construct-side units still fire and are recorded.
- **`xgb_predict` from `TerminalCounter.dump`** (the chromosome-level polyA/TSS prediction run after the last
  block) is keyed with the last block of the chromosome. Its `caller` is `_predict_rows`.
- **Per-block counters** (U20 `counter_add_read_info`) record only the features touched by the read. The
  chromosome-level state is in `counter_dump`.
- **HP tags** are synthesized by `iqcap.slice`, and every capture runs with `--read_group tag:HP`. The grouped
  counters and `read_group` fields are therefore populated (`A`, `B`, `unphased`).
- **Subsampled slices**: `slice_def.subsample < 1` means the slice BAM keeps reads with
  `int(sha1('sub:'+qname)[:16],16)/2**64 < subsample` (all alignments of a read together).

---------------------------------------------------------------------------------------------------
## 6. Changelog
- v1 (2026-10-29): first draft.
- v1 (2026-10-29, additive, requested by iq_types): pinned the dumps settings (ensure_ascii=True);
  AssignParams enums written as int codes; added Store `read_assignment_counts_idx` and
  `construction_assignment_idx`.
- v1 (2026-10-29, additive, doc catch-up): U7 `correct_assigned_read` `in.ref`
  ([RefWin], one merged window per read exon ±32) and U7 `get_assignment_strand` `in.intron_sites`;
  U5 `construct_profiles` `in.read` may be `null` (provisional models / no read in scope). These fields were
  present in every M0a fixture; only the changelog entry was missing.
- v1 (2026-10-29, batch 2, additive; all M0a fixture dirs unchanged, new fields appear only in `b2_*` dirs):
  - U4 `extract_relevant_deletions` `in.pos0` (alignment `reference_start`, 0-based).
  - U5 `match_genomic_features` `in.delta` **semantics fix** (#6): now the comparator's effective delta
    (`comparator.keywords["delta"]`); field name/type unchanged. M0a fixtures keep the old (constructor
    kwarg, 0 at `exon_corrector`) value; runners must keep the `params.json` workaround for non-`b2_` dirs.
  - U13 `construct_fl_isoforms` `in.intron_sites`, U16 `correct_transcript_splice_sites` `in.ref` (#14).
  - `fixture_manifest.json` `slice_def.keep_names` (read-name-restricted slices).
  - New capture dirs `b2_*`, some run with `--report_novel_unspliced true`, `--novel_apa`, or other
    chromosomes (see harness/README.md "Batch 2"); check `params.json` per dir.
