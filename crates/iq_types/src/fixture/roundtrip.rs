//! Schema conformance: re-encode fixture lines through the typed DTOs and check that the
//! bytes are unchanged (U0 pass criterion "Python dump → Rust load → Rust dump →
//! byte-equal JSON"). Key order is part of the contract, so no normalisation is applied.

use super::schema::*;
use super::{pyjson, FixError, FixResult};
use serde::de::DeserializeOwned;
use serde::Serialize;

fn re<I, O>(line: &str) -> FixResult<String>
where
    I: DeserializeOwned + Serialize,
    O: DeserializeOwned + Serialize,
{
    // exceptions (§0.8) have their own out shape
    let raw: Case = serde_json::from_str(line)?;
    if raw.out.get("exception").is_some() {
        let c: Case<I, ExceptionOut> = serde_json::from_str(line)?;
        return Ok(pyjson::to_py_string(&c)?);
    }
    let c: Case<I, O> = serde_json::from_str(line)?;
    Ok(pyjson::to_py_string(&c)?)
}

macro_rules! dispatch {
    ($unit:expr, $fn_:expr, $line:expr; $( ($u:literal, $f:literal, $i:ty, $o:ty) ),+ $(,)?) => {
        match ($unit, $fn_) {
            $( ($u, $f) => re::<$i, $o>($line), )+
            (u, f) => Err(FixError(format!("no DTO for unit {u} fn {f}"))),
        }
    };
}

/// Re-encodes one unit-file line through the DTOs of its `(unit, in.fn)`. The result must
/// equal `line` (without its trailing newline) for a conforming fixture.
pub fn reencode_case_line(line: &str) -> FixResult<String> {
    #[derive(serde::Deserialize)]
    struct Head {
        unit: String,
        #[serde(rename = "in")]
        input: InHead,
    }
    #[derive(serde::Deserialize)]
    struct InHead {
        #[serde(rename = "fn")]
        fn_: String,
    }
    let h: Head = serde_json::from_str(line)?;
    dispatch!(h.unit.as_str(), h.input.fn_.as_str(), line;
        ("U1", "dump_annotation", U1DumpAnnotationIn, U1DumpAnnotationOut),
        ("U1", "get_gene_info_for_region", U1GeneInfoForRegionIn, U1GeneInfoForRegionOut),
        ("U1", "from_models", U1FromModelsIn, IndexOut),
        ("U2", "fetch", U2FetchIn, U2FetchOut),
        ("U3", "split_coverage_regions", U3SplitCoverageIn, U3SplitCoverageOut),
        ("U3", "process_alignments_in_region", U3ProcessRegionIn, U3ProcessRegionOut),
        ("U3", "cluster", U3ClusterIn, U3ClusterOut),
        ("U3", "genes_in_region", U3GenesInRegionIn, GeneIdsOut),
        ("U4", "process_alignment", U4ProcessAlignmentIn, U4ProcessAlignmentOut),
        ("U4", "extract_relevant_deletions", U4ExtractDeletionsIn, U4ExtractDeletionsOut),
        ("U5", "construct_profiles", U5ConstructProfilesIn, U5ConstructProfilesOut),
        ("U5", "match_genomic_features", U5MatchGenomicFeaturesIn, U5MatchGenomicFeaturesOut),
        ("U6", "assign_to_isoform", U6AssignIn, U6AssignOut),
        ("U7", "correct_assigned_read", U7CorrectIn, U7CorrectOut),
        ("U7", "get_assignment_strand", U7StrandIn, U7StrandOut),
        ("U7", "check_antisense", U7AntisenseIn, AssignmentOut),
        ("U8", "collect_block", U8CollectBlockIn, U8CollectBlockOut),
        ("U9", "resolve", U9ResolveIn, U9ResolveOut),
        ("U10a", "detect_peaks_batch", U10aDetectPeaksIn, U10aDetectPeaksOut),
        ("U10a", "xgb_predict", U10aXgbPredictIn, U10aXgbPredictOut),
        ("U10b", "flush", U10bFlushIn, PredictionsOut),
        ("U11", "IntronGraph.__init__", U11GraphIn, U11GraphOut),
        ("U12", "IntronPathStorage.fill", U12FillIn, U12FillOut),
        ("U12", "get_known_spliced_isoforms", U12KnownSplicedIn, U12KnownSplicedOut),
        ("U13", "construct_fl_isoforms", StageStoreIn, StoreOut),
        ("U14", "construct_assignment_based_isoforms", StageStoreIn, StoreOut),
        ("U15a", "pre_filter_transcripts", StageStoreIn, StoreOut),
        ("U15a", "filter_transcripts", StageStoreIn, StoreOut),
        ("U15a", "detect_similar_isoforms", U15aDetectSimilarIn, U15aDetectSimilarOut),
        ("U15b", "correct_novel_transcript_ends", StageModelIn, ModelOut),
        ("U15b", "_mark_terminal_confirmation", StageModelIn, ModelOut),
        ("U15b", "_add_known_alternative_end_models", StageStoreIn, StoreOut),
        ("U15b", "_drop_duplicate_alt_end_models", StageStoreIn, StoreOut),
        ("U15b", "derive_alternative_end_models", StageModelIn, U15bDeriveAltEndOut),
        ("U16", "correct_transcript_splice_sites", StageStoreIn, StoreOut),
        ("U16", "_correct_one_model_splice_sites", U16CorrectOneIn, U16CorrectOneOut),
        ("U17", "assign_reads_to_models", StageStoreIn, StoreOut),
        ("U17", "resolve_read_ownership_conflicts", StageStoreIn, StoreOut),
        ("U18", "join_transcripts", U18JoinIn, ModelsOut),
        ("U19", "build_model_read_assignments", StageStoreIn, U19AssignmentsOut),
        ("U20", "gff_dump", U20GffDumpIn, U20GffDumpOut),
        ("U20", "counter_add_read_info", U20CounterAddIn, U20CounterAddOut),
        ("U20", "counter_dump", U20CounterDumpIn, U20CounterDumpOut),
    )
}

/// Re-encodes one line of a block file (line 0 = header, others = `ReadRecord`).
pub fn reencode_block_line(line: &str, is_header: bool) -> FixResult<String> {
    if is_header {
        Ok(pyjson::to_py_string(&serde_json::from_str::<BlockHeader>(line)?)?)
    } else {
        Ok(pyjson::to_py_string(&serde_json::from_str::<ReadRecordFixture>(line)?)?)
    }
}

/// Re-encodes an `index/<sha1>.json` file content (no trailing newline).
pub fn reencode_index(text: &str) -> FixResult<String> {
    Ok(pyjson::to_py_string(&serde_json::from_str::<IndexDump>(text)?)?)
}

/// `params.json`: sorted keys + trailing newline (schema §0.2).
pub fn reencode_params(text: &str) -> FixResult<String> {
    let v: serde_json::Value = serde_json::from_str(text)?;
    Ok(pyjson::to_py_string_sorted(&v)? + "\n")
}

/// Result of checking one file.
#[derive(Debug, Default, Clone)]
pub struct FileCheck {
    pub lines: u64,
    pub ok: u64,
    /// `(line_no (1-based), error or first differing byte offset)`.
    pub failures: Vec<(u64, String)>,
}

/// First byte offset where `a` and `b` differ, with context.
pub fn first_difference(a: &str, b: &str) -> String {
    let pos = a.bytes().zip(b.bytes()).position(|(x, y)| x != y).unwrap_or(a.len().min(b.len()));
    let ctx = |s: &str| {
        let lo = pos.saturating_sub(40);
        let hi = (pos + 40).min(s.len());
        s.get(lo..hi).unwrap_or("").to_owned()
    };
    format!("byte {pos}: python …{}… vs rust …{}…", ctx(a), ctx(b))
}

/// Checks a unit file, block file, index file or params.json by re-encoding.
pub fn check_file(path: &std::path::Path, max_failures: usize) -> FixResult<FileCheck> {
    check_file_limited(path, max_failures, u64::MAX)
}

/// [`check_file`] reading at most `max_lines` lines of a JSONL file.
pub fn check_file_limited(path: &std::path::Path, max_failures: usize, max_lines: u64) -> FixResult<FileCheck> {
    use std::io::BufRead;
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    let mut fc = FileCheck::default();
    let record = |fc: &mut FileCheck, no: u64, orig: &str, got: FixResult<String>| {
        fc.lines += 1;
        match got {
            Ok(s) if s == orig => fc.ok += 1,
            Ok(s) => {
                if fc.failures.len() < max_failures {
                    fc.failures.push((no, first_difference(orig, &s)));
                }
            }
            Err(e) => {
                if fc.failures.len() < max_failures {
                    fc.failures.push((no, e.0));
                }
            }
        }
    };
    if name.ends_with(".jsonl.gz") || name.ends_with(".jsonl") {
        let is_block = path.components().any(|c| c.as_os_str() == "block");
        for (k, l) in super::open_lines(path)?.lines().enumerate().take(max_lines.try_into().unwrap_or(usize::MAX)) {
            let l = l?;
            if l.is_empty() {
                continue;
            }
            let got = if is_block { reencode_block_line(&l, k == 0) } else { reencode_case_line(&l) };
            record(&mut fc, k as u64 + 1, &l, got);
        }
    } else if name == "params.json" {
        let t = std::fs::read_to_string(path)?;
        record(&mut fc, 1, &t, reencode_params(&t));
    } else if name.ends_with(".json") {
        let t = std::fs::read_to_string(path)?;
        record(&mut fc, 1, &t, reencode_index(&t));
    } else {
        return Err(FixError(format!("{}: not a fixture file", path.display())));
    }
    Ok(fc)
}

/// Every fixture file under `root` that [`check_file`] understands (unit files, block
/// files, index dumps, params.json), sorted. Skips `outputs/`, `slice/` and the manifest.
pub fn fixture_files(root: &std::path::Path) -> Vec<std::path::PathBuf> {
    fn walk(d: &std::path::Path, out: &mut Vec<std::path::PathBuf>) {
        let Ok(rd) = std::fs::read_dir(d) else { return };
        for e in rd.flatten() {
            let p = e.path();
            let n = p.file_name().and_then(|n| n.to_str()).unwrap_or("").to_owned();
            if p.is_dir() {
                if n != "outputs" && n != "slice" {
                    walk(&p, out);
                }
            } else if n.ends_with(".jsonl.gz") || n == "params.json" || (n.ends_with(".json") && p.parent().is_some_and(|q| q.ends_with("index"))) {
                out.push(p);
            }
        }
    }
    let mut v = Vec::new();
    walk(root, &mut v);
    v.sort();
    v
}
