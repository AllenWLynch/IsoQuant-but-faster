//! Effective parameters (the 4.1 `args` subset iqrs uses) and the two built-in parameter
//! sets of design §7.0. There is no option surface beyond these presets; tests may load a
//! 4.1 `params.json` instead ([`Params::from_params_json`]).
//!
//! Source citations are `isoquant.py:N` (repository root) unless a path is given. Values
//! below are those of IsoQuant 4.1.0, commit `cea5337`.

use crate::enums::{AmbiguityResolving, MultimapStrategy, PolyAUsage, StrandnessReporting};
use serde::{Deserialize, Serialize};

/// `AmbiguityResolvingMethod.minimal_score` (`assignment/long_read_assigner.py:54`).
pub const MINIMAL_SCORE: f64 = -0.5;
/// `AmbiguityResolvingMethod.top_scored_factor` (`assignment/long_read_assigner.py:55`).
pub const TOP_SCORED_FACTOR: f64 = 1.5;

/// Data type preset (`--data-type`).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DataType {
    /// `--data_type pacbio_ccs --model_construction_strategy sensitive_pacbio --fl_data`.
    Pacbio,
    /// `--data_type nanopore --model_construction_strategy sensitive_ont`.
    Ont,
}

/// Read-to-isoform matching strategy (`MatchingStrategy` namedtuple, `isoquant.py:1091-1101`)
/// plus the fixed matching constants of `set_matching_options` (`:1110-1127`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MatchingParams {
    pub delta: i32,
    pub max_intron_shift: i32,
    pub max_missed_exon_len: i32,
    pub max_fake_terminal_exon_len: i32,
    pub max_suspicious_intron_abs_len: i32,
    pub max_suspicious_intron_rel_len: f64,
    pub resolve_ambiguous: AmbiguityResolving,
    pub correct_minor_errors: bool,
    /// `:1110`.
    pub minor_exon_extension: i32,
    /// `:1111`.
    pub major_exon_extension: i32,
    /// `:1118`.
    pub min_abs_exon_overlap: i32,
    /// `:1119`.
    pub min_rel_exon_overlap: f64,
    /// `:1120`.
    pub micro_intron_length: i32,
    /// `min(30, max_intron_shift)`, `:1121`.
    pub max_intron_abs_diff: i32,
    /// `:1122`.
    pub max_intron_rel_diff: f64,
    /// `= minor_exon_extension`, `:1123`.
    pub apa_delta: i32,
    /// `:1124`.
    pub minimal_exon_overlap: i32,
    /// `:1125`.
    pub minimal_intron_absence_overlap: i32,
    /// `:1126`.
    pub polya_window: i32,
    /// `:1127`.
    pub polya_fraction: f64,
    /// `--count_exons` (off; exon profiles are not built).
    pub count_exons: bool,
}

/// `SplicSiteCorrectionStrategy` (`isoquant.py:1143-1161`).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SpliceCorrection {
    pub fuzzy_junctions: bool,
    pub intron_shifts: bool,
    pub skipped_exons: bool,
    pub terminal_exons: bool,
    pub fake_terminal_exons: bool,
    pub microintron_retention: bool,
}

/// `ModelConstructionStrategy` (`isoquant.py:1189-1224`), all 19 fields in tuple order.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ModelStrategy {
    pub min_novel_intron_count: u32,
    pub graph_clustering_ratio: f64,
    pub graph_clustering_distance: i32,
    pub min_novel_isolated_intron_abs: u32,
    pub min_novel_isolated_intron_rel: f64,
    pub terminal_position_abs: u32,
    pub terminal_position_rel: f64,
    pub terminal_internal_position_rel: f64,
    pub min_known_count: u32,
    pub min_nonfl_count: u32,
    pub min_novel_count: u32,
    pub min_novel_count_rel: f64,
    pub min_mono_count_rel: f64,
    pub singleton_adjacent_cov: u32,
    pub fl_only: bool,
    /// `novel_monoexonic`: default of `--report_novel_unspliced` (`:1249-1250`).
    pub novel_monoexonic: bool,
    pub require_monointronic_polya: bool,
    pub require_monoexonic_polya: bool,
    pub report_canonical: StrandnessReporting,
}

/// Read filters (`isoquant.py:289-308`).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FilterParams {
    /// `--min_mapq`, default `None` (`:293-294`).
    pub min_mapq: Option<u8>,
    /// `--inconsistent_mapq_cutoff`, default 5 (`:295-297`). KEEP (changelog 3.4.0 #110).
    pub inconsistent_mapq_cutoff: u8,
    /// `--simple_alignments_mapq_cutoff`, default 1 (`:298-300`).
    pub simple_alignments_mapq_cutoff: u8,
    /// `--max_coverage_small_chr`, default 1000000 (`:301-304`).
    pub max_coverage_small_chr: i64,
    /// `--max_coverage_normal_chr`, default -1 (`:305-308`).
    pub max_coverage_normal_chr: i64,
}

/// All effective parameters of one run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Params {
    pub data_type: DataType,
    /// `--fl_data` (pacbio preset only).
    pub fl_data: bool,
    pub matching: MatchingParams,
    pub splice: SpliceCorrection,
    pub model: ModelStrategy,
    pub filter: FilterParams,
    /// `--report_novel_unspliced`; default `model.novel_monoexonic` (`:1249-1250`).
    pub report_novel_unspliced: bool,
    /// `report_canonical_strategy`: `--report_canonical auto` → strategy value (`:1265-1267`).
    pub report_canonical: StrandnessReporting,
    /// `--polya_requirement`, default `auto` (`:368-372`).
    pub polya_requirement: PolyAUsage,
    /// Run-time decision (barrier, §5.4): polyA fraction ≥ `polya_percentage_threshold`
    /// (`dataset_processor.py:283-300`). `false` until decided (`:1037`).
    pub requires_polya_for_construction: bool,
    /// `:1306`.
    pub simple_models_mapq_cutoff: u8,
    /// `:1307`.
    pub polya_percentage_threshold: f64,
    /// `:1308`.
    pub low_polya_percentage_threshold: f64,
    /// `:1292`.
    pub indel_near_splice_site_dist: i32,
    /// `:1293`.
    pub upstream_region_len: i32,
    /// `:1247`.
    pub min_mono_exon_coverage: f64,
    /// `multimap_strategy = "take_best"` (`:1295-1299`).
    pub multimap_strategy: MultimapStrategy,
    /// `use_tss_model = fl_data` (`model_construction/context.py:144-147`).
    pub use_tss_model: bool,
    /// `create_alt_ts_nics = bool(genedb)` (`context.py:141`); an annotation is always given.
    pub create_alt_ts_nics: bool,
    /// `--novel_apa`, default off (`isoquant.py:376-377`).
    pub novel_apa: bool,
    /// `not --no_splice_site_correction`, default on (`:381-383`).
    pub splice_site_correction: bool,
}

impl MatchingParams {
    /// `'precise': MatchingStrategy(4, 30, 50, 20, 0, 0.0, 'monoexon_and_fsm', True)`
    /// (`isoquant.py:1098`); pacbio_ccs default (`:1022`).
    pub fn precise() -> Self {
        Self::from_strategy(4, 30, 50, 20, 0, 0.0, AmbiguityResolving::MonoexonAndFsm, true)
    }
    /// `'default': MatchingStrategy(6, 60, 100, 40, 60, 1.0, 'monoexon_and_fsm', True)`
    /// (`isoquant.py:1099`); nanopore default (`:1022`).
    pub fn default_ont() -> Self {
        Self::from_strategy(6, 60, 100, 40, 60, 1.0, AmbiguityResolving::MonoexonAndFsm, true)
    }

    /// `set_matching_options` (`isoquant.py:1090-1134`) for one strategy tuple.
    pub fn from_strategy(
        delta: i32,
        max_intron_shift: i32,
        max_missed_exon_len: i32,
        max_fake_terminal_exon_len: i32,
        max_suspicious_intron_abs_len: i32,
        max_suspicious_intron_rel_len: f64,
        resolve_ambiguous: AmbiguityResolving,
        correct_minor_errors: bool,
    ) -> Self {
        let minor_exon_extension = 50;
        MatchingParams {
            delta,
            max_intron_shift,
            max_missed_exon_len,
            max_fake_terminal_exon_len,
            max_suspicious_intron_abs_len,
            max_suspicious_intron_rel_len,
            resolve_ambiguous,
            correct_minor_errors,
            minor_exon_extension,
            major_exon_extension: 300,
            min_abs_exon_overlap: 10,
            min_rel_exon_overlap: 0.2,
            micro_intron_length: 50,
            max_intron_abs_diff: max_intron_shift.min(30),
            max_intron_rel_diff: 0.2,
            apa_delta: minor_exon_extension,
            minimal_exon_overlap: 5,
            minimal_intron_absence_overlap: 20,
            polya_window: 16,
            polya_fraction: 0.75,
            count_exons: false,
        }
    }
}

impl SpliceCorrection {
    /// `'default_pacbio': (True, False, True, False, False, True)` (`isoquant.py:1149`).
    pub fn default_pacbio() -> Self {
        SpliceCorrection {
            fuzzy_junctions: true,
            intron_shifts: false,
            skipped_exons: true,
            terminal_exons: false,
            fake_terminal_exons: false,
            microintron_retention: true,
        }
    }
    /// `'default_ont': (True, False, True, False, True, True)` (`isoquant.py:1151`).
    pub fn default_ont() -> Self {
        SpliceCorrection { fake_terminal_exons: true, ..Self::default_pacbio() }
    }
}

impl ModelStrategy {
    /// `'sensitive_pacbio': (1, 0.5, 5, 2, 0.005, 1, 0.01, 0.02, 1, 2, 2, 0.005, 0.001, 100,
    /// False, True, False, False, only_stranded)` (`isoquant.py:1206-1207`).
    pub fn sensitive_pacbio() -> Self {
        ModelStrategy {
            min_novel_intron_count: 1,
            graph_clustering_ratio: 0.5,
            graph_clustering_distance: 5,
            min_novel_isolated_intron_abs: 2,
            min_novel_isolated_intron_rel: 0.005,
            terminal_position_abs: 1,
            terminal_position_rel: 0.01,
            terminal_internal_position_rel: 0.02,
            min_known_count: 1,
            min_nonfl_count: 2,
            min_novel_count: 2,
            min_novel_count_rel: 0.005,
            min_mono_count_rel: 0.001,
            singleton_adjacent_cov: 100,
            fl_only: false,
            novel_monoexonic: true,
            require_monointronic_polya: false,
            require_monoexonic_polya: false,
            report_canonical: StrandnessReporting::OnlyStranded,
        }
    }
    /// `'sensitive_ont': (1, 0.5, 20, 3, 0.005, 1, 0.01, 0.02, 1, 2, 3, 0.005, 0.005, 10,
    /// False, True, False, False, only_stranded)` (`isoquant.py:1216-1217`).
    pub fn sensitive_ont() -> Self {
        ModelStrategy {
            graph_clustering_distance: 20,
            min_novel_isolated_intron_abs: 3,
            min_novel_count: 3,
            min_mono_count_rel: 0.005,
            singleton_adjacent_cov: 10,
            ..Self::sensitive_pacbio()
        }
    }
}

impl FilterParams {
    pub fn defaults() -> Self {
        FilterParams {
            min_mapq: None,
            inconsistent_mapq_cutoff: 5,
            simple_alignments_mapq_cutoff: 1,
            max_coverage_small_chr: 1_000_000,
            max_coverage_normal_chr: -1,
        }
    }
}

impl Params {
    /// Built-in parameter set for a data type (design §7.0).
    pub fn preset(dt: DataType) -> Self {
        match dt {
            DataType::Pacbio => Self::pacbio(),
            DataType::Ont => Self::ont(),
        }
    }

    /// `--data_type pacbio_ccs --model_construction_strategy sensitive_pacbio --fl_data`:
    /// matching `precise` (`isoquant.py:1022`), splice correction `default_pacbio`
    /// (`:1032`), `resolve_ambiguous = monoexon_and_fsm` (`:1036`, fl_data).
    pub fn pacbio() -> Self {
        let model = ModelStrategy::sensitive_pacbio();
        Params {
            data_type: DataType::Pacbio,
            fl_data: true,
            matching: MatchingParams::precise(),
            splice: SpliceCorrection::default_pacbio(),
            report_novel_unspliced: model.novel_monoexonic,
            report_canonical: model.report_canonical,
            model,
            filter: FilterParams::defaults(),
            polya_requirement: PolyAUsage::Auto,
            requires_polya_for_construction: false,
            simple_models_mapq_cutoff: 30,
            polya_percentage_threshold: 0.7,
            low_polya_percentage_threshold: 0.1,
            indel_near_splice_site_dist: 10,
            upstream_region_len: 20,
            min_mono_exon_coverage: 0.75,
            multimap_strategy: MultimapStrategy::TakeBest,
            use_tss_model: true,
            create_alt_ts_nics: true,
            novel_apa: false,
            splice_site_correction: true,
        }
    }

    /// `--data_type nanopore --model_construction_strategy sensitive_ont`: matching
    /// `default` (`isoquant.py:1022`), splice correction `default_ont` (`:1032`),
    /// `resolve_ambiguous` = strategy default `monoexon_and_fsm` (`:1036`, `:1128-1129`).
    pub fn ont() -> Self {
        let model = ModelStrategy::sensitive_ont();
        Params {
            data_type: DataType::Ont,
            fl_data: false,
            matching: MatchingParams::default_ont(),
            splice: SpliceCorrection::default_ont(),
            report_novel_unspliced: model.novel_monoexonic,
            report_canonical: model.report_canonical,
            model,
            use_tss_model: false,
            ..Self::pacbio()
        }
    }

    /// Delta of the single-model assigners of the final counter:
    /// `max(delta, graph_clustering_distance)` (design §7.11).
    pub fn final_single_delta(&self) -> i32 {
        self.matching.delta.max(self.model.graph_clustering_distance)
    }

    /// Builds `Params` from a 4.1 `params.json` (schema §2.3): starts from the preset
    /// matching `data_type`/`fl_data` and overrides every field present in the JSON.
    /// Unknown or missing keys are ignored.
    pub fn from_params_json(j: &serde_json::Map<String, serde_json::Value>) -> Result<Self, String> {
        use serde_json::Value;
        let dt = match j.get("data_type").and_then(Value::as_str) {
            Some("nanopore") => DataType::Ont,
            Some("pacbio_ccs") | None => DataType::Pacbio,
            Some(o) => return Err(format!("unsupported data_type {o:?}")),
        };
        let mut p = Params::preset(dt);
        let i = |k: &str| j.get(k).and_then(Value::as_i64);
        let f = |k: &str| j.get(k).and_then(Value::as_f64);
        let b = |k: &str| j.get(k).and_then(Value::as_bool);
        macro_rules! set {
            ($dst:expr, i, $k:literal) => {
                if let Some(v) = i($k) {
                    $dst = v as _;
                }
            };
            ($dst:expr, f, $k:literal) => {
                if let Some(v) = f($k) {
                    $dst = v;
                }
            };
            ($dst:expr, b, $k:literal) => {
                if let Some(v) = b($k) {
                    $dst = v;
                }
            };
        }
        set!(p.fl_data, b, "fl_data");
        let m = &mut p.matching;
        set!(m.delta, i, "delta");
        set!(m.max_intron_shift, i, "max_intron_shift");
        set!(m.max_missed_exon_len, i, "max_missed_exon_len");
        set!(m.max_fake_terminal_exon_len, i, "max_fake_terminal_exon_len");
        set!(m.max_suspicious_intron_abs_len, i, "max_suspicious_intron_abs_len");
        set!(m.max_suspicious_intron_rel_len, f, "max_suspicious_intron_rel_len");
        set!(m.correct_minor_errors, b, "correct_minor_errors");
        set!(m.minor_exon_extension, i, "minor_exon_extension");
        set!(m.major_exon_extension, i, "major_exon_extension");
        set!(m.min_abs_exon_overlap, i, "min_abs_exon_overlap");
        set!(m.min_rel_exon_overlap, f, "min_rel_exon_overlap");
        set!(m.micro_intron_length, i, "micro_intron_length");
        set!(m.max_intron_abs_diff, i, "max_intron_abs_diff");
        set!(m.max_intron_rel_diff, f, "max_intron_rel_diff");
        set!(m.apa_delta, i, "apa_delta");
        set!(m.minimal_exon_overlap, i, "minimal_exon_overlap");
        set!(m.minimal_intron_absence_overlap, i, "minimal_intron_absence_overlap");
        set!(m.polya_window, i, "polya_window");
        set!(m.polya_fraction, f, "polya_fraction");
        set!(m.count_exons, b, "count_exons");
        if let Some(v) = j.get("resolve_ambiguous") {
            m.resolve_ambiguous = parse_resolve_ambiguous(v)?;
        }
        let s = &mut p.splice;
        set!(s.fuzzy_junctions, b, "correct_fuzzy_junctions");
        set!(s.intron_shifts, b, "correct_intron_shifts");
        set!(s.skipped_exons, b, "correct_skipped_exons");
        set!(s.terminal_exons, b, "correct_terminal_exons");
        set!(s.fake_terminal_exons, b, "correct_fake_terminal_exons");
        set!(s.microintron_retention, b, "correct_microintron_retention");
        let ms = &mut p.model;
        set!(ms.min_novel_intron_count, i, "min_novel_intron_count");
        set!(ms.graph_clustering_ratio, f, "graph_clustering_ratio");
        set!(ms.graph_clustering_distance, i, "graph_clustering_distance");
        set!(ms.min_novel_isolated_intron_abs, i, "min_novel_isolated_intron_abs");
        set!(ms.min_novel_isolated_intron_rel, f, "min_novel_isolated_intron_rel");
        set!(ms.terminal_position_abs, i, "terminal_position_abs");
        set!(ms.terminal_position_rel, f, "terminal_position_rel");
        set!(ms.terminal_internal_position_rel, f, "terminal_internal_position_rel");
        set!(ms.min_known_count, i, "min_known_count");
        set!(ms.min_nonfl_count, i, "min_nonfl_count");
        set!(ms.min_novel_count, i, "min_novel_count");
        set!(ms.min_novel_count_rel, f, "min_novel_count_rel");
        set!(ms.min_mono_count_rel, f, "min_mono_count_rel");
        set!(ms.singleton_adjacent_cov, i, "singleton_adjacent_cov");
        set!(ms.fl_only, b, "fl_only");
        set!(ms.require_monointronic_polya, b, "require_monointronic_polya");
        set!(ms.require_monoexonic_polya, b, "require_monoexonic_polya");
        let fl = &mut p.filter;
        if let Some(v) = j.get("min_mapq") {
            fl.min_mapq = v.as_i64().map(|x| x as u8);
        }
        set!(fl.inconsistent_mapq_cutoff, i, "inconsistent_mapq_cutoff");
        set!(fl.simple_alignments_mapq_cutoff, i, "simple_alignments_mapq_cutoff");
        set!(fl.max_coverage_small_chr, i, "max_coverage_small_chr");
        set!(fl.max_coverage_normal_chr, i, "max_coverage_normal_chr");
        set!(p.report_novel_unspliced, b, "report_novel_unspliced");
        set!(p.requires_polya_for_construction, b, "requires_polya_for_construction");
        set!(p.simple_models_mapq_cutoff, i, "simple_models_mapq_cutoff");
        set!(p.polya_percentage_threshold, f, "polya_percentage_threshold");
        set!(p.low_polya_percentage_threshold, f, "low_polya_percentage_threshold");
        set!(p.indel_near_splice_site_dist, i, "indel_near_splice_site_dist");
        set!(p.upstream_region_len, i, "upstream_region_len");
        set!(p.min_mono_exon_coverage, f, "min_mono_exon_coverage");
        set!(p.novel_apa, b, "novel_apa");
        if let Some(v) = i("report_canonical_strategy") {
            p.report_canonical = u8::try_from(v)
                .ok()
                .and_then(StrandnessReporting::from_code)
                .ok_or_else(|| format!("bad report_canonical_strategy {v}"))?;
            p.model.report_canonical = p.report_canonical;
        }
        if let Some(v) = i("polya_requirement_strategy") {
            p.polya_requirement = u8::try_from(v)
                .ok()
                .and_then(PolyAUsage::from_code)
                .ok_or_else(|| format!("bad polya_requirement_strategy {v}"))?;
        }
        if let Some(v) = i("multimap_strategy") {
            p.multimap_strategy = u8::try_from(v)
                .ok()
                .and_then(MultimapStrategy::from_code)
                .ok_or_else(|| format!("bad multimap_strategy {v}"))?;
        }
        if let Some(v) = b("no_splice_site_correction") {
            p.splice_site_correction = !v;
        }
        p.use_tss_model = p.fl_data;
        Ok(p)
    }
}

/// `resolve_ambiguous` as found in fixtures: the enum code (int) or the member name (str).
pub fn parse_resolve_ambiguous(v: &serde_json::Value) -> Result<AmbiguityResolving, String> {
    if let Some(c) = v.as_i64() {
        i8::try_from(c).ok().and_then(AmbiguityResolving::from_code)
    } else {
        v.as_str().and_then(AmbiguityResolving::from_py_name)
    }
    .ok_or_else(|| format!("bad resolve_ambiguous {v}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn presets_match_design_table() {
        let p = Params::pacbio();
        let m = &p.model;
        assert_eq!((p.matching.delta, p.matching.max_intron_shift, p.matching.max_missed_exon_len), (4, 30, 50));
        assert_eq!(p.matching.max_fake_terminal_exon_len, 20);
        assert_eq!(p.matching.max_intron_abs_diff, 30);
        assert_eq!((m.graph_clustering_distance, m.min_novel_isolated_intron_abs, m.min_novel_count), (5, 2, 2));
        assert_eq!((m.min_mono_count_rel, m.singleton_adjacent_cov), (0.001, 100));
        assert!(p.use_tss_model && p.fl_data && !p.splice.fake_terminal_exons);
        assert_eq!(p.final_single_delta(), 5);

        let o = Params::ont();
        let m = &o.model;
        assert_eq!((o.matching.delta, o.matching.max_intron_shift, o.matching.max_missed_exon_len), (6, 60, 100));
        assert_eq!((o.matching.max_suspicious_intron_abs_len, o.matching.max_suspicious_intron_rel_len), (60, 1.0));
        assert_eq!((m.graph_clustering_distance, m.min_novel_isolated_intron_abs, m.min_novel_count), (20, 3, 3));
        assert_eq!((m.min_mono_count_rel, m.singleton_adjacent_cov), (0.005, 10));
        assert_eq!(m.min_nonfl_count, 2);
        assert!(!o.use_tss_model && !o.fl_data && o.splice.fake_terminal_exons);
        assert_eq!(o.final_single_delta(), 20);
        assert!(o.report_novel_unspliced && p.report_novel_unspliced);
        assert_eq!(o.report_canonical, StrandnessReporting::OnlyStranded);
    }

    #[test]
    fn from_params_json_overrides() {
        let j: serde_json::Value = serde_json::json!({
            "data_type": "nanopore", "delta": 7, "resolve_ambiguous": 20,
            "resolve_ambiguous__name": "monoexon_and_fsm", "report_canonical_strategy": 1,
            "min_mapq": null, "fl_data": false
        });
        let p = Params::from_params_json(j.as_object().unwrap()).unwrap();
        assert_eq!(p.matching.delta, 7);
        assert_eq!(p.report_canonical, StrandnessReporting::OnlyCanonical);
        assert_eq!(p.model.graph_clustering_distance, 20);
        let q = Params::from_params_json(serde_json::json!({"resolve_ambiguous": "all"}).as_object().unwrap()).unwrap();
        assert_eq!(q.matching.resolve_ambiguous, AmbiguityResolving::All);
    }

    #[test]
    fn serde_roundtrip() {
        for p in [Params::pacbio(), Params::ont()] {
            let j = serde_json::to_string(&p).unwrap();
            assert_eq!(serde_json::from_str::<Params>(&j).unwrap(), p);
        }
    }
}
