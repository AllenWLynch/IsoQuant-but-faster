//! U10a: native XGBoost tree evaluation in f32 (design §7.9).
//!
//! Ports what `XGBClassifier.predict` / `predict_proba` / `get_booster().predict(output_margin=True)`
//! compute for IsoQuant's two bundled `binary:logistic` gbtree models (xgboost 3.4.2):
//! - features f64 → f32 (round to nearest); NaN = missing → `default_left`;
//! - a split goes left iff `f32(x) < split_conditions[n]`; a leaf (`left_children == -1`)
//!   contributes `split_conditions[leaf]`;
//! - margin = base margin + leaves in tree order, summed in f32;
//! - probability = xgboost `common::Sigmoid` in f32; label = `p > 0.5`.

use iq_types::serde_json::{self, Value};

/// Number of features in both models (`FEATURE_COLUMNS`).
pub const N_FEATURES: usize = 8;

const NO_CHILD: u32 = u32::MAX;

/// One flattened tree node (absolute indices into `XgbModel::nodes`).
#[derive(Clone, Copy, Debug)]
struct Node {
    /// Split threshold, or the leaf value when `left == NO_CHILD`.
    value: f32,
    left: u32,
    right: u32,
    feature: u8,
    default_left: bool,
}

/// One loaded `binary:logistic` model.
#[derive(Debug)]
pub struct XgbModel {
    nodes: Vec<Node>,
    /// Root node index of each tree, in tree order.
    roots: Vec<u32>,
    base_margin: f32,
}

/// Base margins as xgboost 3.4.2 computes them for the bundled models (`base_score` string →
/// f32 bits). Recomputing `-log(1/b - 1)` with other log routines is off by 1-3 ulp, which
/// flips many margins (design §7.9); these bits were checked against the fixtures' margins.
const KNOWN_BASE_MARGINS: &[(&str, u32)] = &[("[5.7763106E-1]", 0x3ea0_48e8), ("[3.413823E-1]", 0xbf28_3a5a)];

fn field<'a>(v: &'a Value, k: &str) -> Result<&'a Value, String> {
    v.get(k).ok_or_else(|| format!("xgb model: missing field `{k}`"))
}

fn arr<'a>(v: &'a Value, k: &str) -> Result<&'a Vec<Value>, String> {
    field(v, k)?.as_array().ok_or_else(|| format!("xgb model: `{k}` is not an array"))
}

fn ints(v: &Value, k: &str) -> Result<Vec<i64>, String> {
    arr(v, k)?.iter().map(|x| x.as_i64().ok_or_else(|| format!("xgb model: `{k}` has a non-integer"))).collect()
}

impl XgbModel {
    /// `XgbModel::load(json)`: parses an xgboost JSON model (`save_model("*.json")`).
    ///
    /// Only what IsoQuant's models use is supported: gbtree, one target, numeric splits, no
    /// categorical features; anything else is rejected.
    pub fn load(json: &[u8]) -> Result<XgbModel, String> {
        let root: Value = serde_json::from_slice(json).map_err(|e| format!("xgb model: {e}"))?;
        let learner = field(&root, "learner")?;
        let objective = field(field(learner, "objective")?, "name")?.as_str().unwrap_or("");
        if objective != "binary:logistic" {
            return Err(format!("xgb model: unsupported objective {objective:?}"));
        }
        let lmp = field(learner, "learner_model_param")?;
        let base_score = field(lmp, "base_score")?.as_str().ok_or("xgb model: base_score not a string")?;
        let base_margin = match KNOWN_BASE_MARGINS.iter().find(|(s, _)| *s == base_score) {
            Some(&(_, bits)) => f32::from_bits(bits),
            None => {
                // Not one of the bundled models: best effort (may differ by an ulp from xgboost).
                let b: f32 = base_score
                    .trim_matches(|c| c == '[' || c == ']')
                    .parse()
                    .map_err(|e| format!("xgb model: base_score {base_score:?}: {e}"))?;
                -(1.0f32 / b - 1.0).ln()
            }
        };
        let booster = field(learner, "gradient_booster")?;
        if field(booster, "name")?.as_str() != Some("gbtree") {
            return Err("xgb model: only gbtree is supported".into());
        }
        let trees = arr(field(booster, "model")?, "trees")?;
        let mut nodes = Vec::new();
        let mut roots = Vec::with_capacity(trees.len());
        for (t, tree) in trees.iter().enumerate() {
            let lc = ints(tree, "left_children")?;
            let rc = ints(tree, "right_children")?;
            let si = ints(tree, "split_indices")?;
            let dl = ints(tree, "default_left")?;
            let st = ints(tree, "split_type")?;
            let sc: Vec<f32> = arr(tree, "split_conditions")?
                .iter()
                // Every split condition of the bundled models is the shortest f32 repr, and
                // decimal → f64 → f32 equals decimal → f32 for all of them (checked once in Python).
                .map(|x| x.as_f64().map(|v| v as f32).ok_or_else(|| "xgb model: bad split_conditions".to_string()))
                .collect::<Result<_, _>>()?;
            let n = lc.len();
            if [rc.len(), si.len(), dl.len(), st.len(), sc.len()].iter().any(|&l| l != n) || n == 0 {
                return Err(format!("xgb model: tree {t} has inconsistent array lengths"));
            }
            if !arr(tree, "categories")?.is_empty() || st.iter().any(|&s| s != 0) {
                return Err(format!("xgb model: tree {t} has categorical splits"));
            }
            let off = nodes.len() as u32;
            roots.push(off);
            for i in 0..n {
                let leaf = lc[i] == -1;
                let child = |c: i64| -> Result<u32, String> {
                    if leaf {
                        Ok(NO_CHILD)
                    } else if (0..n as i64).contains(&c) {
                        Ok(off + c as u32)
                    } else {
                        Err(format!("xgb model: tree {t} node {i} has a bad child {c}"))
                    }
                };
                let feature = if leaf { 0 } else { si[i] };
                if !(0..N_FEATURES as i64).contains(&feature) {
                    return Err(format!("xgb model: tree {t} node {i} splits on feature {feature}"));
                }
                nodes.push(Node {
                    value: sc[i],
                    left: child(lc[i])?,
                    right: child(rc[i])?,
                    feature: feature as u8,
                    default_left: dl[i] != 0,
                });
            }
        }
        Ok(XgbModel { nodes, roots, base_margin })
    }

    /// Raw margin (f32, trees summed in order from the hard-coded base margin).
    pub fn margin(&self, features: &[f64; N_FEATURES]) -> f32 {
        let x: [f32; N_FEATURES] = features.map(|v| v as f32);
        let mut m = self.base_margin;
        for &r in &self.roots {
            let mut n = self.nodes[r as usize];
            while n.left != NO_CHILD {
                let v = x[n.feature as usize];
                let go_left = if v.is_nan() { n.default_left } else { v < n.value };
                n = self.nodes[if go_left { n.left } else { n.right } as usize];
            }
            m += n.value;
        }
        m
    }

    /// `predict_proba(...)[:, 1]` for one row.
    pub fn proba(&self, features: &[f64; N_FEATURES]) -> f32 {
        sigmoid(self.margin(features))
    }

    /// `predict(...)` for one row: `proba > 0.5`.
    pub fn predict(&self, features: &[f64; N_FEATURES]) -> bool {
        self.proba(features) > 0.5
    }

    /// Number of trees.
    pub fn n_trees(&self) -> usize {
        self.roots.len()
    }
}

/// xgboost 3.4.2 `common::Sigmoid` (`src/common/math.h`), in f32.
pub fn sigmoid(x: f32) -> f32 {
    const EPS: f32 = 1e-16;
    let x = (-x).min(88.7);
    let denom = x.exp() + 1.0 + EPS;
    1.0 / denom
}

/// The two bundled models.
#[derive(Debug)]
pub struct XgbModels {
    pub polya: XgbModel,
    pub tss: XgbModel,
}

impl XgbModels {
    /// Loads the bundled `MODEL_POLYA_JSON` / `MODEL_TSS_JSON`.
    pub fn bundled() -> Result<XgbModels, String> {
        Ok(XgbModels { polya: XgbModel::load(crate::MODEL_POLYA_JSON)?, tss: XgbModel::load(crate::MODEL_TSS_JSON)? })
    }

    /// Model by fixture name (`"polya"` / `"tss"`).
    pub fn by_name(&self, name: &str) -> Option<&XgbModel> {
        match name {
            "polya" => Some(&self.polya),
            "tss" => Some(&self.tss),
            _ => None,
        }
    }
}

/// Process-wide bundled models (parsed once).
pub fn bundled_models() -> &'static XgbModels {
    static M: std::sync::OnceLock<XgbModels> = std::sync::OnceLock::new();
    M.get_or_init(|| XgbModels::bundled().expect("bundled xgboost models parse"))
}

pub mod fixture {
    use super::*;
    use iq_types::fixture::schema::{U10aXgbPredictIn, U10aXgbPredictOut};
    use iq_types::fixture::{Diff, DiffEntry, FixResult, FixtureCtx, PyF64, UnitFixture, Verdict};

    /// Converts a fixture feature row.
    pub(crate) fn row(f: &[PyF64]) -> FixResult<[f64; N_FEATURES]> {
        if f.len() != N_FEATURES {
            return Err(format!("feature row has {} values", f.len()).into());
        }
        Ok(std::array::from_fn(|i| f[i].0))
    }

    /// Judges one model output (margin, proba, label) against 4.1 (§3 U10a pass criteria):
    /// probabilities |Δ| ≤ 1e-6, labels exact unless |margin| < 1e-5. Margins are expected
    /// bit-equal; a margin difference alone (within the proba tolerance) is still reported.
    pub(crate) fn judge_pred(path: &str, got: (f64, f64, i32), want: (f64, f64, i32), out: &mut Vec<DiffEntry>) {
        let (gm, gp, gl) = got;
        let (wm, wp, wl) = want;
        let push = |out: &mut Vec<DiffEntry>, f: &str, g: String, w: String| {
            out.push(DiffEntry { path: format!("{path}.{f}"), got: g, want: w })
        };
        if (gp - wp).abs() > 1e-6 || gp.is_nan() != wp.is_nan() {
            push(out, "proba", gp.to_string(), wp.to_string());
        }
        if gm.to_bits() != wm.to_bits() && (gm - wm).abs() > 1e-6 {
            push(out, "margin", gm.to_string(), wm.to_string());
        }
        if gl != wl && wm.abs() >= 1e-5 {
            push(out, "label", gl.to_string(), wl.to_string());
        }
    }

    pub(crate) fn verdict(entries: Vec<DiffEntry>) -> Verdict {
        if entries.is_empty() {
            Verdict::Exact
        } else {
            Verdict::Mismatch(Diff { entries, ..Default::default() })
        }
    }

    pub struct U10aXgbPredict;
    impl UnitFixture for U10aXgbPredict {
        const UNIT: &'static str = "U10a";
        const FN: &'static str = "xgb_predict";
        type In = U10aXgbPredictIn;
        type Out = U10aXgbPredictOut;

        fn run(input: &Self::In, _ctx: &FixtureCtx) -> FixResult<Self::Out> {
            let model =
                bundled_models().by_name(&input.model).ok_or_else(|| format!("unknown model {:?}", input.model))?;
            let mut out = U10aXgbPredictOut { labels: Vec::new(), proba: Vec::new(), margin: Vec::new() };
            for f in &input.features {
                let m = model.margin(&row(f)?);
                let p = sigmoid(m);
                out.labels.push((p > 0.5) as i32);
                out.proba.push(p as f64);
                out.margin.push(m as f64);
            }
            Ok(out)
        }

        fn compare(got: &Self::Out, want: &Self::Out, _ctx: &FixtureCtx) -> Verdict {
            if got == want {
                return Verdict::Exact;
            }
            let mut e = Vec::new();
            if got.labels.len() != want.labels.len() {
                return Verdict::Mismatch(Diff::error("row count differs"));
            }
            for i in 0..want.labels.len() {
                judge_pred(
                    &format!("[{i}]"),
                    (got.margin[i], got.proba[i], got.labels[i]),
                    (want.margin[i], want.proba[i], want.labels[i]),
                    &mut e,
                );
            }
            verdict(e)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bundled_models_load() {
        let m = XgbModels::bundled().unwrap();
        assert_eq!(m.polya.n_trees(), 417);
        assert_eq!(m.tss.n_trees(), 417);
    }

    #[test]
    fn sigmoid_matches_xgboost_clamp() {
        assert_eq!(sigmoid(0.0), 0.5);
        assert!(sigmoid(-1000.0) > 0.0);
        assert_eq!(sigmoid(1000.0), 1.0);
    }
}
