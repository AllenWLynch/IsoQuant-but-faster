//! Report model (metrics with thresholds, tables) and its markdown / JSON / TSV renderings.

use iq_types::serde_json::{self, json, Map, Value};

/// How a metric's deficit (distance from the ideal value) is judged.
#[derive(Clone, Debug)]
pub struct Metric {
    pub name: String,
    pub section: &'static str,
    pub value: f64,
    pub ref_value: Option<f64>,
    pub test_value: Option<f64>,
    /// Distance from perfect agreement (>= 0); `NaN` when undefined.
    pub deficit: f64,
    /// Maximum allowed deficit; `None` = informational.
    pub threshold: Option<f64>,
    /// `value` and `threshold` printed as: "frac", "pp", "rel", "corr", "n".
    pub unit: &'static str,
    pub desc: String,
    /// Filled by [`Report::judge`].
    pub floor_deficit: Option<f64>,
    pub pass: Option<bool>,
    /// Deviation that explains a failure of this gated metric (`--explained`); a failing
    /// metric with an explanation is reported as EXPLAINED, not FAIL.
    pub explained: Option<String>,
}

impl Metric {
    pub fn new(name: impl Into<String>, section: &'static str, unit: &'static str, value: f64, deficit: f64) -> Metric {
        Metric {
            name: name.into(),
            section,
            value,
            ref_value: None,
            test_value: None,
            deficit,
            threshold: None,
            unit,
            desc: String::new(),
            floor_deficit: None,
            pass: None,
            explained: None,
        }
    }
    pub fn thr(mut self, t: f64) -> Self {
        self.threshold = Some(t);
        self
    }
    pub fn sides(mut self, r: f64, t: f64) -> Self {
        self.ref_value = Some(r);
        self.test_value = Some(t);
        self
    }
    pub fn desc(mut self, d: impl Into<String>) -> Self {
        self.desc = d.into();
        self
    }
    /// `Some("pass" | "fail" | "explained")` for gated metrics.
    pub fn verdict(&self) -> Option<&'static str> {
        self.pass.map(|p| match (p, &self.explained) {
            (true, _) => "pass",
            (false, Some(_)) => "explained",
            (false, None) => "fail",
        })
    }
    fn verdict_md(&self) -> String {
        match (self.verdict(), &self.explained) {
            (Some("pass"), _) => "PASS".into(),
            (Some("explained"), Some(d)) => format!("EXPLAINED ({d})"),
            (Some(_), _) => "**FAIL**".into(),
            (None, _) => String::new(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Table {
    pub section: &'static str,
    pub title: String,
    pub headers: Vec<String>,
    pub rows: Vec<Vec<String>>,
    /// Rows shown in markdown (all rows go to JSON).
    pub md_limit: usize,
}

impl Table {
    pub fn new(section: &'static str, title: impl Into<String>, headers: &[&str]) -> Table {
        Table {
            section,
            title: title.into(),
            headers: headers.iter().map(|s| s.to_string()).collect(),
            rows: Vec::new(),
            md_limit: usize::MAX,
        }
    }
    pub fn row(&mut self, r: Vec<String>) {
        self.rows.push(r);
    }
}

#[derive(Default)]
pub struct Report {
    pub title: String,
    pub header_lines: Vec<String>,
    pub metrics: Vec<Metric>,
    pub tables: Vec<Table>,
    pub warnings: Vec<String>,
    pub floor_source: Option<String>,
    /// `--explained` rules: (metric name or name prefix, deviation id).
    pub explained: Vec<(String, String)>,
}

pub const SECTIONS: [(&str, &str); 9] = [
    ("overview", "Overview"),
    ("chain", "Intron-chain match (read-weighted)"),
    ("counts", "Count agreement over matched models"),
    ("fates", "Read fates"),
    ("tpg", "Transcripts per gene"),
    ("novelty", "Novel / known"),
    ("grouped", "Grouped counts"),
    ("breakdown", "Breakdowns"),
    ("residuals", "Residuals"),
];

pub fn fmt_val(v: f64, unit: &str) -> String {
    if v.is_nan() {
        return "NA".into();
    }
    match unit {
        "frac" => format!("{:.4}", v),
        "corr" => format!("{:.5}", v),
        "pp" => format!("{:+.3} pp", v),
        "rel" => format!("{:+.3}%", v * 100.0),
        "n" => {
            if v.fract() == 0.0 {
                format!("{}", v as i64)
            } else {
                format!("{:.2}", v)
            }
        }
        _ => format!("{:.4}", v),
    }
}

fn fmt_dev(v: f64, unit: &str) -> String {
    if v.is_nan() {
        return "NA".into();
    }
    match unit {
        "pp" => format!("{:.3} pp", v),
        "rel" => format!("{:.3}%", v * 100.0),
        "n" => format!("{}", v),
        _ => format!("{:.5}", v),
    }
}

pub fn jnum(v: f64) -> Value {
    if v.is_finite() {
        json!(v)
    } else {
        Value::Null
    }
}

impl Report {
    /// Sets `pass` on gated metrics. With a noise floor, the allowed deficit is
    /// `threshold + floor deficit` of the same-named metric in the floor report.
    pub fn judge(&mut self, floor: Option<&Value>) {
        let mut fmap = std::collections::HashMap::new();
        if let Some(Value::Array(ms)) = floor.and_then(|f| f.get("metrics")) {
            for m in ms {
                if let (Some(n), Some(d)) =
                    (m.get("name").and_then(Value::as_str), m.get("deficit").and_then(Value::as_f64))
                {
                    fmap.insert(n.to_string(), d);
                }
            }
        }
        for m in &mut self.metrics {
            m.floor_deficit = fmap.get(&m.name).copied();
            if let Some(t) = m.threshold {
                let allowed = t + m.floor_deficit.unwrap_or(0.0);
                m.pass = Some(!m.deficit.is_nan() && m.deficit <= allowed + 1e-12);
                m.explained =
                    self.explained.iter().find(|(k, _)| m.name.starts_with(k.as_str())).map(|(_, d)| d.clone());
            }
        }
    }

    /// (passing, gated) metric counts.
    pub fn n_gated(&self) -> (usize, usize) {
        let g: Vec<_> = self.metrics.iter().filter_map(|m| m.pass).collect();
        (g.iter().filter(|&&p| p).count(), g.len())
    }

    /// Gated metrics that fail but are explained by a named deviation.
    pub fn n_explained(&self) -> usize {
        self.metrics.iter().filter(|m| m.verdict() == Some("explained")).count()
    }

    /// One-line verdict: "15/15 pass" or "14/15 pass, 1 explained (D-M6), 0 fail".
    pub fn gate_line(&self) -> String {
        let (p, n) = self.n_gated();
        let e = self.n_explained();
        if e == 0 {
            return format!("{p}/{n} pass");
        }
        let mut devs: Vec<&str> = self
            .metrics
            .iter()
            .filter(|m| m.verdict() == Some("explained"))
            .filter_map(|m| m.explained.as_deref())
            .collect();
        devs.sort_unstable();
        devs.dedup();
        format!("{p}/{n} pass, {e} explained ({}), {} fail", devs.join(", "), n - p - e)
    }

    pub fn markdown(&self) -> String {
        let mut o = String::new();
        o += &format!("# {}\n\n", self.title);
        for l in &self.header_lines {
            o += &format!("- {l}\n");
        }
        o += &format!(
            "\n**Gated metrics: {}**{}\n\n",
            self.gate_line(),
            match &self.floor_source {
                Some(f) => format!(" (thresholds widened by the noise floor `{f}`)"),
                None => " (absolute thresholds; no noise floor given)".into(),
            }
        );
        if !self.warnings.is_empty() {
            o += "## Warnings\n\n";
            for w in &self.warnings {
                o += &format!("- {w}\n");
            }
            o += "\n";
        }
        for (sec, title) in SECTIONS {
            let ms: Vec<&Metric> = self.metrics.iter().filter(|m| m.section == sec).collect();
            let ts: Vec<&Table> = self.tables.iter().filter(|t| t.section == sec).collect();
            if ms.is_empty() && ts.is_empty() {
                continue;
            }
            o += &format!("## {title}\n\n");
            if !ms.is_empty() {
                o += "| metric | ref | test | value | deviation | threshold | floor | verdict |\n|---|---|---|---|---|---|---|---|\n";
                for m in ms {
                    let opt = |v: Option<f64>| {
                        v.map(|v| fmt_val(v, if m.unit == "pp" || m.unit == "rel" { "frac" } else { m.unit }))
                            .unwrap_or_default()
                    };
                    o += &format!(
                        "| {} | {} | {} | {} | {} | {} | {} | {} |\n",
                        m.name,
                        opt(m.ref_value),
                        opt(m.test_value),
                        fmt_val(m.value, m.unit),
                        fmt_dev(m.deficit, m.unit),
                        m.threshold.map(|t| format!("≤ {}", fmt_dev(t, m.unit))).unwrap_or_else(|| "info".into()),
                        m.floor_deficit.map(|f| fmt_dev(f, m.unit)).unwrap_or_default(),
                        m.verdict_md()
                    );
                }
                o += "\n";
            }
            for t in ts {
                o += &format!("**{}**", t.title);
                if t.rows.len() > t.md_limit {
                    o += &format!(" (top {} of {})", t.md_limit, t.rows.len());
                }
                o += "\n\n";
                o += &format!("| {} |\n|{}\n", t.headers.join(" | "), "---|".repeat(t.headers.len()));
                for r in t.rows.iter().take(t.md_limit) {
                    o += &format!("| {} |\n", r.iter().map(|c| c.replace('|', "\\|")).collect::<Vec<_>>().join(" | "));
                }
                o += "\n";
            }
        }
        o += "## Metric definitions\n\n";
        let mut seen = std::collections::HashSet::new();
        for m in &self.metrics {
            if !m.desc.is_empty() && seen.insert(m.desc.clone()) {
                o += &format!("- `{}`: {}\n", m.name, m.desc);
            }
        }
        o
    }

    pub fn json(&self) -> Value {
        let metrics: Vec<Value> = self
            .metrics
            .iter()
            .map(|m| {
                json!({
                    "name": m.name, "section": m.section, "unit": m.unit,
                    "value": jnum(m.value), "ref": m.ref_value.map(jnum), "test": m.test_value.map(jnum),
                    "deficit": jnum(m.deficit), "threshold": m.threshold, "floor_deficit": m.floor_deficit,
                    "gated": m.threshold.is_some(), "pass": m.pass, "verdict": m.verdict(),
                    "explained_by": if m.verdict() == Some("explained") { m.explained.clone() } else { None },
                    "desc": m.desc,
                })
            })
            .collect();
        let tables: Vec<Value> = self
            .tables
            .iter()
            .map(|t| {
                let rows: Vec<Value> = t
                    .rows
                    .iter()
                    .take(t.md_limit.max(1000))
                    .map(|r| {
                        let mut o = Map::new();
                        for (h, c) in t.headers.iter().zip(r) {
                            o.insert(h.clone(), Value::String(c.clone()));
                        }
                        Value::Object(o)
                    })
                    .collect();
                json!({"section": t.section, "title": t.title, "n_rows": t.rows.len(), "rows": rows})
            })
            .collect();
        let (p, n) = self.n_gated();
        let e = self.n_explained();
        let rules: Vec<Value> = self.explained.iter().map(|(k, d)| json!({"metric": k, "deviation": d})).collect();
        json!({
            "title": self.title, "info": self.header_lines, "warnings": self.warnings,
            "noise_floor": self.floor_source, "gated_pass": p, "gated_total": n, "all_pass": p == n,
            "gated_explained": e, "gated_fail": n - p - e, "all_pass_or_explained": p + e == n,
            "explained_rules": rules,
            "metrics": metrics, "tables": tables,
        })
    }

    pub fn tsv(&self) -> String {
        let mut o = String::from("metric\tvalue\tdeficit\tthreshold\tfloor_deficit\tpass\n");
        let f = |v: Option<f64>| v.map(|v| format!("{v}")).unwrap_or_else(|| "NA".into());
        for m in &self.metrics {
            o += &format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                m.name,
                m.value,
                m.deficit,
                f(m.threshold),
                f(m.floor_deficit),
                match (m.verdict(), &m.explained) {
                    (Some("pass"), _) => "PASS".to_string(),
                    (Some("explained"), Some(d)) => format!("EXPLAINED:{d}"),
                    (Some(_), _) => "FAIL".to_string(),
                    (None, _) => "info".to_string(),
                }
            );
        }
        o
    }
}

pub fn to_pretty(v: &Value) -> String {
    serde_json::to_string_pretty(v).unwrap_or_default()
}
