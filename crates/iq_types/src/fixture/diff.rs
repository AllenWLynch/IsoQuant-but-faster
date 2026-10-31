//! Generic structural diff for fixture outputs: both sides are serialized to
//! `serde_json::Value` and walked together, producing field-level entries such as
//! `.assignment.matches[0].events[1].t: got 23, want 25`.

use serde::Serialize;
use serde_json::Value;
use std::fmt;

/// Options for [`diff_values`].
#[derive(Clone, Debug)]
pub struct DiffOptions {
    /// Floats compare equal if `|a-b| <= abs_tol` or `|a-b| <= rel_tol * max(|a|,|b|)`.
    /// Default 0/0 = exact equality (fixtures are bit-exact unless a unit says otherwise).
    pub abs_tol: f64,
    pub rel_tol: f64,
    /// Maximum entries kept; the rest are counted in `Diff::truncated`.
    pub max_entries: usize,
    /// Compare arrays at these paths (exact path match, e.g. `.isoforms`) as multisets.
    pub unordered_paths: Vec<String>,
}

impl Default for DiffOptions {
    fn default() -> Self {
        DiffOptions { abs_tol: 0.0, rel_tol: 0.0, max_entries: 20, unordered_paths: Vec::new() }
    }
}

/// One differing field.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffEntry {
    pub path: String,
    pub got: String,
    pub want: String,
}

/// A readable mismatch report.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Diff {
    pub entries: Vec<DiffEntry>,
    /// Entries beyond `max_entries`.
    pub truncated: usize,
    /// Set when the case could not be run (panic, error, conversion failure).
    pub error: Option<String>,
}

impl Diff {
    pub fn error(msg: impl Into<String>) -> Self {
        Diff { error: Some(msg.into()), ..Default::default() }
    }
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty() && self.truncated == 0 && self.error.is_none()
    }
    pub fn len(&self) -> usize {
        self.entries.len() + self.truncated
    }
}

impl fmt::Display for Diff {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(e) = &self.error {
            writeln!(f, "  error: {e}")?;
        }
        for e in &self.entries {
            let p = if e.path.is_empty() { "<root>" } else { &e.path };
            writeln!(f, "  {p}: got {}, want {}", e.got, e.want)?;
        }
        if self.truncated > 0 {
            writeln!(f, "  ... and {} more", self.truncated)?;
        }
        Ok(())
    }
}

fn short(v: &Value) -> String {
    let s = super::pyjson::to_py_string(v).unwrap_or_else(|_| v.to_string());
    if s.len() > 120 {
        let mut cut = 117;
        while !s.is_char_boundary(cut) {
            cut -= 1;
        }
        format!("{}...", &s[..cut])
    } else {
        s
    }
}

struct Walker<'o> {
    opts: &'o DiffOptions,
    out: Diff,
}

impl Walker<'_> {
    fn push(&mut self, path: &str, got: String, want: String) {
        if self.out.entries.len() < self.opts.max_entries {
            self.out.entries.push(DiffEntry { path: path.to_string(), got, want });
        } else {
            self.out.truncated += 1;
        }
    }

    fn num_eq(&self, a: &serde_json::Number, b: &serde_json::Number) -> bool {
        if a == b {
            return true;
        }
        match (a.as_f64(), b.as_f64()) {
            (Some(x), Some(y)) if a.is_f64() || b.is_f64() => {
                if x == y {
                    return a.is_f64() == b.is_f64() || (self.opts.abs_tol > 0.0 || self.opts.rel_tol > 0.0);
                }
                let d = (x - y).abs();
                d <= self.opts.abs_tol || d <= self.opts.rel_tol * x.abs().max(y.abs())
            }
            _ => false,
        }
    }

    fn walk(&mut self, path: &mut String, got: &Value, want: &Value) {
        match (got, want) {
            (Value::Object(g), Value::Object(w)) => {
                for (k, wv) in w {
                    let len = path.len();
                    path.push('.');
                    path.push_str(k);
                    match g.get(k) {
                        Some(gv) => self.walk(path, gv, wv),
                        None => self.push(path, "<missing>".into(), short(wv)),
                    }
                    path.truncate(len);
                }
                for (k, gv) in g {
                    if !w.contains_key(k) {
                        let len = path.len();
                        path.push('.');
                        path.push_str(k);
                        self.push(path, short(gv), "<missing>".into());
                        path.truncate(len);
                    }
                }
            }
            (Value::Array(g), Value::Array(w)) => {
                if self.opts.unordered_paths.iter().any(|p| p == path) {
                    let mut gs: Vec<String> = g.iter().map(short).collect();
                    let mut ws: Vec<String> = w.iter().map(short).collect();
                    gs.sort();
                    ws.sort();
                    if gs != ws {
                        self.push(path, short(got), short(want));
                    }
                    return;
                }
                if g.len() != w.len() {
                    self.push(&format!("{path}.len()"), g.len().to_string(), w.len().to_string());
                }
                for (i, (gv, wv)) in g.iter().zip(w).enumerate() {
                    let len = path.len();
                    path.push_str(&format!("[{i}]"));
                    self.walk(path, gv, wv);
                    path.truncate(len);
                }
                for (i, wv) in w.iter().enumerate().skip(g.len()) {
                    self.push(&format!("{path}[{i}]"), "<missing>".into(), short(wv));
                }
                for (i, gv) in g.iter().enumerate().skip(w.len()) {
                    self.push(&format!("{path}[{i}]"), short(gv), "<missing>".into());
                }
            }
            (Value::Number(a), Value::Number(b)) => {
                if !self.num_eq(a, b) {
                    self.push(path, short(got), short(want));
                }
            }
            _ => {
                if got != want {
                    self.push(path, short(got), short(want));
                }
            }
        }
    }
}

/// Field-level diff of two JSON values. `None` if equal under `opts`.
pub fn diff_json(got: &Value, want: &Value, opts: &DiffOptions) -> Option<Diff> {
    let mut w = Walker { opts, out: Diff::default() };
    w.walk(&mut String::new(), got, want);
    if w.out.is_empty() {
        None
    } else {
        Some(w.out)
    }
}

/// Field-level diff of any two serializable values (typically fixture DTOs).
pub fn diff_values<G: Serialize + ?Sized, W: Serialize + ?Sized>(
    got: &G,
    want: &W,
    opts: &DiffOptions,
) -> Option<Diff> {
    match (serde_json::to_value(got), serde_json::to_value(want)) {
        (Ok(g), Ok(w)) => diff_json(&g, &w, opts),
        (Err(e), _) | (_, Err(e)) => Some(Diff::error(format!("serialize for diff: {e}"))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn reports_paths() {
        let g = json!({"a": 1, "m": [{"t": 23, "p": 0.1}, {"t": 1}], "x": "s"});
        let w = json!({"a": 1, "m": [{"t": 25, "p": 0.1}], "y": true});
        let d = diff_json(&g, &w, &DiffOptions::default()).unwrap();
        let paths: Vec<_> = d.entries.iter().map(|e| e.path.as_str()).collect();
        assert_eq!(paths, vec![".m.len()", ".m[0].t", ".m[1]", ".y", ".x"]);
        assert!(d.to_string().contains(".m[0].t: got 23, want 25"));
        assert!(diff_json(&g, &g, &DiffOptions::default()).is_none());
    }

    #[test]
    fn floats_and_tolerance() {
        let g = json!({"f": 1.0000000001, "i": 3});
        let w = json!({"f": 1.0, "i": 3.0});
        let d = diff_json(&g, &w, &DiffOptions::default()).unwrap();
        assert_eq!(d.entries.len(), 2); // int vs float is a difference in exact mode
        let o = DiffOptions { rel_tol: 1e-9, ..Default::default() };
        assert!(diff_json(&g, &w, &o).is_none());
    }

    #[test]
    fn unordered() {
        let o = DiffOptions { unordered_paths: vec![".s".into()], ..Default::default() };
        assert!(diff_json(&json!({"s": [1, 2]}), &json!({"s": [2, 1]}), &o).is_none());
        assert!(diff_json(&json!({"s": [1, 2]}), &json!({"s": [2, 3]}), &o).is_some());
    }

    #[test]
    fn truncation() {
        let g = json!((0..30).collect::<Vec<_>>());
        let w = json!((1..31).collect::<Vec<_>>());
        let d = diff_json(&g, &w, &DiffOptions::default()).unwrap();
        assert_eq!((d.entries.len(), d.truncated), (20, 10));
    }
}
