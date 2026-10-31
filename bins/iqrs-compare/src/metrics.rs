//! Metric computations.

use crate::annot::{LABEL_CLASSES, STRUCT_CLASSES};
use crate::io::{coarse_of, Model, Side, COARSE_FATES, FINE_FATES, NO_MODEL};
use crate::matching::{recip, Matching};
use crate::report::{Metric, Report, Table};
use rustc_hash::FxHashMap;

// ---------------------------------------------------------------- statistics

/// Average ranks (ties share the mean rank), 1-based.
pub fn ranks(v: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].total_cmp(&v[b]));
    let mut r = vec![0.0; v.len()];
    let mut i = 0;
    while i < idx.len() {
        let mut j = i;
        while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
            j += 1;
        }
        let avg = (i + j) as f64 / 2.0 + 1.0;
        for &k in &idx[i..=j] {
            r[k] = avg;
        }
        i = j + 1;
    }
    r
}

/// Pearson correlation; 1.0 for identical vectors (even constant ones), NaN if undefined.
pub fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len();
    if n == 0 {
        return f64::NAN;
    }
    if x == y {
        return 1.0;
    }
    let mx = x.iter().sum::<f64>() / n as f64;
    let my = y.iter().sum::<f64>() / n as f64;
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for i in 0..n {
        let (a, b) = (x[i] - mx, y[i] - my);
        sxy += a * b;
        sxx += a * a;
        syy += b * b;
    }
    if sxx == 0.0 || syy == 0.0 {
        return f64::NAN;
    }
    sxy / (sxx * syy).sqrt()
}

pub fn spearman(x: &[f64], y: &[f64]) -> f64 {
    if x == y && !x.is_empty() {
        return 1.0;
    }
    pearson(&ranks(x), &ranks(y))
}

fn pct(a: f64, b: f64) -> f64 {
    if b > 0.0 {
        a / b
    } else {
        f64::NAN
    }
}
fn f4(v: f64) -> String {
    if v.is_nan() {
        "NA".into()
    } else {
        format!("{v:.4}")
    }
}
fn cnt(v: f64) -> String {
    if v.fract() == 0.0 {
        format!("{}", v as i64)
    } else {
        format!("{v:.2}")
    }
}

// ---------------------------------------------------------------- context

pub struct Ctx<'a> {
    pub r: &'a Side,
    pub t: &'a Side,
    pub mt: &'a Matching,
    pub chroms: &'a [String],
    /// Per global model: locus category index into `cats` (with annotation).
    pub cat: Option<Vec<u8>>,
    pub cats: Vec<String>,
    /// Per global model: structural class (with annotation).
    pub sclass: Option<Vec<u8>>,
    /// `--exclude-nongene-loci`: non-`gene` locus categories are excluded from gated metrics.
    pub exclude_nongene: bool,
    pub top: usize,
}

const SIDES: [&str; 2] = ["ref", "test"];

impl<'a> Ctx<'a> {
    pub fn nr(&self) -> usize {
        self.mt.nr
    }
    pub fn n(&self) -> usize {
        self.mt.comp_of.len()
    }
    pub fn model(&self, g: usize) -> &'a Model {
        if g < self.nr() {
            &self.r.models[g]
        } else {
            &self.t.models[g - self.nr()]
        }
    }
    pub fn range(&self, side: usize) -> std::ops::Range<usize> {
        if side == 0 {
            0..self.nr()
        } else {
            self.nr()..self.n()
        }
    }
    pub fn side(&self, side: usize) -> &'a Side {
        if side == 0 {
            self.r
        } else {
            self.t
        }
    }
    pub fn matched(&self, g: usize) -> bool {
        self.mt.comps[self.mt.comp_of[g] as usize].matched()
    }
    pub fn is_nongene(&self, g: usize) -> bool {
        match &self.cat {
            Some(c) => {
                let n = &self.cats[c[g] as usize];
                n != "gene" && n != "intergenic"
            }
            None => false,
        }
    }
    /// Models that enter gated metrics (all, or only `gene`/intergenic loci with `--exclude-nongene-loci`).
    pub fn gated(&self, g: usize) -> bool {
        !(self.exclude_nongene && self.is_nongene(g))
    }
    fn cat_name(&self, g: usize) -> &str {
        match &self.cat {
            Some(c) => &self.cats[c[g] as usize],
            None => "",
        }
    }
    fn sclass_name(&self, g: usize) -> &'static str {
        match &self.sclass {
            Some(c) => STRUCT_CLASSES[c[g] as usize],
            None => "",
        }
    }

    /// `(weighted fraction, unweighted fraction, matched weight, total weight, matched n, n)`.
    fn frac(&self, side: usize, f: impl Fn(usize) -> bool) -> (f64, f64, f64, f64, usize, usize) {
        let (mut wn, mut wd, mut nm, mut n) = (0.0, 0.0, 0, 0);
        for g in self.range(side) {
            if !f(g) {
                continue;
            }
            let c = self.model(g).count;
            let m = self.matched(g);
            wd += c;
            n += 1;
            if m {
                wn += c;
                nm += 1;
            }
        }
        (pct(wn, wd), pct(nm as f64, n as f64), wn, wd, nm, n)
    }
}

// ---------------------------------------------------------------- sections

pub fn overview(c: &Ctx, rep: &mut Report) {
    let mut t = Table::new("overview", "Inputs", &["", "ref", "test"]);
    let row = |name: &str, f: &dyn Fn(usize) -> String| vec![name.to_string(), f(0), f(1)];
    t.row(row("label", &|s| c.side(s).label.clone()));
    t.row(row("transcript models", &|s| c.side(s).models.len().to_string()));
    t.row(row("  multi-exon", &|s| c.range(s).filter(|&g| c.model(g).exons.len() > 1).count().to_string()));
    t.row(row("  mono-exon", &|s| c.range(s).filter(|&g| c.model(g).is_mono()).count().to_string()));
    t.row(row("  matched (structure)", &|s| c.range(s).filter(|&g| c.matched(g)).count().to_string()));
    t.row(row("  unmatched", &|s| c.range(s).filter(|&g| !c.matched(g)).count().to_string()));
    t.row(row("assigned reads (count rows)", &|s| cnt(c.side(s).stats.assigned)));
    t.row(row("__ambiguous", &|s| cnt(c.side(s).stats.ambiguous)));
    t.row(row("__no_feature", &|s| cnt(c.side(s).stats.no_feature)));
    t.row(row("__not_aligned", &|s| cnt(c.side(s).stats.not_aligned)));
    t.row(row("grouped files", &|s| c.side(s).grouped.len().to_string()));
    t.row(row("per-read file", &|s| {
        c.side(s).reads.as_ref().map(|r| r.path.display().to_string()).unwrap_or_else(|| "-".into())
    }));
    rep.tables.push(t);
    let units = c.mt.comps.iter().filter(|u| u.matched()).count();
    for s in 0..2 {
        let un = c.range(s).filter(|&g| !c.matched(g)).count();
        rep.metrics.push(
            Metric::new(format!("unmatched_models_{}", SIDES[s]), "overview", "n", un as f64, un as f64)
                .desc("models (GTF transcripts) with no structural partner on the other side"),
        );
    }
    rep.metrics.push(Metric::new("matched_units", "overview", "n", units as f64, 0.0).desc(
        "connected components of the match graph containing models of both sides (alt-end models with one chain collapse into one unit)",
    ));
}

pub fn chain(c: &Ctx, rep: &mut Report) {
    let dir = ["ref_to_test", "test_to_ref"];
    for s in 0..2 {
        let (w, u, ..) = c.frac(s, |g| c.gated(g));
        let note = if c.exclude_nongene { " Models at non-`gene` loci excluded." } else { "" };
        rep.metrics.push(Metric::new(format!("chain_match_{}", dir[s]), "chain", "frac", w, 1.0 - w).thr(0.01).desc(format!(
            "sum of counts of {} models whose structure (multi-exon: chr+strand+exact intron chain; mono-exon: same strand, reciprocal overlap >= 0.8) exists on the other side / sum of all {} counts.{note}",
            SIDES[s], SIDES[s]
        )));
        rep.metrics.push(
            Metric::new(format!("chain_match_{}_models", dir[s]), "chain", "frac", u, 1.0 - u)
                .desc("same, unweighted (fraction of models)"),
        );
        let (wm, ..) = c.frac(s, |g| c.gated(g) && c.model(g).exons.len() > 1);
        let (wo, ..) = c.frac(s, |g| c.gated(g) && c.model(g).is_mono());
        rep.metrics.push(
            Metric::new(format!("chain_match_{}_multi", dir[s]), "chain", "frac", wm, 1.0 - wm)
                .desc("read-weighted, multi-exon models only"),
        );
        rep.metrics.push(
            Metric::new(format!("chain_match_{}_mono", dir[s]), "chain", "frac", wo, 1.0 - wo)
                .desc("read-weighted, mono-exon models only"),
        );
        if c.exclude_nongene && c.cat.is_some() {
            let (wa, ..) = c.frac(s, |_| true);
            rep.metrics.push(
                Metric::new(format!("chain_match_{}_all_loci", dir[s]), "chain", "frac", wa, 1.0 - wa)
                    .desc("read-weighted, including non-`gene` loci (not gated)"),
            );
        }
    }
}

/// Matched units admitted to count metrics.
fn unit_gated(c: &Ctx, u: &crate::matching::Comp) -> bool {
    u.members.iter().any(|&g| c.gated(g as usize))
}

pub fn counts(c: &Ctx, rep: &mut Report) {
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for u in &c.mt.comps {
        if u.matched() && unit_gated(c, u) {
            x.push(u.ref_count);
            y.push(u.test_count);
        }
    }
    let sp = spearman(&x, &y);
    let lx: Vec<f64> = x.iter().map(|v| v.ln_1p()).collect();
    let ly: Vec<f64> = y.iter().map(|v| v.ln_1p()).collect();
    let pl = pearson(&lx, &ly);
    let (mn, mx) = x.iter().zip(&y).fold((0.0, 0.0), |(a, b), (p, q)| (a + p.min(*q), b + p.max(*q)));
    let conc = pct(mn, mx);
    rep.metrics.push(Metric::new("spearman_matched", "counts", "corr", sp, 1.0 - sp).thr(0.005).desc(format!(
        "Spearman (average ranks) of (sum ref count, sum test count) over the {} matched units",
        x.len()
    )));
    rep.metrics.push(
        Metric::new("pearson_log1p_matched", "counts", "corr", pl, 1.0 - pl)
            .desc("Pearson of log1p counts over matched units"),
    );
    rep.metrics.push(
        Metric::new("count_concordance_matched", "counts", "frac", conc, 1.0 - conc)
            .desc("sum_u min(ref,test) / sum_u max(ref,test) over matched units"),
    );
}

pub fn fates(c: &Ctx, rep: &mut Report) {
    // count-table fates (always available)
    let both_reads = c.r.reads.is_some() && c.t.reads.is_some();
    let tab = |s: &Side| {
        let st = &s.stats;
        let tot = st.assigned + st.ambiguous + st.no_feature + st.not_aligned;
        [st.assigned / tot, st.ambiguous / tot, st.no_feature / tot, st.not_aligned / tot]
    };
    let (fr, ft) = (tab(c.r), tab(c.t));
    for (k, name) in ["assigned", "ambiguous", "no_feature", "not_aligned"].iter().enumerate() {
        let d = (ft[k] - fr[k]) * 100.0;
        let mut m = Metric::new(format!("fate_table_{name}"), "fates", "pp", d, d.abs()).sides(fr[k], ft[k]).desc(
            "fraction of the counts table: assigned (sum of feature rows), __ambiguous, __no_feature, __not_aligned over their total; test - ref in percentage points",
        );
        if !both_reads {
            m = m.thr(0.5);
        }
        rep.metrics.push(m);
    }
    // per-read files
    let mut ft_tab =
        Table::new("fates", "Per-read fine fates (reads)", &["fate", "ref", "test", "ref frac", "test frac"]);
    let mut fine = [[0u64; 10]; 2];
    let mut ntot = [0u64; 2];
    for s in 0..2 {
        if let Some(r) = &c.side(s).reads {
            for rec in &r.recs {
                fine[s][rec.fine() as usize] += 1;
            }
            ntot[s] = r.recs.len() as u64;
        }
    }
    if c.r.reads.is_none() && c.t.reads.is_none() {
        return;
    }
    for (k, f) in FINE_FATES.iter().enumerate() {
        if fine[0][k] + fine[1][k] == 0 {
            continue;
        }
        ft_tab.row(vec![
            f.to_string(),
            fine[0][k].to_string(),
            fine[1][k].to_string(),
            f4(pct(fine[0][k] as f64, ntot[0] as f64)),
            f4(pct(fine[1][k] as f64, ntot[1] as f64)),
        ]);
    }
    rep.tables.push(ft_tab);
    if !both_reads {
        return;
    }
    let mut coarse = [[0u64; 5]; 2];
    for s in 0..2 {
        for k in 0..10 {
            coarse[s][coarse_of(k as u8)] += fine[s][k];
        }
    }
    for (k, name) in COARSE_FATES.iter().enumerate() {
        let (a, b) = (pct(coarse[0][k] as f64, ntot[0] as f64), pct(coarse[1][k] as f64, ntot[1] as f64));
        let d = (b - a) * 100.0;
        let mut m = Metric::new(format!("fate_reads_{name}"), "fates", "pp", d, d.abs()).sides(a, b).desc(
            "fraction of reads in the per-read file by final fate: unique (unique, unique_minor_difference), ambiguous, inconsistent (inconsistent*), no_feature (noninformative, intergenic, discarded); test - ref in pp",
        );
        if k < 4 {
            m = m.thr(0.5);
        }
        rep.metrics.push(m);
    }
    // confusion matrix joined on read id
    let (rr, tr) = (c.r.reads.as_ref().unwrap(), c.t.reads.as_ref().unwrap());
    let mut conf = [[0u64; 6]; 6];
    let (mut uu, mut uu_same) = (0u64, 0u64);
    let (ra, ta) = (&rr.recs, &tr.recs);
    let (mut i, mut j) = (0, 0);
    while i < ra.len() || j < ta.len() {
        let ord = match (ra.get(i), ta.get(j)) {
            (Some(a), Some(b)) => a.hash().cmp(&b.hash()),
            (Some(_), None) => std::cmp::Ordering::Less,
            _ => std::cmp::Ordering::Greater,
        };
        match ord {
            std::cmp::Ordering::Less => {
                conf[coarse_of(ra[i].fine())][5] += 1;
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                conf[5][coarse_of(ta[j].fine())] += 1;
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                let (a, b) = (ra[i], ta[j]);
                let (ca, cb) = (coarse_of(a.fine()), coarse_of(b.fine()));
                conf[ca][cb] += 1;
                if ca == 0 && cb == 0 && a.model() != NO_MODEL && b.model() != NO_MODEL {
                    uu += 1;
                    if c.mt.comp_of[a.model() as usize] == c.mt.comp_of[c.nr() + b.model() as usize] {
                        uu_same += 1;
                    }
                }
                i += 1;
                j += 1;
            }
        }
    }
    let names: Vec<&str> = COARSE_FATES.iter().copied().chain(["absent"]).collect();
    let mut hdr = vec!["ref \\ test"];
    hdr.extend(names.iter());
    let mut t = Table::new("fates", "Per-read fate confusion matrix (joined on read_id)", &hdr);
    for (i, n) in names.iter().enumerate() {
        let mut row = vec![n.to_string()];
        row.extend(conf[i].iter().map(|v| v.to_string()));
        t.row(row);
    }
    rep.tables.push(t);
    let diag: u64 = (0..6).map(|i| conf[i][i]).sum();
    let all: u64 = conf.iter().flatten().sum();
    let agree = pct(diag as f64, all as f64);
    rep.metrics.push(
        Metric::new("fate_read_agreement", "fates", "frac", agree, 1.0 - agree)
            .desc("fraction of reads (union of both files) with the same coarse fate on both sides"),
    );
    let ma = pct(uu_same as f64, uu as f64);
    rep.metrics.push(
        Metric::new("read_model_agreement", "fates", "frac", ma, 1.0 - ma).desc(
            "among reads unique on both sides, fraction assigned to structurally matched models (same match unit)",
        ),
    );
    let mut t = Table::new("fates", "Per-read file notes", &["", "ref", "test"]);
    t.row(vec!["rows".into(), rr.n_rows.to_string(), tr.n_rows.to_string()]);
    t.row(vec!["reads".into(), rr.recs.len().to_string(), tr.recs.len().to_string()]);
    t.row(vec!["reads at >1 locus (best fate kept)".into(), rr.multi_locus.to_string(), tr.multi_locus.to_string()]);
    rep.tables.push(t);
    if rr.unknown_models + tr.unknown_models > 0 {
        rep.warnings.push(format!(
            "per-read files name models absent from the GTF: ref {}, test {}",
            rr.unknown_models, tr.unknown_models
        ));
    }
}

fn hist_bin(n: usize) -> usize {
    match n {
        0 | 1 => 0,
        2 => 1,
        3..=5 => 2,
        _ => 3,
    }
}
const HIST: [&str; 4] = ["1", "2", "3-5", "6+"];

/// Models per gene ID, and per gene-agnostic locus (same-strand clusters of overlapping spans).
pub fn tpg(c: &Ctx, rep: &mut Report) {
    let mut res = [[0.0f64; 2]; 2]; // [gene|locus][side]
    let mut hist = [[[0usize; 4]; 2]; 2];
    let mut ngroups = [[0usize; 2]; 2];
    for s in 0..2 {
        let mut by_gene: FxHashMap<&str, usize> = FxHashMap::default();
        let mut spans: Vec<(u32, u8, i32, i32)> = Vec::new();
        for g in c.range(s) {
            if !c.gated(g) {
                continue;
            }
            let m = c.model(g);
            *by_gene.entry(m.gene.as_str()).or_default() += 1;
            if !m.exons.is_empty() {
                let (a, b) = m.span();
                spans.push((m.chr, m.strand, a, b));
            }
        }
        for &n in by_gene.values() {
            hist[0][s][hist_bin(n)] += 1;
        }
        ngroups[0][s] = by_gene.len();
        res[0][s] = pct(by_gene.values().sum::<usize>() as f64, by_gene.len() as f64);
        spans.sort_unstable();
        let mut loci: Vec<usize> = Vec::new();
        let mut cur: Option<(u32, u8, i32)> = None;
        for &(chr, st, a, b) in &spans {
            match cur {
                Some((cc, cs, ce)) if cc == chr && cs == st && a <= ce => {
                    *loci.last_mut().unwrap() += 1;
                    cur = Some((cc, cs, ce.max(b)));
                }
                _ => {
                    loci.push(1);
                    cur = Some((chr, st, b));
                }
            }
        }
        for &n in &loci {
            hist[1][s][hist_bin(n)] += 1;
        }
        ngroups[1][s] = loci.len();
        res[1][s] = pct(spans.len() as f64, loci.len() as f64);
    }
    let names = ["tx_per_gene_mean", "tx_per_locus_mean"];
    let descs = [
        "mean models per gene_id in the GTF; value = relative difference (test - ref) / ref",
        "mean models per gene-agnostic locus (same-strand clusters of overlapping model spans; factors out gene-joiner differences); relative difference",
    ];
    for k in 0..2 {
        let d = pct(res[k][1] - res[k][0], res[k][0]);
        rep.metrics
            .push(Metric::new(names[k], "tpg", "rel", d, d.abs()).sides(res[k][0], res[k][1]).thr(0.01).desc(descs[k]));
    }
    let mut t = Table::new("tpg", "Models per gene / locus", &["grouping", "side", "groups", "1", "2", "3-5", "6+"]);
    for (k, kn) in ["gene_id", "locus"].iter().enumerate() {
        for s in 0..2 {
            let mut row = vec![kn.to_string(), SIDES[s].into(), ngroups[k][s].to_string()];
            row.extend((0..4).map(|b| hist[k][s][b].to_string()));
            t.row(row);
        }
    }
    let _ = HIST;
    rep.tables.push(t);
}

pub fn novelty(c: &Ctx, rep: &mut Report) {
    let label: Vec<u8> = (0..c.n()).map(|g| crate::annot::label_class(c.model(g))).collect();
    let mut views: Vec<(&str, &[&str], &[u8], bool)> = Vec::new();
    let has_s = c.sclass.is_some();
    if let Some(sc) = &c.sclass {
        views.push(("struct", &STRUCT_CLASSES, sc.as_slice(), true));
    }
    views.push(("label", &LABEL_CLASSES, label.as_slice(), !has_s));
    let mut t = Table::new(
        "novelty",
        "Model classes",
        &["view", "class", "ref models", "test models", "ref reads", "test reads"],
    );
    for (vname, names, cls, gate) in views {
        let mut nm = [[0usize; 4]; 2];
        let mut nw = [[0.0f64; 4]; 2];
        for s in 0..2 {
            for g in c.range(s) {
                if c.gated(g) {
                    nm[s][cls[g] as usize] += 1;
                    nw[s][cls[g] as usize] += c.model(g).count;
                }
            }
        }
        let tm = [nm[0].iter().sum::<usize>() as f64, nm[1].iter().sum::<usize>() as f64];
        let tw = [nw[0].iter().sum::<f64>(), nw[1].iter().sum::<f64>()];
        for k in 0..4 {
            if nm[0][k] + nm[1][k] == 0 {
                continue;
            }
            t.row(vec![
                vname.into(),
                names[k].into(),
                nm[0][k].to_string(),
                nm[1][k].to_string(),
                cnt(nw[0][k]),
                cnt(nw[1][k]),
            ]);
            for (what, a, b) in [
                ("models", pct(nm[0][k] as f64, tm[0]), pct(nm[1][k] as f64, tm[1])),
                ("reads", pct(nw[0][k], tw[0]), pct(nw[1][k], tw[1])),
            ] {
                let d = (b - a) * 100.0;
                let mut m = Metric::new(format!("{vname}_{}_{what}", names[k]), "novelty", "pp", d, d.abs()).sides(a, b).desc(if vname == "struct" {
                    "fraction of models / of model reads in each structural class vs the annotation (known: chain equals an annotated chain or mono-exon recip. overlap >= 0.8; nic: all introns or all splice sites annotated; nnic: otherwise); test - ref in pp"
                } else {
                    "fraction of models / reads per ID-label class (known: annotation source; nic/nnic: ID suffix); test - ref in pp"
                });
                if gate {
                    m = m.thr(1.0);
                }
                rep.metrics.push(m);
            }
        }
    }
    rep.tables.push(t);
}

pub fn grouped(c: &Ctx, rep: &mut Report) {
    for gr in &c.r.grouped {
        let Some(gt) = c.t.grouped.iter().find(|g| g.key == gr.key) else {
            rep.warnings.push(format!("grouped file {} has no counterpart on the test side", gr.key));
            continue;
        };
        let short = gr
            .key
            .trim_start_matches("discovered_transcript_grouped_")
            .trim_end_matches("_counts.linear.tsv")
            .trim_end_matches("counts.linear.tsv");
        let short = if short.is_empty() { "grouped" } else { short.trim_end_matches('_') };
        // per-side, per-group model counts (global model idx)
        let mut names: Vec<String> = gr.groups.clone();
        for g in &gt.groups {
            if !names.contains(g) {
                names.push(g.clone());
            }
        }
        let gidx = |n: &str| names.iter().position(|x| x == n).unwrap();
        let ng = names.len();
        let mut per: Vec<FxHashMap<u32, f64>> = vec![FxHashMap::default(); ng];
        let mut tot = vec![[0.0f64; 2]; ng];
        let mut unknown = 0usize;
        let mut sum_model = vec![0.0f64; c.n()];
        for (s, gf) in [(0usize, gr), (1usize, gt)] {
            let off = if s == 0 { 0 } else { c.nr() as u32 };
            for &(m, gi, v) in &gf.rows {
                if m == u32::MAX {
                    unknown += 1;
                    continue;
                }
                let g = gidx(&gf.groups[gi as usize]);
                let glob = m + off;
                *per[g].entry(glob).or_default() += v;
                sum_model[glob as usize] += v;
                if c.gated(glob as usize) {
                    tot[g][s] += v;
                }
            }
        }
        if unknown > 0 {
            rep.warnings.push(format!("{}: {unknown} grouped rows name features absent from the GTF", gr.key));
        }
        let bad = (0..c.n()).filter(|&g| (sum_model[g] - c.model(g).count).abs() > 0.01).count();
        if bad > 0 {
            rep.warnings
                .push(format!("{}: {bad} models whose grouped counts do not sum to the ungrouped count", gr.key));
        }
        let mut t = Table::new(
            "grouped",
            format!("Grouped counts: {}", gr.key),
            &[
                "group",
                "ref reads",
                "test reads",
                "rel diff",
                "match ref->test",
                "match test->ref",
                "spearman",
                "units",
            ],
        );
        for (g, name) in names.iter().enumerate() {
            let mut wsum = [[0.0f64; 2]; 2]; // [side][matched?]
            let mut units: FxHashMap<u32, [f64; 2]> = FxHashMap::default();
            for (&m, &v) in &per[g] {
                let gm = m as usize;
                if !c.gated(gm) {
                    continue;
                }
                let s = usize::from(gm >= c.nr());
                wsum[s][usize::from(c.matched(gm))] += v;
                if c.matched(gm) {
                    units.entry(c.mt.comp_of[gm]).or_default()[s] += v;
                }
            }
            // include matched units with zero count in this group on both sides? no: units touched by the group only
            let (x, y): (Vec<f64>, Vec<f64>) = units.values().map(|a| (a[0], a[1])).unzip();
            let sp = spearman(&x, &y);
            let fr = pct(wsum[0][1], wsum[0][0] + wsum[0][1]);
            let ft = pct(wsum[1][1], wsum[1][0] + wsum[1][1]);
            let rd = pct(tot[g][1] - tot[g][0], tot[g][0]);
            t.row(vec![
                name.clone(),
                cnt(tot[g][0]),
                cnt(tot[g][1]),
                f4(rd),
                f4(fr),
                f4(ft),
                f4(sp),
                x.len().to_string(),
            ]);
            let pre = format!("grouped_{short}_{name}");
            rep.metrics.push(
                Metric::new(format!("{pre}_spearman"), "grouped", "corr", sp, 1.0 - sp).thr(0.005).desc(
                    "per group: Spearman over matched units touched by the group (grouped counts summed per unit)",
                ),
            );
            rep.metrics.push(
                Metric::new(format!("{pre}_chain_match_ref_to_test"), "grouped", "frac", fr, 1.0 - fr)
                    .thr(0.01)
                    .desc("per group: read-weighted chain match using grouped counts"),
            );
            rep.metrics
                .push(Metric::new(format!("{pre}_chain_match_test_to_ref"), "grouped", "frac", ft, 1.0 - ft).thr(0.01));
            rep.metrics.push(
                Metric::new(format!("{pre}_total"), "grouped", "rel", rd, rd.abs())
                    .sides(tot[g][0], tot[g][1])
                    .desc("per group: relative difference of total grouped counts"),
            );
        }
        rep.tables.push(t);
    }
    for gt in &c.t.grouped {
        if !c.r.grouped.iter().any(|g| g.key == gt.key) {
            rep.warnings.push(format!("grouped file {} has no counterpart on the ref side", gt.key));
        }
    }
}

pub fn breakdowns(c: &Ctx, rep: &mut Report) {
    let row_for = |name: String, f: &dyn Fn(usize) -> bool| -> Vec<String> {
        let a = c.frac(0, f);
        let b = c.frac(1, f);
        vec![name, a.5.to_string(), b.5.to_string(), cnt(a.3), cnt(b.3), f4(a.0), f4(b.0)]
    };
    let hdr = ["", "ref models", "test models", "ref reads", "test reads", "match ref->test", "match test->ref"];
    if let Some(cat) = &c.cat {
        let mut t = Table::new("breakdown", "By locus category (GFF3 gene type; * = not `gene`/intergenic)", &hdr);
        for (k, n) in c.cats.iter().enumerate() {
            let label = if n != "gene" && n != "intergenic" { format!("{n} *") } else { n.clone() };
            t.row(row_for(label, &|g| cat[g] as usize == k));
        }
        rep.tables.push(t);
    }
    let mut t = Table::new("breakdown", "By chromosome", &hdr);
    for (k, n) in c.chroms.iter().enumerate() {
        t.row(row_for(n.clone(), &|g| c.gated(g) && c.model(g).chr as usize == k));
    }
    t.md_limit = 30;
    rep.tables.push(t);
    let bins = [
        (0.0, 0.0, "0"),
        (1.0, 1.0, "1"),
        (2.0, 4.0, "2-4"),
        (5.0, 19.0, "5-19"),
        (20.0, 99.0, "20-99"),
        (100.0, f64::MAX, "100+"),
    ];
    let mut t = Table::new("breakdown", "By unit depth (max of ref/test reads in the model's match unit)", &hdr);
    for (lo, hi, n) in bins {
        t.row(row_for(n.into(), &|g| {
            let u = &c.mt.comps[c.mt.comp_of[g] as usize];
            let d = u.ref_count.max(u.test_count);
            c.gated(g) && d >= lo && d <= hi
        }));
    }
    rep.tables.push(t);
}

// ---------------------------------------------------------------- residuals

/// Per-side spatial index for nearest-model explanations.
struct SpanIdx {
    v: Vec<(u32, u8, i32, u32)>,
    maxlen: i32,
}
impl SpanIdx {
    fn new(c: &Ctx, side: usize) -> SpanIdx {
        let mut v = Vec::new();
        let mut maxlen = 0;
        for g in c.range(side) {
            let m = c.model(g);
            if m.exons.is_empty() {
                continue;
            }
            let (a, b) = m.span();
            maxlen = maxlen.max(b - a + 1);
            v.push((m.chr, m.strand, a, g as u32));
        }
        v.sort_unstable();
        SpanIdx { v, maxlen }
    }
}

fn shared_introns(a: &Model, b: &Model) -> usize {
    let (x, y): (Vec<_>, Vec<_>) = (a.introns().collect(), b.introns().collect());
    let (mut i, mut j, mut k) = (0, 0, 0);
    while i < x.len() && j < y.len() {
        match x[i].cmp(&y[j]) {
            std::cmp::Ordering::Equal => {
                k += 1;
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
        }
    }
    k
}

fn exon_overlap(a: &Model, b: &Model) -> i64 {
    let (mut i, mut j, mut o) = (0, 0, 0i64);
    while i < a.exons.len() && j < b.exons.len() {
        let (x, y) = (a.exons[i], b.exons[j]);
        let ov = (x.1.min(y.1) - x.0.max(y.0) + 1) as i64;
        if ov > 0 {
            o += ov;
        }
        if x.1 < y.1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    o
}

/// Nearest model on the other side and a short description of how it differs.
fn explain(c: &Ctx, idx: &SpanIdx, g: usize) -> (String, String) {
    let m = c.model(g);
    if m.exons.is_empty() {
        return ("".into(), "no exons".into());
    }
    let (a, b) = m.span();
    let lo = (m.chr, m.strand, a - idx.maxlen, 0u32);
    let i0 = idx.v.partition_point(|x| *x < lo);
    let mut best: Option<(usize, i64, usize)> = None;
    for &(chr, st, s, o) in &idx.v[i0..] {
        if chr != m.chr || st != m.strand || s > b {
            break;
        }
        let om = c.model(o as usize);
        let ov = exon_overlap(m, om);
        if ov <= 0 {
            continue;
        }
        let k = shared_introns(m, om);
        if best.is_none_or(|(bk, bo, _)| (k, ov) > (bk, bo)) {
            best = Some((k, ov, o as usize));
        }
    }
    let Some((k, _, o)) = best else {
        return ("".into(), "no overlapping model on the other side (same strand)".into());
    };
    let om = c.model(o);
    let near = format!("{} ({})", om.id, cnt(om.count));
    let (na, nb) = (m.exons.len() - 1, om.exons.len() - 1);
    let why = match (na, nb) {
        (0, 0) => format!("mono-exon; best reciprocal overlap {:.2} < 0.8", recip(m.exons[0], om.exons[0])),
        (0, _) => format!("mono-exon overlapping a {}-intron model", nb),
        (_, 0) => "multi-exon; other side only has an overlapping mono-exon model".into(),
        _ if k == na && na < nb => format!("sub-chain: other model has {} extra intron(s)", nb - na),
        _ if k == nb && nb < na => format!("super-chain: {} extra intron(s) vs other model", na - nb),
        _ if na == nb => {
            let x: Vec<_> = m.introns().collect();
            let y: Vec<_> = om.introns().collect();
            let d = x.iter().zip(&y).map(|(p, q)| (p.0 - q.0).abs().max((p.1 - q.1).abs())).max().unwrap_or(0);
            format!("same intron count ({na}); {} intron(s) differ, max splice-site shift {d} bp", na - k)
        }
        _ => format!("shares {k}/{na} introns (other model has {nb})"),
    };
    (near, why)
}

/// Builds the residual table (top N in markdown) and returns `mismatches.tsv` content.
pub fn residuals(c: &Ctx, rep: &mut Report) -> String {
    let idx = [SpanIdx::new(c, 0), SpanIdx::new(c, 1)];
    let chroms = c.chroms;
    let locus = |m: &Model| {
        let (a, b) = m.span();
        format!("{}:{}-{}{}", chroms[m.chr as usize], a, b, m.strand as char)
    };
    struct Row {
        weight: f64,
        cols: Vec<String>,
    }
    let mut rows: Vec<Row> = Vec::new();
    for g in 0..c.n() {
        if c.matched(g) {
            continue;
        }
        let s = usize::from(g >= c.nr());
        let m = c.model(g);
        let (near, why) = explain(c, &idx[1 - s], g);
        let (rc, tc) = if s == 0 { (m.count, 0.0) } else { (0.0, m.count) };
        rows.push(Row {
            weight: m.count,
            cols: vec![
                format!("unmatched_{}", SIDES[s]),
                m.id.clone(),
                locus(m),
                m.exons.len().to_string(),
                cnt(rc),
                cnt(tc),
                c.cat_name(g).to_string(),
                c.sclass_name(g).to_string(),
                near,
                why,
            ],
        });
    }
    for u in &c.mt.comps {
        if !u.matched() || u.ref_count == u.test_count {
            continue;
        }
        let ids = |side: usize| {
            let v: Vec<&str> = u
                .members
                .iter()
                .filter(|&&g| usize::from(g as usize >= c.nr()) == side)
                .map(|&g| c.model(g as usize).id.as_str())
                .collect();
            if v.len() > 3 {
                format!("{},... ({} models)", v[..3].join(","), v.len())
            } else {
                v.join(",")
            }
        };
        let g0 = u.members[0] as usize;
        let lr = ((u.test_count + 1.0) / (u.ref_count + 1.0)).log2();
        rows.push(Row {
            weight: (u.test_count - u.ref_count).abs(),
            cols: vec![
                "discordant".into(),
                format!("{} | {}", ids(0), ids(1)),
                locus(c.model(g0)),
                c.model(g0).exons.len().to_string(),
                cnt(u.ref_count),
                cnt(u.test_count),
                c.cat_name(g0).to_string(),
                c.sclass_name(g0).to_string(),
                String::new(),
                format!("log2((test+1)/(ref+1)) = {lr:+.2}"),
            ],
        });
    }
    rows.sort_by(|a, b| b.weight.total_cmp(&a.weight).then_with(|| a.cols[1].cmp(&b.cols[1])));
    let hdr = [
        "kind",
        "model(s)",
        "locus",
        "exons",
        "ref count",
        "test count",
        "category",
        "struct class",
        "nearest other-side model",
        "explanation",
    ];
    let mut t = Table::new("residuals", "Largest unmatched / discordant models (by read weight)", &hdr);
    t.md_limit = c.top;
    let mut tsv = hdr.join("\t") + "\n";
    for r in &rows {
        tsv += &r.cols.join("\t");
        tsv.push('\n');
    }
    t.rows = rows.into_iter().map(|r| r.cols).collect();
    rep.tables.push(t);
    tsv
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn average_ranks() {
        assert_eq!(ranks(&[10.0, 20.0, 20.0, 5.0]), vec![2.0, 3.5, 3.5, 1.0]);
    }

    #[test]
    fn correlations() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(spearman(&x, &[2.0, 4.0, 6.0, 8.0, 100.0]), 1.0);
        assert!((spearman(&x, &[5.0, 4.0, 3.0, 2.0, 1.0]) + 1.0).abs() < 1e-12);
        // scipy.stats.spearmanr([1,2,3,4,5],[5,6,7,8,7]) = 0.8207826816681233
        assert!((spearman(&x, &[5.0, 6.0, 7.0, 8.0, 7.0]) - 0.8207826816681233).abs() < 1e-12);
        assert_eq!(pearson(&[3.0, 3.0], &[3.0, 3.0]), 1.0);
        assert!(pearson(&[3.0, 3.0], &[1.0, 2.0]).is_nan());
    }
}
