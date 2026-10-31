//! U10b: `TerminalCounter.add_read_info/flush/_predict_rows` and the PolyA/TSS subclasses
//! (`terminal_prediction/terminal_counter.py:67-312, 497-552`), ungrouped, per block.
//!
//! 4.1 feeds every read of a block to the ungrouped polyA (always) and TSS (`--fl_data` only)
//! counters, then `flush()` predicts on that block's accumulation (`last_gene_predictions`,
//! consumed by the intron graph, U11) and clears it (`parallel_workers.py:312-330`). The
//! chromosome-wide `dump()` TSV is a side output and is not ported (design §1.3).
//!
//! **Accumulation.** One flat `(slot, pos)` buffer per counter in read order plus a dense
//! `TxIdx → slot` map (slot = first-read order of the transcript = 4.1 dict order). No per-read
//! allocation; the buffers are reused across blocks. Read order is kept (not a histogram)
//! because `np.var`/`stats.skew` run numpy's pairwise sum over `data` in add order, so the
//! float bits of the `var`/`skew` features depend on it ([`read_order_moments`]).
//!
//! **Prediction.** `_predict_rows` is the same pipeline as U10a's `detect_peaks` except for the
//! data order of `var`/`skew`: histogram → `find_peaks`/`peak_widths` → rank by height (D-T1
//! ties by position) → XGBoost → accepted peaks `(tx, start + peak_location)`. Rows are in
//! transcript dict order, ranked order within a transcript; transcripts whose histogram has no
//! peak (unreachable: the zero padding guarantees one) come first, as 4.1 concatenates the
//! zero-peak frame before the real-peak frame.

use crate::peaks::{feature_row, histogram_record, np_sum, rank_peaks, HistRecord, Histogram};
use crate::xgb::{XgbModel, XgbModels};
use iq_types::{AssignmentType, BlockStore, IsoRef, IsoformIndex, Params, Pos, Predictions, ReadView, Strand, TxIdx};

/// Which counter (`PolyACounter` / `TSSCounter`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TerminalKind {
    PolyA,
    Tss,
}

impl TerminalKind {
    /// `_passes_filter` + `_extract_position` (`terminal_counter.py:507-516, 536-545`).
    #[inline]
    pub fn read_end(self, r: &ReadView<'_>) -> Option<Pos> {
        match self {
            TerminalKind::PolyA => {
                // `polyA_found and polya_info`: polya_info is always present on a record.
                if !r.polya_found() {
                    return None;
                }
                let info = r.polya();
                match r.strand() {
                    Strand::Plus if info.ext_a != -1 => Some(info.ext_a),
                    Strand::Minus if info.ext_t != -1 => Some(info.ext_t),
                    _ => None,
                }
            }
            TerminalKind::Tss => match r.strand() {
                // `corrected_exons` non-empty is checked by the caller.
                Strand::Plus => Some(r.start()),
                Strand::Minus => Some(r.end()),
                Strand::Unknown => None,
            },
        }
    }
}

/// `_ACCEPTED_ASSIGNMENT_TYPES` (`terminal_counter.py:47-52`).
#[inline]
fn accepted_type(t: AssignmentType) -> bool {
    matches!(
        t,
        AssignmentType::Unique
            | AssignmentType::UniqueMinorDifference
            | AssignmentType::Inconsistent
            | AssignmentType::InconsistentNonIntronic
    )
}

/// `TerminalCounter.transcripts` for one block: per-transcript read-end positions in add order.
#[derive(Clone, Debug, Default)]
pub struct TerminalAccumulator {
    /// `TxIdx → slot + 1` (0 = not seen in this block); only touched entries are reset.
    slot_of: Vec<u32>,
    /// Transcripts in first-read order.
    txs: Vec<TxIdx>,
    /// `(slot, pos)` in add order.
    ends: Vec<(u32, Pos)>,
    // ---- flush scratch ----
    offs: Vec<u32>,
    fill: Vec<u32>,
    data: Vec<Pos>,
    sq: Vec<f64>,
    hist: Histogram,
}

impl TerminalAccumulator {
    pub fn new() -> Self {
        Self::default()
    }

    /// No transcript accumulated.
    pub fn is_empty(&self) -> bool {
        self.txs.is_empty()
    }

    /// `add_read_info(read_assignment)` (`terminal_counter.py:186-219`), ungrouped.
    #[inline]
    pub fn add_read(&mut self, kind: TerminalKind, r: &ReadView<'_>) {
        if !accepted_type(r.atype()) || r.n_corrected_exons() == 0 {
            return;
        }
        let Some(m) = r.first_match() else { return };
        let Some(pos) = kind.read_end(r) else { return };
        match m.tx() {
            Some(IsoRef::Ref(tx)) => self.add(tx, pos),
            // Reference assignments never carry block-model ids; `Predictions` is keyed by
            // reference transcript.
            Some(IsoRef::Model(_)) | None => {}
        }
    }

    /// Appends one read end of transcript `tx`.
    #[inline]
    pub fn add(&mut self, tx: TxIdx, pos: Pos) {
        let t = tx as usize;
        if t >= self.slot_of.len() {
            self.slot_of.resize(t + 1, 0);
        }
        let mut s = self.slot_of[t];
        if s == 0 {
            self.txs.push(tx);
            s = self.txs.len() as u32;
            self.slot_of[t] = s;
        }
        self.ends.push((s - 1, pos));
    }

    /// Per-transcript positions in add order, transcripts in first-read order (the 4.1
    /// `self.transcripts` dict: `{tx: {'data': [...]}}`).
    pub fn transcripts(&self) -> Vec<(TxIdx, Vec<Pos>)> {
        let mut v: Vec<(TxIdx, Vec<Pos>)> = self.txs.iter().map(|&t| (t, Vec::new())).collect();
        for &(s, p) in &self.ends {
            v[s as usize].1.push(p);
        }
        v
    }

    /// Forgets the accumulation (keeps capacity).
    pub fn clear(&mut self) {
        for &t in &self.txs {
            self.slot_of[t as usize] = 0;
        }
        self.txs.clear();
        self.ends.clear();
    }

    /// `flush()` (`terminal_counter.py:221-243`) via `_predict_rows` (`:280-312`): appends
    /// `last_gene_predictions` to `out` and clears the accumulation.
    pub fn flush_into(&mut self, model: &XgbModel, out: &mut Vec<(TxIdx, Pos)>) {
        let n_tx = self.txs.len();
        if n_tx == 0 {
            return;
        }
        // Stable counting sort of `ends` by slot: data of slot s = data[offs[s]..offs[s+1]].
        self.offs.clear();
        self.offs.resize(n_tx + 1, 0);
        for &(s, _) in &self.ends {
            self.offs[s as usize + 1] += 1;
        }
        for s in 0..n_tx {
            self.offs[s + 1] += self.offs[s];
        }
        self.data.clear();
        self.data.resize(self.ends.len(), 0);
        self.fill.clear();
        self.fill.extend_from_slice(&self.offs[..n_tx]);
        for &(s, p) in &self.ends {
            let f = &mut self.fill[s as usize];
            self.data[*f as usize] = p;
            *f += 1;
        }
        let base = out.len();
        let mut n_zero = 0;
        for s in 0..n_tx {
            let tx = self.txs[s];
            let data = &self.data[self.offs[s] as usize..self.offs[s + 1] as usize];
            let rec = record_for(data, &mut self.hist, &mut self.sq);
            if rec.peaks.is_empty() {
                // Zero-peak fallback (unreachable, ported): the mode, always accepted; this
                // frame precedes the real-peak rows in 4.1's concat.
                out.insert(base + n_zero, (tx, (rec.start as i64 + rec.mode) as Pos));
                n_zero += 1;
                continue;
            }
            for (i, rel) in rank_peaks(&rec) {
                let p = &rec.peaks[i];
                if model.predict(&feature_row(&rec, p, rel)) {
                    out.push((tx, (rec.start as i64 + p.location) as Pos));
                }
            }
        }
        self.clear();
    }
}

/// `_build_peak_dataframe` record for one transcript's `data` (read order, non-empty).
pub(crate) fn record_for(data: &[Pos], hist: &mut Histogram, sq: &mut Vec<f64>) -> HistRecord {
    // Histogram counts are order-free; build them in position order.
    let mut sorted = data.to_vec();
    sorted.sort_unstable();
    hist.clear();
    for &p in &sorted {
        *hist.entry(p).or_insert(0) += 1;
    }
    let mut rec = histogram_record(hist).expect("non-empty transcript data");
    let (var, skew) = read_order_moments(data, sq);
    rec.var = var;
    rec.skew = skew;
    rec
}

/// `np.var(data)` and the counter's skew (`stats.skew(data)` if `len(data) > 3 and var >=
/// 1e-12`, else NaN) with `data` in add order — the same operations as U10a's
/// `histogram_record` (`peaks.rs`), which expands the histogram in dict order instead.
pub fn read_order_moments(data: &[Pos], sq: &mut Vec<f64>) -> (f64, f64) {
    let n = data.len();
    let nf = n as f64;
    // np.mean of an int64 array: float64 sum of integers (exact below 2**53) / n.
    let sum: i128 = data.iter().map(|&p| p as i128).sum();
    let mean = sum as f64 / nf;
    sq.clear();
    sq.extend(data.iter().map(|&p| {
        let d = p as f64 - mean;
        d * d
    }));
    let var = np_sum(sq) / nf;
    if !(n > 3 && var >= 1e-12) {
        return (var, f64::NAN);
    }
    // scipy: m2 = mean(d**2) (== var, same ops), m3 = mean((d*d)*d).
    for (v, &p) in sq.iter_mut().zip(data) {
        *v *= p as f64 - mean;
    }
    let m3 = np_sum(sq) / nf;
    let eps_mean = f64::EPSILON * mean;
    let skew = if var <= eps_mean * eps_mean { f64::NAN } else { m3 / var.powf(1.5) };
    (var, skew)
}

/// `fn block_predictions(&BlockStore, &IsoformIndex, &XgbModels) -> Predictions`: feeds every
/// read of the block (block order) to the ungrouped polyA and, when `params.use_tss_model`
/// (`--fl_data`), TSS counters and flushes them. `_index` is unused: the annotated position it
/// would give only sets the side-output `flag` column.
pub fn block_predictions(store: &BlockStore, _index: &IsoformIndex, models: &XgbModels, params: &Params) -> Predictions {
    let mut polya = TerminalAccumulator::new();
    let mut tss = TerminalAccumulator::new();
    let mut out = Predictions::default();
    block_predictions_with(store, models, params.use_tss_model, &mut polya, &mut tss, &mut out);
    out
}

/// [`block_predictions`] with reusable accumulators (cleared on return) and output.
pub fn block_predictions_with(
    store: &BlockStore,
    models: &XgbModels,
    use_tss: bool,
    polya: &mut TerminalAccumulator,
    tss: &mut TerminalAccumulator,
    out: &mut Predictions,
) {
    out.polya.clear();
    out.tss.clear();
    for r in store.iter() {
        polya.add_read(TerminalKind::PolyA, &r);
        if use_tss {
            tss.add_read(TerminalKind::Tss, &r);
        }
    }
    polya.flush_into(&models.polya, &mut out.polya);
    tss.flush_into(&models.tss, &mut out.tss);
}

pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U10bFlush>();
}

pub mod fixture {
    use super::*;
    use crate::peaks::{local_maxima, raw_peak, HISTOGRAM_PAD, PEAK_DISTANCE};
    use crate::xgb::bundled_models;
    use iq_types::fixture::schema::{PredictionsOut, U10bFlushIn};
    use iq_types::fixture::{Diff, DiffEntry, FixResult, FixtureCtx, UnitFixture, Verdict};
    use iq_types::DevId;
    use std::cell::RefCell;

    /// What the last `run` on this thread accumulated: `(kind, [(tx, data)])`. `compare` is
    /// called right after `run` on the same thread (`UnitFixture::judge`) and needs the
    /// histograms for the D-T1 check, which `PredictionsOut` does not carry.
    type LastRun = (TerminalKind, Vec<(String, Vec<Pos>)>);
    thread_local! {
        static LAST_RUN: RefCell<Option<LastRun>> = const { RefCell::new(None) };
    }

    pub struct U10bFlush;
    impl UnitFixture for U10bFlush {
        const UNIT: &'static str = "U10b";
        const FN: &'static str = "flush";
        const DEVIATIONS: &'static [DevId] = &[DevId::T1];
        type In = U10bFlushIn;
        type Out = PredictionsOut;

        fn run(input: &Self::In, ctx: &FixtureCtx) -> FixResult<Self::Out> {
            LAST_RUN.with(|l| *l.borrow_mut() = None);
            let block = ctx.block(&input.block)?;
            let index = ctx.index(&block.header.index)?;
            let params = ctx.params()?;
            let models = bundled_models();
            let kind = match input.kind.as_str() {
                "polya" => TerminalKind::PolyA,
                "tss" => TerminalKind::Tss,
                k => return Err(format!("unknown kind {k:?}").into()),
            };
            if kind == TerminalKind::Tss && !params.use_tss_model {
                return Err("TSS flush recorded without --fl_data".into());
            }
            // Accumulation check against the recorded `self.transcripts`.
            let mut acc = TerminalAccumulator::new();
            for r in block.store.iter() {
                acc.add_read(kind, &r);
            }
            let got_acc: Vec<(String, Vec<Pos>)> = {
                let ids = ctx.ids();
                acc.transcripts().into_iter().map(|(t, d)| (ids.iso_name(IsoRef::Ref(t)).to_owned(), d)).collect()
            };
            let want_acc: Vec<(String, Vec<Pos>)> =
                input.transcripts.iter().map(|(t, e)| (t.clone(), e.data.clone())).collect();
            if want_acc.is_empty() && !got_acc.is_empty() {
                // The second flush of the chromosome's last block: `dump()` calls `flush()`
                // again on the already-flushed (empty) accumulator (schema §5).
                return Ok(PredictionsOut { predictions: Vec::new() });
            }
            if got_acc != want_acc {
                return Err(format!(
                    "accumulation differs: got {} tx / {} reads, want {} tx / {} reads",
                    got_acc.len(),
                    got_acc.iter().map(|x| x.1.len()).sum::<usize>(),
                    want_acc.len(),
                    want_acc.iter().map(|x| x.1.len()).sum::<usize>()
                )
                .into());
            }
            let preds = block_predictions(&block.store, &index, models, &params);
            let list = match kind {
                TerminalKind::PolyA => preds.polya,
                TerminalKind::Tss => preds.tss,
            };
            LAST_RUN.with(|l| *l.borrow_mut() = Some((kind, got_acc)));
            let ids = ctx.ids();
            Ok(PredictionsOut {
                predictions: list.into_iter().map(|(t, p)| (ids.iso_name(IsoRef::Ref(t)).to_owned(), p)).collect(),
            })
        }

        fn compare(got: &Self::Out, want: &Self::Out, _ctx: &FixtureCtx) -> Verdict {
            if got == want {
                return Verdict::Exact;
            }
            let last = LAST_RUN.with(|l| l.borrow_mut().take());
            let res = match &last {
                Some((kind, txs)) => {
                    let model = match kind {
                        TerminalKind::PolyA => &bundled_models().polya,
                        TerminalKind::Tss => &bundled_models().tss,
                    };
                    tie_explained(&got.predictions, &want.predictions, txs, model)
                }
                None => Err("no accumulation recorded for the D-T1 check".to_owned()),
            };
            match res {
                Ok(()) => Verdict::Explained(DevId::T1),
                Err(msg) => Verdict::Mismatch(Diff {
                    entries: vec![
                        DiffEntry {
                            path: ".predictions".into(),
                            got: format!("{:?}", got.predictions),
                            want: format!("{:?}", want.predictions),
                        },
                        DiffEntry { path: "(D-T1 check)".into(), got: "inconsistent".into(), want: msg },
                    ],
                    truncated: 0,
                    error: None,
                }),
            }
        }
    }

    /// Splits a prediction list into per-transcript runs, in `txs` (dict) order. Errors if a
    /// transcript's rows are not contiguous, or runs are out of dict order / unknown.
    fn per_tx(preds: &[(String, Pos)], txs: &[(String, Vec<Pos>)]) -> Result<Vec<Vec<Pos>>, String> {
        let mut out = vec![Vec::new(); txs.len()];
        let mut cur = 0usize;
        let mut prev: Option<&str> = None;
        for (t, p) in preds {
            if prev != Some(t.as_str()) {
                let k = txs[cur..]
                    .iter()
                    .position(|(x, _)| x == t)
                    .ok_or_else(|| format!("transcript {t} out of dict order or unknown"))?;
                cur += k;
                if !out[cur].is_empty() {
                    return Err(format!("transcript {t} rows not contiguous"));
                }
                prev = Some(t.as_str());
            }
            out[cur].push(*p);
        }
        Ok(out)
    }

    /// D-T1: is `want` what 4.1 produces under some order of equal heights? Per transcript,
    /// either the rows are a permutation of iqrs's rows within equal heights (ranking ties), or
    /// some visiting order of equal-height maxima in the distance filter yields the same set of
    /// accepted peaks, listed height-descending. Peak acceptance does not depend on rank
    /// (rank is not a feature), only on the kept set.
    pub(crate) fn tie_explained(
        got: &[(String, Pos)],
        want: &[(String, Pos)],
        txs: &[(String, Vec<Pos>)],
        model: &XgbModel,
    ) -> Result<(), String> {
        let g = per_tx(got, txs)?;
        let w = per_tx(want, txs)?;
        let mut hist = Histogram::default();
        let mut sq = Vec::new();
        for (k, (tx, data)) in txs.iter().enumerate() {
            if g[k] == w[k] {
                continue;
            }
            let rec = record_for(data, &mut hist, &mut sq);
            let height = |pos: Pos| -> Option<u32> {
                let i = (pos as i64 - rec.start as i64) + HISTOGRAM_PAD as i64;
                rec.padded.get(usize::try_from(i).ok()?).copied()
            };
            let hs: Vec<u32> = w[k].iter().map(|&p| height(p).unwrap_or(0)).collect();
            if hs.windows(2).any(|h| h[1] > h[0]) {
                return Err(format!("{tx}: rows not height-descending"));
            }
            let mut want_set = w[k].clone();
            want_set.sort_unstable();
            if want_set.windows(2).any(|p| p[0] == p[1]) {
                return Err(format!("{tx}: duplicate rows"));
            }
            if !reachable_accepted_sets(&rec, model).contains(&want_set) {
                return Err(format!("{tx}: {:?} not produced by any tie order", w[k]));
            }
        }
        Ok(())
    }

    /// Distinct sorted accepted-position sets over distance-filter tie orders: exhaustive when
    /// the product of tie-group permutations is ≤ 5040, else 5040 seeded random orders.
    fn reachable_accepted_sets(rec: &HistRecord, model: &XgbModel) -> Vec<Vec<Pos>> {
        const CAP: u64 = 5040;
        let x = &rec.padded;
        let xf: Vec<f64> = x.iter().map(|&c| c as f64).collect();
        let cand = local_maxima(x);
        // Tie groups (indices into cand), height descending.
        let mut by_h: Vec<usize> = (0..cand.len()).collect();
        by_h.sort_by(|&a, &b| x[cand[b]].cmp(&x[cand[a]]).then(a.cmp(&b)));
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for i in by_h {
            match groups.last_mut() {
                Some(gr) if x[cand[gr[0]]] == x[cand[i]] => gr.push(i),
                _ => groups.push(vec![i]),
            }
        }
        let fact = |n: usize| (1..=n as u64).try_fold(1u64, |a, b| a.checked_mul(b)).unwrap_or(u64::MAX);
        let total = groups.iter().try_fold(1u64, |a, gr| a.checked_mul(fact(gr.len()))).unwrap_or(u64::MAX);
        let mut sets: Vec<Vec<Pos>> = Vec::new();
        let mut rng: u64 = 0x9e37_79b9_7f4a_7c15;
        let n_iter = total.min(CAP);
        for it in 0..n_iter {
            let mut order = Vec::with_capacity(cand.len());
            let mut code = it;
            for gr in &groups {
                let mut g = gr.clone();
                if total <= CAP {
                    // Mixed-radix (factorial number system) decode of `it`.
                    let mut perm = Vec::with_capacity(g.len());
                    while !g.is_empty() {
                        let r = (code % g.len() as u64) as usize;
                        code /= g.len() as u64;
                        perm.push(g.remove(r));
                    }
                    g = perm;
                } else {
                    for i in (1..g.len()).rev() {
                        rng ^= rng << 13;
                        rng ^= rng >> 7;
                        rng ^= rng << 17;
                        g.swap(i, (rng % (i as u64 + 1)) as usize);
                    }
                }
                order.extend(g);
            }
            let mut keep = vec![true; cand.len()];
            for &j in &order {
                if !keep[j] {
                    continue;
                }
                for (k, kp) in keep.iter_mut().enumerate() {
                    if k != j && cand[k].abs_diff(cand[j]) < PEAK_DISTANCE {
                        *kp = false;
                    }
                }
            }
            let peaks: Vec<_> = cand.iter().zip(&keep).filter(|(_, &k)| k).map(|(&p, _)| raw_peak(&xf, x, p)).collect();
            let r = HistRecord { peaks, ..rec.clone() };
            let mut acc: Vec<Pos> = rank_peaks(&r)
                .into_iter()
                .filter(|&(i, rel)| model.predict(&feature_row(&r, &r.peaks[i], rel)))
                .map(|(i, _)| (r.start as i64 + r.peaks[i].location) as Pos)
                .collect();
            acc.sort_unstable();
            if !sets.contains(&acc) {
                sets.push(acc);
            }
        }
        sets
    }
}

#[cfg(test)]
mod tests {
    use super::fixture::tie_explained;
    use super::*;
    use crate::peaks::{local_maxima, HISTOGRAM_PAD};
    use crate::xgb::bundled_models;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// Clustered read ends around a few sites, in a random (read) order.
    fn random_data(rng: &mut Rng) -> Vec<Pos> {
        let base = 1_000_000 + rng.below(100_000_000) as Pos;
        let n_sites = 1 + rng.below(4);
        let n = 1 + rng.below(300) as usize;
        let sites: Vec<Pos> = (0..n_sites).map(|_| base + rng.below(200) as Pos).collect();
        (0..n)
            .map(|_| sites[rng.below(n_sites) as usize] + rng.below(7) as Pos - 3)
            .collect()
    }

    #[test]
    fn moments_match_u10a_in_dict_order() {
        let mut rng = Rng(0x1234_5678_9abc_def1);
        let mut sq = Vec::new();
        for _ in 0..2000 {
            let data = random_data(&mut rng);
            let hist: Histogram = data.iter().fold(Histogram::default(), |mut h, &p| {
                *h.entry(p).or_insert(0) += 1;
                h
            });
            // U10a expands the histogram in dict order: positions grouped by first occurrence.
            let grouped: Vec<Pos> = hist.iter().flat_map(|(&p, &c)| std::iter::repeat_n(p, c as usize)).collect();
            let rec = histogram_record(&hist).unwrap();
            let (var, skew) = read_order_moments(&grouped, &mut sq);
            assert_eq!(var.to_bits(), rec.var.to_bits());
            assert!(skew.to_bits() == rec.skew.to_bits() || (skew.is_nan() && rec.skew.is_nan()));
        }
    }

    #[test]
    fn accumulator_order_and_reuse() {
        let mut acc = TerminalAccumulator::new();
        acc.add(7, 100);
        acc.add(3, 50);
        acc.add(7, 101);
        acc.add(1000, 5);
        assert_eq!(acc.transcripts(), vec![(7, vec![100, 101]), (3, vec![50]), (1000, vec![5])]);
        let mut out = Vec::new();
        acc.flush_into(&bundled_models().polya, &mut out);
        assert!(acc.is_empty());
        acc.add(3, 9);
        acc.add(7, 8);
        assert_eq!(acc.transcripts(), vec![(3, vec![9]), (7, vec![8])]);
    }

    /// Flush output = per transcript, accepted peaks of `detect_peaks` on its histogram, with
    /// var/skew from read order; equal to `detect_peaks` whenever those bits agree.
    #[test]
    fn flush_matches_detect_peaks() {
        let mut rng = Rng(0xdead_beef_cafe_f00d);
        let model = &bundled_models().polya;
        let (mut same, mut total, mut var_bits_differ) = (0, 0, 0);
        for _ in 0..500 {
            let data = random_data(&mut rng);
            let mut acc = TerminalAccumulator::new();
            for &p in &data {
                acc.add(1, p);
            }
            let mut out = Vec::new();
            acc.flush_into(model, &mut out);
            let hist: Histogram = data.iter().fold(Histogram::default(), |mut h, &p| {
                *h.entry(p).or_insert(0) += 1;
                h
            });
            let dp: Vec<(TxIdx, Pos)> = crate::detect_peaks(&hist, model).iter().map(|p| (1, p.position)).collect();
            total += 1;
            let rec = histogram_record(&hist).unwrap();
            if read_order_moments(&data, &mut Vec::new()).0.to_bits() != rec.var.to_bits() {
                var_bits_differ += 1;
            }
            if out == dp {
                same += 1;
            } else {
                // Only a var/skew bit difference may separate them.
                let rec = histogram_record(&hist).unwrap();
                let (v, s) = read_order_moments(&data, &mut Vec::new());
                assert!(v.to_bits() != rec.var.to_bits() || s.to_bits() != rec.skew.to_bits());
            }
        }
        eprintln!("flush == detect_peaks: {same}/{total}; var bits differ by data order: {var_bits_differ}");
        assert!(same * 10 >= total * 9, "{same}/{total}");
    }

    #[test]
    fn tie_checker_accepts_iqrs_and_rejects_shifts() {
        let mut rng = Rng(0x0bad_cafe_1234_5678);
        let mut rejected = 0;
        for it in 0..1600 {
            let model = if it % 2 == 0 { &bundled_models().polya } else { &bundled_models().tss };
            let data = random_data(&mut rng);
            let txs = vec![("t".to_owned(), data.clone())];
            let mut acc = TerminalAccumulator::new();
            for &p in &data {
                acc.add(0, p);
            }
            let mut out = Vec::new();
            acc.flush_into(model, &mut out);
            let got: Vec<(String, Pos)> = out.iter().map(|&(_, p)| ("t".to_owned(), p)).collect();
            assert!(tie_explained(&got, &got, &txs, model).is_ok());
            let Some(&(_, p0)) = got.first() else { continue };
            let rec = record_for(&data, &mut Histogram::default(), &mut Vec::new());
            let cand = local_maxima(&rec.padded);
            let shifted = p0 + 1;
            let idx = (shifted - rec.start) as usize + HISTOGRAM_PAD;
            if cand.contains(&idx) {
                continue;
            }
            let mut bad = got.clone();
            bad[0].1 = shifted;
            assert!(tie_explained(&got, &bad, &txs, model).is_err());
            // Without ties at all, a dropped row is a different (unreachable) accepted set.
            let hs: Vec<u32> = cand.iter().map(|&c| rec.padded[c]).collect();
            if hs.iter().enumerate().all(|(i, h)| hs[i + 1..].iter().all(|g| g != h)) {
                assert!(tie_explained(&got, &got[1..], &txs, model).is_err());
            }
            rejected += 1;
        }
        assert!(rejected > 100, "{rejected}");
    }
}
