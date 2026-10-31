"""Generates randomized U10a cases from the live IsoQuant 4.1 Python (numpy/scipy/xgboost).

Usage (from iqrs/):
    ../env/bin/python crates/iq_terminal/tests/data/gen_u10a.py [native|stable] [seed] [n_batch] [n_xgb]
Writes crates/iq_terminal/tests/data/u10a_random_<mode>.jsonl.gz in the fixture schema
(harness/iqcap/SCHEMA.md, unit U10a, fns detect_peaks_batch / xgb_predict), computed exactly as
harness/iqcap/hooks.py captures them. `cargo test -p iq_terminal` replays the file; Python is not
needed at test time.
"""
import gzip
import json
import math
import os
import sys

NOSIMD = "AVX512_SPR AVX512_ICL X86_V4 X86_V3"
if (sys.argv[1] if len(sys.argv) > 1 else "native") == "stable" and os.environ.get("NPY_DISABLE_CPU_FEATURES") != NOSIMD:
    # stable mode = numpy's non-SIMD argsort (introsort + insertion sort), which iqrs emulates
    # exactly (D-T1, `iq_terminal::peaks::np_argsort`); the dispatch is fixed at numpy import.
    os.execve(sys.executable, [sys.executable] + sys.argv, dict(os.environ, NPY_DISABLE_CPU_FEATURES=NOSIMD))

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../../../../IsoQuant"))
from isoquant_lib.terminal_prediction import terminal_peaks as TP  # noqa: E402

MODE = sys.argv[1] if len(sys.argv) > 1 else "native"
OUT = os.path.join(HERE, f"u10a_random_{MODE}.jsonl.gz")

if MODE == "stable":
    # No patching: numpy's own non-SIMD argsort (see the re-exec above) is the D-T1 order iqrs
    # emulates, so every case must be bit-exact.
    assert np.argsort(np.zeros(17)).tolist() != list(range(17)), "SIMD dispatch still active?"
elif MODE != "native":
    sys.exit(f"unknown mode {MODE}")


def fval(x):
    x = float(x)
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Infinity" if x > 0 else "-Infinity"
    return x


def dumps(o):
    return json.dumps(o, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


MODELS = {"polya": TP.get_polya_model(), "tss": TP.get_tss_model()}


def batch_case(histograms, model_name):
    """detect_peaks_batch in/out, as hooks.py w_dpb records them."""
    model = MODELS[model_name]
    hs = [None if not h else [[int(p), int(c)] for p, c in h.items()] for h in histograms]
    records, rows = [], []
    for i, counts in enumerate(histograms):
        if not counts:
            records.append(None)
            continue
        r = TP._histogram_record(dict(counts))
        rr = {"start": int(r["start"]), "mode": int(r["mode"]), "var": fval(r["var"]), "skew": fval(r["skew"]),
              "entropy": fval(r["entropy"]), "mean_height": fval(r["mean_height"]),
              "peak_count": int(r["peak_count"]), "histogram": [int(x) for x in r["histogram"]],
              "raw_peaks": None}
        if r["peak_count"]:
            rp = r["raw_peaks"]
            rr["raw_peaks"] = {"peak_location": [int(x) for x in rp["peak_location"]],
                               "peak_width": [fval(x) for x in rp["peak_width"]],
                               "peak_left": [int(x) for x in rp["peak_left"]],
                               "peak_right": [int(x) for x in rp["peak_right"]],
                               "peak_heights": [int(x) for x in rp["peak_heights"]]}
            for pk in TP._rank_peaks(r):
                rows.append((i, [r["var"], r["skew"], r["peak_count"], pk["peak_width"], r["entropy"],
                                 r["mean_height"], pk["peak_heights"], pk["relative_height"]]))
        records.append(rr)
    row_out = []
    if rows:
        feats = pd.DataFrame([x[1] for x in rows], columns=TP.FEATURE_COLUMNS).astype(float)
        lab = model.predict(feats)
        proba = model.predict_proba(feats)[:, 1]
        margin = np.asarray(model.predict(feats, output_margin=True)).ravel()
        for j, (i, _) in enumerate(rows):
            row_out.append({"hist": i, "features": [fval(v) for v in feats.iloc[j].tolist()],
                            "proba": fval(proba[j]), "margin": fval(margin[j]), "label": int(lab[j])})
    res = TP.detect_peaks_batch(histograms, model)
    out = {"records": records, "rows": row_out,
           "peaks": [[[int(p.position), int(p.count), int(p.left), int(p.right)] for p in pl] for pl in res]}
    return {"fn": "detect_peaks_batch", "model": model_name, "histograms": hs}, out


def xgb_case(rows, model_name):
    model = MODELS[model_name]
    feats = pd.DataFrame(rows, columns=TP.FEATURE_COLUMNS).astype(float)
    lab = model.predict(feats)
    proba = model.predict_proba(feats)[:, 1]
    margin = np.asarray(model.predict(feats, output_margin=True)).ravel()
    arr = np.asarray(feats, dtype=float)
    inp = {"fn": "xgb_predict", "model": model_name, "caller": "random",
           "features": [[fval(v) for v in row] for row in arr.tolist()]}
    out = {"labels": [int(v) for v in np.asarray(lab).tolist()], "proba": [fval(v) for v in proba.tolist()],
           "margin": [fval(v) for v in margin.tolist()]}
    return inp, out


def random_histogram(rng):
    """Shapes that stress the numerics: sparse/dense, plateaus, equal heights, big counts,
    large genomic positions, shuffled dict order, zero counts."""
    kind = rng.integers(0, 9)
    base = int(rng.choice([0, 1, 1000, 8_704_428, 248_000_000, int(rng.integers(0, 3e8))]))
    if kind == 0:  # single position
        h = {base: int(rng.integers(1, 50))}
    elif kind == 1:  # few sparse positions with equal counts (ties)
        k = int(rng.integers(2, 12))
        pos = rng.choice(int(rng.integers(k, 300)), size=k, replace=False)
        c = int(rng.integers(1, 4))
        h = {base + int(p): c for p in pos}
    elif kind == 2:  # plateau(s)
        h = {}
        for _ in range(int(rng.integers(1, 4))):
            s = int(rng.integers(0, 200))
            w = int(rng.integers(1, 15))
            c = int(rng.integers(1, 6))
            for p in range(s, s + w):
                h[base + p] = c
    elif kind == 3:  # mixture of gaussians, realistic polyA-like
        n = int(rng.integers(5, 3000))
        centers = rng.integers(0, 500, size=int(rng.integers(1, 4)))
        sd = rng.uniform(0.5, 20)
        vals = np.round(rng.choice(centers, size=n) + rng.normal(0, sd, size=n)).astype(int)
        h = {}
        for v in vals:
            h[base + int(v)] = h.get(base + int(v), 0) + 1
    elif kind == 4:  # wide sparse (long span), many height-1 peaks
        n = int(rng.integers(2, 200))
        span = int(rng.integers(100, 20000))
        vals = rng.integers(0, span, size=n)
        h = {}
        for v in vals:
            h[base + int(v)] = h.get(base + int(v), 0) + 1
    elif kind == 5:  # big counts on a few positions (n > 8192 -> chunking paths)
        k = int(rng.integers(1, 30))
        pos = rng.integers(0, 400, size=k)
        h = {base + int(p): int(rng.integers(1, 5000)) for p in pos}
    elif kind == 6:  # random small dense, periodic -> ties within distance 10
        period = int(rng.integers(2, 15))
        h = {base + i * period: int(rng.choice([1, 2, 2, 3])) for i in range(int(rng.integers(2, 40)))}
    elif kind == 7:  # uniform random counts incl. zeros
        span = int(rng.integers(1, 120))
        h = {base + i: int(rng.integers(0, 4)) for i in range(span)}
        if all(v == 0 for v in h.values()):
            h[base] = 1
    else:  # heavy tails: geometric counts
        span = int(rng.integers(1, 800))
        h = {base + int(i): int(rng.geometric(0.3)) for i in rng.integers(0, span, size=int(rng.integers(1, 200)))}
    items = list(h.items())
    if rng.random() < 0.5:  # dict order is load-bearing for float sums
        rng.shuffle(items)
    return dict((int(p), int(c)) for p, c in items)


def random_rows(rng, model_name, n):
    """Feature rows: realistic ranges, NaNs, and values on / one f32 ulp around split thresholds."""
    thresholds = {}
    booster = MODELS[model_name].get_booster()
    for tree in json.loads(booster.save_raw("json"))["learner"]["gradient_booster"]["model"]["trees"]:
        for lc, si, sc in zip(tree["left_children"], tree["split_indices"], tree["split_conditions"]):
            if lc != -1:
                thresholds.setdefault(si, []).append(sc)
    rows = []
    for _ in range(n):
        r = [rng.uniform(0, 5000) * rng.random() ** 3, rng.normal(0, 2), float(rng.integers(1, 40)),
             rng.uniform(0.5, 60), rng.uniform(0, 6), rng.uniform(0, 30), float(rng.integers(1, 2000)),
             rng.random()]
        for f in range(8):
            u = rng.random()
            if u < 0.35 and thresholds.get(f):
                t = np.float32(rng.choice(thresholds[f]))
                v = rng.integers(0, 4)
                if v == 0:
                    r[f] = float(t)
                elif v == 1:
                    r[f] = float(np.nextafter(t, np.float32(np.inf)))
                elif v == 2:
                    r[f] = float(np.nextafter(t, np.float32(-np.inf)))
                else:  # f64 just below the f32 threshold: rounds to it
                    r[f] = float(np.nextafter(float(t), -np.inf))
            elif u < 0.4:
                r[f] = float("nan")
        rows.append(r)
    return rows


def main():
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261029
    n_batch = int(sys.argv[3]) if len(sys.argv) > 3 else 1500
    n_xgb = int(sys.argv[4]) if len(sys.argv) > 4 else 200
    rng = np.random.default_rng(seed)
    with gzip.open(OUT, "wt") as f:
        for k in range(n_batch):
            model = "polya" if rng.random() < 0.6 else "tss"
            hs = [random_histogram(rng) if rng.random() > 0.05 else None for _ in range(int(rng.integers(1, 4)))]
            inp, out = batch_case(hs, model)
            f.write(dumps({"unit": "U10a", "v": 1, "case": f"random-{MODE}:{seed}:-:-:detect_peaks_batch:n.{k}",
                           "in": inp, "out": out}) + "\n")
        for k in range(n_xgb):
            model = "polya" if k % 2 == 0 else "tss"
            inp, out = xgb_case(random_rows(rng, model, 25), model)
            f.write(dumps({"unit": "U10a", "v": 1, "case": f"random-{MODE}:{seed}:-:-:xgb_predict:n.{k}",
                           "in": inp, "out": out}) + "\n")
    print(OUT, os.path.getsize(OUT))


if __name__ == "__main__":
    main()
