"""python -m iqcap.run --units U1,U4,... --out DIR [--sample-rate R] -- <isoquant args>"""
import argparse
import datetime
import glob
import json
import os
import platform
import shutil
import socket
import subprocess
import sys
import time
from enum import Enum

from . import SCHEMA_VERSION, ISOQUANT_DIR, ensure_isoquant_on_path
from .sink import Sinks, fval

ALL_UNITS = ["U1", "U3", "U4", "U5", "U6", "U7", "U8", "U9", "U10a", "U10b", "U11", "U12", "U13", "U14",
             "U15a", "U15b", "U16", "U17", "U18", "U19", "U20"]


def parse_rates(s):
    rates = {"default": 1.0}
    if not s:
        return rates
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            rates[k.strip()] = float(v)
        else:
            rates["default"] = float(part)
    return rates


def cpu_info():
    d = {"model": None, "avx512f": False, "avx2": False, "flags_hash": None}
    try:
        import hashlib
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name") and d["model"] is None:
                    d["model"] = line.split(":", 1)[1].strip()
                if line.startswith("flags"):
                    flags = line.split(":", 1)[1].split()
                    d["avx512f"] = "avx512f" in flags
                    d["avx2"] = "avx2" in flags
                    d["flags_hash"] = hashlib.sha1(" ".join(flags).encode()).hexdigest()
                    break
    except OSError:
        pass
    return d


def versions():
    v = {"python": platform.python_version()}
    for m in ("numpy", "scipy", "xgboost", "pandas", "pysam", "gffutils", "pyfaidx"):
        try:
            v[m] = __import__(m).__version__
        except Exception:
            v[m] = None
    return v


def isoquant_info():
    d = {"path": ISOQUANT_DIR, "version": None, "commit": None, "dirty": None}
    try:
        d["version"] = open(os.path.join(ISOQUANT_DIR, "VERSION")).read().strip()
    except OSError:
        pass
    try:
        d["commit"] = subprocess.check_output(["git", "-C", ISOQUANT_DIR, "rev-parse", "HEAD"], text=True).strip()
        st = subprocess.check_output(["git", "-C", ISOQUANT_DIR, "status", "--porcelain", "--untracked-files=no"], text=True)
        d["dirty"] = bool(st.strip())
    except Exception:
        pass
    return d


def jsonable_args(args):
    out = {}
    for k in sorted(vars(args)):
        v = getattr(args, k)
        if isinstance(v, Enum):
            out[k] = v.value
            out[k + "__name"] = v.name
        elif isinstance(v, bool) or v is None or isinstance(v, (int, str)):
            out[k] = v
        elif isinstance(v, float):
            out[k] = fval(v)
        elif isinstance(v, (list, tuple)) and all(isinstance(x, (int, str, bool)) or x is None for x in v):
            out[k] = list(v)
        else:
            out[k] = repr(v)
    return out


def finalize_unit_files(out_dir):
    res = {}
    for path in glob.glob(os.path.join(out_dir, "units", "*.jsonl.gz")):
        base = os.path.basename(path)
        parts = base.split(".")
        res.setdefault(parts[0], []).append(path)
    final = {}
    for unit, paths in res.items():
        pid_files = [p for p in paths if len(os.path.basename(p).split(".")) == 4]
        if len(pid_files) == 1 and len(paths) == 1:
            dst = os.path.join(out_dir, "units", unit + ".jsonl.gz")
            os.replace(pid_files[0], dst)
            final[unit] = ["units/" + unit + ".jsonl.gz"]
        else:
            final[unit] = sorted(os.path.relpath(p, out_dir) for p in paths)
    return final


def copy_outputs(iq_out, out_dir):
    dst = os.path.join(out_dir, "outputs")
    os.makedirs(dst, exist_ok=True)
    copied = []
    pats = ["*.transcript_models.gtf", "*.discovered_transcript_counts.tsv", "*.discovered_transcript_grouped*",
            "*.discovered_gene_counts.tsv", "*.discovered_gene_grouped*"]
    for sd in glob.glob(os.path.join(iq_out, "*")):
        if not os.path.isdir(sd):
            continue
        for p in pats:
            for f in glob.glob(os.path.join(sd, p)):
                shutil.copy2(f, dst)
                copied.append(os.path.basename(f))
    return sorted(copied)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--" in argv:
        i = argv.index("--")
        own, iq_args = argv[:i], argv[i + 1:]
    else:
        own, iq_args = argv, []
    ap = argparse.ArgumentParser(prog="iqcap.run")
    ap.add_argument("--units", default="all")
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", default=None)
    ap.add_argument("--slice", default=None, dest="slice_name")
    ap.add_argument("--sample-rate", default="1.0")
    ap.add_argument("--block-rate", type=float, default=1.0)
    ap.add_argument("--block-cap", type=int, default=200000)
    ap.add_argument("--slice-def", default=None, help="JSON file with the slice definition (from iqcap.slice)")
    a = ap.parse_args(own)

    units = ALL_UNITS if a.units == "all" else [u.strip() for u in a.units.split(",") if u.strip()]
    out_dir = os.path.abspath(a.out)
    os.makedirs(out_dir, exist_ok=True)
    sample = a.sample or os.path.basename(os.path.dirname(out_dir))
    slice_name = a.slice_name or os.path.basename(out_dir)
    rates = parse_rates(a.sample_rate)

    ensure_isoquant_on_path()
    import isoquant
    from . import hooks

    sinks = Sinks(out_dir, units, rates, a.block_rate, a.block_cap)
    ctx = hooks.Ctx(sinks, sample, slice_name, out_dir)
    hooks.install(ctx)

    captured = {}
    orig_run = isoquant.run_pipeline

    def run_pipeline(args):
        captured["args"] = args
        return orig_run(args)
    isoquant.run_pipeline = run_pipeline

    t0 = time.time()
    rc = 0
    try:
        isoquant.main(iq_args)
    except SystemExit as e:
        rc = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
    finally:
        hooks._flush_u8() if ctx.u8_pending else None
        sinks.close()
    wall = time.time() - t0

    files = finalize_unit_files(out_dir)
    args = captured.get("args")
    if args is not None:
        from .sink import dumps
        with open(os.path.join(out_dir, "params.json"), "w") as f:
            f.write(dumps(jsonable_args(args), sort_keys=True) + "\n")
    outputs = copy_outputs(args.output, out_dir) if args is not None else []
    slice_def = None
    if a.slice_def and os.path.exists(a.slice_def):
        slice_def = json.load(open(a.slice_def))
    units_meta = {}
    for u in units:
        units_meta[u] = {"files": files.get(u, []), "cases": sum(sinks.counts[u].values()),
                         "by_fn": dict(sinks.counts[u])}
    manifest = {
        "v": SCHEMA_VERSION, "tool": "iqcap",
        "created": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sample": sample, "slice": slice_name,
        "isoquant": isoquant_info(), "versions": versions(), "cpu": cpu_info(), "host": socket.gethostname(),
        "slice_def": slice_def, "isoquant_args": iq_args,
        "iqcap_args": {"units": units, "sample_rate": rates, "block_rate": a.block_rate, "block_cap": a.block_cap},
        "units": units_meta, "blocks": ctx.blocks_meta, "skipped_blocks": ctx.skipped_blocks,
        "outputs": outputs, "exit_code": rc, "wall_seconds": round(wall, 1),
    }
    with open(os.path.join(out_dir, "fixture_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    return rc


if __name__ == "__main__":
    sys.exit(main())
