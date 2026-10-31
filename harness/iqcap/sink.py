"""Fixture sinks: JSON encoding, deterministic sampling, per-PID gzip JSONL files, index side table."""
import gzip
import hashlib
import json
import math
import os
from collections import defaultdict

from . import SCHEMA_VERSION


def fval(x):
    """Float encoding: finite -> float, non-finite -> string (SCHEMA §0.4)."""
    x = float(x)
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Infinity" if x > 0 else "-Infinity"
    return x


def dumps(obj, sort_keys=False):
    """The single JSON writer for all fixture payloads (SCHEMA §0.2)."""
    return json.dumps(obj, separators=(",", ":"), ensure_ascii=True, allow_nan=False, sort_keys=sort_keys)


def keep(key, rate):
    if rate >= 1.0:
        return True
    if rate <= 0.0:
        return False
    h = int(hashlib.sha1(key.encode()).hexdigest()[:16], 16)
    return h / float(1 << 64) < rate


class Sinks:
    def __init__(self, out_dir, units, rates, block_rate=1.0, block_cap=200000):
        self.out_dir = out_dir
        self.units = set(units)
        self.rates = rates  # dict unit -> rate, 'default' key
        self.block_rate = block_rate
        self.block_cap = block_cap
        self.files = {}
        self.counts = defaultdict(lambda: defaultdict(int))
        self.index_written = set()
        os.makedirs(os.path.join(out_dir, "units"), exist_ok=True)
        os.makedirs(os.path.join(out_dir, "index"), exist_ok=True)

    def wants(self, unit):
        return unit in self.units

    def rate(self, unit):
        return self.rates.get(unit, self.rates.get("default", 1.0))

    def _fh(self, unit):
        pid = os.getpid()
        k = (unit, pid)
        fh = self.files.get(k)
        if fh is None:
            path = os.path.join(self.out_dir, "units", "%s.%d.jsonl.gz" % (unit, pid))
            fh = gzip.open(path, "wt", encoding="utf-8", compresslevel=6)
            self.files[k] = fh
        return fh

    def write(self, unit, case, inp, out):
        line = dumps({"unit": unit, "v": SCHEMA_VERSION, "case": case, "in": inp, "out": out})
        self._fh(unit).write(line + "\n")
        self.counts[unit][inp.get("fn", "")] += 1

    def write_index(self, payload_bytes, sha):
        if sha in self.index_written:
            return
        path = os.path.join(self.out_dir, "index", sha + ".json")
        if not os.path.exists(path):
            tmp = path + ".tmp.%d" % os.getpid()
            with open(tmp, "wb") as f:
                f.write(payload_bytes)
            os.replace(tmp, path)
        self.index_written.add(sha)

    def close(self):
        for fh in self.files.values():
            fh.close()
        self.files = {}


class BlockWriter:
    """block/<sample>/<chr>/<block_no>.jsonl.gz"""

    def __init__(self, out_dir, sample, chr_id, block_no):
        d = os.path.join(out_dir, "block", sample, chr_id)
        os.makedirs(d, exist_ok=True)
        self.path = os.path.join(d, "%d.jsonl.gz" % block_no)
        self.fh = gzip.open(self.path, "wt", encoding="utf-8", compresslevel=6)

    def write(self, obj):
        self.fh.write(dumps(obj) + "\n")

    def close(self):
        self.fh.close()
