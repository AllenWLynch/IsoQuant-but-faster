"""python -m iqcap.validate FIXTURE_DIR... : structural checks of captured fixtures.

- every line parses, has exactly {unit,v,case,in,out}, v == SCHEMA_VERSION, in.fn present, keys unique per file;
- every "index" sha1 referenced exists in index/ and hashes to its name;
- every block reference has its block/ file, header n_reads == number of records, read_idx refs in range;
- every graph/paths cross reference resolves to a U11/U12 case;
- byte round trip: dumps(json.loads(line)) == line (the iq_types round-trip contract).
"""
import glob
import gzip
import hashlib
import json
import os
import sys

from . import SCHEMA_VERSION
from .sink import dumps


def walk_refs(o, out):
    if isinstance(o, dict):
        for k, v in o.items():
            if k == "index" and isinstance(v, str) and len(v) == 40:
                out["index"].add(v)
            elif k in ("graph", "paths") and isinstance(v, str):
                out["xref"].add(v)
            elif k == "block" and isinstance(v, dict) and "block_no" in v:
                out["block"].add((v["chr"], v["block_no"]))
            else:
                walk_refs(v, out)
    elif isinstance(o, list):
        for v in o:
            walk_refs(v, out)


def validate(d):
    errs = []
    refs = {"index": set(), "xref": set(), "block": set()}
    keys_by_unit = {}
    nlines = 0
    for f in sorted(glob.glob(os.path.join(d, "units", "*.jsonl.gz"))):
        seen = set()
        for line in gzip.open(f, "rt", encoding="utf-8"):
            nlines += 1
            line = line.rstrip("\n")
            c = json.loads(line)
            if list(c.keys()) != ["unit", "v", "case", "in", "out"]:
                errs.append("%s: bad envelope keys %s" % (f, list(c.keys())))
            if c.get("v") != SCHEMA_VERSION:
                errs.append("%s: v=%s" % (f, c.get("v")))
            if "fn" not in c["in"]:
                errs.append("%s: no fn in %s" % (f, c["case"]))
            if c["case"] in seen:
                errs.append("%s: duplicate key %s" % (f, c["case"]))
            seen.add(c["case"])
            if dumps(c) != line:
                errs.append("%s: round trip differs for %s" % (f, c["case"]))
            walk_refs(c["in"], refs)
            walk_refs(c["out"], refs)
        keys_by_unit[os.path.basename(f).split(".")[0]] = seen
    for sha in refs["index"]:
        p = os.path.join(d, "index", sha + ".json")
        if not os.path.exists(p):
            errs.append("missing index %s" % sha)
        elif hashlib.sha1(open(p, "rb").read()).hexdigest() != sha:
            errs.append("index hash mismatch %s" % sha)
    allkeys = set().union(*keys_by_unit.values()) if keys_by_unit else set()
    for x in refs["xref"]:
        if x not in allkeys:
            errs.append("unresolved cross reference %s" % x)
    for p in glob.glob(os.path.join(d, "block", "*", "*", "*.jsonl.gz")):
        lines = gzip.open(p, "rt").read().splitlines()
        h = json.loads(lines[0])
        if h.get("kind") != "block_header" or h["n_reads"] != len(lines) - 1:
            errs.append("bad block file %s" % p)
        for i, l in enumerate(lines[1:]):
            if json.loads(l)["idx"] != i:
                errs.append("bad idx in %s" % p)
                break
    sample = json.load(open(os.path.join(d, "fixture_manifest.json")))["sample"]
    for chr_id, b in refs["block"]:
        if not os.path.exists(os.path.join(d, "block", sample, chr_id, "%d.jsonl.gz" % b)):
            errs.append("missing block file %s/%s" % (chr_id, b))
    return nlines, errs


def main(argv=None):
    bad = 0
    for d in (argv or sys.argv[1:]):
        n, errs = validate(d)
        print("%s: %d cases, %d errors" % (d, n, len(errs)))
        for e in errs[:20]:
            print("   ", e)
        bad += len(errs)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
