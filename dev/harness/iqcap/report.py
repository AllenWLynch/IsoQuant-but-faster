"""python -m iqcap.report [FIXTURES_DIR] [--md] : cases per unit per sample (and per slice) from the manifests."""
import argparse
import glob
import json
import os
import sys
from collections import defaultdict

UNITS = ["U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8", "U9", "U10a", "U10b", "U11", "U12", "U13", "U14",
         "U15a", "U15b", "U16", "U17", "U18", "U19", "U20"]


def du(path):
    tot = 0
    for root, _, files in os.walk(path):
        for f in files:
            try:
                tot += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return tot


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("root", nargs="?", default=os.path.join(os.path.dirname(__file__), "..", "..", "..", "fixtures"))
    ap.add_argument("--by-slice", action="store_true")
    a = ap.parse_args(argv)
    root = os.path.abspath(a.root)
    table = defaultdict(lambda: defaultdict(int))
    slices = defaultdict(list)
    for m in sorted(glob.glob(os.path.join(root, "*", "*", "fixture_manifest.json"))):
        d = json.load(open(m))
        key = d["sample"] if not a.by_slice else "%s/%s" % (d["sample"], d["slice"])
        slices[d["sample"]].append(d["slice"])
        for u, um in d["units"].items():
            table[key][u] += um["cases"]

    import gzip
    for u2 in glob.glob(os.path.join(root, "*", "u2", "units", "U2.jsonl.gz")):
        s = u2.split(os.sep)[-4]
        table[s if not a.by_slice else s + "/u2"]["U2"] += sum(1 for _ in gzip.open(u2, "rt"))
    ann = os.path.join(root, "annotation", "U1.annotation.jsonl.gz")
    if os.path.exists(ann):
        for c in table:
            table[c]["U1"] += 0
        n_ann = sum(1 for _ in gzip.open(ann, "rt"))
    cols = sorted(table)
    print("| unit | " + " | ".join(cols) + " |")
    print("|---|" + "---|" * len(cols))
    for u in UNITS:
        print("| %s | " % u + " | ".join(str(table[c].get(u, 0)) for c in cols) + " |")
    print()
    for s in sorted(slices):
        print("- %s: slices %s; size %.1f MB" % (s, ", ".join(sorted(slices[s])), du(os.path.join(root, s)) / 1e6))
    if os.path.exists(ann):
        print("- U1 annotation dump: %d chromosome case(s), %s (%.1f MB)" % (n_ann, ann, os.path.getsize(ann) / 1e6))
    print("- total fixtures size: %.1f MB" % (du(root) / 1e6))


if __name__ == "__main__":
    sys.exit(main())
