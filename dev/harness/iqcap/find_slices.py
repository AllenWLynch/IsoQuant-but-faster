"""python -m iqcap.find_slices --bam X.bam --annotation fixtures/annotation/U1.annotation.jsonl.gz --chr 8 --out cand.json

Approximates IsoQuant's pass A (gap-separated clusters bridged by shared genes, then
split_coverage_regions with IsoQuant's own static function) to propose slice regions per design §10.2:
(a) ordinary loci, (b) deepest block, (c) multi-gene both-strand cluster, (d) annotation-free region,
(e) coverage-valley split. Only used to choose slices; the capture itself is authoritative.
"""
import argparse
import bisect
import gzip
import json
import sys
from collections import defaultdict

import pysam

from . import ensure_isoquant_on_path


class FakeStorage:
    def __init__(self, cov, n):
        self.coverage_dict = cov
        self.n = n

    def get_read_count(self):
        return self.n


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--chr", default="8")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    ensure_isoquant_on_path()
    from isoquant_lib.alignment.alignment_processor import AlignmentCollector

    genes = []
    for line in gzip.open(a.annotation, "rt"):
        c = json.loads(line)
        if c["out"]["chr"] == a.chr:
            for g in c["out"]["genes"]:
                genes.append((g["start"], g["end"], g["id"], g["strand"], len(g["transcripts"])))
    genes.sort()
    gstarts = [g[0] for g in genes]
    maxlen = max(g[1] - g[0] for g in genes) + 1

    def genes_in(s, e):  # 1-based closed overlap (approximation of gffutils region)
        i = bisect.bisect_left(gstarts, s - maxlen)
        return [g for g in genes[i:bisect.bisect_right(gstarts, e)] if g[1] >= s and g[0] <= e]

    bam = pysam.AlignmentFile(a.bam, "rb")
    clusters = []
    cur = None
    for r in bam.fetch(a.chr):
        if r.is_unmapped:
            continue
        s, e = r.reference_start, r.reference_end - 1
        if cur is not None and s > cur["e"]:
            pg = set(g[2] for g in genes_in(cur["s"] + 1, cur["e"] + 1))
            ng = set(g[2] for g in genes_in(s + 1, e + 1))
            if not (pg & ng):
                clusters.append(cur)
                cur = None
        if cur is None:
            cur = {"s": s, "e": e, "n": 0, "cov": defaultdict(int), "spliced": 0}
        cur["e"] = max(cur["e"], e)
        cur["n"] += 1
        if r.cigartuples and any(op == 3 for op, _ in r.cigartuples):
            cur["spliced"] += 1
        for b in range(s // 256, e // 256 + 1):
            cur["cov"][b] += 1
    if cur:
        clusters.append(cur)

    res = []
    for c in clusters:
        splits = AlignmentCollector.split_coverage_regions((c["s"], c["e"]), FakeStorage(c["cov"], c["n"]))
        gs = genes_in(c["s"] + 1, c["e"] + 1)
        res.append({"region": "%s:%d-%d" % (a.chr, c["s"] + 1, c["e"] + 1), "len": c["e"] - c["s"] + 1,
                    "n": c["n"], "spliced": c["spliced"], "n_genes": len(gs),
                    "strands": "".join(sorted(set(g[3] for g in gs))), "n_tx": sum(g[4] for g in gs),
                    "n_blocks": len(splits), "max_cov": max(x[1] for x in splits)})
    json.dump(res, open(a.out, "w"), indent=0)
    print("%d clusters" % len(res), file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
