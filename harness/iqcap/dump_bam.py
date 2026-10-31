"""python -m iqcap.dump_bam --bam X.bam --region 8:1001-2000 --out DIR [--sample S] [--tags HP] [--max N]

U2 fixture: pysam walk of a region (as BAMOnlineMerger fetches it), no IsoQuant run.
Writes DIR/units/U2.jsonl.gz (SCHEMA §4 U2).
"""
import argparse
import gzip
import os
import sys

import pysam

from . import SCHEMA_VERSION
from .sink import dumps
from .slice import parse_region


def tag_value(a, t):
    try:
        v, typ = a.get_tag(t, with_value_type=True)
    except KeyError:
        return None
    if isinstance(v, (bytes, bytearray)):
        v = v.decode()
    elif not isinstance(v, (int, float, str)):
        v = list(v) if hasattr(v, "__iter__") else str(v)
    return {"type": typ, "value": v}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="iqcap.dump_bam")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--region", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--tags", default="HP")
    ap.add_argument("--max", type=int, default=0)
    a = ap.parse_args(argv)
    chrom, s, e = parse_region(a.region)
    tags = [t for t in a.tags.split(",") if t]
    os.makedirs(os.path.join(a.out, "units"), exist_ok=True)
    region = "%s:%d-%d" % (chrom, s, e)
    n = 0
    bam = pysam.AlignmentFile(a.bam, "rb")
    with gzip.open(os.path.join(a.out, "units", "U2.jsonl.gz"), "wt", encoding="utf-8") as fh:
        for rec in bam.fetch(chrom, s - 1, e):
            out = {"qname": rec.query_name, "flag": rec.flag, "ref_id": rec.reference_id,
                   "ref_name": rec.reference_name, "pos0": rec.reference_start, "ref_end0": rec.reference_end,
                   "mapq": rec.mapping_quality, "cigar": [[o, l] for o, l in (rec.cigartuples or [])],
                   "seq": rec.query_sequence, "qual_len": len(rec.query_qualities) if rec.query_qualities is not None else 0,
                   "tags": {}}
            for t in tags:
                tv = tag_value(rec, t)
                if tv is not None:
                    out["tags"][t] = tv
            case = {"unit": "U2", "v": SCHEMA_VERSION, "case": "%s:%s:%d" % (a.sample, region, n),
                    "in": {"fn": "fetch", "bam": os.path.abspath(a.bam), "region": region, "n": n}, "out": out}
            fh.write(dumps(case) + "\n")
            n += 1
            if a.max and n >= a.max:
                break
    print("U2: %d records from %s %s" % (n, a.bam, region))


if __name__ == "__main__":
    sys.exit(main())
