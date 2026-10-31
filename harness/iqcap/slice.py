"""python -m iqcap.slice --bam IN.bam --regions 8:1-1000,8:5000-9000 --out OUT.bam [--no-hp]

Region slicing of a coordinate-sorted indexed BAM (like `samtools view -b IN regions`; overlapping or
adjacent regions are merged, every record overlapping a merged region is written once, in input order)
plus deterministic HP tag synthesis: HP:Z = ["A","B","unphased"][int(sha1(qname)[:16],16) % 3]
(existing HP tags are replaced). Writes OUT.bam, OUT.bam.bai and OUT.slice.json (slice definition).
"""
import argparse
import hashlib
import json
import os
import sys

import pysam

HP_VALUES = ("A", "B", "unphased")


def hp_for(qname):
    return HP_VALUES[int(hashlib.sha1(qname.encode()).hexdigest()[:16], 16) % 3]


def parse_region(r):
    chrom, rng = r.rsplit(":", 1)
    s, e = rng.replace(",", "").split("-")
    return chrom, int(s), int(e)


def merge(regions):
    out = []
    for c, s, e in sorted(regions, key=lambda x: (x[0], x[1], x[2])):
        if out and out[-1][0] == c and s <= out[-1][2] + 1:
            out[-1] = (c, out[-1][1], max(out[-1][2], e))
        else:
            out.append((c, s, e))
    return out


def make_slice(bam, regions, out, hp=True, hp_int=False, subsample=1.0, keep_names=None):
    regs = merge([parse_region(r) for r in regions])
    keep = None
    if keep_names:
        with open(keep_names) as f:
            keep = set(l.strip() for l in f if l.strip())
    inb = pysam.AlignmentFile(bam, "rb")
    outb = pysam.AlignmentFile(out, "wb", template=inb)
    n = 0
    seen_end = {}
    for c, s, e in regs:
        prevs = seen_end.setdefault(c, [])
        for a in inb.fetch(c, s - 1, e):
            # skip records already written for an earlier merged region of the same chromosome
            if a.reference_end is not None and any(a.reference_start < pe and a.reference_end > ps for ps, pe in prevs):
                continue
            if subsample < 1.0 and int(hashlib.sha1(("sub:" + a.query_name).encode()).hexdigest()[:16], 16) / 2.0 ** 64 >= subsample:
                continue
            if keep is not None and a.query_name not in keep:
                continue
            if hp:
                v = hp_for(a.query_name)
                if hp_int:
                    a.set_tag("HP", {"A": 1, "B": 2}.get(v, 0), value_type="i")
                else:
                    a.set_tag("HP", v, value_type="Z")
            outb.write(a)
            n += 1
        prevs.append((s - 1, e))
    outb.close()
    inb.close()
    pysam.index(out)
    d = {"bam": os.path.abspath(bam), "regions": ["%s:%d-%d" % r for r in regs], "requested": list(regions),
         "hp_synthesis": ("HP:i A=1,B=2,unphased=0" if hp_int else "HP:Z sha1(qname)[:16] % 3 -> A|B|unphased")
         if hp else None, "slice_bam": os.path.abspath(out), "n_records": n,
         "subsample": subsample, "subsample_rule": "keep iff int(sha1('sub:'+qname)[:16],16)/2**64 < subsample"}
    if keep is not None:
        names = sorted(keep)
        d["keep_names"] = {"file": os.path.abspath(keep_names), "n": len(names),
                           "sha1": hashlib.sha1("\n".join(names).encode()).hexdigest(),
                           "rule": "keep iff qname in file (applied after subsample)"}
    with open(os.path.splitext(out)[0] + ".slice.json", "w") as f:
        json.dump(d, f, indent=1)
    return d


def main(argv=None):
    ap = argparse.ArgumentParser(prog="iqcap.slice")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--regions", required=True, help="comma-separated chr:start-end (1-based closed)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-hp", action="store_true")
    ap.add_argument("--subsample", type=float, default=1.0, help="deterministic read-name subsampling fraction")
    ap.add_argument("--keep-names", default=None, help="file of read names (one per line); only these are kept")
    ap.add_argument("--hp-int", action="store_true", help="HP:i variant (D-O2 test; IsoQuant 4.1 cannot run it)")
    a = ap.parse_args(argv)
    d = make_slice(a.bam, a.regions.split(","), a.out, hp=not a.no_hp, hp_int=a.hp_int, subsample=a.subsample,
                   keep_names=a.keep_names)
    print(json.dumps(d))


if __name__ == "__main__":
    sys.exit(main())
