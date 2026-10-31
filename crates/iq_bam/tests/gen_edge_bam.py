"""Generates tests/data/edge.bam (+.bai) and edge.expected.jsonl.gz: pysam's view of edge-case
records for iq_bam's integration test. Run with the project env:
    ../env/bin/python crates/iq_bam/tests/gen_edge_bam.py
The expected output is a genuine U2 fixture file (SCHEMA §4 U2), tests/data/edge/units/U2.jsonl.gz,
with `in.bam` relative to the crate dir (cargo test's cwd).
"""
import array
import gzip
import json
import os
import random

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data")
IUPAC = "=ACMGRSVTWYHKDBN"


def tag_value(a, t):
    try:
        v, typ = a.get_tag(t, with_value_type=True)
    except KeyError:
        return None
    if isinstance(v, (bytes, bytearray)):
        v = v.decode()
    elif not isinstance(v, (int, float, str)):
        v = str(v)  # array.array: what str() of the tag gives a read grouper
    return {"type": typ, "value": v}


def main():
    rng = random.Random(7)
    os.makedirs(OUT, exist_ok=True)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "c1", "LN": 3_000_000}, {"SN": "c2", "LN": 100_000}]}
    recs = []

    def mk(name, tid, pos, cigar, seq, flag=0, mapq=60, qual=True, tags=()):
        a = pysam.AlignedSegment()
        a.query_name = name
        a.flag = flag
        a.reference_id = tid
        a.reference_start = pos
        a.mapping_quality = mapq
        if cigar is not None:
            a.cigartuples = cigar
        a.query_sequence = seq
        if seq and qual:
            a.query_qualities = array.array("B", [30] * len(seq))
        for t, v, ty in tags:
            if ty is None:
                a.set_tag(t, v)
            else:
                a.set_tag(t, v, value_type=ty)
        recs.append(a)

    # IUPAC bases in every position, odd length
    mk("iupac", 0, 100, [(0, 16)], IUPAC, tags=[("HP", 1, "i")])
    mk("iupac_odd", 0, 101, [(4, 2), (0, 13)], IUPAC[1:], tags=[("HP", 2, "C")])
    mk("hp_s", 0, 102, [(0, 4)], "ACGT", tags=[("HP", -3, "s")])
    mk("hp_c", 0, 102, [(0, 4)], "ACGT", tags=[("HP", -1, "c")])
    mk("hp_I", 0, 102, [(0, 4)], "ACGT", tags=[("HP", 4000000000, "I")])
    mk("hp_S", 0, 102, [(0, 4)], "ACGT", tags=[("HP", 65000, "S")])
    mk("hp_Z", 0, 103, [(0, 4)], "ACGT", tags=[("HP", "2", "Z")])
    mk("hp_A", 0, 103, [(0, 4)], "ACGT", tags=[("HP", "x", "A")])
    mk("hp_f", 0, 103, [(0, 4)], "ACGT", tags=[("HP", 1.5, "f")])
    mk("hp_B", 0, 103, [(0, 4)], "ACGT", tags=[("HP", array.array("B", [1, 2, 3]), None)])
    mk("hp_none", 0, 104, [(0, 4)], "ACGT", tags=[("XX", 1, "i")])
    # no sequence / no quality
    mk("noseq", 0, 105, [(0, 50)], None, flag=256, mapq=0)
    mk("noqual", 0, 105, [(0, 5)], "ACGTA", qual=False)
    # placed unmapped read: bam_endpos = pos + 1
    mk("unmapped_placed", 0, 2000, None, "ACGTACGT", flag=4, mapq=0)
    # mapped, zero reference length CIGAR: endpos = pos + 1
    mk("insonly", 0, 3000, [(4, 3), (1, 5)], "ACGTACGT")
    # long read crossing many 16 kb windows, with an intron
    mk("spliced", 0, 10_000, [(0, 100), (3, 200_000), (0, 100)], "A" * 200, flag=16)
    # >65535 CIGAR ops: written as CG:B:I by htslib
    n = 35_000
    mk("longcigar", 0, 50_000, [(0, 1), (1, 1)] * n, "ACGT" * (n // 2), tags=[("HP", "1", "Z")])
    # bulk: random reads
    for k in range(3000):
        pos = rng.randrange(0, 2_900_000)
        ln = rng.randrange(1, 3000)
        if rng.random() < 0.2:
            cig = [(0, ln // 2), (3, rng.randrange(10, 50_000)), (0, ln - ln // 2)]
        elif rng.random() < 0.1:
            cig = [(4, 10), (0, ln), (2, 3), (0, 7), (5, 20)]
        else:
            cig = [(0, ln)]
        qlen = sum(l for o, l in cig if o in (0, 1, 4, 7, 8))
        seq = "ACGTN"[k % 5] * qlen  # compresses well; base decoding is covered above
        flag = rng.choice([0, 16, 256, 272, 2048, 2064])
        tags = [("HP", rng.choice([1, 2]), "i")] if rng.random() < 0.5 else []
        mk("r%d" % k, 0, pos, cig, seq, flag=flag, mapq=rng.randrange(0, 61), tags=tags)
    for k in range(200):
        mk("c2r%d" % k, 1, rng.randrange(0, 99_000), [(0, 500)], "C" * 500)
    # unplaced unmapped at the end
    mk("unplaced", -1, -1, None, "ACGT", flag=4, mapq=0)

    recs.sort(key=lambda a: (a.reference_id if a.reference_id >= 0 else 1 << 30, a.reference_start))
    path = os.path.join(OUT, "edge.bam")
    with pysam.AlignmentFile(path, "wb", header=header) as f:
        for a in recs:
            f.write(a)
    pysam.index(path)

    regions = [("c1", 100, 110), ("c1", 101, 101), ("c1", 1, 100), ("c1", 1, 101), ("c1", 2000, 2001), ("c1", 2001, 2001),
               ("c1", 2002, 2002), ("c1", 3001, 3001), ("c1", 3000, 3000), ("c1", 3002, 3002),
               ("c1", 150_000, 150_100), ("c1", 210_111, 210_111), ("c1", 50_001, 50_001),
               ("c1", 1, 3_000_000), ("c2", 1, 100_000), ("c2", 99_990, 100_000), ("c1", 16_384, 16_385)]
    for _ in range(150):
        s = rng.randrange(1, 2_999_000)
        regions.append(("c1", s, s + rng.randrange(0, 200_000)))
    bam = pysam.AlignmentFile(path, "rb")
    os.makedirs(os.path.join(OUT, "edge", "units"), exist_ok=True)
    with gzip.open(os.path.join(OUT, "edge", "units", "U2.jsonl.gz"), "wt") as fh:
        for chrom, s, e in regions:
            region = "%s:%d-%d" % (chrom, s, e)
            for n, rec in enumerate(bam.fetch(chrom, s - 1, e)):
                d = {"qname": rec.query_name, "flag": rec.flag, "ref_id": rec.reference_id,
                     "ref_name": rec.reference_name, "pos0": rec.reference_start, "ref_end0": rec.reference_end,
                     "mapq": rec.mapping_quality, "cigar": [[o, l] for o, l in (rec.cigartuples or [])],
                     "seq": rec.query_sequence,
                     "qual_len": len(rec.query_qualities) if rec.query_qualities is not None else 0, "tags": {}}
                tv = tag_value(rec, "HP")
                if tv is not None:
                    d["tags"]["HP"] = tv
                case = {"unit": "U2", "v": 1, "case": "edge:%s:%d" % (region, n),
                        "in": {"fn": "fetch", "bam": "tests/data/edge.bam", "region": region, "n": n}, "out": d}
                fh.write(json.dumps(case, separators=(",", ":")) + "\n")
    print("unmapped", bam.unmapped, "mapped", bam.mapped)
    with open(os.path.join(OUT, "edge.unmapped.txt"), "w") as fh:
        fh.write("%d\n" % bam.unmapped)


if __name__ == "__main__":
    main()
