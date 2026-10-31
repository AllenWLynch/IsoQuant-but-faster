"""pysam/pyfaidx side of `bam_bench check|fasta`, plus a pysam full-fetch timing baseline.

    python bam_check.py regions <bam|fa> <chr> <n> <seed> > regions.tsv   random 1-based regions
    python bam_check.py check <bam> <regions.tsv>
    python bam_check.py fasta <fa> <regions.tsv>
    python bam_check.py time <bam> <chr>
"""
import random
import sys
import time
import zlib

FNV0 = 0xcbf29ce484222325
MASK = (1 << 64) - 1


def fnv(h, bs):
    for b in bs:
        h = ((h ^ b) * 0x100000001b3) & MASK
    return h


def regions(path):
    for line in open(path):
        c, s, e = line.rstrip("\n").split("\t")
        yield c, int(s), int(e)


def main():
    cmd = sys.argv[1]
    if cmd == "regions":
        import pysam
        path, chrom, n, seed = sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
        if path.endswith(".bam"):
            ln = pysam.AlignmentFile(path).get_reference_length(chrom)
        else:
            ln = pysam.FastaFile(path).get_reference_length(chrom)
        rng = random.Random(seed)
        for k in range(n):
            w = rng.choice([1, 2, 59, 60, 61, 120, 1000, 10_000, 100_000, rng.randrange(1, 200_000)])
            s = rng.randrange(1, ln - w + 2)
            if k % 50 == 0:  # past the end: clipped
                s = ln - rng.randrange(0, 100)
            print("%s\t%d\t%d" % (chrom, s, s + w - 1))
    elif cmd == "check":
        import pysam
        bam = pysam.AlignmentFile(sys.argv[2])
        for c, s, e in regions(sys.argv[3]):
            h, n = 0, 0
            for r in bam.fetch(c, s - 1, e):
                n += 1
                # end: RecordLite::end0 (CIGAR ref length; pysam reference_end is None if unmapped)
                end = r.reference_start + sum(l for o, l in (r.cigartuples or []) if o in (0, 2, 3, 7, 8))
                hp = r.get_tag("HP", with_value_type=True) if r.has_tag("HP") else None
                hp = "None" if hp is None else ('Some(Str("%s"))' % hp[0] if hp[1] == "Z" else 'Some(Int(%d))' % hp[0])
                h = zlib.crc32(("%s\t%d\t%d\t%d\t%s\t%s\n" % (r.query_name, r.reference_start, r.flag, end,
                                                          r.query_sequence or "", hp)).encode(), h)
            print("%s:%d-%d\t%d\t%08x" % (c, s, e, n, h))
    elif cmd == "fasta":
        from pyfaidx import Fasta
        fa = Fasta(sys.argv[2])
        for c, s, e in regions(sys.argv[3]):
            seq = str(fa[c][s - 1:e]).encode()
            print("%s:%d-%d\t%016x\t%d" % (c, s, e, fnv(FNV0, seq), len(seq)))
    elif cmd == "time":
        import pysam
        bam = pysam.AlignmentFile(sys.argv[2])
        t = time.time()
        n = acc = 0
        for r in bam.fetch(sys.argv[3]):
            n += 1
            acc += len(r.query_name) + r.flag + r.mapping_quality + (r.reference_end or 0) + len(r.cigartuples or ())
        s = time.time() - t
        print("pysam\tfull\t1\tnoseq\t%d\t%.3f\t%.0f\t-\tpython loop: qname,flag,mapq,reference_end,cigartuples" % (n, s, n / s))
        t = time.time()
        n = 0
        for r in bam.fetch(sys.argv[3]):
            n += len(r.query_sequence or "")
        s = time.time() - t
        print("pysam\tfull\t1\tseq\t-\t%.3f\t-\t-\tquery_sequence only" % s)


if __name__ == "__main__":
    main()
