"""python -m iqcap.verify DIR_A DIR_B : compare two IsoQuant output dirs (all .gtf/.tsv/.tsv.gz), ignoring
the '# Command line:' header line. Exit 1 on any difference."""
import gzip
import os
import sys


def lines(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as f:
        return [l for l in f if not l.startswith("# Command line:")]


def main(argv=None):
    a, b = (argv or sys.argv[1:])[:2]
    fa, fb = set(os.listdir(a)), set(os.listdir(b))
    bad = 0
    for f in sorted(fa | fb):
        if f not in fa or f not in fb:
            print("MISSING", f, "in", a if f not in fa else b)
            bad += 1
            continue
        if lines(os.path.join(a, f)) != lines(os.path.join(b, f)):
            print("DIFF", f)
            bad += 1
        else:
            print("same", f)
    print("files: %d, differing: %d" % (len(fa | fb), bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
