"""Extra U1(b) cases straight from gffutils + 4.1 GeneInfo, plus an exon_id table."""
import sys, os, gzip, random, hashlib
# project root holding IsoQuant/ (4.1 clone), iqrs/ and data/; default: the parent of the iqrs checkout
R = os.environ.get("IQRS_PROJECT", os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../..")))
sys.path.insert(0, R + "/IsoQuant"); sys.path.insert(0, R + "/iqrs/dev/harness")
import gffutils
from isoquant_lib.gene_info import GeneInfo
from isoquant_lib.utils.id_policy import FeatureIdStorage, SimpleIDDistributor
from iqcap.enc import gene_info_dump
from iqcap.sink import dumps
out = sys.argv[1]; n_per_chr = int(sys.argv[2])
os.makedirs(out + "/units", exist_ok=True); os.makedirs(out + "/index", exist_ok=True)
db = gffutils.FeatureDB(R + "/data/ref/isoquant.annotation.db")
random.seed(7)
chrs = [str(i) for i in range(1, 26)]
fh = gzip.open(out + "/units/U1.jsonl.gz", "wt")
k = 0
allg = {}
for g in db.features_of_type("gene"):
    allg.setdefault(g.seqid, []).append(g)
for c in chrs:
    genes = allg[c]
    regions = []
    for _ in range(n_per_chr):
        g = random.choice(genes)
        kind = random.random()
        if kind < 0.3:   # boundary probes: 1bp off either gene end, 0-based tuple as 4.1 passes it
            p = random.choice([g.start - 1, g.start, g.end, g.end + 1, g.start - 2])
            regions.append((p, p + random.choice([0, 0, 5, 300])))
        elif kind < 0.9:
            s = g.start + random.randint(-20000, 20000)
            regions.append((max(s, 0), max(s, 0) + random.choice([1, 500, 5000, 20000, 60000])))
        else:
            regions.append((0, random.choice([1, 5000, 30000])))
    for (s, e) in regions:
        gl = list(db.region(seqid=c, start=s, end=e, featuretype="gene"))
        if gl:
            gl = sorted(gl, key=lambda x: x.start)
            gi = GeneInfo(gl, db, 4)
        else:
            gi = GeneInfo.from_region(c, s, e, 4, None)
        b, sha = gene_info_dump(gi)
        with open(out + "/index/" + sha + ".json", "wb") as f:
            f.write(b)
        inp = {"fn": "get_gene_info_for_region", "chr": c, "region": [s, e], "delta": 4, "needs_reference": False}
        o = {"gene_ids": [g.id for g in gi.gene_db_list], "index": sha, "empty": bool(gi.empty())}
        fh.write(dumps({"unit": "U1", "v": 1, "case": "extra:gen:%s:-:get_gene_info_for_region:n.%d" % (c, k), "in": inp, "out": o}) + "\n")
        k += 1
    print(c, k, file=sys.stderr, flush=True)
fh.close()
