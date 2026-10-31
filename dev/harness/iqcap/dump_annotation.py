"""python -m iqcap.dump_annotation --db isoquant.annotation.db --gff3 annotation.gff3 --out FILE [--chr 8,...]

U1 hook (a): the per-chromosome annotation view of the pipeline's gffutils DB, as IsoQuant 4.1 sees it
(genes = featuretype exactly 'gene'; transcripts = DEFAULT_RNA_TYPES children; exons after the
CDS-as-exon rule). Uses IsoQuant's own GeneInfo code per gene (prepare_profiles=False), no IsoQuant run.
Writes one case per chromosome (SCHEMA §4 U1 a).
"""
import argparse
import gzip
import os
import sys

from . import SCHEMA_VERSION, ensure_isoquant_on_path
from .sink import dumps


def attrs_list(f):
    return [[k, list(v)] for k, v in f.attributes.items()]


def dump_chr(db, chr_id, GeneInfo):
    genes = []
    max_end = list(db.conn.execute("select max(end) from features where seqid = ?", (chr_id,)))[0][0] or 1
    for g in db.features_of_type("gene", limit="%s:1-%d" % (chr_id, max_end), order_by="start"):
        if g.seqid != chr_id:
            continue
        gi = GeneInfo([g], db, delta=0, prepare_profiles=False)
        unordered = [t.id for t in db.children(g, featuretype=GeneInfo.DEFAULT_RNA_TYPES)]
        txs = []
        for t in db.children(g, featuretype=GeneInfo.DEFAULT_RNA_TYPES, order_by="start"):
            txs.append({
                "id": t.id, "type": t.featuretype, "start": t.start, "end": t.end, "strand": t.strand,
                "source": t.source, "attributes": attrs_list(t), "attrs": gi.feature_attributes.get(t.id, ""),
                "exons": [[int(s), int(e)] for s, e in gi.all_isoforms_exons.get(t.id, [])],
                "other_features": [[int(s), int(e), str(ft)] for s, e, ft in gi.other_features.get(t.id, [])],
                "children_order_index": unordered.index(t.id),
            })
        genes.append({"id": g.id, "start": g.start, "end": g.end, "strand": g.strand, "source": g.source,
                      "attributes": attrs_list(g), "attrs": gi.feature_attributes.get(g.id, ""),
                      "transcripts": txs})
    return genes


def main(argv=None):
    ap = argparse.ArgumentParser(prog="iqcap.dump_annotation")
    ap.add_argument("--db", required=True)
    ap.add_argument("--gff3", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chr", default=None, help="comma-separated chromosomes (default: all in the DB)")
    a = ap.parse_args(argv)
    ensure_isoquant_on_path()
    import gffutils
    from isoquant_lib.gene_info import GeneInfo
    db = gffutils.FeatureDB(a.db)
    if a.chr:
        chrs = a.chr.split(",")
    else:
        chrs = [r[0] for r in db.execute("select distinct seqid from features where featuretype='gene'")]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with gzip.open(a.out, "wt", encoding="utf-8") as fh:
        for c in chrs:
            genes = dump_chr(db, c, GeneInfo)
            case = {"unit": "U1", "v": SCHEMA_VERSION, "case": "annotation:%s" % c,
                    "in": {"fn": "dump_annotation", "db": os.path.abspath(a.db),
                           "gff3": os.path.abspath(a.gff3) if a.gff3 else None, "chr": c, "gene_types": "isoquant"},
                    "out": {"chr": c, "genes": genes}}
            fh.write(dumps(case) + "\n")
            print("chr %s: %d genes" % (c, len(genes)), file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
