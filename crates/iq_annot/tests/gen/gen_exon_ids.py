"""Extra U1(b) cases straight from gffutils + 4.1 GeneInfo, plus an exon_id table."""
import sys, os, gzip, random, hashlib
R = "/n/data1/hms/dbmi/park/allen_l/lariat/isoquant_r"
sys.path.insert(0, R + "/IsoQuant"); sys.path.insert(0, R + "/iqrs/harness")
import gffutils
from isoquant_lib.gene_info import GeneInfo
from isoquant_lib.utils.id_policy import FeatureIdStorage, SimpleIDDistributor
from iqcap.enc import gene_info_dump
from iqcap.sink import dumps
out = sys.argv[1]
db = gffutils.FeatureDB(R + "/data/ref/isoquant.annotation.db")
chrs = [str(i) for i in range(1, 26)]
with open(out + "/exon_ids.tsv", "w") as f:
    for c in chrs:
        st = FeatureIdStorage(SimpleIDDistributor(), db, c, "exon")
        for (ch, s, e, strand), v in st.id_dict.items():
            f.write("%s\t%d\t%d\t%s\t%s\n" % (ch, s, e, strand, v))
