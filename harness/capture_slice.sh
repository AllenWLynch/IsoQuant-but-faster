#!/bin/bash
# Capture one slice: BAM slice (+HP synthesis) -> iqcap.run (all units) -> fixtures/<sample>/<slice>/
# Usage: capture_slice.sh SAMPLE SLICE REGIONS DTYPE [SAMPLE_RATE] [SUBSAMPLE] [UNITS]
#   REGIONS: comma-separated 8:s-e (1-based closed); DTYPE: nanopore|pacbio_ccs
#   SAMPLE_RATE: iqcap --sample-rate spec (default 1.0); SUBSAMPLE: read fraction kept (default 1.0)
# Env: STAGE (dir with genome.fa(.fai) + isoquant.annotation.db on local disk; staged to $TMPDIR if unset),
#      NO_CAPTURE=1 runs plain IsoQuant on the same slice; with VERIFY_DIR set, both modes copy all IsoQuant
#      outputs to $VERIFY_DIR/<slice>/{plain,capture} for the no-perturbation check (iqcap.verify).
#      Batch-2 overrides (defaults reproduce the M0a behaviour):
#        NOVEL_UNSPLICED=true|false (default false)  -> --report_novel_unspliced
#        SRC_BAM=path      source BAM to slice (default data/chr8/<SAMPLE>.chr8.bam)
#        KEEP_NAMES=file   only keep these read names in the slice (iqcap.slice --keep-names)
#        EXTRA_IQARGS="..." appended to the IsoQuant args
#        BLOCK_RATE / BLOCK_CAP  -> iqcap --block-rate / --block-cap (defaults: iqcap's)
#SBATCH -t 4:00:00
set -euo pipefail
SAMPLE=$1; SLICE=$2; REGIONS=$3; DTYPE=$4; RATE=${5:-1.0}; SUB=${6:-1.0}; UNITS=${7:-all}
P=/n/data1/hms/dbmi/park/allen_l/lariat/isoquant_r
H=$P/iqrs/harness
PY=$P/env/bin/python
export PATH=$P/env/bin:$PATH
OUT=$P/fixtures/$SAMPLE/$SLICE
WORK=${TMPDIR:-/tmp}/iqcap_${SAMPLE}_${SLICE}_$$
mkdir -p "$WORK" "$OUT/slice"
if [ -z "${STAGE:-}" ]; then
  STAGE=$WORK/stage; mkdir -p $STAGE
  cp $P/data/ref/isoquant.annotation.db $P/data/ref/genome.fa $P/data/ref/genome.fa.fai $STAGE/
fi
if [ "$DTYPE" = pacbio_ccs ]; then STRAT="sensitive_pacbio --fl_data"; else STRAT=sensitive_ont; fi
BAM=$OUT/slice/$SAMPLE.$SLICE.bam
if [ ! -s "$BAM.bai" ]; then
  cd $H && $PY -m iqcap.slice --bam ${SRC_BAM:-$P/data/chr8/$SAMPLE.chr8.bam} --regions "$REGIONS" --out $BAM --subsample $SUB ${KEEP_NAMES:+--keep-names $KEEP_NAMES}
fi
IQARGS="--reference $STAGE/genome.fa --genedb $STAGE/isoquant.annotation.db --bam $BAM --data_type $DTYPE
  --model_construction_strategy $STRAT --report_novel_unspliced ${NOVEL_UNSPLICED:-false} --read_group tag:HP --threads 1 --force ${EXTRA_IQARGS:-}"
cd $H
if [ "${NO_CAPTURE:-0}" = 1 ]; then
  $PY $P/IsoQuant/isoquant.py $IQARGS -o $WORK/iq > $WORK/plain.log 2>&1
  V=${VERIFY_DIR:-$OUT/../_verify}/$SLICE/plain; mkdir -p $V
  cp $WORK/iq/OUT/*.gtf $WORK/iq/OUT/*.tsv $WORK/iq/OUT/*.tsv.gz $V/ 2>/dev/null || true
else
  /usr/bin/time -v -o $OUT/time.txt $PY -m iqcap.run --units $UNITS --out $OUT --sample $SAMPLE --slice $SLICE \
    --sample-rate "$RATE" ${BLOCK_RATE:+--block-rate $BLOCK_RATE} ${BLOCK_CAP:+--block-cap $BLOCK_CAP} --slice-def ${BAM%.bam}.slice.json -- $IQARGS -o $WORK/iq > $OUT/capture.log 2>&1
  if [ -n "${VERIFY_DIR:-}" ]; then
    V=$VERIFY_DIR/$SLICE/capture; mkdir -p $V
    cp $WORK/iq/OUT/*.gtf $WORK/iq/OUT/*.tsv $WORK/iq/OUT/*.tsv.gz $V/ 2>/dev/null || true
  fi
fi
rm -rf "$WORK"
echo "done $SAMPLE $SLICE"
