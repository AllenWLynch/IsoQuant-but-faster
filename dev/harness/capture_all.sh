#!/bin/bash
# The M0a fixture set. Usage: capture_all.sh local|sbatch [filter-regex]
# local: runs entries tagged L on this node (STAGE must point at a local copy of the refs), N in parallel.
# sbatch: submits entries tagged S as batch jobs (each job stages refs to $TMPDIR; extra options in $SBATCH_OPTS).
set -uo pipefail
MODE=$1; FILT=${2:-.}
H=$(cd "$(dirname "$0")" && pwd)
export IQRS_PROJECT=${IQRS_PROJECT:-$(cd "$H/../../.." && pwd)}
SMALL="1.0"
MID="U4=0.25,U5=0.1,U6=0.5,U7=0.25,U20=0.25,default=1"
DEEP="U4=0.05,U5=0.02,U6=0.25,U7=0.05,U20=0.1,default=1"
FULL="U20=0.05,default=1"
# where sample slice regions dtype rates subsample units
TABLE=$(cat <<'T'
L SRR33046254 ordinary 8:8699085-8704430,8:14170567-14183265,8:16758304-16769798,8:22521729-22542545,8:51350453-51367298 nanopore SMALL 1.0 all
L SRR33046254 deepest 8:30699499-30827785 nanopore MID 1.0 all
L SRR33046254 multigene 8:23678436-23882242 nanopore SMALL 1.0 all
L SRR33046254 nogene 8:40502837-40522974,8:50137473-50139985,8:7644321-7658264 nanopore SMALL 1.0 all
L SRR33046254 valley 8:21104496-21142297 nanopore SMALL 1.0 all
L SRR32588720 ordinary 8:7205524-7308373,8:11884671-11988065,8:14058185-14080534,8:23355484-23381548,8:50742975-50839564 pacbio_ccs SMALL 1.0 all
S SRR32588720 deepest 8:47098781-47214606 pacbio_ccs MID 1.0 all
L SRR32588720 multigene 8:20832584-21179323 pacbio_ccs SMALL 1.0 all
L SRR32588720 nogene 8:40502753-40526405,8:28384127-28391287,8:45504296-45506113 pacbio_ccs SMALL 1.0 all
L SRR32588720 valley 8:15269280-15305546 pacbio_ccs SMALL 1.0 all
S SRR32588726 ordinary 8:10304981-10308771,8:16615798-16675464,8:18010568-18036639,8:25034529-25040762,8:52314542-52350684 pacbio_ccs MID 1.0 all
S SRR32588726 deepest 8:13749762-13972736 pacbio_ccs DEEP 0.15 all
S SRR32588726 multigene 8:16469721-16769809 pacbio_ccs MID 1.0 all
L SRR32588726 nogene 8:28551295-28574739,8:34163831-34166057,8:34080370-34134026 pacbio_ccs SMALL 1.0 all
S SRR32588726 valley 8:18449060-18829702 pacbio_ccs MID 1.0 all
S SRR33046254 fullchr8 8:1-54304671 nanopore FULL 1.0 U8,U20
S SRR32588720 fullchr8 8:1-54304671 pacbio_ccs FULL 1.0 U8,U20
T
)
N=${N:-6}
echo "$TABLE" | grep -E "$FILT" | { while read W S SL R D RT SUB U; do
  RATE=${!RT}
  if [ "$MODE" = local ] && [ "$W" = L ]; then
    while [ $(jobs -r | wc -l) -ge $N ]; do sleep 5; done
    echo "local $S $SL"; $H/capture_slice.sh $S $SL $R $D "$RATE" $SUB $U > /dev/null 2>&1 &
  elif [ "$MODE" = sbatch ] && [ "$W" = S ]; then
    mkdir -p $H/logs
    sbatch ${SBATCH_OPTS:--t 4:00:00 -c 2 --mem=16G} -J iqcap_${S}_${SL} -o $H/logs/${S}_${SL}.%j.out \
      $H/capture_slice.sh $S $SL $R $D "$RATE" $SUB $U
  fi
done; wait; }
