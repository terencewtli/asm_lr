#!/bin/bash
#$ -N G03_betabinom_region_chr15
#$ -l h_data=6G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-8:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/G03/G03_betabinom_region_chr15.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/G03/G03_betabinom_region_chr15.$JOB_ID.$TASK_ID
#$ -cwd

# Region-level beta-binomial ASM calling, chr15 (chosen for its density of
# imprinted loci -- SNRPN/PWS-AS region -- as a targeted pilot), across the
# same 8 good donors as G02_betabinom_region_chr1.sh. chr15 pileup input
# already exists for all 8 (W03 stage), so this needs no upstream regen.
# chr15 is ~1/3 the size of chr1 (~102Mb vs ~248Mb autosomal), so expect well
# under the ~16-17 min/sample timed for chr1.

set -euo pipefail
export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/whatshap-env/bin:$PATH"

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
SAMPLES=($(cat $PROJDIR/txt/samples/good_donors.txt))
ID=$SGE_TASK_ID
SAMPLE=${SAMPLES[$((ID - 1))]}

echo "$(date): G03 task $ID -- $SAMPLE chr15"
echo "Hostname: $(hostname)"

OUTFILE=$PROJDIR/tables/dasm_betabinom/$SAMPLE/regions_chr15.tsv
if [ -f "$OUTFILE" ]; then
    echo "$(date): $OUTFILE exists, skipping"
    exit 0
fi

time python3 "$PROJDIR/scripts/wg/good_donors/G02_betabinom_region_chr1.py" "$SAMPLE" chr15

echo "$(date): G03 done -- $SAMPLE chr15"
