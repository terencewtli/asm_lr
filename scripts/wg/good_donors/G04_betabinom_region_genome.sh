#!/bin/bash
#$ -N G04_betabinom_region_genome
#$ -l h_data=6G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-132:1
#$ -tc 40
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/G04/G04_betabinom_region_genome.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/G04/G04_betabinom_region_genome.$JOB_ID.$TASK_ID
#$ -cwd

# Genome-wide (22 autosomes) region-level beta-binomial ASM calling for the
# 6 donors already confirmed sane (excludes HG00344/NA21144 pending their
# dispersion-anomaly root-cause -- see md notes). One task per (donor, chrom),
# txt/samples/six_donor_chrom_tasks.txt = 6 x 22 = 132 tasks.
# Timed on real data: chr1 ~14min/sample, chr15 ~5min/sample (roughly linear
# in chrom size, ~0.05-0.056 min/Mb) -- longest single task (chr1) should
# still finish inside the 1hr budget with room to spare.

set -euo pipefail
export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/whatshap-env/bin:$PATH"

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
TASKS=$PROJDIR/txt/samples/six_donor_chrom_tasks.txt
ID=$SGE_TASK_ID
LINE=$(sed -n "${ID}p" "$TASKS")
SAMPLE=$(echo "$LINE" | cut -f1)
CHROM=$(echo "$LINE" | cut -f2)

echo "$(date): G04 task $ID -- $SAMPLE $CHROM"
echo "Hostname: $(hostname)"

OUTFILE=$PROJDIR/tables/dasm_betabinom/$SAMPLE/regions_${CHROM}.tsv
if [ -f "$OUTFILE" ]; then
    echo "$(date): $OUTFILE exists, skipping"
    exit 0
fi

time python3 "$PROJDIR/scripts/wg/good_donors/G02_betabinom_region_chr1.py" "$SAMPLE" "$CHROM"

echo "$(date): G04 done -- $SAMPLE $CHROM"
