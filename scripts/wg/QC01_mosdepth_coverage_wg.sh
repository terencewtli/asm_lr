#!/bin/bash
#$ -N QC01_mosdepth_coverage_wg
#$ -l h_data=2G,h_rt=4:00:00
#$ -pe shared 4
#$ -t 1-18:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/QC01/QC01_mosdepth_coverage_wg.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/QC01/QC01_mosdepth_coverage_wg.$JOB_ID.$TASK_ID
#$ -cwd

# Real per-sample coverage from the actual haplotagged wg BAMs, to replace
# the stale pre-alignment planning estimates in notebooks/qc/csv/qc_cohort_metadata.csv
# (ont_total_cov_x -- dated 2026-06-13, before these BAMs existed). Runs
# independently of W02/W03/W04 -- only needs the haplotagged BAM, which
# exists for all 18 samples already.
#
# -n/--no-per-base skips base-level output (we only want the summary), -x
# fast-mode skips CIGAR-based mate-overlap correction (recommended for long
# reads anyway). mosdepth is lightweight -- no HP-tag splitting here, this
# is total coverage only; per-haplotype coverage would need a separate,
# more expensive pass filtering by the HP tag.

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate whatshap-env

MOSDEPTH=/u/project/cluo/terencew/programs/mosdepth/v0.3.14/mosdepth

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
OUTDIR=$PROJDIR/qc/mosdepth

SAMPLES=($(cat $PROJDIR/txt/samples/final_samples.txt))
ID=$SGE_TASK_ID
SAMPLE=${SAMPLES[$((ID - 1))]}

BAM=$PROJDIR/bam/haplotagged/${SAMPLE}_wg_haplotagged.bam

mkdir -p "$OUTDIR"

echo "$(date): QC01 task $ID — $SAMPLE"

if [ ! -f "$BAM" ]; then
    echo "ERROR: $BAM not found"
    exit 1
fi

if [ -f "$OUTDIR/${SAMPLE}.mosdepth.summary.txt" ]; then
    echo "$(date): ${SAMPLE}.mosdepth.summary.txt exists, skipping"
    exit 0
fi

time $MOSDEPTH --no-per-base --fast-mode --threads 4 "$OUTDIR/${SAMPLE}" "$BAM"

echo "$(date): $SAMPLE done"
cat "$OUTDIR/${SAMPLE}.mosdepth.summary.txt"
