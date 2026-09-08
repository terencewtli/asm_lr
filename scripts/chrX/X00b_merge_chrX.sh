#!/bin/bash
#$ -N X00b_merge_chrX
#$ -l h_data=2G,h_rt=8:00:00
#$ -pe shared 4
#$ -t 1-18:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/X00b/X00b_merge_chrX.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/X00b/X00b_merge_chrX.$JOB_ID.$TASK_ID
#$ -cwd

# chrX analog of W00b_merge_wg.sh: merge per-uBAM chrX-aligned scratch BAMs
# for one sample into bam/chrX/merged_bam/{SAMPLE}_chrX.bam. Same array
# indexing (final_samples.txt order) as the autosome pipeline, so this can
# run for all 18 even though only the 8 female donors have a real
# haplotype-comparison use for chrX -- male donors' merged chrX BAM is still
# useful as the negative-control / bulk-methylation reference (see
# scripts/chrX/README.md).

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
SCRATCH=/u/project/cluo_scratch/terencew/claude/asm_lr/chrX_align
TMPDIR_MERGE=$PROJDIR/tmp/bam

SAMTOOLS=/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools
THREADS=10

SAMPLES=($(cat $PROJDIR/txt/samples/final_samples.txt))
SAMPLE=${SAMPLES[$((SGE_TASK_ID - 1))]}

echo "$(date): X00b task $SGE_TASK_ID — $SAMPLE"
echo "Hostname: $(hostname)"

OUT_BAM=$PROJDIR/bam/chrX/merged_bam/${SAMPLE}_chrX.bam
mkdir -p "$(dirname "$OUT_BAM")"

if [ -f "$OUT_BAM" ] && [ -s "$OUT_BAM" ]; then
    echo "$(date): $OUT_BAM exists — skipping"
    exit 0
fi

ALIGNED_BAMS=($(ls $SCRATCH/${SAMPLE}_*_chrX.bam 2>/dev/null))
N_RUNS=${#ALIGNED_BAMS[@]}

if [ $N_RUNS -eq 0 ]; then
    echo "ERROR: no scratch chrX BAMs found for $SAMPLE — run X00a first"
    exit 1
fi

EXPECTED=$(awk -F, -v s="$SAMPLE" 'NR>1 && $2==s {count++} END {print count}' \
    $PROJDIR/csv/meta/ubam_manifest.csv)
if [ "$N_RUNS" -ne "$EXPECTED" ]; then
    echo "ERROR: found $N_RUNS scratch chrX BAMs but expected $EXPECTED for $SAMPLE"
    exit 1
fi

echo "$(date): merging $N_RUNS chrX BAMs for $SAMPLE"
mkdir -p "$TMPDIR_MERGE"
MERGE_TMP=$TMPDIR_MERGE/${SAMPLE}_chrX_merged_tmp.bam

if [ $N_RUNS -eq 1 ]; then
    cp "${ALIGNED_BAMS[0]}" "$MERGE_TMP"
else
    time $SAMTOOLS merge -f -@ $THREADS "$MERGE_TMP" "${ALIGNED_BAMS[@]}"
fi
time $SAMTOOLS index "$MERGE_TMP"

TOTAL=$($SAMTOOLS view -c "$MERGE_TMP")
echo "$(date): merged total = $TOTAL reads ($(du -sh $MERGE_TMP | cut -f1))"

mv "$MERGE_TMP"       "$OUT_BAM"
mv "${MERGE_TMP}.bai" "${OUT_BAM}.bai" 2>/dev/null || time $SAMTOOLS index "$OUT_BAM"

for B in "${ALIGNED_BAMS[@]}"; do
    rm -f "$B" "${B}.bai" "${B}.done"
done

echo "$(date): X00b done — $OUT_BAM ($TOTAL reads)"
