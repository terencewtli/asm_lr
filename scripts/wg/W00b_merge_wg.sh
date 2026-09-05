#!/bin/bash
#$ -N W00b_merge_wg
#$ -l h_data=2G,h_rt=8:00:00
#$ -pe shared 4
#$ -t 1-18:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W00b/W00b_merge_wg.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W00b/W00b_merge_wg.$JOB_ID.$TASK_ID
#$ -cwd

# Merge per-ubam scratch BAMs for one sample into the final bam/{SAMPLE}_wg.bam.
# Run after ALL W00a tasks for the sample have completed successfully.
# Array: 1-18 (one per sample, matching order in txt/samples/final_samples.txt).

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
SCRATCH=/u/project/cluo_scratch/terencew/claude/asm_lr/wg_align
TMPDIR_MERGE=$PROJDIR/tmp/bam

SAMTOOLS=/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools
THREADS=10

SAMPLES=($(cat $PROJDIR/txt/samples/final_samples.txt))
SAMPLE=${SAMPLES[$((SGE_TASK_ID - 1))]}

echo "$(date): W00b task $SGE_TASK_ID — $SAMPLE"
echo "Hostname: $(hostname)"

OUT_BAM=$PROJDIR/bam/merged_bam/${SAMPLE}_wg.bam

# Skip if already done
if [ -f "$OUT_BAM" ] && [ -s "$OUT_BAM" ]; then
    echo "$(date): $OUT_BAM exists — skipping"
    exit 0
fi

# Collect per-ubam scratch BAMs for this sample
ALIGNED_BAMS=($(ls $SCRATCH/${SAMPLE}_*_wg.bam 2>/dev/null))
N_RUNS=${#ALIGNED_BAMS[@]}

if [ $N_RUNS -eq 0 ]; then
    echo "ERROR: no scratch BAMs found for $SAMPLE — run W00a first"
    exit 1
fi

# Verify all expected ubams are present
EXPECTED=$(awk -F, -v s="$SAMPLE" 'NR>1 && $2==s {count++} END {print count}' \
    $PROJDIR/csv/meta/ubam_manifest.csv)
if [ "$N_RUNS" -ne "$EXPECTED" ]; then
    echo "ERROR: found $N_RUNS scratch BAMs but expected $EXPECTED for $SAMPLE"
    echo "  found: ${ALIGNED_BAMS[*]}"
    exit 1
fi

echo "$(date): merging $N_RUNS BAMs for $SAMPLE"
for B in "${ALIGNED_BAMS[@]}"; do
    echo "  $(du -sh $B | cut -f1)  $B"
done

mkdir -p "$TMPDIR_MERGE"
MERGE_TMP=$TMPDIR_MERGE/${SAMPLE}_wg_merged_tmp.bam

if [ $N_RUNS -eq 1 ]; then
    cp "${ALIGNED_BAMS[0]}" "$MERGE_TMP"
else
    time $SAMTOOLS merge -f -@ $THREADS "$MERGE_TMP" "${ALIGNED_BAMS[@]}"
    time $SAMTOOLS index "$MERGE_TMP"
fi

TOTAL=$($SAMTOOLS view -c "$MERGE_TMP")
echo "$(date): merged total = $TOTAL reads ($(du -sh $MERGE_TMP | cut -f1))"

MM=$($SAMTOOLS view "$MERGE_TMP" | head -500 | grep -c "MM:Z:" || true)
echo "$(date): MM tag check (first 500 reads): $MM have MM:Z:"
if [ "$MM" -eq 0 ]; then
    echo "ERROR: no MM tags in merged BAM"
    exit 1
fi

mv "$MERGE_TMP"       "$OUT_BAM"
mv "${MERGE_TMP}.bai" "${OUT_BAM}.bai" 2>/dev/null || time $SAMTOOLS index "$OUT_BAM"

# Clean up per-ubam scratch BAMs
for B in "${ALIGNED_BAMS[@]}"; do
    rm -f "$B" "${B}.bai"
done

echo "$(date): W00b done"
echo "  BAM:   $OUT_BAM  ($(du -sh $OUT_BAM | cut -f1))"
echo "  Reads: $TOTAL"
echo "  Next:  qsub scripts/wg/W01_phase_wg.sh"
