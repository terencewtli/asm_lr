#!/bin/bash
#$ -N W00_align_wg
#$ -l h_data=2G,h_rt=24:00:00
#$ -pe shared 16
#$ -t 1-18:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W00_align_wg.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W00_align_wg.$JOB_ID.$TASK_ID
#$ -cwd

# Merge + align ONT uBAMs for one cohort sample to GRCh38 (whole genome).
# Reads pre-downloaded uBAMs from bam/ubam/{SAMPLE}/, aligns each run to scratch,
# merges into a single sorted BAM at bam/{SAMPLE}_wg.bam.
#
# Prerequisites: D01_download_ont_ubams.sh must have completed for this sample.
# Output: bam/{SAMPLE}_wg.bam  (+.bai)  — input for W02_phase_wg.sh
# Skips: if bam/{SAMPLE}_wg.bam already exists and is non-empty.
#
# MM/ML methylation tags are preserved through alignment via:
#   samtools bam2fq -T MM,ML | minimap2 -y  (-y passes query tags to output SAM)

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
MANIFEST=$PROJDIR/csv/meta/cohort_downloads.tsv
SCRATCH=/u/project/cluo_scratch/terencew/claude/asm_lr/wg_align
REF=/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa

MINIMAP2=/u/home/t/terencew/bin/minimap2
SAMTOOLS=/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools
THREADS=10

mkdir -p "$PROJDIR/logs" "$SCRATCH"

# Map task ID → sample (18 Dorado samples only; see txt/samples/final_samples.txt)
SAMPLES=($(cat $PROJDIR/txt/samples/final_samples.txt))
SAMPLE=${SAMPLES[$((SGE_TASK_ID - 1))]}

UBAM_DIR=$PROJDIR/bam/ubam/$SAMPLE
OUT_BAM=$PROJDIR/bam/${SAMPLE}_wg.bam

echo "$(date): W00 task $SGE_TASK_ID — $SAMPLE"
echo "Hostname: $(hostname)"

# Skip if already done
if [ -f "$OUT_BAM" ] && [ -s "$OUT_BAM" ]; then
    echo "$(date): $OUT_BAM exists and non-empty — skipping"
    exit 0
fi

# Verify uBAMs are present
if [ ! -d "$UBAM_DIR" ]; then
    echo "ERROR: $UBAM_DIR not found — run D01_download_ont_ubams.sh first"
    exit 1
fi

UBAMS=($UBAM_DIR/*.bam)
N_RUNS=${#UBAMS[@]}
if [ $N_RUNS -eq 0 ]; then
    echo "ERROR: no .bam files found in $UBAM_DIR"
    exit 1
fi
echo "$(date): found $N_RUNS uBAM file(s) in $UBAM_DIR"

# Align each run to scratch
ALIGNED_BAMS=()
for UBAM in "${UBAMS[@]}"; do
    BNAME=$(basename "$UBAM" .bam)
    ALIGNED=$SCRATCH/${SAMPLE}_${BNAME}_wg.bam

    echo "$(date): aligning $BNAME ($(du -sh $UBAM | cut -f1))..."
    time $SAMTOOLS bam2fq -T MM,ML "$UBAM" \
      | $MINIMAP2 -y -a -x map-ont -t $THREADS "$REF" - \
      | $SAMTOOLS view -bh - \
      | $SAMTOOLS sort -T "${SCRATCH}/sort_${SAMPLE}_${BNAME}" -m 1G -@ $THREADS \
                       -o "$ALIGNED"
    $SAMTOOLS index "$ALIGNED"

    N=$($SAMTOOLS view -c "$ALIGNED")
    echo "$(date): $BNAME → $N reads aligned"
    ALIGNED_BAMS+=("$ALIGNED")
done

# Merge (or rename if only one run)
MERGE_TMP=$SCRATCH/${SAMPLE}_wg_merged_tmp.bam
if [ $N_RUNS -eq 1 ]; then
    mv "${ALIGNED_BAMS[0]}" "$MERGE_TMP"
    mv "${ALIGNED_BAMS[0]}.bai" "${MERGE_TMP}.bai"
else
    echo "$(date): merging $N_RUNS runs..."
    time $SAMTOOLS merge -f -@ $THREADS "$MERGE_TMP" "${ALIGNED_BAMS[@]}"
    $SAMTOOLS index "$MERGE_TMP"
fi

TOTAL=$($SAMTOOLS view -c "$MERGE_TMP")
echo "$(date): merged total = $TOTAL reads ($(du -sh $MERGE_TMP | cut -f1))"

# Verify MM tags survived alignment
MM=$($SAMTOOLS view "$MERGE_TMP" | head -500 | grep -c "MM:Z:" || true)
echo "$(date): MM tag check (first 500 reads): $MM have MM:Z:"
if [ "$MM" -eq 0 ]; then
    echo "ERROR: no MM tags in output — methylation will be missing downstream"
    exit 1
fi

# Move to final location
mv "$MERGE_TMP"       "$OUT_BAM"
mv "${MERGE_TMP}.bai" "${OUT_BAM}.bai"

# Clean up per-run scratch BAMs
for ALIGNED in "${ALIGNED_BAMS[@]}"; do
    rm -f "$ALIGNED" "${ALIGNED}.bai"
done

echo "$(date): W00 done"
echo "  BAM:   $OUT_BAM  ($(du -sh $OUT_BAM | cut -f1))"
echo "  Reads: $TOTAL"
echo "  Next:  qsub scripts/wg/W02_phase_wg.sh (requires dipcall VCF or 1000G VCF)"
