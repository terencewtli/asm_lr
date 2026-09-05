#!/bin/bash
#$ -N W00a_align_ubam
#$ -l h_data=2G,h_rt=24:00:00
#$ -pe shared 10
#$ -t 1-4:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W00a_align_ubam.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W00a_align_ubam.$JOB_ID.$TASK_ID
#$ -cwd

# Align one uBAM to GRCh38 and save sorted BAM to scratch.
# Array: one task per uBAM (72 total, see csv/meta/ubam_manifest.csv).
# Run W00b_merge_wg.sh after all tasks for a sample complete.

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
MANIFEST=$PROJDIR/csv/meta/ubam_manifest_NA21093.csv
SCRATCH=/u/project/cluo_scratch/terencew/claude/asm_lr/wg_align
REF=/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa

MINIMAP2=/u/home/t/terencew/bin/minimap2
SAMTOOLS=/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools
THREADS=10

mkdir -p "$PROJDIR/logs" "$SCRATCH"

# Read task row from manifest (skip header, 1-based)
MANIFEST_ROW=$(awk -F, -v t="$SGE_TASK_ID" 'NR==t+1 {print}' "$MANIFEST")
SAMPLE=$(echo "$MANIFEST_ROW" | cut -d, -f2)
UBAM=$(echo "$MANIFEST_ROW" | cut -d, -f3)
BNAME=$(basename "$UBAM" .bam)

echo "$(date): W00a task $SGE_TASK_ID — $SAMPLE / $BNAME"
echo "Hostname: $(hostname)"

if [ -z "$SAMPLE" ] || [ -z "$UBAM" ]; then
    echo "ERROR: could not parse manifest row for task $SGE_TASK_ID"
    exit 1
fi

ALIGNED=$SCRATCH/${SAMPLE}_${BNAME}_wg.bam

# Skip if already done
if [ -f "$ALIGNED" ] && [ -s "$ALIGNED" ]; then
    echo "$(date): $ALIGNED already exists — skipping"
    N=$($SAMTOOLS view -c "$ALIGNED")
    echo "$(date): existing BAM has $N reads"
    exit 0
fi

if [ ! -f "$UBAM" ]; then
    echo "ERROR: $UBAM not found"
    exit 1
fi

echo "$(date): aligning $(du -sh $UBAM | cut -f1) uBAM..."
time $SAMTOOLS bam2fq -T MM,ML "$UBAM" \
  | $MINIMAP2 -y -a -x map-ont -t $THREADS "$REF" - \
  | $SAMTOOLS view -bh - \
  | $SAMTOOLS sort -T "${SCRATCH}/sort_${BNAME}" -m 1G -@ $THREADS \
                   -o "$ALIGNED"
$SAMTOOLS index "$ALIGNED"

N=$($SAMTOOLS view -c "$ALIGNED")
echo "$(date): aligned $N reads → $ALIGNED ($(du -sh $ALIGNED | cut -f1))"

# Verify MM tags survived alignment
MM=$($SAMTOOLS view "$ALIGNED" | head -500 | grep -c "MM:Z:" || true)
echo "$(date): MM tag check (first 500 reads): $MM have MM:Z:"
if [ "$MM" -eq 0 ]; then
    echo "ERROR: no MM tags in output — check -y flag and MM,ML bam2fq"
    exit 1
fi

echo "$(date): W00a done — run W00b_merge_wg.sh after all ubams for $SAMPLE complete"
