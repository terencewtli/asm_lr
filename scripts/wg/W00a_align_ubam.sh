#!/bin/bash
#$ -N W00a_align_ubam
#$ -l h_data=3G,h_rt=16:00:00
#$ -pe shared 8
#$ -t 1-72:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W00a/W00a_align_ubam.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W00a/W00a_align_ubam.$JOB_ID.$TASK_ID
#$ -cwd

# Align one uBAM to GRCh38 and save sorted BAM to scratch.
# Array: one task per uBAM (72 total, see csv/meta/ubam_manifest.csv).
# Run W00b_merge_wg.sh after all tasks for a sample complete.
#
# h_data bumped 2G->3G on 2026-09-07: several uBAM tasks (across NA20762,
# HG00253, NA19776, NA20870) were repeatedly silently OOM-killed ("Killed"
# in the log, no further output) under the old 2G x 8 = 16GB total budget --
# one surviving task's log showed minimap2 alone hitting Peak RSS 16.57GB,
# i.e. at/over the old ceiling before samtools sort/view even ran alongside
# it. 3G x 8 = 24GB gives real headroom. This was a plain resource-sizing
# bug, not data-quality-specific -- it's plausible it also silently affected
# some of the original 12 unaffected samples' tasks in the very first run
# and just wasn't caught because those samples completed after retries.

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
MANIFEST=$PROJDIR/csv/meta/ubam_manifest.csv
SCRATCH=/u/project/cluo_scratch/terencew/claude/asm_lr/wg_align
REF=/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa

MINIMAP2=/u/home/t/terencew/bin/minimap2
SAMTOOLS=/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools
THREADS=8

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
DONE_MARKER=${ALIGNED}.done

# Skip only if a previous run both finished AND passed the completeness check
# below (marked by DONE_MARKER). A file that merely exists (e.g. left behind
# by a job killed mid-alignment by walltime/OOM) is NOT trusted -- this was
# the root cause of several silently-truncated samples in the 18-donor cohort
# (see md/20260905.progress.md): the old skip-check only tested non-empty,
# so a truncated BAM from a killed job was skipped forever on every retry.
if [ -f "$ALIGNED" ] && [ -f "$DONE_MARKER" ]; then
    echo "$(date): $ALIGNED already verified complete (found $DONE_MARKER) — skipping"
    exit 0
fi

if [ ! -f "$UBAM" ]; then
    echo "ERROR: $UBAM not found"
    exit 1
fi

echo "$(date): counting source reads in $UBAM (for completeness check)..."
time SRC_READS=$($SAMTOOLS view -c "$UBAM")
echo "$(date): source uBAM has $SRC_READS reads"

echo "$(date): aligning $(du -sh $UBAM | cut -f1) uBAM..."
rm -f "$DONE_MARKER"
time $SAMTOOLS bam2fq -T MM,ML "$UBAM" \
  | $MINIMAP2 -y -a -x map-ont -t $THREADS "$REF" - \
  | $SAMTOOLS view -bh - \
  | $SAMTOOLS sort -T "${SCRATCH}/sort_${BNAME}" -m 1G -@ $THREADS \
                   -o "$ALIGNED"
time $SAMTOOLS index "$ALIGNED"

N=$($SAMTOOLS view -c "$ALIGNED")
echo "$(date): aligned $N reads → $ALIGNED ($(du -sh $ALIGNED | cut -f1))"

# Completeness check: aligned read count should be close to source read count
# (some loss from unmapped/filtered reads is normal, but a large shortfall
# means the job was killed mid-pipe or the source was truncated).
PCT=$(awk -v n="$N" -v s="$SRC_READS" 'BEGIN { if (s==0) print 0; else printf "%.1f", n/s*100 }')
echo "$(date): aligned/source read ratio: ${PCT}%"
if awk -v p="$PCT" 'BEGIN { exit !(p < 90) }'; then
    echo "ERROR: only ${PCT}% of source reads made it into the aligned BAM ($N / $SRC_READS)."
    echo "  This BAM is likely truncated (killed mid-alignment, or source uBAM itself"
    echo "  truncated). NOT marking as done -- rerun this task; do not resubmit blindly"
    echo "  without checking why it's short (walltime? scratch disk full? source corrupt?)."
    exit 1
fi

# Verify MM tags survived alignment
MM=$($SAMTOOLS view "$ALIGNED" | head -500 | grep -c "MM:Z:" || true)
echo "$(date): MM tag check (first 500 reads): $MM have MM:Z:"
if [ "$MM" -eq 0 ]; then
    echo "ERROR: no MM tags in output — check -y flag and MM,ML bam2fq"
    exit 1
fi

touch "$DONE_MARKER"
echo "$(date): W00a done — run W00b_merge_wg.sh after all ubams for $SAMPLE complete"
