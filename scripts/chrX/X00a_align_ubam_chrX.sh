#!/bin/bash
#$ -N X00a_align_ubam_chrX
#$ -l h_data=3G,h_rt=16:00:00
#$ -pe shared 8
#$ -t 1-72:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/X00a/X00a_align_ubam_chrX.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/X00a/X00a_align_ubam_chrX.$JOB_ID.$TASK_ID
#$ -cwd

# chrX analog of scripts/wg/W00a_align_ubam.sh: aligns the SAME raw uBAMs
# (csv/meta/ubam_manifest.csv, one task per uBAM, same 72-row manifest --
# chrX reads are a subset of the same raw data, not a separate download)
# against a chrX-ONLY reference instead of the autosome-only one. Aligning
# against a single small chromosome means most reads (everything from other
# chromosomes) simply fail to place -- wasteful but correct, and much faster
# per-uBAM than the full-genome alignment since there's no autosomal target
# to search.
#
# Reference: reference/chrX_ref/GRCh38.chrX.fa (extracted 2026-09-07 via
# `samtools faidx .../GRCh38.fa chrX` from the full-genome reference at
# /u/project/cluo/terencew/reference/hg38_igvf/GRCh38.fa -- NOT the
# project's usual GRCh38.autosome.fa, which has no sex chromosomes at all).
# Deliberately chrX-only, not chrX+chrY: since chrY isn't present, PAR1/PAR2
# reads have nowhere else to map and simply place uniquely on chrX -- no
# PAR-masking logic needed, unlike a reference that includes both.
#
# h_data bumped straight to 3G (matching the 2026-09-07 fix to the autosome
# W00a) to avoid the same OOM failure mode from the start.

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
MANIFEST=$PROJDIR/csv/meta/ubam_manifest.csv
SCRATCH=/u/project/cluo_scratch/terencew/claude/asm_lr/chrX_align
REF=$PROJDIR/reference/chrX_ref/GRCh38.chrX.fa

MINIMAP2=/u/home/t/terencew/bin/minimap2
SAMTOOLS=/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools
THREADS=8

mkdir -p "$PROJDIR/logs/X00a" "$SCRATCH"

MANIFEST_ROW=$(awk -F, -v t="$SGE_TASK_ID" 'NR==t+1 {print}' "$MANIFEST")
SAMPLE=$(echo "$MANIFEST_ROW" | cut -d, -f2)
UBAM=$(echo "$MANIFEST_ROW" | cut -d, -f3)
BNAME=$(basename "$UBAM" .bam)

echo "$(date): X00a task $SGE_TASK_ID — $SAMPLE / $BNAME"
echo "Hostname: $(hostname)"

if [ -z "$SAMPLE" ] || [ -z "$UBAM" ]; then
    echo "ERROR: could not parse manifest row for task $SGE_TASK_ID"
    exit 1
fi

ALIGNED=$SCRATCH/${SAMPLE}_${BNAME}_chrX.bam
DONE_MARKER=${ALIGNED}.done

if [ -f "$ALIGNED" ] && [ -f "$DONE_MARKER" ]; then
    echo "$(date): $ALIGNED already verified complete — skipping"
    exit 0
fi

if [ ! -f "$UBAM" ]; then
    echo "ERROR: $UBAM not found"
    exit 1
fi

echo "$(date): counting source reads in $UBAM..."
time SRC_READS=$($SAMTOOLS view -c "$UBAM")
echo "$(date): source uBAM has $SRC_READS reads (genome-wide; only a chrX-sized fraction will align here)"

echo "$(date): aligning $(du -sh $UBAM | cut -f1) uBAM against chrX-only reference..."
rm -f "$DONE_MARKER"
time $SAMTOOLS bam2fq -T MM,ML "$UBAM" \
  | $MINIMAP2 -y -a -x map-ont -t $THREADS "$REF" - \
  | $SAMTOOLS view -bh - \
  | $SAMTOOLS sort -T "${SCRATCH}/sort_${BNAME}" -m 1G -@ $THREADS \
                   -o "$ALIGNED"
time $SAMTOOLS index "$ALIGNED"

N=$($SAMTOOLS view -c "$ALIGNED")
echo "$(date): aligned $N reads to chrX -> $ALIGNED ($(du -sh $ALIGNED | cut -f1))"
# NOTE: no aligned/source ratio completeness check here, unlike the autosome
# W00a -- a low ratio is EXPECTED (most reads are from other chromosomes,
# not truncated data), so that check would misfire. If N is suspiciously
# close to 0 for every task of a sample, that's worth a manual look, but
# there's no principled per-task threshold the way there is for a
# whole-genome alignment.

MM=$($SAMTOOLS view "$ALIGNED" | head -500 | grep -c "MM:Z:" || true)
echo "$(date): MM tag check (first 500 aligned reads): $MM have MM:Z:"

touch "$DONE_MARKER"
echo "$(date): X00a done — run X00b_merge_chrX.sh after all ubams for $SAMPLE complete"
