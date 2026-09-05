#!/bin/bash
#$ -N W01g_merge_haplotagged_wg
#$ -l h_data=2G,h_rt=4:00:00
#$ -pe shared 10
#$ -t 1-18:1
#$ -hold_jid W01_phase_per_chrom
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W01g/W01g_merge_haplotagged_wg.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W01g/W01g_merge_haplotagged_wg.$JOB_ID.$TASK_ID
#$ -cwd

# Merge per-chromosome haplotagged BAMs into a whole-genome BAM per sample.
# Input:  bam/haplotagged/by_chrom/{SAMPLE}_chr{N}_haplotagged.bam
# Output: bam/haplotagged/{SAMPLE}_wg_haplotagged.bam
#
# Generalized to all 18 samples via final_samples.txt (was previously a
# hardcoded 16-sample list that excluded NA18620/NA21093 pending a separate
# patch -- that patch is applied, all 18 now have complete by_chrom/ input).

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
SAMTOOLS=/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools
INDIR=$PROJDIR/bam/haplotagged/by_chrom
OUTDIR=$PROJDIR/bam/haplotagged

SAMPLES=($(cat $PROJDIR/txt/samples/final_samples.txt))
ID=$SGE_TASK_ID

SAMPLE=${SAMPLES[$((ID - 1))]}

OUT=$OUTDIR/${SAMPLE}_wg_haplotagged.bam

echo "$(date): W01g task $ID — $SAMPLE"
echo "Hostname: $(hostname)"

if [ -f "$OUT" ] && [ -s "$OUT" ]; then
    echo "$(date): $OUT exists, skipping"
    exit 0
fi

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
        chr20 chr21 chr22)

MISSING=0
for CHROM in "${CHROMS[@]}"; do
    f=$INDIR/${SAMPLE}_${CHROM}_haplotagged.bam
    if [ ! -f "$f" ] || [ ! -s "$f" ]; then
        echo "ERROR: missing or empty $f"
        MISSING=$((MISSING + 1))
    fi
done
if [ "$MISSING" -gt 0 ]; then
    echo "ERROR: $MISSING chromosome BAMs missing — aborting"
    exit 1
fi

BAM_LIST=()
for CHROM in "${CHROMS[@]}"; do
    BAM_LIST+=("$INDIR/${SAMPLE}_${CHROM}_haplotagged.bam")
done

echo "$(date): merging 22 chrom BAMs..."
time $SAMTOOLS merge -f -@ 10 "$OUT" "${BAM_LIST[@]}"
time $SAMTOOLS index "$OUT"

TOTAL=$($SAMTOOLS view -c "$OUT")
SIZE=$(du -sh "$OUT" | cut -f1)
echo "$(date): merged $TOTAL reads → $OUT ($SIZE)"

HP=$($SAMTOOLS view "$OUT" | head -1000 | grep -c "HP:i:" || true)
echo "$(date): HP tag check (first 1000 reads): $HP/1000 have HP:i:"
if [ "$HP" -eq 0 ]; then
    echo "ERROR: no HP tags in merged BAM"
    exit 1
fi

echo "$(date): W01g done — $SAMPLE"
