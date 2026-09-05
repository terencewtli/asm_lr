#!/bin/bash
#$ -N W03b_merge_haplotagged_wg
#$ -l h_data=2G,h_rt=2:00:00
#$ -pe shared 10
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W03b_merge_haplotagged_wg.$JOB_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W03b_merge_haplotagged_wg.$JOB_ID
#$ -hold_jid W03_phase_wg
#$ -cwd

# Merge per-chromosome haplotagged BAMs (output of W03) into a single whole-genome BAM.
# Required by W02_modkit_pileup_wg.sh which expects bam/NA21093_wg_haplotagged.bam.

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
SAMPLE=NA21093
SAMTOOLS=/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools
BAMDIR=$PROJDIR/bam
OUT=$BAMDIR/${SAMPLE}_wg_haplotagged.bam

echo "$(date): W03b — merging per-chromosome haplotagged BAMs"

# verify all 22 per-chrom BAMs exist before merging
MISSING=0
for i in $(seq 1 22); do
    f=$BAMDIR/${SAMPLE}_chr${i}_haplotagged.bam
    if [ ! -f "$f" ]; then
        echo "ERROR: missing $f"
        MISSING=$((MISSING + 1))
    fi
done
if [ "$MISSING" -gt 0 ]; then
    echo "ERROR: $MISSING chromosome BAMs missing — recheck W03 job array"
    exit 1
fi

$SAMTOOLS merge -f -@ 10 "$OUT" \
    $(for i in $(seq 1 22); do echo $BAMDIR/${SAMPLE}_chr${i}_haplotagged.bam; done)

$SAMTOOLS index "$OUT"

TOTAL=$($SAMTOOLS view -c "$OUT")
echo "$(date): merged $TOTAL reads → $OUT  ($(du -sh $OUT | cut -f1))"

# verify HP tags present
HP=$($SAMTOOLS view "$OUT" | head -1000 | grep -c "HP:i:" || true)
echo "$(date): HP tag check (first 1000 reads): $HP have HP:i:"
if [ "$HP" -eq 0 ]; then
    echo "ERROR: no HP tags found in merged BAM"
    exit 1
fi

echo "$(date): W03b done"
echo "Next: qsub scripts/qsub/wg/W02_modkit_pileup_wg.sh"
