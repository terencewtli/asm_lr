#!/bin/bash
#$ -N W03_modkit_pileup_wg
#$ -l h_data=8G,h_rt=1:00:00
#$ -pe shared 10
#$ -t 1-22
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W03_modkit_pileup_wg.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W03_modkit_pileup_wg.$JOB_ID.$TASK_ID
#$ -hold_jid W02_phase_wg
#$ -cwd

# HP-stratified modkit pileup — job array, one element per autosome (~20 min each).
# Requires: haplotagged BAMs from W02. Runs after W02_phase_wg.
# Output: pileup/{sample}_{chrom}_hp_pileup_{1,2,ungrouped}.bed.gz (tabix-indexed)

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate whatshap-env

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
SAMPLE=NA21093
CHROM=chr${SGE_TASK_ID}
BAM=$PROJDIR/bam/${SAMPLE}_${CHROM}_haplotagged.bam
REF=/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa
OUTDIR=$PROJDIR/pileup
LOGFILE=$PROJDIR/logs/modkit_pileup_${SAMPLE}_${CHROM}.$JOB_ID.log

mkdir -p $OUTDIR

echo "$(date): W03 modkit pileup — $SAMPLE $CHROM"
modkit pileup \
    --ref $REF \
    --region $CHROM \
    --partition-tag HP \
    --prefix ${OUTDIR}/${SAMPLE}_${CHROM}_hp_pileup \
    --threads 10 \
    --log-filepath $LOGFILE \
    $BAM \
    ${OUTDIR}/${SAMPLE}_${CHROM}_hp_pileup

for f in ${OUTDIR}/${SAMPLE}_${CHROM}_hp_pileup_*.bed.gz; do
    tabix -p bed $f && echo "indexed $f"
done

echo "$(date): $CHROM pileup done"
