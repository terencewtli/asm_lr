#!/bin/bash
#$ -N W04_dasm_loci_wg
#$ -l h_data=16G,h_rt=5:00:00
#$ -pe shared 4
#$ -t 1-22
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W04_dasm_loci_wg.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W04_dasm_loci_wg.$JOB_ID.$TASK_ID
#$ -hold_jid W01_modkit_extract_wg,W03_modkit_pileup_wg
#$ -cwd

# dASM loci — job array per chromosome.
# Calls W04_dasm_loci.py (delta threshold, then Mann-Whitney/AUC at locus level).
# h_data=16G needed to load the full-genome modkit extract TSV per chromosome.
# Runs after W01 (extract) and W03 (pileup), which both depend on W02 (phase).

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate allcools

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
CHROM=chr${SGE_TASK_ID}

echo "$(date): W04 dASM loci — $CHROM"
python $PROJDIR/scripts/wg/W04_dasm_loci.py \
    --chrom  $CHROM \
    --sample NA21093 \
    --projdir $PROJDIR

echo "$(date): $CHROM done"
