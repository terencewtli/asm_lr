#!/bin/bash
#$ -N W07_classify_annotate_wg
#$ -l h_data=8G,h_rt=1:00:00
#$ -pe shared 2
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W07_classify_annotate_wg.$JOB_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W07_classify_annotate_wg.$JOB_ID
#$ -hold_jid W05_merge_wg
#$ -cwd

# Classify dASM loci by SNP distance, annotate with CpG context and ChromHMM.
# Reads dasm_loci_wg.parquet, writes dasm_loci_wg_annotated.parquet.
# Runs after W05_merge_wg.

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate allcools

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
NB_IN=$PROJDIR/notebooks/wg/W07_classify_annotate_wg.ipynb
NB_OUT=$PROJDIR/notebooks/wg/executed/W07_classify_annotate_wg.ipynb

mkdir -p $PROJDIR/notebooks/wg/executed

echo "$(date): W07 classify + annotate"
papermill \
    --kernel allcools \
    --log-output \
    $NB_IN \
    $NB_OUT

echo "$(date): W07 done"
