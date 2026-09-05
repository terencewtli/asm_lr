#!/bin/bash
#$ -N W12_epiallele_entropy
#$ -cwd
#$ -l h_data=16G,h_rt=2:00:00
#$ -pe shared 4
#$ -t 1-22:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W12_epiallele_entropy.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W12_epiallele_entropy.$JOB_ID.$TASK_ID

# Per-chromosome epiallele entropy scoring.
# After all 22 tasks finish, run the merge step:
#   source activate allcools && python scripts/wg/W12_epiallele_entropy.py

. /u/local/Modules/default/init/modules.sh
source /u/home/t/terencew/project-cluo/miniconda3/bin/activate allcools

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
cd "$PROJDIR"
mkdir -p logs

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
        chr20 chr21 chr22)
CHROM=${CHROMS[$((SGE_TASK_ID - 1))]}

echo "Start: $(date)  Chrom: $CHROM"
python -u scripts/wg/W12_epiallele_entropy.py --chrom "$CHROM"
echo "End: $(date)"
