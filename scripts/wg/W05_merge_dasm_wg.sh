#!/bin/bash
#$ -N W05_merge_dasm_wg
#$ -l h_data=4G,h_rt=1:00:00
#$ -pe shared 1
#$ -hold_jid W04_dasm_loci_per_chrom
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W05/W05_merge_dasm_wg.$JOB_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W05/W05_merge_dasm_wg.$JOB_ID
#$ -cwd

# The one merge step in the whole pipeline: per-chromosome Fisher ASM output
# (tables/dasm_fisher/{sample}/{dasm_loci,phasing_results}_{chrom}_fisher.*)
# -> one whole-genome table per sample -> one cohort-wide table. Everything
# upstream of this (W02 extract, W03 pileup, W04 ASM calling) stays
# per-chromosome all the way through; nothing else in the pipeline merges.
#
# Single task, not an array -- concatenating small per-chrom tables is cheap
# (low minutes for the whole cohort), no need for cluster parallelism here.
#
# Output:
#   tables/dasm_fisher/{sample}_wg_dasm_loci_fisher.tsv
#   tables/dasm_fisher/{sample}_wg_phasing_results_fisher.parquet
#   tables/dasm_fisher/cohort_wg_dasm_loci_fisher.parquet       (all samples, + sample column)
#   tables/dasm_fisher/cohort_wg_phasing_results_fisher.parquet (all samples, + sample column)

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate whatshap-env

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
time python3 "$PROJDIR/scripts/wg/W05_merge_dasm_wg.py"
