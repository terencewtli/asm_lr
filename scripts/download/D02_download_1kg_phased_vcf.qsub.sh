#!/bin/bash
#$ -cwd
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/D02_download_1kg_phased_vcf.$JOB_ID
#$ -j y
#$ -N D02_download_1kg_vcf
#$ -l h_data=4G,h_rt=24:00:00
#$ -pe shared 2

echo "Job $JOB_ID started on: " `hostname -s`
echo "Job $JOB_ID started on: " `date`
echo " "

time bash /u/project/cluo/terencew/claude/project_ideas/asm_lr/scripts/download/D02_download_1kg_phased_vcf.sh

echo "Job $JOB_ID done on: " `hostname -s`
echo "Job $JOB_ID done on: " `date`
