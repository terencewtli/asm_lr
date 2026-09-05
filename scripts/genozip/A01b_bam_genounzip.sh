#!/bin/bash
#$ -cwd
#$ -o logs/A01b_genounzip_bam.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A01b_genounzip_bam
#$ -l h_data=1G,h_rt=2:00:00
#$ -pe shared 1
#$ -t 1-20:1

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

PROJDIR=/u/home/t/terencew/project-cluo/igvf/2023_YR2/multiome
cd $PROJDIR

# ID=1
ID=${SGE_TASK_ID}
BAMS=$PROJDIR/txt/genozip/crarc_bam_paths.txt
BAM=$(head -${ID} $BAMS | tail -1)

GENOZIP=/u/home/t/terencew/project-cluo/programs/genozip-linux-x86_64/genounzip
GENO=${BAM}.genozip

time $GENOZIP $GENO -@ 8 --output $BAM

echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "
