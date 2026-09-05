#!/bin/bash
#$ -cwd
#$ -o logs/A01a_bam_genozip.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A01a_bam_genozip
#$ -l h_data=4G,h_rt=24:00:00
#$ -pe shared 6
#$ -t 1-139:1

# ubams here run 2.4-95GB (median ~49GB) — much larger than the igvf multiome
# bams this template was copied from. h_rt=24:00:00 is a first-pass estimate,
# not measured on this data; check `grep -a real logs/A01a_bam_genozip.*` after
# the first few tasks finish and lower it once real per-GB throughput is known,
# so future reruns queue faster.

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
cd $PROJDIR

# ID=1
ID=${SGE_TASK_ID}
BAMS=$PROJDIR/txt/genozip/ubam_paths.txt
BAM=$(head -${ID} $BAMS | tail -1)

GENOZIP=/u/home/t/terencew/project-cluo/programs/genozip-linux-x86_64/genozip

time $GENOZIP $BAM -@ 8

echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "
