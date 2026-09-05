#!/bin/bash
#$ -N W01_modkit_extract_wg
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 10
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W01_modkit_extract_wg.$JOB_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W01_modkit_extract_wg.$JOB_ID
#$ -cwd

# Extract per-read 5mC calls (whole genome, no region filter)
# Output: single TSV covering all chromosomes — downstream scripts filter by chrom
# Runtime est: ~45 min for NA21093 full genome

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate whatshap-env

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
SAMPLE=NA21093
BAM=$PROJDIR/bam/${SAMPLE}_wg.bam
REF=/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa
OUTTSV=$PROJDIR/modkit/${SAMPLE}_wg_modkit_extract.tsv.gz
SCRATCH=/u/project/cluo_scratch/terencew/claude/asm_lr/modkit_extract
LOGFILE=$PROJDIR/logs/modkit_extract_${SAMPLE}_wg.$JOB_ID.log

mkdir -p $PROJDIR/modkit $SCRATCH

TMPOUT=$SCRATCH/${SAMPLE}_wg_modkit_extract_$JOB_ID.tsv

echo "$(date): modkit extract — $SAMPLE whole genome"
modkit extract \
    --ref $REF \
    --threads 10 \
    --log-filepath $LOGFILE \
    $BAM \
    $TMPOUT

echo "$(date): gzipping extract..."
gzip -c $TMPOUT > $OUTTSV
rm -f $TMPOUT

echo "$(date): done — $(zcat $OUTTSV | wc -l) rows"
