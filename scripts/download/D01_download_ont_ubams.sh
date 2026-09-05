#!/bin/bash
#$ -N D01_download_ont_ubams
#$ -cwd
#$ -l h_data=1G,h_rt=8:00:00
#$ -pe shared 1
#$ -t 1-18:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/D01_download_ont_ubams.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/D01_download_ont_ubams.$JOB_ID.$TASK_ID

# Download HPRC ONT uBAMs for the 18 Dorado (R10.4.1 5mCG_5hmCG) cohort members.
# Each task handles one sample; downloads all BAM files listed for that sample
# in the cohort manifest, skipping files that already exist.
#
# uBAMs are unaligned but carry MM/ML methylation tags (5mCG).
# After download, run W00_align_wg.sh to merge + align to hg38.
#
# Requires: aws CLI in igvf_utils conda env (public HPRC bucket, no credentials needed)
# Output:   bam/ubam/{SAMPLE}/{filename}.bam

set -euo pipefail

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate igvf_utils

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
MANIFEST=$PROJDIR/csv/meta/cohort_downloads.tsv
OUTBASE=$PROJDIR/bam/ubam

mkdir -p $PROJDIR/logs $OUTBASE

# ID=10  # hardcoded for testing; production uses SGE_TASK_ID
ID=$SGE_TASK_ID
# Map task ID → sample (18 Dorado samples only; see txt/samples/final_samples.txt)
SAMPLES=($(cat $PROJDIR/txt/samples/final_samples.txt))
SAMPLE=${SAMPLES[$((ID - 1))]}

echo "$(date): Task $ID — sample $SAMPLE"

# Get S3 directory for this sample from manifest
S3DIR=$(awk -v s=$SAMPLE '$1==s {print $7}' $MANIFEST)
if [ -z "$S3DIR" ]; then
    echo "ERROR: sample $SAMPLE not found in manifest"
    exit 1
fi

OUTDIR=$OUTBASE/$SAMPLE
mkdir -p $OUTDIR

echo "$(date): listing $S3DIR"
# List all .bam files in the S3 directory
BAM_KEYS=$(aws s3 ls --no-sign-request "$S3DIR" | awk '{print $4}' | grep '\.bam$' | grep -v '_fail\.bam$' || true)

if [ -z "$BAM_KEYS" ]; then
    echo "ERROR: no BAM files found at $S3DIR"
    exit 1
fi

N_FILES=$(echo "$BAM_KEYS" | wc -l)
echo "$(date): found $N_FILES BAM file(s)"

DOWNLOADED=0
SKIPPED=0
for BNAME in $BAM_KEYS; do
    DEST=$OUTDIR/$BNAME
    S3PATH="${S3DIR}${BNAME}"
    if [ -f "$DEST" ]; then
        echo "$(date): skipping (exists): $BNAME"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    echo "$(date): downloading $BNAME"
    time aws s3 cp --no-sign-request "$S3PATH" "$DEST"
    DOWNLOADED=$((DOWNLOADED + 1))
done

echo "$(date): $SAMPLE done — downloaded=$DOWNLOADED skipped=$SKIPPED"
echo "$(date): files in $OUTDIR:"
ls -lh $OUTDIR/*.bam 2>/dev/null | awk '{print $5, $9}'
