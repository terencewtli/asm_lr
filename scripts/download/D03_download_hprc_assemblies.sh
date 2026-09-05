#!/bin/bash
#$ -N D03_download_hprc_assemblies
#$ -cwd
#$ -l h_data=1G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-18:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/D03_download_hprc_assemblies.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/D03_download_hprc_assemblies.$JOB_ID.$TASK_ID

# Download hap1/hap2 assembly FASTAs for the 18 Dorado panel members from HPRC S3.
# Y1 samples: year1_f1_assembly_v2_genbank (paternal/maternal .fa.gz)
# R2 samples: release2/{SAMPLE}_hap1/hap2_hprc_r2_v1.0.1.fa.gz
#
# Output: assemblies/{SAMPLE}_hap1.fa.gz and assemblies/{SAMPLE}_hap2.fa.gz
# These are inputs for D04_run_dipcall.sh (assembly-backed phasing for all 30 samples).
#
# Prerequisite: D01_download_ont_ubams.sh (igvf_utils env with aws CLI)

set -euo pipefail

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate igvf_utils

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
MANIFEST=$PROJDIR/csv/meta/cohort_downloads.tsv
OUTDIR=$PROJDIR/assemblies
mkdir -p "$OUTDIR" "$PROJDIR/logs"

SAMPLES=($(cat "$PROJDIR/txt/samples/final_samples.txt"))
SAMPLE=${SAMPLES[$((SGE_TASK_ID - 1))]}

HAP1_S3=$(awk -F'\t' -v s="$SAMPLE" '$1==s {print $11}' "$MANIFEST")
HAP2_S3=$(awk -F'\t' -v s="$SAMPLE" '$1==s {print $12}' "$MANIFEST")

OUT_HAP1=$OUTDIR/${SAMPLE}_hap1.fa.gz
OUT_HAP2=$OUTDIR/${SAMPLE}_hap2.fa.gz

echo "$(date): D03 task $SGE_TASK_ID — $SAMPLE"
echo "  hap1: $HAP1_S3"
echo "  hap2: $HAP2_S3"

if [ -z "$HAP1_S3" ] || [ -z "$HAP2_S3" ]; then
    echo "ERROR: assembly S3 paths not found for $SAMPLE"
    exit 1
fi

for HAP in 1 2; do
    if [ "$HAP" -eq 1 ]; then S3_PATH="$HAP1_S3"; OUT="$OUT_HAP1"
    else                       S3_PATH="$HAP2_S3"; OUT="$OUT_HAP2"
    fi

    if [ -f "$OUT" ]; then
        echo "$(date): hap${HAP} exists, skipping: $OUT"
        continue
    fi

    echo "$(date): downloading hap${HAP}..."
    time aws s3 cp --no-sign-request "$S3_PATH" "$OUT"
    SIZE=$(du -sh "$OUT" | cut -f1)
    echo "$(date): hap${HAP} done — $SIZE"
done

echo "$(date): $SAMPLE assemblies complete"
echo "  $(du -sh $OUT_HAP1 | cut -f1)  $OUT_HAP1"
echo "  $(du -sh $OUT_HAP2 | cut -f1)  $OUT_HAP2"
echo "Next: qsub D04_run_dipcall.sh (after all 30 tasks complete)"
