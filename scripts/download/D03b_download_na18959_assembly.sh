#!/bin/bash
#$ -N D03b_download_na18959_assembly
#$ -cwd
#$ -l h_data=2G,h_rt=2:00:00
#$ -pe shared 2
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/D03b_download_na18959_assembly.$JOB_ID
#$ -j y

# Download NA18959 hap1+hap2 assembly FASTAs from HPRC S3.
# NA18959 is the only Dorado sample in the 18-sample panel missing from asm_lr/assemblies/.
# R2 release: s3://human-pangenomics/working/HPRC/NA18959/assemblies/release2/

set -euo pipefail

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate igvf_utils

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
OUTDIR=$PROJDIR/assemblies
mkdir -p "$OUTDIR" "$PROJDIR/logs"

SAMPLE=NA18959
BASE=s3://human-pangenomics/working/HPP/${SAMPLE}/assemblies/release2

echo "Start: $(date)"
echo "Sample: $SAMPLE"

for HAP in hap1 hap2; do
    OUT=$OUTDIR/${SAMPLE}_${HAP}.fa.gz
    if [ -f "$OUT" ] && [ -s "$OUT" ]; then
        echo "$(date): $OUT exists — skipping"
        continue
    fi
    SRC="${BASE}/${SAMPLE}_${HAP}_hprc_r2_v1.0.1.fa.gz"
    echo "$(date): downloading $SRC"
    time aws s3 cp --no-sign-request "$SRC" "$OUT"
    samtools faidx "$OUT" 2>/dev/null || true
    echo "$(date): done — $(du -sh $OUT | cut -f1)"
done

ls -lh $OUTDIR/${SAMPLE}_hap1.fa.gz $OUTDIR/${SAMPLE}_hap2.fa.gz
echo "End: $(date)"
