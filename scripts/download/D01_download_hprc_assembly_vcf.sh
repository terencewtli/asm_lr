#!/bin/bash
#$ -N D01_download_hprc_assembly_vcf
#$ -cwd
#$ -l h_data=2G,h_rt=4:00:00
#$ -pe shared 2
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/D01_download_hprc_assembly_vcf.$JOB_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/D01_download_hprc_assembly_vcf.$JOB_ID

# Download HPRC Year 1 dipcall VCFs for the 8 Y1-assembly samples in the 30-person panel.
# These VCFs provide assembly-backed phasing with no switch errors within contigs (N50 ~40-60 Mb).
# The remaining 22 panel members are R2-only and use the 1000G phased VCF via W02.
#
# Y1 samples with dipcall VCFs: HG00621 HG01258 HG01361 HG01891 HG02572 HG03453 HG00438 HG01928
# Output: vcf/read_based/{SAMPLE}_dip.vcf.gz (+ .tbi)
#
# Usage: qsub D01_download_hprc_assembly_vcf.sh
#        bash D01_download_hprc_assembly_vcf.sh [SAMPLE]  # single sample

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
OUTDIR=$PROJDIR/vcf/read_based
mkdir -p "$OUTDIR" "$PROJDIR/logs"

BASE=https://human-pangenomics.s3.amazonaws.com/working/HPRC
DIPCALL_PATH=assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.1

# All 8 Y1 samples — HG01258/HG01361/HG00621 already downloaded but script is idempotent
Y1_SAMPLES=(HG00621 HG01258 HG01361 HG01891 HG02572 HG03453 HG00438 HG01928)

if [ $# -ge 1 ]; then
    SAMPLES=("$1")
else
    SAMPLES=("${Y1_SAMPLES[@]}")
fi

for SAMPLE in "${SAMPLES[@]}"; do
    VCF_URL="${BASE}/${SAMPLE}/${DIPCALL_PATH}/${SAMPLE}.f1_assembly_v2.dip.vcf.gz"
    OUT="${OUTDIR}/${SAMPLE}_dip.vcf.gz"

    if [ -f "$OUT" ] && [ -f "${OUT}.tbi" ]; then
        echo "$(date) — ${SAMPLE}: already present, skipping"
        continue
    fi

    echo "$(date) — ${SAMPLE}: downloading dip.vcf.gz..."
    time wget -q -c --tries=5 --waitretry=10 -O "$OUT" "$VCF_URL"

    echo "$(date) — ${SAMPLE}: indexing..."
    /u/local/apps/htslib/1.12/gcc-4.8.5/bin/tabix -p vcf "$OUT"

    N=$(zcat "$OUT" | grep -cv "^#" || true)
    echo "$(date) — ${SAMPLE}: ${N} variants  →  $OUT"
done

echo "$(date) — done. Files in $OUTDIR:"
ls -lh "$OUTDIR"/*.vcf.gz 2>/dev/null | awk '{print $5, $9}'
