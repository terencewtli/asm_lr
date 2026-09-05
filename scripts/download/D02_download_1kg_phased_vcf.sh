#!/bin/bash
# Pre-download 1000G phased het SNPs for all 30 panel members.
# Strategy: stream each per-chrom 1000G VCF ONCE, extract all 30 samples,
# then split per-sample locally — avoids 660 independent S3 streams.
#
# Output: vcf/read_based/{SAMPLE}_{CHROM}_het_snps.vcf.gz (+.tbi) — 22 files per sample
# These are consumed directly by W02_phase_wg.sh (which skips generation if files exist).
#
# Run on login node: bash scripts/download/D02_download_1kg_phased_vcf.sh
# Resume-safe: skips chroms where all 30 per-sample files exist; skips individual
# per-sample files that already have a .tbi index.

set -euo pipefail

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate igvf_utils

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
MANIFEST=$PROJDIR/csv/meta/cohort_downloads.tsv
OUTDIR=$PROJDIR/vcf/read_based/by_chrom
TMPDIR=/u/project/cluo_scratch/terencew/claude/asm_lr/vcf_tmp
BCFTOOLS=/u/local/apps/bcftools/1.11/gcc-4.8.5/bin/bcftools
TABIX=/u/local/apps/htslib/1.12/gcc-4.8.5/bin/tabix

mkdir -p "$OUTDIR" "$TMPDIR" "$PROJDIR/logs"

KG_BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"

SAMPLES=($(awk 'NR>1 {print $1}' "$MANIFEST"))
SAMPLES_STR=$(IFS=,; echo "${SAMPLES[*]}")
N_SAMPLES=${#SAMPLES[@]}

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
        chr20 chr21 chr22)

echo "$(date): D02 start — $N_SAMPLES samples, ${#CHROMS[@]} chroms"
echo "Samples: ${SAMPLES[*]}"

for CHROM in "${CHROMS[@]}"; do
    # Check if all per-sample files already exist
    ALL_DONE=1
    for SAMPLE in "${SAMPLES[@]}"; do
        OUT=$OUTDIR/${SAMPLE}_${CHROM}_het_snps.vcf.gz
        if [ ! -f "$OUT" ] || [ ! -f "${OUT}.tbi" ]; then
            ALL_DONE=0
            break
        fi
    done
    if [ "$ALL_DONE" -eq 1 ]; then
        echo "$(date): $CHROM — all $N_SAMPLES samples done, skipping"
        continue
    fi

    TMP_VCF=$TMPDIR/${CHROM}_30samples.vcf.gz
    EBI_URL="${KG_BASE}/1kGP_high_coverage_Illumina.${CHROM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

    if [ ! -f "$TMP_VCF" ] || [ ! -f "${TMP_VCF}.tbi" ]; then
        echo "$(date): $CHROM — streaming 1000G VCF via EBI HTTPS (all 30 samples in one pass)..."
        time curl -s "$EBI_URL" \
          | $BCFTOOLS view -s "$SAMPLES_STR" -O z -o "$TMP_VCF" -
        $TABIX "$TMP_VCF"
        echo "$(date): $CHROM — batch downloaded ($(du -sh $TMP_VCF | cut -f1))"
    else
        echo "$(date): $CHROM — batch VCF cached, splitting only"
    fi

    # Split into per-sample het VCFs (local, fast)
    for SAMPLE in "${SAMPLES[@]}"; do
        OUT=$OUTDIR/${SAMPLE}_${CHROM}_het_snps.vcf.gz
        if [ -f "$OUT" ] && [ -f "${OUT}.tbi" ]; then
            continue
        fi
        $BCFTOOLS view -s "$SAMPLE" -g het -c 1 -O z -o "$OUT" "$TMP_VCF"
        $TABIX "$OUT"
        N=$(zcat "$OUT" | grep -cv "^#" || true)
        echo "  ${SAMPLE} ${CHROM}: $N het SNPs"
    done

    # Remove chrom-level tmp once all samples split
    rm -f "$TMP_VCF" "${TMP_VCF}.tbi"
    echo "$(date): $CHROM — done"
done

echo "$(date): D02 complete. Files in $OUTDIR:"
ls "$OUTDIR"/*.vcf.gz 2>/dev/null | wc -l
