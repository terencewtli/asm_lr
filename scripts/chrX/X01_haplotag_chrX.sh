#!/bin/bash
#$ -N X01_haplotag_chrX
#$ -l h_data=2G,h_rt=8:00:00
#$ -pe shared 4
#$ -t 1-8:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/X01/X01_haplotag_chrX.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/X01/X01_haplotag_chrX.$JOB_ID.$TASK_ID
#$ -cwd

# whatshap haplotag on chrX, FEMALE DONORS ONLY (array 1-8, matching
# txt/samples/female_donors.txt). Non-PAR chrX is hemizygous in males (one
# copy, maternal-origin only) -- there's no haplotype pair to compare, so
# male chrX haplotagging isn't meaningful here and isn't run by this script.
# (Male chrX bulk coverage is still useful as a built-in negative control --
# see scripts/chrX/README.md -- but that just needs the X00b merged BAM,
# not haplotagging.)
#
# *** REQUIRES A HET-SNP CHRX VCF SOURCE THIS SCRIPT DOES NOT YET HAVE. ***
# The autosome pipeline uses the 1000G 3202-sample statistically-phased
# panel (per-chromosome het VCFs already staged under vcf/read_based/).
# For chrX, the project owner has a separately TOPMed-imputed, statistically
# phased chrX VCF already generated -- NOT yet located/wired in here. Set
# CHRX_VCF_SOURCE below (or change the per-sample lookup) before running.
# The whatshap haplotag INVOCATION itself needs no chrX-specific change from
# the autosome version (chrX is genuinely diploid in females outside PAR, so
# no --ploidy flag or other special-casing is needed) -- only the VCF INPUT
# differs.

set -euo pipefail

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate whatshap-env

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
REF=$PROJDIR/reference/chrX_ref/GRCh38.chrX.fa
BCFTOOLS=/u/local/apps/bcftools/1.11/gcc-4.8.5/bin/bcftools
SAMTOOLS=/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools

CHRX_VCF_SOURCE=""   # <-- FILL IN: path (or path template) to the TOPMed-imputed,
                      #     phased chrX het-SNP VCF, per sample or a single
                      #     cohort-wide file. Not yet wired in.

SAMPLES=($(cat "$PROJDIR/txt/samples/female_donors.txt"))
SAMPLE=${SAMPLES[$((SGE_TASK_ID - 1))]}

BAM=$PROJDIR/bam/chrX/merged_bam/${SAMPLE}_chrX.bam
OUTDIR=$PROJDIR/bam/chrX/haplotagged
OUT_TAGGED=$OUTDIR/${SAMPLE}_chrX_haplotagged.bam
mkdir -p "$OUTDIR"

echo "$(date): X01 task $SGE_TASK_ID — $SAMPLE chrX"
echo "Hostname: $(hostname)"

if [ -z "$CHRX_VCF_SOURCE" ]; then
    echo "ERROR: CHRX_VCF_SOURCE is not set -- point this at the TOPMed chrX"
    echo "  phased VCF before running. See this script's header comment."
    exit 1
fi
if [ -f "$OUT_TAGGED" ] && [ -s "$OUT_TAGGED" ]; then
    echo "$(date): $OUT_TAGGED exists — skipping"
    exit 0
fi
if [ ! -f "$BAM" ]; then
    echo "ERROR: $BAM not found — run X00b_merge_chrX.sh first"
    exit 1
fi

HET_VCF=$PROJDIR/vcf/read_based/${SAMPLE}_chrX_het_snps.vcf.gz
if [ ! -f "$HET_VCF" ]; then
    echo "$(date): extracting chrX het SNPs for $SAMPLE from $CHRX_VCF_SOURCE..."
    $BCFTOOLS view -s "$SAMPLE" -r chrX -i 'GT="het"' "$CHRX_VCF_SOURCE" \
        -Oz -o "$HET_VCF"
    tabix -p vcf "$HET_VCF"
fi

# Same haplotag invocation as the autosome pipeline (W01_phase_per_chrom.sh)
# -- no ploidy flag needed, chrX is diploid in females outside PAR.
time whatshap haplotag \
    --reference "$REF" \
    --ignore-read-groups \
    --output "$OUT_TAGGED" \
    "$HET_VCF" \
    "$BAM"
$SAMTOOLS index "$OUT_TAGGED"

echo "$(date): X01 done — $OUT_TAGGED"
