#!/bin/bash
#$ -N P02_phase_readbased_dipcall_vcf
#$ -l h_data=6G,h_rt=8:00:00
#$ -pe shared 4
#$ -t 1-6:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/P02/P02_phase_readbased_dipcall_vcf.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/P02/P02_phase_readbased_dipcall_vcf.$JOB_ID.$TASK_ID
#$ -cwd

# whatshap phase, chr1 only, for the 6-donor phasing-QC pilot -- but using
# the P01-generated assembly-backed dipcall VCF (vcf/read_based/{SAMPLE}_dip.vcf.gz)
# as the het-SNP input instead of the 1000G population-panel VCF that
# W01b_phase_readbased_per_chrom.sh uses. Rationale (see
# scripts/wg/phasing_qc/README.md and md/20260907_phasing_qc_review.md):
#
#   - the panel VCF only carries variants in the 1000G panel, which
#     underrepresents an individual's rare/private heterozygous sites --
#     refining ITS phase with whatshap phase fixes switch errors among
#     already-known common variants, but doesn't add phasing-informative
#     density anywhere the panel itself is sparse.
#   - the dipcall VCF (from P01) is built directly from this donor's own
#     assembled diploid genome, so it carries every het site the donor
#     actually has -- common, rare, and private -- with no population-
#     frequency dependence at all.
#
# Running whatshap phase with the dipcall VCF + the ONT BAM combines that
# higher-density, assembly-backed base phase with the reads' own physical
# linkage -- a strictly higher information-content input than either the
# panel-only or dipcall-only phase alone, and the natural next step to
# compare against panel-only haplotagging rate.
#
# Output: vcf/read_based_phased/{SAMPLE}_chr1_whatshap_phased_dipcall.vcf.gz (+.tbi)
# (Filename kept distinct from W01b's *_whatshap_phased.vcf.gz output, which
# is panel-VCF-based, so both can be compared side by side.)

set -euo pipefail

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate whatshap-env

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
REF=/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa
OUTDIR=$PROJDIR/vcf/read_based_phased
CHROM=chr1

SAMPLES=($(cat "$PROJDIR/txt/samples/phasing_qc_pilot_donors.txt"))
SAMPLE=${SAMPLES[$((SGE_TASK_ID - 1))]}

BAM=$PROJDIR/bam/merged_bam/${SAMPLE}_wg.bam
DIP_VCF=$PROJDIR/vcf/read_based/${SAMPLE}_dip.vcf.gz
OUT=$OUTDIR/${SAMPLE}_${CHROM}_whatshap_phased_dipcall.vcf.gz

mkdir -p "$OUTDIR"

echo "$(date): P02 task $SGE_TASK_ID — $SAMPLE $CHROM"
echo "Hostname: $(hostname)"

if [ -f "$OUT" ] && [ -f "${OUT}.tbi" ]; then
    echo "$(date): $OUT exists, skipping"
    exit 0
fi
if [ ! -f "$BAM" ]; then
    echo "ERROR: $BAM not found — run W00b_merge_wg.sh first"
    exit 1
fi
if [ ! -f "$DIP_VCF" ]; then
    echo "ERROR: $DIP_VCF not found — run P01_run_dipcall_pilot.sh first"
    exit 1
fi

# Restrict the dipcall VCF (genome-wide) to this chromosome's het sites,
# matching W01b's per-chromosome scope.
CHROM_HET_VCF=$OUTDIR/${SAMPLE}_${CHROM}_dip_het.vcf.gz
bcftools view -r "$CHROM" -i 'GT="het"' "$DIP_VCF" -Oz -o "$CHROM_HET_VCF"
tabix -p vcf "$CHROM_HET_VCF"

TMP_OUT=$OUTDIR/${SAMPLE}_${CHROM}_whatshap_phased_dipcall.tmp.vcf

time whatshap phase \
    --reference "$REF" \
    --chromosome "$CHROM" \
    --ignore-read-groups \
    --output "$TMP_OUT" \
    "$CHROM_HET_VCF" \
    "$BAM"

bgzip -f "$TMP_OUT"
mv "${TMP_OUT}.gz" "$OUT"
tabix -p vcf "$OUT"

echo "$(date): P02 done — $SAMPLE $CHROM"
echo "  $OUT"
