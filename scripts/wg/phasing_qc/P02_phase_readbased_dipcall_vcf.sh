#!/bin/bash
#$ -N P02_phase_readbased_dipcall_vcf
#$ -l h_data=6G,h_rt=8:00:00
#$ -pe shared 4
#$ -t 1-132:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/P02/P02_phase_readbased_dipcall_vcf.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/P02/P02_phase_readbased_dipcall_vcf.$JOB_ID.$TASK_ID
#$ -cwd

# whatshap phase, GENOME-WIDE (22 autosomes, not just chr1 -- widened
# 2026-09-07 once the project owner decided this should become the
# production phasing source, not just a pilot QC check), for the 6-donor
# phasing-QC pilot -- using the P01-generated assembly-backed dipcall VCF
# (vcf/read_based/{SAMPLE}_dip.vcf.gz, genome-wide by construction: dipcall
# aligns the whole hap1/hap2 assembly, not a single chromosome) as the
# het-SNP input instead of the 1000G population-panel VCF that
# W01b_phase_readbased_per_chrom.sh uses.
#
# One task per (sample, chrom), txt/samples/phasing_qc_pilot_chrom_tasks.txt
# (6 x 22 = 132 rows, same convention as G04's six_donor_chrom_tasks.txt).
#
# Rationale (see scripts/wg/phasing_qc/README.md and
# md/20260907_phasing_qc_review.md): the dipcall VCF carries every het site
# this donor's own assembly shows, common/rare/private, with no population-
# frequency dependence -- strictly more informative than the 1000G panel,
# which is why this is meant to REPLACE (not just QC-compare against) the
# panel-based haplotagging for donors it's run on.
#
# Output: vcf/read_based_phased/{SAMPLE}_{CHROM}_whatshap_phased_dipcall.vcf.gz (+.tbi)
# Next: re-run haplotag (W01-style) against THIS phased VCF instead of the
# panel VCF, then re-run W02/W03/beta-binomial calling downstream for
# whichever donors this replaces panel-phasing for.

set -euo pipefail

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate whatshap-env

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
REF=/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa
OUTDIR=$PROJDIR/vcf/read_based_phased

TASKS=$PROJDIR/txt/samples/phasing_qc_pilot_chrom_tasks.txt
LINE=$(sed -n "${SGE_TASK_ID}p" "$TASKS")
SAMPLE=$(echo "$LINE" | cut -f1)
CHROM=$(echo "$LINE" | cut -f2)

BAM=$PROJDIR/bam/merged_bam/${SAMPLE}_wg.bam
DIP_VCF=$PROJDIR/vcf/read_based/${SAMPLE}_dip.vcf.gz
OUT=$OUTDIR/${SAMPLE}_${CHROM}_whatshap_phased_dipcall.vcf.gz

mkdir -p "$OUTDIR" "$PROJDIR/logs/P02"

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
