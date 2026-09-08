#!/bin/bash
#$ -N P01_run_dipcall_pilot
#$ -cwd
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 8
#$ -t 1-6:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/P01/P01_run_dipcall_pilot.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/P01/P01_run_dipcall_pilot.$JOB_ID.$TASK_ID

# Pilot-scoped copy of scripts/download/D04_run_dipcall.sh: same dipcall
# pipeline (assembly vs. reference -> phased truth VCF), but restricted to
# the 6-donor phasing-QC pilot (txt/samples/phasing_qc_pilot_donors.txt:
# the 4 donors with anomalous beta-binomial ASM rates -- NA18508, HG01167,
# NA21093, NA18620 -- plus 2 confirmed-sane controls -- NA19700, HG00146)
# instead of qsub'ing D04 against the full 18-donor cohort. Assemblies are
# already downloaded for all 6 (D03 already ran for the full 18-donor set).
#
# Output: vcf/read_based/{SAMPLE}_dip.vcf.gz (+.tbi) -- same output path/
# format as D04, so this can be reused directly by D04-downstream steps or
# by P02_phase_readbased_dipcall_vcf.sh in this directory. If the full-cohort
# D04 run happens later, it will just see these 6 files already present and
# skip them (same existence check).
#
# See scripts/wg/phasing_qc/README.md and md/20260907_phasing_qc_review.md
# for why this pilot exists: neither coverage nor the beta-binomial caller
# switch explains 4 donors' anomalous ASM rates, and this cohort has never
# had a non-simulated (assembly-truth-based) switch-error measurement.

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
ASM_DIR=$PROJDIR/assemblies
OUTDIR=$PROJDIR/vcf/read_based
SCRATCH=/u/project/cluo_scratch/terencew/claude/asm_lr/dipcall_pilot
USERBIN=/u/home/t/terencew/bin
THREADS=8

MINIMAP2=$USERBIN/minimap2
BCFTOOLS=/u/local/apps/bcftools/1.11/gcc-4.8.5/bin/bcftools
TABIX=/u/local/apps/htslib/1.12/gcc-4.8.5/bin/tabix
K8=$USERBIN/k8
PAFTOOLS=$USERBIN/paftools.js
DIPCALL_AUX=$USERBIN/dipcall-aux.js
REF=/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa

mkdir -p "$OUTDIR" "$SCRATCH" "$PROJDIR/logs/P01"

if [ ! -f "$K8" ] || [ ! -f "$PAFTOOLS" ] || [ ! -f "$DIPCALL_AUX" ]; then
    echo "ERROR: k8/paftools.js/dipcall-aux.js not found in $USERBIN -- run"
    echo "scripts/download/D04_run_dipcall.sh's one-time setup first (or this"
    echo "script standalone; it downloads them the same way D04 does)."
    for tool_url in \
        "$K8 https://github.com/attractivechaos/k8/releases/download/v1.2/k8-x86_64-Linux" \
        "$PAFTOOLS https://raw.githubusercontent.com/lh3/minimap2/master/misc/paftools.js" \
        "$DIPCALL_AUX https://raw.githubusercontent.com/lh3/dipcall/master/dipcall-aux.js"; do
        f=$(echo "$tool_url" | cut -d' ' -f1)
        u=$(echo "$tool_url" | cut -d' ' -f2)
        if [ ! -f "$f" ]; then
            echo "$(date): downloading $(basename $f)..."
            wget -q -O "$f" "$u"
            chmod +x "$f"
        fi
    done
fi

SAMPLES=($(cat "$PROJDIR/txt/samples/phasing_qc_pilot_donors.txt"))
SAMPLE=${SAMPLES[$((SGE_TASK_ID - 1))]}

OUT_VCF=$OUTDIR/${SAMPLE}_dip.vcf.gz
if [ -f "$OUT_VCF" ] && [ -f "${OUT_VCF}.tbi" ]; then
    echo "$(date): $SAMPLE — dip.vcf.gz exists, skipping"
    exit 0
fi

HAP1=$ASM_DIR/${SAMPLE}_hap1.fa.gz
HAP2=$ASM_DIR/${SAMPLE}_hap2.fa.gz
if [ ! -f "$HAP1" ] || [ ! -f "$HAP2" ]; then
    echo "ERROR: assemblies missing for $SAMPLE — expected from D03, check $ASM_DIR"
    exit 1
fi

echo "$(date): P01 task $SGE_TASK_ID — $SAMPLE"
echo "Hostname: $(hostname)"
echo "  hap1: $(du -sh $HAP1 | cut -f1)"
echo "  hap2: $(du -sh $HAP2 | cut -f1)"

WORKDIR=$SCRATCH/$SAMPLE
mkdir -p "$WORKDIR"

for HAP in 1 2; do
    if [ "$HAP" -eq 1 ]; then ASM="$HAP1"; else ASM="$HAP2"; fi
    PAF=$WORKDIR/${SAMPLE}_hap${HAP}.paf.gz
    if [ ! -f "$PAF" ]; then
        echo "$(date): aligning hap${HAP}..."
        time $MINIMAP2 -cx asm5 --cs -r2k -t $THREADS "$REF" "$ASM" \
          | sort -k6,6 -k8,8n \
          | bgzip > "$PAF"
    fi
done

for HAP in 1 2; do
    PAF=$WORKDIR/${SAMPLE}_hap${HAP}.paf.gz
    VCF=$WORKDIR/${SAMPLE}_hap${HAP}.vcf
    if [ ! -f "$VCF" ]; then
        echo "$(date): calling variants hap${HAP}..."
        zcat "$PAF" | $K8 "$PAFTOOLS" call -l 5000 -q 5 - > "$VCF"
        N=$(grep -cv "^#" "$VCF" || true)
        echo "$(date): hap${HAP} — $N variant calls"
    fi
done

echo "$(date): merging haplotypes..."
VCF1=$WORKDIR/${SAMPLE}_hap1.vcf
VCF2=$WORKDIR/${SAMPLE}_hap2.vcf
TMP_VCF=$WORKDIR/${SAMPLE}_dip_tmp.vcf.gz

$K8 "$DIPCALL_AUX" merge "$VCF1" "$VCF2" \
  | $BCFTOOLS sort \
  | bgzip > "$TMP_VCF"
$TABIX "$TMP_VCF"

N_HET=$(zcat "$TMP_VCF" | grep -v "^#" | awk '$10 ~ /0[|\/]1|1[|\/]0/' | wc -l || true)
echo "$(date): $N_HET het variants in merged VCF"

mv "$TMP_VCF"       "$OUT_VCF"
mv "${TMP_VCF}.tbi" "${OUT_VCF}.tbi"

rm -rf "$WORKDIR"

echo "$(date): P01 done — $SAMPLE"
echo "  $OUT_VCF  ($(du -sh $OUT_VCF | cut -f1))"
