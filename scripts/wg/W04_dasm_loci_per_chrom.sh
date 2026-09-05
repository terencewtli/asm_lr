#!/bin/bash
#$ -N W04_dasm_loci_per_chrom
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 8
#$ -t 1-396:1
#$ -hold_jid_ad W03_modkit_pileup_per_chrom
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W04/W04_dasm_loci_per_chrom.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W04/W04_dasm_loci_per_chrom.$JOB_ID.$TASK_ID
#$ -cwd

# Fisher-based dASM calling, split out from the pileup step (W03) into its
# own array so it can be re-tuned/rerun cheaply without re-paying the
# multi-hour modkit pileup cost each time. Same 396-row manifest as W03
# (txt/samples/sample_chrom_tasks.txt), so -hold_jid_ad chains task i of this
# array to task i of W03 specifically -- each (sample, chrom) ASM call
# starts as soon as its own pileup is done, not after the whole W03 array.
#
# Per task: call-loci (Fisher + BH-FDR + cluster) -> targeted modkit extract
# (--include-bed, loci only) -> phase-loci (Mann-Whitney U + AUC).
#
# NOT YET TIMED — the NA19700/chr21 pilot never got past the pileup step
# (killed by h_rt before reaching call-loci). h_data=4G/pe shared 8/h_rt=4h
# here are placeholders scaled down from the pileup step on the assumption
# that CpG-level Fisher testing (benchmarked ~750us/call, parallel over
# pe shared workers) and a loci-restricted extract are both much cheaper
# than a full chromosome pileup -- validate against the pilot (once rerun
# with its bgzip fix) before trusting this at full 396-task scale.

ulimit -c 0
source /u/home/t/terencew/project-cluo/miniconda3/bin/activate whatshap-env

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
REF=/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa
SCRIPT=$PROJDIR/scripts/wg/W04_dasm_loci.py
SCRATCH=/u/project/cluo_scratch/terencew/claude/asm_lr/dasm_fisher

TASKS=$PROJDIR/txt/samples/sample_chrom_tasks.txt
LINE=$(sed -n "${SGE_TASK_ID}p" "$TASKS")
SAMPLE=$(echo "$LINE" | cut -f1)
CHROM=$(echo "$LINE" | cut -f2)

BAM=$PROJDIR/bam/haplotagged/${SAMPLE}_wg_haplotagged.bam
PILEUP_DIR=$PROJDIR/pileup
# W04_dasm_loci_fisher.py names its outputs by --chrom only (dasm_loci_{chrom}_fisher.tsv,
# phasing_results_{chrom}_fisher.parquet, etc). With 18 samples sharing chrom names,
# two tasks for the same chromosome (different samples) running concurrently against
# a shared outdir would clobber each other's intermediate files. Giving each task its
# own per-sample outdir avoids that -- no rename-at-the-end needed.
OUTDIR=$PROJDIR/tables/dasm_fisher/$SAMPLE

mkdir -p "$OUTDIR" "$SCRATCH"

echo "$(date): W04 task $SGE_TASK_ID — $SAMPLE $CHROM"
echo "Hostname: $(hostname)"

if [ ! -f "${PILEUP_DIR}/${SAMPLE}_${CHROM}_hp_pileup_1.bed.gz" ]; then
    echo "ERROR: ${PILEUP_DIR}/${SAMPLE}_${CHROM}_hp_pileup_1.bed.gz not found — W03 task for this (sample, chrom) hasn't produced it"
    exit 1
fi

if [ -f "$OUTDIR/phasing_results_${CHROM}_fisher.parquet" ]; then
    echo "$(date): $OUTDIR/phasing_results_${CHROM}_fisher.parquet exists, skipping"
    exit 0
fi

echo "$(date): [1/3] call-loci — $SAMPLE $CHROM (Fisher + BH-FDR + clustering)"
time python3 "$SCRIPT" call-loci \
    --chrom "$CHROM" --sample "$SAMPLE" \
    --pileup-dir "$PILEUP_DIR" --outdir "$OUTDIR" --threads 8
echo "$(date): [1/3] done"

LOCI_BED=$OUTDIR/dasm_loci_${CHROM}_fisher.bed
N_LOCI=$(wc -l < "$LOCI_BED")
echo "$(date): [2/3] modkit extract — $SAMPLE $CHROM restricted to $N_LOCI candidate loci"
if [ "$N_LOCI" -gt 0 ]; then
    TMPOUT=$SCRATCH/${SAMPLE}_${CHROM}_loci_extract_$JOB_ID.tsv
    time modkit extract \
        --ref "$REF" \
        --region "$CHROM" \
        --include-bed "$LOCI_BED" \
        --threads 8 \
        --log-filepath "$PROJDIR/logs/W04/modkit_extract_dasm_${SAMPLE}_${CHROM}.$JOB_ID.log" \
        "$BAM" \
        "$TMPOUT"
    EXIT_CODE=$?
    if [ $EXIT_CODE -ne 0 ]; then
        echo "$(date): ERROR — modkit extract exited $EXIT_CODE for $SAMPLE $CHROM"
        rm -f "$TMPOUT"
        exit $EXIT_CODE
    fi
    time gzip -f "$TMPOUT"
    EXTRACT_PATH="${TMPOUT}.gz"
else
    echo "$(date): no candidate loci — skipping extract"
    EXTRACT_PATH=""
fi
echo "$(date): [2/3] done"

echo "$(date): [3/3] phase-loci — $SAMPLE $CHROM (Mann-Whitney U + AUC)"
if [ -n "$EXTRACT_PATH" ]; then
    time python3 "$SCRIPT" phase-loci \
        --chrom "$CHROM" --sample "$SAMPLE" \
        --extract-path "$EXTRACT_PATH" --bam-path "$BAM" --outdir "$OUTDIR"
    rm -f "$EXTRACT_PATH"
fi
echo "$(date): [3/3] done — $SAMPLE $CHROM complete"
