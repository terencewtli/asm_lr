#!/bin/bash
#$ -N W02_modkit_extract_per_chrom
#$ -l h_data=2G,h_rt=6:00:00
#$ -pe shared 10
#$ -t 1-396:1
#$ -tc 40
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W02/W02_modkit_extract_per_chrom.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W02/W02_modkit_extract_per_chrom.$JOB_ID.$TASK_ID
#$ -cwd

# Canonical per-chromosome modkit extract: one task per (sample, chrom) in
# txt/samples/sample_chrom_tasks.txt (396 = 18 samples x 22 autosomes).
# Replaces the original whole-genome-at-once W02 (moved to
# old/W02_modkit_extract_wg_original_oom.sh) which OOM'd/panicked on 14/18
# samples because it silently gzipped a truncated file on crash and
# requested nowhere near enough memory for a 100-300GB BAM in one shot.
#
# NOT on the ASM-calling critical path: the Fisher caller (W03 pileup + W04)
# reads modkit *pileup* output for CpG-level testing and does its own small
# --include-bed-restricted extract for read-level phase-loci validation.
# This whole-genome-per-read table is only consumed by the retired raw-delta
# caller (old/W04_dasm_loci_rawdelta.py) and QC02_read_length_statistics.ipynb
# (which could be recomputed straight from the BAM more cheaply). Skip this
# job entirely if you don't need either of those.
#
# h_data=2G x pe shared 8 = 16GB per task. A prior run of this stage (job
# 14660928, manually qdel'ed once we switched pipelines rather than left to
# finish) dispatched 43/308 tasks including chr1 for two different samples
# (the worst case by data volume) with no OOM/crash symptoms in any log
# before the qdel -- also, per-chromosome tasks turned out much faster than
# expected (10-25 min for most chroms, chr1 still running at ~23 min), which
# suggests --region genuinely fixes the whole-genome version's blowup rather
# than just delaying it. Not fully proven (chr1 never ran to completion),
# but no evidence so far that 16GB is too low.
#
# Runs uniformly across all 18 samples for consistency, including the 4
# that already have a complete, correct WHOLE-GENOME extract from the
# original run (NA19700, HG00126, NA19776, HG03784) -- those legacy
# modkit/{sample}_wg_modkit_extract.tsv.gz files still exist and are fine to
# use as-is if you don't need the per-chromosome split for them specifically;
# this job does not touch or delete them, it only adds the new
# modkit/{sample}_{chrom}_modkit_extract.tsv.gz files alongside.

ulimit -c 0
source /u/home/t/terencew/project-cluo/miniconda3/bin/activate whatshap-env

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
REF=/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa
SCRATCH=/u/project/cluo_scratch/terencew/claude/asm_lr/modkit_extract

TASKS=$PROJDIR/txt/samples/sample_chrom_tasks.txt
LINE=$(sed -n "${SGE_TASK_ID}p" "$TASKS")
SAMPLE=$(echo "$LINE" | cut -f1)
CHROM=$(echo "$LINE" | cut -f2)

BAM=$PROJDIR/bam/haplotagged/${SAMPLE}_wg_haplotagged.bam
OUTTSV=$PROJDIR/modkit/${SAMPLE}_${CHROM}_modkit_extract.tsv.gz
LOGFILE=$PROJDIR/logs/W02/modkit_extract_${SAMPLE}_${CHROM}.$JOB_ID.log
TMPOUT=$SCRATCH/${SAMPLE}_${CHROM}_modkit_extract_$JOB_ID.tsv

echo "$(date): W02 task $SGE_TASK_ID — $SAMPLE $CHROM"
echo "Hostname: $(hostname)"

if [ ! -f "$BAM" ]; then
    echo "ERROR: $BAM not found"
    exit 1
fi

if [ -f "$OUTTSV" ] && [ -s "$OUTTSV" ]; then
    echo "$(date): $OUTTSV exists, skipping"
    exit 0
fi

mkdir -p "$PROJDIR/modkit" "$SCRATCH"

echo "$(date): modkit extract — $SAMPLE $CHROM"
time modkit extract \
    --ref "$REF" \
    --region "$CHROM" \
    --threads 8 \
    --log-filepath "$LOGFILE" \
    "$BAM" \
    "$TMPOUT"
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
    echo "$(date): ERROR — modkit extract exited $EXIT_CODE for $SAMPLE $CHROM"
    rm -f "$TMPOUT"
    exit $EXIT_CODE
fi

echo "$(date): gzipping extract..."
time gzip -c "$TMPOUT" > "$OUTTSV"
rm -f "$TMPOUT"

echo "$(date): done — $SAMPLE $CHROM — $(zcat "$OUTTSV" | wc -l) rows"
