#!/bin/bash
#$ -N W03_modkit_pileup_per_chrom
#$ -l h_data=2G,h_rt=8:00:00
#$ -pe shared 10
#$ -t 1-396:1
#$ -tc 40
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W03/W03_modkit_pileup_per_chrom.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W03/W03_modkit_pileup_per_chrom.$JOB_ID.$TASK_ID
#$ -cwd

# HP-stratified modkit pileup, generalized from the HG01258-only pilot to the
# full 18-sample x 22-autosome cohort (396 tasks = txt/samples/sample_chrom_tasks.txt).
#
# Runs directly against the merged whole-genome haplotagged BAM
# (bam/haplotagged/{sample}_wg_haplotagged.bam) with --region $CHROM --
# modkit pileup does index-assisted region access, so no per-chromosome BAM
# split is needed (the old per-sample-hardcoded version assumed one existed;
# it doesn't for this cohort).
#
# Timing (from the NA19700/chr21 pilot, job 14650896): ~4.9 kb/s of alignable
# sequence per task at 8-10 threads -> chr21 (~41.6Mb alignable) took ~2.3-2.5h.
# Largest chrom (chr1, ~230Mb alignable) x largest sample (NA21093, 1.55x
# NA19700's data density) could run ~18-20h -- h_rt=24h covers that with margin.
# h_data=6G x 8 slots = 48GB is generous relative to the ~8.6GB maxvmem modkit
# pileup actually used on the pilot (pileup streams in bounded chunks, unlike
# modkit extract's per-read tables -- nowhere near the W02 OOM situation).
#
# -tc 40 caps concurrent running tasks so we don't hit 396 processes doing
# indexed random access into multi-hundred-GB BAMs on the shared filesystem
# at once. Adjust down if this contends badly with other jobs, up if the
# queue has headroom.

ulimit -c 0
source /u/home/t/terencew/project-cluo/miniconda3/bin/activate whatshap-env

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
REF=/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa
OUTDIR=$PROJDIR/pileup

TASKS=$PROJDIR/txt/samples/sample_chrom_tasks.txt
LINE=$(sed -n "${SGE_TASK_ID}p" "$TASKS")
SAMPLE=$(echo "$LINE" | cut -f1)
CHROM=$(echo "$LINE" | cut -f2)

BAM=$PROJDIR/bam/haplotagged/${SAMPLE}_wg_haplotagged.bam
PREFIX=${OUTDIR}/${SAMPLE}_${CHROM}_hp_pileup
LOGFILE=$PROJDIR/logs/W03/modkit_pileup_${SAMPLE}_${CHROM}.$JOB_ID.log

echo "$(date): W03 task $SGE_TASK_ID — $SAMPLE $CHROM"
echo "Hostname: $(hostname)"

if [ ! -f "$BAM" ]; then
    echo "ERROR: $BAM not found — run W01_phase_per_chrom.sh + W01g_merge_haplotagged_wg.sh first"
    exit 1
fi

if [ -f "${PREFIX}_1.bed.gz" ] && [ -f "${PREFIX}_2.bed.gz" ]; then
    echo "$(date): ${PREFIX}_{1,2}.bed.gz exist, skipping"
    exit 0
fi

mkdir -p "$OUTDIR"

echo "$(date): modkit pileup — $SAMPLE $CHROM (HP-partitioned)"
time modkit pileup \
    --ref "$REF" \
    --region "$CHROM" \
    --partition-tag HP \
    --prefix "$PREFIX" \
    --threads 8 \
    --log-filepath "$LOGFILE" \
    "$BAM" \
    "$OUTDIR"
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
    echo "$(date): ERROR — modkit pileup exited $EXIT_CODE for $SAMPLE $CHROM"
    exit $EXIT_CODE
fi

for f in ${PREFIX}_*.bed; do
    time (bgzip -f "$f" && tabix -p bed "${f}.gz") && echo "$(date): indexed ${f}.gz"
done

echo "$(date): $SAMPLE $CHROM pileup done"
