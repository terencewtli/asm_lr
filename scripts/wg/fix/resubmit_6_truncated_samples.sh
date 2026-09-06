#!/bin/bash
# Submission helper (not itself a qsub job) to regenerate the 6 samples whose
# whole-genome BAMs are still built from the truncated W00a alignment bug
# (see md/20260905.progress.md and W00a_align_ubam.sh's fixed skip-check):
#   HG00126  NA20762  HG00253  NA19776  HG03784  NA20870
#
# Raw uBAM recovery (W00c) already completed for these 6 on 2026-09-05 -- this
# script covers everything downstream: W00a (align, now fixed) -> W00b (merge)
# -> W01_phase_per_chrom (haplotag) -> W01g (merge haplotagged) -> W02 (modkit
# extract) -> W03 (modkit pileup).
#
# The scratch alignment dir (/u/project/cluo_scratch/.../wg_align) is already
# empty (scratch appears to purge automatically) -- nothing to clean there.
# But W00b/W01g/W02/W03 all skip on bare file existence, not content, so their
# STALE (pre-fix) outputs for these 6 samples must be removed by hand before
# resubmitting, or every stage will just silently keep the old truncated data.
# Per project convention this script does NOT delete anything itself -- it
# checks for stale files and prints the exact rm commands, then refuses to
# submit until you've run them.
#
# Usage: bash scripts/wg/fix/resubmit_6_truncated_samples.sh

set -euo pipefail
cd /u/project/cluo/terencew/claude/project_ideas/asm_lr

SAMPLES=(HG00126 NA20762 HG00253 NA19776 HG03784 NA20870)

STALE=()
for S in "${SAMPLES[@]}"; do
    for f in bam/merged_bam/${S}_wg.bam bam/merged_bam/${S}_wg.bam.bai \
             bam/haplotagged/${S}_wg_haplotagged.bam bam/haplotagged/${S}_wg_haplotagged.bam.bai; do
        [ -e "$f" ] && STALE+=("$f")
    done
    for f in modkit/${S}_chr*_modkit_extract.tsv.gz pileup/${S}_chr*_hp_pileup_*.bed.gz; do
        [ -e "$f" ] && STALE+=("$f")
    done
done

if [ "${#STALE[@]}" -gt 0 ]; then
    echo "Found ${#STALE[@]} stale (pre-fix) output files for these 6 samples."
    echo "W00b/W01g/W02/W03 all skip on bare existence, so these must be removed"
    echo "first or every stage below will silently keep the old truncated data."
    echo
    echo "Run this yourself, then re-run this script:"
    echo
    echo "  rm -f \\"
    for S in "${SAMPLES[@]}"; do
        echo "    bam/merged_bam/${S}_wg.bam bam/merged_bam/${S}_wg.bam.bai \\"
        echo "    bam/haplotagged/${S}_wg_haplotagged.bam bam/haplotagged/${S}_wg_haplotagged.bam.bai \\"
        echo "    modkit/${S}_chr*_modkit_extract.tsv.gz \\"
        echo "    pileup/${S}_chr*_hp_pileup_*.bed.gz \\"
    done
    echo
    exit 1
fi

echo "No stale output files found for these 6 samples -- proceeding."

# task IDs, confirmed directly against each stage's own manifest/sample list
W00A_POOLS="28-30:1,34-44:1,48-50:1,67-72:1"   # csv/meta/ubam_manifest.csv rows
IDX_18="7,9,10,12,17,18"                        # final_samples.txt line numbers (W00b/W01/W01g)
W0203_POOLS="133-154:1,177-220:1,243-264:1,353-396:1"  # sample_chrom_tasks.txt rows

echo "Submitting W00a_align_ubam.sh (-t $W00A_POOLS)"
JID_A=$(qsub -t $W00A_POOLS scripts/wg/W00a_align_ubam.sh | grep -oE '[0-9]+' | head -1)
echo "  -> job $JID_A"

echo "Submitting W00b_merge_wg.sh (-hold_jid $JID_A, -t $IDX_18)"
JID_B=$(qsub -hold_jid "$JID_A" -t $IDX_18 scripts/wg/W00b_merge_wg.sh | grep -oE '[0-9]+' | head -1)
echo "  -> job $JID_B"

echo "Submitting W01_phase_per_chrom.sh (-hold_jid $JID_B, -t $IDX_18)"
JID_C=$(qsub -hold_jid "$JID_B" -t $IDX_18 scripts/wg/W01_phase_per_chrom.sh | grep -oE '[0-9]+' | head -1)
echo "  -> job $JID_C"

echo "Submitting W01g_merge_haplotagged_wg.sh (-hold_jid $JID_C, -t $IDX_18)"
JID_D=$(qsub -hold_jid "$JID_C" -t $IDX_18 scripts/wg/W01g_merge_haplotagged_wg.sh | grep -oE '[0-9]+' | head -1)
echo "  -> job $JID_D"

echo "Submitting W02_modkit_extract_per_chrom.sh (-hold_jid $JID_D, -t $W0203_POOLS)"
JID_E=$(qsub -hold_jid "$JID_D" -t $W0203_POOLS scripts/wg/W02_modkit_extract_per_chrom.sh | grep -oE '[0-9]+' | head -1)
echo "  -> job $JID_E"

echo "Submitting W03_modkit_pileup_per_chrom.sh (-hold_jid $JID_E, -t $W0203_POOLS)"
JID_F=$(qsub -hold_jid "$JID_E" -t $W0203_POOLS scripts/wg/W03_modkit_pileup_per_chrom.sh | grep -oE '[0-9]+' | head -1)
echo "  -> job $JID_F"

cat <<EOF

Chained: $JID_A (align) -> $JID_B (merge) -> $JID_C (phase/haplotag)
       -> $JID_D (merge haplotagged) -> $JID_E (modkit extract) -> $JID_F (pileup)

Once $JID_F clears, these 6 samples' pileup/ output is fresh and their
G02/G03 (chr1/chr15) beta-binomial results can be regenerated the same way
as the other 8 good donors (they are NOT currently in txt/samples/good_donors.txt
-- add them there, or run G02/G03 directly against these 6, once ready).
EOF
