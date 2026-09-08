#!/bin/bash
# Submission helper (not itself a qsub job) to regenerate the samples whose
# whole-genome BAMs were built from the truncated W00a alignment bug
# (see md/20260905.progress.md and W00a_align_ubam.sh's fixed skip-check).
# Originally 6 samples: HG00126, NA20762, HG00253, NA19776, HG03784, NA20870.
#
# 2026-09-07: HG03784 fully completed and verified this run (W00b log:
# 10,178,154 reads, quickcheck passed, MM tag check passed on 500/500 reads,
# 184G merged BAM -- consistent with the rest of the cohort's 170-330G range,
# not anomalous) -- DROPPED from SAMPLES below. Do not add it back: this
# script's stale-file check below only looks at whether an output file
# exists, not whether it's the correct post-fix version, so re-including
# HG03784 here would make the script demand deletion of its now-correct,
# already-verified, ~16hr-to-regenerate merged BAM. Its downstream stages
# (W01/W01g/W02/W03) are unaffected since those still run as full 1-18
# arrays and will pick it up via their own existing skip-checks.
#
# The other 5 are each partially through W00a as of 2026-09-07 (see
# md/20260907_phasing_qc_review.md or ask Claude to re-diagnose via
# scripts/wg/fix/ubam_alignment_status.sh) -- some uBAMs done and marked
# with a .done marker in scratch, some killed (OOM -- see the W00a memory
# bump below) and never retried, some interrupted by a mid-run cancellation.
# None of this requires re-touching genozip/S3 recovery -- purely a matter
# of letting W00a finish the remaining uBAM alignments.
#
# Raw uBAM recovery (W00c) already completed for these samples on 2026-09-05
# -- this script covers everything downstream: W00a (align, now fixed) ->
# W00b (merge) -> W01_phase_per_chrom (haplotag) -> W01g (merge haplotagged)
# -> W02 (modkit extract) -> W03 (modkit pileup).
#
# The scratch alignment dir (/u/project/cluo_scratch/.../wg_align) is NOT
# fully purged -- as of 2026-09-07 it still holds several already-completed,
# .done-marked per-uBAM files from earlier attempts (confirmed for NA20762,
# HG00253, NA20870); W00a's skip-check only trusts a file paired with its
# own .done marker, so those will be skipped and not wastefully redone.
# But W00b/W01g/W02/W03 all skip on bare file existence, not content, so
# their STALE (pre-fix) outputs for these samples must be removed by hand
# before resubmitting, or every stage will just silently keep the old
# truncated data. Per project convention this script does NOT delete
# anything itself -- it checks for stale files and prints the exact rm
# commands, then refuses to submit until you've run them.
#
# Usage: bash scripts/wg/fix/resubmit_6_truncated_samples.sh

set -euo pipefail
cd /u/project/cluo/terencew/claude/project_ideas/asm_lr

SAMPLES=(HG00126 NA20762 HG00253 NA19776 NA20870)

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

# This SGE build rejects multi-range -t specs ("only allows one range
# specification"), so restricted task-ID lists (28-30:1,34-44:1,... etc.)
# don't work here. Instead, submit each stage's own FULL default array
# (no -t override) and rely on each stage's now-correct skip-check: the
# 12-49 unaffected tasks per stage just skip almost instantly, and only
# the 6 samples whose stale output was just removed actually redo.

echo "Submitting W00a_align_ubam.sh (full array, skip-checked)"
JID_A=$(qsub scripts/wg/W00a_align_ubam.sh | grep -oE '[0-9]+' | head -1)
echo "  -> job $JID_A"

echo "Submitting W00b_merge_wg.sh (-hold_jid $JID_A, full array)"
JID_B=$(qsub -hold_jid "$JID_A" scripts/wg/W00b_merge_wg.sh | grep -oE '[0-9]+' | head -1)
echo "  -> job $JID_B"

echo "Submitting W01_phase_per_chrom.sh (-hold_jid $JID_B, full array)"
JID_C=$(qsub -hold_jid "$JID_B" scripts/wg/W01_phase_per_chrom.sh | grep -oE '[0-9]+' | head -1)
echo "  -> job $JID_C"

echo "Submitting W01g_merge_haplotagged_wg.sh (-hold_jid $JID_C, full array)"
JID_D=$(qsub -hold_jid "$JID_C" scripts/wg/W01g_merge_haplotagged_wg.sh | grep -oE '[0-9]+' | head -1)
echo "  -> job $JID_D"

echo "Submitting W02_modkit_extract_per_chrom.sh (-hold_jid $JID_D, full array)"
JID_E=$(qsub -hold_jid "$JID_D" scripts/wg/W02_modkit_extract_per_chrom.sh | grep -oE '[0-9]+' | head -1)
echo "  -> job $JID_E"

echo "Submitting W03_modkit_pileup_per_chrom.sh (-hold_jid $JID_E, full array)"
JID_F=$(qsub -hold_jid "$JID_E" scripts/wg/W03_modkit_pileup_per_chrom.sh | grep -oE '[0-9]+' | head -1)
echo "  -> job $JID_F"

cat <<EOF

Chained: $JID_A (align) -> $JID_B (merge) -> $JID_C (phase/haplotag)
       -> $JID_D (merge haplotagged) -> $JID_E (modkit extract) -> $JID_F (pileup)

Once $JID_F clears, these 6 samples' pileup/ output is fresh and their
G02/G03 (chr1/chr15) beta-binomial results can be regenerated the same way
as the other 8 good donors (they are NOT currently in txt/samples/good_donors.txt
-- add them there, or run G02/G03 directly against these 6, once ready).
EOF
