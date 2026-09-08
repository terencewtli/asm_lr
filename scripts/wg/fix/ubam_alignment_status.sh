#!/bin/bash
# Diagnostic (not a qsub job): for each given sample, reports how many of its
# uBAMs (per csv/meta/ubam_manifest.csv) have actually finished W00a alignment
# -- i.e. have a .done marker in scratch, which is the only thing W00a's own
# skip-check trusts (a bare .bam file without .done may be a killed/partial
# run, see W00a_align_ubam.sh's skip-check comment).
#
# Written 2026-09-07 after discovering that only 1 of the 6 originally-
# truncated samples (HG03784) had actually made it all the way through W00a
# + W00b on the first resubmit attempt -- the other 5 were each partially
# done (some uBAMs finished and .done-marked, some OOM-killed and never
# retried, some interrupted by that day's job cancellation). This makes that
# per-sample bookkeeping a single command instead of manual log archaeology.
#
# Usage: bash scripts/wg/fix/ubam_alignment_status.sh [SAMPLE ...]
#   (defaults to the 5 remaining samples from resubmit_6_truncated_samples.sh
#   if no samples are given)

set -euo pipefail
cd /u/project/cluo/terencew/claude/project_ideas/asm_lr

MANIFEST=csv/meta/ubam_manifest.csv
SCRATCH=/u/project/cluo_scratch/terencew/claude/asm_lr/wg_align

SAMPLES=("$@")
if [ ${#SAMPLES[@]} -eq 0 ]; then
    SAMPLES=(HG00126 NA20762 HG00253 NA19776 NA20870)
fi

for S in "${SAMPLES[@]}"; do
    echo "=== $S ==="
    MERGED=bam/merged_bam/${S}_wg.bam
    if [ -f "$MERGED" ]; then
        echo "  merged_bam already exists ($(du -sh "$MERGED" | cut -f1)) -- W00a/W00b done, nothing to check"
        continue
    fi

    rows=$(awk -F, -v s="$S" '$2==s {print $1","$3}' "$MANIFEST")
    n_total=0
    n_done=0
    n_missing=0
    n_partial_no_done=0
    while IFS=, read -r task_id ubam_path; do
        [ -z "$task_id" ] && continue
        n_total=$((n_total+1))
        bname=$(basename "$ubam_path" .bam)
        aligned="$SCRATCH/${S}_${bname}_wg.bam"
        done_marker="${aligned}.done"
        if [ -f "$done_marker" ]; then
            n_done=$((n_done+1))
        elif [ -f "$aligned" ]; then
            n_partial_no_done=$((n_partial_no_done+1))
            echo "  task $task_id: partial file, no .done marker -- will be redone from scratch ($(du -sh "$aligned" 2>/dev/null | cut -f1))"
        else
            n_missing=$((n_missing+1))
            echo "  task $task_id: no output at all yet -- needs a fresh W00a run"
        fi
    done <<< "$rows"
    echo "  -> $n_done/$n_total uBAMs done, $n_partial_no_done partial (will be redone), $n_missing not started"
done
