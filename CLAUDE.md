# asm_lr — read README.md first

This project has a `README.md` in this directory. Read it before making changes here,
especially the "Cohort status" section at the top.

**One-line version:** the active whole-genome pipeline runs on **18 donors**
(`txt/samples/final_samples.txt`). That is the full, intentional working cohort — not a
partial or incomplete subset of an original 30-donor plan. The 12 excluded donors
(`txt/samples/excluded_samples.txt`) use different ONT flow cell chemistry (R9.4.1/Guppy
vs. the 18's R10.4.1/Dorado) and are deliberately held out as a separate stratum — see
`md/basecaller_cohort_split.md`.

**Do not write scripts to bring the excluded 12 back into the primary pipeline** without
reading that file first — this has happened before (2026-08-21), was based on a mistaken
"incomplete download" assumption, and was reverted. Seeing "30" referenced in
`scripts/download/D02_*` or `D03_*`, or in `csv/meta/cohort_downloads.tsv`, does not mean
the active cohort should be 30 — those deliberately keep reference data for all 30 in case
of future expansion.

## Check for other sessions' work before trusting local git state

This project spans two git repos (`github/asm_lr` and `github/ont_asm_caller`), and the user
runs multiple Claude Code sessions in parallel against them — one session's local working
copy of a repo can silently fall behind another session's pushed commits. This already
happened once (2026-09-07): a beta-binomial caller calibration fix, pushed to
`github/ont_asm_caller` by a different session, sat unpulled for a full session's worth of
analysis before a `git push` rejection surfaced it.

**At the start of any session that will read or modify code/docs in either repo**, run
`git fetch && git log --oneline HEAD..origin/master` in that repo before trusting its
current state — not just when you're about to push. If it shows anything, read it before
proceeding; it may change your diagnosis.
