# Outstanding items, 2026-09-07 — full punch list

Written at the end of a long day's session covering beta-binomial caller calibration,
truncated-BAM recovery, a finalized locus list, and chrX/XCI scoping. Organized by area, roughly
in priority order within each.

## Cohort completeness (blocking some downstream numbers)

1. **5 of 6 originally-truncated samples still incomplete**: HG00126, NA20762, HG00253,
   NA19776, NA20870. `scripts/wg/fix/ubam_alignment_status.sh` gives current per-uBAM status.
   `W00a_align_ubam.sh`'s OOM fix (h_data 2G→3G) applied but not yet exercised at scale.
2. **NA20762 and NA20870 have no beta-binomial results at all yet** (blocked on #1) — 2 of the
   original 6 flagged donors are still totally unevaluated under the new caller.
3. **mosdepth is stale** for all 5 incomplete samples, and for HG03784 (whose merge just
   finished but hasn't been re-measured) — rerun once W00a→W03 clears, per the "triple-check
   coverage matches expected" ask from earlier today.
4. **Coverage manifest** (`tables/phasing_qc/coverage_expected_vs_measured.tsv`) only cross-checks
   against this project's own recorded values — no independent external source verified.

## Beta-binomial caller calibration (the big open methods question)

5. **Design effect (pseudoreplication) — only one 5Mb pilot window measured, on 2 donors.**
   Real and substantial (both donors), but didn't show the predicted flagged-vs-sane direction.
   Needs: multiple windows across full chromosomes, more donors (all 4 flagged + both zero-hit +
   more of the sane six), ideally weighted toward regions that are actually called significant.
   `github/ont_asm_caller/docs/2026-09-07_next_steps_for_discovery_handoff.md` has the plan.
6. **`estimate_dispersion_trend` (the mean-dependent dispersion fix) not yet wired into
   production calling** — `G02_betabinom_region_chr1.py` still uses the original global
   estimator.
7. **The 2 zero-hit donors (HG00344, NA21144) have no root cause.** Extreme ρ̂ (0.097, 0.116)
   doesn't obviously follow from either documented defect. HG00344's low coverage (51-56x, lowest
   in cohort) is a lead but doesn't extend to NA21144 (79-87x, unremarkable).
8. **CPEL (Jiang et al. 2020) not evaluated** as a possible off-the-shelf alternative to
   continued in-house patching.
9. **Read-level confirmatory pass (`P06_readlevel_confirm.py`) only run for one donor, one
   chromosome (HG00146, chr1), and slow** (a naive per-line Python loop streaming a multi-GB
   `modkit extract` file — tens of minutes for one (sample, chrom)). Needs (a) a faster
   implementation before it can scale to 6 donors × 22 chromosomes, (b) actually running that
   full scope, (c) reporting what fraction of the finalized locus list survives.
10. **The new bulk-%mCG gradient (§7 of `md/20260907_phasing_qc_review.md`) and design effect
    (§6) are hypothesized to share a root cause — untested.** Correlate per-donor design effect
    against per-donor bulk methylation rate once #5 has more data points.

## Finalized locus list (today's discovery-arm deliverable)

11. **The ≥2-of-6-donor, 250bp-slop finalization rule is explicitly arbitrary** — not
    sensitivity-tested against other thresholds (≥3 donors? different slop? requiring
    consistent effect-size sign across supporting donors?). Worth a quick robustness check
    before treating 3,821 as a fixed number in any writeup.
12. **chromHMM enrichment not run — no local reference data.** Needs a decision on relevant cell
    type (blood/PBMC-appropriate, matching the tissue this cohort's ONT data comes from) and an
    external download; not attempted here since guessing a source wasn't appropriate.
13. **Gustafson et al. 2024 (1KGP-ONT)**, repeatedly cited as "the closest published analog," has
    never actually been pulled and read in this project — worth doing before leaning on it
    further as a comparison point.

## chrX / XCI

14. **Pipeline scripts written (`scripts/chrX/X00a`, `X00b`, `X01`) and chrX-only reference
    extracted, but nothing has been run.** `X01_haplotag_chrX.sh`'s `CHRX_VCF_SOURCE` is an
    empty placeholder — needs the TOPMed-imputed chrX phased VCF path.
15. **No pilot scale has been chosen or approved.** Recommended in `scripts/chrX/README.md`:
    start with the 2 confirmed-sane female donors (HG00146, NA21110), not all 8 females or all
    18 samples, given every uBAM realignment is a genuine (if smaller-reference) compute cost.
16. **XCI as a design-effect/read-level-test positive control (the mechanistic tie-in noted in
    `scripts/chrX/README.md`) is conceptual only — not executed.**

## Phasing QC (yesterday's still-open thread, now secondary priority)

17. **P01 (dipcall)/P02 (whatshap phase w/ dipcall VCF), jobs 14697093/14697095, still queued**,
    no results yet.
18. **`whatshap stats` on the panel-phased VCF was confirmed uninformative** (population phasing
    is one giant block per chromosome by construction) — will only become useful once P02's
    read-based-phased output exists; not yet checked.

## Process / documentation

19. **Both git repos are currently in sync** (`asm_lr`@`5ab92d9`, `ont_asm_caller`@`55cce16`) —
    keep it that way; the `CLAUDE.md` git-fetch guardrail exists now specifically to catch this
    going forward, in both repos being touched by parallel sessions.
20. **`doc/20260907.results.docx` is a generated artifact** (via `pandoc` from
    `md/20260907.results.summary.md`) — not committed to either git repo, matching their
    existing convention of keeping heavy/binary/generated files out of the trimmed mirrors.
    Regenerate rather than hand-edit if the source md changes.
