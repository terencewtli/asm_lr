# Beta-binomial re-check of the 4 flagged donors, and why phasing QC needs to happen before more ASM-caller work

**Written 2026-09-07.** Two things: (1) the 4-of-6 previously-Fisher-flagged donors that already have
beta-binomial `chr1` results, checked against the deCODE 0.51% benchmark; (2) a stock-take of what
phasing QC we actually have (short answer: none, for this cohort) and a concrete, non-simulated plan
to get it, prioritized ahead of further ASM-caller iteration. Also nails down coverage/platform,
since it's easy to get wrong: **this is ONT, not PacBio HiFi, and yes, several samples really are
close to 80-110x.**

---

## 1. Does beta-binomial rescue the 4 previously-flagged donors?

Of the original 6 samples flagged anomalous under the retired per-CpG Fisher pipeline
(`NA18508, NA20762, HG01167, NA21093, NA18620, NA20870`), 4 now have beta-binomial region-level
results on chr1 (`G05_betabinom_region_other4_chr1.sh`, ran 2026-09-06). The other 2
(`NA20762`, `NA20870`) are the truncated-BAM-bug samples currently being reprocessed
(`scripts/wg/fix/resubmit_6_truncated_samples.sh`) and have no beta-binomial output yet.

| Donor | Fisher pct sig (old) | Beta-binomial pct sig, chr1 (new) | rho_hat | Verdict |
|---|---|---|---|---|
| NA18508 | 9.84% | **8.70%** | 0.0098 | Still ~17x the deCODE 0.51% benchmark. Unresolved. |
| HG01167 | 3.59% | **5.44%** | 0.0053 | Got *worse* under the new caller. |
| NA21093 | 2.16% | **1.81%** | 0.0065 | Mild improvement, still ~3.5x benchmark. |
| NA18620 | 1.87% | **1.85%** | 0.0065 | Unchanged, still ~3.5x benchmark. |

Compare to the 6 confirmed-sane donors run genome-wide: rho_hat 0.001-0.018, chr1 rates 0.12-0.53%.
**Conclusion: switching statistical methods does not rescue this batch.** That's informative in its
own right — it argues the problem is upstream of the statistical test (data or phasing), not an
artifact of Fisher's known FDR blowup.

**New, unrelated finding surfaced by the same batch of runs:** 2 of the original **8 "good"**
donors — `HG00344` and `NA21144` — show a different failure mode under beta-binomial: rho_hat
0.097 and 0.116 (vs. 0.001-0.018 for the other 6) and **zero significant regions on chr1**. Excluded
from the genome-wide run (`G04_betabinom_region_genome.sh`) pending root-cause. Not yet investigated.

Before this doc, none of the above comparison existed anywhere written down — only the *plan* to
re-run these 4 was documented (`md/20260905.progress.md` §14 item 6, and the `G05` script header).
The numbers above are read directly out of `tables/dasm_betabinom/*/summary_chr1.json`.

---

## 2. Coverage and platform — double-checked, stated plainly

**This is Oxford Nanopore (ONT) long-read sequencing. It is not PacBio HiFi.** `csv/meta/cohort_downloads.tsv`'s
`ont_basecaller` column reads `Dorado 0.6.0 5mCG+5hmCG` (Dorado is ONT's basecaller) for all 18
active donors, and every S3 path is under `raw_data/nanopore/dorado0.6.0_.../`. HiFi is used
elsewhere in this project only as the platform HPRC used to build the **diploid genome assemblies**
(`assemblies/{sample}_hap{1,2}.fa.gz`) that serve as phasing truth — the methylation reads
themselves are 100% ONT.

And yes — several of these samples really are ~80-110x, independently confirmed two ways:

| Sample | Manifest `ont_total_cov_x` | Measured (mosdepth, autosomal total) |
|---|---|---|
| NA18508 | 80 | 88.92x |
| NA18620 | 79 | 88.93x |
| NA21093 | 97 | 109.01x |
| HG01167 | 67 | 75.63x |
| HG00146 (good) | 76 | 76.80x |
| NA19700 (good) | 62 | 69.61x |

(Elsewhere in the full 30-donor HPRC manifest, outside our 18-donor cohort, coverage goes even
higher — HG01361 at 140x, HG00438 at 137x — so 80-110x is not unusual for this dataset.)

This matters for the anomaly: coverage is **not confounded with anomaly status**. NA18620/NA21093
are among the deepest samples in the cohort and are still elevated; HG00146 is comparably deep and
is one of the "good" 6. This was already established in the phase-2 QC writeup
(`md/20260905.progress.md`) via the NA20762/HG00126 near-identical-coverage comparison, and it holds
here too — depth alone still doesn't explain the pattern.

---

## 3. Is this still statistical (population-panel) phasing, or read-based?

**Still statistical, cohort-wide, no exceptions.** Production haplotagging
(`scripts/wg/W01_phase_per_chrom.sh`) runs `whatshap haplotag` against the 1000 Genomes 3202-sample
phased panel (`vcf/read_based/{sample}_{chrom}_het_snps.vcf.gz`, sourced from the EBI 1000G
high-coverage phased-panel release). `whatshap haplotag` only **assigns** ONT reads to one of the two
already-fixed panel haplotypes — it never uses the ONT reads' own multi-SNP linkage to inform or
correct the phase itself.

The hybrid statistical+read-based alternative (`whatshap phase`, `scripts/wg/W01b_phase_readbased_per_chrom.sh`)
exists and is described as ready in its own header, but **has been run exactly once**: one
(sample, chromosome) pair, `HG00146 chr21`, on 2026-09-05
(`logs/W01b/W01b_phase_readbased_per_chrom.14683548.131`). It is not wired into the production
pipeline, and it has not been run for any of the 4 flagged donors.

**So: yes, every donor in the active 18-donor cohort — flagged and "good" alike — is currently
haplotagged against pure population-panel phase, with zero read-based correction.**

---

## 4. Could population-panel phasing explain the anomaly, and how would that even work?

Plausible, not yet confirmed — no switch-error measurement exists for this cohort at all (see §5).
But the mechanism is real and worth spelling out, because "switch errors → false positives" is not
obvious on its own:

**What a switch error does.** A switch error means the phasing flips which physical haplotype is
called "HP1" vs "HP2" at some point along a chromosome — usually where the statistical phasing
algorithm had weak local LD support (sparse heterozygosity, or a region poorly represented in the
reference panel). Downstream of a switch, reads that are truly all from one parental chromosome get
split across both HP labels, and reads from the two true parental chromosomes get merged into one
label.

**Why that manufactures ASM signal rather than just adding noise.** If methylation were pooled at
random with respect to the switch, this would just wash out real signal — average both groups
toward the same value and reduce power, not manufacture false positives. But it doesn't stay random
whenever something *else* correlates with genotype near the switch point: reference-mapping bias,
segmental duplication or repeat-driven allele-specific alignment artifacts, or systematic quality
differences between reads carrying the reference vs. alternate allele. Any of those turns a scrambled
"HP1 vs HP2" grouping into a grouping that is still non-random with respect to some other variable —
which a statistical test (Fisher or beta-binomial, either one) can then call significant. This is the
same class of confound as reference-bias in allele-specific *expression* (the WASP/MBASED literature
already cited in `md/20260905.progress.md`), just not previously ported to long-read methylation ASM.

**Why beta-binomial region-pooling doesn't fix this.** Region pooling helps with per-CpG sparsity and
overdispersion, but a phase-block switch inside a tested region mixes reads from a "correctly phased"
sub-block with reads from a "post-switch" sub-block *within the same region* — pooling more CpGs
under a broken haplotype label doesn't un-break the label. This is consistent with what we actually
saw in §1: the anomaly survived the caller switch.

**Why this could plausibly be donor-specific rather than cohort-wide.** Statistical phasing accuracy
against a population panel is known to vary with how well the panel represents that individual's
ancestry and local haplotype diversity — already flagged in `md/20260905.progress.md` as the
"most promising open lead," specifically for the AFR sample (NA18508). That's a real, literature-backed
axis (population-panel phasing is measurably worse for underrepresented ancestries) that could produce
exactly the kind of per-sample, non-coverage-tracking variance we're seeing. **But this is still a
hypothesis** — nobody has actually measured a switch-error rate for any of these samples, flagged or
good, so we don't know whether the "good" 6-8 donors are just as switch-error-prone but happen not to
trip the statistical test as often, or whether they're genuinely better phased.

---

## 5. What phasing QC do we actually have? (Answer: none, for this cohort.)

Checked directly — no `whatshap stats` output, no phase block N50/NG50, no switch-error rate, no
`whatshap compare` output exists anywhere in the project for any of the 18 active donors. Every
mention of block N50 / switch error in `md/` is a **recommendation for future work**
(`md/20260903_qc_review.md` §5.2-5.3), never something that was run. This is worth saying plainly:
**we've spent weeks iterating on the ASM statistical caller without ever running the phasing QC step
that would tell us whether the input haplotype labels are trustworthy in the first place.** That's
backwards, and it should be fixed before more caller work, since a phasing problem invalidates ASM
calls regardless of which downstream test is used.

### 5a. Phase block size / N50 — no truth data needed, can run today

`whatshap stats` (v2.2, available in the `whatshap-env` conda env) computes per-chromosome block
count, phase block N50/NG50, and fraction of variants phased directly from a single phased VCF —
no comparison file needed:

```
whatshap stats --tsv stats.tsv --block-list blocks.txt --chr-lengths hg38.chrom.sizes \
    vcf/read_based/{sample}_{chrom}_het_snps.vcf.gz
```

This can be run right now, today, against the existing 1000G-panel phased VCFs for all 18 donors —
cheapest possible QC step, zero new upstream infra needed. It won't tell you switch error, but a
donor with unusually short/fragmented phase blocks compared to the rest of the cohort is an
immediate red flag worth checking before trusting its ASM calls, and it costs essentially nothing to
compute.

### 5b. Switch error, non-simulated — the infrastructure already exists and is unused

Switch error needs an independent ground-truth phasing to compare against — you can't get it from
the panel-phased VCF alone. The obvious ground truth here is the **HPRC diploid assembly** for each
donor (`assemblies/{sample}_hap1.fa.gz` / `_hap2.fa.gz`), since assembly-vs-reference alignment gives
a long-range, sequencing-independent phase that has nothing to do with the 1000G panel or ONT reads.

The pipeline to turn that into a real switch-error number **already exists in this repo and has
simply never been run for the 18-donor cohort**:

1. **Assemblies are already downloaded** for all 4 flagged donors (confirmed:
   `assemblies/{NA18508,HG01167,NA21093,NA18620}_hap{1,2}.fa.gz` all present) — and per
   `scripts/download/D03_download_hprc_assemblies.sh`, for the full 18-donor set.
2. **`scripts/download/D04_run_dipcall.sh`** is a ready-to-submit SGE array job, scoped exactly to
   `-t 1-18` (i.e. already targeted at the active cohort, not the excluded 12), that runs
   `minimap2 asm5` + `paftools.js call` + `dipcall-aux.js merge` to turn each donor's two assembly
   haplotypes into `vcf/read_based/{sample}_dip.vcf.gz` — an assembly-backed truth-phased VCF. **This
   has not been run for any of the 18 donors**; the only `*_dip.vcf.gz` files that exist
   (`HG00438, HG00621, HG01258, HG01361, HG01891, HG01928, HG02572, HG03453`) are HPRC's own
   precomputed Year-1 dipcall VCFs for the *excluded* Guppy cohort, downloaded by a different script
   (`D01_download_hprc_assembly_vcf.sh`) — not something this project generated, and not applicable
   to our 18 Dorado donors.
3. Once that truth VCF exists per donor, **`whatshap compare`** directly computes a real,
   non-simulated switch error rate against it:

```
whatshap compare --names panel,truth --tsv-pairwise switch_stats.tsv \
    --switch-error-bed switches.bed \
    vcf/read_based/{sample}_{chrom}_het_snps.vcf.gz \
    vcf/read_based/{sample}_dip.vcf.gz
```

This gives, per sample per chromosome: switch error rate, Hamming distance, and a BED file of exactly
where the switches happen — genuine ground truth, no simulation involved. The same comparison can
later be repeated with the `whatshap phase` (read-based) output in place of the panel VCF, which
would directly answer whether read-based phasing measurably reduces switch error for this cohort —
the comparison `md/20260903_qc_review.md` §5.2 already proposed as a figure.

### 5c. Does this nominate NA18508/HG01167/NA21093/NA18620 as bad, right now?

**We don't know yet — nothing in §5a or §5b has actually been run.** That's the point of this
section: the phasing QC that would answer "are these donors' haplotype labels trustworthy" doesn't
exist yet for any donor, flagged or good. Recommended order of operations, before any further
ASM-caller iteration (including re-running Fisher/beta-binomial on more chromosomes, and before
folding the 2 recovering truncated-BAM samples back in):

1. Run `whatshap stats` (§5a) on all 18 donors' existing panel-phased VCFs — cheap, no new
   dependencies, gives block N50/NG50/phased-fraction today.
2. Submit `D04_run_dipcall.sh` (§5b) for the 18-donor cohort — assemblies are already staged, this
   is a straight resource/time cost, not a design decision.
3. Run `whatshap compare` per sample per chromosome against the resulting dipcall truth — gives the
   actual switch-error rate.
4. Cross-tabulate switch error / block N50 against the beta-binomial anomaly status (the 4 flagged
   donors, the 2 newly-anomalous HG00344/NA21144, and the 6 confirmed-sane donors) — this is the
   test of the hypothesis in §4, and it's a strictly better diagnostic than anything derived from
   ASM candidate rates themselves, since it doesn't depend on which statistical caller is in fashion.

None of steps 1-3 have been submitted as jobs by this session — they're infrastructure that already
exists in the repo, staged and ready, not yet run.

---

## 6. Later the same day: a cross-session discovery changes the leading hypothesis

While pushing `P04_permutation_null.py` to `github/ont_asm_caller`, the push was rejected —
the remote had a commit (`ccabe0f`, 2026-09-05 21:41, from a **different, parallel Claude Code
session** the user ran with Opus 5) that had never been pulled into this working copy. That
commit is a calibration critique of the exact beta-binomial caller every result in §1-§5 (and
in `md/20260905.progress.md`, `md/20260906.progress.md`) was generated with. Full writeup:
`github/ont_asm_caller/docs/2026-09-05_calibration_critique.md`.

**Defect 1 (major): `region.cluster_cpgs` pools raw counts across CpGs in a window as if each
CpG's reads were independent.** An ONT read spans every CpG in a 500bp window, so 8 CpGs × 15
reads is **15 molecules counted 8 times**, not 120 independent observations. Measured on
simulated null data: type I error is fine at co-methylation=0 (what v1's own benchmark used),
but **2.9× inflated at α=0.05 and 11× inflated at α=0.001** at co-methylation=1, worsening with
region size — exactly the small-p regime genome-wide BH-FDR operates in.

**Defect 2: the dispersion estimator (`estimate_dispersion`) is unbiased only at μ=0.5** (again,
what v1's own benchmark used) and is **2-3.5× too low on a realistic bimodal methylome** —
anti-conservative, more false positives, worse the further a locus sits from μ=0.5.

**Why this reframes everything above:** neither defect depends on phasing accuracy. Defect 1's
magnitude is a function of a *donor's own within-read co-methylation structure* (how "blocky"
their real methylation landscape is) and region size — not switch error rate. This directly
answers the objection that "the 6 good donors use the identical phasing method and are fine" —
under this hypothesis they don't need to be phased any better, they just need to have less
extreme co-methylation/region-size behavior for their real ASM signal, or plainly less
degenerate methylation distributions.

**Confirmed real on our own data, not just simulated.** Implemented
`readlevel.design_effect()` from the newly-pulled commit against real `modkit extract` output
(`scripts/wg/phasing_qc/P05_design_effect_real_data.py` — HP tags recovered separately from the
haplotagged BAM since modkit extract itself isn't haplotype-split; design effect measured on
the production caller's own region boundaries, restricted to chr1:1-5,000,000 for tractability,
not genome-wide):

| Donor | n regions | median DE | mean DE | % DE>2 | % DE>5 |
|---|---|---|---|---|---|
| HG00146 (sane) | 300 | 1.39 | 2.09 | 32.0% | 5.7% |
| NA18508 (flagged) | 300 | 1.13 | 1.52 | 18.3% | 1.0% |

(pilot: chr1:1-5,000,000, the production caller's own region boundaries, HP tags recovered from
the haplotagged BAM since modkit extract isn't itself haplotype-split;
`scripts/wg/phasing_qc/P05_design_effect_real_data.py`,
`tables/phasing_qc/design_effect_{sample}_chr1_pilot.tsv`)

**Design effect is confirmed real and substantial on our actual data — both donors sit well
above DE=1** (the "CpGs are independent" assumption v1's caller and its benchmark rely on), so
Defect 1 is not just a theoretical concern here. **But this specific pilot window does not show
the predicted direction**: HG00146 (sane) has the *higher* median/mean design effect of the two,
not NA18508 (flagged) — the opposite of what the "flagged donors have worse co-methylation
structure" hypothesis predicts. Read this as inconclusive, not disconfirming: the window is the
first 5Mb of chr1, subtelomeric, with unusually different HP-tagged read counts between the two
donors in this specific region (18,585 vs. 7,079) — not necessarily representative of the rest
of the chromosome, where the real candidate-ASM-rate divergence actually lives. **Next step:
sample multiple windows spread across the whole chromosome (not just the telomere-adjacent
start) before concluding anything about whether donor-to-donor DE differences track anomaly
status** — this pilot establishes the measurement pipeline works and that DE is real, not yet
the donor-comparison answer.

**Important, humbling correction to §1's permutation-null result:** that bootstrap simulated
fresh BetaBinomial draws directly at the already-*pooled* per-region depths (n1, n2) — i.e. it
inherited the same "N counted once per CpG" assumption Defect 1 says is wrong, and therefore
**cannot detect Defect 1's failure mode** — it can only show the model's own math is
self-consistent *given* that assumption. "0% null for every donor" is still true and still rules
out gross dispersion-fitting failure, but it is no longer a clean "therefore the caller isn't at
fault" verdict.

**Correcting for this when calling ASM — concrete options, in the package already:**
1. **`readlevel.test_region_reads` / `test_region_perm`** (added in the same commit): summarise
   each molecule to one methylation fraction over the region and compare haplotypes'
   read-level distributions directly. Needs no dispersion parameter, handles complete
   separation (imprinted-DMR-strength signal) exactly via exact randomisation, and is
   **correctly calibrated at every co-methylation level tested** (type I 0.042/0.033/0.009 vs.
   nominal 0.05) — at the cost of real power loss when the true design effect is high (which is
   *correct* conservatism, not a bug: at DE≈5.5 there genuinely are ~5× fewer independent
   observations than reads).
2. `dispersion.estimate_dispersion_trend`: mean-dependent ρ(μ), simulation-calibrated instead of
   assuming E[χ²]=1. Better than the global estimator but still overshoots (2-3×) at the
   boundary — a safe (conservative) failure mode, but costs power exactly where imprinting-style
   ASM lives.
3. Use the pooled beta-binomial test only for **screening**, and require read-level
   confirmation (option 1) before calling anything significant — converts the pooled test's
   power advantage into a first pass without inheriting its false-positive blowup.

**Is this "intrinsic to the data" in a way phasing can't fix?** Design effect is a property of
each donor's own methylation processivity/region structure, not of phasing accuracy, so
switching to read-based phasing (P01/P02, §5) won't address it — the two are orthogonal
problems that could both be contributing simultaneously. And design effect is not confined to
the flagged donors: the pilot measurement above finds real design effect **even in a confirmed-
sane donor**, which reframes the open question from "why are 4-6 donors broken" to "how much of
*every* donor's pooled-region ASM call set is inflated, and does that inflation happen to cross
the significance threshold more often for some donors than others" — a matter of degree, not of
some donors having a categorically different data-generating process. That degree-of-inflation
question is exactly what the global bulk-methylation-rate finding (§7) may be tracking.

## 7. Global bulk methylation rate — a new, previously undocumented per-donor gradient

Bulk %mCG (chr1, HP1+HP2 combined) across the 4 flagged + 2 dispersion-broken + 6 sane donors
spans **60.5% (NA21093) to 75.5% (NA19700)** — a ~15-point spread never previously measured.
5 of the 6 problem donors (NA21093, NA21144, HG00344, NA18620, HG01167) occupy the bottom 5 of
12 positions; only NA18508 breaks rank (7th of 12). Table:
`tables/phasing_qc/bulk_methylation_chr1.tsv`.

Plausible unifying mechanism with §6: a donor whose methylation sits closer to the boundaries
(more stably-methylated/unmethylated stretches, lower bulk %mCG) would show both lower bulk
%mCG *and* higher within-read co-methylation (long uniform stretches) — i.e. §6 and §7 may be
two views of the same root cause, not two separate hypotheses. Not yet proven; the design-effect
measurement above, extended across donors and correlated against bulk %mCG, would test this
directly.

## 8. BAM regeneration status (2026-09-07) and script fixes

Of the 6 originally-truncated samples, only **HG03784 fully completed** (verified: 10,178,154
reads, quickcheck + MM-tag checks passed, 184GB — in line with the other 12 successful samples'
170-330GB range). The other 5 are each partially through W00a (some uBAMs `.done`-marked and
complete, some OOM-killed under the old 2GB memory budget and never retried, some interrupted by
a same-day job cancellation). Diagnostic: `scripts/wg/fix/ubam_alignment_status.sh`.

Two fixes applied:
- `scripts/wg/W00a_align_ubam.sh`: `h_data` 2G→3G (16GB→24GB total). One surviving log showed
  minimap2 alone hitting Peak RSS 16.57GB against the old 16GB ceiling — a plain OOM
  resource-sizing bug, unrelated to data quality, likely explains the repeated silent "Killed"
  tasks across NA20762/HG00253/NA19776/NA20870.
- `scripts/wg/fix/resubmit_6_truncated_samples.sh`: dropped HG03784 from `SAMPLES` — its
  stale-file check only tests existence, not correct-version, so leaving HG03784 in would have
  demanded deleting its freshly-completed, ~16hr-to-regenerate merged BAM.

Coverage manifest (self-reported vs. mosdepth-measured, all 18 donors):
`tables/phasing_qc/coverage_expected_vs_measured.tsv`. Note: mosdepth for the 5 still-incomplete
samples (and HG03784, whose merge just finished) is stale and needs rerunning once W00a→W03
clears.

## 9. Checked, not the right paper

`reference/zhuo_2026/zhuo_2026.pdf` (Zhuo et al. 2026, *Genome Research* — "Characterizing
cytosine methylation of polymorphic transposable element insertions using the human pangenome
resources," same Ting Wang lab already flagged as the epigenome-portal incumbent) is a
TE-methylation study using HiFi+WGBS+ONT and HiPhase (HiFi) phasing — not an ASM-calling paper,
no mention of switch error/dispersion/per-sample ASM heterogeneity. Not a useful comparison for
this project's donor-anomaly question. Gustafson et al. 2024 (1KGP-ONT) remains the closer
analog if this needs checking against another group's experience.

---

## 10. Proximal/distal breakdown was wrong — found and fixed a real bug, not just staleness

The 91.9% proximal / 7.4% distal / 0.7% very_distal breakdown reported earlier today (and in
`md/20260907.results.summary.md`) is **wrong**, not just built on stale donor coverage (§ above
already fixed that separately). Directly comparing the original `six_donor_sig_regions.tsv`
against a from-scratch recompute of `nearest_snp_dist` (same method: nearest het-SNP position
per `vcf/read_based/{sample}_{chrom}_het_snps.vcf.gz`, same 200bp/5kb thresholds) found:
**HG00146 chr1 matches exactly between old and new (249/249 rows, identical distances)**, but
every *other* chromosome's category values in the original aggregate are wrong — HG00146's
genome-wide old value was 98.1% proximal, but the corrected per-chromosome values (all 6 donors,
all 22 chroms) range smoothly 22-64% proximal with no chromosome anywhere near 98%. The original
per-chromosome parallel computation (`ProcessPoolExecutor`, one job per (sample, chrom)) has a
real bug for the non-chr1 jobs; not root-caused, but the recompute is verified correct against
chr1 as ground truth and produces a smooth, sane per-chromosome pattern everywhere.

**Corrected numbers** (`tables/asm_analysis/six_donor_sig_regions_COMPLETE.tsv`,
`tables/asm_analysis/finalized_asm_loci.tsv`):

| | proximal (≤200bp) | distal (200bp-5kb) | very distal (>5kb) |
|---|---|---|---|
| Per-donor, unreplicated (41,461 sig. regions) | 44.3% | 48.9% | 6.9% |
| Finalized, ≥2-of-6-donor replicated (4,466 loci) | 55.6% | 41.0% | 3.4% |

Two things follow from this:

1. **The corrected unreplicated split (44%/49%/7%) is close to phase-1's own single-donor number
   (33.4% proximal / 66.6% long-read-only, `md/progress.md`)** — much more consistent with the
   "why long-read" thesis than the wrong 91.9% figure suggested, and consistent with the
   intuition that 30-40kb reads should reach plenty of CpGs far from the nearest phaseable SNP.
2. **Requiring ≥2-donor replication measurably shifts the mix toward proximal** (44%→56%
   proximal, 7%→3% very-distal — very-distal roughly halves). Read as a statistical-power
   effect, not a biological one: a locus anchored by a *nearby* het SNP has a shorter, more
   direct single-molecule linkage to the tested CpG (fewer opportunities for something to go
   wrong over the intervening span), so it's more likely to independently clear significance in
   multiple donors with different het-SNP sets and coverage. Distal/very-distal loci — the class
   that's actually novel to long reads — are disproportionately filtered out by a
   replication-based inclusion criterion.

**Implication for the paper's evidence chain**: don't use the ≥2-donor-replicated set as the
basis for the "long reads reveal otherwise-invisible ASM" claim — that claim should rest on
phase-1's already-validated single-donor ascertainment analysis (which isn't subject to this
replication-power bias). The replicated set is a different, narrower, complementary claim
("these specific loci reproduce across individuals") and should be presented as such, with the
proximal-skew called out explicitly rather than left implicit.

`md/20260907.results.summary.md` and `doc/20260907.results.docx` still have the old, wrong
91.9% figure and need regenerating with the corrected numbers/figures.

## 11. Paper scope decision (discussed at length, not just QC)

Project owner's call, after discussion: **one paper, not two, in a ~2 month timeline.**
Biology-primary framing ("what do we gain from long-read vs. short-read ASM detection"), with
the donor-heterogeneity/caller-calibration work included as supporting validation content
("methods baked into main figures"), not as a second, methods-led paper. Explicitly ruled out
of scope for this paper: chrX/XCI, CPEL evaluation, fully root-causing the two zero-hit donors,
population/ancestry-stratified biology claims (n=6, at most n=18, is underpowered for this —
already directly tested: the one same-ancestry pair among the six replicates at a rate
indistinguishable from cross-ancestry pairs).

**What "baked-in methods" needs to be credible, scoped down from a full methods paper**:
(1) a principled, transparent donor-exclusion criterion (the permutation-null test, already
built and run on 6 donors) rather than fully explaining *why* specific donors fail it;
(2) the read-level confirmatory pass (`P06_readlevel_confirm.py`) run genome-wide on the
retained donors (currently only HG00146 chr1: 248/249, 99.6%, confirmed but not yet a strong
test — see below); (3) the positive controls already in hand (imprinting recovery, meQTL
enrichment, effect-size concordance across donors, §12).

**Decided: switch phasing from statistical (1000G panel) to assembly-backed
(dipcall + `whatshap phase`) as the production standard**, not just a QC side-comparison — the
project owner's reasoning (agreed): the panel VCF is structurally worse than using the HPRC
assemblies directly, given they're already downloaded, and it's not a large compute lift.
Re-submitted genome-wide (not chr1-only pilot): `P01_run_dipcall_pilot.sh` (job 14697667,
already genome-wide by construction — dipcall aligns the whole assembly) →
`P02_phase_readbased_dipcall_vcf.sh` (job 14697668, widened from chr1-only to all 22
chromosomes, 132 tasks = 6 donors x 22 chroms, held on the dipcall job). The prior chr1-only
P01/P02 submission (jobs 14697093/14697095) never actually ran — no logs, no output — likely
lost to queue turnover; this is a fresh resubmission, not a continuation.

**Next step once P02 clears**: re-haplotag the 6 sane donors against the dipcall-phased VCF
(replacing the 1000G panel VCF), then re-run W02/W03/beta-binomial calling downstream for those
donors — this is a real pipeline re-execution, not just a QC comparison, since the decision is
to make this the new production source.

## 12. Big next steps, in priority order (superseding the punch list in `md/20260907_outstanding_items.md` for anything paper-scoped)

1. Regenerate `md/20260907.results.summary.md` / `doc/20260907.results.docx` with the corrected
   proximal/distal numbers (§10) and the replication-bias caveat.
2. Scale `P06_readlevel_confirm.py` to run genome-wide across all 6 sane donors — the real test
   of "do the retained calls hold up," and (separately, lower priority) run it once on a flagged
   donor's candidates as the actual test of the pseudoreplication hypothesis.
3. Let P01/P02 (dipcall + genome-wide whatshap phase) finish, then re-haplotag and re-call the
   6 sane donors against the new phased VCF — this becomes the production locus set the paper
   is built on, superseding the panel-phased results throughout.
4. Once (3) lands, redo the finalized-locus-list build, the four discovery checks, and the
   effect-size-consistency notebook against the dipcall-phased calls — the current versions of
   all of these are built on panel-phased data that's about to be superseded.
5. Nail down a defensible, transparent per-donor QC/exclusion criterion (permutation-null-based)
   as the paper's stated basis for which donors are included, rather than presenting an
   unresolved anomaly.

---

## 13. Session close-out (end of day) — results, decisions, and priorities for whoever picks this up next

**Read this section first if resuming.** Today ran long and touched a lot of moving parts across
this file (§1-13); this section is the synthesis. Also relevant, same date: `md/20260907.results.summary.md`
(has stale numbers, not yet regenerated — see §10), `md/20260907_outstanding_items.md` (written mid-day,
partly superseded by the scope decision in §11 below), and `github/ont_asm_caller/docs/2026-09-07_next_steps_for_discovery_handoff.md`
(the methods-arm handoff, in the *other* repo).

### New results that shape the paper's main flavor

1. **Per-donor proximal/distal split is ~50/50, not 90/10.** Two real bugs were found and fixed
   today (§10): a stale intermediate aggregate missing most of NA18959's chromosomes, and a
   genuine computational bug in the original per-chromosome nearest-het-SNP-distance calculation
   (verified: chr1 matches a from-scratch recompute exactly; every other chromosome's original
   values don't). Corrected, per-donor, no replication requirement, no reliance on phase-1
   (deliberately dropped as an evidentiary basis — different basecaller/chemistry than the main
   cohort, and predates this project's whole calling-methodology evolution):

   | Donor | proximal | distal | very distal |
   |---|---|---|---|
   | HG00146 | 35.2% | 55.5% | 9.3% |
   | HG02392 | 46.9% | 46.6% | 6.5% |
   | NA18959 | 46.4% | 45.8% | 7.8% |
   | NA19682 | 48.7% | 46.5% | 4.7% |
   | NA19700 | 56.0% | 41.6% | 2.3% |
   | NA21110 | 51.2% | 43.2% | 5.6% |

   Consistent across all 6 donors independently — not driven by one outlier. **This is now the
   paper's core "you get a lot more loci with long reads" number**: roughly half of each donor's
   own significant ASM calls sit beyond where short-read phasing can reach (proximal ≤200bp is
   phase-1's own short-read-detectable-proxy definition, so this is directly comparable in kind,
   just measured fresh on the modern cohort/pipeline instead of reused from phase-1).

2. **Imprinting recovery, stratified by cross-donor replication count, is a clean monotonic
   dose-response** — the single best new result of the day:

   | n donors supporting locus | n loci | % that are known imprinted DMRs |
   |---|---|---|
   | 1 | 25,425 | 0.10% |
   | 2 | 3,052 | 0.33% |
   | 3 | 859 | 1.40% |
   | 4 | 339 | 5.60% |
   | 5 | 154 | 18.83% |
   | 6 (all donors) | 62 | **53.23%** |

   Loci replicating in all 6 donors are majority known imprinted DMRs; donor-private loci are
   essentially never imprinted (0.10%). This validates the replication-count axis as tracking
   real biology (imprinting is the one truly population-universal ASM mechanism, and it
   dominates exactly the tier where universal replication would predict it should) — a strong,
   clean, headline-figure-worthy result.

3. **The heavy single-donor concentration (85% of merged loci, 25,425/29,891) is mostly NOT an
   ascertainment artifact.** Sampled 500 singleton loci and checked whether the other 5 donors
   even had a *testable* (adequately-covered, regardless of significance) region there:
   **90.8% were testable in all 5 other donors — they just weren't called significant.** Only
   0.2% were genuinely untestable elsewhere. This redirects the explanation toward genuine
   allele-frequency-limited, genotype-driven ASM (a real causal variant needs to be present, not
   just a nearby het SNP) rather than "other donors never had the opportunity to detect it" —
   consistent with, not contradictory to, real biology, though some residual technical noise
   contribution can't be fully excluded.

4. **Effect-size concordance across donors at shared loci** (§ from earlier today,
   `notebooks/asm_analysis/ASM02_effect_size_consistency.ipynb`): 79.6% of |Δ| variance is
   between-locus (locus identity), not between-donor; pairwise concordance r=0.724 (n=10,133
   donor-pairs). Detected *loci* differ donor-to-donor (consistent with #3 above — different
   donors carry different het SNPs/causal variants), but effect *size* at a shared locus is
   consistent.

**Excitement flagged by the project owner**: the ~50/50 proximal/distal split is a bigger,
more compelling gain than the earlier (wrong) 90/10 split suggested, and opens a genuine new
direction — **investigating specific distal loci for mechanism** (chromatin context, gene
proximity, what's actually going on biologically at the long-read-only class) is now an
explicit next priority, not just a QC-validated count. Not yet scoped or started.

### Decisions made today about the paper itself

- **One paper, ~2 months, biology-primary** ("what do we gain from long-read vs. short-read ASM
  detection"), not two papers and not methods-led. Donor-heterogeneity/caller-calibration work
  is supporting validation content, scoped down to: a transparent donor-exclusion criterion +
  the read-level confirmatory pass + the existing positive controls (§11, unchanged).
- **Explicitly out of scope for this paper**: chrX/XCI, CPEL evaluation, fully root-causing the
  two zero-hit donors (HG00344, NA21144), and any population/ancestry-stratified biology claim
  (n=6, or even n=18, is underpowered — directly tested: the one same-ancestry pair among the
  six replicates indistinguishably from cross-ancestry pairs).
- **Phase-1 (HG01258) dropped as the evidentiary basis** for the long-read-gain claim — wrong
  basecaller/chemistry relative to the main cohort, methodologically superseded. Use the fresh
  per-donor breakdown above instead.
- **Production phasing is switching from statistical (1000G panel) to assembly-backed
  (dipcall + `whatshap phase`)** — not just a QC side-comparison anymore, the project owner's
  call given the assemblies are already available and it's not a large compute lift. **Still
  running as of session close**: `P01_run_dipcall_pilot.sh` (job 14697667, genome-wide, still
  `qw`) → `P02_phase_readbased_dipcall_vcf.sh` (job 14697668, 132 tasks = 6 donors × 22 chroms,
  held, not started). Once this clears: re-haplotag the 6 sane donors against the new phased
  VCF, then re-run W02/W03/beta-binomial calling downstream — **every locus-list/discovery
  number in this file and in `md/20260907.results.summary.md` is about to be superseded by
  that re-call**, not final.

### Open discussion points / priorities, not yet resolved

Roughly in priority order for whoever picks this up:

1. **Regenerate everything once the dipcall re-call (above) lands** — finalized locus list,
   the four discovery checks, the effect-size-consistency notebook, and the corrected
   proximal/distal table are all currently built on soon-to-be-superseded panel-phased data.
2. **Investigate specific distal loci for mechanism** (new, per the project owner's excitement
   above) — not scoped yet. Natural next questions: what genes/regulatory elements do they sit
   in, is there a chromatin-state or repeat-content story (phase-1's W10 already found distal
   and proximal loci functionally indistinguishable on that axis for HG01258 — worth checking
   whether that holds on the modern multi-donor set too), any enrichment pattern.
3. **Sensitivity check on the 250bp merge slop** (0bp / 250bp / 500bp) — flagged twice today,
   still not run. Priority raised given how much is now resting on the replication-count
   figures (imprinting stratification, effect-size concordance).
4. **Scale `P06_readlevel_confirm.py` genome-wide** across the 6 sane donors (currently only
   HG00146 chr1: 248/249 confirmed, not yet a strong test — the real test is running it on a
   flagged donor's candidates, still not done either).
5. **"Testable but not significant" trend-direction check** (proposed today, not run): for
   singleton loci, do the other testable-but-non-significant donors at least trend the same
   direction, even without reaching significance? Would further distinguish "real but
   underpowered elsewhere" from "genuinely private" for the singleton-locus question.
6. Literature check: is the heavy single-donor/"erratic ASM" pattern already reported elsewhere
   (the Gutiérrez-Arcelus-style genotype-vs-donor-specific ASM/ASE literature)? Not searched.
7. A defensible, transparent per-donor QC/exclusion criterion (permutation-null-based) still
   needs to be finalized as the paper's stated basis for donor inclusion.

### For a future Claude Code session picking this up

- **Read order**: this section first, then skim §1-12 of this same file for the day's full
  narrative, then `md/20260907_outstanding_items.md` for anything not superseded by the scope
  decision above, then `github/ont_asm_caller/docs/2026-09-07_next_steps_for_discovery_handoff.md`
  if picking up methods-arm work specifically (separate repo, separate session by design — see
  `CLAUDE.md`'s git-fetch guardrail, added today for exactly this multi-session situation).
- **Check job status first**: `qstat -u terencew` for jobs 14697667/14697668 (dipcall + genome-
  wide read-based phasing) — if they've cleared, item 1 above (the full re-call) is the
  immediate next action, and it will invalidate the specific numbers in §10-13 of this file
  until redone.
- **Today's numbers came with two real, silently-wrong bugs** (§10) caught only by directly
  spot-checking a from-scratch recompute against a known-good slice (chr1) rather than trusting
  an existing aggregate. Worth the same discipline going forward — any inherited intermediate
  table (not freshly regenerated from `tables/dasm_betabinom/*/regions_chr*.tsv`) should be
  treated as suspect until spot-checked, not assumed current.
- `doc/20260907.results.docx` and `md/20260907.results.summary.md` still have the wrong 91.9%
  proximal figure — do not cite them as-is; regenerate after the dipcall re-call, not before
  (regenerating twice would be wasted effort).
