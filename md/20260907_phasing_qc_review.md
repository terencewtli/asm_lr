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
