> **Canonical copy:** `asm_lr_hprc2/reference/meredith_2026_comparison.md` (not a git repo; tracked here).

# asm_lr_hprc2 / ont_asm_caller vs. Meredith et al. 2026 (ASM-LR) — what was and wasn't scooped

**Written 2026-09-15** from a full read of `2024.12.16.628723v2.full.pdf` (56 pp., bioRxiv v2
posted 2026-08-27; text extracted beside it as `.full.txt`).

> Meredith, M.\*, Daida, K.\*, Moller, A.\*, … Paten, B., Blauwendraat, C., Billingsley, K. J.
> *Haplotype-Resolved Long-Read Sequencing in Hundreds of Diverse Brains Identifies Structural
> Variant Impacts on Expression and Allele-Specific Methylation.* bioRxiv 2024.12.16.628723 v2.
> NIH CARD. Not peer reviewed. Tool: https://github.com/molleraj/CARDlongread_allele_specific_QTL

`asm_lr` catalogued **v1** of this preprint on 2026-09-03 (`md/candidates.md`, `md/lit_review.md`),
but only as a *candidate dataset*. The ASM-LR analysis is what v2 adds.

**Naming collision, act on it regardless of scope:** their tool is called **ASM-LR**. This project
is `asm_lr`. Any manuscript, repo or talk using that name will read as theirs.

---

## 1. What ASM-LR actually is

Stated plainly because the abstract does not: **ASM-LR does not call allele-specific methylation in
any individual.** It is a *phased methylation QTL*, a population association with haplotypes as
the unit of observation.

1. **Cohort.** 351 post-mortem frontal cortex samples: NABEC, 205 EUR, **R9.4.1** / Guppy 6.12,
   46×; HBCC, 146 AFR/admixed after dropping 8, **R10.4.1** / Guppy 6.38, 39×. Read N50 about 25 kb.
2. **Phasing.** NAPU pipeline: Shasta + Hapdup *locally* phased assemblies, with Margin harmonising
   variant and methylation haplotypes. Phase-block NG50 is **1.25 Mb (NABEC) and 2.43 Mb (HBCC)**.
   No trio data. Haplotype labels are "arbitrary and assigned locally" (p. 19).
3. **Methylation.** `modkit pileup`, CpG motif, strands combined. The per-haplotype methylation
   fraction is **averaged over a region** (CGI, promoter or gene body) with ≥10 CpGs and ≥5 reads.
4. **Model.** Per variant–region pair, OLS:
   `Methylation_region(haplotype) ~ Allele_variant(haplotype) + age + sex + biobank + 20 PCs`.
   Each individual contributes two rows, one per haplotype. They describe this as "the number of
   observations is effectively doubled".
5. **Variants.** MAF ≥ 5%, HWE p > 1e-6, call rate ≥ 95%, **LD-pruned at r² < 0.3**, within ±500 kb.
6. **Scope.** **171 AD/PD GWAS genes** only; not genome-wide.
7. **Significance.** BH FDR < 0.05. No permutation or empirical null.

**Headline results.** 234 NABEC + 217 HBCC = **451 ASM-LR mQTLs** (Table 3; 228,087 and 1,745,252
tests). With CpGs inside deletions removed, this falls to **64** (54 + 10). SVs are the lead
variant in about 60% of hits and are about 10× as likely as SNVs to be a hit. **132 loci** were
"masked in unphased models". Extrapolating, they estimate about 1,500 extra mQTLs genome-wide.

### 1a. Verified against their code, 2026-09-15

`long_read_QTLs.py` in `molleraj/CARDlongread_allele_specific_QTL` (HEAD, read in full):

- **Model:** `smf.ols("Q(region) ~ variant + covariates", data=merged_data)`. Ordinary least
  squares; **no weights, no sample random effect, no sample fixed effect, no clustered errors.**
- **Rows:** `genetic_data.merge(methylation_data, on=['SAMPLE','HAPLOTYPE']).merge(metadata, on='SAMPLE')`.
  One row per haplotype. The per-person covariates (age, sex, biobank, 20 PCs) are copied onto both
  rows.
- **Haplotype matching is assumed, not checked:** "This is based on the assumption that for each
  sample, genetics H1 and methylation H1 match" (code comment). The regression does not ask whether
  the variant and the region sit in the same phase block.
- **Unphased mode** (`--force_dephased`) sums alleles to 0/1/2 dosage and averages the two haplotypes.
  That is ordinary mQTL, the tensorQTL comparison.
- **Significance:** `statsmodels.multitest.fdrcorrection` (BH) over all variant–region tests pooled.

### 1b. Why they call it "ASM" when no individual is called

It is a naming convention, not a mislabel, but it hides which comparison drives the result.

"Allele-specific methylation" has long meant two different things:

1. **Per-individual allelic imbalance.** One person's two alleles differ at a locus. This is what
   `ont_asm_caller` tests, and what "calling ASM" means in `asm_lr`.
2. **Allele-resolved methylation as a phenotype**, tested for association with genotype across people.
   This is the "haplotype-dependent ASM" / hap-ASM line of work, where the question is "which
   allele carries which methylation level?"

ASM-LR is sense 2. "Allele-specific" describes the **measurement** (methylation measured per
haplotype), not the **inference**, which is a population QTL. Their own results heading says it more
accurately: "haplotype-resolved mQTL".

**What the regression actually compares.** With haplotypes as rows and no per-person term, the slope
mixes two comparisons with no control over their relative weight:

- **Within heterozygotes:** H1 carries the variant and H2 doesn't, in the same person. This is a
  true allelic contrast, the per-individual ASM part. Age, post-mortem interval, cell composition,
  ancestry and trans effects all cancel inside a person.
- **Between homozygotes:** alt/alt people versus ref/ref people. This is ordinary mQTL signal. It
  does not cancel those confounders, which is why the model needs age, sex, biobank and 20 PCs.

Adding a **sample fixed effect** would leave only the within-heterozygote contrast: a pure cis
allelic test that controls confounding by design and needs no PCs. They have the data for that
design and did not fit it. As fitted, the "ASM" label rests on the heterozygotes, while an unknown
share of the power comes from between-person comparisons that a standard mQTL already makes. That
may be one reason unphased tensorQTL found more hits on shared tests (557 vs 132, §2).

## 2. Three things in the paper that cut against its own framing

Worth knowing before citing it, and for a reviewer-proof positioning statement.

- **Phasing did not dominate on shared tests.** Of the 1,014 region–variant pairs tested by both
  methods, **557 were significant only in unphased tensorQTL, and 132 only in ASM-LR** (p. 21). The
  authors attribute the 557 to effects under 1%, but the paper does not show that phasing adds
  power overall. It shows larger effect sizes plus some unique loci.
- **Most raw hits are CpG presence/absence.** Removing deleted CpGs cuts hits 451 → 64 (86%). Their
  flagship large effect (a 524 bp insertion in a chr12 CGI, 20–30% lower methylation) is explained
  in the text as CpGs *moving into the insertion*: a reference-coordinate artifact, not regulation.
  This is the artifact class that `asm_lr`'s `VALIDATION_PLAN.md` D1 (CpG-destroying variants) is
  built to partition.
- **Ancestry and chemistry are confounded.** EUR is R9 and AFR is R10; they say so (p. 26). Every
  cross-cohort comparison mixes the two.

Three methodological properties are **not addressed** in the methods. Each is structural, not a
guess about their results.

- **Within-individual pseudoreplication.** Both haplotypes of one person enter OLS as independent
  observations with identical covariates. Shared trans, cell-composition and post-mortem effects
  make them correlated, so the doubled *n* overstates the effective sample size. This is the same
  class of problem as the read×CpG design effect in `ont_asm_caller`, one level up.
- **No coverage weighting.** A 5-read haplotype average counts the same as a 40-read one. This is
  the depth-heteroscedasticity problem `asm_lr` diagnosed with Fisher's exact test.
- **Phase breaks inside the test window.** Local phase blocks (median 1.25–2.43 Mb) against ±500 kb
  windows and whole gene bodies mean the haplotype↔allele link can break inside a tested pair,
  biasing toward the null.

## 3. Side by side

| | **ASM-LR (Meredith 2026)** | **`ont_asm_caller`** | **`asm_lr_hprc2` plan** |
|---|---|---|---|
| Question | Which common variant shifts haplotype methylation, across people? | Do this person's two haplotypes differ here? | How much single-molecule methylation structure is genetic, and how far does it reach? |
| Unit | haplotype × individual, region average | read×CpG counts pooled per region, per individual (read-level path: molecules) | molecules, per haplotype, per individual |
| Needs a variant? | yes, MAF ≥ 5%, LD-pruned | **no**; the test is genotype-free | classification step uses variants |
| Overdispersion / coverage | none; OLS on averages, ≥5 reads | beta-binomial, depth in the likelihood, global ρ | inherits caller |
| Non-independence | haplotypes treated as independent | design effect measured: median **1.29** at 500 bp (HG00146 chr15) | same |
| Multiple testing | BH over variant–region pairs | BH over regions + Δ ≥ 0.10 | same |
| Tissue / n | brain, 351, R9 + R10 | LCL, 6 QC-pass of 18, R10.4.1 | LCL, 229 with ONT modbed, 202 with per-hap ASE |
| Ancestry | EUR (R9) vs AFR (R10) | 6 donors | 5 superpopulations; chemistry uniformity **unverified** |
| Genome-wide? | **no**; 171 genes | yes | yes |
| Imprinted / PofO ASM | **invisible by design** | detected as ASM; labelled by DMR overlap + cross-donor consistency | a taxonomy class |
| Non-genetic ASM | **invisible by design** | detected | a taxonomy class |
| Rare / private cis variants | excluded (MAF ≥ 5%) | detected per individual | detectable |
| Molecule structure (epialleles, entropy, co-methylation) | discarded (region average) | read-level + pattern paths | core question |
| Distance of genetic influence | not analysed | proximal/distal ~50/50 per donor (`asm_lr` §13) | core question |
| Expression link | unphased gene-level short-read eQTL | — | **haplotype-resolved ASM × ASE**, same individuals |
| SVs | central | none | not in plan |

"Invisible by design" is worth spelling out, because it is the whole positioning argument. In an
association across people, an imprinted DMR's methylation follows parental origin, which is
independent of which allele a person carries at a nearby variant, so there is no genotype
association to find. The same holds for stochastic or non-genetic ASM. A population regression can
recover **only** the genotype-driven classes, and only for common variants.

## 4. Verdict: what was scooped

**Taken:**
- Long-read haplotype-resolved methylation QTL at hundreds-of-genomes scale, with diverse ancestry.
- **SVs as drivers** of allele-specific methylation.
- "Phased reveals mQTLs that unphased misses", as a claim, with a named tool.
- The name.

**Not taken, and structurally unreachable from their design:**
1. **Per-individual ASM calling with calibrated FDR.** ASM-LR has no per-individual output.
2. **The four-class taxonomy's two non-genetic classes** (imprinted/parental, non-genetic), and
   therefore the headline question "what fraction of ASM is genetic?" Their design can only see
   the numerator.
3. **Single-molecule structure.** They average every molecule away before testing.
4. **How far genetic influence reaches**, and the long-read vs short-read phasing gain for ASM.
5. **Haplotype-resolved ASM × ASE in the same individuals.** Their RNA is unphased and gene-level.
6. **Genome-wide.** They tested 171 genes; the genome-wide figure is an extrapolation.
7. **LCL, five superpopulations, one chemistry.** Their ancestry contrast is chemistry-confounded;
   HPRC2's is not, *if* the ONT data are uniform. **Not verified here; check before claiming.**

**Bottom line.** If `asm_lr_hprc2` were executed as "phased mQTL across ancestries", it would be
substantially scooped. As written, its core question and taxonomy are the part of the problem
ASM-LR cannot reach. The paper helps the positioning: it is the citable population-association
half, and the HPRC2 project is the per-individual half that explains what association misses.

## 5. Corrections this comparison forced in `asm_lr_hprc2/md/README.md`

The README says the caller "carries forward … with its two confirmed calibration fixes
(region-pooling pseudoreplication correction, mean-dependent dispersion estimator)". **Neither is in
production.** Verified against the Hoffman2 production script on 2026-09-15:

- `scripts/wg/good_donors/G02_betabinom_region_chr1.py` calls the **global**
  `estimate_dispersion`, not `estimate_dispersion_trend`. The trend estimator exists in the package
  but was never wired in (`asm_lr/md/20260907_outstanding_items.md` item 6 also says so).
- The pooled region test has **no** pseudoreplication correction. The read-level path is a
  *separate* test. What exists is a measurement: DE is small at the 500 bp scale in use (median
  1.29), so pooling is broadly defensible, which is a different claim from "corrected".

Before the caller carries forward, also note that **HPRC2 modbeds are already binarised**. The
confidence filter and call-error estimate used on `modkit extract` (P08/P09; ε ≈ 0.043) cannot
be applied to them.

## 6. Recommended next steps (for the HPRC2 plan only)

1. **Rename** away from ASM-LR / `asm_lr` in anything outward-facing.
2. **Lead with the decomposition** ASM-LR cannot do: per individual, the share of ASM that is
   genotype-driven (proximal / distal), imprinted, or non-genetic, genome-wide, with the
   imprinting dose-response (0.1% at 1 donor → 53% at 6 donors in `asm_lr`) as the positive control.
3. **Promote ASM × ASE** (202 samples) from validation to a main result. It is the cleanest
   claim ASM-LR's data cannot make.
4. **Cite ASM-LR explicitly** as the population-association complement, and pre-empt the obvious
   reviewer question with the §2 facts (557 vs 132; 451 → 64).
5. **Verify HPRC2 ONT chemistry uniformity** before any "unconfounded ancestry" statement.
6. **Fix the README's calibration claim** (§5) before it propagates into a methods section.
