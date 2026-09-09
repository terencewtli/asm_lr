# Literature review: where this project sits in the long-read ASM landscape

**Date:** 2026-09-08
**Purpose:** (1) map the current landscape of long-read allele-specific methylation (ASM);
(2) determine whether anyone has already done this project's specific re-analysis of HPRC ONT
data; (3) adjudicate, claim by claim, the novelty arguments this project rests on.

**Search method and its limits.** ~20 web searches plus direct retrieval of ~15 primary sources,
covering: long-read ASM methods and tools; population-scale long-read methylomes; HPRC/1KGP-ONT
methylation work; short-read ASM atlases; imprinting and parent-of-origin long-read studies;
methylation-haplotype-block / epiallele literature; XCI; and long-read-vs-short-read comparisons.
This is a web-index search, not a systematic PubMed/Scopus query with a registered protocol.
Nature-family full texts (`nature.com`) redirect to an auth wall, so several Nature/Nat Genet/Nat
Commun items below are characterized from PubMed/PMC mirrors, abstracts, or search snippets rather
than from the full text — those are flagged **[abstract-level]**. Anything not flagged was
retrieved directly. Preprint-only items are flagged **[preprint]**.

---

## 1. Bottom line

**Nobody has published the specific thing this project is building.** No genome-wide, unbiased,
cross-donor-replicated ASM catalog on HPRC ONT data exists in the literature or on the preprint
servers as of this search. That much of the concern is unfounded.

**But the reasons the project gave for why it is novel are not all equally strong, and two of the
five are substantially weaker than assumed.** In descending order of how well each holds up:

| # | Novelty claim | Verdict |
|---|---|---|
| 4 | Co-methylation / epiallele structure analyzed *per genetic haplotype* | **Holds.** Closest to unoccupied ground. |
| 5 | No concrete quantification of long-read gain over short-read for ASM | **Holds** — but the project's own estimate is currently unstable (see §7). |
| 2 | No standardized, statistically careful single-sample ASM caller | **Partly holds**, and shrinking fast. Two 2026 methods now occupy adjacent ground. |
| 3 | Field emphasizes imprinting over unbiased ASM discovery | **Largely does not hold.** At least four unbiased genome-wide ASM studies exist. |
| 1 | Nobody has done it on HPRC | **True but weak.** Dataset novelty is the least persuasive kind, and HPRC's own methylation working group is the obvious party to scoop it. |

The correct framing is therefore **not** "first unbiased ASM discovery" and **not** "first on
HPRC." It is: *a calibrated, single-sample ASM caller with an explicit treatment of ONT read
pseudoreplication, used to measure long-range haplotype–epiallele structure at a depth the
population-scale cohorts cannot reach.* Depth and method, not scale or dataset.

The most serious problem found in this review is internal, not external: the project currently
carries **three mutually inconsistent estimates of its own headline number** (§7). That is a
larger threat to publication than anything in the literature.

---

## 2. The two papers that prompted this review

### 2.1 Sigurpalsdottir et al. 2026, *Nature Communications* — the "aging clock" paper

> "Nanopore sequencing identifies parent-of-origin specific age-associated methylation changes at
> imprinted loci in the human genome." doi:10.1038/s41467-026-76240-w
> <https://pmc.ncbi.nlm.nih.gov/articles/PMC13522577/>

| | |
|---|---|
| Cohort | 7,284 Icelandic whole-blood samples (deCODE), ages 1–106; validation n=545 (HerediGene/Intermountain) |
| Platform | ONT PromethION, **R9.4.1** |
| Coverage | **mean 15.9×**, N50 ~18.9 kb |
| Methylation calling | Nanopolish v0.13.3 (log-likelihood ratios) |
| Phasing | Parent-of-origin via long-range phasing against 63,460 Icelandic WGS; 8.96M high-quality variants; Graphtyper |
| Statistics | Linear regression of methylation on age, sex, 10 PCs, **and cis ASM-QTL genotype**; Bonferroni p < 2.8×10⁻⁹ |
| Findings | 20.8% of ~18M CpG units age-associated (92.4% hypomethylating); clock MAE 2.43 yr; **702 CpG units with parent-of-origin-specific age association**, 94.5% within 100 kb of imprinted loci; DIRAS3 paternal hypermethylation with age |
| Data access | Raw data restricted by Icelandic law; summary statistics on Zenodo |

**Why this is less overlapping than it looks.** It is an *epidemiological* study that uses ASM as a
covariate and a stratification, not an ASM-discovery study. Three concrete separations:

- **Coverage.** 15.9× total is ~8× per haplotype. That is at or below the floor for per-locus
  read-level epiallele analysis. This project's donors sit at 51–87× (5–10× more per haplotype).
  Anything requiring the *within-molecule CpG vector* — epialleles, entropy, co-methylation
  blocks, model-order selection — is out of reach at their depth. This is the single cleanest
  differentiation available and should be stated numerically in any paper.
- **Chemistry.** R9.4.1 with Nanopolish, vs. this project's R10.4.1/Dorado with joint
  5mCG_5hmCG. R10 is documented as materially better for methylation calling. Notably, this is
  the same chemistry split the project already made internally (`md/basecaller_cohort_split.md`)
  — the argument for excluding the 12 R9 donors is *the same argument* that separates this work
  from deCODE's, which is a rhetorically useful symmetry.
- **Access.** Their data cannot be re-analyzed by anyone. HPRC's can. That is a genuine, if
  unglamorous, contribution argument.

**Where it genuinely competes.** They *have* a cis-ASM-QTL map at n=7,284 that they use as a
regression covariate. No population-scale ASM claim from n=18 will survive comparison. Do not
frame anything in this project as a population-genetic result.

### 2.2 Akbari et al. 2022, *eLife* 77898 — the PofO imprinting paper

> "Genome-wide detection of imprinted differentially methylated regions using nanopore sequencing."
> <https://elifesciences.org/articles/77898>

| | |
|---|---|
| Cohort | **12 LCLs** from 1KGP + GIAB, all with trio data (incl. NA12878, NA19240); draws on Human Pangenomics/GIAB public nanopore data (HG002, HG005) |
| Method | NanoMethPhase (their own tool, Genome Biol 2021) + Nanopolish; trio-SNV phasing |
| Statistics | **DSS** (dispersion-shrinkage beta-binomial), p<0.001, \|Δ\|>0.20, present in ≥4 cell lines |
| Scale | Phased 26.5M autosomal CpGs = 95% of autosomal methylome |
| Findings | 143 iDMRs: 101 known + **42 novel** (16 germline, 26 somatic); recovered 94% of well-characterized iDMRs |

**The user's read of this paper is correct.** The statistical treatment is DSS applied off the
shelf. There is no per-donor dispersion estimation, no design-effect or pseudoreplication
correction, no coverage-power model, and the cross-sample criterion (≥4 of 12 cell lines) is a
frequency filter, not a calibrated FDR. The replication rule this project uses (≥2 of 6) is the
same species of arbitrary threshold — which is worth noticing, because the criticism cuts both
ways until the sensitivity analysis in `md/20260907_outstanding_items.md` #11 is actually run.

**Scope separation is clean:** they target imprinting (PofO), require trios, and use n=12 LCLs at
unstated (but 2021-era, therefore modest) coverage. This project targets genotype-driven ASM
without trios, at 18 donors, with proximal/distal architecture classification. Different question.

---

## 3. The landscape, in six clusters

### Cluster A — Population-scale long-read haplotype-resolved methylation (the real competition)

| Study | Data | What it does | Distance from this project |
|---|---|---|---|
| **Gustafson et al. 2024**, *Genome Research* 34:2061 — 1KGP-ONT ([link](https://genome.cshlp.org/content/34/11/2061)) | ONT, **37× mean, N50 54 kb**; ~100 samples reported, consortium target **800+** | Variant catalog. Methylation is a *demonstration*: imprinted-locus patterns, XCI skew, some novel DMRs; SNURF-SNRPN haplotype imbalance in 2 samples | **The closest published analog, and it is not close.** Methylation is a figure, not the paper. But the consortium owns the natural follow-up. Repo item #13 — still worth pulling in full. |
| **1KGP-ONT via methylmap** (BMC Bioinformatics 2025) | **226 individuals**, PromethION | Visualization of haplotype-specific methylation across the cohort | Shows the haplotype-resolved 1KGP-ONT methylome is already *public and browsable* at 12× this project's n |
| **Chinese population ONT methylation atlas** [preprint] [abstract-level] ([bioRxiv Apr 2026](https://www.biorxiv.org/content/10.64898/2026.04.20.719515v1.full)) | **106 individuals**, 19 provinces, ONT | Haplotype-resolved atlas; 27.6M CpGs; defines three DMR classes: segment-related (sDMR), **haplotype-based (hDMR)**, population-specific (pDMR); SVs as pervasive methylation covariate; altitude as environmental determinant | **The single most design-similar paper found.** Same architecture: population ONT cohort → haplotype-resolved methylation → DMR taxonomy. Differentiators are ancestry scope (one population), likely coverage, and no proximal/distal SNP-distance architecture. Read this in full before writing. |
| **METAFORA** — Jensen et al. 2026, *medRxiv* [preprint] ([link](https://pmc.ncbi.nlm.nih.gov/articles/PMC13278280/)) | **551 long-read genomes** (UDN/GREGoR), ONT R9/R10 + PacBio Revio, median ~25–35× | Population-scale methylation *outlier* detection: depth-aware beta test at each CpG, change-point segmentation of >36M CpGs into **correlated blocks**, hidden-factor PCA, M-value residualization, winsorized depth-dependent thresholds. Haplotype-specific profiles intersected with outliers; 32,808 outlier regions; **16.4% show clear allelic effects** | **The most direct methods competitor.** It already does depth-awareness, co-methylation block structure, and covariate adjustment at population scale. Different *target* (rare outliers vs. common ASM) but overlapping statistical territory. This one materially weakens novelty claim #2. |
| **Groza et al. 2025**, *Genome Research* ([link](https://pmc.ncbi.nlm.nih.gov/articles/PMC12047246/)) | **435 PacBio HiFi methylomes** (GA4K), 26×, 470 phased assemblies | Pangenome-graph methylome: 14.6M non-reference CpGs (+43.6%); 230,464 SV-associated methylation bins; SVs 17.2× more likely to be mQTLs than SNPs; mean SV-QTL interaction distance **37.8 kb** vs. 19.8 kb for SNP-QTLs | Important prior art for "distance matters." Their 37.8 kb figure is the closest existing analog to this project's distal-ASM argument — cite it, and note it is SV-driven and HiFi, not ONT/SNP |

### Cluster B — HPRC-specific methylation work

| Study | HPRC usage | Verdict |
|---|---|---|
| **Zhuo et al. 2026**, *Genome Research* ([link](https://pmc.ncbi.nlm.nih.gov/articles/PMC13262945/)) | **32 HPRC Year-1 samples** (HiFi methylation) + 5 with ONT R9.4.1 + WGBS | Polymorphic transposable-element insertion methylation, genome-wide, with **haplotype-resolved comparison** on phased HiFi/ONT reads. Validates long-read vs. WGBS concordance: **r = 0.901 (HiFi), 0.925 (ONT)**. Finds TEs mostly methylated; ~3% flanking-methylation spreading within 300 bp. **This is HPRC ASM-adjacent work, scoped to TE insertions.** Useful as a concordance citation; not a scooping risk |
| **HPRC Release 2 Epigenomes** ([AWS Open Data](https://registry.opendata.aws/hprc-epigenome/); [HPRC](https://humanpangenome.org/hprc-data-release-2/)) | 232 individuals; ONT UL for 200+; harmonized HiFi methylation subset curated by an **HPRC methylation working group**; interactive HPRC Epigenome Browser (WashU) | **This is the main scooping risk in the entire review.** A consortium working group with harmonized methylation across all 232 R2 samples is the natural author of "the HPRC ASM atlas." No such paper found yet — but a consortium resource paper would supersede an n=18 external re-analysis on scale alone |
| **HPRC2 pangenome paper** [abstract-level] ([bioRxiv/PubMed 42539208](https://pubmed.ncbi.nlm.nih.gov/42539208/)) | 460 haplotypes, >99% of common variation in All of Us v8 | Assembly/variation paper; no ASM component found |

### Cluster C — Short-read ASM atlases (the "where" this project answers "why" for)

- **Rosenski et al. 2025**, *Nat Commun* 16:2141 ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11897249/)) — 202 WGBS samples, ~40 purified cell types, 135 donors, 150 bp PE, ~984M read pairs/sample. **324,759 bimodal regions; 34,426 SNP-linked (~10%); 460 parental-ASM regions** (45 known iDMRs, 78 near imprinted genes). EM bimodal calling + LLR. Explicitly states the read-length limitation: with ~200 bp reads "it is impossible to assess the presence of distant genetic variants controlling methylation," and therefore that they may underestimate sequence-dependent ASM. **This is the project's strongest citation: the field's best short-read atlas naming exactly the gap this project fills, in the authors' own words.**
- **Abante, Fang, Feinberg & Goutsias 2020 (CPEL)**, *Nat Commun* ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7567826/)) — haplotype-dependent ASM in WGBS via a 1-D Ising model over SNP haplotypes; information-theoretic comparison of methylation *stochasticity* between alleles. **Their stated motivation is verbatim the project's methods argument**: existing approaches "perform statistically independent marginal analysis at individual CpG sites, thus ignoring correlations in the methylation state," or do unreliable joint analysis over four CpGs. Repo item #8 flags CPEL as unevaluated — it is more than an off-the-shelf alternative, it is **prior art for the core statistical claim**, and must be cited and positioned against.
- **Zink et al. 2018**, *Nat Genet* 50:1542 [abstract-level] ([link](https://www.nature.com/articles/s41588-018-0232-7)) — deCODE, 285 methylomes + 11,617 transcriptomes with PofO-phased haplotypes. Imprinted methylation is *continuous*, not binary; polymorphic imprinted methylation at VTRNA2-1, PARD6G, CHRNE. Already used in this project as the 229-DMR validation set.

### Cluster D — ASM tools and callers

| Tool | Year / venue | Approach | Relevance |
|---|---|---|---|
| **NanoMethPhase** | 2021, *Genome Biol* ([link](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02283-5)) | Megabase-scale methylation phasing; DMA module; works at low coverage | The de facto reference implementation. Any new caller is benchmarked against it |
| **nanoASM** — Tian et al. 2026 [preprint] ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13320949/)) | ONT R10.4.1; 37 normal + 15 tumor prostate | Partitions reads by allele, **de novo binary segmentation with metilene**; **406,011 tumor / 391,918 normal ASM events**; methylation entropy; nominates rs6885084 (IRX4), rs4736369 (PSCA); validates vs. TCGA mQTL / GTEx eQTL | **Unbiased genome-wide ASM discovery already exists.** This is the paper that most damages novelty claim #3. Note it does *not* compare long-read vs. short-read detection — leaving claim #5 intact |
| **NANOME** | 2025, *Brief Bioinform*/[PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12236756/) | Nextflow; XGBoost consensus of Megalodon/Nanopolish/DeepSignal; +11% single-base precision, +2.4% single-molecule accuracy, ~200k extra CpGs; haplotype-aware ASM at imprinting controls | Caller-level, not statistics-level. Positions R9-era consensus calling; less relevant post-R10/Dorado |
| **cyberDMR** — Li et al. 2026, *BMC Biology* 24:165 ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13371177/)) | **Weighted beta regression + LRT**; F-like pre-filter for excess intra-group variability; coverage-adaptive smoothing (recovers 38.9% of <5× sites); seed-growing clustering for "local methylation coherence." Benchmarked vs. DSS, metilene, BSmooth, HOME, MethyLasso, DiffMethylTools. **Reports 92 hDMRs from phased GIAB trios** | **Second-most-damaging to claim #2.** A 2026 beta-regression DMR caller that explicitly handles overdispersion, coverage variation, and CpG coherence, and demonstrates haplotype DMRs. The project's beta-binomial caller must be benchmarked against this or the methods claim will not survive review |
| **ASMS** — 2026 ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13326645/), [GitHub](https://github.com/ecmra/asms)) | Finds ASM **without phasing**, from MM/ML tags; scans thousands of loci fast | Orthogonal approach worth a comparison paragraph; a reviewer may ask why phasing is needed at all |
| **MethPhaser** | 2024, *Nat Commun* [abstract-level] ([link](https://www.nature.com/articles/s41467-024-49588-0)) | Uses methylation to *extend* haplotype phasing | Directly relevant to the project's open haplotagging-rate problem (~38% HP-tagged on NA18508/chr20 vs. a ≥60% gate) |
| **Methylation-aware long-read phasing** [preprint] [abstract-level] ([bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.03.11.710820v2.full)) | Reports phase-block N50 increases of **78–151%** at 83.4–98.7% phasing accuracy on ONT R9/R10 | Same — a concrete, citable remedy for the phasing-rate gap |
| Visualization: **methylartist**, **modbamtools**, **modkit**, **methylmap** | 2022–2025 | Plotting/aggregation | Ecosystem context |

### Cluster E — Methylation haplotype blocks, epialleles, single-molecule structure

- **Guo et al. 2017**, *Nat Genet* [abstract-level] ([PubMed](https://pubmed.ncbi.nlm.nih.gov/28263317/)) — methylation haplotype blocks (MHBs) defined by LD (r²) of epialleles; used for tissue deconvolution and plasma tissue-of-origin. The founding paper of the co-methylation-block idea. **Short-read, bulk, not haplotype-stratified.**
- **MHB landscape in human tissues and preimplantation embryos**, 2023 ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10760529/)) — regulatory elements defined by co-methylation patterns.
- **Genome-wide single-molecule analysis of long-read DNA methylation**, *PLoS Genetics* 2023 ([link](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010958)) — heterogeneous single-molecule patterns at heterochromatin reflecting **nucleosome organization**; nucleosome-scale periodicity in methylation correlating with accessibility. The best existing treatment of within-molecule structure — and it is *not* haplotype-stratified.
- **Tian et al. 2025**, *HGG Advances* 7:100532 ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12642126/)) — nanopore adaptive sampling on 22Rv1 + public cancer/normal cell lines. Identifies methylation haplotype blocks via LD r², does de novo DMR detection between haplotypes, finds **6,390 haplotype-specific DMRs in 22Rv1**, shows cancer genomes have **smaller methylation blocks** than normal, and states **nanopore captures co-methylation better than bisulfite**. Also improves CTCF binding prediction by combining methylation + motif. **The nearest neighbor to novelty claim #4 — but in cell lines, targeted (adaptive sampling), n≈5, and it does not stratify block structure by genetic haplotype across donors.**
- Metric vocabulary already established and to be adopted rather than reinvented: methylation frequency, **methylation entropy**, **epipolymorphism**, **methylation haplotype load**, per-read Shannon entropy, read transition statistics.
- Note against the project's framing: at least one benchmark reports that **high correlation between neighboring CpG sites is the exception, not the rule.** If co-methylation is weaker than assumed genome-wide, both the pseudoreplication concern *and* the epiallele-structure payoff are region-dependent, not global. Measure it, don't assume it.

### Cluster F — Adjacent applications (context, not competition)

- **ASD parent-of-origin methylomes**, *Sci Adv* 2026 [abstract-level] ([link](https://www.science.org/doi/10.1126/sciadv.aee4069)) — PacBio HiFi, **124 individuals (31 quartets)**, phased methylomes; 114 paternal / 106 maternal DMCs, 45/46 DMRs, 2,425/2,693 methylation outliers.
- **APOE haplotype-resolved methylation**, *npj Dementia* 2026 [abstract-level] ([link](https://www.nature.com/articles/s44400-026-00094-8)) — 332 postmortem brains, two ancestries, 18 novel APOE-associated DMCs. Locus-scoped.
- **Voronina et al. 2025**, *IJMS* ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12525282/)) — **single genome**, PBMC, ONT at **112×**, trio-binned + computational phasing, Megalodon + NanoMethPhase, coverage titration 12×–112×. **308,673 ASM sites (0.5% of 58.5M positions)**; SNV density 4–10× elevated near unmethylated ASM sites; ASM enriched at bivalent promoters/enhancers; **45% of ASM sites overlap known mQTLs**. Closest single-sample deep-ONT analog; note the 0.5% rate is strikingly close to this project's 0.39–0.54% and deCODE's 0.51% benchmark — a useful convergent-validity citation.
- **XCI:** adaptive-sampling XCI skew (*Genome Research* 2024, [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610576/)); haplotype-resolved XCI in AIFM1/PDHA1 carriers (2026, [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13347936/)); lineage-specific XCI escape and skew [preprint] ([bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.08.26.739472v1.full)). Escapees show reduced promoter-CGI methylation — the orthogonal validation the project's chrX arm would rely on. **The XCI niche is more crowded than the ASM niche; treat chrX as a positive control, not a chapter.**
- **Non-Mendelian inheritance of ASM**, *Nat Genet* 2026 [abstract-level] ([link](https://www.nature.com/articles/s41588-026-02603-0)) — genome-wide long-read ASM framework, but **in mice**. ~93% of autosomal epigenetic inheritance Mendelian via cis-meQTL; 522 non-Mendelian cases; emergent epialleles; paramutation at Capn11. Not a human competitor, but it establishes "genome-wide long-read ASM framework" as a Nature Genetics-caliber framing.

---

## 4. Claim-by-claim adjudication

### Claim 1 — "Nobody's done it on HPRC data." **True, but the weakest argument available.**

Confirmed: no genome-wide unbiased ASM catalog on HPRC ONT exists. The only HPRC methylation
paper found (Zhuo et al. 2026) is TE-scoped and mostly HiFi.

Three corrections to the reasoning behind the claim:

1. **"HPRC is higher coverage" needs checking.** HPRC R2's stated ONT target is ~30× ultra-long
   (>100 kb); 1KGP-ONT (Gustafson) reports **37× mean with N50 54 kb**. On the published targets,
   1KGP-ONT is the *deeper* resource. This project's own donors measure 51–87× — but that is a
   property of *this multi-run merge*, not of HPRC generically. The defensible sentence is
   "our merged per-donor coverage of 51–87× exceeds published population ONT cohorts
   (deCODE 15.9×, 1KGP-ONT 37×)", not "HPRC is higher-coverage."
2. **"Diverse lines" does not differentiate from 1KGP-ONT.** HPRC's samples are drawn largely
   from the same 1KGP LCL panel, and 1KGP-ONT is at 226+ individuals and heading to 800. On
   both n and ancestry breadth, 1KGP-ONT dominates. (The project already knows AFR is thin at
   n=2 — `md/20260904_cohort_provenance.md`.)
3. **LCL provenance is an unaddressed reviewer objection.** *[Flagged for verification, not
   confirmed in this search:]* HPRC/1KGP material is LCL-derived. EBV transformation, culture
   drift, and clonal expansion all produce monoallelic-looking methylation and inflated epiallele
   structure that is not germline-driven. Rosenski used **purified primary cell types**; deCODE
   used **whole blood**. Both will be held up against an LCL-based ASM catalog. This needs an
   explicit paragraph and ideally a positive control (e.g. do imprinted DMRs behave normally?
   they do — 48.9% Zink recovery — which is reassuring but does not address clonality at
   non-imprinted loci).

**Recommendation:** demote the HPRC argument to a data-availability/reproducibility point
("uniquely re-analyzable at this depth"), and stop leading with it.

### Claim 2 — "No standardized, statistically careful single-sample ASM calling." **Partly holds; eroding.**

The user's reading of the field is substantially right. Concrete evidence *for* the claim found
in this search:

- Akbari et al. use off-the-shelf DSS with hard thresholds and a ≥4-of-12 frequency filter.
- A nanopore methylation benchmark states plainly that **"in the absence of statistically robust
  DMR-finding methods for nanopore data"** the authors used permissive thresholds and accepted
  that many lower-ranked DMRs are likely false positives.
- nanoASM uses metilene binary segmentation — a change-point method with no dispersion model.
- No paper found reports a **design-effect / pseudoreplication measurement for ONT read-level
  methylation.** This project has measured one on real data (`P05_design_effect_real_data.py`).
  **That specific empirical result appears genuinely unpublished and is the most defensible
  methods contribution in the project.**

But three pieces of prior art constrain how the claim can be phrased:

- **CPEL (Abante et al. 2020)** already made the "marginal per-CpG analysis ignores methylation
  correlation" argument, and solved it with an Ising model — for WGBS. The project cannot claim
  the *insight*; it can claim the *ONT-specific instantiation and empirical quantification*.
- **METAFORA (2026)** already does depth-aware beta testing, correlated-block segmentation, and
  hidden-factor covariate adjustment at n=551 with haplotype-specific output.
- **cyberDMR (2026)** already does beta regression + LRT with coverage-adaptive smoothing, an
  explicit intra-group-variability filter, and demonstrated hDMRs on GIAB trios.

**Revised, survivable phrasing:** *"a haplotype-contrast ASM caller for ONT that estimates
per-donor dispersion and corrects for the read-level design effect induced by multi-CpG reads —
neither of which is handled by existing region-level DMR methods (DSS, metilene, cyberDMR) or
population-outlier frameworks (METAFORA)."* And then **benchmark against cyberDMR and DSS**,
because a reviewer will ask.

### Claim 3 — "Emphasis is on imprinting, not unbiased discovery." **Largely does not hold.**

Unbiased genome-wide ASM discovery has been done, repeatedly:

- nanoASM 2026: 406,011 / 391,918 ASM events, de novo, genome-wide, ONT R10.4.1.
- Voronina et al. 2025: 308,673 ASM sites, single genome, ONT 112×, unbiased.
- Rosenski et al. 2025: 324,759 bimodal regions, 34,426 SNP-linked, unbiased (short-read).
- Chinese ONT atlas 2026: hDMRs called genome-wide across 106 individuals.
- METAFORA 2026: 32,808 outlier regions, 16.4% with allelic effects.

The imprinting-heavy impression is real for the *long-read* subfield specifically (Akbari, Zink,
Sigurpalsdottir, the ASD paper) — but "unbiased ASM discovery" is not by itself an open niche.

**What is still open, and should be the actual claim:** unbiased discovery *plus* (a) explicit
proximal/distal/very-distal architecture classification by SNP distance, (b) cross-donor
replication with a calibrated rather than arbitrary threshold, and (c) per-locus read-level
confirmation. No paper found does all three.

### Claim 4 — "Methylation haplotype structure per genetic haplotype." **Holds. Strongest claim.**

Search found no study that computes co-methylation block structure, epiallele entropy, or
epipolymorphism **separately within each genetic haplotype, across a multi-donor cohort.** The
nearest work and why each falls short:

- Guo 2017 MHBs and the 2023 MHB landscape: bulk, short-read, no haplotype stratification.
- PLoS Genetics 2023 single-molecule heterogeneity: long-read, genuinely single-molecule, but not
  haplotype-stratified — it partitions by chromatin context, not by allele.
- Tian 2025 (HGG Adv): does LD-based MHBs *and* haplotype DMRs on nanopore — closest — but n≈5,
  cell lines, adaptive sampling, and blocks and haplotypes are analyzed as separate axes rather
  than blocks *conditioned on* haplotype.
- nanoASM: computes methylation entropy, but by sample/tumor status, not per haplotype.

The scientific question — *does the two-haplotype contrast differ in co-methylation architecture,
not just in mean, and is that difference genotype-linked?* — is well-posed, requires exactly the
coverage this project has, and is unclaimed. The existing epiallele scaffolding
(`md/epiallele.results.md`, `W12_epiallele_entropy.py`) is the right asset; it is currently
chr19/HG01258 only, on the *excluded* R9 pilot sample, and has no `scripts/wg/` equivalent. **This
is the highest-value unfinished work in the repo.**

Caveat to respect: the benchmark finding that strong neighboring-CpG correlation is the exception
means the effect may be confined to specific region classes. Frame as "where co-methylation
exists, is it haplotype-structured?" rather than assuming it is everywhere.

### Claim 5 — "No satisfactory quantification of long-read gain over short-read." **Holds — and it is the most valuable gap found.**

Searched specifically for a head-to-head quantification and found none. What exists instead is a
uniform pattern of **qualitative assertion**: "long reads improve ASM analyses," "allele-biased
patterns are more clearly detected," "obtaining genome-wide haplotyped methylomes with short reads
remains challenging." Rosenski et al. state the limitation precisely and cannot measure it
(10% of bimodal regions SNP-linked; the rest "may feature ASM with distant sequence determinants
that cannot be captured with short-read sequencing"). nanoASM does not compare to short reads.
Zhuo et al. compare *concordance* (r = 0.90–0.93) but not *ASM detection yield*.

The nearest quantitative anchors that do exist, and which any such claim must engage:

- Groza et al.: SV-QTL mean interaction distance **37.8 kb** vs. SNP-QTL **19.8 kb**.
- Rosenski: **~10%** of bimodal regions SNP-linkable at 200 bp.
- A co-methylation report: most co-methylation is **≤200 bp**, with only **1–3%** at ≥1,000 bp —
  a number that, if it generalizes, *caps* how much a purely co-methylation-based long-read
  advantage can be.

**This is publishable on its own** — but see §7 first, because the project does not currently
have a stable version of this number.

---

## 5. Scooping risk, ranked

1. **HPRC methylation working group / HPRC Epigenome resource paper.** Harmonized methylation
   across 232 R2 samples already exists as a data product. A consortium ASM/methylation resource
   paper would supersede an n=18 external re-analysis on scale, and consortium papers are not
   telegraphed in advance. *Mitigation: do not compete on scale. Compete on the caller, the
   design-effect result, and per-haplotype epiallele structure — none of which a resource paper
   would do well.*
2. **1KGP-ONT consortium methylation paper.** 226+ haplotype-resolved methylomes are already
   public and browsable; the consortium is heading to 800. The obvious next paper is exactly
   "haplotype-resolved methylation across 1KGP-ONT." *Mitigation: same. Also consider framing
   this cohort as the deep-coverage complement to their breadth.*
3. **The Chinese ONT atlas group extending to ASM architecture.** They already have the cohort,
   the haplotype-resolved pipeline, and an hDMR class defined. Adding SNP-distance classification
   is a small increment for them.
4. **nanoASM authors generalizing off prostate.** They have a working unbiased genome-wide ASM
   pipeline and would need only a public cohort.
5. **METAFORA/cyberDMR authors turning their statistics on common ASM.** Both have the
   statistical machinery; neither has aimed it at haplotype-contrast ASM discovery yet.

---

## 6. What a reviewer will attack

- **"Why n=18 when 1KGP-ONT has 226 public and deCODE has 7,284?"** Needs a one-sentence answer
  in the abstract: depth per haplotype, R10.4.1 chemistry, and single-molecule analyses that
  cannot be done at 16–37×. Lead with the number.
- **"Your data are LCLs."** Needs an explicit paragraph and, ideally, a clonality control.
- **"Your replication threshold (≥2 of 6, 250 bp slop) is arbitrary."** It is — the repo says so
  (`20260907_outstanding_items.md` #11). Akbari's ≥4-of-12 is equally arbitrary, but that defence
  works only once the sensitivity analysis exists.
- **"Six donors, not eighteen."** Two zero-hit donors are unexplained (HG00344, NA21144 — extreme
  ρ̂ 0.097/0.116), five samples are still incomplete, and one (NA18959) had no confirmed raw-data
  directory. A methods paper that drops a third of its cohort for unexplained reasons will be
  asked why. **This is arguably a bigger obstacle than any competing paper.**
- **"You claim a statistical contribution; how does it compare to cyberDMR / DSS / CPEL /
  METAFORA?"** Currently unanswerable. Needs a benchmark table.
- **"38% haplotagging rate."** Below the project's own ≥60% gate. MethPhaser and the 2026
  methylation-aware phasing preprint are the citable remedies.

---

## 7. The internal problem, which is larger than the external one

The project currently contains **three incompatible answers to its own headline question**:

| Source | Cohort | Proximal | Distal | Very distal | Non-genetic |
|---|---|---|---|---|---|
| `progress.md` (W07, phase 1) | HG01258, N=1, ~2.9×/hap, R9 | 6.1% | 10.6% | 1.5% | **81.8%** |
| `20260906.progress.md` (chr1, beta-binomial) | 6 donors, R10 | **50.5%** | **46.0%** | 3.5% | n/a |
| `20260907.results.summary.md` (finalized ≥2/6) | 6 donors, genome-wide | **91.3%** | **8.2%** | 0.4% | n/a |

The long-read-advantage claim reads as "3× more genotype-driven ASM" in the first, "distal is
nearly half of all ASM" in the second, and "distal is 8%" in the third. Those are three different
papers. The 2026-09-07 note that the finalized breakdown "matches the pre-finalization ASM01
breakdown (91.9/7.4/0.7) closely" is true *within* that analysis but does not reconcile it with
the 50.5/46.0/3.5 from the previous day, which the same doc attributes to a corrected bug.

Compounding this: the multiplier has **never been measured against actual short-read data.** W11
(bismark + short-read phasing baseline) is proposed in `progress.md` and never run. The current
number is a ≤200 bp *proxy* for short-read detectability. Given that claim #5 is the most valuable
open gap in the literature, shipping a proxy-based estimate into a field that has never measured
this would invite exactly the criticism the paper is trying to level at everyone else.

**Before any writing starts, one number has to be settled and defended.** Recommended order:

1. Reconcile the three breakdowns and document which is authoritative and why.
2. Run W11 (or an equivalent: down-sample the ONT reads to 200 bp-equivalent phasing windows and
   re-call, which is cheaper than a bismark pipeline and arguably a cleaner controlled comparison
   since it holds the caller, cohort, and coverage fixed and varies only the phasing distance).
   The in-silico version is the stronger experiment — same molecules, one variable.
3. Only then quantify the multiplier, with a confidence interval across donors.

---

## 8. Recommended repositioning

**Do not write:** "the first unbiased genome-wide ASM catalog," or "the first ASM analysis of
HPRC." Both are contestable and neither is the strength.

**Write instead**, in this order:

1. **Methods contribution.** A calibrated single-sample haplotype-contrast ASM caller for ONT
   with per-donor dispersion estimation and an explicit read-level design-effect correction —
   with the *measurement* of that design effect on real data as a standalone result, since no
   published equivalent was found. Benchmark vs. DSS, cyberDMR, and NanoMethPhase's DMA.
2. **The controlled long-read-gain experiment.** Same reads, same caller, phasing window varied
   from 200 bp to 50 kb. This directly answers the question Rosenski et al. pose and cannot
   answer. Expect the honest answer to be smaller than 3×.
3. **Per-haplotype epiallele architecture.** Co-methylation block structure, entropy, and
   epipolymorphism computed within each genetic haplotype — the least-occupied ground found in
   this review, and the analysis that most requires 51–87× coverage.
4. **QC/reporting heterogeneity across HPRC lines** as a resource contribution: the truncation
   bug, the dispersion anomalies (ρ̂ 0.097/0.116), the coverage manifest, the haplotagging-rate
   audit. This is real and useful and nobody publishes it — but it is a *discussion section or a
   supplementary resource*, not a chapter. It does not carry a paper.
5. **Validation** (already in hand and genuinely good): 48.9% Zink iDMR recovery at a ~30%
   background; 2.6× monocyte-meQTL enrichment; ASM rate 0.39–0.54% converging with deCODE's 0.51%
   and Voronina's 0.5%.

**Venue implication:** this is a *Genome Research* / *Genome Biology* methods-plus-resource paper,
not a *Nature Genetics* discovery paper. The comparison set (Gustafson, Groza, Zhuo, Tian all in
Genome Research; cyberDMR in BMC Biology; nanoASM heading to a methods venue) supports that
placement. Framed as a methods paper, n=18 is not a weakness — deCODE's 7,284 genomes are
irrelevant to whether the caller is calibrated.

**Overall verdict on whether the idea survives:** yes, but not in its current framing, and not
before the three-way inconsistency in §7 is resolved and the cohort is either restored to 18 or
the reduction to 6 is justified. The literature is denser than assumed, and two of the five
novelty pillars do not bear weight — but the two that do (per-haplotype epiallele structure; the
first real measurement of long-read ASM gain) are each strong enough to anchor a paper, and the
methods work is genuinely ahead of what the field currently uses.

---

## 9. Gaps in this review

- Nature-family full texts were not accessible; Zink 2018, MethPhaser, the *Nat Genet* 2026
  non-Mendelian paper, the *Sci Adv* ASD paper, and the *npj Dementia* APOE paper are
  characterized at abstract level only.
- The Chinese ONT atlas (the most design-similar paper) was only reachable at abstract level due
  to repeated bioRxiv rate-limiting. **It should be retrieved and read in full before writing.**
- Gustafson et al. 2024 was still not read in full (repo item #13 remains open).
- No systematic PubMed/Scopus query with recorded search strings; no conference-abstract search
  (ASHG/AGBT 2025–2026), where a competing HPRC ASM analysis would most likely surface first.
- LCL provenance of HPRC ONT material was not verified against a primary source in this session.

## 10. References

Retrieved in full unless marked. `[A]` = abstract-level only, `[P]` = preprint.

1. Sigurpalsdottir et al. 2026. Nanopore sequencing identifies parent-of-origin specific age-associated methylation changes at imprinted loci. *Nat Commun*. doi:10.1038/s41467-026-76240-w — <https://pmc.ncbi.nlm.nih.gov/articles/PMC13522577/>
2. Akbari et al. 2022. Genome-wide detection of imprinted DMRs using nanopore sequencing. *eLife* 11:e77898 — <https://elifesciences.org/articles/77898>
3. Rosenski et al. 2025. Atlas of imprinted and allele-specific DNA methylation in the human body. *Nat Commun* 16:2141 — <https://pmc.ncbi.nlm.nih.gov/articles/PMC11897249/>
4. Gustafson et al. 2024. High-coverage nanopore sequencing of 1000 Genomes Project samples. *Genome Res* 34:2061 — <https://genome.cshlp.org/content/34/11/2061>
5. Zink et al. 2018. Insights into imprinting from parent-of-origin phased methylomes and transcriptomes. *Nat Genet* 50:1542 `[A]` — <https://www.nature.com/articles/s41588-018-0232-7>
6. Abante, Fang, Feinberg, Goutsias 2020. Detection of haplotype-dependent allele-specific DNA methylation in WGBS data (CPEL). *Nat Commun* — <https://pmc.ncbi.nlm.nih.gov/articles/PMC7567826/>
7. Akbari et al. 2021. Megabase-scale methylation phasing using nanopore long reads and NanoMethPhase. *Genome Biol* — <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02283-5>
8. Tian et al. 2026. nanoASM: long-read allele-specific DNA methylation profiling in human prostate tissues `[P]` — <https://pmc.ncbi.nlm.nih.gov/articles/PMC13320949/>
9. Tian et al. 2025. Fine mapping regulatory variants by characterizing native CpG methylation with nanopore long-read sequencing. *HGG Adv* 7:100532 — <https://pmc.ncbi.nlm.nih.gov/articles/PMC12642126/>
10. Jensen et al. 2026. METAFORA: population-scale detection of methylation outliers from long-read genome sequencing. *medRxiv* `[P]` — <https://pmc.ncbi.nlm.nih.gov/articles/PMC13278280/>
11. Li et al. 2026. cyberDMR: accurate and robust identification of DMRs from WGS-derived methylomes. *BMC Biol* 24:165 — <https://pmc.ncbi.nlm.nih.gov/articles/PMC13371177/>
12. Groza et al. 2025. Expanded methylome and QTL detection by long-read profiling of personal DNA. *Genome Res* — <https://pmc.ncbi.nlm.nih.gov/articles/PMC12047246/>
13. Zhuo et al. 2026. Characterizing cytosine methylation of polymorphic transposable element insertions using the human pangenome resources. *Genome Res* — <https://pmc.ncbi.nlm.nih.gov/articles/PMC13262945/>
14. NANOME 2025. Nextflow pipeline for haplotype-aware allele-specific consensus DNA methylation detection — <https://pmc.ncbi.nlm.nih.gov/articles/PMC12236756/>
15. ASMS 2026. Finding allele-specific methylation in human genomes without phasing — <https://pmc.ncbi.nlm.nih.gov/articles/PMC13326645/> · <https://github.com/ecmra/asms>
16. A comprehensive DNA methylation atlas for the Chinese population through nanopore long-read sequencing of 106 individuals. *bioRxiv* 2026 `[P]` `[A]` — <https://www.biorxiv.org/content/10.64898/2026.04.20.719515v1.full>
17. MethPhaser 2024. Methylation-based long-read haplotype phasing of human genomes. *Nat Commun* `[A]` — <https://www.nature.com/articles/s41467-024-49588-0>
18. Methylation-aware long-read phasing significantly improves genome-wide haplotype reconstruction. *bioRxiv* 2026 `[P]` `[A]` — <https://www.biorxiv.org/content/10.64898/2026.03.11.710820v2.full>
19. Guo et al. 2017. Identification of methylation haplotype blocks aids in deconvolution of heterogeneous tissue samples. *Nat Genet* `[A]` — <https://pubmed.ncbi.nlm.nih.gov/28263317/>
20. A DNA methylation haplotype block landscape in human tissues and preimplantation embryos. 2023 — <https://pmc.ncbi.nlm.nih.gov/articles/PMC10760529/>
21. Genome-wide single-molecule analysis of long-read DNA methylation reveals heterogeneous patterns at heterochromatin that reflect nucleosome organisation. *PLoS Genet* 2023 — <https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010958>
22. Voronina et al. 2025. Interplay of genetic variants and allele-specific methylation in a single human genome. *IJMS* — <https://pmc.ncbi.nlm.nih.gov/articles/PMC12525282/>
23. Haplotype-resolved methylomes reveal parent-of-origin DNA methylation imbalance in autism spectrum disorder. *Sci Adv* 2026 `[A]` — <https://www.science.org/doi/10.1126/sciadv.aee4069>
24. Haplotype-resolved DNA methylation at the APOE locus. *npj Dementia* 2026 `[A]` — <https://www.nature.com/articles/s44400-026-00094-8>
25. Allele-specific methylation uncovers non-Mendelian inheritance (mouse). *Nat Genet* 2026 `[A]` — <https://www.nature.com/articles/s41588-026-02603-0>
26. Measuring X-chromosome inactivation skew for X-linked diseases with adaptive nanopore sequencing. *Genome Res* 2024 — <https://pmc.ncbi.nlm.nih.gov/articles/PMC11610576/>
27. Nanopore-based haplotype-resolved XCI for severity assessment in X-linked disorders. 2026 — <https://pmc.ncbi.nlm.nih.gov/articles/PMC13347936/>
28. Lineage-specific X chromosome inactivation escape and skew. *bioRxiv* 2026 `[P]` `[A]` — <https://www.biorxiv.org/content/10.64898/2026.08.26.739472v1.full>
29. HPRC Data Release 2 — <https://humanpangenome.org/hprc-data-release-2/> · Epigenomes of HPRC Release 2 (AWS Open Data) — <https://registry.opendata.aws/hprc-epigenome/>
30. HPRC2: a human pangenome reference with near-complete coverage of common genetic variation. 2026 `[A]` — <https://pubmed.ncbi.nlm.nih.gov/42539208/>
31. Methylartist: tools for visualizing modified bases from nanopore sequence data. *Bioinformatics* 2022 `[A]` — <https://academic.oup.com/bioinformatics/article/38/11/3109/6575433>
32. DNA methylation-calling tools for Oxford Nanopore sequencing: a survey and human epigenome-wide evaluation. 2021 — <https://pmc.ncbi.nlm.nih.gov/articles/PMC8524990/>
33. Methylmap: visualization of modified nucleotides for large cohort sizes. *BMC Bioinformatics* 2025 `[A]` — <https://link.springer.com/article/10.1186/s12859-025-06106-3>
