---
Assessment against: Rosenski et al. (2025) "Atlas of imprinted and allele-specific DNA methylation
in the human body." Nature Communications 16:2141.
Updated: 2026-05-12
---

## Assessment of Original Proposal

The Rosenski atlas directly validates the core premise of this project. Using deep WGBS across 39
cell types, they identified 324,759 bimodal methylation loci — but could associate only ~34,000
(~10%) with a nearby SNP. The remaining ~290,000 bimodal regions are unexplained: the authors
explicitly acknowledge that "short-read sequencing (~200 bp) cannot capture distant regulatory
variants" and that it is impossible to distinguish sequence-dependent from stochastic methylation
without informative SNPs. Long-read sequencing is the direct answer to this limitation.

The proposal is well-motivated. The main gaps relative to the article:
- Does not leverage the atlas as a pre-computed target set (saves de novo calling)
- The epiallele clustering framing is underspecified — the Rosenski EM approach already handles
  bimodal calling; the proposal needs to articulate what long reads add *beyond* bimodality
- No treatment of parental/imprinted ASM, which the article elevates as a major category
- Long-range phasing — the key technical contribution of long reads — is not named explicitly

---

## Background

Allele-specific methylation (ASM) is traditionally defined using heterozygous SNPs to phase
sequencing reads and compare methylation levels between haplotypes. This captures genotype-driven
cis effects but reduces methylation to a single per-allele summary statistic.

Epigenomic data at single-molecule resolution reveals richer structure: individual reads contain
multi-CpG methylation patterns (epialleles) reflecting coordinated states along a molecule. These
arise from genetic variation, cell-state heterogeneity, or stochastic regulation.

The Rosenski et al. (2025) atlas provides the most comprehensive characterization of bimodal
methylation to date (325K loci, 39 cell types, deep WGBS), but is fundamentally limited by
short reads:
- SNP-ASM detection is restricted to SNPs within ~200 bp of the bimodal region
- ~90% of bimodal loci have no associated SNP explanation
- The method cannot distinguish stochastic from genotype-driven bimodality at unexplained loci
- Multi-allelic (>2 epiallele) patterns are forced into a bimodal EM model

Long-read sequencing (nanopore) addresses all four limitations simultaneously by jointly observing
phased SNPs and dense CpG patterns along kilobase-scale molecules.

---

## Core Question

To what extent is single-molecule methylation structure explained by genetic variation — and how
far does that genetic influence extend along the chromosome?

The Rosenski atlas defines *where* bimodality exists. This project asks *why*.

---

## Proposed Project

### Goal

Use long-read sequencing to resolve the genetic architecture of bimodal methylation loci,
specifically to:
1. Classify unexplained bimodal loci as genotype-driven (distal cis effects), parental, or
   stochastically regulated
2. Characterize long-range haplotype–epiallele coupling beyond the 200 bp short-read limit
3. Detect multi-allelic and continuous epiallele structures invisible to WGBS bimodal calling

---

### Approach

#### 1. Target the atlas, don't re-derive it

Use the 324,759 bimodal regions from Rosenski et al. as the input locus set. Prioritize:
- The ~290K bimodal loci with no associated SNP (highest information gain)
- The 347 novel parental-ASM loci (functionally uncharacterized)
- Known imprinted loci with tissue-specific escape (e.g., IGF2, IGF2R, CHD7)

This eliminates the de novo discovery step and focuses effort on interpretation.

#### 2. Long-range phasing

For each bimodal locus, extend SNP-based phasing to a ±5–10 kb window using nanopore reads
(vs. 200 bp for WGBS). For each read overlapping a bimodal locus:
- Assign reads to haplotypes using heterozygous SNPs anywhere in the long-read window
- Extract per-read CpG methylation pattern (epiallele)
- Compute haplotype–epiallele alignment: mutual information or AUC of haplotype predicting
  methylation state

This directly tests whether short-read SNP-ASM underestimates the fraction of genotype-driven
bimodal loci.

#### 3. Classify loci by genetic architecture

For each bimodal locus, classify into one of four regimes:

| Class | Signature |
|---|---|
| Proximal genotype-driven | SNP within 200 bp predicts epiallele (replicates Rosenski) |
| Distal genotype-driven | SNP at 200 bp–10 kb predicts epiallele (new with long reads) |
| Parental (imprinted) | Epiallele correlates with parent-of-origin, not SNP identity |
| Non-genetic | No haplotype–epiallele alignment; entropy is symmetric across haplotypes |

For the parental class: trio-sequenced samples (child + parents) allow direct parent-of-origin
phasing without relying on allele-flipping across families (as Rosenski required).

#### 4. Epiallele complexity beyond bimodality

The Rosenski EM model forces bimodality (2 components). Long reads can test this:
- Apply a model-order selection criterion (BIC or likelihood ratio) to choose between 1, 2,
  or ≥3 epiallele components per locus
- Characterize the CpG pattern of each component (not just mean methylation)
- Identify loci with genuinely continuous or multi-allelic methylation distributions

This is only possible with long reads, where the full within-read pattern is observable.

#### 5. Case studies with functional grounding

Draw from the atlas's specific findings to anchor the analysis:

- **IGF2 hepatocyte escape**: Rosenski identified two liver-specific unmethylated enhancers
  enabling biallelic IGF2 expression. Long reads in hepatocyte-enriched material can test whether
  the enhancer methylation state is phased to the normally-imprinted allele or represents a
  bona fide epiallele switch.
- **CHD7 intron / CHARGE syndrome**: A novel maternal-methylated iDMR adjacent to a CTCF site.
  Long reads span the CTCF binding site and the iDMR simultaneously, allowing direct chromatin
  loop–methylation co-occurrence testing (if CTCF ChIP or CUT&RUN is available).
- **IGF2R colon-specific biallelic methylation**: This locus flips from monoallelic to biallelic
  methylation in colon epithelium. Long reads can resolve whether this is a continuous graded
  shift or a discrete epiallele switch.

---

### New Ideas (Beyond Original Proposal)

#### A. Quantify short-read underestimation of ASM

Compare:
- Fraction of bimodal loci with proximal SNP-ASM in matched short-read data
- Fraction explained when extending to a ±5 kb phasing window with long reads

This produces a concrete estimate of how much ASM the field has missed due to read-length
constraints — a publishable result in itself.

#### B. Epiallele entropy as a cell-type marker

Rosenski et al. showed bimodal loci are largely shared across cell types (65% present in ≥10
samples). But epiallele entropy — the balance between the two allelic states — may be
cell-type-specific even when bimodality is shared. Map epiallele entropy across available
cell-type-matched public long-read datasets (e.g., ONT data in ENCODE or HPRC).

#### C. Long-range regulatory variant discovery

For loci where phasing identifies a distal SNP as the driver, test whether the SNP overlaps:
- Cell-type-specific enhancers (H3K27ac, ATAC peaks)
- CTCF binding sites (known to anchor methylation boundaries)
- Known eQTLs (GTEx)

This directly links the methylation architecture to regulatory function and generates candidate
causal variants for unexplained ASM loci.

#### D. Multi-CpG pattern clustering (beyond mean methylation)

The Rosenski EM model uses a single parameter per epiallele (mean methylation). Long reads allow
clustering reads by their full CpG vector (e.g., {M, U, M, M, U} per read). Use:
- UMAP or PCA on per-read methylation vectors within each bimodal locus
- Gaussian mixture models on the full pattern (not just the mean)

This may reveal that some "bimodal" loci are actually structured mixtures of distinct methylation
programs that average to look bimodal.

---

### Feasibility

- **Minimum viable dataset**: 1–2 deeply sequenced nanopore genomes (ideally a trio) mapped to
  the Rosenski bimodal locus set
- **Public data option**: HPRC long-read data includes trio-sequenced samples with existing
  methylation calls (e.g., Dorado modkit output)
- **Core analysis scope**: Long-range phasing + locus classification is implementable in ~3
  months; case studies and cell-type entropy analysis add another 2 months
- **No de novo peak calling needed**: Rosenski locus coordinates are public (GEO GSE186458)

---

### Expected Contributions

1. A genome-wide classification of bimodal methylation loci by genetic architecture (proximal,
   distal, parental, non-genetic) — the first such map at long-read resolution
2. A quantitative estimate of ASM underestimation by short-read methods
3. A catalog of multi-allelic and pattern-structured loci invisible to WGBS bimodal calling
4. Mechanistic grounding for specific disease-relevant imprinted loci (IGF2, CHD7, IGF2R)

---

### Relation to Rosenski et al.

This project is complementary, not redundant. Rosenski provides the *where* (bimodal locus atlas,
39 cell types, deep coverage). This project provides the *why* (genetic architecture, long-range
phasing, epiallele pattern structure) for the ~90% of bimodal loci they could not explain.

---

## Key Changes from Original Proposal

1. **Use the Rosenski atlas as input, not a starting point for de novo discovery** — 325K
   pre-called bimodal loci are a better target set than anything called from scratch. This
   eliminates the detection step and focuses effort on interpretation.

2. **Name long-range phasing explicitly as the core contribution** — extending SNP-ASM phasing
   from 200 bp (WGBS limit) to ±5–10 kb is the key technical advance. The original proposal did
   not foreground this.

3. **Four-class locus taxonomy** — proximal genotype-driven, distal genotype-driven, parental,
   non-genetic. Rosenski's unexplained ~290K bimodal loci are the primary targets for
   reclassification.

4. **Quantify short-read ASM underestimation** — comparing proximal vs. extended phasing
   fractions produces a concrete, publishable estimate of what WGBS misses.

5. **Epiallele entropy as a cell-type marker** — even when bimodality is shared across cell
   types, the balance between allelic states may vary. Testable with public ONT datasets.

6. **Multi-CpG pattern clustering beyond bimodal EM** — long reads support model-order selection
   (1, 2, or ≥3 epiallele components) and full-pattern clustering, detecting structure that
   WGBS forces into a bimodal model.

7. **Long-range regulatory variant discovery** — distal SNPs driving ASM can be linked to
   enhancers, CTCF sites, and eQTLs, grounding methylation architecture in function.

8. **Concrete case studies from the paper** — IGF2 hepatocyte escape, CHD7/CHARGE syndrome
   maternal iDMR, IGF2R colon biallelic flip — provide functional grounding and disease
   relevance the original proposal lacked.
