# The HPRC Epigenome resource (232 samples): what it contains and what it changes

**Date:** 2026-09-08. Addendum to [`20260908_lit_review.md`](20260908_lit_review.md).
**Resource:** <https://epigenome.humanpangenome.org/> (Ting Wang lab, WashU) ·
bucket `s3://hprc-epigenome` (us-east-2, public, no-sign-request) ·
[AWS Open Data registry](https://registry.opendata.aws/hprc-epigenome/)

Everything below was verified by direct S3 listing and byte-range reads of the actual files on
2026-09-08, not from the portal's documentation.

---

## 1. Headline

**221 individuals have both haplotype-resolved ONT methylation and 1000 Genomes phased
genotypes.** Not 18. The genotypes are from the *same 3202-sample phased panel this project's
W01 already uses*.

This resolves the single biggest structural weakness identified in the literature review — that
n=18 is indefensible next to 1KGP-ONT's 226 and deCODE's 7,284 — and it makes an
ancestry-stratified, genotype-linked analysis supportable for the first time.

It also raises the scooping risk (§6), and it does **not** remove the calibration work that
currently blocks the project (§7).

---

## 2. What is actually in the bucket

232 sample directories, 11,355 objects. Per-sample contents, with the count of samples that
have each file:

| File | n/232 | What it is |
|---|---|---|
| `methylation.ONT.hap{1,2}.modbed.gz` | **229** | **Read-level** ONT methylation, per haplotype |
| `methylation.PacBio.hap{1,2}.methylc.gz` | 221 | Aggregate per-CpG PacBio methylation, per haplotype |
| `expression.{plus,minus}.hap{1,2}.{combined,specific}.bw` | 202 | **Allele-specific expression** (Kinnex RNA), stranded, per haplotype |
| `hg38.hic`, `chm13.hic` + insulation/boundary bedgraphs | 214 | Hi-C contact maps and TAD boundaries |
| `hap{1,2}_vs_{hg38,chm13}.gz` | 232 | Assembly-to-reference alignments (projection) |
| `CGI.bed.gz`, `RepeatMasker.hap{1,2}.bb`, `HMMFlagger.{ONT,PacBio}.bed.gz` | 232 | CpG islands, repeats, assembly-reliability flags |
| `hap{1,2}.refbed.gz` | 232 | Gene annotation on assembly coordinates |

Missing ONT methylation: `HG002`, `HG02622`, `NA19087`.

### The ONT methylation is read-level — this matters enormously

```
HG00146#1#CM090013.1  29  8878  a3c68340-2c7f-448b-bfb6-0d01508ab901  0  -  <mod offsets>  <unmod offsets>
```

Column 4 is an **ONT read UUID**. Each line is one molecule, with its methylated and unmethylated
CpG offsets. Epiallele structure, co-methylation, entropy, per-read confirmation, and the
design-effect/pseudoreplication analysis are all computable directly from these files. ~2.6 GB
per haplotype per sample → **~1.2 TB for the full read-level ONT set.**

### The PacBio methylation is *not* count data

```
HG00146#1#CM090013.1  9439  9440  CG  0.953  +  4
```

A methylation value of 0.953 at coverage 4 is not a ratio — it is a **model score** (pb-CpG-tools
model mode). **The beta-binomial caller cannot consume these files**, because there is no
(n_modified, n_total) pair to put in a binomial likelihood. Treat PacBio here as an orthogonal
concordance check only, not as an inference substrate.

### Coordinates are per-sample assembly, not GRCh38

Contigs are PanSN-named: `SAMPLE#HAPLOTYPE#CONTIG` (e.g. `HG00146#1#CM090013.1`). The
`hap1_vs_hg38.gz` files are WashU-browser `genomealign` records, not standard liftOver chains:

```
chr1  10336  10800  id:1,genomealign:{chr:"HG00146#1#JBHIKO010000053.1",start:84129873,stop:84130268,strand:"-",targetseq:"..."}
```

**Consequence:** to link a 1KGP genotype (GRCh38) to a methylation value (assembly coordinates)
you must project across 464 haplotype assemblies. Regions that fail to align do not fail at
random — they are exactly the repetitive, structurally variable regions where long-read ASM is
most interesting. This is a real, non-trivial source of ascertainment bias.

---

## 3. Genotype availability — verified

Intersected the 232 sample IDs against the 1KGP 3202-sample high-coverage pedigree
(`1kGP.3202_samples.pedigree_info.txt`, `20130606_g1k_3202_samples_ped_population.txt`):

- **223 of 232 (96.1%)** are in the 1KGP 3202 panel → phased genotypes directly available, from
  the panel `scripts/download/D02_*` already fetches.
- 9 are not: `HG002`, `HG005` (GIAB — have their own truth sets), `HG01123`, `HG02109`,
  `HG02486`, `HG02559`, `HG03471`, `HG06807`, `NA21309`.
- **Analyzable cohort = 221** (both ONT haplotype modbeds *and* 3202-panel genotypes).
  Sample list saved during this analysis; regenerate with the S3 listing + pedigree intersect.

### Ancestry composition of the n=221 analyzable cohort

| Superpopulation | n |
|---|---|
| AFR | **63** |
| EAS | 49 |
| AMR | 43 |
| SAS | 36 |
| EUR | 30 |

**25 distinct 1KGP populations.** Sex: 113 male, **108 female**.

Two things to notice. First, this fixes the documented AFR-thin problem
(`20260904_cohort_provenance.md`: AFR at n=2, no ancestry-stratified claim supportable). AFR at
n=63 is supportable. Second, the cohort is **AFR-majority and EUR-minority** — the inverse of
almost every published methylation cohort, including deCODE (Icelandic) and Blueprint. That is a
genuine and citable strength, not a footnote.

### Relatedness

**Zero parent–offspring pairs within the 232.** Checked every sample's father/mother ID against
the set; no first-degree relatives are both present. Cryptic relatedness and genotype PCs still
need checking before any population-level claim, but there is no trio structure to unwind.

### Power

At n=221, expected heterozygote carriers per site:

| MAF | Expected hets |
|---|---|
| 0.05 | 21 |
| 0.10 | 40 |
| 0.20 | 71 |
| 0.50 | 110 |

**This is well-powered for ASM, and the reason is structural, not just sample size.** ASM is a
*within-individual* contrast: each heterozygote carries both alleles in the same nucleus, same
cell population, same batch, same environment. The between-individual confounders that force
conventional meQTL studies to n≈1,000+ (cell composition, age, batch, environment) are
differenced out. 71 heterozygotes contributing paired allelic observations is a strong test.
For comparison, Blueprint's per-cell-type meQTL analyses ran at roughly this n and yielded tens
of thousands of meQTLs.

Conventional population meQTL at n=221 is more modest — adequate for large-effect cis signals,
underpowered for trans and for small effects. **Frame the design as ASM-first with meQTL as
corroboration, not as a meQTL study.**

---

## 4. What this unlocks that was not previously available

1. **Ancestry-stratified ASM.** Supportable at AFR 63 / EAS 49 / AMR 43 / SAS 36 / EUR 30.
2. **ASM × allele-specific expression in the same individuals** (202 samples with per-haplotype
   RNA). Does haplotype-biased methylation predict haplotype-biased expression, on the same
   molecules-worth of the same donor? No paper found in the literature review does this at scale
   in a long-read cohort. **This is a stronger novelty claim than anything in the current plan.**
3. **ASM × 3D genome** (214 samples with Hi-C and insulation scores). Distal ASM-driving SNPs
   and TAD boundaries in matched samples directly tests the "does the distal SNP contact the
   CpG" question that the proximal/distal taxonomy raises but cannot answer.
4. **XCI at n=108 females**, rather than the 2-donor pilot in `scripts/chrX/README.md`.
5. **Allele frequency as a covariate.** With genotypes, ASM can be tested against MAF, LD
   structure, and derived-allele status — actual population genetics, not just cataloguing.

---

## 5. The critical decision: use this resource's files, or keep the existing pipeline?

**Recommendation: keep the existing pipeline for the ASM calling; use the resource for the sample
manifest and the other assays.**

Reasoning:

- The pipeline already aligns HPRC ONT uBAMs to GRCh38 and haplotags with WhatsHap against the
  3202 panel. That produces **GRCh38 coordinates and genotype-linked haplotype labels in one
  step** — precisely the two things the resource's files do *not* give you.
- The resource's `hap1`/`hap2` are **assembly haplotypes**: arbitrary labels, not parent-of-origin,
  and with no inherent link to which allele each carries. Recovering that link requires
  assembly-based variant calls (`D03`/`D04`/dipcall — already written but only run as a pilot)
  plus projection through the genomealign files. That is strictly more work than the existing
  route, and it introduces the alignment-dropout bias described in §2.
- The one thing the resource offers that the pipeline does not is **the other assays** (ASE,
  Hi-C, gene annotation) already harmonized per haplotype. Those are worth taking as-is.

The honest cost accounting for scaling the existing pipeline from 18 to ~221 donors: 221 uBAM
downloads, alignments, merges, haplotagging, `modkit extract`/`pileup`, at roughly 12× the
current compute and storage. That is a large but ordinary job — **and it should not start until
§7 is resolved.**

---

## 6. Scooping risk, revised upward

The literature review ranked "HPRC methylation working group publishes first" as risk #1. This
resource confirms that risk and sharpens it: the Wang lab has **already assembled, harmonized,
and published** per-haplotype methylation, per-haplotype expression, and Hi-C for 232 samples.
Every input for "allele-specific methylation and expression across the human pangenome" is in
their hands, released as of 2025-11-30.

They may well publish only a resource/browser paper — that is the most common outcome for a data
portal, and a resource paper typically does not do calibrated per-locus ASM inference. But
planning should assume the descriptive layer ("here is haplotype-biased methylation across 232
diverse genomes, here is its overlap with ASE") is not available to claim.

**What stays defensible regardless of what they publish:** the read-level statistics
(design-effect/pseudoreplication correction, per-donor dispersion, calibrated per-locus FDR), the
controlled phasing-window experiment quantifying long-read gain, and per-genetic-haplotype
epiallele architecture. None of those is a resource-paper analysis. Concentrate there.

---

## 7. What has to happen before scaling

Scaling 18 → 221 multiplies every unresolved calibration problem by twelve. In priority order:

1. **Resolve the three-way inconsistency** in the proximal/distal breakdown
   (`20260908_lit_review.md` §7). Running 221 donors through a caller whose headline output has
   three incompatible values produces 221 donors' worth of ambiguity.
2. **Root-cause the zero-hit donors** (HG00344, NA21144; ρ̂ 0.097/0.116). At n=18 these are two
   anomalies. At n=221 the same failure mode is ~25 donors, and a reviewer will read it as
   uncontrolled.
3. **Fix the haplotagging rate** (~38% vs. the ≥60% gate). MethPhaser and the 2026
   methylation-aware phasing preprint are the citable remedies. At scale this is the difference
   between a usable and an unusable cohort.
4. **Then** pilot the scale-up on ~40 donors balanced across the five superpopulations before
   committing to all 221.

A caveat that does not go away with n: these are **LCL-derived samples**. At n=221 with matched
ASE, LCL clonality becomes testable rather than merely admitted — clonal expansion should produce
correlated allelic skew across methylation *and* expression, which is a checkable signature. Do
that check; it converts a reviewer objection into a result.

---

## 8. Revised verdict

The literature review concluded the project should be repositioned as a methods-plus-depth paper
because n=18 could not support population claims. **That constraint is gone.** With n=221,
genotypes, 25 populations, AFR-majority composition, 108 females, and matched allele-specific
expression and Hi-C, a genuinely population-oriented project is well-powered and well-founded.

The recommended shape changes accordingly:

1. Calibrated read-level ASM caller with explicit design-effect correction (unchanged — still
   the methods contribution, and now applied at a scale that justifies it).
2. Genome-wide ASM across 221 diverse genomes, with proximal/distal architecture, tested against
   allele frequency and LD.
3. **ASM × ASE concordance in matched individuals** — promoted to a primary chapter.
4. Ancestry-stratified ASM, which is now supportable.
5. The controlled phasing-window experiment quantifying long-read gain (unchanged — still the
   clearest open gap in the literature).

XCI at n=108 and ASM × Hi-C are strong secondary chapters; the QC-heterogeneity material remains
a discussion section.

**Reproducibility note.** All numbers in this document were derived from the public bucket and
the public 1KGP pedigree; the derivation is a bucket listing plus an ID intersect and is
reproducible in a few minutes with `curl` alone.

---

## 9. Addendum: access policy, what HPRC has actually done, and how little you need to download

Added the same day, after the questions "am I scooping a consortium?" and "do I need all the
BAMs?".

### 9.1 The data-use policy explicitly permits this

From the [HPRC Data Use and Publication Protocol](https://humanpangenome.org/publication_protocol/)
and [Data Use](https://humanpangenome.org/data-use/) pages:

- HPRC data "are in the public domain," not patented or copyrighted, and donors "have given broad
  informed consent for re-use."
- "Publicly released" = data that is either published, **or** part of an unpublished data freeze
  more than one year old.
- **"Researchers are encouraged to publish publicly released data without contacting the HPRC
  directly."**
- The notify-first requirement applies only to genome-wide publication on data that is *not yet*
  publicly released; in that case HPRC asks for consortium-banner authorship after steering
  committee review.
- Users must cite the latest integrated HPRC publication and the accessions used, and remain
  bound by any sample-specific restrictions in the HPRC data use table.

Release 2 assemblies are public, the epigenome resource is published on AWS Open Data with a
prescribed citation format, and the HPRC2 preprint (July 2026) exists as the integrated
publication to cite. This is not scooping; it is the intended use.

**Recommended anyway: email the resource contact (`dli23@wustl.edu`, Ting Wang lab) before
committing months of compute.** This is the dominant strategy regardless of the answer. Either
they say they are not pursuing ASM — which removes the largest risk in this project at the cost
of one email — or they are, and the conversation becomes a collaboration instead of a race lost
silently. The downside is negligible because the policy already grants permission.

### 9.2 What HPRC has and has not done

Searched again specifically for consortium ASM/meQTL work. Findings:

- **No HPRC ASM or meQTL paper exists.** The HPRC-adjacent methylation work remains Zhuo et al.
  2026 (*Genome Res*), scoped to polymorphic transposable-element insertions.
- **The Fiber-seq worry is partly misplaced.** Enumerating all 11,355 objects in the bucket, the
  per-sample assay files are: ONT methylation, PacBio methylation, per-haplotype expression,
  Hi-C, CGI, RepeatMasker, HMMFlagger, gene annotation, and assembly-to-reference alignments.
  **No Fiber-seq / FIRE / ATAC files are present.** The AWS registry's "chromatin accessibility"
  phrasing appears to refer to the Hi-C tracks and to the browser's broader scope.
- **Fiber-seq haplotype work exists but elsewhere.** "A haplotype-resolved view of human gene
  regulation" (Stergachis lab, bioRxiv 2024) maps haplotype-selective chromatin accessibility
  (1,231 genome-wide-significant elements) using FIRE — on GM12878, COLO829, fibroblast lines,
  CD8+ T cells and primary tissue, **not HPRC samples**, and it explicitly does **not** perform
  ASM or meQTL discovery even though it measures CpG methylation alongside.

That last point is worth reading twice. The nearest thing to this project's design that has been
published is single-modality (accessibility only), on a handful of cell lines, and it left
methylation on the table.

### 9.3 A methylation-only paper does not look thin — a shallow one does

The instinct that "if I only do methylation it'll look weird next to their Hi-C and RNA" conflates
two different paper shapes. Multi-omic integration across every available assay is a *consortium
resource* paper. A single-modality paper carried by methodological depth is a different and
equally standard shape — and it is what the comparison set actually looks like: the FIRE paper is
accessibility-only, Rosenski is methylation-only, nanoASM is methylation-only.

The correct response to a data-rich neighbourhood is not to add modalities until the scope becomes
unmanageable. It is to go deeper on one. Adding **ASE as a single validation chapter** is cheap
(202 samples, bigwigs, no new pipeline) and buys most of the integrative credit; Hi-C is a third
chapter only if the distal-ASM result actually needs it.

### 9.4 You do not need the BAMs

**Verified working, remotely, with no bulk download.** Every methylation and alignment file in the
bucket is tabix-indexed and readable over HTTPS by `tabix` directly:

```
B=https://hprc-epigenome.s3.us-east-2.amazonaws.com/samples/HG00146

# step 1 — project an hg38 region to this sample's assembly coordinates
tabix $B/hap1_vs_hg38.gz chr7:130490000-130500000
#   -> chr7 129808123 130591387 id:33866,genomealign:{chr:"HG00146#1#CM090016.1",
#      start:130815525,stop:131598785,strand:"+",targetseq:"..."}

# step 2 — pull read-level methylation there
tabix $B/methylation.ONT.hap1.modbed.gz "HG00146#1#CM090016.1:131490000-131500000"
#   -> 30 reads, 1.3 seconds
```

What this means for scope:

| Need | Source | BAM required? |
|---|---|---|
| Per-read CpG patterns (epialleles, entropy, co-methylation) | `methylation.ONT.hap*.modbed.gz` — read UUIDs + per-CpG offsets | **No** |
| Per-CpG per-haplotype counts for the beta-binomial caller | derive by aggregating the same modbed | **No** |
| Haplotype assignment of reads | already done, assembly-based | **No** |
| Design-effect / pseudoreplication measurement | read-level modbed | **No** |
| Het sites and which haplotype carries which allele | dipcall on the two assemblies (`D03`/`D04` already written), or 1KGP panel + projection | No — assemblies, not BAMs |
| Coordinate projection to hg38 | `hap*_vs_{hg38,chm13}.gz`, tabix-indexed | **No** |
| Assembly-reliability masking | `HMMFlagger.{ONT,PacBio}.bed.gz` | **No** |

Storage: full download of read-level ONT methylation is ~1.2 TB (2.6 GB × 2 × 229), versus roughly
11–22 TB for the equivalent BAMs. With region-restricted remote queries it is closer to **zero** —
stream per sample, compute, discard.

**Two things you genuinely give up by not using BAMs**, and both should be stated as limitations
rather than discovered later:

1. **Per-call modification probabilities are gone.** The modbed encodes modified and unmodified
   CpG positions as offset lists — already binarized at somebody else's threshold. There is no ML
   tag to filter on and no way to propagate call uncertainty into the beta-binomial likelihood.
   For a project whose contribution is statistical calibration, that is a real constraint. It may
   justify keeping BAMs for the read-level Tier 2 subset only.
2. **No re-alignment or re-basecalling.** The haplotype assignment, the basecaller, and the
   modification model are fixed to HPRC's choices. This is mostly a benefit (standardization
   across 229 donors is exactly the thing the current 18-donor cohort lacks), but it removes the
   R10.4.1/Dorado-vs-R9.4.1 control this project currently exercises deliberately.

**Also note what assembly-based haplotype assignment fixes for free:** the ~38% WhatsHap
haplotagging rate, which is currently one of the top three blockers, is not a problem in this data
— reads were assigned to haplotypes by mapping to the diploid assembly, not by spanning a
phased het SNP. That removes an entire open item.

### 9.5 Unexpected observation: these are ultra-long reads

In the MEST-region query above, the 6 reads inspected had alignment spans of **130–278 kb**, one
carrying 2,111 methylated and 711 unmethylated CpG calls on a single molecule. HPRC ONT is
explicitly the "Ultralong" product (>100 kb), which is consistent.

Two consequences, both flagged as *observed on 6 reads in one window of one sample — compute the
cohort-wide read-length distribution before relying on this*:

- `README.md` currently states "Long reads span 5–20 kb." If the cohort-wide distribution matches
  what was observed, that line understates the data by an order of magnitude and should be fixed.
- The proximal/distal/very-distal taxonomy is capped at 50 kb, which was calibrated to 5–20 kb
  reads. Molecules spanning 130–278 kb make SNP–CpG linkage observable at distances the current
  classification cannot express, and they differentiate this data from both 1KGP-ONT (N50 54 kb)
  and deCODE (N50 19 kb) on the one axis this project's core claim depends on.

### 9.6 A tiered scope that keeps the project manageable

The worry that n=221 makes the project unmanageable is legitimate, but n=221 does not have to mean
221 donors through every analysis. Availability of a cohort does not obligate its use — but with
the data this easy to query, "why only 30?" now needs an answer, and a tiered design is that
answer:

| Tier | n | Analyses | Cost |
|---|---|---|---|
| 1 | ~221 (all) | Aggregate per-CpG per-haplotype ASM; proximal/distal architecture; ancestry stratification; allele-frequency and LD tests | Streaming remote queries; no BAMs; modest |
| 2 | ~40–60, balanced across the 5 superpopulations | Read-level epiallele structure, co-methylation, entropy, design-effect calibration, per-locus confirmation | Read-level; the expensive tier; possibly with BAMs for call probabilities |
| 3 | ~202 (ASE available) | ASM × allele-specific expression concordance | Cheap — bigwigs, no new pipeline |

This is defensible as written: the expensive analysis runs where depth is needed, the cheap
analysis runs everywhere, and no reviewer asks why the cohort was truncated.

**The sequencing constraint from §7 still holds.** None of this starts before the three-way
proximal/distal inconsistency is resolved and the zero-hit-donor failure mode is root-caused —
those are caller problems, and they get 12× worse at scale, not better.
