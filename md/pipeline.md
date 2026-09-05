# ASM Long-Read Pipeline

**Sample**: HG01258 (1000 Genomes, HPRC trio)  
**Data**: Nanopore ONT, Guppy 6.5.7 basecaller, MM/ML 5mCG methylation tags  
**Reference**: GRCh38 autosomal

---

## Scientific framing

The central question is: **to what extent does single-molecule methylation structure arise from genetic variation?**

Long reads jointly observe phased heterozygous SNPs and dense CpG methylation on the same molecule. This lets you:
1. Assign reads to haplotypes via phased SNPs
2. Measure methylation per CpG per haplotype
3. Test whether the two haplotypes differ at each locus (allele-specific methylation, ASM)
4. Classify ASM by how far its driving SNP is: proximal (<200 bp, detectable by short reads) vs. distal (200 bp–50 kb, long-read-only) vs. non-genetic (no phased SNP linkage at any distance)
5. Characterize single-molecule methylation *patterns* (epialleles) beyond mean differences — are individual molecules uniformly intermediate, or do they fall into discrete all-methylated/all-unmethylated states?

This is not "conditioning on genetic ASM and asking if epialleles differ." It is the reverse: measure all haplotype-stratified methylation structure first, then ask how much of it is explained by genetics. Genetic ASM is a subset of the output, not a filter applied at the start.

The key innovation over short-read ASM studies: short reads can only detect SNP-CpG associations within ~200 bp (one read). Long reads span 5–20 kb, phasing methylation over distances that are invisible to Illumina.

---

## Pipeline overview

```
ONT BAM (HPRC S3)
    │
    ├── A: Chr19 proof-of-concept (low coverage; framework validation)
    │
    └── W: Whole-genome production (all 22 autosomes)
         │
         ├─ STAGE 1: Align + phase
         │       BAM → minimap2 → haplotagged BAM
         │       1000G VCF → WhatsHap → phased VCF
         │
         ├─ STAGE 2: Methylation extraction
         │       modkit extract (per-read per-CpG) + pileup (aggregate)
         │       stratified by HP:1 / HP:2 tags
         │
         ├─ STAGE 3: ASM detection
         │       per-CpG per-haplotype methylation → |Δ| ≥ 0.4 threshold
         │       cluster adjacent CpGs → dASM loci
         │       AUC test: does HP tag predict methylation at each locus?
         │
         ├─ STAGE 4: Classification + annotation
         │       nearest phased SNP distance → proximal/distal/non-genetic
         │       CpG context (island/shore/shelf/open-sea)
         │       ChromHMM broad state (Vu 2022 universal)
         │       Rosenski bimodal atlas overlap
         │
         ├─ STAGE 5: Validation
         │       Blueprint mono meQTL overlap (external population-level signal)
         │       Imprinting DMR overlap (for non-genetic class)
         │
         └─ STAGE 6: Epiallele complexity (planned)
                 GMM k-order per locus (k=1 uniform / k=2 bimodal / k≥3 multimodal)
                 Shannon entropy of read methylation distribution
                 Cell-type specificity (requires multi-sample)
```

---

## Stage 1: Alignment and phasing

### A00a / W00 — Align ONT BAM

**Input**: HG01258 unaligned BAMs from HPRC S3 (3 runs: 1, 3, 5; ~15 GB each, public)  
**Tool**: minimap2 v2.24, `map-ont` preset, `--MD -y` (preserve MM/ML methylation tags)  
**Output**: `bam/HG01258_wg.bam` — sorted, indexed, merged across runs  
**Status**: Done. BAM header (@PG lines) confirms runs 1, 3, and 5 were all aligned and merged.

Coverage note: The merged BAM is 48 GB. HP-stratified pileup depth is ~2.9x per haplotype (chr22). This is the real working coverage — reads must span a phased het-SNP to be haplotagged, so haplotagged depth is always a fraction of total ONT depth.

### A01a — Get phased VCF

**Input**: 1000 Genomes NYGC high-coverage phased VCF (3202 samples, streamed via bcftools from HTTPS)  
**Tool**: bcftools view → filter to HG01258 heterozygous SNPs only  
**Output**: `HG01258_wg_het_snps.vcf.gz` (~50K SNPs per chromosome)  
**Status**: Done (chr19); WG version uses same logic per chromosome

### A03a / W03 — WhatsHap haplotagging

**Input**: Aligned BAM + het SNP VCF  
**Tool**: WhatsHap `haplotag` — assigns each ONT read to HP:1 or HP:2 based on overlapping phased SNPs  
**Output**: `HG01258_wg_haplotagged.bam` — reads tagged with `HP:i:1` or `HP:i:2`  
**Status**: Done (chr19): 6,247 reads tagged (3,143 HP:1, 3,104 HP:2)

---

## Stage 2: Methylation extraction

### A02a / W01–W02 — modkit extract + pileup

**Input**: Haplotagged BAM + reference FASTA  
**Tool**: `modkit extract` (per-read per-CpG calls with HP tag preserved) + `modkit pileup` (aggregate bedMethyl per strand)  
**Output**:
- `modkit_extract.tsv.gz` — one row per read×CpG: position, probability, HP tag, strand
- `pileup.bed.gz` — aggregate bedMethyl: position, coverage, mod_fraction per strand  
**Status**: Done (chr19); WG not yet run

The HP tag in the extract output is the key link: it lets you compute per-CpG methylation separately for each haplotype.

---

## Stage 3: ASM detection

### W04/W05 (notebook: W05_dasm_loci.ipynb, run as job array per chromosome)

**Input**: HP-stratified modkit pileup per chromosome  
**Logic**:
1. Filter CpGs to ≥3 reads per haplotype
2. Compute `mean_meth_hap1` and `mean_meth_hap2` per CpG
3. Flag CpGs where `|mean_meth_hap1 - mean_meth_hap2| ≥ 0.4`
4. Cluster adjacent flagged CpGs within 200 bp into loci
5. For each locus: AUC test — can HP tag predict read methylation (hap1 vs hap2)?
6. Compute `nearest_snp_dist`: distance from locus midpoint to nearest phased heterozygous SNP

**Output**: `tables/dasm_loci_chr[1-22].tsv` + `phasing_results_chr[1-22].parquet`  
**Key columns**: `locus_id`, `chrom`, `locus_start/end`, `n_cpgs`, `mean_abs_delta`, `n_reads_hap1/2`, `mean_meth_hap1/2`, `auc`, `pval`, `is_phased`, `nearest_snp_dist`, `in_rosenski`  
**Status**: Done (all 22 chr); merged into `tables/dasm_loci_wg.parquet`

**WG results** (304,514 loci total):
- AUC ≥ 0.7 + p < 0.05 ("phased"): 55,315 (18.2%)
- Mean |Δ|: 0.557 (phased), 0.450 (non-phased)
- Mean AUC: 0.82 (phased), 0.62 (non-phased)

### W06 — Merge chromosomes (W06_merge_chromosomes.ipynb)

**Input**: Per-chromosome TSVs and parquets  
**Output**: `tables/dasm_loci_wg.parquet` (304,514 rows × 15 columns)  
**Status**: Done

---

## Stage 4: Classification and annotation

### W07 — Classify and annotate whole genome (W07_classify_annotate_wg.ipynb)

**Input**: `dasm_loci_wg.parquet`  
**Steps**:

**4a. Classification by SNP distance**

| Class | Criterion | n | % |
|---|---|---|---|
| proximal-genotype | is_phased + nearest_snp_dist ≤ 200 bp | 18,498 | 6.1% |
| distal-genotype | is_phased + 200 bp < dist ≤ 5 kb | 32,326 | 10.6% |
| very-distal-genotype | is_phased + 5 kb < dist ≤ 50 kb | 4,490 | 1.5% |
| non-genetic | not is_phased | 249,199 | 81.8% |

Long reads detect 3× more genotype-driven ASM than short reads could (36,816 long-read-only vs 18,498 short-read-detectable).

**4b. CpG island context** (ref: `cpgIslandExt.hg38.bed`, 31,448 islands)  
- island (direct overlap) / shore (≤2 kb) / shelf (2–4 kb) / open_sea (>4 kb)  
- Key finding: genotype-driven ASM is depleted from CpG islands (OR=0.71, p=5×10⁻¹⁶) relative to non-genetic. Most ASM (~90%) occurs in open sea.

**4c. ChromHMM broad state** (ref: Vu 2022 universal states, `vu.states.bed.gz`)  
- Proximal-genotype loci depleted at Active TSS (OR=0.60, FDR=3×10⁻⁴) and exons (OR=0.75) — consistent with purifying selection at constrained elements.

**4d. Rosenski atlas overlap** (ref: Rosenski 2025 bimodal CpG atlas)  
- ~95% of all loci overlap Rosenski regardless of class — the atlas captures both genetic and non-genetic ASM with equal sensitivity.

**Output**: `tables/dasm_loci_wg_annotated.parquet` (adds `classification`, `cpg_context`, `chromhmm_broad`)  
**Status**: Done

---

## Stage 5: Validation

### W08 — Blueprint meQTL overlap (W08_meqtl_overlap.ipynb)

**Scientific question**: Do dASM loci classified as genotype-driven overlap with known population-level meQTL CpGs from an independent bulk study (Blueprint / Bonder 2017)?

This is the key external validation: if proximal-genotype ASM loci are enriched for Blueprint meQTL signal, it confirms the phasing correctly identifies cis-genetic methylation effects. The distal classes should show intermediate enrichment, quantifying the added meQTL discovery value of long reads.

**Input**:
- `dasm_loci_wg_annotated.parquet`
- Blueprint mono meQTL summary stats (`data/raw/meqtl_mono.fdr05.tsv.gz`, 17 GB, all tested pairs)
- 450K probe positions hg38 (`data/450k_probes_hg38.bed`, 473,699 probes — generated from R `IlluminaHumanMethylation450kanno.ilmn12.hg19` + pyliftover)

**Preprocessing** (W08a_extract_meqtl_cpgs.sh, job 13608677, running):
- Stream 17 GB file; extract CpG probe IDs where corrected pval < 0.05
- Output: `data/meqtl_mono_sig_cpgs.txt` (unique significant CpG probe IDs)

**Analysis**:
1. 450K coverage: what fraction of dASM loci sit on a 450K-assayable CpG? (baseline: ~4% overall, flat across classes — expected given 90% open-sea loci)
2. meQTL enrichment: among 450K-covered loci, is there enrichment for meQTL signal by class? Fisher's exact vs non-genetic reference.
3. Effect concordance: does |ASM delta| correlate with |meQTL beta| at overlapping loci?

**Status**: 450K coverage done; meQTL enrichment awaiting W08a job completion (~2–3 hr)

### W09 (planned) — Imprinting overlap

**Scientific question**: What is the non-genetic class? Likely a mixture of: (a) imprinted loci, (b) stochastic/noise ASM, (c) loci where heterozygosity is present but read length didn't reach a phased SNP.

**Input**: Zink et al. 2018 (Nat Genetics) Supplementary Table 4 — ~4,000 Icelandic ASM loci with imprinting flags, already hg38. Download Supplementary Table 4 from the paper.

**Analysis**: Overlap non-genetic dASM loci with known imprinted DMRs. Expected: non-genetic loci enriched at imprinted regions; proximal-genotype loci depleted.

---

## Stage 6: Epiallele complexity (planned, chr19 framework exists)

Beyond mean haplotype differences, individual reads carry multi-CpG methylation *patterns*. A locus could be:
- **k=1** (uniform): all reads ~50% methylated — likely artifact or cell-type mixture, not true ASM
- **k=2** (bimodal): reads cluster into fully methylated vs. fully unmethylated — canonical ASM, consistent with Rosenski bimodal CpG framework
- **k≥3** (multimodal): discrete but more than two states — rare, could indicate multi-allelic regulation or cell-type admixture

### A07 / WG equivalent — GMM complexity per locus

**Input**: `epiallele_store` (read × CpG matrix per locus, HP labels)  
**Tool**: Gaussian Mixture Model (k=1,2,3), BIC model selection  
**Output**: Per-locus `k_optimal`, mixture weights, per-read assignment  
**Status**: Chr19 framework implemented (A07_epiallele_complexity.ipynb); WG run pending

### A09 / WG equivalent — Entropy analysis

Shannon entropy H = -p·log₂(p) - (1-p)·log₂(1-p), where p = fraction of methylated reads at the locus.

- H≈1 (balanced 50:50) → classic heterozygous ASM
- H≈0 (monoallelic) → imprinting-like (one haplotype fully silenced)
- Intermediate H with no genetic link → stochastic/cell-type-driven

Cross-sample entropy comparison (requires multiple individuals or tissues) can identify tissue-specific vs. constitutive epiallele structure.

**Status**: Chr19 framework implemented (A09_entropy_analysis.ipynb); WG run pending

---

## Key output files

| File | Description | Status |
|---|---|---|
| `tables/dasm_loci_wg.parquet` | 304,514 dASM loci, 15 columns | Done |
| `tables/dasm_loci_wg_annotated.parquet` | + classification, CpG context, ChromHMM | Done |
| `tables/W08_meqtl_overlap_summary.csv` | meQTL enrichment by class | Ready — re-run W08 |
| `figures/W07_*.pdf` | Classification and annotation figures | Done |
| `figures/W08_*.pdf` | meQTL overlap figures | Partial |

---

## Naming conventions

- `A` series: chr19 proof-of-concept (low coverage, framework validation)
- `W` series: whole-genome production analysis
- `W0N_*.sh`: qsub job scripts (STAGE 1–3, computationally heavy)
- `WNN_*.ipynb`: analysis notebooks (STAGE 4–6, interactive)
- `WNNa_*.sh`: preprocessing helper jobs (e.g., W08a extracts meQTL CpGs)
