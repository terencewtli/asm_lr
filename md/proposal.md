⸻
Background
Allele-specific methylation (ASM) is traditionally defined using heterozygous SNPs to phase sequencing reads and compare methylation levels between haplotypes. This framework captures genotype-driven cis effects on DNA methylation, but it reduces methylation to a single summary statistic per allele.
In parallel, epigenomic data (especially single-molecule or long-read sequencing) reveals richer structure: individual reads often contain multi-CpG methylation patterns (“epialleles”) that reflect coordinated methylation states along a molecule. These patterns can arise from:
* genetic variation (cis effects),
* cell-state heterogeneity,
* or stochastic epigenetic regulation.
Current approaches treat these two views separately:
* ASM studies focus on haplotype-level methylation differences
* epiallele studies focus on within-read methylation structure, usually without genotype context
As a result, there is no systematic framework linking:
genotype-defined ASM ↔ single-molecule methylation structure (epialleles)
Long-read sequencing (e.g., nanopore) makes this link possible by jointly observing:
* phased SNPs (haplotypes)
* dense CpG methylation patterns along the same reads
⸻
Core Idea
The central question is:
To what extent is single-molecule methylation structure (epialleles) explained by genetic variation (haplotypes)?
This reframes ASM as a special case of a broader phenomenon:
methylation pattern structure on individual DNA molecules
⸻
Proposed Project
Goal
Quantitatively decompose DNA methylation structure into:
* genetic (haplotype-driven) components
* non-genetic (epigenetic or stochastic) components
⸻
Approach
Using long-read sequencing data (e.g., nanopore):
1. Identify heterozygous SNPs to phase reads into haplotypes
2. Extract CpG methylation patterns per read (epialleles)
3. Cluster reads into epiallele states based on methylation pattern similarity
4. For each genomic locus:
    * compute classical ASM (haplotype methylation differences)
    * compute epiallele structure (entropy / clustering)
    * quantify alignment between them (e.g., mutual information or predictive accuracy)
⸻
Outputs
A genome-wide classification of loci into:
* genotype-driven methylation structure (strong haplotype–epiallele alignment)
* mixed control loci (partial alignment)
* genotype-independent epialleles (structure without genetic explanation)
Plus illustrative case studies of each regime.
⸻
Feasibility
* Single individual or small cohort (1–3 genomes) is sufficient for the core analysis
* Population-scale data (e.g., HPRC) is optional and mainly useful for generalization
* A ~5-month project is feasible if scoped around:
    * careful locus filtering
    * simple, robust metrics (entropy, clustering, MI)
⸻
Expected Contribution
A conceptual and quantitative framework that:
* connects ASM and epiallele structure at single-molecule resolution
* reinterprets “methylation heterogeneity” as a mixture of genetic and non-genetic signals
* provides a new lens for long-read epigenomics beyond average methylation levels
