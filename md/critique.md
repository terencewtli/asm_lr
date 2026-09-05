What the pipeline does

Single individual (HG01258), ONT whole-genome → haplotagged BAM → haplotype-stratified methylation pileups → dASM loci
called at |Δ| ≥ 0.3, AUC ≥ 0.70, p < 0.05 → classified by distance to nearest phased het SNP → annotated with CpG
context, ChromHMM, Rosenski bimodal atlas, Blueprint meQTLs.

---
Validity concerns

The biggest problem: slide 5 claims meQTL validation hasn't been run. W08 executed with HAS_MEQTL=False —
W08a_extract_meqtl_cpgs.sh hasn't been run, so meqtl_mono_sig_cpgs.txt is missing. The "6× enriched for Blueprint
meQTLs" claim on slide 5 and the summary slide are unsupported by any computed number. This needs to be fixed before
sharing the slides.

The "3× more" claim is a proxy, not an empirical test. The comparison is 55,315 phased vs 18,498 proximal (<200 bp). The
200 bp cutoff approximates short-read detectability, but you never actually ran short-read phasing on the same sample.
The correct test would be to run bismark + a short-read phasing tool (MethHap, biscuit, or phASER) on HG01258 and
compare yields directly. Without that, the multiplier is an assumption baked into the threshold.

"Phased" and "near a SNP" are not independent. Loci are classified non-genetic because they're far from phased het SNPs.
But haplotagging (WhatsHap) itself requires SNPs — loci distant from SNPs may just be in low-heterozygosity regions or
poorly phased blocks, not genuinely genotype-independent. The 81.8% non-genetic rate likely conflates true non-genetic
regulation with phasing dropout.

Systematic HP1/HP2 read imbalance. In chr1, hap1 consistently has fewer reads than hap2 (e.g., 8–9 vs 13). This is
visible in the W05 output across many loci. A systematic phasing asymmetry could bias AUC calling. Worth checking
whether this is uniform across chromosomes or localized.

The ultra-distal class has 1 locus. This is an artifact, not a finding — probably a VCF edge case or a phasing block
that runs unusually far. Should be collapsed into very-distal or excluded.

The epiallele store is built but never used. The proposal's core idea — linking genotype ASM to single-molecule
methylation structure — is never executed. W05 builds a full per-read methylation matrix for every locus, but W07/W08
only use per-read mean methylation. This is a significant gap between the proposal and current results.

Rosenski overlap at 94.5% is validating but also redundant. It confirms you're finding real variable methylation, but it
doesn't validate the phasing or distance classification specifically. More informative would be the meQTL enrichment
(when it's actually run).

---
What's solid

- Effect sizes (median |Δ| = 0.557, AUC ≈ 0.82) are genuinely consistent across all distance classes — this is the
cleanest result in the deck. No attenuation with distance is a real and interesting finding.
- ChromHMM and CpG context distributions look coherent — distal loci are in the same regulatory contexts as proximal,
which supports the functional equivalence claim.
- The pipeline itself is clean — vectorized epiallele store in W05, dual threshold phasing, proper Fisher's exact with
FDR correction.

---
Next steps

Immediate (before sharing slides):
1. Run W08a_extract_meqtl_cpgs.sh and execute W08 with meQTL data — slide 5 is currently unsupported
2. Clarify the "3× more" claim as "3× by 200 bp proxy" not empirical comparison

Short-term:
3. Multi-donor replication — even 3 additional HPRC ONT samples. The multiplier is meaningless from N=1
4. Investigate HP1/HP2 read imbalance across chromosomes — could be a WhatsHap phasing artifact
5. Stratify all results by n_cpgs ≥ 2 — single-CpG loci (the median) are noisy; the multi-CpG subset is the stronger
story

Longer-term (connecting to the proposal):
6. Actually use the epiallele matrices — compute per-locus methylation entropy or epiallele clustering and test whether
genotype-driven loci have lower entropy (more bimodal) than non-genetic loci. This is the scientifically novel part of
the proposal that's currently missing
7. For distal loci, distinguish phasing-tag SNP distance from causal SNP distance — comparing against meQTL SNP
positions would help here
