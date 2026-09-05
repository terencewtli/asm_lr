# Epiallele Entropy Analysis — W12 Results

*Notebook: W12_epiallele_entropy.executed.local.ipynb | Run: 2026-06-12 (local, single sample)*

## Dataset

304,514 methylation loci genome-wide (all 22 autosomes), single sample pilot. Loci classified by proximity to nearest heterozygous SNP:

| Class | n | % |
|---|---|---|
| non-genetic | 249,199 | 82% |
| distal-genotype | 32,326 | 11% |
| proximal-genotype | 18,498 | 6% |
| very-distal-genotype | 4,490 | 1% |

## Metrics

Five metrics computed per locus from the read-level methylation matrix (reads × CpGs):

- **per_cpg_entropy**: Binary entropy of average methylation at each CpG, averaged across CpGs (aggregate; captures between-haplotype bimodality)
- **within_hp_entropy**: Same entropy computed within HP1 and HP2 separately, then averaged (captures within-haplotype read variability — determinism)
- **hp_sep**: |mean_HP1 − mean_HP2| — between-haplotype methylation difference (dASM effect size)
- **within_hp_var**: Variance of per-read methylation means within each haplotype, averaged across HP1/HP2
- **bimodality_coef**: (skewness² + 1) / (kurtosis + 3) — standard bimodality index applied to per-read means

## Main Results (MWU, BH FDR)

| Class | hp_sep | within_hp_entropy | within_hp_var | per_cpg_entropy | bimodality_coef |
|---|---|---|---|---|---|
| non-genetic (ref) | 0.048 | 0.265 | 0.0188 | 0.284 | 0.497 |
| proximal-genotype | **0.161** ↑ p→0 | 0.254 ↓ p<1e-54 | 0.0162 ↓ p<1e-86 | 0.292 ↑ p<1e-26 | 0.512 ↑ p<1e-60 |
| distal-genotype | **0.162** ↑ p→0 | 0.248 ↓ p<1e-243 | 0.0170 ↓ p<1e-59 | 0.282 (ns) | 0.512 ↑ p<1e-114 |
| very-distal-genotype | **0.168** ↑ p→0 | 0.246 ↓ p<1e-41 | 0.0174 ↓ p<1e-8 | 0.280 (ns) | 0.509 ↑ p<1e-10 |

## Key Biological Findings

### 1. Genotype-driven ASM is deterministic

`within_hp_entropy` and `within_hp_var` are significantly **lower** for all genotype classes. Reads on HP1 look alike; reads on HP2 look alike. In contrast, non-genetic methylation variation shows higher within-haplotype read-to-read variability (stochastic epivariation). This determinism vs stochasticity axis is the most biologically interpretable finding: local genotype constrains the methylation state, while non-genetic methylation is a more probabilistic process.

### 2. hp_sep does not decay with distance

Between-haplotype separation stays high across all distance classes (0.160–0.168), and does not decrease as distance from the nearest het SNP increases. Combined with the fact that `per_cpg_entropy` loses significance at distal/very-distal classes while `within_hp_entropy` remains significant, this suggests:
- The aggregate bimodality (captured by `per_cpg_entropy`) weakens at distance — fewer reads show a clean two-group split across all CpGs simultaneously
- But the *mean* methylation difference between haplotypes persists

### 3. The distance gradient in within_hp_entropy is counterintuitive

Within-haplotype entropy decreases monotonically with distance from nearest SNP:
- proximal: 0.254 → distal: 0.248 → very-distal: 0.246

You might expect the genotype effect to weaken with distance and within-haplotype reads to become noisier. Instead they become more homogeneous. Working hypotheses:
- Very-distal-genotype loci may be enriched for large epigenetic domains (imprinting control regions, TAD-bounded allelic domains, CTCF loops) where methylation state is propagated coherently far from the causal SNP
- Selection bias: loci that pass the dasm_loci threshold at long distance may be the most extreme cases (very high hp_sep + very clean within-haplotype reads), not typical distal effects

**Follow-up**: intersect very-distal-genotype loci with known imprinting control regions, CTCF sites, and TAD boundaries.

## Novelty Assessment

This is an interesting decomposition. Most dASM papers report only the between-haplotype mean difference (analogous to `hp_sep`). Separating:

- **hp_sep**: is there dASM? (between-haplotype)
- **within_hp_entropy / within_hp_var**: is the dASM deterministic or stochastic? (within-haplotype)
- **asymmetric HP variance** (not yet computed): is one haplotype more variable than the other?

is a more complete picture of epiallele structure. The finding that genetic proximity → determinism is testable across the 30-sample cohort and would be a strong component of the paper.

## Planned Follow-ups

### Near-term
1. **Run W12 across all 30 samples** — replicate distributions; identify loci with consistent vs individual-specific ASM
2. **Intersect very-distal-genotype with known biology** — ICRs, CTCF peaks, TADs, imprinted genes
3. **Asymmetric HP variance metric** — compute `|var(HP1_read_means) - var(HP2_read_means)|` per locus; expected to be highest at imprinted loci (one allele deterministic, other variable)

### Population-level questions
4. **Population variation in epiallele structure** — across 30 individuals, does the within_hp_entropy at ASM loci correlate with local heterozygosity or haplotype diversity?
5. **Shared vs private ASM** — loci where hp_sep is high in most individuals (constitutive ASM) vs rare individuals (private epivariation)

### Asymmetric HP variance (proposed)
The current `within_hp_var` averages `var(HP1_read_means)` and `var(HP2_read_means)`. This discards the asymmetry between haplotypes. Computing `|var(HP1) - var(HP2)|` would specifically detect loci where one haplotype is fixed and the other is variable — the signature of imprinting (one allele locked by parent-of-origin imprint, other allele displaying somatic epivariation). This metric should be highest at known imprinting control regions, providing a built-in positive control.
