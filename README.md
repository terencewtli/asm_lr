# asm_lr — Allele-Specific Methylation from Long Reads

## Overview

Long-read (ONT) sequencing jointly observes phased heterozygous SNPs and dense
CpG methylation on the same molecule. This lets us:

1. Assign reads to haplotypes via phased SNPs
2. Measure methylation per CpG per haplotype
3. Test whether the two haplotypes differ at each locus (allele-specific
   methylation, ASM)
4. Classify ASM by how far its driving SNP is: proximal (<200 bp, detectable
   by short reads) vs. distal (200 bp–50 kb, long-read-only) vs. non-genetic
   (no phased SNP linkage at any distance)
5. Characterize single-molecule methylation *patterns* (epialleles) beyond
   mean differences

The central question: **to what extent does single-molecule methylation
structure arise from genetic variation?** This reframes ASM as a special case
of a broader phenomenon — methylation pattern structure on individual DNA
molecules — rather than a haplotype-difference summary statistic. See
[`md/updated_proposal.md`](md/updated_proposal.md) for the full scientific framing.

Short reads can only detect SNP-CpG associations within ~200 bp (one read).
Long reads span 5–20 kb, phasing methylation over distances invisible to
Illumina.

## Cohort

**18 donors** (`txt/samples/final_samples.txt`) — the whole-genome
replication cohort — are the working cohort for this project. They are not a
partial subset of a larger planned cohort; the scope was deliberately
narrowed from an original 30-donor design to these 18 for a documented
scientific reason (not just convenience):

- All 18 are basecalled with **Dorado 0.6.0** (one sample 0.5.3) on
  **R10.4.1** chemistry, with joint `5mCG_5hmCG` modification calling.
- A separate batch of 12 donors (`txt/samples/excluded_samples.txt`,
  including the original HG01258 chr19 pilot) uses **Guppy 6.5.7** on
  **R9.4.1** chemistry with 5mC-only calling. Different pore chemistry and a
  different modification-calling model both bear directly on ASM calling, so
  pooling the two batches would let basecaller identity masquerade as
  biological signal. See [`md/basecaller_cohort_split.md`](md/basecaller_cohort_split.md)
  for the full rationale.
- Full provenance of how the 30-donor pool was originally assembled and why
  it narrowed to 18 is reconstructed in
  [`md/20260904_cohort_provenance.md`](md/20260904_cohort_provenance.md),
  including ancestry composition (AFR is thin at n=2 — no ancestry-stratified
  claim is supportable at that count) and a documented path to a balanced
  38-donor expansion later, if needed. Scaling is strictly additive: nothing
  computed on the 18 is invalidated by adding donors under the same
  chemistry/model.
- Open item: `NA18959` (one of the 18) has no confirmed raw-data directory
  under HPRC `working/` as of the last check — resolve before treating the
  cohort as a firm n=18.

## Pipeline architecture

Whole-genome production pipeline, orchestrated via qsub (SGE array jobs) —
`scripts/wg/`, prefixed `W`/`QC`. (An earlier Snakemake-based orchestration
was tried and abandoned; the pipeline is qsub-only now.)

| Stage | Script | Does | Status |
|---|---|---|---|
| Coverage QC | `QC01_mosdepth_coverage_wg.sh` | mosdepth coverage per sample | — |
| W00a | `W00a_align_ubam.sh` | minimap2 align ONT uBAM (`map-ont`, `--MD -y` to preserve MM/ML tags) | Done, all 18 |
| W00b | `W00b_merge_wg.sh` | merge per-run BAMs into one whole-genome BAM per sample | Done, all 18 |
| W01 | `W01_phase_per_chrom.sh` | **WhatsHap `haplotag`** per chromosome against het SNPs pulled from the **1000 Genomes 3202 phased panel** (read-based haplotagging, not statistical SHAPEIT4 phasing) | Done, all 18 — see caveat below |
| W01g | `W01g_merge_haplotagged_wg.sh` | merge per-chrom haplotagged BAMs into one per sample | Done, all 18 |
| W02 | `W02_modkit_extract_per_chrom.sh` | `modkit extract` — per-read per-CpG calls, HP tag preserved | Mostly done |
| W03 | `W03_modkit_pileup_per_chrom.sh` | `modkit pileup` — HP-stratified aggregate bedMethyl | 392/396 sample×chrom cells done |
| W04 | `W04_dasm_loci_per_chrom.sh` + `W04_dasm_loci.py` | Fisher-exact + BH-FDR dASM calling per (sample, chrom), then Mann-Whitney U / AUC phasing test | **Pilot in progress** — see below |
| W05 | `W05_merge_dasm_wg.sh` / `.py` | merge per-chrom dASM results genome-wide | Pending W04 |

`scripts/download/` (`D01`–`D04`) fetches HPRC ONT uBAMs, 1000G phased VCFs,
and HPRC assemblies (the latter needed for `dipcall` assembly-backed
phasing, see caveat below). `scripts/genozip/` compresses/decompresses BAMs
for archival. `scripts/NA21093/` keeps only `W12_epiallele_entropy.{sh,py}` —
the epiallele-complexity analysis behind [`md/epiallele.results.md`](md/epiallele.results.md),
for which no equivalent exists yet in `scripts/wg/`; the rest of that earlier
single-donor development pass was dropped as superseded by `scripts/wg/`.

Chr19-scope proof-of-concept results for the original HG01258 pilot sample
(SHAPEIT4-phased, Guppy/R9.4.1 — now in the excluded batch) are documented in
[`md/pipeline.md`](md/pipeline.md), including the full downstream
classification/annotation/validation pipeline (Stages 3–6: dASM detection,
CpG-context and ChromHMM annotation, meQTL/imprinting validation, planned
epiallele-complexity analysis). That pipeline design carries forward to the
18-donor whole-genome run above.

## Open items

- **Haplotagging rate needs systematic validation across the 18-donor
  cohort.** W01 already uses WhatsHap read-based phasing (not the SHAPEIT4
  approach criticized for the old HG01258 pilot in
  [`md/20260903_qc_review.md`](md/20260903_qc_review.md) — that review is
  about the older single-donor pilot and a different phasing approach; don't
  conflate it with the current W01 pipeline, though its QC methodology notes
  and priority list remain broadly relevant). A spot check on one sample/chr
  (NA18508, chr20) showed ~38% of reads receiving an HP tag — better than the
  old pilot's ~18%, but still under the ≥60% QC gate this project targets.
  Worth checking systematically before trusting phased/haplotype-stratified
  results.
- **Fisher-based dASM calling (W04) has never been run at scale.** A small
  timing pilot (a handful of tasks, not the full 396-task array) was just
  submitted, since its resource requests (`h_data`/`h_rt`) were unvalidated
  placeholders and the one earlier pilot attempt died during the upstream
  modkit pileup step before ever reaching the Fisher call. Full-scale timing
  and submission are still pending.

## Other documentation

`md/` is kept intentionally lean — methodology and results documentation only, no
working notes or chat-style critique threads:

- [`progress.md`](md/progress.md) — dated status log of pipeline results
- [`epiallele.results.md`](md/epiallele.results.md) — chr19 epiallele-complexity results (HG01258 pilot)
- [`pipeline.md`](md/pipeline.md) — full pipeline architecture and POC results
- [`basecaller_cohort_split.md`](md/basecaller_cohort_split.md) — why the cohort is 18, not 30
- [`20260904_cohort_provenance.md`](md/20260904_cohort_provenance.md) — cohort selection methodology
- [`20260903_qc_review.md`](md/20260903_qc_review.md) — QC findings and phasing-validity review
- [`updated_proposal.md`](md/updated_proposal.md) — scientific framing and motivation
