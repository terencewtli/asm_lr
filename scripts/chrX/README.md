# chrX pipeline — status, design decisions, and what's still missing

**Written 2026-09-07. Scripts here are ready to review but NOT yet run at any scale** — see
"What's missing before this can execute" below.

## Why chrX, and why it's not just "port the autosome pipeline"

X-chromosome inactivation (XCI) is a natural special case of ASM: female cells randomly
silence one parental X per cell and densely methylate its CpG islands. Framed the same way as
imprinting — a parent-of-origin-driven methylation asymmetry, independent of any local genetic
variant — except (a) it's *usually* not systematically skewed toward one specific parent at
the population level (unlike imprinting, which always picks the same parent), and (b) it's
**cell-autonomous, not individual-autonomous**: each cell independently chooses which X to
silence. That has a real consequence for haplotype-phased bulk data:

- **Random (unskewed) XCI, phased by parental origin**: pooling reads from many cells, each
  haplotype (paternal-origin X, maternal-origin X) is a roughly 50/50 mixture of
  active-in-this-cell and inactive-in-this-cell molecules. The *average* methylation per
  haplotype comes out similar for both — **no population-level ASM signal**, even though every
  single cell has a real, strong, binary XCI pattern. The signal shows up as **within-haplotype,
  between-read bimodality** instead (some HP1 reads high, some HP1 reads low, because they came
  from different cells) — not as a between-haplotype mean difference.
- **Skewed XCI** (one parental X preferentially inactivated across most cells — happens with
  age, selection, or by chance) *does* produce a real between-haplotype mean difference, i.e.
  genuine population-level ASM, in exactly the donors where it's skewed.
- **Males have no chrX ASM at all outside PAR** (hemizygous, one copy, maternal-origin only) —
  a clean, free negative control that autosomal loci don't have.
- This also makes chrX a natural **positive control for the design-effect/read-level question**
  from `md/20260907_phasing_qc_review.md` §6: XCI is a genuine, strong, real source of
  within-molecule co-methylation (every CpG on a read reflects that one cell's XCI state), so
  it's a good place to check whether `readlevel.test_region_reads` correctly detects real
  bimodal structure that the pooled beta-binomial test would be maximally vulnerable to.

## Design decisions

1. **Reference: chrX only, not chrX+chrY.** Extracted from
   `/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.fa` (which has both) via
   `samtools faidx ... chrX` → `reference/chrX_ref/GRCh38.chrX.fa`. With no chrY in the
   reference, PAR1/PAR2 reads have nowhere else to map and simply place uniquely on chrX — no
   PAR-masking logic needed, which would otherwise be required (chrY normally needs its PAR
   copy masked with N's to avoid multi-mapping ambiguity with chrX's PAR).
2. **Female donors only for haplotagging** (`txt/samples/female_donors.txt`, 8 of 18 donors:
   NA18508, HG00146, HG00344, HG00253, NA19776, NA21110, NA21144, HG03784). Non-PAR chrX in
   males is hemizygous — there's no haplotype pair to phase or compare. Male chrX is still
   useful (bulk coverage as the negative control above), just doesn't need `X01_haplotag_chrX.sh`.
3. **The `whatshap haplotag` invocation itself needs no chrX-specific change** — no `--ploidy`
   flag, nothing PAR-specific — because chrX is genuinely diploid in females outside PAR, same
   as any autosome. The only thing that changes is the **input VCF**: chrX isn't in the 1000G
   panel VCF set the autosome pipeline uses (or if it technically is, it wasn't the source this
   project standardized on) — the project owner has a separate, already-generated
   **TOPMed-imputed, statistically phased chrX VCF** for exactly this purpose.
4. **X00a aligns the full raw uBAM against the chrX-only reference** (same
   `csv/meta/ubam_manifest.csv`, same 72 tasks) rather than trying to extract chrX reads from
   the existing whole-genome BAMs — those were aligned against `GRCh38.autosome.fa`, which has
   no sex chromosomes at all, so there are no chrX-mapped reads to extract from them. This is a
   genuine realignment from raw data, not a cheap slice of existing output.

## What's missing before this can execute

- **`CHRX_VCF_SOURCE` in `X01_haplotag_chrX.sh` is an empty placeholder.** Point it at the
  TOPMed chrX phased VCF (single file or per-sample) before running X01. X00a/X00b need no such
  input and are ready to submit as-is.
- **Not yet run at any scale.** Aligning against a small reference is faster per-uBAM than the
  full autosome alignment, but this is still a genuine realignment of every raw uBAM for every
  sample being run — a real compute commitment, not free. Recommend piloting on 1-2 female
  donors (e.g. the two among the confirmed-sane six: HG00146, NA21110) before committing to all
  8 female donors or the full 18-donor set (10 of which are male and only get bulk-coverage
  value, not haplotype comparison).
