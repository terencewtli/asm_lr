#!/usr/bin/env python3
"""
The one merge step in the pipeline: per-chromosome Fisher ASM output ->
whole-genome per sample -> cohort-wide. Everything upstream (W02, W03, W04)
stays per-chromosome; this is the only place results get concatenated.
"""

import pandas as pd

PROJDIR = '/u/project/cluo/terencew/claude/project_ideas/asm_lr'
SAMPLES = open(f'{PROJDIR}/txt/samples/final_samples.txt').read().split()
CHROMS = [f'chr{i}' for i in range(1, 23)]
OUTDIR = f'{PROJDIR}/tables/dasm_fisher'

cohort_loci = []
cohort_phasing = []

for sample in SAMPLES:
    sample_dir = f'{OUTDIR}/{sample}'

    loci_parts = []
    for chrom in CHROMS:
        f = f'{sample_dir}/dasm_loci_{chrom}_fisher.tsv'
        try:
            loci_parts.append(pd.read_csv(f, sep='\t'))
        except FileNotFoundError:
            print(f'  missing: {f}')
    phasing_parts = []
    for chrom in CHROMS:
        f = f'{sample_dir}/phasing_results_{chrom}_fisher.parquet'
        try:
            phasing_parts.append(pd.read_parquet(f))
        except FileNotFoundError:
            print(f'  missing: {f}')

    if not loci_parts and not phasing_parts:
        print(f'{sample}: no chromosome output found, skipping')
        continue

    loci_df = pd.concat(loci_parts, ignore_index=True) if loci_parts else pd.DataFrame()
    phasing_df = pd.concat(phasing_parts, ignore_index=True) if phasing_parts else pd.DataFrame()

    loci_df.to_csv(f'{OUTDIR}/{sample}_wg_dasm_loci_fisher.tsv', sep='\t', index=False)
    phasing_df.to_parquet(f'{OUTDIR}/{sample}_wg_phasing_results_fisher.parquet', index=False)
    print(f'{sample}: {len(loci_df):,} loci across {len(loci_parts)}/22 chroms, '
          f'{len(phasing_df):,} phasing rows across {len(phasing_parts)}/22 chroms')

    if len(loci_df):
        loci_df.insert(0, 'sample', sample)
        cohort_loci.append(loci_df)
    if len(phasing_df):
        phasing_df.insert(0, 'sample', sample)
        cohort_phasing.append(phasing_df)

if cohort_loci:
    pd.concat(cohort_loci, ignore_index=True).to_parquet(
        f'{OUTDIR}/cohort_wg_dasm_loci_fisher.parquet', index=False)
if cohort_phasing:
    pd.concat(cohort_phasing, ignore_index=True).to_parquet(
        f'{OUTDIR}/cohort_wg_phasing_results_fisher.parquet', index=False)

print(f'\ncohort: {len(cohort_loci)}/{len(SAMPLES)} samples with loci, '
      f'{len(cohort_phasing)}/{len(SAMPLES)} samples with phasing results')
