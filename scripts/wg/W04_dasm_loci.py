#!/usr/bin/env python3
"""
dASM locus calling -- Fisher's exact replacement for the |delta| >= 0.3
raw-delta threshold in W04_dasm_loci.py.

Stage 1 (CpG level) used to be an effect-size-only cutoff: flag any CpG
where |meth_hp1 - meth_hp2| >= 0.3, regardless of depth. A CpG at 5x/5x
and one at 200x/200x got the same treatment, so low-depth noise and
genuine signal were indistinguishable. This replaces it with a per-CpG
Fisher's exact test on the [n_mod, n_canonical] x [hp1, hp2] contingency
table (from modkit pileup --partition-tag HP bedMethyl counts) + BH-FDR
correction across all tested CpGs, with a small effect-size floor
(DELTA_FLOOR) kept only to drop FDR-significant-but-biologically-trivial
deltas at very high depth.

Stage 2 (locus level, Mann-Whitney U + AUC on per-read methylation) is
unchanged from W04_dasm_loci.py -- it was already a real statistical test.

Two subcommands, run in sequence by the qsub wrapper:
  call-loci   pileup HP1/HP2 bedMethyl -> Fisher+FDR -> cluster -> loci BED
  phase-loci  modkit-extract (--include-bed loci BED) + haplotagged BAM
              -> per-locus Mann-Whitney U / AUC -> final tables

Splitting it this way means the expensive per-read extract only has to
cover the (small) set of candidate loci, not the whole genome -- see
--include-bed in `modkit extract --help`.
"""

import argparse
import json
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pysam
from scipy.stats import fisher_exact, mannwhitneyu, false_discovery_control
from sklearn.metrics import roc_auc_score

# ── thresholds ────────────────────────────────────────────────────────────
MIN_COV_PER_HAP = 5
FDR_THRESH      = 0.05  # BH-FDR on the per-CpG Fisher p-values
DELTA_FLOOR     = 0.10  # effect-size floor kept alongside FDR significance
GAP_MERGE_BP    = 500
FLANK           = 500
MIN_HAP_READS   = 3
AUC_THRESH      = 0.70
PVAL_THRESH     = 0.05  # locus-level Mann-Whitney threshold (unchanged)

BEDMETHYL_COLS = [
    'chrom', 'start', 'end', 'mod_code', 'depth', 'strand',
    'start2', 'end2', 'color', 'depth2', 'pct_mod',
    'n_mod', 'n_canonical', 'n_other_mod', 'n_delete',
    'n_fail', 'n_diff', 'n_no_call',
]


def load_pileup(path):
    df = pd.read_csv(
        path, sep='\t', header=None, names=BEDMETHYL_COLS,
        usecols=['chrom', 'start', 'end', 'mod_code', 'n_mod', 'n_canonical'],
    )
    return df[df['mod_code'] == 'm'].drop(columns='mod_code')


def _fisher_row(args):
    n_mod1, n_can1, n_mod2, n_can2 = args
    _, p = fisher_exact([[n_mod1, n_can1], [n_mod2, n_can2]])
    return p


def call_cpgs(merged, threads):
    tables = list(zip(
        merged['n_mod_hp1'], merged['n_canonical_hp1'],
        merged['n_mod_hp2'], merged['n_canonical_hp2'],
    ))
    print(f'  running Fisher\'s exact on {len(tables):,} CpGs ({threads} procs)...')
    if threads > 1 and len(tables) > 10_000:
        with Pool(threads) as pool:
            pvals = pool.map(_fisher_row, tables, chunksize=2000)
    else:
        pvals = [_fisher_row(t) for t in tables]
    merged['pval'] = pvals
    merged['qval'] = false_discovery_control(merged['pval'], method='bh')
    return merged


def cluster_cpgs(dasm_cpgs, chrom, gap_bp):
    loci = []
    cur_start = cur_end = None
    cur_cpgs = []
    for _, row in dasm_cpgs.iterrows():
        if cur_start is None or row['start'] - cur_end > gap_bp:
            if cur_start is not None:
                loci.append(_locus_record(chrom, cur_start, cur_end, cur_cpgs))
            cur_start, cur_end, cur_cpgs = int(row['start']), int(row['end']), [row]
        else:
            cur_end = int(row['end'])
            cur_cpgs.append(row)
    if cur_start is not None:
        loci.append(_locus_record(chrom, cur_start, cur_end, cur_cpgs))
    return loci


def _locus_record(chrom, start, end, cpgs):
    directions = [1 if c['meth_hp1'] > c['meth_hp2'] else -1 for c in cpgs]
    frac_concordant = abs(sum(directions)) / len(directions)
    return {
        'chrom':            chrom,
        'start':            start,
        'end':              end,
        'locus_id':         f'{chrom}:{start}-{end}',
        'n_cpgs':           len(cpgs),
        'mean_abs_delta':   float(np.mean([abs(c['delta']) for c in cpgs])),
        'mean_meth_hp1':    float(np.mean([c['meth_hp1']  for c in cpgs])),
        'mean_meth_hp2':    float(np.mean([c['meth_hp2']  for c in cpgs])),
        'min_qval':         float(np.min([c['qval'] for c in cpgs])),
        'frac_concordant':  frac_concordant,
    }


def cmd_call_loci(args):
    print(f'{args.sample}  {args.chrom}  call-loci')
    os.makedirs(args.outdir, exist_ok=True)

    hp1 = load_pileup(f'{args.pileup_dir}/{args.sample}_{args.chrom}_hp_pileup_1.bed.gz')
    hp2 = load_pileup(f'{args.pileup_dir}/{args.sample}_{args.chrom}_hp_pileup_2.bed.gz')
    assert (hp1['chrom'] == args.chrom).all() and (hp2['chrom'] == args.chrom).all(), \
        'pileup file contains rows outside the expected chromosome'
    print(f'  HP1: {len(hp1):,}  HP2: {len(hp2):,}')

    merged = hp1.merge(hp2, on=['chrom', 'start', 'end'], suffixes=('_hp1', '_hp2'))
    merged['cov_hp1'] = merged['n_mod_hp1'] + merged['n_canonical_hp1']
    merged['cov_hp2'] = merged['n_mod_hp2'] + merged['n_canonical_hp2']
    merged = merged[(merged['cov_hp1'] >= MIN_COV_PER_HAP) &
                     (merged['cov_hp2'] >= MIN_COV_PER_HAP)].copy()
    merged['meth_hp1']  = merged['n_mod_hp1'] / merged['cov_hp1']
    merged['meth_hp2']  = merged['n_mod_hp2'] / merged['cov_hp2']
    merged['delta']     = merged['meth_hp1'] - merged['meth_hp2']
    print(f'  CpGs covered >= {MIN_COV_PER_HAP}x on both haps: {len(merged):,}')

    merged = call_cpgs(merged, args.threads)
    dasm_cpgs = merged[(merged['qval'] < FDR_THRESH) &
                        (merged['delta'].abs() >= DELTA_FLOOR)]
    dasm_cpgs = dasm_cpgs.sort_values('start').reset_index(drop=True)
    print(f'  dASM CpGs (FDR<{FDR_THRESH} & |delta|>={DELTA_FLOOR}): {len(dasm_cpgs):,} / {len(merged):,} tested')

    loci = cluster_cpgs(dasm_cpgs, args.chrom, GAP_MERGE_BP)
    loci_df = pd.DataFrame(loci) if loci else pd.DataFrame(columns=[
        'chrom', 'start', 'end', 'locus_id', 'n_cpgs', 'mean_abs_delta',
        'mean_meth_hp1', 'mean_meth_hp2', 'min_qval', 'frac_concordant'])
    print(f'  dASM loci: {len(loci_df):,}')

    loci_df.to_csv(f'{args.outdir}/dasm_loci_{args.chrom}_fisher.tsv', sep='\t', index=False)

    with open(f'{args.outdir}/dasm_loci_{args.chrom}_fisher.bed', 'w') as f:
        for _, r in loci_df.iterrows():
            f.write(f"{r['chrom']}\t{max(0, int(r['start']) - FLANK)}\t{int(r['end']) + FLANK}\t{r['locus_id']}\n")

    summary = {
        'chrom': args.chrom, 'sample': args.sample,
        'n_cpgs_covered': int(len(merged)), 'n_dasm_cpgs': int(len(dasm_cpgs)),
        'n_dasm_loci': int(len(loci_df)),
    }
    with open(f'{args.outdir}/dasm_call_summary_{args.chrom}_fisher.json', 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


def load_extract_for_loci(path):
    cols   = ['read_id', 'ref_position', 'chrom', 'mod_qual', 'mod_code']
    dtypes = {'mod_qual': 'float32', 'ref_position': 'int32'}
    df = pd.read_csv(path, sep='\t', usecols=cols, dtype=dtypes, low_memory=False)
    df = df[(df['ref_position'] >= 0) & (df['mod_code'] == 'm')].copy()
    df['is_methylated'] = (df['mod_qual'] >= 0.5).astype(np.float32)
    return df


def build_epiallele_store(df_ex, loci_df, bam_path, flank):
    store = {}
    df_sorted = df_ex.sort_values('ref_position').reset_index(drop=True)
    pos_arr   = df_sorted['ref_position'].values

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for _, locus in loci_df.iterrows():
            lo = int(locus['start']) - flank
            hi = int(locus['end'])   + flank
            i0 = int(np.searchsorted(pos_arr, lo,    side='left'))
            i1 = int(np.searchsorted(pos_arr, hi + 1, side='left'))
            if i1 <= i0:
                continue
            region   = df_sorted.iloc[i0:i1]
            cpg_pos  = np.sort(region['ref_position'].unique())
            read_ids = region['read_id'].unique().tolist()
            if not len(read_ids) or not len(cpg_pos):
                continue

            hp_map = {}
            for read in bam.fetch(locus['chrom'], max(0, lo), hi):
                if read.has_tag('HP'):
                    hp_map[read.query_name] = int(read.get_tag('HP'))

            c_idx = {p: i for i, p in enumerate(cpg_pos)}
            r_idx = {r: i for i, r in enumerate(read_ids)}
            mat   = np.full((len(read_ids), len(cpg_pos)), np.nan, dtype=np.float32)
            ri_arr = region['read_id'].map(r_idx).values
            ci_arr = region['ref_position'].map(c_idx).values
            mat[ri_arr, ci_arr] = region['is_methylated'].values

            hp_labels = np.full(len(read_ids), np.nan)
            for j, rid in enumerate(read_ids):
                h = hp_map.get(rid)
                if h is not None:
                    hp_labels[j] = float(h)

            store[locus['locus_id']] = {
                'hp_labels':     hp_labels,
                'per_read_mean': np.nanmean(mat, axis=1),
                'locus':         locus.to_dict(),
            }
    return store


def compute_phasing(store, min_reads, auc_thresh, pval_thresh):
    results = []
    for lid, v in store.items():
        locus     = v['locus']
        hp_labels = v['hp_labels']
        per_read  = v['per_read_mean']

        hap1 = per_read[(hp_labels == 1) & ~np.isnan(per_read)]
        hap2 = per_read[(hp_labels == 2) & ~np.isnan(per_read)]

        auc, pval = np.nan, np.nan
        read_level_support = False
        if len(hap1) >= min_reads and len(hap2) >= min_reads:
            _, pval = mannwhitneyu(hap1, hap2, alternative='two-sided')
            try:
                auc = roc_auc_score([1] * len(hap1) + [0] * len(hap2), list(hap1) + list(hap2))
                auc = max(auc, 1 - auc)
            except Exception:
                pass
            read_level_support = (not np.isnan(auc)) and (auc >= auc_thresh) and (pval < pval_thresh)

        results.append({
            'locus_id': lid, **{k: locus[k] for k in
                ('chrom', 'start', 'end', 'n_cpgs', 'mean_abs_delta', 'min_qval')},
            'n_reads_hap1': len(hap1), 'n_reads_hap2': len(hap2),
            'mean_meth_hap1': float(np.mean(hap1)) if len(hap1) else np.nan,
            'mean_meth_hap2': float(np.mean(hap2)) if len(hap2) else np.nan,
            'auc': float(auc) if not np.isnan(auc) else np.nan,
            'pval': float(pval) if not np.isnan(pval) else np.nan,
            'read_level_support': read_level_support,
        })
    return pd.DataFrame(results)


def cmd_phase_loci(args):
    print(f'{args.sample}  {args.chrom}  phase-loci')
    loci_df = pd.read_csv(f'{args.outdir}/dasm_loci_{args.chrom}_fisher.tsv', sep='\t')
    if loci_df.empty:
        print('  no candidate loci — writing empty phasing output')
        pd.DataFrame().to_parquet(f'{args.outdir}/phasing_results_{args.chrom}_fisher.parquet')
        return

    print('  loading targeted modkit extract (--include-bed restricted)...')
    df_ex = load_extract_for_loci(args.extract_path)
    print(f'  extract rows: {len(df_ex):,}')

    print('  building epiallele store + HP labels from BAM...')
    store = build_epiallele_store(df_ex, loci_df, args.bam_path, FLANK)
    print(f'  store: {len(store):,} / {len(loci_df):,} loci')

    phasing_df = compute_phasing(store, MIN_HAP_READS, AUC_THRESH, PVAL_THRESH)
    n_phased = phasing_df['read_level_support'].sum() if len(phasing_df) else 0
    print(f'  read-level phased: {n_phased:,} / {len(phasing_df):,}')

    phasing_df.to_parquet(f'{args.outdir}/phasing_results_{args.chrom}_fisher.parquet', index=False)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)

    common = dict(required=True)
    p1 = sub.add_parser('call-loci')
    p1.add_argument('--chrom', **common)
    p1.add_argument('--sample', **common)
    p1.add_argument('--pileup-dir', **common)
    p1.add_argument('--outdir', **common)
    p1.add_argument('--threads', type=int, default=10)
    p1.set_defaults(func=cmd_call_loci)

    p2 = sub.add_parser('phase-loci')
    p2.add_argument('--chrom', **common)
    p2.add_argument('--sample', **common)
    p2.add_argument('--extract-path', **common)
    p2.add_argument('--bam-path', **common)
    p2.add_argument('--outdir', **common)
    p2.set_defaults(func=cmd_phase_loci)

    args = ap.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
