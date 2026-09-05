#!/usr/bin/env python3
"""
Per-chromosome dASM locus calling.

Replaces the papermill-executed W05_dasm_loci.ipynb notebook.
Stage 1: delta threshold (|HP1 - HP2| >= 0.3) — effect-size filter at CpG level.
Stage 2: Mann-Whitney U + AUC on per-read methylation — statistical filter at locus level.

Outputs:
  tables/phasing_results_{chrom}.parquet
  tables/dasm_loci_{chrom}.tsv
  tables/phasing_summary_{chrom}.json
"""

import argparse
import json
import os
import numpy as np
import pandas as pd
import pysam
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

# ── thresholds ────────────────────────────────────────────────────────────────
MIN_COV_PER_HAP = 5
DELTA_THRESH    = 0.3   # |HP1_meth - HP2_meth| threshold at CpG level
GAP_MERGE_BP    = 500
METH_THRESH     = 0.5   # mod_qual ≥ 0.5 → methylated (modkit extract)
FLANK           = 500
MIN_HAP_READS   = 3
AUC_THRESH      = 0.70
PVAL_THRESH     = 0.05  # Mann-Whitney p-value threshold at locus level

BEDMETHYL_COLS = [
    'chrom', 'start', 'end', 'mod_code', 'depth', 'strand',
    'start2', 'end2', 'color', 'depth2', 'pct_mod',
    'n_mod', 'n_canonical', 'n_other_mod', 'n_delete',
    'n_fail', 'n_diff', 'n_no_call',
]


def load_pileup(path):
    return pd.read_csv(
        path, sep='\t', header=None, names=BEDMETHYL_COLS,
        usecols=['chrom', 'start', 'end', 'depth', 'pct_mod'],
    )


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
        'mean_abs_delta':   float(np.mean([c['abs_delta'] for c in cpgs])),
        'mean_meth_hp1':    float(np.mean([c['meth_hp1']  for c in cpgs])),
        'mean_meth_hp2':    float(np.mean([c['meth_hp2']  for c in cpgs])),
        'frac_concordant':  frac_concordant,
    }


def load_modkit_extract(path, chrom, meth_thresh):
    cols   = ['read_id', 'ref_position', 'chrom', 'mod_qual', 'mod_code']
    dtypes = {'mod_qual': 'float32', 'ref_position': 'int32'}
    chunks = []
    for chunk in pd.read_csv(path, sep='\t', usecols=cols, dtype=dtypes,
                             low_memory=False, chunksize=1_000_000):
        sub = chunk[chunk['chrom'] == chrom]
        if len(sub):
            chunks.append(sub)
    df = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=cols)
    df = df[(df['ref_position'] >= 0) & (df['mod_code'] == 'm')].copy()
    df['is_methylated'] = (df['mod_qual'] >= meth_thresh).astype(np.float32)
    return df


def load_hp_tags(bam_path, chrom):
    hp_map = {}
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam.fetch(chrom):
            if read.has_tag('HP'):
                hp_map[read.query_name] = int(read.get_tag('HP'))
    return hp_map


def build_epiallele_store(df_ex, loci_df, hp_map, flank):
    store = {}
    df_sorted = df_ex.sort_values('ref_position').reset_index(drop=True)
    pos_arr   = df_sorted['ref_position'].values

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

        c_idx = {p: i for i, p in enumerate(cpg_pos)}
        r_idx = {r: i for i, r in enumerate(read_ids)}
        mat   = np.full((len(read_ids), len(cpg_pos)), np.nan, dtype=np.float32)
        hp_labels = np.full(len(read_ids), np.nan)

        ri_arr = region['read_id'].map(r_idx).values
        ci_arr = region['ref_position'].map(c_idx).values
        mat[ri_arr, ci_arr] = region['is_methylated'].values

        for j, rid in enumerate(read_ids):
            h = hp_map.get(rid)
            if h is not None:
                hp_labels[j] = float(h)

        store[locus['locus_id']] = {
            'matrix':        mat,
            'cpg_positions': cpg_pos,
            'read_ids':      read_ids,
            'hp_labels':     hp_labels,
            'per_read_mean': np.nanmean(mat, axis=1),
            'locus':         locus.to_dict(),
        }
    return store


def load_het_snps(vcf_path, chrom):
    rows = []
    vcf = pysam.VariantFile(vcf_path)
    for rec in vcf.fetch(chrom):
        if len(rec.alleles) != 2:
            continue
        s  = list(rec.samples.values())[0]
        gt = s['GT']
        if None in gt or set(gt) in ({0}, {1}):
            continue
        rows.append({'pos0': rec.start, 'phased': s.phased})
    vcf.close()
    return pd.DataFrame(rows)


def compute_phasing(loci_df, store, snp_pos_arr, min_reads, auc_thresh, pval_thresh):
    results = []
    for lid, v in store.items():
        locus    = v['locus']
        hp_labels = v['hp_labels']
        per_read  = v['per_read_mean']

        hap1 = per_read[(hp_labels == 1) & ~np.isnan(per_read)]
        hap2 = per_read[(hp_labels == 2) & ~np.isnan(per_read)]

        dists = np.abs(snp_pos_arr - (locus['start'] + locus['end']) // 2) \
                if len(snp_pos_arr) > 0 else np.array([np.nan])
        nearest_snp_dist = int(dists.min()) if len(dists) > 0 else np.nan

        auc, pval = np.nan, np.nan
        read_level_support = False
        if len(hap1) >= min_reads and len(hap2) >= min_reads:
            _, pval = mannwhitneyu(hap1, hap2, alternative='two-sided')
            try:
                auc = roc_auc_score([1] * len(hap1) + [0] * len(hap2),
                                    list(hap1) + list(hap2))
                auc = max(auc, 1 - auc)
            except Exception:
                pass
            read_level_support = (not np.isnan(auc)) and (auc >= auc_thresh) and (pval < pval_thresh)

        results.append({
            'locus_id':         lid,
            'chrom':            locus['chrom'],
            'locus_start':      locus['start'],
            'locus_end':        locus['end'],
            'n_cpgs':           locus['n_cpgs'],
            'mean_abs_delta':   locus['mean_abs_delta'],
            'n_reads_hap1':     len(hap1),
            'n_reads_hap2':     len(hap2),
            'mean_meth_hap1':   float(np.mean(hap1)) if len(hap1) else np.nan,
            'mean_meth_hap2':   float(np.mean(hap2)) if len(hap2) else np.nan,
            'auc':              float(auc)  if not np.isnan(auc)  else np.nan,
            'pval':             float(pval) if not np.isnan(pval) else np.nan,
            'read_level_support':        read_level_support,
            'nearest_snp_dist': nearest_snp_dist,
        })
    return pd.DataFrame(results)


def add_rosenski_overlap(phasing_df, projdir, chrom):
    def _load(fname):
        p = f'{projdir}/data/{fname}'
        df = pd.read_csv(p, sep='\t', header=None, compression='gzip',
                         usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
        return df[df['chrom'] == chrom].reset_index(drop=True)

    b205 = _load('GSE186458_blocks.s205.bed.gz')
    b207 = _load('GSE186458_blocks.s207.hg38.bed.gz')

    def overlaps(qs, qe, ref):
        return ((ref['start'] < qe) & (ref['end'] > qs)).any()

    phasing_df['in_rosenski'] = phasing_df.apply(
        lambda r: overlaps(r['locus_start'], r['locus_end'], b205) or
                  overlaps(r['locus_start'], r['locus_end'], b207),
        axis=1,
    )
    return phasing_df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--chrom',   required=True)
    ap.add_argument('--sample',  default='NA21093')
    ap.add_argument('--projdir', default='/u/project/cluo/terencew/claude/project_ideas/asm_lr')
    ap.add_argument('--outdir',  default=None,
                    help='Output directory (default: {projdir}/tables)')
    args = ap.parse_args()

    chrom   = args.chrom
    sample  = args.sample
    projdir = args.projdir

    print(f'{sample}  {chrom}')
    out = args.outdir or f'{projdir}/tables'
    os.makedirs(out, exist_ok=True)

    pileup_dir   = f'{projdir}/pileup'
    extract_path = f'{projdir}/modkit/{sample}_wg_modkit_extract.tsv.gz'
    bam_path     = f'{projdir}/bam/{sample}_{chrom}_haplotagged.bam'
    vcf_path     = f'{projdir}/vcf/{sample}_{chrom}_phased.vcf.gz'

    # ── step 1: load pileups ─────────────────────────────────────────────────
    print('loading pileups...')
    hp1 = load_pileup(f'{pileup_dir}/{sample}_{chrom}_hp_pileup_1.bed.gz')
    hp2 = load_pileup(f'{pileup_dir}/{sample}_{chrom}_hp_pileup_2.bed.gz')
    print(f'  HP1: {len(hp1):,}  HP2: {len(hp2):,}')

    merged = hp1.merge(hp2, on=['chrom', 'start', 'end'], suffixes=('_hp1', '_hp2'))
    mask   = (merged['depth_hp1'] >= MIN_COV_PER_HAP) & (merged['depth_hp2'] >= MIN_COV_PER_HAP)
    merged = merged[mask].copy()
    merged['meth_hp1']  = merged['pct_mod_hp1'] / 100
    merged['meth_hp2']  = merged['pct_mod_hp2'] / 100
    merged['delta']     = merged['meth_hp1'] - merged['meth_hp2']
    merged['abs_delta'] = merged['delta'].abs()
    print(f'  CpGs covered ≥{MIN_COV_PER_HAP}x on both haps: {len(merged):,}')

    # ── step 2: delta filter ─────────────────────────────────────────────────
    dasm_cpgs = merged[merged['abs_delta'] >= DELTA_THRESH].sort_values('start').reset_index(drop=True)
    print(f'  dASM CpGs (|Δ|≥{DELTA_THRESH}): {len(dasm_cpgs):,}')

    # ── step 3: cluster ───────────────────────────────────────────────────────
    loci    = cluster_cpgs(dasm_cpgs, chrom, GAP_MERGE_BP)
    loci_df = pd.DataFrame(loci)
    if loci_df.empty:
        print('  no dASM loci found — writing empty outputs')
        loci_df = pd.DataFrame(columns=[
            'chrom', 'start', 'end', 'locus_id', 'n_cpgs',
            'mean_abs_delta', 'mean_meth_hp1', 'mean_meth_hp2',
        ])
    else:
        n_conc = (loci_df['frac_concordant'] == 1.0).sum()
        print(f'  dASM loci: {len(loci_df):,}  (median {loci_df["n_cpgs"].median():.0f} CpGs/locus)')
        print(f'  concordant (fc=1): {n_conc:,} / {len(loci_df):,}  ({n_conc/len(loci_df):.1%})')

    # ── step 4: modkit extract ────────────────────────────────────────────────
    print('loading modkit extract (chunked)...')
    df_ex = load_modkit_extract(extract_path, chrom, METH_THRESH)
    print(f'  extract rows: {len(df_ex):,}')

    # ── step 5: HP tags ────────────────────────────────────────────────────────
    print('reading HP tags from haplotagged BAM...')
    hp_map = load_hp_tags(bam_path, chrom)
    df_ex['hp'] = df_ex['read_id'].map(hp_map).astype('Int8')
    print(f'  reads with HP tag: {len(hp_map):,}')

    # ── step 6: epiallele store ────────────────────────────────────────────────
    if loci_df.empty:
        store = {}
    else:
        print('building epiallele store...')
        store = build_epiallele_store(df_ex, loci_df, hp_map, FLANK)
        print(f'  store: {len(store):,} loci')

    # ── step 7: het SNPs ────────────────────────────────────────────────────────
    df_snps  = load_het_snps(vcf_path, chrom)
    snp_pos  = df_snps['pos0'].values if len(df_snps) else np.array([])
    print(f'  het SNPs: {len(df_snps):,}')

    # ── step 8: per-locus AUC phasing ──────────────────────────────────────────
    print('computing per-locus AUC...')
    phasing_df = compute_phasing(loci_df, store, snp_pos,
                                 MIN_HAP_READS, AUC_THRESH, PVAL_THRESH)
    n_phased = phasing_df['read_level_support'].sum() if len(phasing_df) else 0
    print(f'  phased: {n_phased:,} / {len(phasing_df):,}')

    # ── step 9: Rosenski overlap ────────────────────────────────────────────────
    if len(phasing_df):
        phasing_df = add_rosenski_overlap(phasing_df, projdir, chrom)
        n_ros = phasing_df['in_rosenski'].sum()
        print(f'  in Rosenski: {n_ros:,}')
    else:
        phasing_df['in_rosenski'] = False

    # ── save ────────────────────────────────────────────────────────────────────
    phasing_df.to_parquet(f'{out}/phasing_results_{chrom}.parquet', index=False)
    loci_df.to_csv(f'{out}/dasm_loci_{chrom}.tsv', sep='\t', index=False)

    summary = {
        'chrom':              chrom,
        'sample':             sample,
        'n_cpgs_covered':     int(len(merged)),
        'n_dasm_cpgs':        int(len(dasm_cpgs)),
        'n_dasm_loci':        int(len(loci_df)),
        'n_phased':           int(n_phased),
        'n_in_rosenski':      int(phasing_df['in_rosenski'].sum()) if len(phasing_df) else 0,
    }
    with open(f'{out}/phasing_summary_{chrom}.json', 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
