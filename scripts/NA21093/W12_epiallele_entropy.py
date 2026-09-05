"""W12: Epiallele entropy analysis.

For each dASM locus, compute single-molecule methylation structure metrics
from the per-read × per-CpG matrices in the epiallele stores, then test
whether genotype-driven loci (proximal/distal/very-distal) have different
entropy/bimodality structure than non-genetic loci.

Run modes:
  --chrom chrN   score one chromosome → tables/W12_metrics_chrN.csv
  (no --chrom)   merge all per-chrom CSVs → final tables/W12_epiallele_metrics.csv
                 + MWU stats + figures

Inputs:
  tables/dasm_loci_wg_annotated.parquet
  /u/project/cluo_scratch/terencew/claude/asm_lr/epiallele_store_chr*.pkl

Outputs:
  tables/W12_metrics_chrN.csv        — per-chrom metrics (intermediate)
  tables/W12_epiallele_metrics.csv   — merged per-locus metrics + classification
  tables/W12_mwu_stats.tsv           — MWU test results with FDR correction
  figures/W12_entropy_by_class.pdf
  figures/W12_hpsep_vs_entropy.pdf
  figures/W12_bimodality_cdf.pdf
"""
import argparse
import os
import pickle
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

projdir  = '/u/project/cluo/terencew/claude/project_ideas/asm_lr'
scratch  = '/u/project/cluo_scratch/terencew/claude/asm_lr'
tabdir   = f'{projdir}/tables'
figdir   = f'{projdir}/figures'
os.makedirs(figdir, exist_ok=True)

CHROMS = [f'chr{i}' for i in range(1, 23)]

CLS_ORDER  = ['proximal-genotype', 'distal-genotype', 'very-distal-genotype', 'non-genetic']
CLS_COLORS = {
    'proximal-genotype':    '#d62728',
    'distal-genotype':      '#ff7f0e',
    'very-distal-genotype': '#9467bd',
    'non-genetic':          '#aec7e8',
}

def ts():
    return datetime.now().strftime('[%H:%M:%S]')


# ---------------------------------------------------------------------------
# Metric functions
# ---------------------------------------------------------------------------

def per_cpg_entropy(matrix):
    # NaN = read doesn't cover that CpG — use nanmean so missing positions
    # don't propagate through entropy calculation.
    p = np.nanmean(matrix, axis=0)
    valid = ~np.isnan(p)
    if valid.sum() == 0:
        return np.nan
    p = np.clip(p[valid], 1e-6, 1 - 1e-6)
    h = -p * np.log2(p) - (1 - p) * np.log2(1 - p)
    return float(h.mean())


def per_read_entropy(matrix):
    # Flatten and drop NaN (positions not covered by a given read).
    flat = matrix.ravel()
    flat = flat[~np.isnan(flat)]
    if len(flat) == 0:
        return np.nan
    p = np.clip(flat, 1e-6, 1 - 1e-6)
    h = -p * np.log2(p) - (1 - p) * np.log2(1 - p)
    return float(h.mean())


def bimodality_coef(arr):
    if len(arr) < 4:
        return np.nan
    s = stats.skew(arr)
    k = stats.kurtosis(arr)
    return float((s ** 2 + 1) / (k + 3))


def compute_locus_metrics(entry):
    mat     = entry['matrix']        # (n_reads, n_cpgs) float32, NaN = no coverage
    hp      = entry['hp_labels']     # (n_reads,) float64: 1.0 / 2.0 / nan
    pr_mean = entry['per_read_mean'] # (n_reads,) float32

    if mat.shape[0] < 4 or mat.shape[1] < 1:
        return None

    hp1 = hp == 1
    hp2 = hp == 2

    result = {
        'n_reads':          int(len(hp)),
        'n_hp1':            int(hp1.sum()),
        'n_hp2':            int(hp2.sum()),
        'per_cpg_entropy':  per_cpg_entropy(mat),
        'per_read_entropy': per_read_entropy(mat),
        'bimodality_coef':  bimodality_coef(pr_mean),
    }

    if hp1.sum() >= 2 and hp2.sum() >= 2:
        # nanmean per haplotype so missing-coverage NaNs don't poison the score
        m1 = float(np.nanmean(mat[hp1]))
        m2 = float(np.nanmean(mat[hp2]))
        result['hp_sep'] = abs(m1 - m2)
        result['within_hp_var'] = float(np.mean([
            np.nanmean(mat[hp1], axis=1).var(),
            np.nanmean(mat[hp2], axis=1).var(),
        ]))
    else:
        result['hp_sep']        = np.nan
        result['within_hp_var'] = np.nan

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_chrom(chrom):
    """Score one chromosome and write tables/W12_metrics_{chrom}.csv."""
    store_path = f'{scratch}/epiallele_store_{chrom}.pkl'
    if not os.path.exists(store_path):
        print(f'{ts()} ERROR: missing {store_path}')
        return

    print(f'{ts()} Loading {chrom} store...')
    with open(store_path, 'rb') as fh:
        store = pickle.load(fh)
    print(f'  {len(store):,} loci in store')

    rows = []
    n_ok = n_skip = 0
    for locus_id, entry in store.items():
        m = compute_locus_metrics(entry)
        if m is not None:
            m['locus_id'] = locus_id
            rows.append(m)
            n_ok += 1
        else:
            n_skip += 1

    print(f'{ts()} {chrom}: {n_ok:,} computed  {n_skip:,} skipped (< 4 reads or 0 CpGs)')

    out = f'{tabdir}/W12_metrics_{chrom}.csv'
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f'{ts()} Saved {out}')


def run_merge():
    """Merge per-chrom CSVs, run MWU tests, produce figures."""
    print(f'{ts()} Loading annotated loci...')
    df = pd.read_parquet(f'{tabdir}/dasm_loci_wg_annotated.parquet')
    print(f'  {len(df):,} loci  |  classes: {df["classification"].value_counts().to_dict()}')

    print(f'{ts()} Merging per-chrom metrics...')
    parts = []
    for chrom in CHROMS:
        p = f'{tabdir}/W12_metrics_{chrom}.csv'
        if not os.path.exists(p):
            print(f'  WARNING: missing {p}, skipping')
            continue
        parts.append(pd.read_csv(p))
    metrics_df = pd.concat(parts, ignore_index=True)
    print(f'{ts()} Metrics for {len(metrics_df):,} loci across {len(parts)} chroms')

    merged = df.merge(metrics_df, on='locus_id', how='left')
    for c in CLS_ORDER:
        n = merged[merged['classification'] == c]['per_cpg_entropy'].notna().sum()
        total_c = (merged['classification'] == c).sum()
        print(f'  {c}: {n:,}/{total_c:,} with entropy')

    # --- MWU tests --------------------------------------------------------
    print(f'{ts()} Running Mann-Whitney U tests...')
    METRICS = ['per_cpg_entropy', 'per_read_entropy', 'hp_sep', 'within_hp_var', 'bimodality_coef']
    ref = merged[merged['classification'] == 'non-genetic']
    test_rows = []
    for cls in ['proximal-genotype', 'distal-genotype', 'very-distal-genotype']:
        sub = merged[merged['classification'] == cls]
        for m in METRICS:
            a = sub[m].dropna()
            b = ref[m].dropna()
            if len(a) < 5 or len(b) < 5:
                continue
            stat, p = mannwhitneyu(a, b, alternative='two-sided')
            test_rows.append({
                'class':      cls,
                'metric':     m,
                'n_cls':      len(a),
                'n_ref':      len(b),
                'median_cls': a.median(),
                'median_ref': b.median(),
                'direction':  'higher' if a.median() > b.median() else 'lower',
                'pval':       p,
            })
    stat_df = pd.DataFrame(test_rows)
    _, stat_df['padj'], _, _ = multipletests(stat_df['pval'], method='fdr_bh')
    stat_df = stat_df.sort_values(['metric', 'class'])
    print(stat_df[['class', 'metric', 'median_cls', 'median_ref', 'direction', 'padj']].to_string(index=False))
    stat_df.to_csv(f'{tabdir}/W12_mwu_stats.tsv', sep='\t', index=False)

    # --- figures ----------------------------------------------------------
    print(f'{ts()} Plotting...')
    plot_cls = [c for c in CLS_ORDER if c in merged['classification'].unique()]

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, metric, label in zip(
        axes,
        ['per_cpg_entropy', 'per_read_entropy'],
        ['Per-CpG entropy (bits)', 'Per-read entropy (bits)'],
    ):
        data_by_cls = [merged[merged['classification'] == c][metric].dropna().values for c in plot_cls]
        valid = [(d, c) for d, c in zip(data_by_cls, plot_cls) if len(d) > 1]
        if not valid:
            continue
        data_v, cls_v = zip(*valid)
        parts_v = ax.violinplot(list(data_v), positions=range(len(cls_v)),
                                showmedians=True, showextrema=False)
        for pc, c in zip(parts_v['bodies'], cls_v):
            pc.set_facecolor(CLS_COLORS[c]); pc.set_alpha(0.7)
        parts_v['cmedians'].set_color('black'); parts_v['cmedians'].set_linewidth(1.5)
        ax.set_xticks(range(len(cls_v)))
        ax.set_xticklabels([c.replace('-', '\n') for c in cls_v], fontsize=7)
        ax.set_ylabel(label)
    plt.tight_layout()
    plt.savefig(f'{figdir}/W12_entropy_by_class.pdf', bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(5, 4))
    for cls in plot_cls:
        sub = merged[merged['classification'] == cls][['hp_sep', 'per_read_entropy']].dropna()
        if cls == 'non-genetic':
            sub = sub.sample(min(10_000, len(sub)), random_state=1)
        ax.scatter(sub['hp_sep'], sub['per_read_entropy'],
                   s=1, alpha=0.3, label=cls, color=CLS_COLORS[cls])
    ax.set_xlabel('HP separation (|mean_HP1 − mean_HP2|)')
    ax.set_ylabel('Per-read entropy (bits)')
    ax.legend(markerscale=6, fontsize=7)
    plt.tight_layout()
    plt.savefig(f'{figdir}/W12_hpsep_vs_entropy.pdf', bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(5, 4))
    for cls in plot_cls:
        vals = np.sort(merged[merged['classification'] == cls]['bimodality_coef'].dropna())
        if len(vals) == 0:
            continue
        ax.plot(vals, np.linspace(0, 1, len(vals)), label=cls, color=CLS_COLORS[cls])
    ax.axvline(5 / 9, color='k', linestyle='--', linewidth=0.8, label='threshold (5/9)')
    ax.set_xlabel('Bimodality coefficient'); ax.set_ylabel('CDF'); ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(f'{figdir}/W12_bimodality_cdf.pdf', bbox_inches='tight')
    plt.close()

    # --- save & summary ---------------------------------------------------
    out_cols = ['locus_id', 'classification', 'per_cpg_entropy', 'per_read_entropy',
                'hp_sep', 'within_hp_var', 'bimodality_coef']
    merged[out_cols].to_csv(f'{tabdir}/W12_epiallele_metrics.csv', index=False)
    print(f'{ts()} Saved {tabdir}/W12_epiallele_metrics.csv')

    print(f'\n{ts()} Summary:')
    for cls in CLS_ORDER:
        sub = merged[merged['classification'] == cls]
        print(f'  {cls:30s}  entropy={sub["per_cpg_entropy"].median():.3f}'
              f'  hp_sep={sub["hp_sep"].median():.3f}'
              f'  bimodality={sub["bimodality_coef"].median():.3f}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--chrom', default=None,
                    help='Score one chromosome (e.g. chr1). Omit to run final merge.')
    args = ap.parse_args()

    print(f'{ts()} Start W12_epiallele_entropy.py')
    if args.chrom:
        run_chrom(args.chrom)
    else:
        run_merge()
    print(f'{ts()} Done')
