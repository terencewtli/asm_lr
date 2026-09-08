#!/usr/bin/env python3
"""
Build a finalized, cross-donor-replicated ASM locus list from the 6
confirmed-sane donors' significant beta-binomial regions, then run the
quick discovery checks (proximal/distal, effect size, imprinting overlap,
450K/meQTL overlap) against it.

Finalization rule (arbitrary, as instructed -- a first pass, not a
methodologically final choice): merge each donor's significant intervals
genome-wide (bedtools merge, -d 250 slop for cross-donor boundary jitter),
keep merged loci supported by >=2 of the 6 donors.

Input: tables/asm_analysis/six_donor_sig_regions.tsv (38,245 rows, 6 donors,
       all 22 autosomes, already has category/abs_delta from ASM01 notebook)
Output: tables/asm_analysis/finalized_asm_loci.tsv
        tables/phasing_qc/*.png figures
"""
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

PROJ = "/u/project/cluo/terencew/claude/project_ideas/asm_lr"
FIGDIR = f"{PROJ}/tables/phasing_qc"
BEDTOOLS = "/u/home/t/terencew/bin/bedtools"

sig = pd.read_csv(f"{PROJ}/tables/asm_analysis/six_donor_sig_regions.tsv", sep="\t")
print(f"loaded {len(sig)} significant regions across {sig.donor.nunique()} donors")

# --- merge across donors, genome-wide, per chromosome, with 250bp slop ---
bed_in = sig[["chrom", "start", "end", "donor"]].sort_values(["chrom", "start"])
with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
    bed_in.to_csv(f, sep="\t", header=False, index=False)
    bed_path = f.name

cmd = f"sort -k1,1 -k2,2n {bed_path} | {BEDTOOLS} merge -d 250 -c 4 -o distinct -i -"
out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
rows = []
for line in out.stdout.strip().split("\n"):
    chrom, start, end, donors = line.split("\t")
    donor_list = donors.split(",")
    rows.append((chrom, int(start), int(end), len(donor_list), donors))
merged = pd.DataFrame(rows, columns=["chrom", "start", "end", "n_donors", "donors"])
print(f"merged into {len(merged)} loci genome-wide")

finalized = merged[merged.n_donors >= 2].copy()
finalized["mid"] = (finalized.start + finalized.end) // 2
print(f"finalized (>=2/6 donors): {len(finalized)} loci")
finalized.to_csv(f"{PROJ}/tables/asm_analysis/finalized_asm_loci.tsv", sep="\t", index=False)

# --- attach best-available per-locus stats (delta, category) by picking the
# supporting donor row with the smallest qval at each finalized locus ---
sig_idx = sig.set_index(["chrom"])
attach_rows = []
for r in finalized.itertuples():
    cand = sig[(sig.chrom == r.chrom) & (sig.end >= r.start - 250) & (sig.start <= r.end + 250)]
    if cand.empty:
        attach_rows.append((np.nan, "NA", np.nan))
        continue
    best = cand.loc[cand.qval.idxmin()]
    attach_rows.append((best.abs_delta, best.category, best.nearest_snp_dist))
finalized[["best_abs_delta", "category", "nearest_snp_dist"]] = attach_rows

# --- proximal/distal breakdown ---
cat_counts = finalized["category"].value_counts()
print("\ncategory breakdown, finalized loci:")
print(cat_counts)

fig, ax = plt.subplots(figsize=(4, 4))
cat_counts.reindex(["proximal", "distal", "very_distal"]).plot.bar(ax=ax, color=["#4C72B0", "#DD8452", "#55A868"])
ax.set_ylabel("n loci")
ax.set_title(f"Finalized ASM loci by category (n={len(finalized)})")
plt.tight_layout()
plt.savefig(f"{FIGDIR}/fig_category_breakdown.png", dpi=150)
plt.close()

# --- effect size ---
fig, ax = plt.subplots(figsize=(4, 4))
finalized["best_abs_delta"].dropna().plot.hist(bins=30, ax=ax, color="#4C72B0")
ax.set_xlabel("|delta| (best supporting donor)")
ax.set_title("Effect size, finalized loci")
plt.tight_layout()
plt.savefig(f"{FIGDIR}/fig_effect_size.png", dpi=150)
plt.close()

# --- n_donors replication histogram ---
fig, ax = plt.subplots(figsize=(4, 4))
merged["n_donors"].value_counts().sort_index().plot.bar(ax=ax, color="#55A868")
ax.set_xlabel("n donors (of 6) supporting locus")
ax.set_ylabel("n merged loci (all, pre-filter)")
ax.set_title("Cross-donor replication, all merged loci")
plt.tight_layout()
plt.savefig(f"{FIGDIR}/fig_replication_hist.png", dpi=150)
plt.close()

# --- imprinting overlap, restricted to finalized loci ---
zink = pd.read_csv(f"{PROJ}/data/zink_2018_supp5_pofo_dmrs.csv", skiprows=2)
def overlaps_finalized(row, fin):
    sub = fin[fin.chrom == row.Chrom]
    return ((sub.start <= row.peakStop) & (sub.end >= row.peakStart)).any()
zink["hit_by_finalized"] = zink.apply(lambda r: overlaps_finalized(r, finalized), axis=1)
n_dmr_hit = zink["hit_by_finalized"].sum()
print(f"\nimprinted DMRs (n={len(zink)}) recovered by finalized list: {n_dmr_hit} ({100*n_dmr_hit/len(zink):.1f}%)")

# --- 450K probe + meQTL overlap ---
probes = pd.read_csv(f"{PROJ}/data/450k_probes_hg38.bed", sep="\t", header=None,
                      names=["chrom", "start", "end", "probe_id"])
meqtl_cpgs = set(open(f"{PROJ}/data/meqtl_mono_sig_cpgs.txt").read().split())

def nearest_probe(row, probes, window=250):
    sub = probes[(probes.chrom == row.chrom) & (probes.start >= row.start - window) & (probes.start <= row.end + window)]
    if sub.empty:
        return None
    return sub.iloc[(sub.start - row.mid).abs().argsort()].iloc[0].probe_id

finalized["nearest_450k_probe"] = finalized.apply(lambda r: nearest_probe(r, probes), axis=1)
finalized["has_450k_probe"] = finalized["nearest_450k_probe"].notna()
finalized["has_meqtl"] = finalized["nearest_450k_probe"].isin(meqtl_cpgs)
print(f"\nfinalized loci with a 450K probe within 250bp: {finalized.has_450k_probe.sum()} ({100*finalized.has_450k_probe.mean():.1f}%)")
print(f"finalized loci with a 450K probe that's also a monocyte meQTL: {finalized.has_meqtl.sum()} "
      f"({100*finalized.has_meqtl.sum()/max(finalized.has_450k_probe.sum(),1):.1f}% of probed loci)")
# background rate: fraction of ALL 450k probes genome-wide that are meqtl-sig
bg_rate = len(meqtl_cpgs) / len(probes)
print(f"background rate (meqtl_cpgs / all 450k probes genome-wide): {100*bg_rate:.1f}%")

fig, ax = plt.subplots(figsize=(4, 4))
vals = [100*bg_rate, 100*finalized.has_meqtl.sum()/max(finalized.has_450k_probe.sum(),1)]
ax.bar(["genome-wide\nbackground", "finalized ASM\nloci (probed)"], vals, color=["#999999", "#C44E52"])
ax.set_ylabel("% probes that are\nmonocyte meQTL-significant")
ax.set_title("meQTL enrichment at finalized loci")
plt.tight_layout()
plt.savefig(f"{FIGDIR}/fig_meqtl_enrichment.png", dpi=150)
plt.close()

finalized.to_csv(f"{PROJ}/tables/asm_analysis/finalized_asm_loci.tsv", sep="\t", index=False)
print(f"\nwrote tables/asm_analysis/finalized_asm_loci.tsv ({len(finalized)} loci)")
print("wrote figures to tables/phasing_qc/fig_*.png")
