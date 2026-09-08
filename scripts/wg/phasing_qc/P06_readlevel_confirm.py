#!/usr/bin/env python3
"""
Read-level confirmatory pass: re-test the pooled beta-binomial caller's
already-significant regions using readlevel.test_region_reads (Welch +
Mann-Whitney on per-read methylation fractions -- the molecule is the unit
of observation, so no dispersion parameter is needed and there is no
region-pooling pseudoreplication to worry about). This is a CONFIRMATION
step on a pre-screened candidate set, not a genome-wide test -- see
docs/2026-09-05_calibration_critique.md / 2026-09-07_next_steps... in
github/ont_asm_caller.

Input:
  - tables/dasm_betabinom/{sample}/regions_{chrom}.tsv (significant==True rows)
  - bam/haplotagged/by_chrom/{sample}_{chrom}_haplotagged.bam (for HP tags,
    restricted via -L to just the candidate regions -- fast, indexed)
  - modkit/{sample}_{chrom}_modkit_extract.tsv.gz (per-read, per-CpG calls;
    not haplotype-split, not position-indexed -- streamed once per chrom)

Output: tables/phasing_qc/readlevel_confirm_{sample}_{chrom}.tsv
        (one row per candidate region: pooled-test qval/delta alongside the
        read-level test's pval/delta/method, plus BH-FDR within this
        confirmatory candidate set)

Usage: python3 P06_readlevel_confirm.py <SAMPLE> <CHROM>
"""
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, "/u/project/cluo/terencew/claude/project_ideas/asm_lr/github/ont_asm_caller")
from ont_asm_caller.readlevel import test_region_reads

PROJ = "/u/project/cluo/terencew/claude/project_ideas/asm_lr"
SAMTOOLS = "/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools"
MOD_QUAL_THRESH = 0.5
MIN_CALLS = 1
MIN_READS = 3


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def get_hp_tags_for_regions(sample, chrom, regions):
    bam = f"{PROJ}/bam/haplotagged/by_chrom/{sample}_{chrom}_haplotagged.bam"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as bf:
        for r in regions.itertuples():
            bf.write(f"{chrom}\t{r.start}\t{r.end}\n")
        bed_path = bf.name

    cmd = [SAMTOOLS, "view", "-L", bed_path, bam]
    hp = {}
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        qname = fields[0]
        for f in fields[11:]:
            if f.startswith("HP:i:"):
                hp[qname] = int(f[5:])
                break
    proc.wait()
    return hp


def stream_calls_into_regions(sample, chrom, region_starts, region_ends):
    """Single pass over modkit extract; assign each mod_code=='m' call to the
    (unique, non-overlapping, sorted) region it falls in via bisect, and
    accumulate (read_id, pos, called) triples per region index."""
    path = f"{PROJ}/modkit/{sample}_{chrom}_modkit_extract.tsv.gz"
    buckets = [[] for _ in range(len(region_starts))]
    cmd = f"zcat {path} | awk 'BEGIN{{FS=\"\\t\"}} NR>1 && $12==\"m\" && $3>=0 {{print $1\"\\t\"$3\"\\t\"$11}}'"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    starts = region_starts
    ends = region_ends
    n_regions = len(starts)
    for line in proc.stdout:
        rid, pos_s, qual_s = line.rstrip("\n").split("\t")
        pos = int(pos_s)
        idx = np.searchsorted(starts, pos, side="right") - 1
        if idx < 0 or idx >= n_regions:
            continue
        if pos > ends[idx]:
            continue
        buckets[idx].append((rid, pos, float(qual_s) >= MOD_QUAL_THRESH))
    proc.wait()
    return buckets


def main():
    sample, chrom = sys.argv[1], sys.argv[2]

    regions = pd.read_csv(f"{PROJ}/tables/dasm_betabinom/{sample}/regions_{chrom}.tsv", sep="\t")
    regions = regions[regions.significant == True].copy()
    regions = regions.sort_values("start").reset_index(drop=True)
    print(f"{sample} {chrom}: {len(regions)} pooled-test-significant regions to confirm", flush=True)
    if len(regions) == 0:
        return

    print(f"{sample} {chrom}: recovering HP tags for candidate regions (samtools -L)...", flush=True)
    hp = get_hp_tags_for_regions(sample, chrom, regions)
    print(f"  {len(hp):,} reads with an HP tag overlapping a candidate region", flush=True)

    print(f"{sample} {chrom}: streaming modkit extract, assigning calls to regions...", flush=True)
    starts = regions["start"].to_numpy()
    ends = regions["end"].to_numpy()
    buckets = stream_calls_into_regions(sample, chrom, starts, ends)
    print(f"  done streaming", flush=True)

    results = []
    for i, r in enumerate(regions.itertuples()):
        rows = buckets[i]
        if not rows:
            continue
        df = pd.DataFrame(rows, columns=["read_id", "pos", "called"])
        df["hp"] = df["read_id"].map(hp)
        df = df.dropna(subset=["hp"])
        if df.empty:
            continue
        df["hp"] = df["hp"].astype(int)
        m = {}
        for hp_val in (1, 2):
            sub = df[df.hp == hp_val]
            if sub.empty:
                m[hp_val] = None
                continue
            piv = sub.pivot_table(index="read_id", columns="pos", values="called", aggfunc="first")
            m[hp_val] = piv.to_numpy(dtype=float)
        if m[1] is None or m[2] is None:
            continue
        res = test_region_reads(chrom, r.start, r.end, m[1], m[2], min_calls=MIN_CALLS, min_reads=MIN_READS)
        if res is None:
            continue
        results.append({
            "chrom": chrom, "start": r.start, "end": r.end,
            "pooled_qval": r.qval, "pooled_delta": r.delta,
            "readlevel_pval": res.pval, "readlevel_delta": res.delta,
            "readlevel_method": res.method,
            "n_reads1": res.n_reads1, "n_reads2": res.n_reads2,
        })
        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(regions)} candidates tested", flush=True)

    out = pd.DataFrame(results)
    if len(out):
        out["readlevel_qval"] = bh_fdr(out["readlevel_pval"].to_numpy())
        out["confirmed"] = out["readlevel_qval"] < 0.05
    outpath = f"{PROJ}/tables/phasing_qc/readlevel_confirm_{sample}_{chrom}.tsv"
    out.to_csv(outpath, sep="\t", index=False)

    print(f"\n{sample} {chrom}: {len(out)}/{len(regions)} candidates had usable read-level data")
    if len(out):
        print(f"  confirmed (readlevel BH-FDR<0.05): {out['confirmed'].sum()} ({100*out['confirmed'].mean():.1f}%)")
    print(f"wrote {outpath}")


if __name__ == "__main__":
    main()
