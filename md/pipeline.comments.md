### Comments on Steps 3–4: W03_phase_wg + W05_dasm_loci (dASM calling)

---

#### Step 3 ordering — your instinct is correct

Yes, W02 (pileup) must run after W03 (phasing) because pileup needs the
haplotagged BAM. The numbering is slightly misleading. The actual dependency
graph is:

```
W00 (align) ──► W01 (extract)          ──► W05 (dasm)
            └─► W03 (phase/haplotag) ──► W02 (pileup) ──► W05 (dasm)
```

W01 and W03 can run in parallel after W00. W02 waits for W03. W05 waits for
both W01 and W02.

---

#### dASM calling — what the pipeline actually does

Your concern is valid, but there's more statistical rigor here than just a
delta threshold. The calling happens in two stages:

**Stage 1 — delta filter at the CpG level (notebook cell 3):**
- Merges HP1 and HP2 pileup BEDs
- Requires ≥5 reads per haplotype at each CpG
- Calls a CpG as candidate dASM if |HP1_pct_mod − HP2_pct_mod| ≥ 30%
- Clusters adjacent candidate CpGs within 500 bp → loci

**Stage 2 — statistical test at the locus level (notebook cell 8):**
- For each locus, pulls per-read methylation from the modkit extract
- Assigns HP1/HP2 to each read via the HP tag
- Tests whether HP1 reads have different mean methylation than HP2 reads
  using Mann-Whitney U (two-sided)
- Computes an AUC: how well does per-read mean methylation predict haplotype?
- A locus passes only if: AUC ≥ 0.70 AND p < 0.05

So the final output (`is_phased = True`) requires both a biological signal (delta)
and a statistical test (Mann-Whitney). The loci in the output parquet are the
statistically validated subset.

---

#### What's legitimately fragile — and what Fisher's exact would fix

Your instinct about Fisher's exact is correct, but it addresses a different
problem than what the Mann-Whitney does. Here's the gap:

**The coverage-delta interaction at Stage 1 is the real fragility.** At 5x
per haplotype (the minimum), a 30% delta requires something like 2/5 vs 0/5
modified reads — that's a delta of 0.4, which passes the threshold, but with
only 5 observations this is a very noisy estimate. The pileup fractions at
low coverage have huge uncertainty that the delta threshold ignores entirely.

Fisher's exact at the CpG level would fix this properly:

```
           HP1   HP2
modified    m1    m2
canonical   c1    c2
```

Fisher's exact tests whether the ratio of modified:canonical reads differs
between haplotypes, and crucially it naturally accounts for coverage — a
3/5 vs 0/5 split gives p ~0.08 (not significant), whereas 30/50 vs 0/50
gives p ~1e-10 (very significant). The current delta threshold treats both
the same.

**The second gap is multiple testing correction on the Mann-Whitney p-values.**
With ~50,000+ loci genome-wide, `p < 0.05` at the locus level without
BH/Bonferroni correction means ~2,500 false positives by chance. The
AUC ≥ 0.70 filter mitigates this somewhat (it's a pretty strict effect-size
requirement), but it's not a substitute for MTC.

**Summary of gaps:**
1. Stage 1 delta filter ignores coverage uncertainty — Fisher's exact per CpG
   would make coverage-dependent calls
2. Stage 2 Mann-Whitney has no genome-wide multiple testing correction
3. Min coverage of 5x per haplotype is quite low (~20–25x per haplotype is
   where the estimates become reliable)

---

#### Are the per-chromosome notebooks just the same thing parallelized?

Yes exactly. `W05_dasm_loci.ipynb` is a single parameterized notebook run
via papermill with `chrom` as the injected parameter. The script W04_dasm_loci_wg.sh
(confusingly, "W04" in the script name vs "W05" in the notebook name — a
numbering artifact) runs it as a 22-element array job, one task per autosome.
Each task calls papermill with `-p chrom chrN`.

The per-chromosome executed notebooks in `notebooks/wg/executed/` are the
outputs of those papermill runs. The `tables/dasm_loci_chrN.tsv` files are
the per-chromosome outputs, and W05_merge_wg.sh concatenates them into
`tables/dasm_loci_wg.parquet`.

You don't need to refactor this — it's already an array job. The only thing
worth changing is whether you want to add Fisher's exact at Stage 1 and BH
correction at Stage 2.

---

### Comments on Step 2: W01_modkit_extract_wg.sh + W02_modkit_pileup_wg.sh

---

#### extract vs pileup — what each produces

Your description is right. They're two different views of the same underlying
modification calls:

**`modkit extract`** → one row per (read × CpG position)

Each row answers: "for this specific read, at this specific CpG on the genome,
was the cytosine methylated?" The output has columns like: read_id, chrom,
ref_position, strand, call_code (m = 5mC, h = 5hmC, . = canonical), score
(0–255 probability), and `query_kmer` (more below). This is the raw
per-molecule view — you can see the methylation pattern of individual DNA
molecules, which is what makes epiallele analysis possible.

**`modkit pileup`** → one row per (CpG position × haplotype)

Aggregates across all reads at each position. Output columns are: chrom, start,
end, coverage, fraction_modified, n_modified, n_canonical, n_other_modified.
With `--partition-tag HP` it produces separate BED files for HP1 and HP2.
This is the summarized view used for per-locus methylation comparisons —
e.g., "at this CpG, haplotype 1 is 80% methylated and haplotype 2 is 10%."

So: extract is input to the epiallele notebooks (need single-molecule patterns),
pileup is input to the ASM analysis (need per-haplotype methylation fractions).

---

#### What is query_kmer?

It's the sequence context around the modified base in the **read itself**
(not the reference). Specifically, it's a short k-mer (typically 9–11 bp)
centered on the C being called, taken from the basecalled read sequence.

ONT modification calling is context-dependent — the electrical signal from a
5mC is slightly different depending on neighboring bases, and the model that
converts signal to modification probability was trained on specific k-mers. The
`query_kmer` lets you see what context the model actually saw for each call.

In practice you'd use it to:
- Filter unreliable calls in low-complexity contexts (e.g., runs of the same
  base that confuse the signal model)
- Verify that calls at CpG positions actually have a CpG in the kmer (a sanity
  check for indel-shifted positions)
- Debug unexpected methylation patterns in specific sequence contexts

For most analyses you can ignore it and just use the score threshold.

---

#### Why unhaplotagged BAM for extract but haplotagged BAM for pileup?

This is the most important design choice in this step.

`modkit extract` runs on the **unhaplotagged** BAM. This is fine because the
extract output keeps the `read_id` for every row. Haplotype information can be
joined in later by matching read IDs against the haplotagging output. The
downstream epiallele notebooks do exactly this — they assign haplotypes at
analysis time rather than at extraction time.

`modkit pileup` runs on the **haplotagged** BAM with `--partition-tag HP`.
The HP tag (HP:i:1 or HP:i:2) is what WhatsHap writes onto each read during
haplotagging. Modkit reads that tag and routes each read's calls to the
appropriate output BED before aggregating. This is the only way to get
haplotype-stratified per-position methylation fractions — once you collapse
reads into a position-level count, you've destroyed the per-read identity and
there's no way to recover which haplotype contributed which calls.

The actual execution order in the pipeline reflects this dependency:

```
W00 (align)  →  W01 extract   (unhaplotagged — can run immediately)
W00 (align)  →  W03 phase/haplotag  →  W02 pileup  (needs HP tags)
```

W01 and W03 can be submitted in parallel after W00. W02 must wait for W03.
This is why the script numbering looks slightly out of order (W01 before W02
but W02 depends on W03).

---

### Comments on Step 1: W00_align_wg.sh

---

**Your understanding is correct:** the script downloads three unaligned BAMs from
HPRC, converts to FASTQ (preserving methylation tags), remaps to hg38, and merges.

---

#### Why remap? Aren't the HPRC BAMs already aligned?

The files downloaded here (`...pass.bam`) are **unaligned BAMs** — they're in
BAM format purely as a container for the methylation tags, not because they've
been mapped to a genome. Guppy outputs basecalled reads in BAM format with
two critical tags:

- `MM:Z:` — modification positions relative to each read (e.g., which C's are 5mC)
- `ML:B:` — per-modification probability scores (0–255)

BAM was chosen over FASTQ because FASTQ has no standardized way to carry
per-base modification calls. So the "BAM" you're downloading is really a fancy
FASTQ with methylation attached.

The remap step in the script does:
```
samtools bam2fq -T MM,ML  →  minimap2 -y  →  sorted BAM
```

`-T MM,ML` promotes the methylation tags into the FASTQ header line.
`-y` tells minimap2 to copy those tags verbatim into its output BAM.
The modification probabilities from Guppy are thus preserved in the final
aligned BAM and carry through to all downstream steps (modkit, phasing, etc.).

**The script actually verifies this** — lines 81–86 check that MM tags are
present in the merged BAM and abort if not.

---

#### Why runs 1, 3, 5? What happened to 2 and 4?

HPRC samples were sequenced across multiple PromethION flow cells, numbered
sequentially as runs were submitted. Runs 2 and 4 for HG01258 likely exist
but belong to different individuals or were run interleaved — the HPRC
distributed sequencing across multiple centers and the run numbering reflects
flow cell order in the lab's tracking system, not a guaranteed contiguous
sequence per sample. You'd want to check the HPRC data portal for HG01258
to confirm whether 2 and 4 exist and if so why they're excluded — it could be
failed flow cells, lower quality basecalls, or simply not yet released.

Having 3 runs is normal. A single PromethION flow cell typically yields
~15–20x human genome coverage, so 3 runs gives you ~45–60x total, which is
the depth this project targets.

---

#### Coverage, read length, error rate

For HPRC HG01258 with Guppy 6 `sup` (super-accuracy) basecalling:

| Metric | Typical value |
|---|---|
| Coverage per run | ~15–20x |
| Total (3 runs) | ~45–60x |
| Per-haplotype (after phasing) | ~20–30x |
| Median read length | ~15–25 kb |
| Read N50 | ~20–30 kb |
| Raw base error rate (guppy sup) | ~1–2% |
| Dominant error type | homopolymer insertions/deletions |

The `450bps` in the filename is the translocation speed — at this speed you
get longer reads on average than at 260bps but slightly higher error rate.
Guppy `sup` basecalling dramatically reduces errors compared to `fast` mode
at the cost of compute (~10x slower). The reads here were already basecalled
at `sup` quality, which is why they're usable as-is.

At 45–60x total coverage, the downstream haplotagging step assigns most reads
to one of the two haplotypes, giving you ~20–30x per haplotype — sufficient
for reliable per-haplotype methylation calls at most CpGs.

---

#### CpG methylation coverage in long-read data

Unlike bisulfite sequencing (short-read) where you need separate libraries
and struggle with conversion bias, ONT measures methylation from the raw
electrical signal for every read that passes over a CpG. This means:

- **Coverage is the same as sequencing depth** — a CpG covered by N reads
  has N methylation calls. There's no "bisulfite conversion efficiency" or
  "C→T ambiguity" problem.
- **Single-molecule resolution** — you know the methylation state of every
  CpG on every individual read, not just aggregate frequencies. This is what
  enables haplotype-resolved methylation and epiallele analysis.
- At 45–60x, >95% of autosomal CpGs should be covered by ≥5 reads.
  The mode of the coverage distribution tracks the sequencing depth.
- Repetitive and heterochromatic regions (centromeres, pericentromeric
  satellites, telomeres) are still mostly inaccessible because the reads
  can't be uniquely aligned — same problem as short reads, though long reads
  do span more repetitive elements than short reads.
- After haplotagging, per-haplotype CpG coverage is ~half the total. The
  downstream modkit and epiallele notebooks work on haplotagged reads, so
  they're operating at ~20–30x per haplotype.

---

#### Practical implication for this project

The whole point of using ONT rather than WGBS here is that you get methylation
**and** phasing information from the same molecule. Short reads can't do this:
a bisulfite-converted short read has no SNP context to phase it. ONT reads are
long enough to span a CpG of interest **and** a nearby heterozygous SNP in the
same read, which is what makes haplotype-resolved methylation possible without
needing a separate phasing experiment. That's the core experimental advantage
this entire pipeline is built on.
