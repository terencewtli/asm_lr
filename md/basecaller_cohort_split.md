# Basecaller cohort split — Dorado vs Guppy

**Decision**: the 30-donor HPRC cohort splits into two basecaller batches. Only
the 18 Dorado/R10.4.1 donors are the primary whole-genome replication cohort.
The 12 Guppy/R9.4.1 donors (including the HG01258 pilot) are intentionally
held out of that cohort and are not to be pooled with it.

## Cohort composition

Source: `csv/meta/hprc_r2_ont.csv`, `basecaller` / `basecaller_model` columns.

| Batch | N | Basecaller | Chemistry | Modification model |
|---|---|---|---|---|
| Primary (WG replication) | 18 | Dorado 0.6.0 (one sample: 0.5.3) | R10.4.1 | `5mCG_5hmCG` — joint 5mC+5hmC calling |
| Held out | 12 | Guppy 6.5.7 | R9.4.1 | `5mc_cg` — 5mC only |

- Primary list: `txt/samples/final_samples.txt`
- Held-out list: `txt/samples/excluded_samples.txt`
- Of the 12 held out, `HG01258` is the phase-1 chr19 pilot (`old/`, `scripts/HG01258/`,
  publication-ready at chr19 scope). The other 11 have raw ONT data downloaded
  but no alignment/phasing has been run.

## Why the split

Two independent differences between the batches both bear directly on
methylation-based ASM calling, which is the core measurement of this project:

1. **Different modification-calling models, not just software versions.**
   Guppy 6.5.7 here was run with a 5mCG-only config
   (`dna_r9.4.1_450bps_modbases_5mc_cg_sup_prom`). Dorado 0.6.0 uses a joint
   `5mCG_5hmCG` model. If "modified C" is ever treated as one class downstream,
   Dorado samples carry 5hmC signal that is structurally invisible in the
   Guppy samples — a systematic per-sample confound rather than noise.
   `notebooks/qc/QC01_cohort_overview.ipynb` (cell 9) already flags this via
   `has_5hmcg`.
2. **Different pore chemistry.** R10.4.1 has materially higher raw per-read
   accuracy than R9.4.1, especially in homopolymers (R10.4 modal accuracy
   >99.1%; R9.4.1 lower, with worse resolution on 4–9bp homopolymers) —
   [Nakagawa et al., benchmarking R10.4 vs R9.4.1 flow cells](https://pmc.ncbi.nlm.nih.gov/articles/PMC10070092/).
   Base-calling accuracy is not independent of modification-calling accuracy,
   so this compounds point 1 rather than being a separate, ignorable effect.
   A recent review of ONT methylation callers explicitly cautions that
   concordance figures from R9.4.1/Guppy pipelines are not transferable to
   R10.4.1/Dorado pipelines — [Coverage-aware evaluation of ONT methylation
   callers](https://www.sciencedirect.com/science/article/pii/S135503062600095X).

Net: pooling the two batches into one "cohort" would let basecaller identity
masquerade as biological signal in any cross-donor methylation comparison.

## History (for anyone re-deriving this)

- **2026-06-12**: `txt/samples/excluded_samples.txt` / `final_samples.txt`
  created, splitting exactly along the Dorado/Guppy line above. No comment was
  written explaining the split at the time — this doc backfills that.
- **2026-08-21**: a restart of the project (via a fresh Claude Code session
  that didn't have this context) generated `W00c_align_merge_remaining.sh`,
  `W01i_phase_remaining.sh`, `W01j_merge_haplotagged_remaining.sh`, and
  `txt/samples/remaining_samples.txt`, intending to bring 11 of the 12
  held-out donors (all Guppy) back into the primary whole-genome pipeline.
  The stated rationale in those scripts was purely operational ("never got
  past raw-data download"), with no mention of basecaller/chemistry — i.e.
  the exclusion was reintroduced as an accident of not having this context,
  not a reasoned reversal. None of those jobs were ever submitted (no logs
  exist), and **the scripts were deleted on 2026-09-01** before they could
  run and silently repool the two batches.

## What this means going forward

- Treat the 18-donor Dorado/R10.4.1 cohort as the whole-genome replication
  cohort for all pooled/population-level ASM and epiallele analysis.
- The 12 Guppy/R9.4.1 donors are out of scope for that analysis. If they're
  ever analyzed, they should be a separate, clearly-labeled stratum — not
  merged rows in the same table — and any modification-calling comparison
  across the two would need model harmonization (e.g. collapsing Dorado's
  5hmC calls into 5mC, or restricting to a shared 5mC-only representation)
  before combining.
- Raw ONT uBAMs for all 12 Guppy donors are downloaded (~1.9TB total in
  `bam/ubam/{sample}/`). To make the exclusion structurally obvious as early
  as possible in the pipeline (rather than relying on a sample list a future
  session might not read), the raw data for the 11 non-pilot Guppy donors can
  be deleted — see terminal output from 2026-09-01 for the exact `rm` command,
  run manually by the user. `HG01258`'s raw data should stay, since it backs
  the existing chr19 pilot results.
