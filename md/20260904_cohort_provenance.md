# Cohort provenance and selection basis

**Written 2026-09-04.** Reconstructed from `csv/hprc/hprc_r2_ont.csv` and direct
`s3://human-pangenomics` LIST calls, because no selection script, notebook, or note recorded
how the 30-donor cohort was built. Everything below is re-derivable with
`scripts/W15_survey_rebasecalled.py`.

This exists so the Methods section rests on a documented basis rather than a later
reconstruction, and so the 18-donor decision can be defended on the correct grounds.

---

## 1. Source data

| | |
|---|---|
| Bucket | `s3://human-pangenomics/working/HPRC/{sample}/raw_data/nanopore/` |
| Access | public; no credentials (`aws s3 cp --no-sign-request`, or plain HTTPS) |
| Local manifest | `csv/hprc/hprc_r2_ont.csv` — 232 samples, 929 ONT file records |
| Genotypes | 1000G 2504 high-coverage phased panel (all 30 rows of `cohort_downloads.tsv`) |
| Assemblies | `csv/hprc/hprc_r2_assemblies.csv` |

**`hprc_r2_ont.csv` is a point-in-time snapshot and is stale.** It reports the 12 excluded
donors as R9.4.1-only; they now also have R10.4.1 data (§5). Never use it alone to decide
cohort membership — query S3.

**There is no public POD5 or FAST5 in this bucket** (checked across all seven top-level
prefixes, 2026-09-04). Re-basecalling from raw signal is not possible from public data.
HPRC retains signal internally — visible in the `@PG` line of published BAMs, which
reference `Year4_s3_pod5.list_*` paths — but does not distribute it.

Note for Methods: `working/` is HPRC's pre-release area, not a frozen release. The file
list can change. Cite as "HPRC working data, accessed 2026-09-04."

## 2. The pool: R10.4.1 availability is strongly confounded with ancestry

Of 232 samples with ONT data, 73 have R10.4.1 and 159 are R9.4.1-only. The split is not
uniform across ancestry:

| superpop | R10.4.1 | R9.4.1 only | total | % R10.4.1 |
|---|---|---|---|---|
| EUR | 21 | 9 | 30 | **70%** |
| EAS | 18 | 32 | 50 | 36% |
| AFR | 20 | 44 | 64 | 31% |
| SAS | 9 | 27 | 36 | 25% |
| AMR | 3 | 40 | 43 | **7%** |
| **all** | **73** | **159** | **232** | **31%** |

**This single table explains the cohort's shape.** Any design balanced on ancestry, drawn
from the full pool without chemistry as a criterion, is automatically unbalanced on
chemistry — because the two are correlated in the source data.

## 3. Why the original cohort is 30

30 = 6 per superpopulation × 5, with 15F/15M. That is a **balanced-design cap, not an
availability ceiling.** Two lines of evidence:

- Exactly 6/6/6/6/6 and exactly 15F/15M is a deliberate draw, not what happened to exist.
- The cohort contains **6 AMR donors, but only 3 AMR samples in all of HPRC have R10.4.1**.
  Six AMR could not have been drawn from an R10-filtered pool. Selection therefore ran over
  the full 232-sample pool with chemistry not applied as a filter.

Ruled out as the binding constraint (both checked against the 56-sample pool of §6):
- assembly availability — 56/56 have an R2 assembly;
- ONT coverage — 56/56 ≥ 40×, 55/56 ≥ 50×, median 77×.

Residual uncertainty: if the original search predates mid-2024 or used a Release 1
manifest, some samples were not yet public. 18 of the 56 carry 2024 run dates, but the
majority date from 2023, so this cannot be the main explanation.

## 4. Why 18, and why that is a chemistry decision

The 18/12 split within the cohort is **perfectly confounded with flow cell chemistry**:

| | chemistry | basecaller | n |
|---|---|---|---|
| included | R10.4.1 | dorado (17 × 0.6.0, 1 × 0.5.3) | 18 |
| excluded | R9.4.1 | guppy 6.5.7 | 12 |

R9.4.1 and R10.4.1 are different pores — single vs dual constriction, ~5–6 nt vs ~9–10 nt
sensing region. **No basecaller converts one into the other.** Re-basecalling the R9 samples
with dorado would harmonise the software and leave the physics untouched.

This matters specifically for this project, beyond generic batch effect: §6.2 of
`20260903_qc_review.md` concerns sequence-context calibration of the modification caller.
The context window differs between chemistries *by construction*, so pooling R9 and R10 in
an allelic analysis would confound the exact artifact that section exists to control.

**The exclusion is therefore correct science, not convenience, and should be described that
way.** The "purely for convenience/simplicity" framing in earlier notes understates it and
invites a reviewer question that has a good answer.

## 5. What the 18 costs, and what changed

Cost — the exclusion is not ancestry-neutral:

| | AFR | AMR | EAS | EUR | SAS |
|---|---|---|---|---|---|
| full 30 | 6 | 6 | 6 | 6 | 6 |
| **included 18** | **2** | 3 | 3 | 5 | 5 |
| excluded 12 | 4 | 3 | 3 | 1 | 1 |

AFR drops to 2. AFR genomes carry the most heterozygosity, hence the most phasable het SNPs
and the highest yield of ascertainable ASM loci per donor — so the loss falls on the
best-powered samples. **With n=2 AFR, no ancestry-stratified claim is supportable**, and
`doc/intro.md` must not claim one (see §8).

Changed — between Aug and Dec 2024, HPRC published **new R10.4.1 ultra-long runs with
`5mCG_5hmCG` calls for all 12 excluded donors** (e.g.
`08_26_24_R1041_UL_HPRC_HG01891_1_dorado1.0.0_sup5.0.0_5mCG_5hmCG.bam`; 12/12, ~1.7 TB).
These are at sup5.0.0, not the sup4.3.0 used by the current 18, so adopting them would
reintroduce a model-version mismatch — a milder form of the problem the exclusion avoided.

## 6. The available pool at the current model

Filtering the bucket for R10.4.1 **and** `5mCG_5hmCG` mods **and** sup4.3.0 — i.e. exactly
what the current 18 already are — yields **56 samples**:

| superpop | in pool | in cohort | available to add |
|---|---|---|---|
| AFR | 15 | 2 | 13 |
| AMR | 3 | 3 | **0** |
| EAS | 8 | 2 | 6 |
| EUR | 21 | 5 | 16 |
| SAS | 9 | 5 | 4 |
| **total** | **56** | **17** | **39** |

(17, not 18: NA18959 has no data under `working/HPRC/` — see §9.)

A balanced expansion capped by availability gives **38 donors** — AFR 9, AMR 3, EAS 8,
EUR 9, SAS 9 — all R10.4.1 UL, all dorado sup4.3.0 `5mCG_5hmCG`. One chemistry, one model,
no harmonisation argument required, and **no basecalling**: cost is transfer plus alignment.
21 new samples, ~4.9 TB.

**AMR is hard-capped at 3.** All three are already in the cohort. This is a property of
HPRC's sequencing rollout (§2), not of any choice made here, and should be stated as a
limitation rather than designed around.

Artifacts:
- `csv/W15_survey_all_r10.tsv` — full survey, 133 basecall directories over 73 samples
- `csv/W15_proposed_expansion.tsv` — the 21 proposed additions
- `csv/W15_download_urls.tsv` — 223 BAM files, 56 samples, 12.4 TB, with `s3://` and
  `https://` per file (154 files / 8.6 TB are outside the current cohort)
- `scripts/W16_fetch_expansion.sh` — download → verify mod tags → align → merge → index

Scaling later is strictly additive. Same bucket, same chemistry, same model, so nothing
computed on the 18 is invalidated by adding donors.

## 7. Verification traps

1. **`sup5.0.0/` vs `sup5.0.0_5mCG_5hmCG/`.** HPRC publishes both under near-identical
   names. Verified directly against BAM records (2026-09-04): the suffixed directory gives
   `dorado basecaller "sup@v5.0.0,5mCG_5hmCG"` with `MM:Z`/`ML:B` on 200/200 reads; the
   un-suffixed gives `dorado basecaller "sup@v5.0.0"` with **zero** mod tags. A no-mods BAM
   runs the full downstream pipeline to empty methylation output **without erroring
   anywhere**. `W16` aborts unless ≥900/1000 reads carry `MM:Z`.
2. **`minimap2 -y` is mandatory.** Without it, MM/ML tags are dropped at alignment and every
   methylation call is silently lost.
3. **Reads are ultra-long.** The manifest reports read N50 74,092 bp for a NA21144 run;
   `SQK-ULK114`, `_UL_` prep. Quote the measured N50 from the full BAM, not an adjective.
4. **The autosome-only reference.** Alignment used `GRCh38.autosome.fa` (22 contigs, no
   chrX/Y/M), which blocks the chrX/XCI control in `W10_control_loci.sh`.

## 8. Consequences for the manuscript

- `doc/intro.md` currently says "18 **globally diverse** donors" and claims this is "the
  first study that leverages long-read data to probe patterns of ASM **across diverse
  ancestries**." **Both are unsupportable at AFR n=2** and must be softened, or the cohort
  expanded first.
- The exclusion rationale in Methods should be stated as chemistry (§4), not convenience.
- The AMR cap of 3 belongs in Limitations.

## 9. Open items

- **NA18959** — in the current 18, `status=new`, basecalled with Dorado 0.5.3 (the only
  non-0.6.0 sample), and **has no directory under `working/HPRC/`**. If it was never
  downloaded, the working cohort is 17. Resolve before analysis, not after.
- Coverage: the 56-sample pool has median 77×, min 43×; the current Dorado-18 has median
  78×, min 51×. Per-haplotype depth ~35–40×, which sets the entropy ceiling discussed in
  `20260903_qc_review.md` §7.3 — those numbers were written against ~30× and should be
  restated.
- `csv/cohort_downloads.tsv` (30 rows) is current. `md/cohort_downloads.tsv` and
  `csv/old/cohort_downloads.tsv` are stale 10-row versions overlapping the current cohort by
  6 and 8 samples. Delete or date-prefix them; they have already caused one incorrect
  conclusion in this project.
