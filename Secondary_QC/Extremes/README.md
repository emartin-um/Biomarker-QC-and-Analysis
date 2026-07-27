# Extreme-Value Screen (Secondary QC)

Flags individual samples whose Alamar CNS-panel profile is **unusually extreme**
— on a single marker or across a **disease signature** — so a clinician can
review the participant more closely (possible ALS/motor-neuron disease,
Parkinson–Lewy synucleinopathy, early-onset / pre-symptomatic AD, a rare
brain-derived pTau inflation, an acute inflammatory illness, or low amyloid).

These are **relative (NPQ) protein levels, not a diagnostic assay** — every flag
is a *prompt to look*, not a diagnosis.

## Files

| File | What it is |
|------|------------|
| `extreme_helpers.R` | All logic: robust/INT scoring, disease-panel definitions, shortlist builder, one-call `run_extreme_pipeline()`. I/O-free. |
| `intensity_columns.R` | Whole-sample intensity + global-outlier direction (added 2026-07-27). Base R, file-free, so the backfill and the pipeline share one definition. |
| `backfill_intensity_columns.R` | Adds those columns to already-written outputs without re-running the screen (base R — no dplyr needed). |
| `test_intensity_columns.R` | Acceptance test for both. Run it in a QC run directory after the screen; exits non-zero on failure. |
| `run_extreme_screen.R` | Headless driver — pools U19 + OtherProjects, runs the pipeline, writes the CSVs, prints a tuning summary. |
| `Extreme_Value_Screen.Rmd` | Clinician-facing HTML report (same pipeline + plots/tables). |
| `output_files/` | CSV outputs (below). |

## How to run

```r
# from this directory
Rscript run_extreme_screen.R                       # CSVs + console summary
rmarkdown::render("Extreme_Value_Screen.Rmd",
  output_file = "Extreme_Value_Screen_<run>.html") # report + CSVs
```

Inputs are the run's split datasets `../../datasets/preprocessed_data_{U19,OtherProjects}.rds`
(list with `$data` + `$biomarker_cols`). Both are pooled into one reference
distribution and tagged by `dataset`.

## Outputs (`output_files/`)

- **`clinician_review_shortlist.csv`** — start here. Flagged samples, ranked
  HIGH→MEDIUM, with a plain-language `review_reason`, `review_caveat`,
  `priority`, `pathology_axes`, panel scores, IDs, run/bay, age, dx, APOE.
- `extreme_sample_master.csv` — every sample with extreme counts, all panel
  scores, flags (annotated; `review_reason` is NA for non-flagged).
- `extreme_marker_long.csv` — one row per (sample, marker) extreme event
  (`INT`, `robust_z`, `value`, `percentile`, `direction`, hemolysis flag).
- `extreme_reference_stats.csv` — per-marker median/MAD/SD/quantiles + caution
  flags (hemolysis-sensitive, genetic, tight-core distribution).
- `panel_top_<PANEL>.csv` — top samples per disease signature.

## Method (short)

- **Screening metric = rank-based inverse-normal (INT):** `qnorm((rank-0.5)/n)`
  per marker, so "extreme" is the same empirical rarity (`|INT| ≥ 2.75` ≈
  top/bottom 0.3%) for every marker — robust to tight cores / genotype
  multimodality that break a raw z-score. Magnitude is reported as `robust_z`
  (median/MAD) and the raw value. **Both tails** are screened.
- **Disease panels** score the mean oriented INT across signature markers and
  count markers in the tail; a panel fires when the whole signature is elevated
  **or** ≥3 of its markers spike (catches a diluted "storm").
- **Shortlist priority** is *specialness*, not just magnitude: HIGH = mixed
  pathology (≥2 distinct axes), a very strong single signature, a striking lone
  marker, an inflammation storm, or **discordance** (NCI / young carrying a
  pathology signature). Expected disease signal (e.g. AD with high tau) is MEDIUM.
- **Excluded:** `APOE4`, `APOE` (genotype-determined; APOE4 also Bay1-confounded).
  **Caveats inline:** hemolysis (HBA1 / RBC markers), global outliers (≥8 markers
  → possibly a bad sample), the `20251124-1407_Bay1` APOE4 run.

### Disease panels

`ALS_MND` (NEFL, NEFH, pTDP43_409, TARDBP) · `Synucleinopathy` (Oligo_SNCA,
pSNCA_129, DDC, PARK7) · `AD_tau` (pTau 217/181/231, BD_pTau_217, MAPT) ·
`Amyloid_low` (low A42/A40 ratio, low A42) · `BD_pTau_inflation` (BD_pTau
181/217/231, BD_MAPT) · `Neuroinflammation` (CRP, IL6, TNF, IL1B, SAA1, CXCL10,
GDF15, IL18, CHI3L1) · `Neurodegeneration` (GFAP, NEFL, NEFH, GDF15).

Thresholds are parameters (`thr_int`, `panel_score_thr`, …) in the driver and the
Rmd `params:` block — change them to widen or tighten the net.


---

## Whole-sample intensity, and the direction of a global outlier (added 2026-07-27)

Two columns groups the screen was missing. Full rationale:
`Secondary_QC/Batch_Effects/2026_Batch_Investigation/QC_PIPELINE_RECOMMENDATIONS_ADDENDUM_2026-07-27.md`
§9–§10 (in the QC run directory).

### `mean_INT` — the level question

The screen builds the full INT matrix and then only ever asks **count** questions
of it (`n_extreme_high`, `n_extreme_low`, `max_abs_INT`). It never asks how high
the sample sits **overall**. So a sample that is moderately high on *all* markers
at once is invisible to it: mean `max_abs_INT` is 2.50 in the middle decile of
background intensity and only 2.96 in the brightest, because no single marker
ever reaches the tail.

| column | what it is |
|---|---|
| `mean_INT` | mean INT across all usable markers — "how high does this sample read across the panel" |
| `mean_INT_z` | the same, standardised across samples |
| `intensity_flag` | `HIGH` / `LOW` beyond ±2.5 SD, else empty |

It correlates **r = 0.952** with the per-sample background factor derived
independently in `AD_Prediction_MCI_OD`, and it is **additive**: of the 52 samples
it flags on the 50-plate run, **31 are not global outliers**.

> ⚠️ **Informational, not a drop rule.** Background intensity is a continuous
> nuisance covariate affecting every sample, not a population of bad wells. No
> outlier screen recovers it — the best whole-panel statistic at a +2 SD bar
> flags ~2% of wells. The axis belongs in the analysis model. These columns exist
> to make it *visible*.

### `global_outlier_dir` — which tail

`global_outlier = n_extreme_total >= 8` merges two opposite populations, because
its components move in opposite directions with intensity (across background
deciles `n_extreme_high` runs 0.09 → 1.94 while `n_extreme_low` runs 1.08 → 0.18).
On the 50-plate run the 55 outliers split **26 high-driven / 28 low-driven /
1 mixed**, sitting at `mean_INT_z` **+1.85 vs −2.30** — a bright sample and a dim
sample carrying one label.

| column | what it is |
|---|---|
| `global_outlier_dir` | `HIGH` / `LOW` / `MIXED` for global outliers, empty otherwise |
| `global_outlier_high`, `global_outlier_low` | booleans for filtering |

**`global_outlier` itself is unchanged**, so existing consumers are unaffected.
Note it gates on the **total** — a sample with 5 high and 4 low is an outlier with
neither tail reaching 8 — so the split is "which tail dominates", *not* a union of
two bar tests. Implementing it as a union would silently drop those samples.

### Consuming these downstream

Join **`extreme_sample_master.csv`** (one row per sample), never
`clinician_review_shortlist.csv` (only samples selected for review). A left-join
against the shortlist silently turns "not reviewed" into "clean" — that bug cost
a downstream eligibility gate 10 of its 28 global outliers in July 2026.
