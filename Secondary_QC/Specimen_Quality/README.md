# Specimen Quality (Secondary QC)

**The third quality axis: properties of the patient and the tube.**

The pipeline already surfaces two technical axes — **read depth**, a site property,
and **detectability**, a blank-well property. Neither describes the specimen. This
module reports the per-specimen axis as a standing panel, so differential
collection and patient physiology across the design are visible *before* anyone
models anything.

| axis | what it is | typical structure |
|---|---|---|
| `sequencing_depth` | how deeply the well was sequenced (`log2(reads)`) | strongly site-structured (~22% between sites) |
| `cystatin_C` | CST3, a renal-function surrogate | barely site-structured (~4%) — patient physiology |
| `hemolysis_HBA1` | HBA1, the red-cell lysis marker | moderately site-structured (~14%) |
| `panel_level` | `mean_INT` — how high the sample reads across the whole panel | ~7% |

> ### 🚫 Descriptive only
> This module defines **no drop rule and no gate**. In particular **cystatin C must
> never be added as a covariate** to remove the effect: it is itself raised in
> clinical AD, so conditioning on it deletes disease signal along with the
> nuisance — MAPT loses **25%** of its effect, pTau-181 17%, pTau-217 9%.
> `log2(reads)` costs pTau-217 2% and remains mandatory. See
> `DOCS/ANALYSIS_PLAN_batch_adjustment.md` §3a.

## Running it

```bash
Rscript run_specimen_quality.R
```

Base R only — no dplyr — so it runs wherever the QC runs do. From the module
directory, after the datasets exist. Two optional inputs improve it:

- `Primary_QC/input_files/Sample_QC.csv` gives read depth (joined on
  `SAMPLE_ALIQUOT` = `Sample Name`). Without it the depth panels are skipped.
- `Secondary_QC/Extremes/output_files/extreme_sample_master.csv` gives
  `mean_INT`. Without it the panel-level axis is skipped.

Both degrade with a note rather than failing.

## Files

| file | what it is |
|---|---|
| `specimen_quality_helpers.R` | all logic — variance-between, within-group retention, co-clustering, sign balance. I/O-free. |
| `run_specimen_quality.R` | headless driver; writes the CSVs and prints a summary |
| `output_files/` | the panels below |

| output | what it answers |
|---|---|
| `specimen_axis_by_site.csv` | the distribution of each axis at each site (median, IQR, SD, n) |
| `specimen_axis_by_ancestry.csv` | the same, by ancestry |
| `specimen_axis_variance.csv` | how much of each axis lies between Site / Run / Bay / Ancestry / diagnosis |
| `specimen_axis_cocluster.csv` | do the axes single out the **same** sites? |
| `specimen_axis_within_site.csv` | does each axis survive comparing patients **within** one site? |
| `sign_balance_report.csv` | which composite indices track the sample's overall level |

## The two questions that matter

**"Is it all a site effect?"** Two outputs answer this and they should be read
together. `specimen_axis_cocluster.csv` asks whether the axes pick out the same
sites — on the 50-plate run, depth and cystatin C correlate **−0.15** across site
medians, so they are *different* sites, not one story. And
`specimen_axis_within_site.csv` asks whether the links survive within a site:
cystatin C retains **106%** of its relationship with panel level, depth **93%**.
If these were collection artifacts, comparing patients handled at the same site
would have destroyed them.

**"Which of our scores are exposed?"** `sign_balance_report.csv` computes
`Σw / Σ|w|` for every composite index the pipeline ships — 0 means the weights
self-cancel, 1 means fully exposed to a common shift. On the current definitions
**all seven sit at ±1.00**, including the hemolysis index. The generalised rule
(2026-07-27 addendum §6) is that *any* index built from same-signed markers tracks
the sample's overall level, whatever the markers are — and you can check it before
touching data.

## Where this came from

Implements the depth-by-site panel recommended on 2026-07-13 (§2, never built
until now) and the 2026-07-27 addendum §12 (cystatin C distributions) and §6 (the
sign-balance guard). The companion additions to the Extremes screen — `mean_INT`,
the global-outlier direction split, and the panel scores read against the sample's
own level — are in `Secondary_QC/Extremes/`.
