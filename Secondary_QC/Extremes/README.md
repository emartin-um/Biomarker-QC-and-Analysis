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
