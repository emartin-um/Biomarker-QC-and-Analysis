# Batch Effects Analysis

⚠️ **Written March 2026 — not yet tested on real data. Run and validate before using results.**

Assesses Run- and Bay-level technical batch effects in the merged biomarker dataset.
Runs after `Metadata_Merge` — requires `merged_combined_post_QC.csv`.

## Targeted investigations

- **`investigate_20251124_Bay1.R`** — deep-dive on the APOE4 assay anomaly on run
  `20251124-1407_Bay1` (first surfaced by the APOE secondary QC). Finding: the **APOE4
  target alone** is broken on ~**32% of Bay1 wells** (scrambled, not signal-loss;
  homozygous ε4 crash to ~0), total APOE and the other ~130 proteins intact; not a sample
  sub-batch or ancestry effect. **Recommendation: drop APOE4 / APOE4−APOE for the whole
  bay.** Outputs in `output_files/` (`Bay1_APOE4_anomaly_findings.md`,
  `bay1_apoe4_sample_deviations.csv`, `bay1_apoe4_anomaly.png`,
  `bay1_panel_iqr_ratio.txt`). Genotypes are unaffected — unrelated to the WGS plate swap.

## Usage

1. Open `Batch_Effects_Analysis.Rmd` in RStudio (create `.Rproj` via File → New Project → Existing Directory if needed)
2. Input auto-copied from most recent `Metadata_Merge/output_files_*/merged_combined_post_QC.csv`
3. Knit; outputs saved to `output_files/`
4. Add dataset-specific observations to `notes.md` (gitignored)

## Input File

**Default:** `merged_combined_post_QC.csv` — has separate `Run` and `Bay` columns, includes all biomarkers (standard + low-abundance).

**Fallback:** `merged_standard_post_QC.csv` — standard biomarkers only; `Run`/`Bay` parsed from combined `RUN` column.

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `VP_CORES` | 4 | Cores for variancePartition |
| `VP_MIN_SAMPLES` | 10 | Min samples per Run to include in VP |
| `ICC_MIN_REPS` | 2 | Min replicate pairs for ICC |
| `FDR_THRESHOLD` | 0.05 | FDR cutoff for batch effect significance |

## Sections

| Section | What it does |
|---------|-------------|
| C. Plate Composition | Race/sex/age/site/Dx distribution per Run and Bay |
| C.7 Confounding | Cramér's V between batch variables and covariates |
| D. Variance Partitioning | `variancePartition` with Run and Bay as random effects |
| E. ICC Analysis | Within-Run and cross-Run ICC from HIHG replicates |
| F. Formal Test | Per-protein Type III F-test for Run and Bay terms |
| G. Summary | Automated recommendation for correction approach |

## Output Files

| File | Contents |
|------|----------|
| `batch_covariate_confounding.csv` | Cramér's V for each batch × covariate pair |
| `variance_partition_results.csv` | Per-protein % variance for each source |
| `icc_results.csv` | Per-protein ICC within-Run and cross-Run |
| `batch_effect_test_results.csv` | p and FDR for Run and Bay terms per protein |

## Dependencies

```r
install.packages(c("tidyverse", "lme4", "irr", "kableExtra", "patchwork", "car"))
BiocManager::install(c("variancePartition", "BiocParallel"))
```

## Workflow Position

```
Metadata_Merge
    └── output_files_*/merged_combined_post_QC.csv
            │
            ▼ (auto-copied)
Batch_Effects/input_files/
            │
            ▼ Batch_Effects_Analysis.Rmd
Batch_Effects/output_files/
    ├── batch_covariate_confounding.csv
    ├── variance_partition_results.csv
    ├── icc_results.csv
    └── batch_effect_test_results.csv
```

## Decision Framework

| Evidence | Recommendation |
|----------|---------------|
| No significant proteins, median Run var < 5% | Include Run as fixed covariate only |
| Moderate effects, no high confounding | Run as random effect (lme4) or ComBat |
| High confounding (Cramér's V > 0.5) | Random effects model preferred; ComBat with caution |
| Substantial effects | ComBat with biological covariates as protected |
