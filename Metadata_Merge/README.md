# Metadata Merge

Merges post-QC biomarker data with sample metadata, applies exclusion filters, derives covariate groupings (diagnosis, ancestry, APOE status with WGS correction), and produces analysis-ready datasets including a `preprocessed_data.rds` for downstream pipelines.

## Usage

### Recommended: Run the Pipeline in RStudio

Open `Metadata_Merge_Pipeline.qmd` in RStudio and click **Render** (or press Ctrl/Cmd+Shift+K).

The pipeline runs in 9 steps and generates an HTML report with all decisions documented.

### Render Parameters

Control pipeline behavior via the YAML `params:` block or at render time:

| Parameter | Default | Description |
|---|---|---|
| `filter_only` | `false` | Skip merge/QC; re-apply filters to existing merged files |
| `no_copy` | `false` | Skip auto-copy of NPQ files from Primary_QC |
| `include_derived_covariates` | `true` | Run Step 6 (CDX_collapsed, Ancestry, APOE4 variables) |
| `include_derived_biomarkers` | `false` | Run Step 7 (APOE4_minus_APOE, ABeta42_minus_ABeta40) |
| `filter_na_covars` | `true` | Exclude rows where CDX_collapsed or Ancestry is NA after grouping |
| `output_dir` | `"output_files"` | Output directory for all pipeline files (created if absent) |
| `review_dir` | `"review"` | Dir holding `CDX_review.csv` / `Ancestry_review.csv`. If **both** exist, Step 6b is **data-driven** from their `FINAL_*` columns; otherwise the built-in config is used |

To override at render time:
```r
quarto::quarto_render("Metadata_Merge_Pipeline.qmd",
                      execute_params = list(filter_only = TRUE,
                                            include_derived_biomarkers = TRUE))
```

### Filter-Only Mode

To re-run filtering after manually editing the exclusion report:
```r
quarto::quarto_render("Metadata_Merge_Pipeline.qmd",
                      execute_params = list(filter_only = TRUE))
```

---

## Scripts

### `Metadata_Merge_Pipeline.qmd`
Main pipeline document. See **Pipeline Steps** below.

### `combine_metadata.R`
Standalone helper script for combining multiple metadata source files before a pipeline run. Reads source CSVs/XLSXs from `Datasets/`, adds a `metadata_source` column to each, combines with `bind_rows()`, and writes the result to `Datasets/`. Run this from the `Metadata_Merge/` directory in RStudio, then copy the output into `input_files/` before rendering the pipeline. Do not overwrite the source files.

### `covariate_explorer.R`
Utility library sourced automatically by the pipeline. Provides config-driven functions for exploring covariate distributions and applying grouping decisions. Key functions:

| Function | Purpose |
|---|---|
| `explore_covariate(data, var)` | Frequency table for one variable |
| `explore_combinations(data, ...)` | All unique combinations of multiple variables |
| `make_category_config(...)` | Define a single-column recoding rule |
| `make_grouping_config(...)` | Define a multi-column grouping rule |
| `summarize_category_decisions(data, config)` | Preview what a category config will do (no data change) |
| `summarize_grouping_decisions(data, config)` | Preview what a grouping config will do (no data change) |
| `apply_category_config(data, config)` | Apply recoding and drop excluded rows |
| `apply_grouping_config(data, config)` | Apply multi-column grouping and optionally drop unmatched rows |
| `validate_groups(data, cols, min_n)` | Check final group sizes against a minimum threshold |
| `print_config(config)` | Pretty-print a config for review |

---

## Pipeline Steps

### Step 1 — Input Files
Verifies NPQ files exist in `input_files/`; auto-copies from `Primary_QC/output_files/` if absent.

### Step 2 — Merge Data
Merges biomarker data with metadata. `APOE.geno` is normalized at this step from raw format (e.g., `"E3    E4"`) to standard `"a/b"` digit format (e.g., `"3/4"`).

### Step 3 — Quality Control
Identifies: duplicate samples (keeps most recent run), implausible values (age, BMI, weight, height), illogical Race/Ethnicity combinations (flagged for review, not auto-excluded).

### Step 4 — Apply Exclusion Filters
Removes auto-excluded samples; writes preliminary filtered CSVs. Exclusion report can be manually edited before re-running in filter-only mode. Note: `filtered_standard_post_QC.csv` and `filtered_low_post_QC.csv` are overwritten again at Step 6d after covariate-driven exclusions are applied.

### Step 5 — Exclusion Summary
Reports sample counts before/after filtering.

### Step 6 — Covariate Groupings *(optional: `include_derived_covariates`)*

Config-driven, designed to be re-configured for each new dataset. Four sub-steps:

- **6a Explore**: Normalizes the `Group` column (actual `NA` → string `"NA"`) so ancestry rules match datasets that lack a Group assignment, then prints CDX distribution and Group × Race × Ethnicity combinations to inform grouping decisions.
- **6b Config**: Defines `dx_config` (diagnosis grouping) and `ancestry_config` (ancestry grouping). **Data-driven (2026-06):** if `<review_dir>/CDX_review.csv` and `Ancestry_review.csv` exist, the configs are built from their `FINAL_*` columns (edit those CSVs to regroup — blank = exclude); otherwise the in-document built-in block is used (**edit that block for a new dataset**). See *Adapting for a New Dataset*.
- **6c Preview**: Shows decision tables and any unmatched rows *before* applying anything — no data is modified.
- **6d Apply**: Applies configs, joins WGS APOE genotypes, computes APOE4 carrier variables, validates group sizes, saves updated CSV. Also re-filters `filtered_standard_post_QC.csv` and `filtered_low_post_QC.csv` to match the final sample set (overwriting the Step 4 versions), so all three filtered files share the same sample universe. **The WGS join reads `APOE_WGS` from the auto-detected metadata** (the same file used for the merge); **no WGS plate-swap correction is applied** — the swap-affected AD-Hispanic samples are masked to Sanger via `../../../May_2026_WGS_QC/WGS_chr1_corrected_idmap_2026Jun11.csv` (see the *WGS plate-swap* note below); `APOE.geno_final = coalesce(APOE_WGS_norm, APOE.geno)`.

Derived columns produced:

| Column | Type | Description |
|---|---|---|
| `CDX_collapsed` | character | Collapsed diagnosis: NCI / MCI / AD / Dementia_Other |
| `Ancestry` | character | Grouped ancestry: AA / AFDC / HI_WH / HI_MU / HI_BL |
| `APOE_WGS` | character | Raw WGS genotype (e.g., `"34"`); NA if not genotyped |
| `APOE_WGS_norm` | character | WGS genotype normalized to `"a/b"` format (e.g., `"3/4"`) |
| `APOE_WGS_changed` | logical | TRUE if WGS genotype differs from original `APOE.geno` (**12** samples after the plate-swap correction; was 137 before) |
| `APOE.geno_final` | character | Best available genotype: WGS where available, otherwise original |
| `APOE4_carrier` | logical | TRUE if any E4 allele in `APOE.geno_final` |
| `APOE4_count` | integer | Number of E4 alleles in `APOE.geno_final` (0 / 1 / 2) |

> **Downstream analyses should use `APOE.geno_final`** rather than `APOE.geno` when genotype accuracy matters. The original `APOE.geno` is retained for comparison with prior results.

#### WGS plate-swap — correction REVERTED (updated 2026-06-12)

The May-2026 WGS delivery had a physical plate mix-up on the AD-Hispanic set that attached the wrong WGS APOE genotype to ~282 Hispanic samples. A **2026-06-08 geometric correction** (chrX fingerprint; `apply_wgs_plate_swap_correction.R` → `..._plateswapcorrected_2026Jun.csv` + audit log) was later found **incorrect** and has been **reverted**. The more complete chromosome-1 analysis (`../../../May_2026_WGS_QC/WGS_Plate_Swap_Summary_2026Jun11.md`) shows a different geometry and identifies the **WGS itself** as the mislabeled platform; the candidate relabel is **not yet applied** (pending USUHS confirmation).

- The two June-8 files are archived in `Datasets/archive/` — **do not use.**
- **Step 6d applies no correction.** It reads `APOE_WGS` from the auto-detected metadata, then **blanks it for the swap-affected samples** (`WGS_chr1_corrected_idmap_2026Jun11.csv$SAMPLE_aslabeled`) so `APOE.geno_final` falls back to **Sanger** for them.
- Treat `APOE.geno_final` of the affected AD-Hispanic samples as **provisional** until the chr1 relabel is confirmed and applied.

> **⚠ Downstream propagation.** `APOE_WGS` / `APOE.geno_final` / `APOE4_*` for the affected samples differ from any output built on the (now-reverted) 2026-06-08 plate-swap-corrected data. Re-run downstream pipelines (e.g. `APOE_Amyloid_Tau/`, `MCI_Progression_Analysis/`) against current outputs. Earlier outputs are **preserved, not deleted** — see *Report / output provenance*.

When `filter_na_covars = true` (default), rows where CDX_collapsed or Ancestry is NA are excluded.

### Step 7 — Derived Biomarkers *(optional: `include_derived_biomarkers`, default FALSE)*

Adds log2 difference columns: `APOE4_minus_APOE` and `ABeta42_minus_ABeta40`.

### Step 8 — Longitudinal Handling (Cross-Sectional vs Longitudinal Datasets)

Resolves **subject-level longitudinal repeats** — the same subject (`Record_ID`) sampled at
different visits (different `SAMPLE`). (This is distinct from the Step 3 SAMPLE-level dedup, which
collapses *same-sample reruns* — the same physical sample assayed in multiple runs, identical
metadata — keeping the most recent run.)

Adds four per-subject flag columns, then writes **two datasets** so each analysis uses the right one:

| Column | Type | Description |
|---|---|---|
| `n_visits` | integer | Number of visits (distinct samples) for this subject |
| `is_repeat_subject` | logical | TRUE if `n_visits > 1` |
| `subject_visit_n` | integer | Visit index, 1 = earliest by `age_at_subject` |
| `is_most_recent` | logical | TRUE for the latest visit (the cross-sectional pick) |

| Output | Contents | Use for |
|---|---|---|
| `filtered_combined_post_QC.csv` | **DEFAULT — cross-sectional**: one row per subject (the most recent visit) | All cross-sectional analyses (no double-counting) |
| `filtered_combined_longitudinal.csv` | **Longitudinal-only**: every visit of subjects with ≥2 visits | Trajectory / conversion analyses |

The report prints a longitudinal-subjects table with a neutral `cdx_changed` flag (diagnosis change
across visits is surfaced, not treated as an error). `preprocessed_data.rds` (Step 9) is built from
the cross-sectional default.

### Step 9 — Save RDS

Writes `preprocessed_data.rds` to `output_dir` (default: `output_files/`):
```r
list(
  data           = filtered_combined,   # cross-sectional default (one row per subject) after Step 8
  biomarker_cols = biomarker_cols       # character vector of NPQ biomarker column names
)
```
This RDS is the primary input for downstream analysis pipelines (e.g., APOE_Amyloid_Tau). For
longitudinal analyses, read `filtered/filtered_combined_longitudinal.csv` instead.

---

## Input Files

Place in `input_files/` (not tracked in git):

| File | Description | Source |
|---|---|---|
| `NPQ_post_QC.csv` | Standard post-QC biomarker data | Primary_QC output |
| `NPQ_Low_post_QC.csv` | Low detectability biomarker data | Primary_QC output |
| `*metadata*.csv` | Sample metadata (auto-detected: most recently modified file with `metadata` in the name). Current: `U19_otherprojects_metadata_combined_APOE_WGS_n2871_plateswapcorrected_2026Jun.csv` (also the WGS-join source in Step 6d). Superseded files live in `input_files/old/`. | External / `Datasets/` |

---

## Output Files

Generated in `output_dir` (default: `output_files/`, not tracked in git):

| File | Description |
|---|---|
| `merged_standard_post_QC.csv` | Merged standard biomarkers + metadata |
| `merged_low_post_QC.csv` | Merged low biomarkers + metadata |
| `merged_combined_post_QC.csv` | All biomarkers + metadata |
| `sample_exclusion_report.csv` | Exclusion flags and reasons per sample |
| `samples_not_in_metadata.csv` | NPQ samples missing from metadata (if any) |
| `filtered/filtered_standard_post_QC.csv` | Filtered standard biomarkers (all sample rows) |
| `filtered/filtered_low_post_QC.csv` | Filtered low biomarkers (all sample rows) |
| `filtered/filtered_combined_post_QC.csv` | **DEFAULT cross-sectional** — combined biomarkers + derived/flag columns, **one row per subject** (most recent visit) |
| `filtered/filtered_combined_longitudinal.csv` | **Longitudinal-only** — every visit of subjects with ≥2 visits |
| `preprocessed_data.rds` | R list: `$data` (cross-sectional default) + `$biomarker_cols` — primary downstream input |

### QC / Error reports & mapping keys (`reports/`)

A dedicated QC-reporting section writes per-issue reports (with metadata context, for sending to the contributing groups) plus shareable mapping keys to `<output_dir>/reports/`. Checks run on the full batch present in metadata (so issues show even for samples later excluded). Each report is written as an **`.xlsx`** so APOE genotypes (`3/4`, `2/3`, …) stay **text** and aren't auto-converted to dates when opened in Excel (uses the `writexl` package; falls back to `.csv` if it isn't installed).

| File | Contents |
|---|---|
| `mapping_key_CDX.csv` | CDX raw → `CDX_collapsed` (+ n, kept/EXCLUDED) — share for review |
| `mapping_key_Ancestry.csv` | Group × Race × Ethnicity → Ancestry (+ n, `unexpected_race_eth`) — share for review |
| `error_implausible_values.csv` | out-of-range age/BMI/weight/height + metadata |
| `error_unexpected_race_ethnicity.csv` | flagged combos (formerly "illogical") + metadata + mapped Ancestry |
| `missing_data_by_sample.csv` | samples missing **critical** fields (CDX, Race, Ethnicity, sex, age, APOE-any) |
| `missing_data_summary.csv` | n + % missing per field, incl. APOE three ways (Sanger / WGS / any) |
| `error_APOE_integrity.csv` | Sanger-vs-WGS discordance + invalid genotype format |
| `error_clinical_consistency.csv` | AOO > age, Years_Onset, CDX vs Case_Control conflicts |
| `check_CDX_vs_CaseControl_crosstab.csv` | CDX × Case_Control cross-tab for eyeballing |
| `error_categorical_unexpected_levels.csv` | unexpected sex / Ethnicity / Race values |
| `error_duplicate_SAMPLE_ALIQUOT.csv` | non-unique aliquots (should be empty) |

NPQ samples lacking metadata are **not** duplicated into `reports/` — they remain in `<output_dir>/samples_not_in_metadata.csv` (written in Step 2) and are referenced from the closing section. `sample_exclusion_report.csv` is **internal** (drives filtering) and is not a review file.

The report **ends with a "Files to Send for Review & Verification" section** that lists, with ⚠/✓ flags and row counts, exactly which files to send out (mapping keys + non-empty error reports + `samples_not_in_metadata.csv`) versus keep internal.

> Terminology: the Step-3 Race/Ethnicity flag is **"unexpected"** (was "illogical"); the `Ancestry_review.csv` flag column is `unexpected_race_eth`. Step 6c prints an explicit **Group × Race × Ethnicity → Ancestry** mapping table alongside the CDX one.

---

## Workflow

```
Primary_QC/output_files/
    ├── NPQ_[date]_post_QC.csv
    └── NPQ_[date]_Low_post_QC.csv
            │
            ▼ (auto-copied if input_files/ is empty)

Metadata_Merge/input_files/
    ├── NPQ_post_QC.csv
    ├── NPQ_Low_post_QC.csv
    └── *metadata*.csv  (auto-detected)
            │
            ▼ Metadata_Merge_Pipeline.qmd (Steps 1–5)

Metadata_Merge/{output_dir}/            ← default: output_files/
    ├── merged_*.csv
    └── sample_exclusion_report.csv  ← Review/edit if needed
            │
            ▼ Step 6 (covariate groupings) + Step 7 (derived biomarkers)

Metadata_Merge/{output_dir}/filtered/
    └── (combined + CDX_collapsed, Ancestry, APOE.geno_final, APOE_WGS_changed, APOE4_carrier, APOE4_count)
            │
            ▼ Step 8 (longitudinal flags + split)

Metadata_Merge/{output_dir}/filtered/
    ├── filtered_combined_post_QC.csv       ← DEFAULT cross-sectional (one row/subject, most recent)
    └── filtered_combined_longitudinal.csv  ← longitudinal-only (≥2-visit subjects, all visits)
            │
            ▼ Step 9

Metadata_Merge/{output_dir}/
    └── preprocessed_data.rds  ← downstream analysis input (cross-sectional default)
```

---

## Exclusion Logic

### Auto-Excluded (Steps 3–4)
- **Same-sample reruns** (same `SAMPLE`, multiple `RUN`s, identical metadata): older runs excluded,
  most recent run kept. These are assay re-runs, **not** the study's designed technical replicates.
- **Implausible values**: Age, BMI, weight, height outside plausible ranges

### Subject-Level Longitudinal Repeats (Step 8 — not excluded, split into two datasets)
- The same subject (`Record_ID`) sampled at different visits (different `SAMPLE`) is **kept**. The
  cross-sectional default keeps the most recent visit per subject; all visits go to the separate
  longitudinal dataset. See **Step 8** above.

### Flagged for Review — Not Auto-Excluded (Step 3)
- **Illogical Race/Ethnicity combinations**: Flagged with note; manually set `exclude = TRUE` in the exclusion report if needed

### Covariate-Driven Exclusions (Step 6, when `filter_na_covars = true`)
- CDX values not mapped by `dx_config` (e.g., `"0"`, `"Other"`, `"Insufficient Data"`)
- Group/Race/Ethnicity combinations not matched by any `ancestry_config` rule

---

## Manual Exclusion Editing

1. Run the full pipeline first
2. Open `output_files/sample_exclusion_report.csv`
3. Change `exclude` column: `TRUE` to keep a flagged sample, `FALSE` to exclude one
4. Re-render in filter-only mode:
   ```r
   quarto::quarto_render("Metadata_Merge_Pipeline.qmd",
                         execute_params = list(filter_only = TRUE))
   ```

---

## Adapting for a New Dataset

**Two ways to set the Step 6b groupings:**

1. **Data-driven (recommended, 2026-06).** Provide `<review_dir>/CDX_review.csv` and `Ancestry_review.csv`, each with a `FINAL_*` decision column (`FINAL_CDX_collapsed` / `FINAL_Ancestry`; blank = exclude). Step 6b builds the configs from them — edit the CSVs and re-render to regroup, no code change. The Step 6a HTML lists every CDX value and Group×Race×Ethnicity combo so you can see what needs a decision; Step 6c previews the result. (This is how the 8-plate June-2026 run was done — see `QC_Runs/QC_8plate_10_06_26/Metadata_Merge/review_2026Jun/`.)
2. **Built-in config (fallback).** If those CSVs are absent, edit the `dx_config` / `ancestry_config` block in Step 6b directly.

### Report / output provenance

Rendered reports and `output_files_*` directories are **never deleted** across runs — each records the version it was produced on:

| Report / output | Produced by | Metadata |
|---|---|---|
| `Metadata_Merge_Pipeline.html` (2026-06-08) | pre-2026-06-12 qmd | `..._n2871_plateswapcorrected_2026Jun` (now archived) |
| `output_files_n2851_2026Jun_WGScorrected/` | 2026-06-08 qmd | plate-swap-corrected n2871 (now reverted) |
| `output_files_n2851_2026May/`, `output_files_n2850_*/`, `output_files_og/` | earlier qmd versions | earlier metadata |
| `QC_Runs/QC_8plate_10_06_26/.../output_files_8plate_2026Jun/` | 2026-06-12 qmd (current) | `n4328_2026Jun12` (N = 512) |

When re-rendering, **always pass a fresh `output_dir`** so a prior run's outputs are not overwritten. The CHANGELOG comment at the top of the `.qmd` tracks code versions.

---

## Dependencies

- tidyverse
- knitr
- kableExtra
- writexl (writes the QC reports as Excel-safe `.xlsx` so genotypes aren't date-converted; pipeline falls back to `.csv` if not installed)
