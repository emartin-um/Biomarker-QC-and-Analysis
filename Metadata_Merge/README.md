# Metadata Merge

Merges post-QC biomarker data with sample metadata, applies exclusion filters, derives covariate groupings (diagnosis, ancestry, APOE status), and produces analysis-ready datasets including a `preprocessed_data.rds` for downstream pipelines.

## Usage

### Recommended: Run the Pipeline in RStudio

Open `Metadata_Merge_Pipeline.qmd` in RStudio and click **Render** (or press Ctrl/Cmd+Shift+K).

The pipeline runs in 8 steps and generates an HTML report with all decisions documented.

### Render Parameters

Control pipeline behavior via the YAML `params:` block or at render time:

| Parameter | Default | Description |
|---|---|---|
| `filter_only` | `false` | Skip merge/QC; re-apply filters to existing merged files |
| `no_copy` | `false` | Skip auto-copy of NPQ files from Primary_QC |
| `include_derived_covariates` | `true` | Run Step 6 (CDX_collapsed, Ancestry, APOE4 variables) |
| `include_derived_biomarkers` | `false` | Run Step 7 (APOE4_minus_APOE, ABeta42_minus_ABeta40) |
| `filter_na_covars` | `true` | Exclude rows where CDX_collapsed or Ancestry is NA after grouping |

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
Removes auto-excluded samples; writes filtered CSVs. Exclusion report can be manually edited before re-running in filter-only mode.

### Step 5 — Exclusion Summary
Reports sample counts before/after filtering.

### Step 6 — Covariate Groupings *(optional: `include_derived_covariates`)*

Config-driven, designed to be re-configured for each new dataset. Four sub-steps:

- **6a Explore**: Prints CDX distribution and Group × Race × Ethnicity combinations to inform grouping decisions.
- **6b Config**: User-editable block defining `dx_config` (diagnosis grouping) and `ancestry_config` (ancestry grouping) using `make_category_config()` / `make_grouping_config()`. **Edit this block when processing a new dataset.**
- **6c Preview**: Shows decision tables and any unmatched rows *before* applying anything — no data is modified.
- **6d Apply**: Applies configs, computes APOE4 carrier variables, validates group sizes, saves updated CSV.

Derived columns produced:

| Column | Type | Description |
|---|---|---|
| `CDX_collapsed` | character | Collapsed diagnosis: NCI / MCI / AD / Dementia_Other |
| `Ancestry` | character | Grouped ancestry: AA / AFDC / HI_WH / HI_MU / HI_BL |
| `APOE4_carrier` | logical | TRUE if any E4 allele present |
| `APOE4_count` | integer | Number of E4 alleles (0 / 1 / 2) |

When `filter_na_covars = true` (default), rows where CDX_collapsed or Ancestry is NA are excluded.

### Step 7 — Derived Biomarkers *(optional: `include_derived_biomarkers`, default FALSE)*

Adds log2 difference columns: `APOE4_minus_APOE` and `ABeta42_minus_ABeta40`.

### Step 8 — Save RDS

Writes `output_files/preprocessed_data.rds`:
```r
list(
  data           = filtered_combined,   # all rows/columns after steps 1–7
  biomarker_cols = biomarker_cols       # character vector of NPQ biomarker column names
)
```
This RDS is the primary input for downstream analysis pipelines (e.g., APOE_Amyloid_Tau).

---

## Input Files

Place in `input_files/` (not tracked in git):

| File | Description | Source |
|---|---|---|
| `NPQ_post_QC.csv` | Standard post-QC biomarker data | Primary_QC output |
| `NPQ_Low_post_QC.csv` | Low detectability biomarker data | Primary_QC output |
| `*metadata*.csv` | Sample metadata (auto-detected by filename) | External |

---

## Output Files

Generated in `output_files/` (not tracked in git):

| File | Description |
|---|---|
| `merged_standard_post_QC.csv` | Merged standard biomarkers + metadata |
| `merged_low_post_QC.csv` | Merged low biomarkers + metadata |
| `merged_combined_post_QC.csv` | All biomarkers + metadata |
| `sample_exclusion_report.csv` | Exclusion flags and reasons per sample |
| `samples_not_in_metadata.csv` | NPQ samples missing from metadata (if any) |
| `filtered/filtered_standard_post_QC.csv` | Filtered standard biomarkers |
| `filtered/filtered_low_post_QC.csv` | Filtered low biomarkers |
| `filtered/filtered_combined_post_QC.csv` | Filtered combined (+ derived columns if enabled) |
| `preprocessed_data.rds` | R list: `$data` + `$biomarker_cols` — primary downstream input |

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

Metadata_Merge/output_files/
    ├── merged_*.csv
    └── sample_exclusion_report.csv  ← Review/edit if needed
            │
            ▼ Step 6 (covariate groupings) + Step 7 (derived biomarkers)

Metadata_Merge/output_files/filtered/
    └── filtered_combined_post_QC.csv  ← with CDX_collapsed, Ancestry, APOE4_carrier, APOE4_count
            │
            ▼ Step 8

Metadata_Merge/output_files/
    └── preprocessed_data.rds  ← downstream analysis input
```

---

## Exclusion Logic

### Auto-Excluded (Steps 3–4)
- **Duplicate samples**: Older runs excluded, most recent kept
- **Implausible values**: Age, BMI, weight, height outside plausible ranges

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

The only block that requires editing for a new dataset is **Step 6b** in `Metadata_Merge_Pipeline.qmd`. Update `dx_config` and `ancestry_config` to reflect the new dataset's CDX levels and Group/Race/Ethnicity structure. Use the Step 6a exploration output and Step 6c decision previews to verify before committing.

---

## Dependencies

- tidyverse
- knitr
- kableExtra
