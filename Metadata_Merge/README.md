# Metadata Merge

Merges post-QC biomarker data with sample metadata, applies exclusion filters, derives covariate groupings (diagnosis, ancestry, APOE status with WGS correction), and produces analysis-ready datasets including a `preprocessed_data.rds` for downstream pipelines.

## Usage

### Recommended: Run the Pipeline in RStudio

Open `Metadata_Merge_Pipeline.Rmd` in RStudio and click **Knit** (or press Ctrl/Cmd+Shift+K).

The pipeline runs in 9 steps and generates an HTML report with all decisions documented. The YAML
`knit:` function **auto-stamps the report filename** with the run name from `params$output_dir`
(e.g. `Metadata_Merge_Pipeline_n2851_2026May.html`), so successive runs don't overwrite each other.
Set `output_dir` (and other params) in the YAML header for each run.

### Render Parameters

Control pipeline behavior via the YAML `params:` block or at render time:

| Parameter | Default | Description |
|---|---|---|
| `filter_only` | `false` | Skip merge/QC; re-apply filters to existing merged files |
| `no_copy` | `false` | Skip auto-copy of NPQ files from Primary_QC |
| `include_derived_covariates` | `true` | Run Step 6 (CDX_collapsed, Ancestry, APOE4 variables) |
| `include_derived_biomarkers` | `false` | Run Step 7 (APOE4_minus_APOE, ABeta42_minus_ABeta40) |
| `filter_na_covars` | `true` | Exclude rows where CDX_collapsed or Ancestry is NA after grouping |
| `input_dir` | `"input_files"` | Where NPQ + metadata are read from. Override to run from an alternate location (e.g. the synthetic `example/` run) without touching the real `input_files/` |
| `output_dir` | `"output_files"` | Output directory for all pipeline files (created if absent) |
| `review_dir` | `"review"` | Dir holding `CDX_review.csv` / `Ancestry_review.csv` (Step 6b data-driven config) **and** `Duplicate_review.csv` (Step 3 lab duplicate decisions, `FINAL_keep`) |
| `known_repeats_file` | `null` | Path to the lab's `Dup_sample_or_ind_*.xlsx` (duplicate + longitudinal cross-check). `null` = auto-detect in `input_files/` then `../../../Datasets/` |

To override at render time (or render headlessly):
```r
rmarkdown::render("Metadata_Merge_Pipeline.Rmd",
                  params = list(filter_only = TRUE, include_derived_biomarkers = TRUE))
```
The `knit:` function in the YAML names the HTML from `params$output_dir`. When calling
`rmarkdown::render()` directly with a one-off `output_dir`, pass `output_file = "..."` to control the
name (the `knit:` hook is only honored by the RStudio **Knit** button, which uses the YAML defaults).

### Filter-Only Mode

To re-run filtering after manually editing the exclusion report (set `filter_only: true` in the YAML
and Knit, or):
```r
rmarkdown::render("Metadata_Merge_Pipeline.Rmd", params = list(filter_only = TRUE))
```

---

## Scripts

### `Metadata_Merge_Pipeline.Rmd`
Main pipeline document (R Markdown; converted from Quarto `.qmd` 2026-06-15 for repo consistency). See **Pipeline Steps** below.

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

### `Cohort_Accounting/`
Per-well exit reasons and the standing cohort-loss accounting — implements
`QC_PIPELINE_RECOMMENDATIONS_2026-07-13.md` §1 and §4. Run it **after** the
pipeline; it reads the pipeline's outputs, removes nothing and changes nothing.

Every attempted patient well gets exactly one exit reason, so **"missing from the
analysis set" stops reading as "failed QC"** — on the 50-plate run those are
11.9% and 88.1% of the loss respectively. `Cohort_Accounting/README.md` has the
full taxonomy and the verified chain; the short version is that its
`per_plate_QC_and_cohort_summary.csv` carries the QC and cohort columns side by
side and never sums them, and supersedes reading Primary QC's
`per_plate_QC_summary.csv` on its own.

```bash
cd Cohort_Accounting && Rscript run_cohort_accounting.R && Rscript test_cohort_accounting.R
```

---

## Pipeline Steps

### Step 1 — Input Files
Verifies NPQ files exist in `input_files/`; auto-copies from `Primary_QC/output_files/` if absent.

### Step 2 — Merge Data
Merges biomarker data with metadata. `APOE.geno` is normalized at this step from raw format (e.g., `"E3    E4"`) to standard `"a/b"` digit format (e.g., `"3/4"`).

### Step 2b — Known Duplicate / Longitudinal Cross-Check *(optional)*
Loads the lab's independent list (`Dup_sample_or_ind_*.xlsx`; auto-detected in `input_files/` then `../../../Datasets/`, or set via `known_repeats_file`) **and** computes metadata-level detection. The current file (`Dup_sample_or_ind_Jun2026.xlsx`) holds **62** duplicate samples and **21** longitudinal subjects — and the n4328 metadata's own structure reproduces exactly those counts, so the list serves as confirmation rather than the sole source.

**Detection is metadata-driven and cross-batch.** Duplicates (same `SAMPLE` on >1 `RUN`) and longitudinal subjects (same `Record_ID` on >1 `SAMPLE`) are detected from the **full metadata file**, which spans all batches. So a duplicate is flagged **even when its two runs' biomarkers were assayed in different batches** — this is a standard check, independent of which NPQ batch is being merged. It works the same whether the pipeline is run per-batch or on a combined dataset.

### Step 3 — Quality Control
Identifies: duplicate samples, implausible values (age, BMI, weight, height), illogical Race/Ethnicity combinations (flagged for review, not auto-excluded).

**Duplicate handling (updated 2026-06).** Detection is metadata-level (above); the **keep/exclude** decision applies only to runs that have biomarkers in the *current* merge:
- **Resolvable only when ≥2 runs are present here:** in a per-batch merge where just one run of a duplicate is present, nothing is excluded (the pair is resolved later in a combined run). In a combined merge, both runs are present and one is kept.
- **Provisional default:** the most-recent `RUN` (among runs with biomarkers here) is kept; the other(s) excluded with reason `"Duplicate SAMPLE … - provisional (most-recent RUN; pending lab review)"`.
- **Lab override:** if `<review_dir>/Duplicate_review.csv` exists, the run marked `FINAL_keep = "keep"` wins (others in that group excluded); groups with no `"keep"` fall back to the provisional rule. Mirrors the `{CDX,Ancestry}_review.csv` data-driven pattern — non-blocking, fully auditable.
- **Cross-check:** metadata-detected duplicates are reconciled against the lab list — *concordant*, *on-list-but-not-in-this-metadata*, and *in-metadata-but-not-on-list* (printed for verification).
- **Outputs:** `reports/Duplicate_review.csv` (lab-editable decision template, one row per run **across all batches**, with `has_biomarkers_here` / `n_runs_with_biomarkers_here`) and `filtered/filtered_combined_duplicates.csv` (every duplicate run that has biomarkers here — the non-longitudinal replicate feed for `Secondary_QC/Replicate_Analysis`; usable pairs appear only when ≥2 runs are present, i.e. a combined merge).

> **Deciding which run to keep:** run `Replicate_Analysis` on `filtered_combined_duplicates.csv` (from a combined merge) to compare biomarker values across the two plates, review with the lab, then record the choice in `Duplicate_review.csv` (`FINAL_keep = "keep"` on the winning run), copy it to the review dir, and re-run. Until then the most-recent run is used provisionally and downstream analyses are unaffected.

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

#### WGS plate-swap — correction REVERTED; mask-to-Sanger (updated 2026-06-15)

> **2026-06-23 update (50-plate combined / `n4328_2026Jun22` metadata).** The mask-to-Sanger idmap step
> below is **superseded** for runs on the Jun22 metadata: the swap is now corrected **upstream**, with the
> MIDI re-map baked into the metadata as the `gp_APOE` column (`APOE_WGS = coalesce(gp_APOE, APOE_WGS)`).
> Step 6d then applies a plain `APOE.geno_final = coalesce(APOE_WGS_norm, APOE.geno)` — **no idmap read, no
> masking** (verified: `APOE.geno_final == coalesce(APOE_WGS_norm, APOE.geno)`, 0 mismatches across 3552).
> Affected samples now use the **corrected WGS** call, not a Sanger fallback. The narrative below is history.

The May-2026 WGS delivery had a sample mix-up on the AD-Hispanic set that attached the wrong WGS APOE genotype to ~282 Hispanic samples. A **2026-06-08 geometric correction** (chrX fingerprint; the old `apply_wgs_plate_swap_correction.R` → `..._plateswapcorrected_2026Jun.csv`) was later found **incorrect**, reverted, and the script + its output **deleted (2026-06-15)**. The more complete chromosome-1 analysis (`../../../May_2026_WGS_QC/WGS_Plate_Swap_Summary_2026Jun11.md`) identifies the **WGS itself** as the mislabeled platform (MIDI re-array geometry). **Current pipeline behavior:** Step 6d **masks the ~282 affected samples to Sanger** via `WGS_chr1_corrected_idmap_2026Jun11.csv` — no relabel applied. The latest candidate re-map is the standalone `QC_Runs/QC_ALZ123_repro_2026Jun/WGS_Midi_Remap/remap_wgs_midi_plateswap.R` (WGS-vs-Sanger agreement 47%→99.6%) — a **PROPOSAL pending USUHS confirmation**, kept with its run (not in the repo) and not wired into the merge. For the affected samples, prefer Sanger (`APOE.geno`).

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

If the lab's known-repeats list is loaded (Step 2b), Step 8 also reconciles its Record_ID-based
detection against the list and writes `reports/Longitudinal_review.csv` (concordant / known-not-detected
— e.g. a `Pending` second draw / detected-not-listed).

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
| `*metadata*.csv` | Sample metadata (auto-detected: most recently modified file with `metadata` in the name); it is also the WGS-join source in Step 6d. Supply per run (e.g. the clean `U19_otherprojects_metadata_combined_APOE_WGS_n2871_May2026.csv`, or `n4328_*` for the 8-plate batch). The repo no longer ships a metadata file here — runs pull it from `Datasets/`. | External / `Datasets/` |

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
| `filtered/filtered_combined_duplicates.csv` | **Duplicate replicate feed** — **both** runs of every duplicate `SAMPLE` (biomarkers retained, `dup_group` + `provisional_keep`), for cross-plate concordance in `Replicate_Analysis` |
| `preprocessed_data.rds` | R list: `$data` (cross-sectional default) + `$biomarker_cols` — primary downstream input |

### QC / Error reports & mapping keys (`reports/`)

A dedicated QC-reporting section writes per-issue reports (with metadata context, for sending to the contributing groups) plus shareable mapping keys to `<output_dir>/reports/`. Checks run on the full batch present in metadata (so issues show even for samples later excluded). **Every report — including `Duplicate_review` and `Longitudinal_review` — is written as `.xlsx`** so APOE genotypes (`3/4`, `2/3`, …) stay **text** and aren't auto-converted to dates in Excel (uses `writexl`; only if it isn't installed do a couple fall back to `.csv`). Per-sample reports carry `metadata_source` (provenance) and `ID2` alongside the other identifiers.

**By-sample + summary views:** `QC_error_report_by_sample.xlsx` has **one row per flagged sample** — a logical flag per check (`implausible_value`, `unexpected_race_ethnicity`, `missing_critical_field`, `APOE_integrity`, `clinical_inconsistency`, `categorical_bad_level`, `duplicate_SAMPLE_ALIQUOT`) plus `n_issues` and a concatenated `issues` string. `QC_error_summary.xlsx` gives per-check flagged counts + batch totals (samples in batch / with ≥1 issue). Both are also sheets in the combined workbook.

**One-file option:** `QC_error_report_combined.xlsx` bundles every error/issue report into a single workbook — sheet 1 `SUMMARY` (row count + `review`/`clean` status per check), then one sheet per check — so the whole QC picture can be sent as one attachment (the per-issue files below are still written too). Mapping keys and the crosstab are excluded (they're reference keys, not error reports). Without `writexl`, falls back to `QC_error_report_combined.csv` (+ `_summary.csv`).

**Folder layout** (`<output_dir>/reports/`): overviews + mapping keys live at the **top level**; the detailed per-issue reports are grouped into **`errors/`**, **`missing_data/`**, and **`duplicates/`** sub-folders.

| File | Contents |
|---|---|
| `QC_error_report_combined.xlsx` | **all error/issue reports in one workbook** (SUMMARY + one sheet per check, incl. by-sample + summary) — sendable as a single file |
| `QC_error_report_by_sample.xlsx` | **QC errors by sample** — one row per flagged sample (per-check flags + `n_issues` + `issues`) |
| `QC_error_summary.xlsx` | per-check flagged counts + batch totals |
| `mapping_key_CDX.xlsx` | CDX raw → `CDX_collapsed` (+ n, kept/EXCLUDED) |
| `mapping_key_Ancestry.xlsx` | Group × Race × Ethnicity → Ancestry (+ n, `unexpected_race_eth`) |
| `check_CDX_vs_CaseControl_crosstab.xlsx` | CDX × Case_Control cross-tab |
| **`errors/`** | **— data-error checks —** |
| `errors/error_implausible_values.xlsx` | out-of-range age/BMI/weight/height + metadata |
| `errors/error_unexpected_race_ethnicity.xlsx` | flagged combos + metadata + mapped Ancestry |
| `errors/error_APOE_integrity.xlsx` | Sanger-vs-WGS discordance + invalid genotype format |
| `errors/error_clinical_consistency.xlsx` | AOO > age, Years_Onset, CDX vs Case_Control conflicts |
| `errors/error_categorical_unexpected_levels.xlsx` | unexpected sex / Ethnicity / Race values |
| `errors/error_duplicate_SAMPLE_ALIQUOT.xlsx` | non-unique aliquots (should be empty) |
| **`missing_data/`** | **— missingness —** |
| `missing_data/missing_data_by_sample.xlsx` | samples missing **critical** fields (CDX, Race, Ethnicity, sex, age, APOE-any) |
| `missing_data/missing_data_summary.xlsx` | n + % missing per field, incl. APOE three ways (Sanger / WGS / any) |
| **`duplicates/`** | **— repeat-sample reviews —** |
| `duplicates/Duplicate_review.xlsx` | **lab-editable** — one row per run of each duplicate; mark the run to keep with `FINAL_keep = "keep"`, copy to `<review_dir>/`, re-run (the override read also accepts `.csv`) |
| `duplicates/Longitudinal_review.xlsx` | longitudinal reconciliation vs the lab's known list (concordant / known-not-detected / detected-not-listed) |

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
            ▼ Metadata_Merge_Pipeline.Rmd (Steps 1–5)

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
- **Same-sample duplicates** (same `SAMPLE`, multiple `RUN`s, identical metadata — cross-plate
  replicates): one run kept per sample, the rest excluded. The kept run is the most recent
  **provisionally** (overridable by the lab via `<review_dir>/Duplicate_review.csv`). **Both** runs
  are preserved in `filtered/filtered_combined_duplicates.csv` for the replicate sub-analysis — the
  exclusion only affects the cross-sectional analysis set. See **Step 3** above.
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

> **Note — a CDX value absent from `CDX_review.csv` entirely is dropped too,** by
> the same `is.na()` test as a deliberately excluded one. `apply_category_config()`
> warns about these "unaccounted" levels, but the Step 6d chunk sets
> `warning: false`, so in the rendered report a new CDX spelling drops its samples
> with no visible signal. `Cohort_Accounting/` names any such value in its console
> summary.

### Counting what these exclusions cost
This section says *which* rules drop rows; it does not say **how many, whose, or
compared with what**. That is `Cohort_Accounting/` — attempted → in-cohort with
one exit reason per well and cuts by year, site and ancestry. On the 50-plate run
the covariate-driven exclusions above are the single largest category of loss in
the whole pipeline (361 wells on diagnosis, against 77 for all of Primary QC),
which is not visible from the exclusion logic alone.

---

## Manual Exclusion Editing

1. Run the full pipeline first
2. Open `output_files/sample_exclusion_report.csv`
3. Change the `exclude` column: **`TRUE` drops the sample, `FALSE` keeps it.** Rows merely flagged
   for review are written with `exclude = FALSE` (kept) — set one to `TRUE` to act on the flag, and
   leave it `FALSE` to keep the sample. Step 5 collects the `exclude == TRUE` rows into
   `samples_to_exclude` and filters them out of the merged datasets.
4. Re-render in filter-only mode:
   ```r
   rmarkdown::render("Metadata_Merge_Pipeline.Rmd",
                     params = list(filter_only = TRUE))
   ```

---

## Adapting for a New Dataset

**Two ways to set the Step 6b groupings:**

1. **Data-driven (recommended, 2026-06).** Provide `<review_dir>/CDX_review.csv` and `Ancestry_review.csv`, each with a `FINAL_*` decision column (`FINAL_CDX_collapsed` / `FINAL_Ancestry`; blank = exclude). Step 6b builds the configs from them — edit the CSVs and re-render to regroup, no code change. The Step 6a HTML lists every CDX value and Group×Race×Ethnicity combo so you can see what needs a decision; Step 6c previews the result. (This is how the 8-plate June-2026 run was done — see `QC_Runs/QC_8plate_10_06_26/Metadata_Merge/review_2026Jun/`.)
2. **Built-in config (fallback).** If those CSVs are absent, edit the `dx_config` / `ancestry_config` block in Step 6b directly.

### Report / output provenance

Current runs live in `QC_Runs/`; the repo keeps code only. Old repo outputs were cleaned up **2026-06-15**:

| Report / output | Status (2026-06-15) |
|---|---|
| `output_files_n2851_2026Jun_WGScorrected/` + the `apply_wgs_plate_swap_correction.R` script + `..._n2871_plateswapcorrected_2026Jun.csv` | **DELETED** — the refuted 2026-06-08 plate-swap correction (superseded by the mask-to-Sanger Step 6d; the latest candidate re-map is the standalone `QC_Runs/QC_ALZ123_repro_2026Jun/WGS_Midi_Remap/remap_wgs_midi_plateswap.R`) |
| `output_files_n2851_2026May/`, `output_files_n2850_*/`, `..._n2871_2026May_STALE_457wgs/`, `output_files_og/` | **ARCHIVED** in `_archive_stale_outputs_2026-06-15/` (gitignored; see `MANIFEST_pre_move.txt`) — superseded; the n2851 merge is reproduced by `QC_Runs/QC_ALZ123_repro_2026Jun/` |
| `QC_Runs/QC_ALZ123_repro_2026Jun/.../output_files_ALZ123_repro_2026Jun/` | **CURRENT** ALZ123 run — `n2871_May2026` metadata, mask-to-Sanger (merged 2825 → N = 2615) |
| `QC_Runs/QC_8plate_10_06_26/.../output_files_8plate_2026Jun/` | **CURRENT** 8-plate run — `n4328_2026Jun12` metadata (merged 632 → N = 512) |

When re-rendering, **always pass a fresh `output_dir`** so a prior run's outputs are not overwritten (the `knit:` function names the HTML from `output_dir`, so a fresh dir also gives a fresh report name). The CHANGELOG comment at the top of the `.Rmd` tracks code versions.

---

## Dependencies

- tidyverse
- knitr
- kableExtra
- writexl (writes the QC reports as Excel-safe `.xlsx` so genotypes aren't date-converted; pipeline falls back to `.csv` if not installed)
