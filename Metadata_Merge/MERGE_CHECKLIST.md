# Pre-Merge Checklist: `dev/metadata-merge` → `main`

---

## 1. Render Tests

- [ ] **Full default render** — open `Metadata_Merge_Pipeline.qmd` with all default params and render to HTML. Confirm no errors, report generates cleanly through all 8 steps.
- [ ] **filter_only mode** — re-render with `filter_only: true` after editing `sample_exclusion_report.csv`. Confirm filters apply and the re-generated CSVs and RDS reflect the changes.
- [ ] **`include_derived_covariates: false`** — render with this set to false. Confirm Step 6 is skipped and `filtered_combined_post_QC.csv` / RDS contain **no** `CDX_collapsed`, `Ancestry`, `APOE4_carrier`, or `APOE4_count` columns.
- [ ] **`filter_na_covars: false`** — render with this set to false. Confirm rows with unmatched ancestry are retained in the output with `Ancestry = NA` rather than being dropped.

---

## 2. APOE.geno Normalization

- [ ] After a full run, inspect `filtered_combined$APOE.geno` (or spot-check the output CSV). Values should be in `"3/4"` format — no `"E3    E4"` style strings and no NAs introduced by the regex that weren't already NA in the raw metadata.
- [ ] Verify the two `str_extract` patterns used for normalization (`"\\d"` for first allele, `"\\d$"` for second) handle all real genotype strings in the metadata — including edge cases like `"E2    E3"` or `"E4    E4"`.

---

## 3. Output File Correctness

- [ ] **`filtered_combined_post_QC.csv`** — confirm it has `CDX_collapsed`, `Ancestry`, `APOE4_carrier` (logical TRUE/FALSE), `APOE4_count` (integer 0/1/2) columns.
- [ ] **`preprocessed_data.rds`** — load it and confirm structure: `list(data = ..., biomarker_cols = ...)`.
- [ ] **`biomarker_cols`** — inspect `dat_list$biomarker_cols`. Verify it contains only NPQ biomarker columns — not metadata columns (age, sex, Site, etc.), not pipeline-derived columns (`CDX_collapsed`, `Ancestry`, `APOE4_carrier`, `APOE4_count`), and not the calculated biomarkers from Step 7 unless `include_derived_biomarkers` is true.
- [ ] **`filtered_standard_post_QC.csv` and `filtered_low_post_QC.csv`** — still present and non-empty.

---

## 4. MCI Progression Analysis Compatibility

This is the biggest potential conflict. `MCI_Data_Preprocessing.Rmd` derives `CDX_collapsed` and `Ancestry` itself. Now that Metadata_Merge writes them into the RDS, there are two possibilities when MCI preprocessing runs:

- [ ] **Open `MCI_Data_Preprocessing.Rmd`** and check whether its CDX/Ancestry derivation code will **overwrite** the columns already in the RDS, or whether it checks for column existence first.
  - If it overwrites: verify the grouping logic is identical to what's in Metadata_Merge Step 6b — they must produce the same category labels.
  - If it's redundant: consider removing the derivation from MCI preprocessing so the RDS is the single source of truth (but that's a future task, not a blocker for this merge).
- [ ] Confirm `APOE4_carrier` arrives as **logical** in the RDS. `MCI_Data_Preprocessing.Rmd` includes it in `export_covariate_cols` — verify `sum(APOE4_carrier)` and logical operations work after loading.
- [ ] Confirm `dat_list$biomarker_cols` is compatible with how `MCI_Progression_Analysis.Rmd` uses it to select biomarker columns from `dat_list$data`.

---

## 5. APOE Amyloid/Tau Pipeline Compatibility

- [ ] With the normalized `"3/4"` APOE.geno format in the RDS, confirm `str_extract(APOE.geno, "^\\d")` in `01_APOE_data_prep.Rmd` correctly extracts the first allele digit (it should, since `"3/4"` starts with a digit).
- [ ] Run through the first chunk of `01_APOE_data_prep.Rmd` and check `table(dat$geno_cat, useNA = "ifany")` — should have non-zero counts in E2E3, E3E3, E3E4, etc., and zero NAs (or only NAs for samples that had NA genotype in the original metadata).

---

## 6. Primary_QC Handoff (Upstream Compatibility)

- [ ] Confirm the auto-copy logic still works: when `input_files/` is empty, the pipeline correctly finds and copies `NPQ_*_post_QC.csv` files from `Primary_QC/output_files/`.
- [ ] Confirm metadata auto-detection (`*metadata*.csv`) still picks up the metadata file correctly.

---

## 7. covariate_explorer.R Duplication

`Metadata_Merge/covariate_explorer.R` and `MCI_Progression_Analysis/covariate_explorer.R` are now independent copies.

- [ ] Confirm both files are currently identical:
  ```bash
  diff Metadata_Merge/covariate_explorer.R \
       In_Development/MCI_Progression_Analysis/covariate_explorer.R
  ```
- [ ] Decide on a drift strategy: either accept that they can evolve independently (and document it), or add a note in both READMEs that they must be kept in sync manually when one is updated.

---

## 8. Legacy File Cleanup

- [ ] Confirm `Metadata_Merge/old/` is excluded by `.gitignore` — the repo-level `.gitignore` has `/old/`, so this should already be handled. Verify with `git status`.
- [ ] Confirm the three deleted legacy scripts (`File_Merge.qmd`, `Apply_Filters.qmd`, `Add_Derived_Biomarkers.qmd`) are not referenced anywhere else in the repo.

---

## 9. README Accuracy

- [ ] Read through the new `README.md` against the actual pipeline. Specifically verify the Step 6b "adapting for a new dataset" instructions are accurate.
- [ ] Confirm the output files table includes `preprocessed_data.rds`.

---

## Summary: Highest-Risk Items

| Risk | Severity | What to check |
|---|---|---|
| MCI preprocessing re-derives CDX_collapsed/Ancestry, possibly with different labels | **High** | Compare grouping logic between Metadata_Merge Step 6b config and `MCI_Data_Preprocessing.Rmd` |
| `biomarker_cols` in RDS missing or over-inclusive | **High** | Inspect `dat_list$biomarker_cols` after a real run |
| APOE.geno normalization introduces NAs or wrong format | **Medium** | Spot-check 10–20 rows of output CSV |
| covariate_explorer.R copies diverge in future | **Low** (now) | Diff immediately; document sync expectation |
