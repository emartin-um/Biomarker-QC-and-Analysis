# Metadata_Merge — synthetic example run

A tiny, **fully synthetic** dataset (no real/PHI data — all IDs, ages, genotypes fabricated) so the
Metadata_Merge pipeline can be run and verified end-to-end without touching real data. Useful for
smoke-testing a code change, onboarding, or seeing what the outputs look like.

The data deliberately includes a **combined-merge scenario** that exercises the duplicate/longitudinal
features: **3 duplicate samples** (same `SAMPLE` on two `RUN`s, both present → resolvable, replicate
pairs) and **2 longitudinal subjects** (same `Record_ID`, two `SAMPLE`s → two visits).

## Run it (from the `Metadata_Merge/` directory)

```r
# 1. generate the synthetic inputs into example/input_files/
source("example/make_example_inputs.R")        # or: Rscript example/make_example_inputs.R

# 2. render the pipeline pointed at the example inputs/outputs (does NOT touch the real input_files/)
rmarkdown::render(
  "Metadata_Merge_Pipeline.Rmd",
  params = list(
    input_dir                  = "example/input_files",
    output_dir                 = "example/output_files",
    include_derived_covariates = FALSE,   # smoke test (see note below)
    no_copy                    = TRUE,
    known_repeats_file         = "example/input_files/Dup_sample_or_ind_EXAMPLE.xlsx"
  ),
  output_file = "Metadata_Merge_Pipeline_EXAMPLE.html",
  output_dir  = "example"
)
```

Outputs land only in `example/output_files/` and `example/*.html` (both gitignored). Expected:
cross-sectional **13** rows, longitudinal-only **4** rows / **2** subjects, duplicates feed **6** rows /
**3** pairs, plus `reports/Duplicate_review.csv`, `reports/Longitudinal_review.csv`, and
`preprocessed_data.rds`.

## Notes

- The committed artifact is the **generator** (`make_example_inputs.R`) + this README; the generated
  CSV/XLSX and all outputs are gitignored (regenerate with step 1). Nothing here is real data.
- **`include_derived_covariates = FALSE`** keeps the example self-contained and robust — it exercises the
  merge, the metadata-driven duplicate/longitudinal cross-check (Step 3 + Step 8), filtering, the
  longitudinal split, and the RDS save. The covariate-grouping path (Step 6) needs ancestry/CDX values
  that match the built-in config and is covered by real-data runs; the combined `QC_error_report` is
  produced only when covariates are on.
