# combine_metadata.R
# Combine U19 and OtherProjects metadata for n2850 pipeline run.
#
# Run from the Metadata_Merge/ directory (e.g., set working directory in RStudio
# to the Metadata_Merge/ folder, or source this script from there).
#
# Inputs  (Alamar_QC_Protocol/Datasets/ — outside the repo, do NOT overwrite):
#   U19_Alamar_metadata_2026Jan16.csv
#   OtherProjects_metadata_03052026.xlsx
#
# Output  (also written to Datasets/, then manually copied to input_files/):
#   U19_OtherProjects_metadata_combined_n2850_2026Mar.csv

library(tidyverse)
library(readxl)

# ── Paths ────────────────────────────────────────────────────────────────────
# Adjust DATASETS_DIR if running from a different working directory.
DATASETS_DIR <- "../../../Datasets"

meta_orig_path  <- file.path(DATASETS_DIR, "U19_Alamar_metadata_2026Jan16.csv")
meta_new_path   <- file.path(DATASETS_DIR, "OtherProjects_metadata_03052026.xlsx")
combined_path   <- file.path(DATASETS_DIR, "U19_OtherProjects_metadata_combined_n2850_2026Mar.csv")

stopifnot(file.exists(meta_orig_path))
stopifnot(file.exists(meta_new_path))

# ── Load ─────────────────────────────────────────────────────────────────────
meta_orig <- read_csv(meta_orig_path, show_col_types = FALSE) |>
  mutate(metadata_source = "U19_2026Jan16")

meta_new <- read_xlsx(meta_new_path) |>
  mutate(metadata_source = "OtherProjects_03052026")

cat("U19 rows:           ", nrow(meta_orig), "\n")
cat("OtherProjects rows: ", nrow(meta_new),  "\n")

# ── Combine ───────────────────────────────────────────────────────────────────
meta_combined <- bind_rows(meta_orig, meta_new)

cat("Combined rows:      ", nrow(meta_combined), "\n")

# ── Save ─────────────────────────────────────────────────────────────────────
write_csv(meta_combined, combined_path)
cat("\nSaved:", combined_path, "\n")

# ── Next step ────────────────────────────────────────────────────────────────
# Copy the combined CSV to Metadata_Merge/input_files/ so the pipeline
# auto-detects it (find_metadata_file() picks the most recently modified
# file with 'metadata' in the name).
#
# fs::file_copy(combined_path, "input_files/U19_OtherProjects_metadata_combined_n2850_2026Mar.csv")
# or manually copy in Finder / terminal.
