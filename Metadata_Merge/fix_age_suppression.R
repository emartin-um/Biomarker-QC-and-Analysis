# fix_age_suppression.R
# In the Apr27 U19 metadata update, 60 samples with age >= 90 had age_at_subject
# set to NA (HIPAA suppression). This script sets those values to 90 so they
# can still be used in age-adjusted analyses.
#
# Run from the Metadata_Merge/ directory after combine_metadata.R has produced
# the new combined CSV.
#
# Input/output: the combined metadata CSV (reads and overwrites in place).

library(readr)

# ── Config ────────────────────────────────────────────────────────────────────
DATASETS_DIR   <- "../../../Datasets"
combined_file  <- file.path(DATASETS_DIR, "U19_Alamar_metadata_combined_n2850_2026Apr27.csv")

# SAMPLE_ALIQUOT IDs whose age_at_subject was suppressed (set to NA) in the
# Apr27 U19 metadata update. All had age >= 90 in the Jan16 version.
AGE_SUPPRESSED_IDS <- c(
  "202304238-22",
  "202304269-06",
  "202314866-06",
  "202314889-06",
  "202316227-06",
  "202403987-06",
  "202408029-30",
  "202408031-37",
  "202408041-38",
  "202408115-37",
  "202409026-30",
  "202412465-38",
  "202412537-33",
  "202415247-40",
  "202415250-40",
  "202415539-39",
  "202415593-40",
  "202415794-39",
  "202419918-32",
  "202421805-29",
  "202421852-39",
  "202421915-20",
  "202422027-41",
  "202423752-43",
  "202423753-30",
  "202423975-22",
  "202424023-35",
  "202427484-23",
  "202427595-50",
  "202427623-48",
  "202428044-32",
  "202428091-25",
  "202428976-44",
  "202434371-37",
  "202434566-35",
  "202434735-30",
  "202436508-44",
  "202436954-39",
  "202504558-05",
  "202514219-11",
  "202514246-11",
  "202514319-11",
  "202514507-11",
  "202514518-11",
  "202514564-11",
  "202514580-11",
  "202514581-11",
  "202514588-11",
  "202514589-11",
  "202514592-11",
  "202514647-11",
  "202514648-11",
  "202514655-11",
  "202514657-11",
  "202514716-11",
  "202514744-11",
  "202514754-11",
  "202514756-17",
  "202514974-11",
  "202515017-11"
)

# ── Load ──────────────────────────────────────────────────────────────────────
meta <- read_csv(combined_file, show_col_types = FALSE)
cat("Loaded:", nrow(meta), "rows from", combined_file, "\n")

# ── Patch ─────────────────────────────────────────────────────────────────────
target_rows <- meta$SAMPLE_ALIQUOT %in% AGE_SUPPRESSED_IDS

# Confirm these rows currently have NA age (sanity check)
na_before <- sum(is.na(meta$age_at_subject[target_rows]))
cat("Samples matched:", sum(target_rows), " | Currently NA age:", na_before, "\n")

meta$age_at_subject[target_rows] <- 90

na_after <- sum(is.na(meta$age_at_subject[target_rows]))
cat("After fix — NA age in target rows:", na_after, "\n")
cat("Overall NA age remaining:", sum(is.na(meta$age_at_subject)), "\n")

# ── Save ──────────────────────────────────────────────────────────────────────
write_csv(meta, combined_file)
cat("Saved:", combined_file, "\n")
