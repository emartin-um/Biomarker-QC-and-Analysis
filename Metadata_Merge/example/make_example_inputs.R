#!/usr/bin/env Rscript
# make_example_inputs.R — generate a tiny SYNTHETIC dataset so the Metadata_Merge pipeline can be
# run/verified end-to-end without any real (PHI) data. All IDs, ages, genotypes are fabricated.
#
# Writes into ./input_files (relative to this script's dir):
#   - NPQ_post_QC.csv, NPQ_Low_post_QC.csv   (SampleID + biomarker columns; samples as rows)
#   - U19_example_metadata_synthetic.csv     (the ~24-column metadata the pipeline expects)
#   - Dup_sample_or_ind_EXAMPLE.xlsx         (known duplicate/longitudinal cross-check list)
#
# The data deliberately includes a COMBINED-merge scenario so the duplicate/longitudinal feature
# fully exercises: 3 duplicate samples (same SAMPLE on two RUNs, BOTH runs present) and 2 longitudinal
# subjects (same Record_ID, two SAMPLEs). See ../example/README.md and the repo README "Running a QC batch".
#
# Usage:  Rscript make_example_inputs.R       (run from Metadata_Merge/example/)

suppressWarnings(suppressMessages({ library(readr); library(dplyr); library(tibble) }))
set.seed(42)

here <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(here) == 0 || !nzchar(here)) here <- "."
indir <- file.path(here, "input_files")
dir.create(indir, showWarnings = FALSE, recursive = TRUE)

# ── Sample design ───────────────────────────────────────────────────────────────
RUN1 <- "20260101-1000_Bay1"; RUN2 <- "20260201-1000_Bay2"; RUN3 <- "20260115-1000_Bay3"

# singletons (one row each)
singles <- tribble(
  ~SAMPLE, ~RUN,  ~Record_ID, ~Group, ~Race, ~Ethnicity, ~sex, ~age, ~CDX,  ~APOE,
  3001, RUN1, "S3001", "HI", "WH", "HI", "F", 71, "NCI", "3/3",
  3002, RUN1, "S3002", "AA", "BL", "NH", "M", 68, "MCI", "3/4",
  3003, RUN2, "S3003", "NA", "WH", "NH", "F", 80, "AD",  "4/4",
  3004, RUN2, "S3004", "HI", "BL", "HI", "M", 63, "NCI", "2/3",
  3005, RUN3, "S3005", "HI", "WH", "HI", "F", 75, "MCI", "3/4",
  3006, RUN3, "S3006", "AA", "BL", "NH", "M", 59, "NCI", "3/3",
  3007, RUN1, "S3007", "NA", "WH", "NH", "F", 84, "AD",  "3/4",
  3008, RUN2, "S3008", "HI", "WH", "HI", "M", 66, "MCI", "2/4"
) %>% mutate(aliquot = "01")

# 3 DUPLICATES: same SAMPLE on two RUNs (both present in NPQ -> resolvable, replicate pairs)
dup_ids <- tibble(SAMPLE = c(1001, 1002, 1003),
                  Record_ID = c("S1001","S1002","S1003"),
                  Group = c("HI","AA","HI"), Race = c("WH","BL","BL"),
                  Ethnicity = c("HI","NH","HI"), sex = c("F","M","F"),
                  age = c(72, 70, 65), CDX = c("MCI","NCI","AD"), APOE = c("3/4","3/3","4/4"))
dups <- bind_rows(
  dup_ids %>% mutate(RUN = RUN1, aliquot = "01"),
  dup_ids %>% mutate(RUN = RUN2, aliquot = "02")   # second run = more recent
)

# 2 LONGITUDINAL subjects: same Record_ID, two different SAMPLEs (two visits)
longi <- tribble(
  ~SAMPLE, ~RUN,  ~Record_ID, ~Group, ~Race, ~Ethnicity, ~sex, ~age, ~CDX, ~APOE,
  2001, RUN1, "SUBJ_L1", "HI", "WH", "HI", "F", 70, "NCI", "3/4",
  2002, RUN3, "SUBJ_L1", "HI", "WH", "HI", "F", 72, "MCI", "3/4",   # same subject, later visit
  2003, RUN1, "SUBJ_L2", "AA", "BL", "NH", "M", 66, "MCI", "3/3",
  2004, RUN2, "SUBJ_L2", "AA", "BL", "NH", "M", 69, "AD",  "3/3"
) %>% mutate(aliquot = "01")

meta <- bind_rows(singles, dups, longi) %>%
  mutate(
    SAMPLE_ALIQUOT = paste0(SAMPLE, "-", aliquot),
    Case_Control   = if_else(CDX == "NCI", "Control", "Case"),
    AOO            = if_else(CDX == "AD", age - 3, NA_real_),
    Years_Onset    = if_else(CDX == "AD", 3, NA_real_),
    height_inches  = round(runif(n(), 60, 72)),
    weight_lb      = round(runif(n(), 120, 220)),
    BMI            = round(703 * weight_lb / height_inches^2, 1),
    Site           = sample(c("CU","IHG","NC"), n(), replace = TRUE),
    `Country/State`= "Example State",
    country_of_birth = "EX",
    ID2            = paste0("EX-", SAMPLE_ALIQUOT),
    APOE_WGS       = NA_character_,
    metadata_source = "EXAMPLE_synthetic",
    age_at_subject = age
  ) %>%
  select(metadata_source, RUN, SAMPLE_ALIQUOT, SAMPLE, Record_ID, Group, Ethnicity, Race, sex,
         age_at_subject, CDX, Case_Control, AOO, Years_Onset, APOE, height_inches, weight_lb, BMI,
         Site, `Country/State`, country_of_birth, ID2, APOE_WGS)

write_csv(meta, file.path(indir, "U19_example_metadata_synthetic.csv"))

# ── Biomarker matrices (SampleID + biomarker columns; one row per SAMPLE_ALIQUOT) ──
ids <- meta$SAMPLE_ALIQUOT
mk_npq <- function(bm_names) {
  m <- as.data.frame(matrix(round(rnorm(length(ids) * length(bm_names), 5, 1.5), 3),
                            nrow = length(ids), dimnames = list(NULL, bm_names)))
  bind_cols(tibble(SampleID = ids), m)
}
write_csv(mk_npq(paste0("BM", 1:8)),  file.path(indir, "NPQ_post_QC.csv"))
write_csv(mk_npq(paste0("LOW", 1:3)), file.path(indir, "NPQ_Low_post_QC.csv"))

# ── Known duplicate / longitudinal list (mirrors the metadata; for the cross-check) ──
if (requireNamespace("writexl", quietly = TRUE)) {
  dup_sheet <- meta %>% group_by(SAMPLE) %>% filter(dplyr::n_distinct(RUN) > 1) %>% ungroup() %>%
    transmute(metadata_source, RUN, SAMPLE_ALIQUOT, SAMPLE, ResultQC = "PASS", Record_ID,
              Group, Ethnicity, Race, sex, age_at_subject, CDX, Case_Control, APOE, ID2)
  long_sheet <- meta %>% group_by(Record_ID) %>% filter(dplyr::n_distinct(SAMPLE) > 1) %>% ungroup() %>%
    transmute(metadata_source, RUN, SAMPLE_ALIQUOT, SAMPLE, ResultQC = "PASS", Record_ID,
              Group, Ethnicity, Race, sex, age_at_subject, CDX, Case_Control, APOE, ID2)
  writexl::write_xlsx(list("Duplicate sample" = dup_sheet, "Longitudinal sample" = long_sheet),
                      file.path(indir, "Dup_sample_or_ind_EXAMPLE.xlsx"))
}

cat("Wrote example inputs to", normalizePath(indir), ":\n")
cat(" -", nrow(meta), "metadata rows;",
    dplyr::n_distinct(meta$SAMPLE), "distinct SAMPLE;",
    "3 duplicates (2 runs each), 2 longitudinal subjects (2 visits each).\n")
