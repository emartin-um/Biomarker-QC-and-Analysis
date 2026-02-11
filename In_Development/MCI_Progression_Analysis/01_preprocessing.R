################################################################################
# 01_preprocessing.R
# Functions for sample-level preprocessing (APOE4 status, missing data, filtering)
#
# NOTE: Diagnosis and ancestry groupings are now handled by covariate_explorer.R
# via run_covariate_setup.R. This script handles the remaining preprocessing
# steps that operate on the already-grouped data.
################################################################################

library(tidyverse)

#' Filter minimum sample size per group
#'
#' @param dat_list Output from load_alamar_data
#' @param min_n Minimum samples per Diagnosis x Ancestry group
#' @param grouping_vars Variables to group by (default: c("Ancestry", "CDX_collapsed"))
#' @return Filtered dat_list
filter_min_samples <- function(dat_list,
                               min_n = 5,
                               grouping_vars = c("Ancestry", "CDX_collapsed")) {

  dat <- dat_list$data

  # Check group sizes
  group_counts <- dat %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarize(n = n(), .groups = "drop")

  message("Group sizes:")
  print(group_counts)

  # Filter
  dat <- dat %>%
    group_by(across(all_of(grouping_vars))) %>%
    filter(n() >= min_n) %>%
    ungroup()

  message(sprintf("\nRetained %d samples after filtering (min_n=%d)",
                  nrow(dat), min_n))

  dat_list$data <- dat
  return(dat_list)
}

#' Create APOE4 carrier status
#'
#' @param dat_list Output from load_alamar_data
#' @return Modified dat_list with APOE4_carrier column
create_apoe4_status <- function(dat_list) {

  dat <- dat_list$data

  # Extract APOE4 status from APOE.geno
  # Format appears to be "E2    E3" or similar
  dat <- dat %>%
    mutate(
      APOE4_carrier = case_when(
        grepl("4", APOE.geno) ~ "APOE4+",
        !is.na(APOE.geno) ~ "APOE4-",
        TRUE ~ NA_character_
      ),
      APOE4_count = case_when(
        grepl("4.*4", APOE.geno) ~ 2,
        grepl("4", APOE.geno) ~ 1,
        !is.na(APOE.geno) ~ 0,
        TRUE ~ NA_real_
      )
    )

  message("APOE4 carrier status:")
  print(table(dat$APOE4_carrier, useNA = "ifany"))

  dat_list$data <- dat
  return(dat_list)
}

#' Remove samples with excessive missing biomarkers
#'
#' @param dat_list Output from load_alamar_data
#' @param max_missing_pct Maximum percent missing biomarkers (default 20%)
#' @return Filtered dat_list
filter_missing_biomarkers <- function(dat_list, max_missing_pct = 0.2) {

  biomarker_mat <- get_biomarker_matrix(dat_list)

  # Calculate percent missing per sample
  pct_missing <- rowSums(is.na(biomarker_mat)) / ncol(biomarker_mat)

  message(sprintf("Missing biomarker distribution (min=%.1f%%, max=%.1f%%, median=%.1f%%)",
                  min(pct_missing)*100, max(pct_missing)*100, median(pct_missing)*100))

  # Filter
  keep_samples <- pct_missing <= max_missing_pct
  dat_list$data <- dat_list$data[keep_samples, ]

  message(sprintf("Retained %d/%d samples with <%.0f%% missing biomarkers",
                  sum(keep_samples), length(keep_samples), max_missing_pct*100))

  return(dat_list)
}

#' Standard preprocessing pipeline
#'
#' Loads pre-processed data (with CDX_collapsed and Ancestry already defined
#' by run_covariate_setup.R) and applies remaining preprocessing steps.
#'
#' @param processed_rds_path Path to processed_dat_list.rds from run_covariate_setup.R
#' @param max_missing_pct Maximum percent missing biomarkers (default 20%)
#' @param min_n Minimum samples per group (default 5)
#' @return Preprocessed dat_list
preprocess_standard <- function(processed_rds_path,
                                max_missing_pct = 0.2,
                                min_n = 5) {

  dat_list <- readRDS(processed_rds_path)

  # Validate that covariate setup has been run
  required_cols <- c("CDX_collapsed", "Ancestry")
  missing_cols <- setdiff(required_cols, names(dat_list$data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing columns: %s\nRun run_covariate_setup.R first to define groupings.",
      paste(missing_cols, collapse = ", ")
    ))
  }

  message("=== Starting preprocessing pipeline ===\n")
  message(sprintf("Loaded %d samples with CDX_collapsed and Ancestry defined.",
                  nrow(dat_list$data)))

  dat_list <- dat_list %>%
    create_apoe4_status() %>%
    filter_missing_biomarkers(max_missing_pct = max_missing_pct) %>%
    filter_min_samples(min_n = min_n)

  message("\n=== Preprocessing complete ===")
  message(sprintf("Final dataset: %d samples x %d biomarkers",
                  nrow(dat_list$data), length(dat_list$biomarker_cols)))

  return(dat_list)
}
