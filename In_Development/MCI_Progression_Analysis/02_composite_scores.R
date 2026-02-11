################################################################################
# 02_composite_scores.R
# Create multi-biomarker composite scores
################################################################################

library(tidyverse)

#' Calculate ATN framework scores
#'
#' Biomarker data are NPQ (log2-transformed) values. Scores are computed as
#' z-scored averages so that markers on different scales contribute equally.
#'
#' A_score: Amyloid. Uses A_42 - A_40 (log2 difference = log2(Abeta42/Abeta40)).
#'   Higher values = higher Abeta42/40 ratio = LESS amyloid pathology.
#'   Sign is preserved so that clinical interpretation is natural:
#'   lower A_score = more pathology.
#'
#' T_score: Tau. Mean of z-scored phospho-tau markers (BD and peripheral).
#'   Higher = more tau pathology.
#'
#' N_score: Neurodegeneration. Mean of z-scored neurodegeneration markers
#'   (NEFL, NEFH, GFAP, NRGN, SNAP25). Higher = more neurodegeneration.
#'
#' @param dat_list Output from load_alamar_data
#' @return dat_list with ATN scores added
calculate_atn_scores <- function(dat_list) {

  dat <- dat_list$data

  # A: Amyloid (log2 difference = log2 of ratio)
  # Higher A_42 - A_40 = higher Abeta42/40 = less pathology
  # Guard against A_40 == 0 (invalid on log2 scale)
  if (all(c("A_42", "A_40") %in% names(dat))) {
    dat <- dat %>%
      mutate(
        A_score = case_when(
          A_40 == 0 | is.na(A_40) | is.na(A_42) ~ NA_real_,
          TRUE ~ A_42 - A_40
        )
      )
    n_na <- sum(is.na(dat$A_score))
    if (n_na > 0) message(sprintf("  A_score: %d samples set to NA (A_40 = 0 or missing)", n_na))
  } else if ("ABeta42_minus_ABeta40" %in% names(dat)) {
    # Fallback: use pre-computed column if available
    dat$A_score <- dat$ABeta42_minus_ABeta40
  } else {
    warning("Cannot compute A_score: need A_42 and A_40 (or ABeta42_minus_ABeta40)")
  }

  # T: Tau -- z-score each marker then average (markers have different NPQ scales)
  tau_markers <- c("BD_pTau_181", "BD_pTau_217", "BD_pTau_231",
                   "pTau_181", "pTau_217", "pTau_231")
  available_tau <- intersect(tau_markers, names(dat))

  if (length(available_tau) > 0) {
    tau_z <- dat %>%
      select(all_of(available_tau)) %>%
      mutate(across(everything(), ~as.numeric(scale(.x))))
    dat$T_score <- rowMeans(tau_z, na.rm = TRUE)
    message(sprintf("  T_score: mean of %d z-scored tau markers", length(available_tau)))
  }

  # N: Neurodegeneration -- z-score each marker then average
  neurodegen_markers <- c("NFL", "NEFL", "NEFH", "GFAP", "NRGN", "SNAP25")
  available_neurodegen <- intersect(neurodegen_markers, names(dat))

  if (length(available_neurodegen) > 0) {
    neuro_z <- dat %>%
      select(all_of(available_neurodegen)) %>%
      mutate(across(everything(), ~as.numeric(scale(.x))))
    dat$N_score <- rowMeans(neuro_z, na.rm = TRUE)
    message(sprintf("  N_score: mean of %d z-scored neurodegeneration markers", length(available_neurodegen)))
  }

  message("Created ATN scores")

  dat_list$data <- dat
  return(dat_list)
}

#' Calculate brain-derived enrichment scores
#' 
#' @param dat_list Output from load_alamar_data
#' @return dat_list with BD enrichment scores
calculate_bd_enrichment <- function(dat_list) {
  
  dat <- dat_list$data
  
  # BD-tau / peripheral tau ratios
  if (all(c("BD_pTau_217", "pTau_217") %in% names(dat))) {
    dat <- dat %>%
      mutate(BD_pTau217_ratio = BD_pTau_217 - pTau_217)
  }
  
  if (all(c("BD_pTau_181", "pTau_181") %in% names(dat))) {
    dat <- dat %>%
      mutate(BD_pTau181_ratio = BD_pTau_181 - pTau_181)
  }
  
  if (all(c("BD_pTau_231", "pTau_231") %in% names(dat))) {
    dat <- dat %>%
      mutate(BD_pTau231_ratio = BD_pTau_231 - pTau_231)
  }
  
  message("Created brain-derived enrichment scores")
  
  dat_list$data <- dat
  return(dat_list)
}

#' Calculate inflammation scores via PCA
#' 
#' @param dat_list Output from load_alamar_data
#' @param n_pcs Number of principal components to extract
#' @return dat_list with inflammation PC scores
calculate_inflammation_pcs <- function(dat_list, n_pcs = 3) {
  
  dat <- dat_list$data
  
  # Identify cytokine/chemokine markers
  inflammation_markers <- c(
    "CCL2", "CCL3", "CCL4", "CCL11", "CCL13", "CCL17", "CCL22", "CCL26",
    "CXCL1", "CXCL8", "CXCL10", "CX3CL1",
    "IL1B", "IL2", "IL4", "IL5", "IL6", "IL7", "IL8", "IL9", "IL10",
    "IL12A", "IL12B", "IL12p70", "IL13", "IL15", "IL16", "IL17A", "IL18", "IL33",
    "TNF", "IFNG", "CSF2",
    "CRP", "TREM1", "TREM2"
  )
  
  available_inflam <- intersect(inflammation_markers, names(dat))
  message(sprintf("Found %d/%d inflammation markers", 
                  length(available_inflam), length(inflammation_markers)))
  
  if (length(available_inflam) < 5) {
    warning("Too few inflammation markers for PCA")
    return(dat_list)
  }
  
  # Extract inflammation data, remove missing
  inflam_data <- dat %>%
    select(all_of(available_inflam))
  
  complete_cases <- complete.cases(inflam_data)
  
  # Run PCA
  pca_result <- prcomp(inflam_data[complete_cases, ], 
                       scale. = TRUE, center = TRUE)
  
  # Extract PCs
  pc_scores <- pca_result$x[, 1:min(n_pcs, ncol(pca_result$x))]
  
  # Add back to data
  pc_df <- as.data.frame(pc_scores)
  names(pc_df) <- paste0("Inflammation_PC", 1:ncol(pc_df))
  
  # Initialize PC columns with NA
  for (col in names(pc_df)) {
    dat[[col]] <- NA_real_
  }
  
  # Fill in complete cases
  dat[complete_cases, names(pc_df)] <- pc_df
  
  message(sprintf("Created %d inflammation PCs (variance explained: %.1f%%)",
                  ncol(pc_df), 
                  sum(pca_result$sdev[1:ncol(pc_df)]^2) / sum(pca_result$sdev^2) * 100))
  
  dat_list$data <- dat
  return(dat_list)
}

#' Calculate multi-proteinopathy burden
#' 
#' @param dat_list Output from load_alamar_data
#' @return dat_list with proteinopathy scores
calculate_proteinopathy_burden <- function(dat_list) {
  
  dat <- dat_list$data
  
  # Standardize key pathology markers
  proteinopathy_markers <- list(
    amyloid = "A_score",
    tau = "BD_pTau_217",
    tdp43 = "pTDP43_409",
    synuclein = "Oligo_SNCA"
  )
  
  # Z-score each available marker
  for (marker_name in names(proteinopathy_markers)) {
    marker_col <- proteinopathy_markers[[marker_name]]
    
    if (marker_col %in% names(dat)) {
      z_col <- paste0(marker_name, "_z")
      dat[[z_col]] <- scale(dat[[marker_col]])[,1]
    }
  }
  
  # Sum z-scores (higher = more pathology)
  z_cols <- paste0(names(proteinopathy_markers), "_z")
  available_z <- intersect(z_cols, names(dat))
  
  if (length(available_z) > 0) {
    dat <- dat %>%
      mutate(
        Proteinopathy_burden = rowMeans(select(., all_of(available_z)), na.rm = TRUE)
      )
    
    message(sprintf("Created proteinopathy burden score from %d markers", 
                    length(available_z)))
  }
  
  dat_list$data <- dat
  return(dat_list)
}

#' Create all composite scores
#' 
#' @param dat_list Output from load_alamar_data
#' @return dat_list with all composite scores
create_all_composite_scores <- function(dat_list) {
  
  message("\n=== Creating composite scores ===")
  
  dat_list <- dat_list %>%
    calculate_atn_scores() %>%
    calculate_bd_enrichment() %>%
    calculate_inflammation_pcs(n_pcs = 3) %>%
    calculate_proteinopathy_burden()
  
  message("=== Composite scores complete ===\n")
  
  return(dat_list)
}
