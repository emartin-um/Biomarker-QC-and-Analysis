################################################################################
# 03_elastic_net.R
# Elastic net models for NCI vs AD classification (with covariate support)
################################################################################

library(tidyverse)
library(glmnet)
library(pROC)

#' Run elastic net NCI vs AD by ancestry
#'
#' @param dat_list Output from load_alamar_data
#' @param ancestry_group Ancestry to analyze ("All" or specific group)
#' @param alpha Elastic net mixing (1=LASSO, 0=Ridge)
#' @param nfolds Cross-validation folds
#' @param min_per_class Minimum samples per class (default 10)
#' @param include_covariates Character vector of covariate columns to include
#'   as predictors alongside biomarkers. Numeric covariates are used as-is;
#'   categorical covariates (character/factor) are dummy-coded.
#'   Default NULL = biomarkers only.
#' @return List with model results (or NULL if insufficient data)
run_elastic_net_by_ancestry <- function(dat_list,
                                        ancestry_group = "All",
                                        alpha = 0.5,
                                        nfolds = 10,
                                        min_per_class = 10,
                                        include_covariates = NULL) {

  dat <- dat_list$data

  # Filter to NCI vs AD
  dat_filtered <- dat %>%
    filter(CDX_collapsed %in% c("NCI", "AD"))

  # Filter ancestry if specified
  if (ancestry_group != "All") {
    dat_filtered <- dat_filtered %>%
      filter(Ancestry == ancestry_group)
  }

  # Check class sizes
  n_nci <- sum(dat_filtered$CDX_collapsed == "NCI")
  n_ad <- sum(dat_filtered$CDX_collapsed == "AD")

  message(sprintf("\n=== Elastic net for %s ===", ancestry_group))
  message(sprintf("N samples: %d (NCI=%d, AD=%d)", nrow(dat_filtered), n_nci, n_ad))

  # Check if we have enough samples
  if (n_nci < min_per_class || n_ad < min_per_class) {
    warning(sprintf("Insufficient samples for %s (need >=%d per class). Skipping.",
                    ancestry_group, min_per_class))
    return(NULL)
  }

  # Prepare biomarker matrix
  biomarker_mat <- dat_filtered %>%
    select(all_of(dat_list$biomarker_cols)) %>%
    as.matrix()

  # Remove columns with too many missing or zero variance
  col_keep <- colSums(is.na(biomarker_mat)) < nrow(biomarker_mat) * 0.5 &
    apply(biomarker_mat, 2, sd, na.rm = TRUE) > 0
  biomarker_mat <- biomarker_mat[, col_keep]

  # Track which biomarker columns survived filtering
  kept_biomarker_cols <- colnames(biomarker_mat)

  # Impute missing with column median
  for (j in 1:ncol(biomarker_mat)) {
    biomarker_mat[is.na(biomarker_mat[,j]), j] <- median(biomarker_mat[,j], na.rm = TRUE)
  }

  # Build covariate matrix if requested
  covariate_names <- character(0)
  if (!is.null(include_covariates) && length(include_covariates) > 0) {
    avail_covs <- intersect(include_covariates, names(dat_filtered))
    if (length(avail_covs) > 0) {
      cov_df <- dat_filtered[, avail_covs, drop = FALSE]

      # Dummy-code categorical variables; handle NAs and whitespace
      cov_mat_list <- list()
      for (cov in avail_covs) {
        if (is.numeric(cov_df[[cov]])) {
          vals <- cov_df[[cov]]
          vals[is.na(vals)] <- median(vals, na.rm = TRUE)
          cov_mat_list[[cov]] <- matrix(vals, ncol = 1, dimnames = list(NULL, cov))
        } else {
          # Normalize whitespace in string values (e.g., "3    3" → "3_3")
          cov_df[[cov]] <- gsub("\\s+", "_", trimws(cov_df[[cov]]))
          # Dummy code (drop first level)
          lvls <- sort(unique(na.omit(cov_df[[cov]])))
          if (length(lvls) >= 2) {
            for (k in 2:length(lvls)) {
              dummy_name <- paste0(cov, "_", lvls[k])
              cov_mat_list[[dummy_name]] <- matrix(
                as.numeric(!is.na(cov_df[[cov]]) & cov_df[[cov]] == lvls[k]),
                ncol = 1, dimnames = list(NULL, dummy_name)
              )
            }
          }
        }
      }

      if (length(cov_mat_list) > 0) {
        cov_mat <- do.call(cbind, cov_mat_list)
        covariate_names <- colnames(cov_mat)
        biomarker_mat <- cbind(biomarker_mat, cov_mat)
        message(sprintf("Added %d covariate features: %s",
                        length(covariate_names), paste(covariate_names, collapse = ", ")))
      }
    }
  }

  y <- as.numeric(dat_filtered$CDX_collapsed == "AD")

  # Run cross-validated elastic net
  set.seed(123)
  cv_fit <- cv.glmnet(biomarker_mat, y,
                      family = "binomial",
                      alpha = alpha,
                      nfolds = nfolds,
                      type.measure = "auc")

  # Get coefficients at lambda.min
  coefs <- coef(cv_fit, s = "lambda.min")
  coefs_df <- data.frame(
    biomarker = rownames(coefs),
    coefficient = as.numeric(coefs),
    stringsAsFactors = FALSE
  ) %>%
    filter(biomarker != "(Intercept)", coefficient != 0) %>%
    mutate(is_covariate = biomarker %in% covariate_names) %>%
    arrange(desc(abs(coefficient)))

  # Predictions
  pred_probs <- predict(cv_fit, newx = biomarker_mat, s = "lambda.min", type = "response")[,1]

  # ROC
  roc_obj <- roc(y, pred_probs, quiet = TRUE)

  message(sprintf("AUC: %.3f (95%% CI: %.3f-%.3f)",
                  as.numeric(auc(roc_obj)),
                  ci.auc(roc_obj)[1],
                  ci.auc(roc_obj)[3]))
  n_bio_selected <- sum(!coefs_df$is_covariate)
  n_cov_selected <- sum(coefs_df$is_covariate)
  message(sprintf("Selected %d biomarkers + %d covariates out of %d features",
                  n_bio_selected, n_cov_selected, ncol(biomarker_mat)))

  # Return results
  list(
    ancestry = ancestry_group,
    cv_fit = cv_fit,
    coefficients = coefs_df,
    predictions = pred_probs,
    truth = y,
    roc = roc_obj,
    auc = auc(roc_obj),
    auc_ci = ci.auc(roc_obj),
    n_samples = nrow(biomarker_mat),
    n_features_total = ncol(biomarker_mat),
    n_biomarkers_selected = n_bio_selected,
    n_covariates_selected = n_cov_selected,
    covariate_names = covariate_names,
    kept_biomarker_cols = kept_biomarker_cols,
    sample_ids = dat_filtered$Record_ID,
    include_covariates = include_covariates
  )
}

#' Run elastic net for all ancestries
#'
#' @param dat_list Output from load_alamar_data
#' @param ancestries Vector of ancestries to analyze
#' @param alpha Elastic net mixing
#' @param min_per_class Minimum samples per class
#' @param include_covariates Character vector of covariates to include (passed through)
#' @return List of results by ancestry (skips groups with insufficient data)
run_elastic_net_all_ancestries <- function(dat_list,
                                           ancestries = c("All", "Hispanic",
                                                         "African_American", "African"),
                                           alpha = 0.5,
                                           min_per_class = 10,
                                           include_covariates = NULL) {

  results <- list()

  for (anc in ancestries) {
    result <- tryCatch({
      run_elastic_net_by_ancestry(dat_list, anc, alpha = alpha,
                                  min_per_class = min_per_class,
                                  include_covariates = include_covariates)
    }, error = function(e) {
      warning(sprintf("Error for %s: %s", anc, e$message))
      NULL
    })

    # Only store non-NULL results
    if (!is.null(result)) {
      results[[anc]] <- result
    }
  }

  if (length(results) == 0) {
    stop("No ancestry groups had sufficient samples for analysis")
  }

  return(results)
}

#' Compare biomarkers across ancestries
#'
#' @param elastic_net_results Output from run_elastic_net_all_ancestries
#' @return Dataframe comparing selected biomarkers
compare_biomarkers_across_ancestries <- function(elastic_net_results) {

  # Extract top biomarkers from each ancestry
  all_biomarkers <- unique(unlist(lapply(elastic_net_results, function(x) {
    if (!is.null(x$coefficients)) x$coefficients$biomarker else NULL
  })))

  if (length(all_biomarkers) == 0) {
    warning("No biomarkers selected in any ancestry")
    return(data.frame())
  }

  # Create comparison matrix
  comparison <- data.frame(biomarker = all_biomarkers, stringsAsFactors = FALSE)

  for (anc in names(elastic_net_results)) {
    if (!is.null(elastic_net_results[[anc]]$coefficients)) {
      coefs <- elastic_net_results[[anc]]$coefficients
      comparison[[paste0(anc, "_coef")]] <- coefs$coefficient[match(all_biomarkers, coefs$biomarker)]
    }
  }

  # Count how many ancestries selected each biomarker
  coef_cols <- grep("_coef$", names(comparison), value = TRUE)
  if (length(coef_cols) > 0) {
    comparison$n_ancestries_selected <- rowSums(!is.na(comparison[, coef_cols, drop = FALSE]))
  } else {
    comparison$n_ancestries_selected <- 0
  }

  comparison <- comparison %>%
    arrange(desc(n_ancestries_selected))

  return(comparison)
}

#' Apply trained model to MCI samples
#'
#' @param dat_list Output from load_alamar_data
#' @param elastic_net_result Result from run_elastic_net_by_ancestry
#' @param ancestry_group Ancestry to apply model to
#' @return Dataframe with MCI predictions (or NULL if no model)
predict_mci_risk <- function(dat_list, elastic_net_result, ancestry_group = "All") {

  if (is.null(elastic_net_result)) {
    warning(sprintf("No model available for %s - cannot predict MCI risk", ancestry_group))
    return(NULL)
  }

  dat <- dat_list$data

  # Filter to MCI
  dat_mci <- dat %>%
    filter(CDX_collapsed == "MCI")

  if (ancestry_group != "All") {
    dat_mci <- dat_mci %>%
      filter(Ancestry == ancestry_group)
  }

  if (nrow(dat_mci) == 0) {
    warning(sprintf("No MCI samples for %s", ancestry_group))
    return(NULL)
  }

  message(sprintf("\nPredicting AD risk for %d MCI samples (%s)",
                  nrow(dat_mci), ancestry_group))

  # Prepare biomarker matrix -- must match the columns used during training
  # Use the kept_biomarker_cols from the model to align columns
  if (!is.null(elastic_net_result$kept_biomarker_cols)) {
    train_bio_cols <- elastic_net_result$kept_biomarker_cols
  } else {
    # Fallback for older results without kept_biomarker_cols
    train_bio_cols <- dat_list$biomarker_cols
  }

  biomarker_mat <- dat_mci %>%
    select(all_of(train_bio_cols)) %>%
    as.matrix()

  # Impute missing with column median
  for (j in 1:ncol(biomarker_mat)) {
    if (sum(is.na(biomarker_mat[,j])) > 0) {
      biomarker_mat[is.na(biomarker_mat[,j]), j] <- median(biomarker_mat[,j], na.rm = TRUE)
    }
  }

  # Add covariate columns if the model was trained with covariates
  if (!is.null(elastic_net_result$include_covariates) &&
      length(elastic_net_result$covariate_names) > 0) {

    avail_covs <- intersect(elastic_net_result$include_covariates, names(dat_mci))
    if (length(avail_covs) > 0) {
      cov_df <- dat_mci[, avail_covs, drop = FALSE]

      cov_mat_list <- list()
      for (cov in avail_covs) {
        if (is.numeric(cov_df[[cov]])) {
          vals <- cov_df[[cov]]
          vals[is.na(vals)] <- median(vals, na.rm = TRUE)
          cov_mat_list[[cov]] <- matrix(vals, ncol = 1, dimnames = list(NULL, cov))
        } else {
          # Normalize whitespace (must match training preprocessing)
          cov_df[[cov]] <- gsub("\\s+", "_", trimws(cov_df[[cov]]))
          lvls <- sort(unique(na.omit(cov_df[[cov]])))
          if (length(lvls) >= 2) {
            for (k in 2:length(lvls)) {
              dummy_name <- paste0(cov, "_", lvls[k])
              cov_mat_list[[dummy_name]] <- matrix(
                as.numeric(!is.na(cov_df[[cov]]) & cov_df[[cov]] == lvls[k]),
                ncol = 1, dimnames = list(NULL, dummy_name)
              )
            }
          }
        }
      }

      if (length(cov_mat_list) > 0) {
        cov_mat <- do.call(cbind, cov_mat_list)
        # Align to the covariate columns used in training
        train_cov_names <- elastic_net_result$covariate_names
        # Add any missing covariate columns as zeros
        for (cn in train_cov_names) {
          if (!cn %in% colnames(cov_mat)) {
            cov_mat <- cbind(cov_mat, matrix(0, nrow = nrow(cov_mat), ncol = 1,
                                              dimnames = list(NULL, cn)))
          }
        }
        cov_mat <- cov_mat[, train_cov_names, drop = FALSE]
        biomarker_mat <- cbind(biomarker_mat, cov_mat)
        message(sprintf("  Added %d covariate features for prediction", length(train_cov_names)))
      }
    }
  }

  # Predict
  pred_probs <- predict(elastic_net_result$cv_fit,
                       newx = biomarker_mat,
                       s = "lambda.min",
                       type = "response")[,1]

  # Create results dataframe
  results <- data.frame(
    Record_ID = dat_mci$Record_ID,
    Ancestry = dat_mci$Ancestry,
    age = dat_mci$age_at_subject,
    APOE4_carrier = dat_mci$APOE4_carrier,
    AD_risk_score = pred_probs
  )

  # Categorize risk
  results <- results %>%
    mutate(
      risk_category = cut(AD_risk_score,
                         breaks = c(0, 0.33, 0.66, 1),
                         labels = c("Low", "Moderate", "High"))
    )

  message(sprintf("Risk distribution: Low=%d, Moderate=%d, High=%d",
                  sum(results$risk_category == "Low"),
                  sum(results$risk_category == "Moderate"),
                  sum(results$risk_category == "High")))

  return(results)
}

#' Run elastic net for MCI+AD vs NCI (supervised feature selection)
#'
#' This analysis identifies biomarkers and covariates that best separate
#' MCI+AD (combined) from NCI. It is independent of the ATN framework and
#' answers: which features make MCI look most like AD and least like NCI?
#'
#' @param dat_list Output from load_alamar_data (preprocessed)
#' @param alpha Elastic net mixing (1=LASSO, 0=Ridge)
#' @param nfolds CV folds
#' @param include_covariates Character vector of covariates to include
#' @return List with model results
run_mci_ad_vs_nci <- function(dat_list,
                               alpha = 0.5,
                               nfolds = 10,
                               include_covariates = NULL) {

  dat <- dat_list$data

  # Create binary outcome: MCI+AD = 1, NCI = 0
  dat_filtered <- dat %>%
    filter(CDX_collapsed %in% c("NCI", "MCI", "AD")) %>%
    mutate(outcome = as.numeric(CDX_collapsed %in% c("MCI", "AD")))

  n_nci <- sum(dat_filtered$outcome == 0)
  n_mci_ad <- sum(dat_filtered$outcome == 1)

  message(sprintf("\n=== Elastic net: MCI+AD vs NCI ==="))
  message(sprintf("N samples: %d (NCI=%d, MCI+AD=%d)", nrow(dat_filtered), n_nci, n_mci_ad))

  # Prepare biomarker matrix
  biomarker_mat <- dat_filtered %>%
    select(all_of(dat_list$biomarker_cols)) %>%
    as.matrix()

  col_keep <- colSums(is.na(biomarker_mat)) < nrow(biomarker_mat) * 0.5 &
    apply(biomarker_mat, 2, sd, na.rm = TRUE) > 0
  biomarker_mat <- biomarker_mat[, col_keep]

  for (j in 1:ncol(biomarker_mat)) {
    biomarker_mat[is.na(biomarker_mat[,j]), j] <- median(biomarker_mat[,j], na.rm = TRUE)
  }

  # Build covariate matrix
  covariate_names <- character(0)
  if (!is.null(include_covariates) && length(include_covariates) > 0) {
    avail_covs <- intersect(include_covariates, names(dat_filtered))
    if (length(avail_covs) > 0) {
      cov_df <- dat_filtered[, avail_covs, drop = FALSE]
      cov_mat_list <- list()
      for (cov in avail_covs) {
        if (is.numeric(cov_df[[cov]])) {
          vals <- cov_df[[cov]]
          vals[is.na(vals)] <- median(vals, na.rm = TRUE)
          cov_mat_list[[cov]] <- matrix(vals, ncol = 1, dimnames = list(NULL, cov))
        } else {
          # Normalize whitespace in string values
          cov_df[[cov]] <- gsub("\\s+", "_", trimws(cov_df[[cov]]))
          lvls <- sort(unique(na.omit(cov_df[[cov]])))
          if (length(lvls) >= 2) {
            for (k in 2:length(lvls)) {
              dummy_name <- paste0(cov, "_", lvls[k])
              cov_mat_list[[dummy_name]] <- matrix(
                as.numeric(!is.na(cov_df[[cov]]) & cov_df[[cov]] == lvls[k]),
                ncol = 1, dimnames = list(NULL, dummy_name)
              )
            }
          }
        }
      }
      if (length(cov_mat_list) > 0) {
        cov_mat <- do.call(cbind, cov_mat_list)
        covariate_names <- colnames(cov_mat)
        biomarker_mat <- cbind(biomarker_mat, cov_mat)
        message(sprintf("Added %d covariate features: %s",
                        length(covariate_names), paste(covariate_names, collapse = ", ")))
      }
    }
  }

  y <- dat_filtered$outcome

  set.seed(123)
  cv_fit <- cv.glmnet(biomarker_mat, y,
                      family = "binomial",
                      alpha = alpha,
                      nfolds = nfolds,
                      type.measure = "auc")

  coefs <- coef(cv_fit, s = "lambda.min")
  coefs_df <- data.frame(
    feature = rownames(coefs),
    coefficient = as.numeric(coefs),
    stringsAsFactors = FALSE
  ) %>%
    filter(feature != "(Intercept)", coefficient != 0) %>%
    mutate(is_covariate = feature %in% covariate_names) %>%
    arrange(desc(abs(coefficient)))

  pred_probs <- predict(cv_fit, newx = biomarker_mat, s = "lambda.min", type = "response")[,1]
  roc_obj <- roc(y, pred_probs, quiet = TRUE)

  message(sprintf("AUC: %.3f (95%% CI: %.3f-%.3f)",
                  as.numeric(auc(roc_obj)),
                  ci.auc(roc_obj)[1],
                  ci.auc(roc_obj)[3]))

  # Per-diagnosis predicted scores for interpretation
  score_by_dx <- dat_filtered %>%
    mutate(pred_score = pred_probs) %>%
    group_by(CDX_collapsed) %>%
    summarize(
      n = n(),
      mean_score = mean(pred_score),
      sd_score = sd(pred_score),
      median_score = median(pred_score),
      .groups = "drop"
    )

  list(
    cv_fit = cv_fit,
    coefficients = coefs_df,
    predictions = pred_probs,
    truth = y,
    diagnosis = dat_filtered$CDX_collapsed,
    roc = roc_obj,
    auc = auc(roc_obj),
    auc_ci = ci.auc(roc_obj),
    n_samples = nrow(biomarker_mat),
    n_features_total = ncol(biomarker_mat),
    score_by_diagnosis = score_by_dx,
    covariate_names = covariate_names
  )
}
