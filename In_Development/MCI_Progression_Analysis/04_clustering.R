################################################################################
# 04_clustering.R
# MCI clustering and subtyping
################################################################################

library(tidyverse)
library(cluster)
library(factoextra)

#' Cluster MCI samples by biomarkers (and optionally covariates)
#'
#' @param dat_list Output from load_alamar_data
#' @param n_clusters Number of clusters (or "auto" for gap statistic)
#' @param biomarker_subset Optional vector of biomarkers to use (NULL = use top varying)
#' @param top_n If biomarker_subset is NULL, use top N most variable biomarkers
#' @param include_covariates Character vector of covariate column names to include
#'   in clustering (e.g., c("age_at_subject", "sex", "APOE4_carrier")).
#'   Numeric covariates are median-imputed; categorical covariates are dummy-coded.
#'   All features are equally weighted after scaling.
#' @return List with clustering results
cluster_mci_samples <- function(dat_list,
                                n_clusters = "auto",
                                biomarker_subset = NULL,
                                top_n = 50,
                                include_covariates = NULL) {

  dat <- dat_list$data

  # Filter to MCI only
  dat_mci <- dat %>%
    filter(CDX_collapsed == "MCI")

  cat(sprintf("\n=== Clustering %d MCI samples ===", nrow(dat_mci)))

  # Select biomarkers
  if (is.null(biomarker_subset)) {
    # Use top N most variable biomarkers
    biomarker_mat <- dat_mci %>%
      select(all_of(dat_list$biomarker_cols)) %>%
      as.matrix()

    biomarker_sd <- apply(biomarker_mat, 2, sd, na.rm = TRUE)
    top_biomarkers <- names(sort(biomarker_sd, decreasing = TRUE)[1:min(top_n, length(biomarker_sd))])

   # cat(sprintf("\nUsing top %d most variable biomarkers", length(top_biomarkers)))
  } else {
    top_biomarkers <- biomarker_subset
   # cat(sprintf("\nUsing %d specified biomarkers", length(top_biomarkers)))
  }

  # Extract biomarker matrix
  biomarker_mat <- dat_mci %>%
    select(all_of(top_biomarkers)) %>%
    as.matrix()

  # Remove samples with too many missing
  complete_enough <- rowSums(is.na(biomarker_mat)) < ncol(biomarker_mat) * 0.3
  biomarker_mat <- biomarker_mat[complete_enough, ]
  dat_mci <- dat_mci[complete_enough, ]

  # Impute remaining missing with column median
  for (j in 1:ncol(biomarker_mat)) {
    if (sum(is.na(biomarker_mat[,j])) > 0) {
      biomarker_mat[is.na(biomarker_mat[,j]), j] <- median(biomarker_mat[,j], na.rm = TRUE)
    }
  }

  # Add covariate features if requested (same pattern as 03_elastic_net.R)
  covariate_names <- character(0)
  if (!is.null(include_covariates) && length(include_covariates) > 0) {
    avail_covs <- intersect(include_covariates, names(dat_mci))
    if (length(avail_covs) > 0) {
      cov_df <- dat_mci[, avail_covs, drop = FALSE]
      cov_mat_list <- list()

      for (cov in avail_covs) {
        if (is.numeric(cov_df[[cov]])) {
          # Numeric: median imputation
          vals <- cov_df[[cov]]
          vals[is.na(vals)] <- median(vals, na.rm = TRUE)
          cov_mat_list[[cov]] <- matrix(vals, ncol = 1, dimnames = list(NULL, cov))
        } else {
          # Categorical: dummy coding (drop first level)
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
        cat(sprintf("\nAdded %d covariate feature(s) for clustering: %s\n",
                    length(covariate_names), paste(covariate_names, collapse = ", ")))
      }
    }
  }

  # Track all feature names (biomarkers + covariates) to keep in sync with matrix
  all_feature_names <- colnames(biomarker_mat)

  # Scale
  biomarker_mat_scaled <- scale(biomarker_mat)

  # Remove zero-variance columns (scale() produces NaN for these)
  nan_cols <- apply(biomarker_mat_scaled, 2, function(x) any(is.nan(x)))
  if (any(nan_cols)) {
    removed_names <- all_feature_names[nan_cols]
    cat(sprintf("\nRemoved %d zero-variance feature(s) after scaling: %s\n",
                sum(nan_cols), paste(removed_names, collapse = ", ")))
    biomarker_mat_scaled <- biomarker_mat_scaled[, !nan_cols, drop = FALSE]
    all_feature_names <- all_feature_names[!nan_cols]
    # Update biomarker and covariate name vectors
    top_biomarkers <- intersect(top_biomarkers, all_feature_names)
    covariate_names <- intersect(covariate_names, all_feature_names)
  }

  # Store scaling parameters for downstream use
  scale_center <- attr(biomarker_mat_scaled, "scaled:center")
  scale_sd <- attr(biomarker_mat_scaled, "scaled:scale")

  # Biomarker-only scaling params (for AD distance calculation, which
  # only has biomarker columns, not covariate columns)
  bio_idx <- which(all_feature_names %in% top_biomarkers)
  scale_center_bio <- scale_center[bio_idx]
  scale_sd_bio <- scale_sd[bio_idx]

  # Determine optimal number of clusters
  auto_k <- identical(n_clusters, "auto")
  if (auto_k) {
    #cat("\nDetermining optimal number of clusters...")

    # Gap statistic
    gap_stat <- clusGap(biomarker_mat_scaled,
                       FUN = kmeans,
                       K.max = 6,
                       B = 50,
                       nstart = 25,
                       iter.max = 100)

    n_clusters <- maxSE(gap_stat$Tab[,"gap"], gap_stat$Tab[,"SE.sim"])
    # Ensure at least 2 clusters (silhouette undefined for k=1)
    n_clusters <- max(n_clusters, 2)
    #cat(sprintf("\nOptimal clusters by gap statistic: %d", n_clusters))
  }

  # Perform clustering
  set.seed(123)
  km_result <- kmeans(biomarker_mat_scaled,
                     centers = n_clusters,
                     nstart = 50,
                     iter.max = 100)

  # Silhouette
  d <- dist(biomarker_mat_scaled)
  sil <- silhouette(km_result$cluster, d)
  sil_mat <- as.matrix(sil)  # coerce from silhouette class to plain matrix
  avg_sil <- mean(sil_mat[, "sil_width"])

 # cat(sprintf("\nK-means clustering: %d clusters, avg silhouette = %.3f",
 #                n_clusters, avg_sil))

  # Add cluster assignments to data
  dat_mci$cluster <- km_result$cluster

  # Cluster sizes
  #cat("\nCluster sizes:")
  #knitr::kable(as.data.frame(table(dat_mci$cluster)))

  # Return results
  list(
    data = dat_mci,
    cluster_assignments = km_result$cluster,
    centers = km_result$centers,
    silhouette = sil,
    avg_silhouette = avg_sil,
    n_clusters = n_clusters,
    biomarkers_used = top_biomarkers,
    covariate_names = covariate_names,
    biomarker_mat_scaled = biomarker_mat_scaled,
    scale_center = scale_center,
    scale_sd = scale_sd,
    scale_center_bio = scale_center_bio,
    scale_sd_bio = scale_sd_bio,
    gap_stat = if (auto_k) gap_stat else NULL
  )
}

#' Profile clusters by biomarkers and covariates
#'
#' @param cluster_result Output from cluster_mci_samples
#' @param dat_list Original dat_list for biomarker names
#' @return List with profile dataframes (biomarker, age, apoe, ancestry)
profile_mci_clusters <- function(cluster_result, dat_list) {

  dat_mci <- cluster_result$data

  cat("\n=== Profiling MCI clusters ===")

  # Biomarker profiles (mean by cluster)
  biomarker_profiles <- dat_mci %>%
    select(cluster, all_of(cluster_result$biomarkers_used)) %>%
    group_by(cluster) %>%
    summarize(across(everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop")

  # Covariate profiles
  age_profile <- NULL
  if ("age_at_subject" %in% names(dat_mci)) {
    age_profile <- dat_mci %>%
      group_by(cluster) %>%
      summarize(n = n(),
                mean_age = round(mean(age_at_subject, na.rm = TRUE), 1),
                sd_age = round(sd(age_at_subject, na.rm = TRUE), 1),
                .groups = "drop")
  }

  apoe_profile <- NULL
  if ("APOE4_carrier" %in% names(dat_mci)) {
    apoe_profile <- dat_mci %>%
      count(cluster, APOE4_carrier) %>%
      group_by(cluster) %>%
      mutate(pct = round(n / sum(n) * 100, 1)) %>%
      ungroup()
  }

  ancestry_profile <- NULL
  if ("Ancestry" %in% names(dat_mci)) {
    ancestry_profile <- dat_mci %>%
      count(cluster, Ancestry) %>%
      group_by(cluster) %>%
      mutate(pct = round(n / sum(n) * 100, 1)) %>%
      ungroup()
  }

  list(
    biomarker_profiles = biomarker_profiles,
    age_profile = age_profile,
    apoe_profile = apoe_profile,
    ancestry_profile = ancestry_profile
  )
}

#' Calculate distance of MCI clusters to AD centroid
#'
#' @param dat_list Output from load_alamar_data (with MCI cluster assignments)
#' @param cluster_result Output from cluster_mci_samples
#' @param ancestry_group Ancestry for AD centroid ("All" or specific)
#' @return Dataframe with distances
calculate_cluster_ad_distance <- function(dat_list, cluster_result, ancestry_group = "All") {

  dat <- dat_list$data
  dat_mci <- cluster_result$data

  # Get AD samples
  dat_ad <- dat %>%
    filter(CDX_collapsed == "AD")

  if (ancestry_group != "All") {
    dat_ad <- dat_ad %>%
      filter(Ancestry == ancestry_group)
  }

  cat(sprintf("\nCalculating distance to AD centroid (%s, N=%d)",
                  ancestry_group, nrow(dat_ad)))

  # Extract biomarkers used in clustering
  biomarkers <- cluster_result$biomarkers_used

  # AD centroid
  ad_biomarkers <- dat_ad %>%
    select(all_of(biomarkers)) %>%
    as.matrix()

  # Remove missing
  ad_complete <- complete.cases(ad_biomarkers)
  ad_biomarkers <- ad_biomarkers[ad_complete, ]

  # Scale AD biomarkers using biomarker-only MCI scaling parameters
  # (scale_center_bio/scale_sd_bio exclude covariate columns that may be in the full params)
  if (!is.null(cluster_result$scale_center_bio) && !is.null(cluster_result$scale_sd_bio)) {
    ad_biomarkers <- scale(ad_biomarkers,
                           center = cluster_result$scale_center_bio,
                           scale = cluster_result$scale_sd_bio)
  } else if (!is.null(cluster_result$scale_center) && !is.null(cluster_result$scale_sd)) {
    # Fallback for results without covariate support
    ad_biomarkers <- scale(ad_biomarkers,
                           center = cluster_result$scale_center,
                           scale = cluster_result$scale_sd)
  }

  ad_centroid <- colMeans(ad_biomarkers, na.rm = TRUE)

  # MCI cluster centroids — extract only biomarker columns
  # (cluster centers may include covariate columns if covariates were used)
  biomarker_col_idx <- which(colnames(cluster_result$centers) %in% biomarkers)
  cluster_centroids <- cluster_result$centers[, biomarker_col_idx, drop = FALSE]

  # Calculate Euclidean distance of each cluster centroid to AD centroid
  distances <- apply(cluster_centroids, 1, function(cluster_center) {
    sqrt(sum((cluster_center - ad_centroid)^2))
  })

  distance_df <- data.frame(
    cluster = 1:cluster_result$n_clusters,
    distance_to_AD = distances
  ) %>%
    arrange(distance_to_AD)

  cat("Cluster distances to AD centroid (closer = more AD-like):")
  print(distance_df)

  distance_df
}

#' Dimensional reduction for visualization (UMAP/PCA)
#'
#' Computes PCA (or UMAP) on the MCI biomarker data and optionally projects
#' NCI and AD samples into the same PC space for comparison.
#'
#' @param cluster_result Output from cluster_mci_samples
#' @param dat_list Full dat_list (needed to project NCI/AD)
#' @param method "pca" or "umap"
#' @param n_pcs Number of PC dimensions to store (default 5)
#' @param project_all If TRUE, also project NCI and AD into the MCI PCA space
#' @return cluster_result with dim reduction coords; if project_all,
#'   also includes $all_projected with all samples in the same space
add_dim_reduction <- function(cluster_result, dat_list = NULL,
                               method = "pca", n_pcs = 5,
                               project_all = FALSE) {

  dat_mci <- cluster_result$data
  biomarkers <- cluster_result$biomarkers_used

  # Extract biomarker matrix
  biomarker_mat <- dat_mci %>%
    select(all_of(biomarkers)) %>%
    as.matrix()

  # Impute any remaining NAs
  for (j in 1:ncol(biomarker_mat)) {
    if (sum(is.na(biomarker_mat[,j])) > 0) {
      biomarker_mat[is.na(biomarker_mat[,j]), j] <- median(biomarker_mat[,j], na.rm = TRUE)
    }
  }

  # Scale
  biomarker_mat_scaled <- scale(biomarker_mat)

  # Remove zero-variance columns (scale() produces NaN for these)
  nan_cols <- apply(biomarker_mat_scaled, 2, function(x) any(is.nan(x)))
  if (any(nan_cols)) {
    cat(sprintf("\nRemoved %d zero-variance biomarker(s) in dim reduction\n", sum(nan_cols)))
    biomarker_mat_scaled <- biomarker_mat_scaled[, !nan_cols, drop = FALSE]
    biomarkers <- biomarkers[!nan_cols]
  }

  center_vals <- attr(biomarker_mat_scaled, "scaled:center")
  scale_vals  <- attr(biomarker_mat_scaled, "scaled:scale")

  if (method == "pca") {
    pca_result <- prcomp(biomarker_mat_scaled, center = FALSE, scale. = FALSE)

    n_pcs_avail <- min(n_pcs, ncol(pca_result$x))
    for (k in 1:n_pcs_avail) {
      dat_mci[[paste0("PC", k)]] <- pca_result$x[, k]
    }
    dat_mci$dim1 <- pca_result$x[, 1]
    dat_mci$dim2 <- pca_result$x[, 2]

    var_explained <- summary(pca_result)$importance[2, 1:min(n_pcs_avail, ncol(pca_result$x))] * 100
    cat(sprintf("PCA: %s variance explained",
                    paste(sprintf("PC%d=%.1f%%", 1:length(var_explained), var_explained),
                          collapse = ", ")))

    cluster_result$pca_result <- pca_result
    cluster_result$var_explained <- var_explained

    # Project NCI and AD into the MCI PCA space
    if (project_all && !is.null(dat_list)) {
      dat_all <- dat_list$data %>%
        filter(CDX_collapsed %in% c("NCI", "MCI", "AD"))

      all_bio_mat <- dat_all %>%
        select(all_of(biomarkers)) %>%
        as.matrix()

      # Impute
      for (j in 1:ncol(all_bio_mat)) {
        if (sum(is.na(all_bio_mat[,j])) > 0) {
          all_bio_mat[is.na(all_bio_mat[,j]), j] <- median(all_bio_mat[,j], na.rm = TRUE)
        }
      }

      # Scale using MCI parameters
      all_scaled <- scale(all_bio_mat, center = center_vals, scale = scale_vals)

      # Project
      all_pcs <- all_scaled %*% pca_result$rotation

      all_projected <- dat_all %>%
        select(Record_ID, CDX_collapsed, Ancestry,
               any_of(c("age_at_subject", "sex", "APOE4_carrier")))
      for (k in 1:n_pcs_avail) {
        all_projected[[paste0("PC", k)]] <- all_pcs[, k]
      }
      all_projected$dim1 <- all_pcs[, 1]
      all_projected$dim2 <- all_pcs[, 2]

      # Add MCI cluster labels where applicable
      all_projected$cluster <- NA_integer_
      mci_idx <- which(all_projected$CDX_collapsed == "MCI")
      if (length(mci_idx) > 0) {
        # Match by Record_ID
        cluster_lookup <- setNames(dat_mci$cluster, dat_mci$Record_ID)
        all_projected$cluster[mci_idx] <- cluster_lookup[all_projected$Record_ID[mci_idx]]
      }

      cluster_result$all_projected <- all_projected
      cat(sprintf("Projected %d NCI + MCI + AD samples into MCI PCA space",
                      nrow(all_projected)))
    }

  } else if (method == "umap") {
    if (!requireNamespace("umap", quietly = TRUE)) {
      stop("umap package required for UMAP")
    }

    umap_result <- umap::umap(biomarker_mat_scaled)

    dat_mci$dim1 <- umap_result$layout[, 1]
    dat_mci$dim2 <- umap_result$layout[, 2]

    cat("UMAP coordinates calculated")
  }

  cluster_result$data <- dat_mci
  cluster_result$dim_reduction_method <- method

  return(cluster_result)
}


#' Supervised PCA: define PCA space from AD+NCI, project MCI
#'
#' Instead of computing PCA on MCI (which may produce axes dominated by
#' ancestry or noise), this function defines PCA axes from the known
#' disease contrast (AD + NCI samples). MCI samples are then projected
#' into this reference space, giving each PC direct clinical meaning:
#' position along each axis reflects similarity to AD vs NCI.
#'
#' @param dat_list Full dat_list with all samples
#' @param cluster_result Output from cluster_mci_samples (for biomarker list and cluster labels)
#' @param n_pcs Number of PCs to retain (default 10)
#' @param n_clusters Number of clusters for MCI in projected space ("auto" or integer)
#' @return List with supervised PCA results, including projected coordinates,
#'   cluster assignments, and per-PC AD-direction analysis
supervised_pca_clustering <- function(dat_list, cluster_result,
                                      n_pcs = 10, n_clusters = "auto") {

  biomarkers <- cluster_result$biomarkers_used
  dat <- dat_list$data

  # ---- 1. Build reference PCA from AD + NCI ----
  dat_ref <- dat %>%
    filter(CDX_collapsed %in% c("NCI", "AD"))

  ref_mat <- dat_ref %>%
    select(all_of(biomarkers)) %>%
    as.matrix()

  # Impute missing with column median
  for (j in 1:ncol(ref_mat)) {
    na_idx <- is.na(ref_mat[, j])
    if (any(na_idx)) ref_mat[na_idx, j] <- median(ref_mat[, j], na.rm = TRUE)
  }

  # Scale reference data
  ref_scaled <- scale(ref_mat)

  # Remove zero-variance columns
  nan_cols <- apply(ref_scaled, 2, function(x) any(is.nan(x)))
  if (any(nan_cols)) {
    ref_scaled <- ref_scaled[, !nan_cols, drop = FALSE]
    biomarkers_kept <- biomarkers[!nan_cols]
  } else {
    biomarkers_kept <- biomarkers
  }

  ref_center <- attr(ref_scaled, "scaled:center")
  ref_scale  <- attr(ref_scaled, "scaled:scale")

  # PCA on reference
  pca_ref <- prcomp(ref_scaled, center = FALSE, scale. = FALSE)
  n_pcs_avail <- min(n_pcs, ncol(pca_ref$x))

  var_explained <- summary(pca_ref)$importance[2, 1:n_pcs_avail] * 100

  cat(sprintf("\n=== Supervised PCA (AD+NCI reference, N=%d) ===\n", nrow(dat_ref)))
  cat(sprintf("Variance explained: %s\n",
              paste(sprintf("PC%d=%.1f%%", 1:length(var_explained), var_explained),
                    collapse = ", ")))

  # ---- 2. Project all groups (NCI, MCI, AD) into reference PCA space ----
  dat_all <- dat %>%
    filter(CDX_collapsed %in% c("NCI", "MCI", "AD"))

  all_mat <- dat_all %>%
    select(all_of(biomarkers_kept)) %>%
    as.matrix()

  # Impute
  for (j in 1:ncol(all_mat)) {
    na_idx <- is.na(all_mat[, j])
    if (any(na_idx)) all_mat[na_idx, j] <- median(all_mat[, j], na.rm = TRUE)
  }

  # Scale using reference parameters
  all_scaled <- scale(all_mat, center = ref_center, scale = ref_scale)

  # Project
  all_pcs <- all_scaled %*% pca_ref$rotation[, 1:n_pcs_avail]

  projected <- dat_all %>%
    select(Record_ID, CDX_collapsed, any_of(c("Ancestry", "age_at_subject",
                                               "sex", "APOE4_carrier")))

  for (k in 1:n_pcs_avail) {
    projected[[paste0("PC", k)]] <- all_pcs[, k]
  }

  # Add cluster labels for MCI samples from the original clustering
  projected$cluster <- NA_integer_
  mci_idx <- which(projected$CDX_collapsed == "MCI")
  if (length(mci_idx) > 0) {
    cluster_lookup <- setNames(cluster_result$data$cluster,
                                cluster_result$data$Record_ID)
    projected$cluster[mci_idx] <- cluster_lookup[projected$Record_ID[mci_idx]]
  }

  # ---- 3. Re-cluster MCI in supervised PCA space ----
  mci_projected <- projected %>% filter(CDX_collapsed == "MCI")
  mci_pc_mat <- as.matrix(mci_projected[, paste0("PC", 1:n_pcs_avail)])

  # Remove rows with NaN (shouldn't happen, but safety)
  valid_rows <- complete.cases(mci_pc_mat)
  mci_pc_mat <- mci_pc_mat[valid_rows, ]
  mci_projected <- mci_projected[valid_rows, ]

  if (identical(n_clusters, "auto")) {
    gap_stat <- clusGap(mci_pc_mat, FUN = kmeans, K.max = 6, B = 50, nstart = 25, iter.max = 100)
    opt_k <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"])
    opt_k <- max(opt_k, 2)
  } else {
    opt_k <- n_clusters
  }

  set.seed(123)
  km_sup <- kmeans(mci_pc_mat, centers = opt_k, nstart = 50, iter.max = 100)
  mci_projected$supervised_cluster <- km_sup$cluster

  cat(sprintf("MCI re-clustered in supervised PCA space: %d clusters\n", opt_k))

  # ---- 4. Per-PC AD-direction analysis ----
  # For each PC, compute mean position of NCI, AD, and each MCI cluster
  nci_means <- projected %>%
    filter(CDX_collapsed == "NCI") %>%
    summarize(across(starts_with("PC"), ~mean(.x, na.rm = TRUE))) %>%
    unlist()

  ad_means <- projected %>%
    filter(CDX_collapsed == "AD") %>%
    summarize(across(starts_with("PC"), ~mean(.x, na.rm = TRUE))) %>%
    unlist()

  # AD direction on each PC
  ad_direction <- ad_means - nci_means  # positive = AD side

  # Where does each MCI cluster fall relative to NCI?
  cluster_pc_means <- mci_projected %>%
    group_by(supervised_cluster) %>%
    summarize(across(starts_with("PC"), ~mean(.x, na.rm = TRUE)), .groups = "drop")

  # Convert to fraction of NCI-to-AD distance
  pc_cols <- paste0("PC", 1:n_pcs_avail)
  ad_fraction <- cluster_pc_means
  for (pc in pc_cols) {
    nci_to_ad <- ad_means[pc] - nci_means[pc]
    if (abs(nci_to_ad) > 1e-10) {
      ad_fraction[[pc]] <- (cluster_pc_means[[pc]] - nci_means[pc]) / nci_to_ad
    } else {
      ad_fraction[[pc]] <- NA_real_
    }
  }

  cat("\nMCI cluster positions (fraction of NCI-to-AD distance, 0=NCI, 1=AD):\n")

  list(
    projected = projected,
    mci_projected = mci_projected,
    pca_ref = pca_ref,
    var_explained = var_explained,
    n_pcs = n_pcs_avail,
    supervised_clusters = km_sup,
    n_clusters = opt_k,
    nci_means = nci_means,
    ad_means = ad_means,
    ad_direction = ad_direction,
    ad_fraction = ad_fraction,
    ref_center = ref_center,
    ref_scale = ref_scale,
    biomarkers_used = biomarkers_kept
  )
}
