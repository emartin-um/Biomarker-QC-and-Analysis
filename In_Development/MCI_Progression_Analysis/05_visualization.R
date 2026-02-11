################################################################################
# 05_visualization.R
# Plotting utilities
################################################################################

library(tidyverse)
library(ggplot2)
library(patchwork)
library(pheatmap)

#' Plot ROC curves by ancestry
#'
#' @param elastic_net_results Output from run_elastic_net_all_ancestries
#' @return ggplot object
plot_roc_by_ancestry <- function(elastic_net_results) {

  # Combine ROC data
  roc_data <- map_dfr(names(elastic_net_results), function(anc) {
    if (!is.null(elastic_net_results[[anc]]$roc)) {
      roc_obj <- elastic_net_results[[anc]]$roc

      data.frame(
        ancestry = anc,
        fpr = 1 - roc_obj$specificities,
        tpr = roc_obj$sensitivities,
        auc = as.numeric(elastic_net_results[[anc]]$auc)
      )
    }
  })

  # Create ancestry labels with AUC
  auc_labels <- roc_data %>%
    distinct(ancestry, auc) %>%
    mutate(label = sprintf("%s (AUC=%.3f)", ancestry, auc))
  roc_data <- roc_data %>%
    left_join(auc_labels %>% select(ancestry, label), by = "ancestry")

  # Plot
  p <- ggplot(roc_data, aes(x = fpr, y = tpr, color = label)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    facet_wrap(~label, nrow = 2) +
    labs(
      title = "ROC Curves: NCI vs AD Classification",
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  return(p)
}

#' Plot top biomarkers by ancestry (with covariate flag)
#'
#' @param elastic_net_results Output from run_elastic_net_all_ancestries
#' @param top_n Number of top biomarkers per ancestry
#' @return ggplot object
plot_top_biomarkers <- function(elastic_net_results, top_n = 20) {

  # Extract top biomarkers
  biomarker_data <- map_dfr(names(elastic_net_results), function(anc) {
    if (!is.null(elastic_net_results[[anc]]$coefficients)) {
      elastic_net_results[[anc]]$coefficients %>%
        slice_max(abs(coefficient), n = top_n) %>%
        mutate(ancestry = anc)
    }
  })

  if (nrow(biomarker_data) == 0) return(NULL)

  # Color by direction, shape by covariate
  if ("is_covariate" %in% names(biomarker_data)) {
    biomarker_data <- biomarker_data %>%
      mutate(feature_type = ifelse(is_covariate, "Covariate", "Biomarker"))
  } else {
    biomarker_data$feature_type <- "Biomarker"
  }

  # Plot
  p <- ggplot(biomarker_data, aes(x = reorder(biomarker, coefficient),
                                   y = coefficient,
                                   fill = coefficient > 0,
                                   alpha = feature_type)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ancestry, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                      guide = "none") +
    scale_alpha_manual(values = c("Biomarker" = 1, "Covariate" = 0.5),
                       name = "Feature Type") +
    labs(
      title = sprintf("Top %d Features per Ancestry", top_n),
      subtitle = "Red = higher in AD | Blue = higher in NCI | Faded = covariate",
      x = "Feature",
      y = "Elastic Net Coefficient"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  return(p)
}

#' Plot MCI cluster visualization
#'
#' @param cluster_result Output from cluster_mci_samples (with dim reduction)
#' @param color_by Variable to color points by (default: "cluster")
#' @return ggplot object
plot_mci_clusters <- function(cluster_result, color_by = "cluster") {

  dat_mci <- cluster_result$data

  if (!all(c("dim1", "dim2") %in% names(dat_mci))) {
    stop("Run add_dim_reduction first")
  }

  method_label <- ifelse(cluster_result$dim_reduction_method == "pca",
                        "Principal Component", "UMAP")

  # Add variance explained to axis labels if PCA
  if (!is.null(cluster_result$var_explained) &&
      cluster_result$dim_reduction_method == "pca") {
    xlab <- sprintf("PC1 (%.1f%%)", cluster_result$var_explained[1])
    ylab <- sprintf("PC2 (%.1f%%)", cluster_result$var_explained[2])
  } else {
    xlab <- paste(method_label, "1")
    ylab <- paste(method_label, "2")
  }

  p <- ggplot(dat_mci, aes(x = dim1, y = dim2, color = as.factor(get(color_by)))) +
    geom_point(size = 2, alpha = 0.7) +
    labs(
      title = sprintf("MCI Clusters (%s)", method_label),
      x = xlab,
      y = ylab,
      color = color_by
    ) +
    theme_bw() +
    theme(legend.position = "right")

  return(p)
}

#' Plot all diagnostic groups in PCA space (NCI, MCI, AD)
#'
#' @param cluster_result Output from add_dim_reduction with project_all=TRUE
#' @return ggplot object
plot_all_diagnoses_pca <- function(cluster_result) {

  if (is.null(cluster_result$all_projected)) {
    stop("Run add_dim_reduction with project_all=TRUE and dat_list first")
  }

  dat_all <- cluster_result$all_projected

  # Add variance explained to axis labels
  if (!is.null(cluster_result$var_explained)) {
    xlab <- sprintf("PC1 (%.1f%%)", cluster_result$var_explained[1])
    ylab <- sprintf("PC2 (%.1f%%)", cluster_result$var_explained[2])
  } else {
    xlab <- "PC1"
    ylab <- "PC2"
  }

  # Order factors
  dat_all$CDX_collapsed <- factor(dat_all$CDX_collapsed, levels = c("NCI", "MCI", "AD"))

  p <- ggplot(dat_all, aes(x = dim1, y = dim2, color = CDX_collapsed)) +
    geom_point(size = 1.5, alpha = 0.5) +
    stat_ellipse(level = 0.68, linewidth = 1, linetype = "solid") +
    scale_color_manual(values = c("NCI" = "#4DAF4A", "MCI" = "#FF7F00", "AD" = "#E41A1C")) +
    labs(
      title = "All Diagnostic Groups in MCI PCA Space",
      subtitle = "PCA trained on MCI samples; NCI and AD projected into same space",
      x = xlab, y = ylab,
      color = "Diagnosis"
    ) +
    theme_bw() +
    theme(legend.position = "right")

  return(p)
}

#' Plot higher PC pairs
#'
#' @param cluster_result Output from add_dim_reduction
#' @param pc_x PC number for x-axis
#' @param pc_y PC number for y-axis
#' @param color_by Column to color by
#' @return ggplot object
plot_pc_pair <- function(cluster_result, pc_x = 1, pc_y = 2, color_by = "cluster") {

  dat <- cluster_result$data
  pc_x_col <- paste0("PC", pc_x)
  pc_y_col <- paste0("PC", pc_y)

  if (!all(c(pc_x_col, pc_y_col) %in% names(dat))) {
    stop(sprintf("PC%d and/or PC%d not available. Increase n_pcs in add_dim_reduction.", pc_x, pc_y))
  }

  xlab <- sprintf("PC%d", pc_x)
  ylab <- sprintf("PC%d", pc_y)
  if (!is.null(cluster_result$var_explained)) {
    if (pc_x <= length(cluster_result$var_explained))
      xlab <- sprintf("PC%d (%.1f%%)", pc_x, cluster_result$var_explained[pc_x])
    if (pc_y <= length(cluster_result$var_explained))
      ylab <- sprintf("PC%d (%.1f%%)", pc_y, cluster_result$var_explained[pc_y])
  }

  p <- ggplot(dat, aes(x = .data[[pc_x_col]], y = .data[[pc_y_col]],
                        color = as.factor(.data[[color_by]]))) +
    geom_point(size = 2, alpha = 0.7) +
    labs(title = sprintf("MCI Clusters: %s vs %s", xlab, ylab),
         x = xlab, y = ylab, color = color_by) +
    theme_bw() +
    theme(legend.position = "right")

  return(p)
}

#' Plot cluster profiles heatmap
#'
#' @param cluster_result Output from cluster_mci_samples
#' @param top_n Number of top varying biomarkers to show
#' @return pheatmap object
plot_cluster_heatmap <- function(cluster_result, top_n = 30) {

  dat_mci <- cluster_result$data
  biomarkers <- cluster_result$biomarkers_used

  # Calculate cluster means
  cluster_means <- dat_mci %>%
    select(cluster, all_of(biomarkers)) %>%
    group_by(cluster) %>%
    summarize(across(everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop")

  # Select top varying biomarkers across clusters
  cluster_mat <- as.matrix(cluster_means[, -1])
  rownames(cluster_mat) <- paste0("Cluster ", cluster_means$cluster)

  biomarker_var <- apply(cluster_mat, 2, var, na.rm = TRUE)
  top_biomarkers <- names(sort(biomarker_var, decreasing = TRUE)[1:min(top_n, length(biomarker_var))])

  # Plot
  pheatmap(
    t(cluster_mat[, top_biomarkers]),
    scale = "row",
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    main = sprintf("Cluster Biomarker Profiles (Top %d)", length(top_biomarkers)),
    fontsize = 8
  )
}

#' Plot MCI risk scores
#'
#' @param mci_predictions Output from predict_mci_risk
#' @param group_by Variable to group by (default: "Ancestry")
#' @return ggplot object
plot_mci_risk_scores <- function(mci_predictions, group_by = "Ancestry") {

  p <- ggplot(mci_predictions, aes(x = get(group_by), y = AD_risk_score,
                                   fill = get(group_by))) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    labs(
      title = "AD Risk Scores for MCI Samples",
      x = group_by,
      y = "Predicted AD Risk Score",
      fill = group_by
    ) +
    theme_bw() +
    theme(legend.position = "none")

  return(p)
}

#' Create composite summary plot
#'
#' @param elastic_net_results Elastic net results
#' @param cluster_result Clustering results
#' @param mci_predictions MCI predictions
#' @return Combined plot
create_summary_plot <- function(elastic_net_results, cluster_result, mci_predictions) {

  p1 <- plot_roc_by_ancestry(elastic_net_results)
  p2 <- plot_mci_clusters(cluster_result)
  p3 <- plot_mci_risk_scores(mci_predictions)

  combined <- (p1 | p2) / p3

  combined + plot_annotation(
    title = "Biomarker Analysis Summary",
    subtitle = sprintf("N samples = %d", nrow(cluster_result$data))
  )
}
