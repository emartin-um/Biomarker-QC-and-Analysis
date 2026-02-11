################################################################################
# QUICK REFERENCE GUIDE
################################################################################

# This file contains copy-paste code snippets for common analysis tasks

################################################################################
# SETUP
################################################################################

# Load all modules
source("00_data_loading.R")
source("01_preprocessing.R")
source("02_composite_scores.R")
source("03_elastic_net.R")
source("04_clustering.R")
source("05_visualization.R")
library(tidyverse)

# Standard data loading
dat_list <- preprocess_standard("input_files/filtered_combined_post_QC.csv")

################################################################################
# COMMON TASK 1: Find best biomarkers for NCI vs AD
################################################################################

# Run for all ancestries
results <- run_elastic_net_all_ancestries(dat_list, alpha = 0.5)

# Get universal biomarkers (selected in all groups)
universal <- compare_biomarkers_across_ancestries(results) %>%
  filter(n_ancestries_selected >= 3)
print(universal)

################################################################################
# COMMON TASK 2: Identify MCI subtypes
################################################################################

# Cluster MCI (auto k)
clusters <- cluster_mci_samples(dat_list, n_clusters = "auto", top_n = 50)

# Profile clusters
profiles <- profile_mci_clusters(clusters, dat_list)

# Which cluster is closest to AD?
distances <- calculate_cluster_ad_distance(dat_list, clusters, "All")
print(distances)  # Lower distance = more AD-like

################################################################################
# COMMON TASK 3: Score MCI samples for AD risk
################################################################################

# Train on all ancestries
model <- run_elastic_net_by_ancestry(dat_list, "All", alpha = 0.5)

# Score MCI
mci_scores <- predict_mci_risk(dat_list, model, "All")

# High-risk MCI
high_risk <- mci_scores %>% 
  filter(risk_category == "High") %>%
  arrange(desc(AD_risk_score))

write_csv(high_risk, "high_risk_mci.csv")

################################################################################
# COMMON TASK 4: Test specific biomarker hypothesis
################################################################################

# Example: Do BD-tau/peripheral tau ratios predict diagnosis?
dat_list <- calculate_bd_enrichment(dat_list)

dat_list$data %>%
  filter(CDX_collapsed %in% c("NCI", "AD")) %>%
  filter(!is.na(BD_pTau217_ratio)) %>%
  ggplot(aes(x = CDX_collapsed, y = BD_pTau217_ratio)) +
  geom_boxplot() +
  facet_wrap(~Ancestry) +
  theme_bw()

# Statistical test
wilcox.test(
  BD_pTau217_ratio ~ CDX_collapsed,
  data = dat_list$data %>% filter(CDX_collapsed %in% c("NCI", "AD"))
)

################################################################################
# COMMON TASK 5: Compare ancestry-specific vs universal model
################################################################################

# Train ancestry-specific models
model_hisp <- run_elastic_net_by_ancestry(dat_list, "Hispanic", alpha = 0.5)
model_aa <- run_elastic_net_by_ancestry(dat_list, "African_American", alpha = 0.5)
model_all <- run_elastic_net_by_ancestry(dat_list, "All", alpha = 0.5)

# Compare performance
data.frame(
  Model = c("Hispanic-specific", "AA-specific", "Universal"),
  AUC = c(model_hisp$auc, model_aa$auc, model_all$auc),
  N_biomarkers = c(model_hisp$n_biomarkers_selected,
                   model_aa$n_biomarkers_selected,
                   model_all$n_biomarkers_selected)
)

################################################################################
# COMMON TASK 6: Create ATN profiles
################################################################################

dat_list <- calculate_atn_scores(dat_list)

# Visualize ATN by diagnosis
dat_list$data %>%
  filter(CDX_collapsed %in% c("NCI", "MCI", "AD")) %>%
  select(CDX_collapsed, A_score, T_score, N_score) %>%
  pivot_longer(cols = c(A_score, T_score, N_score), 
               names_to = "Score", values_to = "Value") %>%
  ggplot(aes(x = CDX_collapsed, y = Value, fill = Score)) +
  geom_boxplot() +
  facet_wrap(~Score, scales = "free_y") +
  theme_bw()

################################################################################
# COMMON TASK 7: Examine biomarker correlations within groups
################################################################################

# Correlation matrix for NCI
nci_data <- dat_list$data %>%
  filter(CDX_collapsed == "NCI") %>%
  select(all_of(dat_list$biomarker_cols))

cor_matrix <- cor(nci_data, use = "pairwise.complete.obs")

# Plot heatmap
library(pheatmap)
pheatmap(cor_matrix, 
         show_rownames = FALSE, 
         show_colnames = FALSE,
         main = "Biomarker Correlations in NCI")

################################################################################
# COMMON TASK 8: Export results for external validation
################################################################################

# Export MCI predictions with biomarker data
mci_export <- dat_list$data %>%
  filter(CDX_collapsed == "MCI") %>%
  left_join(mci_scores, by = "Record_ID") %>%
  select(Record_ID, Ancestry, age_at_subject, APOE4_carrier,
         AD_risk_score, risk_category,
         BD_pTau_217, NFL, GFAP, A_42, A_40)

write_csv(mci_export, "mci_predictions_with_biomarkers.csv")

################################################################################
# COMMON TASK 9: Sample size calculations for stratified analysis
################################################################################

# Check if we have enough samples for stratified analysis
dat_list$data %>%
  count(Ancestry, CDX_collapsed, APOE4_carrier) %>%
  pivot_wider(names_from = CDX_collapsed, values_from = n, values_fill = 0)

################################################################################
# COMMON TASK 10: Quick quality check
################################################################################

# Missing data by diagnosis
missing_summary <- dat_list$data %>%
  select(CDX_collapsed, all_of(dat_list$biomarker_cols)) %>%
  group_by(CDX_collapsed) %>%
  summarize(
    n_samples = n(),
    pct_missing = mean(is.na(across(all_of(dat_list$biomarker_cols)))) * 100
  )

print(missing_summary)

# Outlier detection (samples with extreme biomarker values)
biomarker_mat <- get_biomarker_matrix(dat_list)
biomarker_mat_scaled <- scale(biomarker_mat)

extreme_values <- rowSums(abs(biomarker_mat_scaled) > 5, na.rm = TRUE)
outlier_samples <- which(extreme_values > 10)  # >10 biomarkers with z>5

if (length(outlier_samples) > 0) {
  message(sprintf("Warning: %d samples with >10 extreme biomarker values", 
                  length(outlier_samples)))
  print(dat_list$data[outlier_samples, c("Record_ID", "CDX_collapsed", "Ancestry")])
}

################################################################################
# END QUICK REFERENCE
################################################################################
