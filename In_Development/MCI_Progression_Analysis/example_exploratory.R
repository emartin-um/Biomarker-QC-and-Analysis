################################################################################
# example_exploratory.R
# Examples of using modules for exploratory analysis
################################################################################

# Load all modules
source("00_data_loading.R")
source("01_preprocessing.R")
source("02_composite_scores.R")
source("03_elastic_net.R")
source("04_clustering.R")
source("05_visualization.R")

library(tidyverse)

################################################################################
# EXAMPLE 1: Quick data exploration
################################################################################

# Load data
data_path <- "input_files/filtered_combined_post_QC.csv"
dat_list <- load_alamar_data(data_path)

# Check structure
glimpse(dat_list$data)

# What diagnoses do we have?
table(dat_list$data$CDX)

# What ancestries?
table(dat_list$data$Group, dat_list$data$Ethnicity)

################################################################################
# EXAMPLE 2: Custom preprocessing
################################################################################

# Load and apply custom filters
dat_list <- load_alamar_data(data_path) %>%
  collapse_diagnoses(keep_diagnoses = c("NCI", "MCI", "AD"), 
                     collapse_other = FALSE) %>%  # Drop "Other"
  define_ancestry() %>%
  create_apoe4_status() %>%
  filter_missing_biomarkers(max_missing_pct = 0.15) %>%  # Stricter
  filter_min_samples(min_n = 10)  # Require ≥10 per group

# Check final sample sizes
dat_list$data %>%
  count(Ancestry, CDX_collapsed)

################################################################################
# EXAMPLE 3: Test specific biomarker panel
################################################################################

# Define your own biomarker panel
tau_panel <- c("BD_pTau_181", "BD_pTau_217", "BD_pTau_231",
               "pTau_181", "pTau_217", "pTau_231", "BD_MAPT", "MAPT")

neurodegeneration_panel <- c("NFL", "NEFL", "NEFH", "GFAP", "NRGN", 
                             "SNAP25", "VGF", "NPTX1", "NPTX2")

inflammation_panel <- c("IL6", "IL1B", "TNF", "CRP", "TREM2", "GFAP",
                       "CCL2", "CXCL10")

# Cluster MCI using only tau panel
cluster_tau <- cluster_mci_samples(
  dat_list,
  n_clusters = 3,  # Force 3 clusters
  biomarker_subset = tau_panel
)

profile_mci_clusters(cluster_tau, dat_list)

################################################################################
# EXAMPLE 4: Compare model performance across ancestries
################################################################################

# Run elastic net for each ancestry
results_hispanic <- run_elastic_net_by_ancestry(dat_list, "Hispanic", alpha = 0.5)
results_aa <- run_elastic_net_by_ancestry(dat_list, "African_American", alpha = 0.5)
results_african <- run_elastic_net_by_ancestry(dat_list, "African", alpha = 0.5)

# Compare AUCs
data.frame(
  Ancestry = c("Hispanic", "African_American", "African"),
  AUC = c(results_hispanic$auc, results_aa$auc, results_african$auc),
  N_samples = c(results_hispanic$n_samples, results_aa$n_samples, results_african$n_samples),
  N_biomarkers = c(results_hispanic$n_biomarkers_selected, 
                   results_aa$n_biomarkers_selected,
                   results_african$n_biomarkers_selected)
)

################################################################################
# EXAMPLE 5: Examine specific biomarker patterns
################################################################################

# Add composite scores
dat_list <- create_all_composite_scores(dat_list)

# Plot ATN scores by diagnosis and ancestry
dat_list$data %>%
  filter(CDX_collapsed %in% c("NCI", "MCI", "AD")) %>%
  ggplot(aes(x = CDX_collapsed, y = A_score, fill = Ancestry)) +
  geom_boxplot() +
  labs(title = "Amyloid Score by Diagnosis and Ancestry",
       x = "Diagnosis", y = "A score (-log10 Aβ42/40)") +
  theme_bw()

# Compare BD-tau enrichment
dat_list$data %>%
  filter(CDX_collapsed %in% c("NCI", "MCI", "AD")) %>%
  filter(!is.na(BD_pTau217_ratio)) %>%
  ggplot(aes(x = Ancestry, y = BD_pTau217_ratio, fill = CDX_collapsed)) +
  geom_boxplot() +
  labs(title = "Brain-Derived pTau-217 Enrichment",
       y = "BD-pTau217 / Peripheral pTau-217 (log2)") +
  theme_bw()

################################################################################
# EXAMPLE 6: Identify high-risk MCI samples
################################################################################

# Train model on all ancestries
model_all <- run_elastic_net_by_ancestry(dat_list, "All", alpha = 0.5)

# Predict MCI risk
mci_risk <- predict_mci_risk(dat_list, model_all, "All")

# Identify high-risk MCI
high_risk_mci <- mci_risk %>%
  filter(risk_category == "High") %>%
  arrange(desc(AD_risk_score))

message(sprintf("\nFound %d high-risk MCI samples", nrow(high_risk_mci)))
print(head(high_risk_mci, 10))

# What biomarkers drive their risk?
top_biomarkers <- model_all$coefficients %>%
  slice_max(abs(coefficient), n = 10)

message("\nTop 10 biomarkers driving risk:")
print(top_biomarkers)

################################################################################
# EXAMPLE 7: Stratified analysis by APOE4
################################################################################

# Split by APOE4 status
dat_list$data %>%
  filter(CDX_collapsed == "MCI", !is.na(APOE4_carrier)) %>%
  count(APOE4_carrier, Ancestry)

# Do APOE4+ MCI have different biomarker profiles?
# Run clustering separately for APOE4+ and APOE4-
dat_apoe4_pos <- dat_list
dat_apoe4_pos$data <- dat_list$data %>%
  filter(APOE4_carrier == "APOE4+")

dat_apoe4_neg <- dat_list
dat_apoe4_neg$data <- dat_list$data %>%
  filter(APOE4_carrier == "APOE4-")

if (nrow(dat_apoe4_pos$data %>% filter(CDX_collapsed == "MCI")) >= 20) {
  cluster_apoe4_pos <- cluster_mci_samples(dat_apoe4_pos, n_clusters = 3, top_n = 30)
  message("\nAPOE4+ MCI clusters:")
  print(table(cluster_apoe4_pos$cluster_assignments))
}

if (nrow(dat_apoe4_neg$data %>% filter(CDX_collapsed == "MCI")) >= 20) {
  cluster_apoe4_neg <- cluster_mci_samples(dat_apoe4_neg, n_clusters = 3, top_n = 30)
  message("\nAPOE4- MCI clusters:")
  print(table(cluster_apoe4_neg$cluster_assignments))
}

################################################################################
# EXAMPLE 8: Create custom composite score
################################################################################

# Create a vascular dysfunction score
vascular_markers <- c("VEGFA", "VEGFD", "FLT1", "KDR", "TEK", 
                     "VCAM1", "ICAM1", "CD106")

dat_list$data <- dat_list$data %>%
  mutate(
    Vascular_score = rowMeans(select(., any_of(vascular_markers)), na.rm = TRUE)
  )

# Test if vascular score differs by diagnosis
dat_list$data %>%
  filter(CDX_collapsed %in% c("NCI", "MCI", "AD")) %>%
  ggplot(aes(x = CDX_collapsed, y = Vascular_score, fill = CDX_collapsed)) +
  geom_boxplot() +
  facet_wrap(~Ancestry) +
  labs(title = "Vascular Dysfunction Score",
       x = "Diagnosis", y = "Mean Vascular Biomarker Level") +
  theme_bw()

message("\n=== Exploratory analysis examples complete ===")
