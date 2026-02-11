# MCI Progression Analysis Pipeline

## Overview

R/RMarkdown pipeline for analyzing MCI (Mild Cognitive Impairment) progression risk using Alamar biomarker platform data (~130 biomarkers, ~2100 samples across Hispanic, African American, and African populations).

**Goal:** Identify MCI samples likely to progress to AD using cross-sectional biomarker data, discover MCI subtypes, and compare biomarker signatures across ancestry groups.

## Quick Start

**Primary entry point:** Open and knit `MCI_Progression_Analysis.Rmd` in RStudio.

```r
# Knit the full pipeline to HTML
rmarkdown::render("MCI_Progression_Analysis.Rmd")
```

The Rmd sources all modular R scripts automatically and produces a self-contained HTML report with all analyses, tables, and figures.

## Pipeline Structure

```
MCI_Progression_Analysis/
├── MCI_Progression_Analysis.Rmd   # Main pipeline (knit to HTML)
├── 00_data_loading.R              # Load data, define biomarker columns
├── 01_preprocessing.R             # APOE4 status, filtering, cleaning
├── 02_composite_scores.R          # ATN, BD-enrichment, inflammation PCs
├── 03_elastic_net.R               # NCI vs AD, MCI+AD vs NCI, MCI risk
├── 04_clustering.R                # K-means, supervised PCA, profiling
├── 05_visualization.R             # Plotting functions
├── covariate_explorer.R           # Config-driven diagnosis/ancestry grouping
├── QUICK_REFERENCE.R              # Copy-paste commands for interactive use
├── example_exploratory.R          # Tutorial for extended analyses
├── input_files/
│   └── filtered_combined_post_QC_mod.csv
└── output_files/                  # Created automatically
```

## Rmd Sections

| Section | Description |
|---------|-------------|
| A. Overview | Setup, load libraries, source scripts |
| B. Covariate Exploration | CDX/ancestry distributions, config-driven groupings |
| C. Preprocessing | APOE4 status, missing data filter, min sample filter |
| D. Composite Scores | ATN, BD enrichment, inflammation PCs, proteinopathy, AD-likeness score |
| E. Elastic Net (NCI vs AD) | Per-ancestry elastic net with covariates, ROC, top features |
| F. MCI Risk Prediction | Apply NCI/AD model to MCI (pan-ancestry + per-ancestry) |
| G. MCI+AD vs NCI | Supervised feature selection for cognitive impairment |
| H. MCI Clustering | K-means with covariates, PCA, supervised PCA, heatmap |
| I. Summary | Summary statistics |

Supervised analyses (E, F, G) are grouped before unsupervised discovery (H).

## Key Functions

### Composite Scores (`02_composite_scores.R`)
- `calculate_atn_scores()` — Amyloid/Tau/Neurodegeneration framework
- `calculate_bd_enrichment()` — Brain-derived vs peripheral ratios
- `calculate_inflammation_pcs(n_pcs)` — Cytokine principal components
- `calculate_proteinopathy_burden()` — Multi-pathology score

### Elastic Net (`03_elastic_net.R`)
- `run_elastic_net_all_ancestries(dat_list, ancestries, alpha, min_per_class, include_covariates)` — NCI vs AD for all ancestries
- `predict_mci_risk(dat_list, elastic_net_result, ancestry_group)` — Apply model to MCI samples
- `run_mci_ad_vs_nci(dat_list, alpha, nfolds, include_covariates)` — MCI+AD vs NCI classification
- `compare_biomarkers_across_ancestries(results)` — Universal vs ancestry-specific markers

### Clustering (`04_clustering.R`)
- `cluster_mci_samples(dat_list, n_clusters, biomarker_subset, top_n, include_covariates)` — K-means with optional covariates (dummy-coded categoricals, median-imputed numerics)
- `profile_mci_clusters(cluster_result, dat_list)` — Characterize clusters by biomarkers and covariates
- `calculate_cluster_ad_distance(dat_list, cluster_result, ancestry_group)` — Euclidean distance to AD centroid
- `add_dim_reduction(cluster_result, dat_list, method, n_pcs, project_all)` — PCA/UMAP with NCI/AD projection
- `supervised_pca_clustering(dat_list, cluster_result, n_pcs, n_clusters)` — PCA axes from AD+NCI reference, project MCI into disease-defined space

### Visualization (`05_visualization.R`)
- `plot_roc_by_ancestry()` — Faceted ROC curves with AUC
- `plot_top_biomarkers(top_n)` — Elastic net coefficients (covariates shown faded)
- `plot_mci_clusters(color_by)` — PCA scatter plot
- `plot_all_diagnoses_pca()` — NCI/MCI/AD overlay with ellipses
- `plot_cluster_heatmap(top_n)` — Biomarker profiles by cluster

## Output Files (`output_files/`)

| File | Description |
|------|-------------|
| `preprocessed_data.rds` | Cleaned data after preprocessing |
| `data_with_composite_scores.rds` | Data with ATN + BD + inflammation + proteinopathy scores |
| `ATN_scores_with_covariates.csv` | Composite scores with covariates (for external use) |
| `elastic_net_results.rds` | Per-ancestry elastic net models |
| `biomarkers_comparison_across_ancestries.csv` | Cross-ancestry feature comparison |
| `mci_clustering_results.rds` | Clustering results (before dim reduction) |
| `mci_clustering_with_dimred.rds` | Clustering results with PCA |
| `mci_cluster_profiles.csv` | Biomarker profiles by cluster |
| `cluster_distances_to_AD.csv` | Cluster centroid distances to AD |
| `mci_ad_risk_predictions.csv` | MCI risk predictions |
| `mci_ad_vs_nci_results.rds` | MCI+AD vs NCI elastic net results |
| `supervised_pca_results.rds` | Supervised PCA (AD+NCI reference space) |
| `analysis_summary.txt` | Summary statistics |

## Data Requirements

**Input file:** `input_files/filtered_combined_post_QC_mod.csv`

**Required columns:**
- Biomarkers: Columns 26+ (~130 biomarkers, log2-transformed NPQ values)
- Diagnosis: `CDX` (NCI, MCI, AD, Other)
- Demographics: `Group`, `Race`, `Ethnicity`, `age_at_subject`, `sex`
- Genetics: `APOE.geno`

## Dependencies

```r
install.packages(c("tidyverse", "glmnet", "pROC", "cluster", "factoextra",
                   "ggplot2", "patchwork", "pheatmap", "knitr", "rmarkdown"))
```

## Notes

- All biomarker data is log2-transformed (NPQ) — differences = log2 fold-changes
- Missing data handled via median imputation
- Random seed set to 123 for reproducibility
- Elastic net uses alpha=0.5 (elastic net mix of LASSO + Ridge)
- See `HANDOVER.md` for detailed technical documentation and design decisions
