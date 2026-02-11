# Handover Document: MCI Progression Analysis Pipeline

## Project Location

```
/Users/emartin1/Library/CloudStorage/Box-Box/AD/BiomarkerProject2024/
  NewAlamarDataAugust2025/Alamar_QC_Protocol/R_Code_Dev/
  Alamar_Biomarker_QC_Repo/In_Development/MCI_Progression_Analysis/
```

## Purpose

R/RMarkdown pipeline for analyzing MCI (Mild Cognitive Impairment) progression risk using Alamar biomarker platform data (~130 biomarkers, ~2100 samples). The pipeline classifies NCI vs AD, predicts AD risk for MCI patients, identifies features separating cognitively impaired from normal, and clusters MCI samples into subtypes with both unsupervised and supervised PCA approaches.

---

## Architecture Overview

**Main pipeline:** `MCI_Progression_Analysis.Rmd` — sections A through I, knits to HTML.

**Modular R scripts** (sourced by the Rmd):

| File | Role | Status |
|------|------|--------|
| `00_data_loading.R` | Load CSV, define biomarker/covariate columns | Working |
| `covariate_explorer.R` | Config-driven diagnosis + ancestry grouping utilities | Working |
| `01_preprocessing.R` | APOE4 status, missing data filtering, min sample filtering | Working |
| `02_composite_scores.R` | ATN, BD enrichment, inflammation PCs, proteinopathy burden | Working |
| `03_elastic_net.R` | Elastic net classification (NCI vs AD, MCI+AD vs NCI, MCI risk) | Working |
| `04_clustering.R` | K-means clustering, profiling, PCA/UMAP, supervised PCA | Working |
| `05_visualization.R` | All plotting functions | Working |

**Helper files** (not sourced by Rmd):

| File | Role |
|------|------|
| `QUICK_REFERENCE.R` | Copy-paste command reference for interactive use |
| `example_exploratory.R` | Tutorial for extended analyses |

**Data flow:** CSV → `load_alamar_data()` → apply dx/ancestry groupings → preprocess → composite scores → elastic net → MCI risk prediction → MCI+AD vs NCI → clustering → summary

---

## Pipeline Sections in the Rmd

Sections are ordered so that supervised analyses (elastic net) are grouped together (E → F → G), followed by unsupervised discovery (H), then summary (I). This reflects the logical dependency: elastic net models must exist before MCI risk prediction can use them.

| Section | What it does | Status |
|---------|-------------|--------|
| A. Overview | Title, setup, load libraries, source scripts | Working |
| B. Covariate Exploration | CDX/Group/Race/Ethnicity distributions, config-driven groupings | Working |
| C. Preprocessing | APOE4 status, missing biomarker filter, min sample filter | Working |
| D. Composite Scores | ATN, BD enrichment, inflammation PCs, proteinopathy, export CSV, AD-likeness | Working |
| E. Elastic Net (NCI vs AD) | Per-ancestry elastic net with covariates, ROC, top features | Working |
| F. MCI Risk Prediction | Apply NCI/AD model to MCI samples (pan-ancestry + per-ancestry) | Working |
| G. MCI+AD vs NCI | Supervised feature selection for cognitive impairment | Working |
| H. MCI Clustering | K-means with covariates, PCA, supervised PCA, heatmap | Working |
| I. Summary | Summary stats table | Working |

---

## What Has Been Completed

### Session 1 (original development)

1. **Table rendering** — all chunks use `results='asis'` + `kable()` for clean HTML tables.
2. **ATN score export** — Section D.8 exports `ATN_scores_with_covariates.csv`.
3. **Covariates in elastic net** — `03_elastic_net.R` accepts `include_covariates` parameter with dummy-coded categoricals and median-imputed numerics.
4. **More PCs + NCI/AD overlay** — `add_dim_reduction()` stores 5 PCs, projects NCI/AD into MCI PCA space.
5. **MCI+AD vs NCI feature selection** — `run_mci_ad_vs_nci()` function in `03_elastic_net.R`.

### Session 2 (bug fixes + new features)

6. **Fixed elastic net crash** — NA APOE.geno values producing NAs in dummy coding; whitespace normalization; `message()` → `warning()` in error handler. Applied in 3 functions in `03_elastic_net.R`.

7. **Fixed clustering bugs (5 fixes in `04_clustering.R`):**
   - NaN guard after `scale()` — removes zero-variance columns instead of crashing `kmeans()`
   - `gap_stat` storage — uses `auto_k` boolean flag instead of checking `is.character(n_clusters)` after reassignment
   - Scale mismatch in AD distance — stores `scale_center`/`scale_sd` and uses them to scale AD biomarkers before distance calculation
   - Covariate dimension mismatch — stores `scale_center_bio`/`scale_sd_bio` (biomarker-only scaling params) so `calculate_cluster_ad_distance()` works when covariates are included in clustering
   - `clusGap` convergence — added `iter.max = 100` to internal kmeans calls (default was 10)

8. **Section reordering** — moved clustering from F to H, grouped supervised analyses (E → F → G). Added methodology descriptions to each section.

9. **Composite score patterns (D.9)** — violin + jitter plots of all composite scores across NCI/MCI/AD/Dementia_Other.

10. **AD-likeness score (D.10)** — z-scores composite scores relative to NCI distribution (sign-flipped A_score so positive = more AD-like), creates per-sample AD-likeness index.

11. **Ancestry-specific MCI predictions (F.4)** — applies per-ancestry elastic net models to their matching MCI samples (in addition to the pan-ancestry "All" model in F.1-F.3).

12. **Covariates in clustering** — `cluster_mci_samples()` now accepts `include_covariates` parameter; dummy-codes categoricals, median-imputes numerics, equally weights all features after scaling. The Rmd passes `c("age_at_subject", "sex", "APOE4_carrier")`.

13. **PC-ancestry correlation (H.6)** — eta-squared (ANOVA) table showing how much variance in each PC is explained by ancestry.

14. **Supervised PCA (H.7)** — new `supervised_pca_clustering()` function: defines PCA axes from AD+NCI reference samples, projects MCI into disease-defined space. Includes scatter plot, AD-fraction table (position along NCI-to-AD axis), and ancestry eta-squared comparison with unsupervised PCA.

15. **kable() inside if-blocks fix** — wrapped `kable()` calls in `print()` for H.6 and H.7.3 (kable doesn't auto-print inside control flow in RMarkdown).

---

## Key Data Facts

- **Input file:** `input_files/filtered_combined_post_QC_mod.csv`
- **Samples after preprocessing:** ~2100 (exact count depends on missing data filtering)
- **Biomarkers:** ~130 (columns 26 onward in the CSV)
- **All data is log2-transformed (NPQ)** — differences = log2 fold-changes
- **Diagnosis groups (CDX_collapsed):** NCI, MCI, AD, Dementia_Other
- **Ancestry groups:** HI_WH, HI_MU, HI_BL, AA, AFDC
- **APOE.geno:** values like "3    3", "3    4" — 206 NAs — whitespace is normalized in elastic net and clustering but NOT in raw data
- **APOE4_carrier:** derived from APOE.geno in `01_preprocessing.R`

### Sample counts (NCI vs AD by ancestry, approximate):
| Ancestry | NCI | AD | Total |
|----------|-----|----|-------|
| All | 1026 | 302 | 1328 |
| AA | 455 | 65 | 520 |
| HI_MU | 249 | 87 | 336 |
| HI_WH | 141 | 55 | 196 |
| AFDC | 136 | 86 | 222 |
| HI_BL | 37 | 7 | 44 |

HI_BL has too few AD samples (7 < min_per_class=10), so elastic net correctly skips it.

---

## File-by-File Reference

### `00_data_loading.R`
- `load_alamar_data(filepath)` → returns `dat_list` with `$data`, `$biomarker_cols`, `$covariate_cols`, `$n_biomarkers`, `$n_samples`
- Biomarker columns defined as columns 26 onward (hardcoded index)

### `covariate_explorer.R`
- `apply_category_config(dat, config)` — applies keep/collapse/exclude for CDX
- `apply_grouping_config(dat, config)` — applies rule-based ancestry grouping
- `summarize_category_decisions()`, `summarize_grouping_decisions()` — for reporting

### `01_preprocessing.R`
- `create_apoe4_status()` — creates `APOE4_carrier` ("APOE4+"/"APOE4-"/NA) and `APOE4_count` (0/1/2/NA)
- `filter_missing_biomarkers(max_missing_pct = 0.2)` — removes samples with >20% missing biomarkers
- `filter_min_samples(min_n = 5)` — removes Ancestry x CDX groups with <5 samples

### `02_composite_scores.R`
- `calculate_atn_scores()` — A_score (A_42 - A_40), T_score (z-scored pTau mean), N_score (z-scored neurodegeneration mean)
- `calculate_bd_enrichment()` — BD minus peripheral pTau (log2 differences)
- `calculate_inflammation_pcs(n_pcs = 3)` — PCA on cytokine/chemokine markers
- `calculate_proteinopathy_burden()` — z-scored multi-protein index

### `03_elastic_net.R`
- `run_elastic_net_by_ancestry(dat_list, ancestry_group, alpha, nfolds, min_per_class, include_covariates)` — fits elastic net NCI vs AD for one ancestry
- `run_elastic_net_all_ancestries(dat_list, ancestries, alpha, min_per_class, include_covariates)` — runs above for all ancestries, skips failures
- `compare_biomarkers_across_ancestries(elastic_net_results)` — cross-ancestry comparison table
- `predict_mci_risk(dat_list, elastic_net_result, ancestry_group)` — applies NCI/AD model to MCI samples
- `run_mci_ad_vs_nci(dat_list, alpha, nfolds, include_covariates)` — elastic net for MCI+AD vs NCI

### `04_clustering.R`
- `cluster_mci_samples(dat_list, n_clusters, biomarker_subset, top_n, include_covariates)` — K-means on MCI samples. Selects top-N most variable biomarkers, optionally adds covariates (dummy-coded categoricals, median-imputed numerics). Returns cluster assignments, centers, silhouette, scaling params (`scale_center_bio`/`scale_sd_bio` for biomarker-only subset).
- `profile_mci_clusters(cluster_result, dat_list)` — returns list with biomarker/age/apoe/ancestry profiles
- `calculate_cluster_ad_distance(dat_list, cluster_result, ancestry_group)` — Euclidean distance of cluster centroids to AD centroid (uses biomarker-only scaling params)
- `add_dim_reduction(cluster_result, dat_list, method, n_pcs, project_all)` — PCA (or UMAP), stores PC1-5, optionally projects NCI/AD into MCI PCA space
- `supervised_pca_clustering(dat_list, cluster_result, n_pcs, n_clusters)` — defines PCA from AD+NCI reference, projects all groups, re-clusters MCI, computes per-PC AD-direction analysis

### `05_visualization.R`
- `plot_roc_by_ancestry()` — faceted ROC curves with AUC labels
- `plot_top_biomarkers(top_n)` — bar chart of elastic net coefficients, covariates shown faded
- `plot_mci_clusters(color_by)` — scatter plot in PC1/PC2 space
- `plot_all_diagnoses_pca()` — NCI/MCI/AD in MCI PCA space with ellipses
- `plot_pc_pair(pc_x, pc_y, color_by)` — arbitrary PC pairs
- `plot_cluster_heatmap(top_n)` — pheatmap of cluster biomarker profiles
- `plot_mci_risk_scores(group_by)` — boxplot of AD risk scores
- `create_summary_plot()` — composite patchwork plot

---

## Output Files (written to `output_files/`)

| File | Description |
|------|-------------|
| `preprocessed_data.rds` | dat_list after preprocessing |
| `data_with_composite_scores.rds` | dat_list with ATN + BD + inflammation + proteinopathy scores |
| `ATN_scores_with_covariates.csv` | Composite scores with all covariates for external use |
| `elastic_net_results.rds` | Per-ancestry elastic net results |
| `biomarkers_comparison_across_ancestries.csv` | Cross-ancestry feature selection comparison |
| `mci_clustering_results.rds` | Clustering results (before dim reduction) |
| `mci_clustering_with_dimred.rds` | Clustering results with PCA |
| `mci_cluster_profiles.csv` | Biomarker profiles by cluster |
| `cluster_distances_to_AD.csv` | Cluster centroid distances to AD centroid |
| `mci_ad_risk_predictions.csv` | MCI risk predictions from All model |
| `mci_ad_vs_nci_results.rds` | MCI+AD vs NCI elastic net results |
| `supervised_pca_results.rds` | Supervised PCA results (AD+NCI reference space) |
| `analysis_summary.txt` | Summary statistics |

---

## Design Decisions

- **Section ordering:** Supervised analyses (E, F, G) grouped before unsupervised (H) because elastic net models must exist before MCI risk prediction. Clustering is independent and comes last.
- **Clustering is independent of elastic net** — uses top 50 most variable biomarkers across MCI samples. This is intentional unsupervised discovery.
- **Covariates in clustering are equally weighted** — after scaling, age/sex/APOE4 are treated the same as biomarkers. This means a 1-SD change in age has the same clustering influence as a 1-SD change in any biomarker.
- **Supervised PCA complements unsupervised PCA** — unsupervised PCA (on MCI only) may produce axes dominated by ancestry or noise. Supervised PCA (trained on AD+NCI) defines axes with direct clinical meaning, so each PC reflects disease-related variation.
- **PCA is trained on MCI only (unsupervised), NCI/AD projected in** — preserves MCI-specific variance structure while showing where other diagnoses fall.
- **Elastic net uses alpha=0.5** (equal mix of LASSO and Ridge) — balances sparsity and grouped selection.
- **Covariate dummy coding drops first level** — standard approach to avoid multicollinearity.
- **NA APOE.geno → all dummy columns = 0** — means "not known to be any specific genotype," handled correctly by elastic net and clustering.
- **ATN A_score = A_42 - A_40 (not ratio)** — on log2 scale, difference equals log2(ratio), avoiding division issues.
- **T_score and N_score are z-scored before averaging** — because individual markers have very different NPQ ranges.

---

## Possible Future Work

1. **Covariates in PCA** — currently PCA (unsupervised) uses only biomarkers. Could optionally include covariates in PCA the same way clustering does.
2. **PCA-vs-kmeans comparison** — compare cluster assignments from unsupervised vs supervised PCA to see if disease-defined axes change subtype definitions.
3. **MCI+AD vs NCI stratified by ancestry** — currently Section G runs pan-ancestry only. Could add per-ancestry MCI+AD vs NCI models.
4. **Longitudinal validation** — when follow-up data arrives, validate MCI risk scores against actual progression, calibrate ancestry-specific thresholds.
