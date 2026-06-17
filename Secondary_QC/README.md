# Secondary QC

Secondary quality control analyses for Alamar biomarker data. These analyses are optional but recommended for comprehensive QC assessment.

## Contents

### Replicate_Analysis/
Analysis of technical replicate samples (HIHG, IPC, SC, NC) for assessing assay reproducibility and concordance. Uses raw NPQ data directly from Primary_QC.

**Scripts:**
- `QC_Pipeline_Replicates.Rmd` - CV analysis and replicate concordance

**Key outputs:**
- CV statistics by biomarker at plate, Run, Bay, and overall levels
- Concordance reports
- Replicate pair comparisons
- mCherry (IC) diagnostics

### Hemolysis_Check/
Analysis of hemolysis markers (HBA1, PGK1, MDH1, SOD1, ENO2) to identify pre-analytical sample handling issues vs biological hemolysis effects.

**Scripts:**
- `Hemolysis_by_Site.Rmd` - Comprehensive hemolysis analysis

**Key outputs:**
- Hemolysis index by site
- Sex differences (G6PD-related effects)
- Covariate associations
- Outlier detection flags

### APOE/
Analysis of APOE biomarker levels (APOE4 − APOE derived metric) stratified by APOE genotype to validate assay performance and identify samples with discordant biomarker/genotype results.

**Scripts:**
- `APOE_geno_protein.qmd` - APOE genotype vs biomarker analysis

**Key outputs:**
- Boxplots of APOE4−APOE by genotype
- Extreme outliers by genotype (3×IQR threshold)
- `APOE4_minus_APOE_extreme_outliers_by_genotype.csv`

### Batch_Effects/ ⚠️ NEW — NEEDS TESTING
Assesses Run- and Bay-level technical batch effects in the merged dataset. Quantifies
how much variance in each protein is explained by batch vs. biological covariates,
and makes a data-driven recommendation for whether batch correction is needed.

Requires `Metadata_Merge` output (`merged_combined_post_QC.csv`).

**Scripts:**
- `Batch_Effects_Analysis.Rmd` - Full batch effect assessment

**Key outputs:**
- `batch_covariate_confounding.csv` — Cramér's V between Run/Bay and Race/sex/age/Site/Dx
- `variance_partition_results.csv` — Per-protein % variance attributed to Run, Bay, biology
- `icc_results.csv` — Per-protein ICC within-Run (across Bays) and cross-Run
- `batch_effect_test_results.csv` — Per-protein FDR-corrected F-test for Run and Bay terms
- Summary section with automated correction recommendation

**Required packages (beyond base tidyverse):**
```r
BiocManager::install(c("variancePartition", "BiocParallel"))
install.packages(c("lme4", "irr", "patchwork", "car"))
```

> **TODO:** Script was written March 2026 and has not yet been knitted on real data.
> Run on current dataset, check output, and add notes to `notes.md`.

---

## Workflow

```
Primary_QC
├── Replicate_Analysis      (raw NPQ; run directly after Primary_QC)
└── Metadata_Merge
    ├── Hemolysis_Check     (requires merged data)
    ├── APOE                (requires merged data with genotype info)
    └── Batch_Effects       (requires merged_combined_post_QC.csv)  ← NEW
```

## Getting Started

1. Open the relevant `.Rproj` file in RStudio
2. Input files are auto-copied from upstream outputs if not present
3. Run the Rmd scripts; outputs will be saved to `output_files/`
