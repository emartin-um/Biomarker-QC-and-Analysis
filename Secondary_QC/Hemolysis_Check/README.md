# Hemolysis Check

Analysis of hemolysis markers to distinguish pre-analytical sample handling artifacts from biological hemolysis effects.

## Overview

Hemolysis (rupture of red blood cells) can artificially elevate certain biomarkers in plasma samples. This analysis examines five intracellular proteins that indicate RBC lysis:

- **HBA1** (Hemoglobin Subunit Alpha 1): Gold-standard hemolysis indicator
- **PGK1** (Phosphoglycerate Kinase 1): Glycolytic enzyme from RBCs
- **MDH1** (Malate Dehydrogenase 1): Citric acid cycle enzyme
- **SOD1** (Superoxide Dismutase 1): Antioxidant enzyme abundant in RBCs
- **ENO2** (Enolase 2): Present in RBCs, used as hemolysis marker in clinical labs

## Scripts

### Hemolysis_by_Site.Rmd

Comprehensive hemolysis analysis including:

1. **Hemolysis Index Creation** - Mean of available hemolysis markers
2. **Site Effects** - ANOVA testing for site differences
3. **Individual Markers** - Per-marker analysis by site
4. **Correlation Analysis** - Within-site marker correlations
5. **Biomarker Confounding** - Test if site differences are hemolysis artifacts
6. **Site × Hemolysis Interactions** - Site-specific hemolysis effects
7. **Sex Differences** - G6PD deficiency effects (X-linked, more common in males)
8. **Covariate Associations** - RUN, Race, age, CDX, AOO, BMI correlations
9. **Outlier Detection** - Overall and site-specific outlier flagging
10. **Technical Replicate Variance** - Compare sample variance to technical baseline

## Data Sources

### Input
- `input_files/filtered_combined_post_QC.csv` (primary)
- Falls back to `../Metadata_Merge/output_files/filtered/filtered_combined_post_QC.csv`

> **Heads up (Metadata_Merge Step 8):** `filtered_combined_post_QC.csv` is now **cross-sectional —
> one row per subject** (most recent visit). Hemolysis QC is *per-sample*; if you want to screen
> **every assayed sample** (including earlier visits of repeat subjects), read
> `merged_combined_post_QC.csv` or `filtered_combined_longitudinal.csv` instead.

### Technical Replicates
- `../Primary_QC/output_files/hihg_reps_NPQ.csv` - HIHG replicate samples for baseline variance

## Output Files

- `output_files/hemo_data_output.csv` - Hemolysis index with metadata
- `output_files/hemo_data_with_outliers.csv` - Data with outlier flags
- `output_files/biomarker_hemolysis_assessment.csv` - Per-biomarker hemolysis assessment

## Interpretation Guide

### Pre-analytical vs Biological Hemolysis

**Pre-analytical (sample handling) indicators:**
- Site differences in hemolysis with no sex effect
- High inter-marker correlations (coordinated release)
- Hemolysis correlates with RUN but not biological covariates

**Biological hemolysis indicators:**
- Sex differences (males > females, consistent with G6PD deficiency)
- Consistent covariate associations across sites
- Sample variance exceeds technical baseline

### Outlier Thresholds

- **Mean + 2 SD**: ~5% expected flagged (sensitive)
- **Mean + 3 SD**: ~0.3% expected flagged (specific)
- **Site-specific**: Accounts for site-level differences

## Usage

1. Open `Hemolysis_Check.Rproj` in RStudio
2. Ensure data exists in `input_files/` or that Metadata_Merge has been run
3. Knit `Hemolysis_by_Site.Rmd`
4. Review HTML output and CSV files

## Notes

The `notes.md` file contains analysis-specific notes and is included in the HTML output but excluded from git (see `.gitignore`).

## Required Packages

```r
library(tidyverse)
library(broom)
library(car)
library(corrplot)
library(knitr)
library(kableExtra)
library(pheatmap)
```
