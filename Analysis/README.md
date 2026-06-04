# Analysis

Statistical analyses for biomarker association testing with clinical diagnosis (CDX), covariate effects, and interaction analyses.

> **Note:** These scripts now use the unified data pipeline from `Metadata_Merge/` and shared functions from `Utilities/shared_functions.R`.

## Usage

1. Open `Analysis.Rproj` in RStudio
2. **Preferred:** Run `Metadata_Merge/` pipeline first, then copy output to `input_files/`
3. Alternatively, place raw input files in `input_files/` (scripts have fallback loading)
4. Run analysis scripts in order (univariate first, then interaction)
5. Outputs will be saved to `output_files/`

## Scripts

### Assoc_Analysis_univariate_HighBiomarkers_12_19_25.Rmd

Univariate association analysis for high-quality biomarkers.

**Purpose:** Tests associations between biomarkers and covariates:
- Biomarker clustering by correlation
- Covariate effects (sex, age, race/ethnicity)
- Plate/batch effects (Run, Bay)
- Group heatmaps and distribution plots

**Key Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `FDR` | 0.05 | False discovery rate threshold |
| `n_groups` | 8 | Number of biomarker clusters |

**Input Files (in order of preference):**
1. `input_files/filtered_combined_post_QC.csv` - Pre-merged data from Metadata_Merge
2. `../Metadata_Merge/output_files/filtered/filtered_combined_post_QC.csv` - Direct fallback
3. Manual merge from `NPQ_[date]_post_QC.csv` + `U19_Alamar_metadata_[date].csv` (legacy)

> As of Metadata_Merge **Step 8**, `filtered_combined_post_QC.csv` is **cross-sectional — one row
> per subject** (most recent visit), which is correct for association analyses. For longitudinal
> work use `../Metadata_Merge/output_files/filtered/filtered_combined_longitudinal.csv`.

**Output Files:**
- `biomarker_groups.csv` - Cluster assignments for biomarkers
- `meta_plus_race_dx.csv` - Merged biomarker + metadata for downstream analysis
- Correlation heatmaps and association plots

**Shared Functions Used:**
- `test_covariate_associations()` - Association testing with FDR correction
- `plot_group_heatmaps()` - Correlation heatmaps by biomarker group
- `summarize_df_pretty()`, `summarize_counts()`, `protect_name()`

---

### Assoc_Analysis_Interaction_HighBiomarkers_12_31_25_1.Rmd

Interaction analysis for CDX associations with covariate adjustment.

**Purpose:** Tests CDX (clinical diagnosis) associations with:
- Main effects adjusted for covariates (sex, age, race/ethnicity)
- CDX x covariate interactions (CDX:age, CDX:sex, CDX:race)
- Outlier detection and handling
- Estimated marginal means with confidence intervals

**Key Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `FDR` | 0.05 | False discovery rate threshold |
| `min_group_n` | 3 | Minimum samples per CDX group |
| `remove_outliers` | FALSE | Whether to remove outliers |
| `outlier_method` | "iqr" | Outlier detection method |
| `outlier_threshold` | 4 | IQR multiplier for outliers |

**Input Files (in order of preference):**
1. `input_files/meta_plus_race_dx.csv` - From univariate analysis
2. `input_files/filtered_combined_post_QC.csv` - Pre-merged from Metadata_Merge
3. `../Metadata_Merge/output_files/filtered/filtered_combined_post_QC.csv` - Direct fallback
4. `biomarker_groups.csv` - Cluster assignments (optional, can skip clustering)

**Output Files:**
- `outlier_details.csv` - Identified outlier samples
- Association summary tables
- Interaction plots (CDX x age, CDX x sex, CDX x race)

## Input Files

Place your input data in `input_files/` (not tracked in git).

**Preferred:** Use pre-merged data from Metadata_Merge pipeline:
- `filtered_combined_post_QC.csv` - Combined biomarker + metadata

**Alternative:** Raw input files for manual merge:
- Post-QC NPQ data from `Primary_QC/output_files/`
- U19 sample metadata with phenotype information
- Biomarker group assignments (optional)

## Output Files

Generated outputs are saved to `output_files/` (not tracked in git):
- Association test results
- Interaction analysis summaries
- Publication-ready plots

## Dependencies

See main repository README for full package list. Key packages:
- `tidyverse` - Data manipulation
- `emmeans` - Estimated marginal means
- `car` - Type II ANOVA
- `broom` - Model tidying
- `pheatmap`, `caret` - Clustering
- `kableExtra` - Table formatting
- `gridExtra`, `patchwork` - Plot composition

## Workflow

Recommended analysis order:

```
Metadata_Merge/
    |
    +-- Run File_Merge.qmd (merge NPQ + metadata)
    +-- Run Apply_Filters.qmd (apply QC filters)
    |
    v
Metadata_Merge/output_files/filtered/
    |
    +-- filtered_combined_post_QC.csv
    |
    v
Analysis/input_files/ (copy or symlink)
    |
    v
[1. Run Assoc_Analysis_univariate.Rmd]
    |
    +-- biomarker_groups.csv
    +-- meta_plus_race_dx.csv
    |
    v
[2. Run Assoc_Analysis_Interaction.Rmd]
    |
    v
Analysis/output_files/
```

**Alternative (Legacy) Workflow:**
```
Primary_QC/output_files/ + External metadata
    |
    v
Analysis/input_files/
    |
    +-- NPQ_post_QC.csv
    +-- U19_Alamar_metadata.csv
    |
    v
[Scripts merge data internally]
```

## Analysis Notes

### CDX Groups
- **NCI**: No Cognitive Impairment (reference group)
- **MCI**: Mild Cognitive Impairment
- **AD**: Alzheimer's Disease

### Race/Ethnicity Groups
- **AA**: African American
- **AFDC**: African Descent Countries
- **WH_HI**: White Hispanic
- **MU_HI**: Mixed/Unknown Hispanic

### Key Findings to Check
- CCL3 shows polymorphic distribution across populations
- Some biomarkers show age x CDX interactions
- Sex differences may be present in some markers (e.g., GFAP, FCN2)
