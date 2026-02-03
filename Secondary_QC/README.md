# Secondary QC

Secondary quality control analyses for Alamar biomarker data. These analyses are optional but recommended for comprehensive QC assessment.

## Contents

### Replicate_Analysis/
Analysis of technical replicate samples for assessing assay reproducibility and concordance.

**Scripts:**
- `QC_Pipeline_Replicates.Rmd` - CV analysis and replicate concordance

**Key outputs:**
- CV statistics by biomarker
- Concordance reports
- Replicate pair comparisons

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
- `APOE_exam.qmd` - APOE genotype vs biomarker analysis

**Key outputs:**
- Boxplots of APOE4−APOE by genotype
- Extreme outliers by genotype (3×IQR threshold)
- `APOE4_minus_APOE_extreme_outliers_by_genotype.csv`

## Workflow

These analyses can be run after Primary_QC:

```
Primary_QC
├── Replicate_Analysis (can run directly after Primary_QC)
└── Metadata_Merge
    ├── Hemolysis_Check (requires merged data)
    └── APOE (requires merged data with genotype info)
```

## Getting Started

1. Open the relevant `.Rproj` file in RStudio
2. Input files are auto-copied from upstream outputs if not present
3. Run the Rmd scripts; outputs will be saved to `output_files/`
