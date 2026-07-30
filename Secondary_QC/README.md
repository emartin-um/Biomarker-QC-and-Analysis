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
Per-batch APOE QC: validates the APOE assay against genotype and surfaces discordances.

**Canonical script:** `APOE_perbatch_QC.Rmd` (params `merged_csv`, `output_dir`). Three sections:
- **A — APOE4−APOE protein-vs-genotype outliers.** Per run-bay × genotype 3×IQR banding with an **auto
  high-spread override** (2026-06-20): a plate whose within-genotype IQR ≫ the typical plate's is re-scored
  against a **clean reference band** (from normal-spread plates) so a scrambled plate's off-band samples are
  flagged, not hidden; normal-plate band widths are floored at the clean reference so tight clusters don't
  over-flag. → `APOE4_minus_APOE_extreme_outliers.xlsx`, `Outliers_Ranked.xlsx`, `high_spread_plates.xlsx`.
- **B — Sanger-vs-WGS concordance** (WGS-default; disagreements reported). → `APOE_Sanger_vs_WGS_discordant.xlsx`.
- **C — per-run APOE4 assay screen** (off-band % per run-bay; flags Bay1-type failures). → `APOE4_per_run_QC.xlsx`.

`APOE_geno_protein.qmd` is the **archived** May-2026 investigation record (`_archive_APOE_pre_perbatch_2026Jun16/`),
not the current pipeline. See `README.Rmd` for full details.

### Batch_Effects/ — variance partition · ICC · submission-year batch structure
Assesses Run- and Bay-level technical batch effects in the merged dataset. Quantifies
how much variance in each protein is explained by batch vs. biological covariates,
and makes a data-driven recommendation for whether batch correction is needed.

Requires `Metadata_Merge` output (`merged_combined_post_QC.csv`).

**Scripts:**
- `Batch_Effects_Analysis.Rmd` - Full batch effect assessment

**Key outputs:**
- `batch_covariate_confounding.csv` — Cramér's V between Run/Bay/**year** and the **derived** covariates (Ancestry, CDX_collapsed, sex, age, Site)
- `variance_partition_results.csv` — Per-protein % variance: Run, Bay, Ancestry, CDX_collapsed, sex, Site, age
- `icc_results.csv` — Per-protein ICC within-Run (across Bays) and cross-Run (from replicate samples)
- `batch_effect_test_results.csv` — Per-protein FDR-corrected F-test for Run and Bay terms
- `year_effect_per_protein.csv` + `hihg_pca_by_year.png` — **submission-year (2025 vs 2026) structure**: HIHG-control PCA/clustering + per-protein year shift
- Summary section ending with a concrete **what-to-adjust-for** recommendation (year, Site)

**Required packages (beyond base tidyverse):**
```r
BiocManager::install(c("variancePartition", "BiocParallel"))
install.packages(c("lme4", "irr", "patchwork", "car"))
```

> **Status (2026-06-23):** tested on the 50-plate combined run. Reads `merged_combined` (keeps the
> replicate samples the ICC needs) and **joins the derived covariates** (Ancestry, CDX_collapsed) from
> `filtered_combined`; `variancePartition` models categoricals as random effects; a submission-year
> (2025 vs 2026) clustering section + adjust-for recommendation were added. On this dataset the dominant
> batch axis is the **submission year**, confounded with **Site** (Cramér's V 0.55).

### Triage_Metadata/ — do the QC decisions know who the specimen came from?
Primary QC triages and flags on assay statistics only, and is agnostic to phenotype
by design. This module asks whether those decisions land evenly across **clinical
collection site**, diagnosis, ancestry, sex and age. Base R only; **diagnostic only**.

Requires `Metadata_Merge` output plus `Primary_QC/output_files/well_level_QC.csv`
(older runs are reconstructed automatically).

**Scripts:** `run_triage_metadata.R`, `test_triage_metadata.R`

**Key outputs:**
- `mutual_adjustment.csv` — each label given the others **and read depth**: which one is actually carrying an association, and which was a proxy
- `metadata_association_tests.csv` — per axis × variable: two permutation statistics, within-plate versions, conditioned deviances, and a **minimum detectable effect on every null**
- `triage_rates_by_metadata.csv` — per-group rates with Wilson CIs, groups below the size floor kept and marked
- `site_detail.csv`, `design_confounding.csv`, `balance_*.csv`, `unlabelled_wells.csv`

> **Status (2026-07-30):** on the 50-plate run the read flag's apparent site and
> diagnosis associations are **entirely read depth** (site keeps 13% of its
> deviance after adjustment), while the **IC axis** carries mutually independent
> associations with site (p = 0.046), ancestry (p = 0.016) and age (p = 0.030).
> It also **corrects** `Triage_Review` §5: the 40%-triage site is 2 wells of 5.

> ⚠️ **The other modules below are documented; `Extremes/`, `Specimen_Quality/`,
> `Detectability_LOD/` and `Triage_Review/` are not yet listed here.** See each
> module's own README.

---

## Workflow

```
Primary_QC
├── Replicate_Analysis      (raw NPQ; run directly after Primary_QC)
├── Triage_Review           (audits the Primary QC screens themselves)
└── Metadata_Merge
    ├── Hemolysis_Check     (requires merged data)
    ├── APOE                (requires merged data with genotype info)
    ├── Batch_Effects       (requires merged_combined_post_QC.csv)
    └── Triage_Metadata     (QC decisions vs specimen metadata)  ← NEW
```

## Getting Started

1. Open the relevant `.Rproj` file in RStudio
2. Input files are auto-copied from upstream outputs if not present
3. Run the Rmd scripts; outputs will be saved to `output_files/`
