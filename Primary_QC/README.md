# Primary QC

Primary quality control analyses for Alamar biomarker data. This pipeline processes raw NPQ data, applies QC filters, and produces cleaned datasets for downstream analysis.

## Usage

1. Open `Primary_QC.Rproj` in RStudio
2. Place input data files in `input_files/`
3. Run analysis scripts
4. Outputs will be saved to `output_files/`

## Scripts

### QC_Pipeline_Primary_Alamar.Rmd

Main QC pipeline for Alamar fluid biomarker data.

> **Report filename (auto-stamped):** the YAML `knit:` function names the rendered HTML after the NPQ input dataset — e.g. `QC_Pipeline_Primary_Alamar_NPQ_20251220.html` — so each run's report is unique and doesn't overwrite the previous one. Just Knit normally (RStudio's Knit button honors this); no need to rename the `.Rmd`. Falls back to a date stamp if no NPQ file is in `input_files/`.

> **Robust to Alamar export drift:** the `Sample_QC.csv` header row is **auto-detected** (no hardcoded `skip=`), and the run stops with a clear message if expected columns (`plateID`, `Sample Name`, `QC Status`) are missing. Value-format drift is handled: `IC Median`/`Detectability` are coerced from percent strings (`"95.4%"`) or numerics, `IC Reads`/`Reads` tolerate thousands separators, and `biomarker_detectability.csv` is rescaled to the 0–100 scale the thresholds assume if an export delivers 0–1 proportions. `QC Status` matching is case-insensitive.

**Purpose:** Performs comprehensive quality control including:
- Sample and biomarker filtering based on read counts and detectability
- Identification of high-read/high-detectability biomarker subsets
- LOD (Limit of Detection) analysis
- Plate/batch effect visualization
- NPQ data normalization and export

**Key Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `PCA_SD` | 5 | Number of SDs from mean to define PCA outliers |
| `min_detectability` | 50 | Minimum detectability (%) to retain biomarker in post-QC dataset |
| `min_detectability_qc` | 98 | Minimum detectability (%) to retain biomarker in QC Biomarker dataset |
| `read_count_threshold` | 500 | Minimum mean/median raw reads per sample for QC |
| `corr_thresh` | 0.4 | Maximum correlation coefficient allowed between biomarkers in QC Biomarker dataset |
| `samp_out_thresh` | 1.5 | Multiplier for IQR to determine if a sample is an outlier from the median NPQ for a biomarker|
| `read_out_thresh` | 4 | Multiplier for IQR to determine if a sample is an outlier from the median for Total Reads|
| `FDR_threshold` | 0.01 | Significance threshold for outlier burden |

**Required Input Files:**
- Raw NPQ data files (CSV format) from Alamar platform
- Sample metadata with Run/Bay information

**Output Files:**
- `NPQ_[date]_post_QC.csv` - Cleaned NPQ data for high-quality samples/biomarkers
- `biomarker_summary.csv` - QC statistics for each biomarker
- `sample_summary.csv` - QC statistics for each sample
- Various QC plots (PNG/PDF)

## Input Files

Place your input data in `input_files/` (not tracked in git):
- Raw NPQ data and raw count data from Alamar runs
- Biomarker annotation and detectability data
- Alamar Sample QC
- List of any non-Alamar technical controls

## Output Files

Generated outputs are saved to `output_files/` (not tracked in git):
- Filtered NPQ datasets
- QC summary files
- Diagnostic plots

## Dependencies

See main repository README for full package list. Key packages:
- `tidyverse`, `readr` - Data manipulation
- `pheatmap` - Heatmap visualization
- `kableExtra` - Table formatting
