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
- **Whole-panel intensity** as a covariate (§D.0) and a **self-audit of the QC screens** (§D.5)

---

## The pipeline audits itself (added 2026-07-30)

Implements the recommendations in `Secondary_QC/Triage_Review/`. **No sample's
fate changed** — triage and flag membership are identical to the previous
version. What changed is what the pipeline measures, and what it says about
itself.

### §D.0 — `mean_INT`, measured but never gated

Every screen in the pipeline asks a **count** question ("on how many markers is
this sample extreme?") or a **distance** question ("how far from the plate
centre?"). None asks a **level** question. So a sample that reads moderately high
on *all* markers at once is invisible to every one of them: the outlier-burden
test counts markers outside a per-plate IQR fence, and a uniform shift never
makes any single marker extreme. That is structural, not a threshold that can be
loosened.

`mean_INT` is that missing level statistic — a per-marker rank-based
inverse-normal transform, averaged per sample, so every marker contributes on the
same footing.

Two design choices worth knowing:

* It is computed **before triage**, so removed wells land on the same scale as
  retained ones. Run downstream it could only ever see the survivors, and removed
  wells could not be placed on the axis at all.
* It is **reported, never gated.** Background intensity is a continuous nuisance
  covariate affecting every sample, not a population of bad wells — gating on it
  would delete real variation. It belongs in the analysis model. Same definition
  as `Secondary_QC/Extremes/intensity_columns.R`; there is one definition of this
  column, deliberately.

### §D.5 — the screen audit

Diagnostic only; drops nothing. Asks of each gate: does it fire at all, is it
independent of the others, where does its cut sit on its own statistic, and is it
keying on plate position rather than on the specimen.

Four things it surfaces that a pass/fail count hides:

| | |
|---|---|
| **INERT gates** | A gate firing on **zero** wells has not been shown to be safe — it has not been *tested*. Headroom to each cut is printed beside the count, so "correctly found nothing" and "cannot fire" can be told apart. |
| **The read gate is one-sided** | Its lower bound is not plate-relative on any real run and **cannot be made so**: within-plate read spread is wide enough that any fence avoiding a mass of false flags falls below the physical floor. §D.2.4 records the alternatives that were tested and why each fails. The gate is *described* as a plate-relative high flag plus an absolute floor, which is what it has always been; behaviour is unchanged. |
| **The IC gate is a vendor metric re-derived** | We test `IC Reads` against 0.6–1.4× the plate median; Alamar's `IC Median` is the same check on the same quantity. §D.2.3 measures the agreement each run and returns a verdict. A useful cross-check — **not** an independent in-house line of evidence, and methods text should not present it as one. |
| **`FDR_threshold` is not a dial** | The p-value is Poisson-binomial on an **integer** anchor count, so the FDR takes only a handful of distinct values. "FDR < x" really means "k or more anchors out". §D.2.2 prints every setting that exists; between adjacent rows there is nothing to tune. |

**Per-axis triage flags.** `samples_to_triage.csv` carries `tri_IC`, `tri_PCA`,
`tri_burden`, `tri_bad_data` and `n_axes`. A report that assigns one reason per
well by precedence is a valid *partition* — the loss must sum — but its per-axis
numbers are **not axis totals**, and the gap can be large. Use the booleans for
"how many wells did this screen catch"; use precedence only when the parts must
sum to the whole. §D.5 prints both side by side.

**Plate position — and the other meaning of "site".** §D.5.4 tests whether a
*physical well position* is triaged more often than chance across plates, by
permutation within plate-bay (a chi-square is wrong here: triage counts spread
over ~84–96 positions leave most cells at 0 or 1). This is invisible to every
other screen, because they all compare a well only to its own plate — a bad
position is a layout or hardware problem, not a specimen problem.

The spatial tests are **layout-agnostic and name no region in advance**: a row
gradient, a column gradient, an edge-vs-interior contrast, and a descriptive
report of the worst quadrant. To *test* a specific region, set `prior_region` in
§A.3 and run on the **next** dataset — testing a region chosen by looking at the
same data is a search result dressed as a p-value, and the pipeline will not do
it for you.

> **Findings live in the rendered report, not here.** Each run's HTML
> (`QC_Pipeline_Primary_Alamar_<NPQ>.html`) carries §D.0 and §D.5 with that run's
> numbers.

> **Clustering by clinical collection site is a different question**, and it
> cannot be asked here: Primary QC is agnostic to phenotype and covariates by
> design, and the site label only exists after the metadata merge. That is what
> `well_level_QC.csv` is for — one row per well, every gate flag, keyed by
> `SampleID`, so the downstream join is one line. See
> `Secondary_QC/Triage_Review/`.

**Key Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `PCA_SD` | 5 | Number of SDs from mean to define PCA outliers |
| `min_detectability` | 50 | Minimum detectability (%) to retain biomarker in post-QC dataset |
| `min_detectability_qc` | 98 | Minimum detectability (%) to retain biomarker in QC Biomarker dataset |
| `read_count_threshold` | 500 | Minimum mean/median raw reads per sample for QC |
| `corr_thresh` | 0.4 | Maximum correlation coefficient allowed between biomarkers in QC Biomarker dataset |
| `scaffold_mode` | `"frozen"` | QC Biomarker Set: `"frozen"` (fixed 21-anchor scaffold + per-plate tripwire, default since 2026-06-20) or `"rederive"` (per-run derive + de-correlate). Only the triage anchors change; post-QC biomarker content is identical. |
| `samp_out_thresh` | 1.5 | Multiplier for IQR to determine if a sample is an outlier from the median NPQ for a biomarker|
| `read_out_thresh` | 4 | Multiplier for IQR to determine if a sample is a **high**-read outlier from the plate median |
| `read_floor_abs` | 500000 | Absolute low-read floor — the **only** lower read bound. See §D.2.4 for why a plate-relative one cannot work |
| `FDR_threshold` | 0.01 | Significance threshold for outlier burden. Quantised — see §D.2.2 |
| `intensity_sd_flag` | 2.5 | SDs of `mean_INT_z` for the intensity flag. **Informational only — never triages** |
| `ic_median_thresh` | 0.40 | Alamar's own `IC Median` bar, used in §D.2.3 to test whether our IC gate re-derives it |
| `n_perm` | 10000 | Permutation replicates for the plate-position tests in §D.5.4 |

**Required Input Files:**
- Raw NPQ data files (CSV format) from Alamar platform
- Sample metadata with Run/Bay information
- `Sample_QC.csv` must carry a well position (`AUTO_WELLPOSITION`, or `wellRow` +
  `wellCol`) for the §D.5.4 position diagnostics. Without it the run continues
  and those diagnostics are skipped with a warning.

**Output Files:**
- `NPQ_[date]_post_QC.csv` - Cleaned NPQ data for high-quality samples/biomarkers
- `biomarker_summary.csv` - QC statistics for each biomarker
- `sample_summary.csv` - QC statistics for each sample
- Various QC plots (PNG/PDF)

**Audit and covariate outputs** (§D.0 and §D.5 — none of these drops a sample):

| file | what it answers |
|---|---|
| `well_level_QC.csv` | **one row per well**: every gate flag, `mean_INT`, position, reads, IC. The join key for anything downstream — above all clustering by *clinical collection site* |
| `whole_panel_intensity.csv` | `mean_INT` / `mean_INT_z` for all wells, triaged included |
| `qc_gate_inventory.csv` | does each gate fire, or is it **INERT** (untested, not proven safe), and how much headroom to its cut |
| `qc_gate_duplication.csv` | is our IC gate the vendor's `IC Median` re-derived (Jaccard + verdict) |
| `qc_gate_overlap.csv` | pairwise Jaccard between screens — which are one screen reported twice |
| `triage_axis_overlap.csv` | true per-axis totals vs what a one-reason-per-well report shows |
| `triage_axis_combinations.csv` | the full Venn, flattened |
| `burden_operating_points.csv` | the complete list of `FDR_threshold` settings that exist |
| `gate_cut_neighbourhood.csv` | nearest retained and nearest removed well on each axis |
| `read_gate_geometry.csv` | per plate-bay, whether the lower read bound is plate-relative or the absolute floor |
| `qc_axis_attribution.csv` | each triage class on the intensity axis (d, AUC) — **per class, never pooled** |
| `well_position_effects.csv` / `_readflag.csv` | per-position rates + permutation p |
| `position_map_IC.csv` | median IC recovery by plate position |
| `ic_row_gradient.csv` | IC recovery down the plate, per plate-bay |

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
