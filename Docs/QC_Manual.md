# Alamar Biomarker QC — Operator Manual

How to run the U19 Alamar biomarker Primary QC pipeline end to end, starting from a GitHub
pull. Written so a new operator can follow it on **their own computer** without prior context.

This manual lives in `Docs/` inside the `Biomarker-QC-and-Analysis` repo (so it travels with a
clone), with a Word/PDF export for non-coders.

> ⚠️ **The GitHub repo is PUBLIC.** Screenshots marked `📷 [SCREENSHOT]` show patient sample IDs
> (export tables, QC details), so the image files are **not** committed here — they are kept in the
> internal Box working folder
> `Alamar_QC_Protocol/Documentation_in_Progress/QC_Manual_Cowork_06_2026/screenshots/`. Reference
> them from there; never push patient-identifiable images to the public repo.

> Everything here runs in your own **local** working folder. You don't need access to anyone
> else's Box/Drive. (Where completed runs eventually get centralized for the group is a separate,
> later decision.)

> **Layout note (this group's actual setup):** runs live under
> `Alamar_QC_Protocol/QC_Runs/`, created with the in-repo **`scaffold_qc_run.R`** (see §1.3) —
> not the `~/AlamarQC` example path used below for a generic fresh operator.

---

## Section 1 — Get the code and set up a run

### 1.0 What you need once

- **Git** (`git --version`; macOS: install Xcode Command Line Tools or `brew install git`).
- **R ≥ 4.0** and **RStudio**.
- **R packages** — see the install block at the bottom of the repo `README.md`.
- **GitHub access** to `emartin-um/Biomarker-QC-and-Analysis` via an SSH key (preferred) or an
  HTTPS credential in the macOS keychain. Do **not** put a token in the clone URL.
- **Alamar account** with access to ARGO Command Center and the NULISA Analysis Software.

### 1.1 Make a working folder and clone the repo (once)

Pick any local folder you control — here we use `~/AlamarQC`. Clone the repo inside it:

```bash
mkdir -p ~/AlamarQC && cd ~/AlamarQC
git clone git@github.com:emartin-um/Biomarker-QC-and-Analysis.git Alamar_Biomarker_QC_Repo
# HTTPS alternative (authenticates via keychain, no token in the URL):
# git clone https://github.com/emartin-um/Biomarker-QC-and-Analysis.git Alamar_Biomarker_QC_Repo
```

The clone is your **canonical code** — the source of truth. You may pull updates here, but you
do **not** run QC inside it; runs happen in separate folders (1.3), which keeps it clean.

> 📷 [SCREENSHOT — Terminal: successful clone output]

### 1.2 One-time folder setup

In your working folder, create two siblings of the clone:

```bash
cd ~/AlamarQC
mkdir -p QC_Runs Archive
```

Layout on your machine:

```
~/AlamarQC/
├── Alamar_Biomarker_QC_Repo/   ← canonical clone (code only, kept clean)
├── QC_Runs/                    ← your per-batch runs live here
└── Archive/                    ← zipped completed runs
```

### 1.3 Create a run folder for this batch

Each run gets its own folder under `QC_Runs/`, holding a **frozen copy of the repo code** so the
run can never overwrite the repo or another run. Create it with the in-repo helper
**`scaffold_qc_run.R`**, run from the repo root. It copies code only (no data/outputs), stamps
`RUN_INFO.txt` with the repo git rev for provenance, and points Metadata_Merge's `output_dir` at
this run:

```bash
cd ~/AlamarQC/Alamar_Biomarker_QC_Repo   # the canonical clone
# Name the run; pass the modules you'll use (Utilities is needed for prepare_inputs.R):
Rscript scaffold_qc_run.R "QC_7plate_06_17_26" "Primary_QC,Metadata_Merge,Secondary_QC,Utilities"
```

This creates `QC_Runs/QC_7plate_06_17_26/` with each module's code + an empty `input_files/`. (An
older draft, `new_qc_run.R`, used a `git worktree` instead; we use `scaffold_qc_run.R` because a
copy is safe inside the Box-synced tree — a nested worktree `.git` can fight Box sync — and because
it never needs a published tag.) A run folder looks like:

```
QC_Runs/QC_7plate_06_17_26/
├── RUN_INFO.txt                         ← repo git rev this run was frozen from
├── Primary_QC/
│   ├── QC_Pipeline_Primary_Alamar.Rmd   ← frozen copy of the repo code
│   ├── input_files/                      ← Alamar download + prepared CSVs go HERE
│   └── output_files/                     ← results land here (created on knit)
└── Utilities/
    └── prepare_inputs.R                  ← builds the input CSVs (§2.5)
```

> Note: the **Primary QC** Rmd is self-contained (it does not `source()` anything from
> `Utilities/`). `Utilities/` is copied into the run only because **`prepare_inputs.R`** (§2.5)
> lives there. To pick up later repo code changes, re-scaffold a fresh run rather than editing the
> frozen copy in place.

> 📷 [SCREENSHOT — Finder: the new QC_Runs/<batch>/Primary_QC structure]

### 1.4 Point the Alamar download at this run's `input_files/`

Downloads go **into the run folder, not `~/Downloads`** — so the raw export sits next to the QC
that uses it:

1. Chrome → **Settings → Downloads**.
2. Set **Location** to this run's input folder, e.g.
   `~/AlamarQC/QC_Runs/my_batch_2026-06/Primary_QC/input_files`
3. Or turn **on** "Ask where to save each file before downloading" and pick it per download.

> 📷 [SCREENSHOT — Chrome download settings pointed at the run's input_files]

---

## Section 2 — Export the data from Alamar and build the input files

Two Alamar web apps are involved:

- **ARGO Command Center** — `commandcenter.alamarbio.com`. Lists the instrument **runs** (one
  per Bay = one plate). Used to confirm which plates are ready.
- **NULISA Analysis Software (NAS)** — `analysis.alamarbio.com`. Combines selected runs into a
  named **Analysis** and produces the report (`NPQ_Values.xlsx` + `QC_Report.html`).

The pipeline reads **six** CSVs from `input_files/`. Five come from one export
(`NPQ_Values.xlsx`) via `prepare_inputs.R`; the sixth (`Replicate_list.csv`) is analyst-made.

### 2.1 Identify the completed plates (ARGO Command Center)

Sign in to `commandcenter.alamarbio.com`; you land on the **RUNS** tab. Each row is one Bay run
= **one plate**. A plate is ready only when **both `Run Status` and `Analysis Status` show
`Completed`** — `Pending` analysis is not ready, `Cancelled` is excluded. (The AD assay is
`CNS Disease Panel 120 (V2)`, project `ALZ`.)

Example — the June 2026 batch is these **8** plates:

| Run | Bays |
|---|---|
| 20260518-1415 | Bay1, Bay2, Bay3 |
| 20260514-1222 | Bay1, Bay2, Bay3 |
| 20260513-1625 | Bay1 |
| 20260512-1714 | Bay1 |

(The 20260519-1245 Bays show `Analysis Status = Pending`, so they wait for a later batch.)

> 📷 [SCREENSHOT — Command Center RUNS list, the 8 plates with both statuses Completed]

### 2.2 Create the analysis (NULISA Analysis Software)

Sign in to `analysis.alamarbio.com`. The **Analysis List** shows past analyses; **a first-time
user's list is empty** — that's expected.

1. Click the blue **+ New Analysis** button. The **START A NEW ANALYSIS** panel opens.
2. Enter an **Analysis Name** (required). It becomes part of the download folder name
   (`YYYYMMDD_<name>_NULISA_Analysis_Software_Report`), so make it meaningful, e.g.
   `ALZ_8plate_2026Jun`.
3. Choose the **Import from ACC** tab — this pulls runs straight from ARGO Command Center
   (ACC). (The other tab, *Upload Local Files*, is for files already on disk.)
4. Under **ARGO Command Center Projects**, expand the project (e.g. **ALZ**) and tick the **8
   completed runs** from 2.1. The list is paginated (e.g. "1–10 of 46 rows"); use the filter
   row (Project / Run / Assay / Created on) to narrow it. Each ticked run appears in **Run
   Preview** with its plate map and well counts (samples + Negative / Inter-Plate / Sample
   controls) — a good per-plate sanity check.
5. Click **Load Selected (N)** at the bottom right. The button shows the count, so confirm it
   reads **8** before clicking.

> ⚙️ After loading: how the analysis is generated and the report downloaded — to be filled in
> from the live capture.
>
> 📷 [SCREENSHOT — New Analysis panel: name + Import from ACC tab + project list]
> 📷 [SCREENSHOT — the 8 runs selected, shown in Run Preview]

### 2.3 Open the analysis and download the report

After **Load Selected**, the analysis is created and appears at the top of the **Analysis
List**; the **Analysis Preview** panel summarizes it — Runs (8), **Targets (131)**, **QC
Warnings**, **Wells (768** = samples + Negative/Inter-Plate/Sample controls**)**, Covariates
(`plateID`). Sanity check: plates × 96 = wells (8 × 96 = 768). Click **View Details** to open it.

The Details view opens on the **QUALITY CONTROL** tab (sidebar: Summary / Visualization / Assay
Precision / Detectability; tabs: Quality Control / Heatmap / Statistics). The Run/Sample table
shows per-plate QC — all 8 runs should read **Run QC = Passed**, alongside per-plate Warnings,
IC/IPC CV, Run Detectability, and Reads. This is **Alamar's own QC** (captured in
`QC_Report.html`); our pipeline runs its own QC downstream, so a high Alamar warning count on a
plate is informational, not a stop.

Click **DOWNLOAD REPORT** (top right). (The separate **Download QC Tables** button exports
Alamar's per-run/per-sample QC tables — useful reference, not a pipeline input.)

**DOWNLOAD REPORT** opens an export-configuration page with three columns:

- **Main Report** — the plots/tables that build `QC_Report.html` (Plate Layout, Read Summary,
  Heatmap, Quality Control, Detectability, Intra/Interplate Normalization, Box plots,
  Correlation, Clustering, PCA, Batch Effect). Leave **all checked** (default).
- **Statistical Reports** — Differential Expression / Pathway Enrichment. Greyed out until a
  statistical model is run; **not needed** for QC, leave off.
- **Data Files** — the pipeline needs only **NPQ Data File** (it produces `NPQ_Values.xlsx`).
  **XML file(s)** (raw per-run XML) and **NPQ Data Long Format** (~100k-row tidy table) are
  **not** read by the pipeline; unchecking them keeps the export small. `prepare_inputs.R` uses
  four sheets — `NPQ`, `Raw-counts`, `Annotation-Targets`, `Sample QC` — so confirm
  `NPQ_Values.xlsx` has those tabs. (Verified: dropping **XML** leaves `NPQ_Values.xlsx` intact
  at ~1.9 MB. **NPQ Long Format** is also unused by the current pipeline — see the refactor note
  at the end of Section 2.)

⚠️ **Then set both `Samples` and `Targets` to `All`** (the default is "Selected", which exports
only the highlighted sample/target — not the full dataset). With both on `All` the row estimate
jumps to the full count (≈100,608 = 131 targets × 768 wells). Click **Export**.

The export downloads as a **`.zip`** that unzips to a report folder named
`YYYYMMDD_<name>_NULISA_Analysis_Software_Report/`. The key file is:

- `NPQ_Values.xlsx` — sheets: `NPQ`, `Raw-counts`, `Annotation-Targets`, `Annotation-Samples`,
  `Run QC`, `Sample QC`. This is the only file the pipeline reads.
- `QC_Report.html` — Alamar's own QC summary (keep for reference; not consumed by the pipeline).

Because you set Chrome's download location in **1.4**, the zip lands in the run's
`input_files/`.

> 📷 [SCREENSHOT — Details view: the DOWNLOAD REPORT button, top right]
> 📷 [SCREENSHOT — export config: Data Files checked, Samples = All, Targets = All, Export]
> 📷 [SCREENSHOT — the downloaded report folder inside input_files/]

### 2.4 Unzip and confirm what landed

The report arrives as a **`.zip`** (e.g. `20260610_ALZ_8_plate_10_06_2026_..._Report.zip`,
~26 MB). Unzip it in place. With **all** Data Files boxes checked it expands to:

```
input_files/
└── 20260610_<name>_NULISA_Analysis_Software_Report/
    ├── <8 per-run>.xml        ← raw, 1.5 MB each — NOT used by the pipeline
    ├── NPQ_Long_Format.csv    ← ~33 MB — NOT used by the pipeline
    ├── NPQ_Values.xlsx        ← the file the pipeline needs
    └── QC_Report.html         ← Alamar QC summary (reference)
```

If you exported **NPQ Data File only**, the folder is just `NPQ_Values.xlsx` (+ `QC_Report.html`
if Main Report was kept) — a fraction of the size.

> 📷 [SCREENSHOT — Finder: the downloaded .zip, then the unzipped report folder contents]

### 2.5 Build the five derived CSVs with `prepare_inputs.R`

`prepare_inputs.R` reads `NPQ_Values.xlsx` and writes five files into `input_files/`,
replacing the old manual "Save-As each sheet in Excel" step:

| Output | Source sheet | Notes |
|---|---|---|
| `NPQ_<date>.csv` | `NPQ` | biomarkers × samples |
| `Raw_counts.csv` | `Raw-counts` | biomarkers × samples |
| `Annotation_Targets.csv` | `Annotation-Targets` | target annotation |
| `Sample_QC.csv` | `Sample QC` | header rows kept (pipeline reads it with `skip = 12`) |
| `biomarker_detectability.csv` | derived | `targetDetectability`×100 pivoted by plate + `Overall` |

In RStudio, with the working directory set to the run's `Primary_QC/` (open
`Primary_QC.Rproj`), run in the console:

`prepare_inputs.R` is in the repo `Utilities/` (the scaffold copies it into each run). In RStudio,
with the working directory set to the run's `Primary_QC/` (open `Primary_QC.Rproj`), run in the
console:

```r
source("../Utilities/prepare_inputs.R")
prepare_inputs(
  xlsx    = "input_files/YYYYMMDD_<name>_NULISA_Analysis_Software_Report/NPQ_Values.xlsx",
  out_dir = "input_files"
)
```

The date in `NPQ_<date>.csv` is taken from the report folder's leading `YYYYMMDD`; pass
`npq_label = "YYYYMMDD"` to set it explicitly.

> **New-export gotchas (seen on the 2026-06 batches — handle these):**
> - **Loose files, no report folder.** When you export *NPQ Data File only*, the zip may unzip to
>   a bare `NPQ_Values.xlsx` + `QC_Report.html` (no `YYYYMMDD_<name>_…_Report/` wrapper). Then
>   there's no folder date to infer, so **pass `npq_label = "YYYYMMDD"`** (or recreate the dated
>   folder and put the two files in it, which also preserves provenance).
> - **`APOE`/`CRP` detectability arrives as the literal text `"NA"`.** These relative-quant targets
>   have no detectability, and some exports write `"NA"` (not an empty cell), which makes the whole
>   `targetDetectability` column character. `prepare_inputs.R` coerces with `as.numeric()` so they
>   become real `NA` (carried through faithfully) — if you see a `* 100` "non-numeric argument"
>   error, you're on an old copy of the script; re-scaffold from the current repo.

> 📷 [SCREENSHOT — RStudio console: prepare_inputs() output]

### 2.6 Add the replicate list (analyst-maintained)

`prepare_inputs.R` does **not** create `Replicate_list.csv` — it is the list of HIHG replicate
sample IDs for this batch, which you maintain. Create `input_files/Replicate_list.csv` as one
column of sample IDs under a `SampleID` header. `prepare_inputs.R` only *warns* if it is missing,
but the **Primary QC pipeline reads it directly and will error without it** — so it is required.

In practice the HIHG replicates are the **`202505211-*`** wells present in the batch (the repeated
QC sample, one ID per replicate well). Each batch carries its own slice of that series (ALZ123 used
`-31…`, the 8-plate used `-321…-338`, the 7-plate used `-343…-361`). You can list them straight
from the NPQ file:

```r
npq  <- readr::read_csv("input_files/NPQ_<date>.csv", n_max = 1, show_col_types = FALSE)
reps <- sort(grep("^202505211", names(npq), value = TRUE))
writeLines(c("SampleID", reps), "input_files/Replicate_list.csv")
```

Confirm the resulting set is the complete intended replicate list for the batch before knitting.

### 2.7 Verify before running QC

`input_files/` should now hold the six inputs plus the raw export:

```
input_files/
├── YYYYMMDD_<name>_NULISA_Analysis_Software_Report/   ← raw export (kept)
├── NPQ_YYYYMMDD.csv
├── Raw_counts.csv
├── Annotation_Targets.csv
├── Sample_QC.csv
├── biomarker_detectability.csv
└── Replicate_list.csv
```

Quick sanity checks: `NPQ` and `Raw_counts` have a `targetName` column plus one column per
sample; `Sample_QC.csv` still has its leading header rows (the pipeline skips 12); and
`biomarker_detectability.csv` has `Biomarker`, `Plate_01…`, and `Overall`.

> **Expected blanks:** `APOE` and `CRP` have **blank** detectability and LOD — they are
> relative-quant targets (`Curve_Quant = "R"`) for which Alamar computes neither. This is
> faithful to the export, not an extraction error; the pipeline retains them. (Heads-up: the
> detectability filter does `min(..., na.rm = TRUE)`, which returns `Inf` for an all-blank row,
> so these pass the filter via `Inf` and emit a harmless warning — a candidate for an explicit
> guard later.)

> **Refactor note (later):** the NAS **NPQ Long Format** export (one tidy row per
> sample × target, ~100k rows) could replace the wide multi-sheet inputs. Confirmed (Jun 2026)
> to carry the per-observation fields — `PlateID`, `SampleType`, `Target`, `NPQ`,
> `UnnormalizedCount` (raw count), `LOD`, `Sample_QC_Detectability`, and per-metric
> `Sample_QC_*_Status` flags. So a single read could derive wide NPQ/Raw_counts, sample
> metadata (Run/Bay/matrix), control flags, and even the Alamar QC-warning droplist — retiring
> `prepare_inputs.R`'s multi-sheet extraction and the brittle `Sample_QC` `skip = 12`.
> Per-target-per-plate detectability would be recomputed as `mean(NPQ ≥ LOD)` per (PlateID,
> Target) and validated against the workbook's `targetDetectability`. A bounded refactor of the
> pipeline's wide-format reads — worth doing at the next internals pass, not a routine-run
> change.

---

## Section 3 — Configure, knit, and read the QC report

### 3.1 Open the run (not the repo)

Open `QC_Runs/<batch>/Primary_QC/Primary_QC.Rproj` in RStudio. The `.Rproj` sets the working
directory to this run's `Primary_QC/`, so `input_files/` and `../Utilities/shared_functions.R`
resolve and outputs go to this run's `output_files/`. **Do not open the repo's `Primary_QC.Rproj`
— knitting there reads the wrong inputs and writes into the repo.**

> A run folder needs the code (`Primary_QC/*.Rmd` + `Utilities/`) present. Cleanest is a
> `git worktree` of a tagged commit (carries provenance). But a worktree only sees *committed*
> code — if you have **uncommitted** changes you want in this run (e.g. a just-applied guard),
> copy the current `Primary_QC/QC_Pipeline_Primary_Alamar.Rmd`, `Primary_QC.Rproj`, and
> `Utilities/` into the run folder instead.

> 📷 [SCREENSHOT — RStudio with the run's Primary_QC.Rproj open]

### 3.2 Set inputs and thresholds (Section A.3)

In the Rmd's **A.3 Initial Setup** chunk:

- `NPQ_data <- "NPQ_YYYYMMDD.csv"` — your NPQ file in `input_files/`.
- Thresholds (defaults are fine for most runs):

| Param | Default | Meaning |
|---|---|---|
| `PCA_SD` | 5 | SDs from mean to call a PCA outlier |
| `min_detectability` | 50 | min % detectability to keep a biomarker in the post-QC set |
| `min_detectability_qc` | 98 | min % detectability for the QC Biomarker Set |
| `corr_thresh` | 0.4 | max correlation allowed in the QC Biomarker Set |
| `read_count_threshold` | 500 | min mean/median raw reads |
| `samp_out_thresh` | 1.5 | IQR multiplier for sample NPQ outliers |
| `read_out_thresh` | 4 | IQR multiplier for read outliers |
| `FDR_threshold` | 0.01 | outlier-burden significance |

Also set the YAML `title` and the output HTML name for the batch.

> 📷 [SCREENSHOT — Section A.3 with NPQ_data set]

### 3.3 Knit

Click **Knit** (or `rmarkdown::render("QC_Pipeline_Primary_Alamar.Rmd")`). It runs ~5–15 min;
the HTML report opens when done, and CSV outputs are written to `output_files/`.

> 📷 [SCREENSHOT — the knit progress / finished HTML report]

### 3.4 What the pipeline does (so the report makes sense)

1. Builds a **QC Biomarker Set** — biomarkers with ≥ `min_detectability_qc` (98%) detectability on
   every plate and pairwise correlation ≤ `corr_thresh`. Relative-quant targets (APOE, CRP) are
   retained and flagged, not detectability-filtered.
2. Uses that set to find **sample outliers** via PCA (`PCA_SD`) and outlier burden
   (`FDR_threshold`) → **triage** (removed) and **flagged** (kept but noted).
3. Reintegrates all biomarkers, splitting **normal** vs **low-detectability / low-count**
   (`min_detectability` 50%, `read_count_threshold` 500).
4. Writes post-QC datasets split by biomarker group and by triaged / non-triaged samples.

### 3.5 Read the report and key outputs

Open the HTML and check, in order:

| Priority | File | Check |
|---|---|---|
| HIGH | `samples_to_triage.csv` | which samples were removed, and why |
| HIGH | `NPQ_YYYYMMDD_post_QC.csv` | your main analysis dataset |
| MED | `flagged_read_outliers.csv` | samples kept but worth a look |
| MED | `biomarkers_low_detectability.csv` | biomarkers separated out |
| LOW | `Alamar_QC_warn.csv` | Alamar's own warnings |

Control / replicate pulls are written separately: `hihg_reps_NPQ.csv`, `nc_reps_NPQ.csv`,
`sc_reps_NPQ.csv`, `ipc_reps_NPQ.csv`. Confirm the **relative-quant note** (APOE, CRP) appears in
the report and that both are present in the post-QC data.

> 📷 [SCREENSHOT — the QC report: summary + triage section]

### 3.6 Output files reference

- **Main analysis:** `NPQ_YYYYMMDD_post_QC.csv` (QC'd samples, normal biomarkers);
  `NPQ_YYYYMMDD_Low_post_QC.csv` (low-detectability/count biomarkers — consider a
  detected/not-detected analysis); the `*_triage.csv` versions hold the removed samples for
  sensitivity checks.
- **Ancillary** (fixed names, regenerated every run): `QC_Biomarker_Set.csv`,
  `sample_pca_coordinates_round{1,2}.csv`, `upper/lower_outliers_*.csv`, `lod_matrix.csv`,
  `lod_reads_matrix.csv`, `biomarkers_low_counts.csv`, `bad_samples.csv`, `raw_df.csv`.

### 3.7 Archive the completed run

When satisfied, zip the run's `input_files/` + `output_files/` + rendered HTML + the Rmd into
`Archive/<batch>_<date>.zip`, with a `RUN_INFO.txt` recording the code commit/tag. (See the
migration/organization plan.) The next ~10 plates are the natural test: a second operator
follows this manual end to end.
