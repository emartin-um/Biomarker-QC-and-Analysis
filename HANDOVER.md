# HANDOVER — Alamar Biomarker QC (repo + QC_Runs)

**Date:** 2026-06-15  ·  **Branch:** `metadata-update-n2850-OtherProjects-03052026`
**Last commit:** `c9726a8` (pushed to `origin`/github.com/emartin-um/Biomarker-QC-and-Analysis)

This covers the **repo** (`R_Code_Dev/Alamar_Biomarker_QC_Repo/`) and the **`QC_Runs/`** working copies
(`Alamar_QC_Protocol/QC_Runs/`, outside the repo). Read alongside the per-module READMEs and
`QC_Runs/README.txt`.

---

## 1. Model (how this is organized)

- The **repo holds canonical CODE only** — `.gitignore` excludes all `input_files/`, `output_files*/`,
  `*.csv/xlsx/rds/html/RData`, and `_archive*/`/`_backup*/`. **0 data files are tracked** (no PHI in git).
- **Real batches run in their own folder** under `QC_Runs/QC_<batch>_<date>/`, created with
  `Rscript scaffold_qc_run.R "QC_<batch>_<date>"` (copies code only, freezes the repo git rev in
  `RUN_INFO.txt`). Each run's data/outputs are gitignored. **Never render in the repo itself.**
- Run-copies differ from the repo **only** in 2 run-specific YAML lines (`output_dir`, `review_dir`;
  Primary also: title + `NPQ_data`). Re-scaffold or re-copy a file to pick up repo fixes.

## 2. Pipelines

**Primary_QC** (`QC_Pipeline_Primary_Alamar.Rmd`) — sample/biomarker QC. Knit in RStudio; the `knit:`
function auto-names the HTML after the NPQ dataset. Robust to Alamar export drift (Sample_QC header
auto-detected; `IC Median`/`Detectability` parsed from `%` strings or numerics; detectability rescaled
to 0–100 if delivered as 0–1).

**Metadata_Merge** (`Metadata_Merge_Pipeline.Rmd`) — merges post-QC NPQ with metadata. Knit (auto-named
HTML). Key behavior:
- **Duplicates/longitudinal:** metadata-driven cross-check vs the lab list (`Datasets/Dup_sample_or_ind_*.xlsx`).
  Duplicates keep the most-recent run **provisionally**, overridable by the lab via
  `<review_dir>/Duplicate_review.xlsx` (`FINAL_keep="keep"`). Both runs of each duplicate →
  `filtered/filtered_combined_duplicates.csv` (replicate feed). Longitudinal split → `filtered_combined_longitudinal.csv`.
- **WGS APOE:** Step 6d masks the ~282 swap-affected AD-Hispanic samples to **Sanger** (see §5);
  `APOE.geno_final = coalesce(WGS, Sanger)` otherwise.
- **Reports** (`<output_dir>/reports/`, all `.xlsx`): top level = `QC_error_report_combined` / `_by_sample`
  / `_summary` + mapping keys + crosstab; sub-folders `errors/`, `missing_data/`, `duplicates/`
  (Duplicate_review + Longitudinal_review).
- Params: `input_dir`, `output_dir`, `review_dir`, `known_repeats_file`, `filter_only`, `include_derived_*`.
- **Test without real data:** `Metadata_Merge/example/` (run `make_example_inputs.R`, then render with
  `input_dir/output_dir` params — covariates-off smoke test).
- Helper scripts (repo): `combine_metadata.R`, `covariate_explorer.R`, `fix_age_suppression.R`.

## 3. QC_Runs — current runs

| Run | What | Status |
|---|---|---|
| `QC_8plate_10_06_26/` | Production 8-plate batch. Primary: `NPQ_20260610` → 655×111. Metadata: `n4328_2026Jun12`, review_2026Jun (data-driven), merged 632 → **512** cross-sectional. | Rendered 2026-06-15 (consistent w/ 06-10). `Metadata_Merge_RUNBOOK.md`. |
| `QC_ALZ123_repro_2026Jun/` | Reproduction of the original ALZ123 dataset. Primary: `NPQ_20251220` → 2891×103 (identical to prior). Metadata: clean `n2871_May2026`, merged 2825 → **2615**. Plus standalone WGS re-map (§5). | Rendered 2026-06-15. `RUNBOOK.md`, `RUN_INFO.txt`. |

Both runs verified consistent with prior results; pre-rerun snapshots kept (`_results_pre_rerun_*`,
`_backup_code_pre_Rmd_*`).

## 4. What changed this session (committed `c9726a8`)

Pipeline `.qmd`→`.Rmd`; duplicate/longitudinal QC; reports → all-`.xlsx` + combined + by-sample +
summary, organized into sub-folders; Primary_QC export-drift robustness; `scaffold_qc_run.R`;
synthetic `example/`; removed the refuted WGS plate-swap correction (script + outputs + metadata);
archived stale repo outputs into `Metadata_Merge/_archive_stale_outputs_2026-06-15/`; docs updated.

## 5. WGS plate-swap status (IMPORTANT)

The May-2026 WGS delivery mislabeled ~282 AD-Hispanic samples. **The 2026-06-08 geometric correction
was REFUTED, reverted, and deleted.** Current understanding (2026-06-11, chr1 concordance): the **WGS**
is the mislabeled platform (MIDI re-array). **Canonical pipeline = mask affected samples to Sanger** (no
relabel). A candidate re-map exists — `QC_Runs/QC_ALZ123_repro_2026Jun/WGS_Midi_Remap/remap_wgs_midi_plateswap.R`
(WGS-vs-Sanger agreement 47%→99.6%) — kept **with its run, not the repo**, as a **PROPOSAL pending USUHS
(Clifton) confirmation**. For affected samples, **prefer Sanger (`APOE.geno`)**. See
`../../May_2026_WGS_QC/WGS_Plate_Swap_Summary_2026Jun11.md`.

---

## 6. PENDING / NEXT STEPS  ⟵ (reminders)

- [ ] **Commit the rest** (deferred by request 2026-06-15). The following pre-existing, uncommitted
  material was **NOT** in commit `c9726a8` — review and commit (or gitignore) separately:
  - Modified (pre-session): `Metadata_Merge/combine_metadata.R`, `Secondary_QC/APOE/APOE_geno_protein.qmd`,
    `Secondary_QC/Hemolysis_Check/Hemolysis_by_Site.Rmd`
  - Untracked: `Secondary_QC/APOE/{investigate_discordance,november_plate_investigation,random_mixup_test,biomarker_batch_run_investigation}.R`,
    5× `Secondary_QC/APOE/*.png`, `Secondary_QC/APOE/output_files_pre_plateswap_2026Jun/` (figures/data),
    `Secondary_QC/Batch_Effects/` (module + a `.Rbak` backup + data dirs),
    `Secondary_QC/Hemolysis_Check/{Hemolysis_Collaborator_Summary,Hemolysis_Internal_Notes}.md`
  - ⚠ Decide per-file: the PNGs + `output_files_pre_plateswap/` figures are output artifacts; the
    `.Rbak` is a backup; `Hemolysis_Internal_Notes.md` may be internal — consider excluding these.
  - [ ] Also commit **this `HANDOVER.md`** (left uncommitted with the above).
- [ ] **WGS re-map** awaits **USUHS (Clifton)** midi-plate confirmation before applying. Until then,
  affected samples use Sanger.
- [ ] **Lab decisions** pending: fill `reports/duplicates/Duplicate_review.xlsx` (`FINAL_keep`) to pick
  which duplicate run to keep, after comparing biomarker values via `Replicate_Analysis` on
  `filtered_combined_duplicates.csv`.
- [ ] **Secondary QC** (Replicate_Analysis, etc.) deferred until the next combined/8-plate batch.
- [ ] Optional: generalize `QC_Runs/QC_8plate_10_06_26/prepare_inputs.R` into a repo utility; make the
  `example/` work covariates-on (generator currently writes `Group="NA"` which readr reads as NA).

## 7. References
- Repo READMEs: top `README.md`; `Primary_QC/README.md`; `Metadata_Merge/README.md` (+ `example/README.md`);
  `Secondary_QC/README.md` + subdir READMEs.
- `QC_Runs/README.txt` (run log) + each run's `RUNBOOK.md`.
- WGS investigation: `../../May_2026_WGS_QC/`.
