#!/usr/bin/env Rscript
library(tidyverse)
library(readxl)

# Build Primary_QC input_files/ from one Alamar NAS NPQ_Values.xlsx, replacing the
# manual "Save-As each sheet" step. Replicate_list.csv is analyst-made (HIHG
# replicate IDs) and is only checked, never generated.

prepare_inputs <- function(xlsx, out_dir = "input_files", npq_label = NULL) {

  stopifnot(file.exists(xlsx))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(npq_label)) npq_label <- str_extract(basename(dirname(xlsx)), "^[0-9]{8}")
  if (is.na(npq_label))   stop('Could not infer date; pass npq_label = "YYYYMMDD".')

  npq <- read_excel(xlsx, "NPQ", .name_repair = "minimal")
  npq_file <- sprintf("NPQ_%s.csv", npq_label)
  write_csv(npq, file.path(out_dir, npq_file))
  write_csv(read_excel(xlsx, "Raw-counts", .name_repair = "minimal"),
            file.path(out_dir, "Raw_counts.csv"))
  at <- read_excel(xlsx, "Annotation-Targets", .name_repair = "minimal")
  write_csv(at, file.path(out_dir, "Annotation_Targets.csv"))

  # Keep the 12 preamble rows so the pipeline's read_csv(skip = 12) still lands on the header
  read_excel(xlsx, "Sample QC", col_names = FALSE, .name_repair = "minimal") |>
    write_csv(file.path(out_dir, "Sample_QC.csv"), col_names = FALSE)

  # biomarker_detectability: targetDetectability*100, pivoted by plateID, plates sorted and
  # relabeled Plate_NN, rows in NPQ order, Overall = row mean (verified exact vs ALZ123).
  # as.numeric() guards relative-quant targets (APOE, CRP): some exports write a literal "NA"
  # for their blank detectability, which types the whole column as character -> "* 100" would
  # error. Coercing turns those into real NA (carried through, faithful to the export).
  det <- at |>
    distinct(targetName, plateID, .keep_all = TRUE) |>
    transmute(targetName, plateID,
              pct = round(suppressWarnings(as.numeric(targetDetectability)) * 100, 1)) |>
    pivot_wider(names_from = plateID, values_from = pct)
  plates <- sort(setdiff(names(det), "targetName"))
  det <- det |>
    select(targetName, all_of(plates)) |>
    arrange(match(targetName, npq[[1]])) |>
    rename_with(~ sprintf("Plate_%02d", seq_along(plates)), all_of(plates)) |>
    rename(Biomarker = targetName)
  det$Overall <- round(rowMeans(select(det, starts_with("Plate_")), na.rm = TRUE), 1)
  write_csv(det, file.path(out_dir, "biomarker_detectability.csv"))

  if (!file.exists(file.path(out_dir, "Replicate_list.csv")))
    warning("Replicate_list.csv missing in ", out_dir,
            " — analyst-created HIHG replicate IDs; add it before running QC.")

  message("Wrote ", npq_file,
          " + Raw_counts, Annotation_Targets, Sample_QC, biomarker_detectability -> ", out_dir)
  invisible(out_dir)
}

# Example:
# prepare_inputs(
#   "QC_Runs/8plate_2026-06/Primary_QC/input_files/<NAS_report>/NPQ_Values.xlsx",
#   out_dir = "QC_Runs/8plate_2026-06/Primary_QC/input_files")
