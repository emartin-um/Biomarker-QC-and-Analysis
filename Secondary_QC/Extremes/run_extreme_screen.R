#!/usr/bin/env Rscript
# =============================================================================
# run_extreme_screen.R
# -----------------------------------------------------------------------------
# Driver for the Extreme-Value Screen (Secondary QC). Loads the U19 and
# OtherProjects preprocessed datasets, pools them for a single robust reference
# distribution, flags per-marker extremes (both tails) and disease-signature
# panels, and writes:
#
#   output_files/extreme_reference_stats.csv   per-marker median/MAD/quantiles
#   output_files/extreme_marker_long.csv       one row per (sample,marker) extreme
#   output_files/extreme_sample_master.csv     one row per sample (all flags)
#   output_files/clinician_review_shortlist.csv ranked cases + plain-language why
#   output_files/panel_top_<PANEL>.csv         top samples per disease panel
#
# Run headless:  Rscript run_extreme_screen.R
# (run from the Extremes/ module directory)
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

# --- paths (module dir is .../Secondary_QC/Extremes; run root is ../../) ------
HERE     <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                     error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
RUN_ROOT <- normalizePath(file.path(HERE, "..", ".."), mustWork = FALSE)
DATA_DIR <- file.path(RUN_ROOT, "datasets")
OUT_DIR  <- file.path(HERE, "output_files")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(HERE, "intensity_columns.R"))   # must precede the helpers
source(file.path(HERE, "extreme_helpers.R"))

# --- thresholds (tunable; INT = rank-based inverse-normal, ~standard-normal) --
THR_INT         <- 2.75 # |INT| for a single-marker extreme event (~top/bot 0.3%)
PANEL_HIT       <- 1.8  # oriented INT for a marker to count as a panel "hit" (~top 3.6%)
PANEL_SCORE_THR <- 2.0  # mean oriented INT to call a panel signature
PANEL_MIN_HIT   <- 2    # >= this many markers in the tail for a signature
LONE_INT_THR    <- 3.0  # |INT| on one disease marker -> shortlist on its own (~top 0.13%)

# --- load + pool -------------------------------------------------------------
u19o <- readRDS(file.path(DATA_DIR, "preprocessed_data_U19.rds"))
otho <- readRDS(file.path(DATA_DIR, "preprocessed_data_OtherProjects.rds"))

u19 <- u19o$data; u19$dataset <- "U19"
oth <- otho$data; oth$dataset <- "OtherProjects"
common <- intersect(names(u19), names(oth))
dat <- bind_rows(u19[, common], oth[, common])
biomarker_cols <- u19o$biomarker_cols
cat(sprintf("Pooled cohort: %d samples (U19 %d + Other %d); %d biomarker cols\n",
            nrow(dat), nrow(u19), nrow(oth), length(biomarker_cols)))

stopifnot(!any(duplicated(dat$SAMPLE)))   # SAMPLE is the unique key

# --- run the shared compute pipeline -----------------------------------------
res <- run_extreme_pipeline(dat, biomarker_cols, thr_int = THR_INT, panel_hit = PANEL_HIT,
                            panel_score_thr = PANEL_SCORE_THR, panel_min_hit = PANEL_MIN_HIT,
                            lone_int_thr = LONE_INT_THR)
mk <- res$markers; ref <- res$ref; long_ext <- res$long_ext
master <- res$master; short <- res$short
cat(sprintf("Usable markers: %d / %d (excluded: %s)\n",
            length(mk), length(biomarker_cols), paste(EXCLUDE_MARKERS, collapse = ", ")))

# --- write outputs -----------------------------------------------------------
write.csv(ref,      file.path(OUT_DIR, "extreme_reference_stats.csv"),    row.names = FALSE)
write.csv(long_ext, file.path(OUT_DIR, "extreme_marker_long.csv"),        row.names = FALSE)
write.csv(master,   file.path(OUT_DIR, "extreme_sample_master.csv"),      row.names = FALSE)
write.csv(short,    file.path(OUT_DIR, "clinician_review_shortlist.csv"), row.names = FALSE)

# Per-panel top lists (samples with a computable score, ranked by score).
for (pn in names(PANELS)) {
  sc <- paste0(pn, "_score")
  if (!sc %in% names(master)) next
  tp <- master[!is.na(master[[sc]]), ]
  tp <- tp[order(-tp[[sc]]), ]
  keep <- intersect(c("SAMPLE","dataset","CDX_collapsed","age_at_subject","sex",
                      sc, paste0(pn, "_nhit"), "top_extreme","HBA1_INT","frac_hemo_high",
                      "RUN","ID2"), names(tp))
  write.csv(head(tp[, keep], 40), file.path(OUT_DIR, sprintf("panel_top_%s.csv", pn)),
            row.names = FALSE)
}

# --- console summary (for threshold tuning) ----------------------------------
cat("\n================ SUMMARY ================\n")
cat(sprintf("Extreme events (|INT|>=%g): %d across %d distinct samples\n",
            THR_INT, nrow(long_ext), length(unique(long_ext$SAMPLE))))
cat(sprintf("Samples with >=1 extreme: %d (%.1f%%)\n",
            sum(master$n_extreme_total > 0), 100*mean(master$n_extreme_total > 0)))
cat(sprintf("Global outliers (>=8 markers extreme): %d  (high-driven %d, low-driven %d, mixed %d)\n",
            sum(master$global_outlier), sum(master$global_outlier_high),
            sum(master$global_outlier_low),
            sum(master$global_outlier & master$global_outlier_dir == "MIXED")))
cat(sprintf("Whole-sample intensity (mean_INT) beyond +/-2.5 SD: %d HIGH, %d LOW  [informational]\n",
            sum(master$intensity_flag == "HIGH", na.rm = TRUE),
            sum(master$intensity_flag == "LOW",  na.rm = TRUE)))
cat(sprintf("Clinician shortlist: %d samples (HIGH=%d, MEDIUM=%d)\n",
            nrow(short), sum(short$priority == "HIGH"), sum(short$priority == "MEDIUM")))
cat("\n-- shortlist by dataset --\n");      print(table(short$dataset))
cat("\n-- shortlist by CDX --\n");          print(table(short$CDX_collapsed, useNA = "ifany"))
cat("\n-- shortlist by priority x dataset --\n"); print(table(short$priority, short$dataset))
cat("\n-- panel signature hits ((score>=thr & nhit>=", PANEL_MIN_HIT, ") OR nhit>=3) --\n", sep = "")
for (pn in names(PANELS)) {
  sc <- master[[paste0(pn, "_score")]]; nh <- master[[paste0(pn, "_nhit")]]
  thr_pn <- panel_thr(pn, PANEL_SCORE_THR)
  hit <- (!is.na(sc) & sc >= thr_pn & nh >= PANEL_MIN_HIT) | (!is.na(nh) & nh >= 3)
  cat(sprintf("  %-20s %d samples (thr %.1f)\n", pn, sum(hit), thr_pn))
}
cat("\n-- top 15 most-extreme markers by frequency of extreme events --\n")
print(head(sort(table(long_ext$Biomarker), decreasing = TRUE), 15))
cat("\nOutputs written to: ", OUT_DIR, "\n")
