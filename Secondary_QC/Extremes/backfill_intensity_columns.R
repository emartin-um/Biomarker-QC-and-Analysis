#!/usr/bin/env Rscript
# =============================================================================
# backfill_intensity_columns.R  ·  2026-07-27
# -----------------------------------------------------------------------------
# Adds the 2026-07-27 columns (mean_INT + global-outlier direction) to the
# ALREADY-WRITTEN outputs, without re-running the whole screen.
#
# WHY THIS EXISTS: run_extreme_screen.R requires dplyr, which is not installed on
# every machine that needs these outputs current. This script is BASE R ONLY. It
# rebuilds the INT matrix exactly as the screen does — same pooled cohort, same
# usable-marker rule, same int_scores() formula — and calls the SAME
# add_intensity_columns() the pipeline calls, so the two cannot drift.
#
# It is a backfill, not a replacement: once the screen is next run in an
# environment with dplyr, it emits these columns itself and this script is a
# no-op. Every other column is left untouched, and that is asserted.
#
# Run from the Extremes/ module directory:  Rscript backfill_intensity_columns.R
# =============================================================================

HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
RUN_ROOT <- normalizePath(file.path(HERE, "..", ".."), mustWork = FALSE)
DATA_DIR <- file.path(RUN_ROOT, "datasets")
OUT_DIR  <- file.path(HERE, "output_files")

source(file.path(HERE, "intensity_columns.R"))   # the single definition

## --- the pieces of extreme_helpers.R we need, restated here ONLY because that
## --- file cannot be sourced without dplyr. These three must stay in step with it.
EXCLUDE_MARKERS <- c("APOE4", "APOE")            # extreme_helpers.R:56
int_scores <- function(v) {                      # extreme_helpers.R:206-213
  v <- suppressWarnings(as.numeric(v)); ok <- !is.na(v)
  out <- rep(NA_real_, length(v)); n <- sum(ok)
  if (n >= 2) out[ok] <- stats::qnorm((rank(v[ok], ties.method = "average") - 0.5) / n)
  out
}
usable_markers <- function(df, biomarker_cols) { # extreme_helpers.R:171-177
  keep <- setdiff(intersect(biomarker_cols, names(df)), EXCLUDE_MARKERS)
  keep[vapply(keep, function(b) {
    v <- suppressWarnings(as.numeric(df[[b]]))
    sum(!is.na(v)) >= 20 && stats::mad(v, na.rm = TRUE) > 0
  }, logical(1))]
}

## --- rebuild the pooled cohort exactly as run_extreme_screen.R does ----------
u19o <- readRDS(file.path(DATA_DIR, "preprocessed_data_U19.rds"))
otho <- readRDS(file.path(DATA_DIR, "preprocessed_data_OtherProjects.rds"))
u19 <- as.data.frame(u19o$data); u19$dataset <- "U19"
oth <- as.data.frame(otho$data); oth$dataset <- "OtherProjects"
common <- intersect(names(u19), names(oth))
dat <- rbind(u19[, common], oth[, common])          # = dplyr::bind_rows on shared cols
stopifnot(!any(duplicated(dat$SAMPLE)))

if (all(c("A_42", "A_40") %in% names(dat))) dat$A42_A40_ratio <- dat$A_42 - dat$A_40
mk <- usable_markers(dat, u19o$biomarker_cols)
if ("A42_A40_ratio" %in% names(dat)) mk <- c(mk, "A42_A40_ratio")
int <- vapply(mk, function(b) int_scores(dat[[b]]), numeric(nrow(dat)))
colnames(int) <- mk
cat(sprintf("Pooled %d samples x %d usable markers\n", nrow(dat), length(mk)))

## --- update each output in place --------------------------------------------
update_one <- function(fname) {
  path <- file.path(OUT_DIR, fname)
  if (!file.exists(path)) { cat("  skip (absent):", fname, "\n"); return(invisible(NULL)) }
  m <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  before <- m

  i <- match(m$SAMPLE, dat$SAMPLE)
  if (any(is.na(i)))
    stop(sprintf("%s: %d rows have no match in the pooled cohort", fname, sum(is.na(i))))

  m <- add_intensity_columns(m, int[i, , drop = FALSE])

  # assert we changed nothing else
  added <- setdiff(names(m), names(before))
  kept  <- intersect(names(before), names(m))
  bad   <- kept[!vapply(kept, function(v) isTRUE(all.equal(before[[v]], m[[v]])), logical(1))]
  if (length(bad)) stop(sprintf("%s: pre-existing column(s) changed: %s", fname,
                                paste(bad, collapse = ", ")))

  write.csv(m, path, row.names = FALSE)
  cat(sprintf("  %-32s %4d rows  +%d cols (%s)\n", fname, nrow(m), length(added),
              paste(added, collapse = ", ")))
  invisible(m)
}

cat("Updating:\n")
master <- update_one("extreme_sample_master.csv")
short  <- update_one("clinician_review_shortlist.csv")

## --- summary -----------------------------------------------------------------
if (!is.null(master)) {
  cat("\n================ SUMMARY ================\n")
  cat(sprintf("Global outliers: %d  (high-driven %d, low-driven %d, mixed %d)\n",
              sum(master$global_outlier), sum(master$global_outlier_high),
              sum(master$global_outlier_low),
              sum(master$global_outlier & master$global_outlier_dir == "MIXED")))
  cat(sprintf("  mean_INT_z among high-driven: %+.2f   low-driven: %+.2f\n",
              mean(master$mean_INT_z[master$global_outlier_high], na.rm = TRUE),
              mean(master$mean_INT_z[master$global_outlier_low],  na.rm = TRUE)))
  cat(sprintf("Whole-panel intensity beyond +/-2.5 SD: %d HIGH, %d LOW  [informational]\n",
              sum(master$intensity_flag == "HIGH", na.rm = TRUE),
              sum(master$intensity_flag == "LOW",  na.rm = TRUE)))
  ov <- table(intensity = master$intensity_flag != "",
              global_outlier = master$global_outlier)
  cat("\nOverlap with the existing global-outlier flag (rows = intensity-flagged):\n")
  print(ov)
  cat("\n  -> the intensity flag is additive: it is not a rediscovery of global_outlier.\n")
}
cat("\nNOTE: `review_caveat` is NOT refreshed here. This script adds columns; it does not\n")
cat("      re-run annotate_master(), so the caveat TEXT still reads the pre-2026-07-27 wording\n")
cat("      (no direction, no intensity note). The new columns above are authoritative in the\n")
cat("      meantime; the text catches up on the next full run of run_extreme_screen.R.\n")
cat("\nDone. Re-running run_extreme_screen.R (needs dplyr) reproduces these columns natively.\n")
