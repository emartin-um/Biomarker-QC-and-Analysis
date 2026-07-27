#!/usr/bin/env Rscript
# =============================================================================
# test_intensity_columns.R  ·  2026-07-27
# -----------------------------------------------------------------------------
# Acceptance test for the intensity columns (addendum §9-§10). Base R only, so
# it runs anywhere. Exits non-zero on the first failure.
#
#   Rscript test_intensity_columns.R                    # test the current outputs
#   Rscript test_intensity_columns.R <baseline_dir>     # also diff against a baseline
#
# The baseline defaults to output_files/_pre_intensity_cols_2026Jul27/ if present.
#
# WHY A TEST AND NOT JUST A DIFF: output_files/ is gitignored, so checking out
# main does NOT revert the CSVs. Branch-switching cannot tell you whether the
# change was safe — only a content check can.
# =============================================================================

HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
OUT  <- file.path(HERE, "output_files")
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args)) args[1] else file.path(OUT, "_pre_intensity_cols_2026Jul27")

NEW_COLS <- c("mean_INT", "mean_INT_z", "intensity_flag",
              "global_outlier_dir", "global_outlier_high", "global_outlier_low")

fails <- 0L
ok <- function(label, cond, detail = "") {
  pass <- isTRUE(cond)
  if (!pass) fails <<- fails + 1L
  cat(sprintf("  [%s] %s%s\n", if (pass) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0("  -- ", detail) else ""))
}

# In the source repo there is no output_files/ at all — code and docs only, no
# data. That is not a failure; there is simply nothing to check yet. Only treat a
# missing FILE as a failure when the directory exists (i.e. a run that should
# have produced it).
if (!dir.exists(OUT)) {
  cat("\nNo output_files/ here — this is the source repo, not a QC run.\n")
  cat("Nothing to validate. Run this inside a QC run directory after the screen.\n")
  quit(status = 0L)
}

for (f in c("extreme_sample_master.csv", "clinician_review_shortlist.csv")) {
  cat("\n===", f, "===\n")
  p <- file.path(OUT, f)
  if (!file.exists(p)) { ok(paste("exists:", f), FALSE); next }
  m <- read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)

  # --- 1. the new columns are present ---------------------------------------
  ok("all six new columns present", all(NEW_COLS %in% names(m)),
     paste("missing:", paste(setdiff(NEW_COLS, names(m)), collapse = ", ")))
  if (!all(NEW_COLS %in% names(m))) next

  # --- 2. internal consistency of the direction split -----------------------
  go <- m$global_outlier
  ok("direction is set on every global outlier, and only on them",
     all(nzchar(m$global_outlier_dir[go])) && all(!nzchar(m$global_outlier_dir[!go])))
  ok("global_outlier_high/low are subsets of global_outlier",
     all(!m$global_outlier_high | go) && all(!m$global_outlier_low | go))
  ok("high + low + mixed accounts for every global outlier",
     sum(m$global_outlier_high) + sum(m$global_outlier_low) +
       sum(go & m$global_outlier_dir == "MIXED") == sum(go),
     sprintf("%d + %d + %d vs %d", sum(m$global_outlier_high), sum(m$global_outlier_low),
             sum(go & m$global_outlier_dir == "MIXED"), sum(go)))
  ok("HIGH means the high tail dominates",
     all(m$n_extreme_high[m$global_outlier_high] > m$n_extreme_low[m$global_outlier_high]))
  ok("LOW means the low tail dominates",
     all(m$n_extreme_low[m$global_outlier_low] > m$n_extreme_high[m$global_outlier_low]))
  ok("high and low are mutually exclusive", !any(m$global_outlier_high & m$global_outlier_low))

  # --- 3. the direction split is not cosmetic -------------------------------
  if (sum(m$global_outlier_high) > 1 && sum(m$global_outlier_low) > 1) {
    hz <- mean(m$mean_INT_z[m$global_outlier_high], na.rm = TRUE)
    lz <- mean(m$mean_INT_z[m$global_outlier_low],  na.rm = TRUE)
    ok("high-driven outliers really are brighter than low-driven", hz > lz + 1,
       sprintf("mean_INT_z %+.2f vs %+.2f", hz, lz))
  }

  # --- 4. intensity_flag agrees with mean_INT_z -----------------------------
  ok("intensity_flag matches mean_INT_z at +/-2.5",
     all((m$intensity_flag == "HIGH") == (m$mean_INT_z >=  2.5), na.rm = TRUE) &&
     all((m$intensity_flag == "LOW")  == (m$mean_INT_z <= -2.5), na.rm = TRUE))
  ok("mean_INT_z is standardised (mean ~0, sd ~1)",
     abs(mean(m$mean_INT_z, na.rm = TRUE)) < 0.05 &&
       abs(stats::sd(m$mean_INT_z, na.rm = TRUE) - 1) < 0.05)

  # --- 5. the flag is additive, not a rediscovery ---------------------------
  if (f == "extreme_sample_master.csv") {
    flagged <- !is.na(m$intensity_flag) & nzchar(m$intensity_flag)
    ok("intensity flag adds cases the global-outlier flag misses",
       sum(flagged & !go) > 0,
       sprintf("%d intensity-flagged, %d of them not global outliers",
               sum(flagged), sum(flagged & !go)))
  }

  # --- 6. nothing pre-existing moved ----------------------------------------
  b <- file.path(BASE, f)
  if (file.exists(b)) {
    o <- read.csv(b, stringsAsFactors = FALSE, check.names = FALSE)
    ok("row count unchanged vs baseline", nrow(o) == nrow(m),
       sprintf("%d -> %d", nrow(o), nrow(m)))
    shared <- intersect(names(o), names(m))
    moved <- shared[!vapply(shared, function(v) isTRUE(all.equal(o[[v]], m[[v]])), logical(1))]
    ok("no pre-existing column changed vs baseline", length(moved) == 0,
       if (length(moved)) paste("changed:", paste(moved, collapse = ", ")) else "")
    ok("the only added columns are the six", setequal(setdiff(names(m), names(o)), NEW_COLS),
       paste("added:", paste(setdiff(names(m), names(o)), collapse = ", ")))
  } else {
    cat("  [skip] no baseline at ", b, "\n", sep = "")
  }
}

cat(sprintf("\n%s  (%d failure%s)\n",
            if (fails == 0L) "ALL CHECKS PASSED" else "FAILED",
            fails, if (fails == 1L) "" else "s"))
quit(status = if (fails == 0L) 0L else 1L)
