# =============================================================================
# intensity_columns.R  ·  added 2026-07-27
# -----------------------------------------------------------------------------
# Two additions to the per-sample master, both recommended in
#   Batch_Effects/2026_Batch_Investigation/QC_PIPELINE_RECOMMENDATIONS_ADDENDUM_2026-07-27.md
# (§9 and §10). Pure base R and file-free, so the headless backfill script can
# reuse it without dplyr — there must be exactly ONE definition of these columns.
#
# 1. WHOLE-SAMPLE INTENSITY (`mean_INT`).
#    The screen already builds the full INT matrix and then only ever asks COUNT
#    questions of it (n_extreme_high/low/total, max_abs_INT). It never asks a
#    LEVEL question. A sample that is moderately high on ALL markers at once
#    therefore stays invisible: verified 2026-07-27, mean max_abs_INT is 2.50 in
#    the middle decile of background intensity and only 2.96 in the brightest, so
#    no single marker ever reaches the tail.
#    `mean_INT` is that missing level question. It correlates r = 0.952 with the
#    per-sample background factor G derived independently in AD_Prediction_MCI_OD,
#    and separates its ALL3-flagged subjects from unflagged at AUC 0.860.
#
#    ⚠ INFORMATIONAL, NOT A DROP RULE. Background intensity is a continuous
#    nuisance covariate affecting every sample, not a population of bad wells.
#    It belongs in the analysis model. The flag exists to make the axis visible.
#
# 2. DIRECTION OF THE GLOBAL-OUTLIER CALL (`global_outlier_dir`).
#    `global_outlier = n_extreme_total >= 8` merges two opposite populations,
#    because its two components move in opposite directions with intensity:
#    across deciles of G, n_extreme_high runs 0.09 -> 1.94 while n_extreme_low
#    runs 1.08 -> 0.18. Of the 28 global outliers in the U19 analysis frame, 18
#    are high-driven (mean G +1.12) and 9 low-driven (mean G -0.94) — a bright
#    sample and a dim sample carrying one label.
#    `global_outlier` itself is UNCHANGED (downstream consumers depend on it);
#    the direction is added alongside.
# =============================================================================

#' Add whole-sample intensity and global-outlier direction to the master.
#'
#' @param master  per-sample master data frame; must already carry
#'                `n_extreme_high`, `n_extreme_low` and `global_outlier`.
#' @param int     the samples x markers INT matrix used to build `master`,
#'                row-aligned with it.
#' @param sd_mult flag threshold for `intensity_flag`, in SDs of `mean_INT`.
#' @return `master` with mean_INT, mean_INT_z, intensity_flag,
#'         global_outlier_dir, global_outlier_high, global_outlier_low.
add_intensity_columns <- function(master, int, sd_mult = 2.5) {
  need <- c("n_extreme_high", "n_extreme_low", "global_outlier")
  miss <- setdiff(need, names(master))
  if (length(miss))
    stop("add_intensity_columns(): master is missing ", paste(miss, collapse = ", "))
  if (nrow(int) != nrow(master))
    stop(sprintf("add_intensity_columns(): int has %d rows, master has %d — not aligned",
                 nrow(int), nrow(master)))

  # --- 1. whole-sample intensity ---------------------------------------------
  # Mean over the INT matrix. Every marker is on the same standard-normal rank
  # scale, so the mean is directly interpretable as "how high does this sample
  # read across the panel" and is comparable between samples.
  mi <- rowMeans(int, na.rm = TRUE)
  mi[!is.finite(mi)] <- NA_real_
  s  <- stats::sd(mi, na.rm = TRUE)
  z  <- if (is.finite(s) && s > 0) (mi - mean(mi, na.rm = TRUE)) / s else rep(NA_real_, length(mi))

  # Derive the flag from the ROUNDED z, i.e. from the value that actually gets
  # written. Flagging on the full-precision z left one boundary sample (true z
  # -2.4996, written -2.50) flagged by `mean_INT_z <= -2.5` but not by
  # `intensity_flag == "LOW"` — two columns in the same file disagreeing about
  # the same sample. Caught by test_intensity_columns.R on 2026-07-27.
  z <- round(z, 2)

  master$mean_INT       <- round(mi, 3)
  master$mean_INT_z     <- z
  master$intensity_flag <- ifelse(is.na(z), NA_character_,
                           ifelse(z >=  sd_mult, "HIGH",
                           ifelse(z <= -sd_mult, "LOW", "")))

  # --- 2. direction of the global-outlier call --------------------------------
  # NB `global_outlier` gates on the TOTAL, so a sample can be a global outlier
  # with neither tail reaching 8 on its own (e.g. 5 high + 4 low). Direction is
  # therefore "which tail dominates", not "which tail crossed a bar".
  go  <- !is.na(master$global_outlier) & master$global_outlier
  dir <- ifelse(!go, "",
         ifelse(master$n_extreme_high > master$n_extreme_low,  "HIGH",
         ifelse(master$n_extreme_low  > master$n_extreme_high, "LOW", "MIXED")))

  master$global_outlier_dir  <- dir
  master$global_outlier_high <- go & dir == "HIGH"
  master$global_outlier_low  <- go & dir == "LOW"

  # --- 3. panel scores read against the sample's own level (addendum §11) -----
  # Every pathology panel scores as the MEAN ORIENTED INT across its markers, and
  # six of the seven orient every marker +1 (the seventh orients both -1). Sign
  # balance is therefore ±1.00 for all seven — the maximum possible exposure to a
  # common shift. Measured against the background factor: AD_tau 0.64,
  # BD_pTau_inflation 0.67, Neuroinflammation 0.64, Neurodegeneration 0.52.
  #
  # Because they all load the same way, a bright sample drifts up on ALL of them
  # together — so apparent MULTI-AXIS pathology is the specific pattern this
  # artifact produces, and n_pathology_axes is the field most exposed to it.
  # Subtracting the sample's own level asks what a reviewer actually wants to
  # know: is this axis raised RELATIVE TO THIS SPECIMEN?
  #
  # Raw scores are left untouched; these are additional columns.
  #
  # `best_panel_score` is deliberately EXCLUDED: it is a max across panels, not a
  # panel, so `best_panel_score - mean_INT` is not the same thing as "the best of
  # the adjusted scores" and would quietly mislead. Take the max over the *_adj
  # columns if you want that.
  sc <- setdiff(grep("_score$", names(master), value = TRUE), "best_panel_score")
  sc <- sc[vapply(master[sc], is.numeric, logical(1))]
  for (s in sc) master[[paste0(s, "_adj")]] <- round(master[[s]] - mi, 3)

  master
}
