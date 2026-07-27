# =============================================================================
# specimen_quality_helpers.R  ·  added 2026-07-27
# -----------------------------------------------------------------------------
# The THIRD quality axis: per-specimen properties. Base R only, I/O-free.
#
# The pipeline already surfaces two technical axes — read depth (a SITE property)
# and detectability (a BLANK-WELL property). Neither describes the patient. This
# module makes the per-specimen axis visible as a standing panel:
#
#   sequencing depth   how deeply this well was sequenced        (site-structured)
#   cystatin C (CST3)  a renal-function surrogate                (patient physiology)
#   whole-panel level  how high this sample reads across the panel
#
# WHY IT MATTERS. Renal impairment raises many plasma proteins at once, so a site
# or ancestry group carrying more renal burden looks proteomically "sicker" on any
# same-signed score. That is a design-level confounder you want to SEE before
# modelling, not discover afterwards.
#
# 🚫 THIS MODULE IS DESCRIPTIVE. It defines no drop rule and no gate. In
# particular CST3 must NEVER be added as a covariate to remove the effect: it is
# itself raised in clinical AD, so conditioning on it deletes disease signal
# (MAPT loses 25% of its effect, pTau-181 17%, pTau-217 9%; log2(reads) costs
# pTau-217 2%). See DOCS/ANALYSIS_PLAN_batch_adjustment.md §3a.
#
# Recommendations implemented here: 2026-07-13 §2 (depth-by-site, never built),
# and the 2026-07-27 addendum §12 (CST3 distributions) and §6 (sign-balance guard).
# =============================================================================

## --- small helpers -----------------------------------------------------------
.num <- function(v) suppressWarnings(as.numeric(gsub("[%,]", "", as.character(v))))

#' % of a variable's variance that lies between the levels of a grouping factor.
#' Plain one-way R^2 — the same quantity the batch-effect module reports as ICC.
var_between <- function(y, g) {
  ok <- !is.na(y) & !is.na(g)
  if (sum(ok) < 3 || length(unique(g[ok])) < 2) return(NA_real_)
  summary(stats::lm(y[ok] ~ factor(g[ok])))$r.squared
}

#' Does `x` still explain `y` once `g` is partialled out of both?
#' If a relationship is really a site effect, the within-site version collapses.
within_group_r2 <- function(y, x, g) {
  ok <- !is.na(y) & !is.na(x) & !is.na(g)
  if (sum(ok) < 10 || length(unique(g[ok])) < 2) return(c(overall = NA_real_, within = NA_real_))
  yr <- stats::residuals(stats::lm(y[ok] ~ factor(g[ok])))
  xr <- stats::residuals(stats::lm(x[ok] ~ factor(g[ok])))
  c(overall = summary(stats::lm(y[ok] ~ x[ok]))$r.squared,
    within  = summary(stats::lm(yr ~ xr))$r.squared)
}

## --- 1. the per-specimen axes ------------------------------------------------

#' Assemble the axes from a prepared data frame.
#' @param d  one row per sample; needs Site and, where available, l2reads / CST3.
#' @param intensity  optional per-sample whole-panel level (e.g. mean_INT from the
#'   Extremes screen). Pass NULL if the screen has not been run.
specimen_axes <- function(d, intensity = NULL) {
  ax <- list()
  if ("l2reads" %in% names(d)) ax[["sequencing_depth"]] <- .num(d$l2reads)
  if ("CST3"    %in% names(d)) ax[["cystatin_C"]]       <- .num(d$CST3)
  if ("HBA1"    %in% names(d)) ax[["hemolysis_HBA1"]]   <- .num(d$HBA1)
  if (!is.null(intensity))     ax[["panel_level"]]      <- .num(intensity)
  if (!length(ax)) stop("specimen_axes(): none of l2reads / CST3 / HBA1 / intensity present")
  ax
}

#' Per-site (or per-ancestry, per-run…) distribution of each axis.
#' Distributions, not just means — a chronically low-input site is a shape, not a point.
axis_by_group <- function(ax, g, group_name = "Site", min_n = 5) {
  g <- as.character(g)
  do.call(rbind, lapply(names(ax), function(a) {
    y <- ax[[a]]
    sp <- split(y, g)
    sp <- sp[vapply(sp, function(v) sum(!is.na(v)) >= min_n, logical(1))]
    if (!length(sp)) return(NULL)
    data.frame(
      axis      = a,
      group_var = group_name,
      group     = names(sp),
      n         = vapply(sp, function(v) sum(!is.na(v)), integer(1)),
      median    = round(vapply(sp, stats::median, numeric(1), na.rm = TRUE), 4),
      q25       = round(vapply(sp, function(v) stats::quantile(v, .25, na.rm = TRUE), numeric(1)), 4),
      q75       = round(vapply(sp, function(v) stats::quantile(v, .75, na.rm = TRUE), numeric(1)), 4),
      sd        = round(vapply(sp, stats::sd, numeric(1), na.rm = TRUE), 4),
      stringsAsFactors = FALSE, row.names = NULL)
  }))
}

#' How much of each axis lies between the levels of each design variable.
axis_variance <- function(ax, d, vars = c("Site", "Run", "Bay", "Ancestry", "CDX_collapsed")) {
  vars <- intersect(vars, names(d))
  do.call(rbind, lapply(names(ax), function(a) {
    r <- vapply(vars, function(v) var_between(ax[[a]], d[[v]]), numeric(1))
    data.frame(axis = a, variable = vars, pct_between = round(100 * r, 2),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
}

#' Do the axes single out the SAME groups? Correlation of their group medians.
#' Depth and cystatin C being uncorrelated across sites is the finding that stops
#' "it is all one site effect" — they are different sites.
axis_cocluster <- function(ax, g, min_n = 20) {
  g <- as.character(g)
  keep <- names(which(table(g) >= min_n))
  if (length(keep) < 4) return(NULL)
  M <- vapply(ax, function(y) vapply(keep, function(k)
    stats::median(y[g == k], na.rm = TRUE), numeric(1)), numeric(length(keep)))
  if (is.null(dim(M))) return(NULL)
  cm <- suppressWarnings(stats::cor(M, use = "pairwise.complete.obs"))
  out <- expand.grid(axis_a = rownames(cm), axis_b = colnames(cm),
                     stringsAsFactors = FALSE)
  out <- out[out$axis_a < out$axis_b, , drop = FALSE]
  out$r_across_groups <- round(mapply(function(a, b) cm[a, b], out$axis_a, out$axis_b), 3)
  out$n_groups <- length(keep)
  out[order(-abs(out$r_across_groups)), , drop = FALSE]
}

#' Does each axis still track the panel level WITHIN a group?
#' Retention near 100% means the link is per-specimen, not a group artifact.
axis_within_group <- function(ax, g, target = "panel_level") {
  if (!target %in% names(ax)) return(NULL)
  y <- ax[[target]]
  others <- setdiff(names(ax), target)
  do.call(rbind, lapply(others, function(a) {
    r <- within_group_r2(y, ax[[a]], g)
    data.frame(axis = a, target = target,
               pct_overall = round(100 * r["overall"], 2),
               pct_within  = round(100 * r["within"], 2),
               pct_retained = round(100 * r["within"] / r["overall"], 1),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
}

## --- 2. the sign-balance guard (2026-07-27 addendum, §6 amendment) -----------

#' Exposure of a weighted score to a common additive shift on the z-scale.
#'
#' Any score of the form w1*z1 + w2*z2 + ... moves under "every marker up by 1"
#' by exactly sum(w). Normalising by sum(|w|) gives a unit-free coefficient:
#'   0 = weights fully self-cancelling, immune
#'   1 = weights all one sign, fully exposed
#'
#' The original §6 rule was "any index built from DEPTH-loaded markers tracks
#' depth". The stronger rule: any index built from SAME-SIGNED markers tracks the
#' sample's overall level, whatever the markers are — and you can check it before
#' touching data. All seven Extremes pathology panels sit at 1.00.
sign_balance <- function(w) {
  w <- w[is.finite(w)]
  if (!length(w) || sum(abs(w)) == 0) return(NA_real_)
  sum(w) / sum(abs(w))
}

#' Screen a list of named weight vectors and flag the exposed ones.
sign_balance_report <- function(weight_list, warn_at = 0.5) {
  do.call(rbind, lapply(names(weight_list), function(nm) {
    w <- weight_list[[nm]]
    b <- sign_balance(w)
    data.frame(score = nm, n_markers = length(w),
               sum_w = round(sum(w), 4), sum_abs_w = round(sum(abs(w)), 4),
               sign_balance = round(b, 3),
               verdict = if (is.na(b)) NA_character_
                         else if (abs(b) >= warn_at)
                           "EXPOSED - report against the sample's own level"
                         else "largely self-cancelling",
               stringsAsFactors = FALSE, row.names = NULL)
  }))
}
