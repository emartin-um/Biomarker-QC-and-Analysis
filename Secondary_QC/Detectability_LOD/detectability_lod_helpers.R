# =============================================================================
# detectability_lod_helpers.R  ·  added 2026-07-27
# -----------------------------------------------------------------------------
# Implements 2026-07-13 recommendations §2 (detectability-by-bay panel) and §8
# (variance-stabilised LOD). Base R only, I/O-free.
#
# WHAT DETECTABILITY ACTUALLY MEASURES, and why that matters:
#
#   detectability = (# targets with NPQ > targetLOD_NPQ[target, bay]) / n_targets
#   targetLOD_NPQ = log2( mean(2^NC) + 3 * SD(2^NC) )
#
# and the NC term comes from just **four blank wells per plate-bay**. So a bay
# "fails" detectability when ITS BLANKS READ HIGH, not when its samples fail. On
# 20260519-1245 Bay2 and 20260514-1222 Bay3 the patients score 0-5% below a
# normal LOD. That makes detectability an ASSAY-BACKGROUND signal, not a
# sample-quality one — an analysis and Alamar question, never a re-run trigger.
#
# TWO STRUCTURAL FAILURE MODES, both addressed here:
#
#   1. `mean + 3*SD` on n = 4 is fragile. The 3*SD term inherits the sampling
#      error of a four-point SD, so ONE high blank can swing a bay's LOD — and
#      its detectability rate — by tens of points. `lod_table()` therefore also
#      returns a variance-stabilised LOD that keeps the per-bay MEAN (which is
#      the real background signal) and replaces the per-bay SD with the target's
#      typical across-bay SD.
#
#   2. The gate is a knife-edge. Detectability is quantised at k/n_targets and
#      the 0.90 line falls between k = 116 (0.8992, FAIL) and k = 117 (0.9070,
#      PASS) — it bisects the mode, so a tiny background shift flips a large
#      share of a bay's wells.
#
# 🚫 NEVER GATE OR DROP ON DETECTABILITY, and never gate a marker on per-bay LOD
# detection status without first checking that marker's blank trend across the
# batch. In 2026, 26% of ALL wells fall below their own bay's pTau-217 LOD (2025:
# ~1-2%) purely because the blank floor rose, while pTau-217 QUANTITATION is
# unchanged. Use continuous NPQ.
# =============================================================================

.lod_from_nc <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 2) return(c(mean = NA_real_, sd = NA_real_, n = length(v)))
  c(mean = mean(v), sd = stats::sd(v), n = length(v))
}

#' Per plate-bay x target LOD, as shipped and variance-stabilised.
#'
#' @param nc  blank-well NPQ: one row per NC well, with a `plate_bay` column and
#'   one numeric column per target.
#' @param markers  target column names.
#' @param k  SD multiplier (Alamar uses 3).
#' @return data frame: plate_bay, target, nc_n, nc_mean_linear, nc_sd_linear,
#'   lod_asis, lod_stabilised, sd_typical, delta.
lod_table <- function(nc, markers, k = 3) {
  pb <- as.character(nc$plate_bay)
  out <- do.call(rbind, lapply(unique(pb), function(p) {
    idx <- which(pb == p)
    st  <- vapply(markers, function(t) .lod_from_nc(2^nc[[t]][idx]), numeric(3))
    data.frame(plate_bay = p, target = markers,
               nc_n = st["n", ], nc_mean_linear = st["mean", ], nc_sd_linear = st["sd", ],
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  out$lod_asis <- log2(out$nc_mean_linear + k * out$nc_sd_linear)

  # Variance stabilisation: keep the per-bay MEAN (the real background signal),
  # replace the four-point SD with the target's TYPICAL across-bay SD. The median
  # is used rather than the mean so one pathological bay cannot inflate every
  # other bay's LOD.
  typ <- tapply(out$nc_sd_linear, out$target, stats::median, na.rm = TRUE)
  out$sd_typical     <- as.numeric(typ[out$target])
  out$lod_stabilised <- log2(out$nc_mean_linear + k * out$sd_typical)
  out$delta          <- out$lod_stabilised - out$lod_asis
  out
}

#' Leave-one-blank-out: how far can ONE well move a bay's LOD?
#' This is the direct measure of the n = 4 fragility.
lod_leave_one_out <- function(nc, markers, k = 3) {
  pb <- as.character(nc$plate_bay)
  do.call(rbind, lapply(unique(pb), function(p) {
    idx <- which(pb == p)
    if (length(idx) < 3) return(NULL)
    rng <- vapply(markers, function(t) {
      v <- 2^nc[[t]][idx]; v <- v[is.finite(v)]
      if (length(v) < 3) return(NA_real_)
      l <- vapply(seq_along(v), function(i) {
        vv <- v[-i]; log2(mean(vv) + k * stats::sd(vv)) }, numeric(1))
      max(l) - min(l)
    }, numeric(1))
    data.frame(plate_bay = p, target = markers, loo_lod_range = rng,
               stringsAsFactors = FALSE, row.names = NULL)
  }))
}

#' Per-bay blank background: how high do this bay's blanks read overall?
#' Flags a bay whose background sits > `z_flag` SD above the study norm — the
#' §2 recommendation. This is the quantity that a "low detectability" bay is
#' really reporting.
nc_background <- function(nc, markers, z_flag = 2) {
  pb <- as.character(nc$plate_bay)
  # centre each target across bays first, so the summary is not dominated by
  # markers that are simply abundant.
  M <- vapply(markers, function(t) nc[[t]], numeric(nrow(nc)))
  Z <- scale(M)
  lvl <- tapply(rowMeans(Z, na.rm = TRUE), pb, mean, na.rm = TRUE)
  z <- (lvl - mean(lvl, na.rm = TRUE)) / stats::sd(lvl, na.rm = TRUE)
  data.frame(plate_bay = names(lvl),
             nc_background = round(as.numeric(lvl), 4),
             nc_background_z = round(as.numeric(z), 2),
             flag_high_background = as.numeric(z) > z_flag,
             stringsAsFactors = FALSE, row.names = NULL)
}

#' Recompute per-sample detectability against a supplied LOD table.
#'
#' @param npq  matrix/data.frame, rows = targets, cols = samples (Alamar layout).
#' @param sample_bay  named vector mapping sample column name -> plate_bay.
#' @param lod  data frame with plate_bay, target and the LOD column to use.
#' @param lod_col  which LOD column.
detectability_from_lod <- function(npq, sample_bay, lod, lod_col = "lod_asis") {
  tg <- rownames(npq)
  key <- paste(lod$plate_bay, lod$target)
  lv  <- lod[[lod_col]]
  samples <- intersect(colnames(npq), names(sample_bay))
  det <- vapply(samples, function(s) {
    b <- sample_bay[[s]]
    thr <- lv[match(paste(b, tg), key)]
    v <- suppressWarnings(as.numeric(npq[, s]))
    ok <- is.finite(v) & is.finite(thr)
    if (!any(ok)) return(NA_real_)
    mean(v[ok] > thr[ok])
  }, numeric(1))
  data.frame(sample = samples, plate_bay = unname(sample_bay[samples]),
             detectability = round(det, 4), stringsAsFactors = FALSE, row.names = NULL)
}

#' Per-marker share of wells falling below their own bay's LOD, split by a
#' grouping (typically year or lot). This is the pTau-217 check: a marker's
#' detection CALLS can die batch-wide from a rising blank floor while its
#' quantitation is untouched.
marker_below_lod <- function(npq, sample_bay, lod, group, lod_col = "lod_asis") {
  tg <- rownames(npq)
  key <- paste(lod$plate_bay, lod$target); lv <- lod[[lod_col]]
  samples <- intersect(colnames(npq), names(sample_bay))
  g <- group[samples]
  res <- lapply(unique(g[!is.na(g)]), function(gg) {
    ss <- samples[which(g == gg)]
    below <- vapply(tg, function(t) {
      thr <- lv[match(paste(sample_bay[ss], t), key)]
      v <- suppressWarnings(as.numeric(npq[t, ss]))
      ok <- is.finite(v) & is.finite(thr)
      if (!any(ok)) return(NA_real_) else mean(v[ok] <= thr[ok])
    }, numeric(1))
    data.frame(group = gg, target = tg, pct_below_lod = round(100 * below, 2),
               stringsAsFactors = FALSE, row.names = NULL)
  })
  do.call(rbind, res)
}
