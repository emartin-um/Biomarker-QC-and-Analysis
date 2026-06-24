# =============================================================================
# extreme_helpers.R
# -----------------------------------------------------------------------------
# Core, I/O-free functions + disease-signature panel definitions for the
# Extreme-Value Screen (Secondary QC). Sourced by run_extreme_screen.R and by
# Extreme_Value_Screen.Rmd.
#
# Purpose: identify individual samples whose biomarker profile is *super
# extreme* on one marker, or coordinated across a disease signature, so a
# clinician may want to look at the participant more closely (possible ALS,
# Parkinson/DLB, early-onset / aggressive AD, an acute inflammatory illness,
# etc.). Values are log2 NPQ (Alamar NULISAseq CNS panel), normalised per-batch,
# so a single pooled reference distribution across the combined cohort is used.
#
# --- Why rank-based INT, not raw robust-z, is the screening metric -----------
# A plain robust z = (x - median)/mad over-flags markers with a *tight core and
# heavy tails* or *genotype-driven multimodality* (e.g. APOE4 mad/sd = 0.13 ->
# 35% of samples "|z|>5"; PTN 0.16 -> 26%). That is a distribution-shape
# artefact, not biology. We therefore screen on the rank-based inverse-normal
# transform (INT):
#       INT = qnorm( (rank(x) - 0.5) / n_nonNA )
# which maps every marker to a standard-normal *rank* scale, so "extreme" means
# the same empirical rarity for every marker (|INT|>=2.75 ~ top/bottom 0.3%).
# Magnitude is preserved for the reader by reporting robust_z, the raw value and
# the percentile alongside the INT. Both tails are screened - an unusually LOW
# value (e.g. Abeta42 with amyloid deposition) can also be informative.
#
# Panels score the *coordinated* signal: mean oriented INT across the signature
# markers (orientation = +1 if HIGH is the disease direction, -1 if LOW), and
# require >= panel_min_hit markers individually in the tail, so one flier cannot
# manufacture a signature.
#
# Sample key: SAMPLE (accession; unique per row). ID2 is NA for the Other
# cohort, so it is carried only as context, never as a join key.
# Nothing here reads or writes files; callers supply data frames.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a   # null-coalesce (base R lacks it pre-4.4)

# -----------------------------------------------------------------------------
# 1. Marker housekeeping.
# -----------------------------------------------------------------------------

# Hemolysis-sensitive (12) - extremes here are often an RBC-leakage artefact
# rather than biology. Source: biomarker notes Hemolysis_Sensitive == TRUE.
HEMOLYSIS_MARKERS <- c("HBA1","PRDX6","MDH1","ARSA","SOD1","SNCA",
                       "PGK1","IL16","S100A12","ENO2","S100B","CXCL8")

# Genotype-/assay-determined - removed from the screen entirely. APOE4 is also
# confounded by the Bay1 assay anomaly; both APOE4 and total APOE are dominated
# by APOE genotype, so an "extreme" there reflects e4 dosage, not illness.
EXCLUDE_MARKERS <- c("APOE4","APOE")

# Known functional-variant / bimodal markers - kept (a high value can be real)
# but tagged 'caution' in the reference table so reviewers discount a lone hit.
# CHIT1: 24-bp dup loss-of-function; ANXA5/ARSA/FCN2/PARK7: notes list variants.
GENETIC_CAUTION <- c("CHIT1","ANXA5","ARSA","FCN2","PARK7")

# Run-bay with the known APOE4 assay anomaly (APOE4 channel scrambled).
APOE4_BAY1_RUN <- "20251124-1407_Bay1"

# -----------------------------------------------------------------------------
# 2. Disease-signature panels.
#    Each entry: named numeric vector  marker -> orientation (+1 HIGH = disease,
#    -1 LOW = disease).  Only markers present in the data are used.
#    'min_markers' = non-missing panel markers needed to compute a score.
# -----------------------------------------------------------------------------
PANELS <- list(

  # ALS / FTLD - TDP-43 proteinopathy. Kept TDP-43-SPECIFIC: pTDP43_409 + TARDBP
  # are the defining proteinopathy. Neurofilament (NEFL/NEFH) is deliberately
  # NOT here - it is non-specific (high in advanced AD, and partly batch-driven
  # on the reads-high 20251125-1332 plate) and lives in the Neurodegeneration
  # panel instead. A true ALS case lights up BOTH (TDP-43 + neurofilament) ->
  # mixed-axis -> HIGH. SOD1 excluded (familial-ALS gene, hemolysis-confounded).
  # NOTE: pTDP-43/TARDBP are intracellular, so a hemolysed sample (high HBA1_INT)
  # can elevate them non-specifically - the hemolysis caveat applies.
  # score_thr override: pTDP43_409 & TARDBP are highly correlated (r~0.87), so a
  # 2-marker mean over-fires at the default cut; 2.5 keeps the genuinely extreme.
  # solo: a SAFETY NET so a true ALS case cannot slip through - very high
  # neurofilament (NEFL, the ALS hallmark) fires this panel on its own even
  # without TDP-43. NEFL is non-specific (also high in advanced AD), so the bar
  # is the extreme tail (oriented INT >= 3.0, ~top 0.1%) and it is labelled as
  # non-specific. NEFL otherwise lives in the Neurodegeneration panel.
  ALS_MND = list(
    markers = c(pTDP43_409 = 1, TARDBP = 1),
    min_markers = 2, score_thr = 2.5,
    solo = c(NEFL = 3.0),
    solo_note = "Possible ALS / aggressive neurodegeneration - very high neurofilament (non-specific)",
    label = "ALS / FTLD (TDP-43 proteinopathy): high pTDP-43 + TDP-43"
  ),

  # Synucleinopathy: Parkinson's disease / dementia with Lewy bodies.
  # Total SNCA is deliberately omitted - plasma a-synuclein is largely RBC-
  # derived (hemolysis-sensitive), so it is reported separately, not scored.
  # Oligo_SNCA / pSNCA_129 are the pathological species.
  Synucleinopathy = list(
    markers = c(Oligo_SNCA = 1, pSNCA_129 = 1, DDC = 1, PARK7 = 1),
    min_markers = 3,
    label = "Synucleinopathy (PD / DLB): high oligomeric / pSer129 a-synuclein"
  ),

  # AD tau pathology - extreme tau, esp. notable when the participant is young
  # (early-onset) or clinically unimpaired (pre-symptomatic / discordant).
  AD_tau = list(
    markers = c(pTau_217 = 1, pTau_181 = 1, pTau_231 = 1, BD_pTau_217 = 1, MAPT = 1),
    min_markers = 3,
    label = "AD tau pathology: high phospho-tau (flag if young / NCI)"
  ),

  # Low amyloid - a low A42/A40 ratio is the classic brain-amyloid-deposition
  # signal (A42 falls faster than A40). The derived A42_A40_ratio (log2) is
  # added by the driver; low A42 itself is kept as a corroborating marker.
  Amyloid_low = list(
    markers = c(A42_A40_ratio = -1, A_42 = -1),
    min_markers = 1,
    label = "Low amyloid (low A42/A40 ratio): possible deposition"
  ),

  # Brain-derived (exosome) pTau rare inflation - specifically requested.
  BD_pTau_inflation = list(
    markers = c(BD_pTau_181 = 1, BD_pTau_217 = 1, BD_pTau_231 = 1, BD_MAPT = 1),
    min_markers = 2,
    label = "Rare BD-pTau (brain-derived/exosome) inflation"
  ),

  # Neuroinflammation 'storm' - many acute-phase / pro-inflammatory markers up
  # together can flag an intercurrent infection / systemic illness at draw.
  Neuroinflammation = list(
    markers = c(CRP = 1, IL6 = 1, TNF = 1, IL1B = 1, SAA1 = 1,
                CXCL10 = 1, GDF15 = 1, IL18 = 1, CHI3L1 = 1),
    min_markers = 5,
    label = "Acute inflammatory state: many cytokines / acute-phase high"
  ),

  # Global neurodegeneration axis (astroglial + axonal + mitochondrial stress).
  Neurodegeneration = list(
    markers = c(GFAP = 1, NEFL = 1, NEFH = 1, GDF15 = 1),
    min_markers = 3,
    label = "Global neurodegeneration: high GFAP + neurofilament"
  )
)

# Panels that represent a neurodegenerative pathology (used by the shortlist's
# age/dx-discordance amplifier; Neuroinflammation excluded there).
PATHOLOGY_PANELS <- c("ALS_MND","Synucleinopathy","AD_tau",
                      "BD_pTau_inflation","Amyloid_low","Neurodegeneration")

# Distinct pathology AXES - overlapping panels collapse to one axis so that the
# same biology is not double-counted as "two pathologies". Co-occurrence of >=2
# *distinct* axes is the clinically interesting mixed-pathology signal (e.g.
# AD tau + Lewy synuclein, or AD + TDP-43). Neurodegeneration (NEFL/NEFH/GFAP)
# is non-specific and overlaps MND, so it is not its own axis.
PATHOLOGY_AXES <- list(
  TAU       = c("AD_tau", "BD_pTau_inflation"),
  AMYLOID   = c("Amyloid_low"),
  SYNUCLEIN = c("Synucleinopathy"),
  MND       = c("ALS_MND")   # TDP-43 proteinopathy or very-high neurofilament
)

# -----------------------------------------------------------------------------
# 3. Per-marker reference stats and the three aligned matrices
#    (INT for screening; robust-z and raw value for magnitude/reporting).
# -----------------------------------------------------------------------------

#' Markers to screen: numeric, varying, not genotype-excluded.
usable_markers <- function(df, biomarker_cols) {
  keep <- setdiff(intersect(biomarker_cols, names(df)), EXCLUDE_MARKERS)
  keep[vapply(keep, function(b) {
    v <- suppressWarnings(as.numeric(df[[b]]))
    sum(!is.na(v)) >= 20 && stats::mad(v, na.rm = TRUE) > 0
  }, logical(1))]
}

#' Reference table: median/MAD/SD/quantiles + caution flags per marker.
reference_stats <- function(df, markers) {
  do.call(rbind, lapply(markers, function(b) {
    v   <- suppressWarnings(as.numeric(df[[b]]))
    m   <- stats::mad(v, na.rm = TRUE)
    s   <- stats::sd(v, na.rm = TRUE)
    rat <- ifelse(s > 0, m / s, NA_real_)
    data.frame(
      Biomarker     = b,
      n             = sum(!is.na(v)),
      median        = round(stats::median(v, na.rm = TRUE), 3),
      mad           = round(m, 3),
      sd            = round(s, 3),
      mad_sd_ratio  = round(rat, 2),                 # <0.4 => tight core / multimodal
      q005          = round(unname(stats::quantile(v, 0.005, na.rm = TRUE)), 3),
      q995          = round(unname(stats::quantile(v, 0.995, na.rm = TRUE)), 3),
      min           = round(min(v, na.rm = TRUE), 3),
      max           = round(max(v, na.rm = TRUE), 3),
      hemolysis_sensitive = b %in% HEMOLYSIS_MARKERS,
      genetic_caution     = b %in% GENETIC_CAUTION,
      distribution_caution = !is.na(rat) && rat < 0.4,
      row.names = NULL, stringsAsFactors = FALSE
    )
  }))
}

#' Rank-based inverse-normal scores for one vector (NA-safe).
int_scores <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  ok <- !is.na(v)
  out <- rep(NA_real_, length(v))
  n <- sum(ok)
  if (n >= 2) out[ok] <- stats::qnorm((rank(v[ok], ties.method = "average") - 0.5) / n)
  out
}

#' Matrix builders (samples x markers), all aligned by column order = markers.
int_matrix    <- function(df, markers) { m <- vapply(markers, function(b) int_scores(df[[b]]), numeric(nrow(df))); colnames(m) <- markers; m }
robust_z_matrix <- function(df, ref) {
  markers <- ref$Biomarker
  z <- vapply(seq_along(markers), function(i)
    (suppressWarnings(as.numeric(df[[markers[i]]])) - ref$median[i]) / ref$mad[i], numeric(nrow(df)))
  colnames(z) <- markers; z
}
value_matrix  <- function(df, markers) { m <- vapply(markers, function(b) suppressWarnings(as.numeric(df[[b]])), numeric(nrow(df))); colnames(m) <- markers; m }
percentile_matrix <- function(df, markers) {
  p <- vapply(markers, function(b) {
    v <- suppressWarnings(as.numeric(df[[b]]))
    100 * (rank(v, na.last = "keep", ties.method = "average") - 0.5) / sum(!is.na(v))
  }, numeric(nrow(df)))
  colnames(p) <- markers; p
}

# -----------------------------------------------------------------------------
# 4. Long table of per-(sample, marker) extreme events (gated on |INT|).
# -----------------------------------------------------------------------------

#' @param int,zmat,pct,values  aligned matrices (samples x markers)
#' @param ids   SAMPLE ids (length nrow); @param thr_int  |INT| gate
build_long_extremes <- function(int, zmat, pct, values, ids, thr_int = 2.75) {
  idx <- which(abs(int) >= thr_int, arr.ind = TRUE)
  empty <- data.frame(SAMPLE = character(), Biomarker = character(), INT = numeric(),
                      robust_z = numeric(), direction = character(), value = numeric(),
                      percentile = numeric(), hemolysis_sensitive = logical(),
                      stringsAsFactors = FALSE)
  if (nrow(idx) == 0) return(empty)
  markers <- colnames(int)
  out <- data.frame(
    SAMPLE     = ids[idx[, "row"]],
    Biomarker  = markers[idx[, "col"]],
    INT        = round(int[idx], 2),
    robust_z   = round(zmat[idx], 2),
    direction  = ifelse(int[idx] > 0, "HIGH", "LOW"),
    value      = round(values[idx], 3),
    percentile = round(pct[idx], 2),
    stringsAsFactors = FALSE
  )
  out$hemolysis_sensitive <- out$Biomarker %in% HEMOLYSIS_MARKERS
  out[order(-abs(out$INT), -abs(out$robust_z)), ]
}

# -----------------------------------------------------------------------------
# 5. Panel scores (coordinated multi-marker signal), on the INT scale.
# -----------------------------------------------------------------------------

#' '<panel>_score' = mean oriented INT across present markers (>= min_markers
#' non-missing); '<panel>_nhit' = count of panel markers with oriented INT >=
#' panel_hit (i.e. individually in the disease-direction tail).
score_panels <- function(int, panels = PANELS, panel_hit = 1.8) {
  present <- colnames(int)
  res <- list()
  for (nm in names(panels)) {
    spec <- panels[[nm]]
    mk   <- intersect(names(spec$markers), present)
    if (length(mk) == 0) next
    io   <- sweep(int[, mk, drop = FALSE], 2, spec$markers[mk], `*`)   # oriented INT
    n_ok <- rowSums(!is.na(io))
    res[[paste0(nm, "_score")]] <- round(ifelse(n_ok >= spec$min_markers, rowMeans(io, na.rm = TRUE), NA_real_), 2)
    res[[paste0(nm, "_nhit")]]  <- rowSums(io >= panel_hit, na.rm = TRUE)
  }
  as.data.frame(res, stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# 6. Per-sample master row: extreme counts, top markers, panel scores, flags.
# -----------------------------------------------------------------------------

#' @param int  INT matrix; @param meta context df aligned row-wise (needs SAMPLE)
build_master <- function(int, meta, panel_scores, thr_int = 2.75) {
  markers <- colnames(int)
  n_hi <- rowSums(int >=  thr_int, na.rm = TRUE)
  n_lo <- rowSums(int <= -thr_int, na.rm = TRUE)

  # Top-3 extreme markers (by |INT|) as a readable string e.g. "NEFL(+3.4); ...".
  top_str <- vapply(seq_len(nrow(int)), function(i) {
    zi <- int[i, ]; zi[is.na(zi)] <- 0
    o  <- order(abs(zi), decreasing = TRUE)
    o  <- o[abs(zi[o]) >= thr_int]
    if (length(o) == 0) return("")
    paste(sprintf("%s(%+.1f)", markers[head(o, 3)], zi[head(o, 3)]), collapse = "; ")
  }, character(1))

  hemo_cols <- intersect(markers, HEMOLYSIS_MARKERS)
  hemo_hi   <- rowSums(int[, hemo_cols, drop = FALSE] >= thr_int, na.rm = TRUE)
  frac_hemo <- ifelse(n_hi > 0, hemo_hi / n_hi, 0)
  max_abs   <- apply(abs(int), 1, max, na.rm = TRUE); max_abs[!is.finite(max_abs)] <- 0

  master <- data.frame(
    SAMPLE          = meta$SAMPLE,
    n_extreme_high  = n_hi,
    n_extreme_low   = n_lo,
    n_extreme_total = n_hi + n_lo,
    max_abs_INT     = round(max_abs, 2),
    top_extreme     = top_str,
    HBA1_INT        = if ("HBA1" %in% markers) round(int[, "HBA1"], 2) else NA_real_,
    frac_hemo_high  = round(frac_hemo, 2),
    stringsAsFactors = FALSE
  )
  master <- cbind(master, panel_scores)

  # solo-trigger marker INTs (so annotate_master can fire a panel on one marker)
  for (sm in setdiff(solo_markers(), "HBA1"))
    master[[paste0(sm, "_INT")]] <- if (sm %in% markers) round(int[, sm], 2) else NA_real_

  ctx <- intersect(c("dataset","metadata_source","CDX_collapsed","age_at_subject",
                     "sex","RUN","Ancestry","APOE.geno_final","ID2","SAMPLE_ALIQUOT",
                     "Record_ID","Site","Years_Onset"), names(meta))
  master <- cbind(master, meta[, ctx, drop = FALSE])

  master$apoe4_bay1_flag <- if ("RUN" %in% names(meta)) meta$RUN == APOE4_BAY1_RUN else FALSE
  master$age_young        <- !is.na(master$age_at_subject) & master$age_at_subject < 60
  master$global_outlier   <- master$n_extreme_total >= 8    # very sick OR bad sample
  master
}

# -----------------------------------------------------------------------------
# 7. Clinician shortlist - the most reviewable cases + plain-language reason.
# -----------------------------------------------------------------------------

#' Per-panel score threshold (a panel may override the global default).
panel_thr <- function(pn, default, panels = PANELS) {
  st <- panels[[pn]]$score_thr
  if (is.null(st)) default else st
}

#' Markers used as panel "solo" triggers (their oriented INT is needed in the
#' master so annotate_master can fire a panel on a single very-extreme marker).
solo_markers <- function(panels = PANELS) {
  unique(unlist(lapply(panels, function(p) names(p$solo))))
}

#' marker -> disease orientation (+1 HIGH is disease, -1 LOW is disease) across
#' all panels (used to make the lone-marker rule direction-aware).
disease_direction_map <- function(panels = PANELS) {
  m <- list()
  for (p in panels) for (mk in names(p$markers)) if (is.null(m[[mk]])) m[[mk]] <- unname(p$markers[mk])
  unlist(m)
}

#' Qualifies a sample on:
#'   (a) a coordinated panel signal: mean oriented INT >= panel_score_thr AND
#'       >= panel_min_hit panel markers in the disease-direction tail; or
#'   (b) a very rare lone disease marker IN THE DISEASE DIRECTION (oriented INT
#'       >= lone_int_thr), not hemolysis-driven.
#' Age/dx discordance is layered on as an amplifier. Output is ranked by a
#' magnitude-aware priority (NOT the INT cap), with an explicit priority tier.
#' Annotate EVERY sample with review_reason / priority / pathology axes / caveats
#' (no filtering) - the report uses this full table; build_shortlist() filters it.
annotate_master <- function(master, long_ext, panels = PANELS, panel_score_thr = 2.0,
                            panel_min_hit = 2, panel_strong_hit = 3, lone_int_thr = 3.0,
                            present_markers = NULL) {

  panel_names <- names(panels)
  ddir        <- disease_direction_map(panels)       # marker -> +/-1
  # Effective min-hit per panel cannot exceed how many of its markers are present
  # (so a 1-2 marker panel like Amyloid can still trigger).
  panel_emin <- vapply(panels, function(p) {
    np <- if (is.null(present_markers)) length(p$markers) else length(intersect(names(p$markers), present_markers))
    min(panel_min_hit, max(np, 1))
  }, numeric(1))

  reason <- character(nrow(master)); n_axes <- integer(nrow(master))
  best_panel <- numeric(nrow(master)); lone_strength <- numeric(nrow(master))
  max_nhit <- integer(nrow(master)); axes_str <- character(nrow(master))

  for (i in seq_len(nrow(master))) {
    r <- character(0); trig_set <- character(0); bp <- 0; ls <- 0; mnh <- 0

    # (a) coordinated panel signatures. Trigger if the whole signature is
    #     elevated (mean score AND >=min_hit in tail) OR enough individual
    #     markers spike (>=strong_hit) - the latter catches broad 'storm'
    #     panels (e.g. inflammation) where a subset spikes but the mean dilutes.
    for (pn in panel_names) {
      sc <- master[[paste0(pn, "_score")]][i]; nh <- master[[paste0(pn, "_nhit")]][i]
      thr_pn     <- panel_thr(pn, panel_score_thr, panels)
      hit_mean   <- !is.na(sc) && sc >= thr_pn && nh >= panel_emin[[pn]]
      hit_strong <- !is.na(nh) && nh >= panel_strong_hit
      if (hit_mean || hit_strong) {
        sc_disp <- ifelse(is.na(sc) || abs(sc) < 0.05, 0, sc)
        r <- c(r, sprintf("%s [%d markers in tail, mean INT=%.1f]", panels[[pn]]$label, nh, sc_disp))
        trig_set <- c(trig_set, pn)
        bp  <- max(bp, sc_disp, ifelse(hit_strong, panel_score_thr, 0))  # strong-hit counts as >=thr
        mnh <- max(mnh, ifelse(is.na(nh), 0, nh))
      }

      # solo trigger: one very-extreme marker (in the disease direction) fires the
      # panel on its own - a safety net so a hallmark single marker is not missed.
      if (!is.null(panels[[pn]]$solo)) {
        for (sm in names(panels[[pn]]$solo)) {
          col <- paste0(sm, "_INT")
          val <- if (col %in% names(master)) master[[col]][i] else NA_real_
          ori <- if (!is.na(val) && sm %in% names(ddir)) val * ddir[[sm]] else NA_real_
          if (!is.na(ori) && ori >= panels[[pn]]$solo[[sm]]) {
            note <- panels[[pn]]$solo_note %||% sprintf("%s (single-marker)", panels[[pn]]$label)
            r <- c(r, sprintf("%s [%s INT=%+.1f]", note, sm, val))
            if (!pn %in% trig_set) trig_set <- c(trig_set, pn)
            bp <- max(bp, panel_thr(pn, panel_score_thr, panels))
          }
        }
      }
    }

    # Distinct pathology axes triggered (collapses overlapping tau panels etc.).
    axes_hit <- names(PATHOLOGY_AXES)[vapply(PATHOLOGY_AXES,
                  function(p) any(p %in% trig_set), logical(1))]
    has_path <- length(intersect(trig_set, PATHOLOGY_PANELS)) > 0

    # (b) lone disease marker, disease-direction only, not hemolysis-sensitive
    le <- long_ext[long_ext$SAMPLE == master$SAMPLE[i] &
                   long_ext$Biomarker %in% names(ddir) &
                   !long_ext$hemolysis_sensitive, , drop = FALSE]
    if (nrow(le) > 0) {
      le$oriented_int <- le$INT * ddir[le$Biomarker]
      le <- le[le$oriented_int >= lone_int_thr, , drop = FALSE]
      if (nrow(le) > 0) {
        le <- le[order(-le$oriented_int, -abs(le$robust_z)), ][1, ]
        r  <- c(r, sprintf("Extreme %s %s (INT=%+.1f, robust_z=%+.1f, %.2f pctile)",
                           le$Biomarker, le$direction, le$INT, le$robust_z, le$percentile))
        ls <- abs(le$robust_z)
      }
    }

    if (length(r) == 0) { reason[i] <- NA_character_; next }

    # (b2) mixed pathology across >=2 distinct axes - clinically notable.
    if (length(axes_hit) >= 2)
      r <- c(r, sprintf("MIXED pathology: %s axes co-elevated", paste(axes_hit, collapse = "+")))

    # (c) discordance amplifier (only when a pathology signature is present)
    if (has_path) {
      if (isTRUE(master$age_young[i]))
        r <- c(r, sprintf("YOUNG (age %.0f) with pathology signature - possible early-onset",
                          master$age_at_subject[i]))
      if (!is.na(master$CDX_collapsed[i]) && master$CDX_collapsed[i] == "NCI")
        r <- c(r, "Clinically NCI but carries a pathology signature - discordant")
    }
    reason[i] <- paste(r, collapse = " | "); n_axes[i] <- length(axes_hit)
    best_panel[i] <- bp; lone_strength[i] <- ls; max_nhit[i] <- mnh
    axes_str[i] <- paste(axes_hit, collapse = "+")
  }

  master$review_reason    <- reason
  master$n_pathology_axes <- n_axes
  master$pathology_axes   <- axes_str
  master$best_panel_score <- round(best_panel, 2)
  master$max_panel_nhit   <- max_nhit
  master$lone_robustz     <- round(lone_strength, 2)

  # Priority tier - HIGH = clinically *special* (discordant, mixed, or super-
  # extreme), not merely the expected disease signal (e.g. AD with high tau).
  discordant <- (master$age_young | (!is.na(master$CDX_collapsed) & master$CDX_collapsed == "NCI"))
  has_path_v <- master$n_pathology_axes >= 1
  master$priority <- with(master, ifelse(
    n_pathology_axes >= 2 |           # mixed pathology
    best_panel_score >= 2.8 |          # a single, very strong signature
    lone_robustz >= 6 |                # a striking lone marker
    max_panel_nhit >= 5 |              # inflammation storm / many markers
    (discordant & has_path_v),         # NCI / young carrying any pathology
    "HIGH", "MEDIUM"))

  # Caveats (informative, not disqualifying) ---------------------------------
  cav <- rep("", nrow(master))
  cav <- ifelse((master$frac_hemo_high >= 0.5 & master$n_extreme_high >= 3) |
                  (!is.na(master$HBA1_INT) & master$HBA1_INT >= 2.75),
                "CAVEAT: high hemolysis (HBA1/RBC markers) - discount leakage-sensitive hits", cav)
  cav <- ifelse(master$global_outlier,
                trimws(paste(cav, "| CAVEAT: global outlier (>=8 markers) - check sample quality")), cav)
  cav <- ifelse(master$apoe4_bay1_flag,
                trimws(paste(cav, "| NOTE: Bay1 APOE4-anomaly run")), cav)
  master$review_caveat <- trimws(gsub("^\\| ", "", cav))
  master
}

#' Filter the annotated master to flagged samples, ranked HIGH-first.
build_shortlist <- function(master, long_ext, ...) {
  m <- if ("review_reason" %in% names(master)) master else annotate_master(master, long_ext, ...)
  short <- m[!is.na(m$review_reason), , drop = FALSE]
  short[order(factor(short$priority, levels = c("HIGH","MEDIUM")),
              -short$n_pathology_axes, -short$best_panel_score, -short$lone_robustz), ]
}

# -----------------------------------------------------------------------------
# 8. One-call pipeline - the single compute path shared by the headless driver
#    and the Rmd report. Takes a pooled data frame (must have SAMPLE + the
#    biomarker columns + context cols) and returns every artefact.
# -----------------------------------------------------------------------------
run_extreme_pipeline <- function(dat, biomarker_cols,
                                 thr_int = 2.75, panel_hit = 1.8, panel_score_thr = 2.0,
                                 panel_min_hit = 2, lone_int_thr = 3.0) {
  if (all(c("A_42","A_40") %in% names(dat))) dat$A42_A40_ratio <- dat$A_42 - dat$A_40
  mk <- usable_markers(dat, biomarker_cols)
  if ("A42_A40_ratio" %in% names(dat)) mk <- c(mk, "A42_A40_ratio")

  ref  <- reference_stats(dat, mk)
  imat <- int_matrix(dat, mk)
  zmat <- robust_z_matrix(dat, ref)
  vmat <- value_matrix(dat, mk)
  pmat <- percentile_matrix(dat, mk)

  long_ext <- build_long_extremes(imat, zmat, pmat, vmat, ids = dat$SAMPLE, thr_int = thr_int)
  psc      <- score_panels(imat, PANELS, panel_hit = panel_hit)
  master0  <- build_master(imat, dat, psc, thr_int = thr_int)
  master   <- annotate_master(master0, long_ext, PANELS, panel_score_thr = panel_score_thr,
                              panel_min_hit = panel_min_hit, lone_int_thr = lone_int_thr,
                              present_markers = mk)
  short    <- build_shortlist(master, long_ext)   # already annotated -> just filter/sort

  list(markers = mk, ref = ref, int = imat, z = zmat, long_ext = long_ext,
       panel_scores = psc, master = master, short = short, dat = dat,
       params = list(thr_int = thr_int, panel_hit = panel_hit, panel_score_thr = panel_score_thr,
                     panel_min_hit = panel_min_hit, lone_int_thr = lone_int_thr))
}
