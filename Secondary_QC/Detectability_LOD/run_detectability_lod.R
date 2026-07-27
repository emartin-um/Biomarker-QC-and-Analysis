#!/usr/bin/env Rscript
# =============================================================================
# run_detectability_lod.R  ·  Secondary QC — the blank-well axis
# -----------------------------------------------------------------------------
# Standing detectability-by-bay panel + the variance-stabilised LOD.
# 2026-07-13 recommendations §2 and §8. BASE R ONLY.
#
# Run from the module directory:  Rscript run_detectability_lod.R
#
# Writes to output_files/:
#   nc_background_by_bay.csv    each bay's blank background, z vs the study norm
#   lod_comparison.csv          per bay x target: LOD as shipped vs stabilised
#   lod_leave_one_out.csv       how far ONE blank well can move a bay's LOD
#   detectability_by_bay.csv    detectability under both LODs, per bay
#   marker_below_lod_by_lot.csv per marker, % of wells under their own bay's LOD
# =============================================================================

HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
RUN_ROOT <- normalizePath(file.path(HERE, "..", ".."), mustWork = FALSE)
OUT_DIR  <- file.path(HERE, "output_files")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(HERE, "detectability_lod_helpers.R"))

REP <- file.path(RUN_ROOT, "Secondary_QC", "Replicate_Analysis")
NC_PATH  <- file.path(REP, "output_files", "nc_reps_NPQ.csv")
NPQ_PATH <- file.path(REP, "input_files", "NPQ_20260623.csv")
ANN_PATH <- file.path(REP, "input_files", "Annotation_Targets.csv")
SQC_PATH <- file.path(RUN_ROOT, "Primary_QC", "input_files", "Sample_QC.csv")

need <- c(NC_PATH, SQC_PATH)
if (!all(file.exists(need))) {
  cat("Missing required input(s):\n"); cat(paste0("  ", need[!file.exists(need)], "\n"))
  cat("\nThis module needs the Replicate_Analysis NC export and Primary_QC's Sample_QC.\n")
  quit(status = 0L)
}

## --- blanks ------------------------------------------------------------------
nc <- utils::read.csv(NC_PATH, stringsAsFactors = FALSE, check.names = FALSE)
nc$plate_bay <- paste0(nc$Run, "_", nc$Bay)
markers <- setdiff(names(nc), c("SampleID", "Run", "Bay", "plate_bay"))
markers <- markers[vapply(nc[markers], is.numeric, logical(1))]
cat(sprintf("Blank (NC) wells: %d across %d plate-bays; %d targets (%.1f blanks per bay)\n",
            nrow(nc), length(unique(nc$plate_bay)), length(markers),
            nrow(nc) / length(unique(nc$plate_bay))))

## --- LOD: as shipped vs variance-stabilised ---------------------------------
lod <- lod_table(nc, markers)
loo <- lod_leave_one_out(nc, markers)
bg  <- nc_background(nc, markers)

## --- optional cross-check against Alamar's own shipped LOD ------------------
if (file.exists(ANN_PATH)) {
  a <- utils::read.csv(ANN_PATH, stringsAsFactors = FALSE)
  a$plate_bay <- sub("_CNS.*$", "", a$plateID)
  m <- a$targetLOD_NPQ[match(paste(lod$plate_bay, lod$target),
                             paste(a$plate_bay, a$targetName))]
  d <- abs(lod$lod_asis - m)
  cat(sprintf("Cross-check vs Alamar's shipped targetLOD_NPQ: median |diff| %.6f over %d pairs\n",
              stats::median(d, na.rm = TRUE), sum(is.finite(d))))
}

## --- detectability under each LOD -------------------------------------------
det_bay <- NULL; below <- NULL
if (file.exists(NPQ_PATH)) {
  npqdf <- utils::read.csv(NPQ_PATH, stringsAsFactors = FALSE, check.names = FALSE)
  rn <- npqdf[[1]]; npq <- as.matrix(npqdf[, -1, drop = FALSE]); rownames(npq) <- rn

  q <- utils::read.csv(SQC_PATH, skip = 12, check.names = FALSE, stringsAsFactors = FALSE)
  q <- q[q[["Sample Type"]] == "Sample", ]
  sb <- sub("_CNS.*$", "", q[["plateID"]]); names(sb) <- q[["Sample Name"]]
  sb <- sb[intersect(names(sb), colnames(npq))]
  cat(sprintf("Patient wells mapped to a plate-bay: %d\n", length(sb)))

  d_asis <- detectability_from_lod(npq, sb, lod, "lod_asis")
  d_stab <- detectability_from_lod(npq, sb, lod, "lod_stabilised")
  d_asis$detectability_stabilised <- d_stab$detectability[match(d_asis$sample, d_stab$sample)]

  det_bay <- do.call(rbind, lapply(split(d_asis, d_asis$plate_bay), function(x)
    data.frame(plate_bay = x$plate_bay[1], n_wells = nrow(x),
               detectability_asis = round(mean(x$detectability, na.rm = TRUE), 4),
               detectability_stabilised = round(mean(x$detectability_stabilised, na.rm = TRUE), 4),
               pct_below_090_asis = round(100 * mean(x$detectability < 0.90, na.rm = TRUE), 1),
               pct_below_090_stabilised = round(100 * mean(x$detectability_stabilised < 0.90, na.rm = TRUE), 1),
               stringsAsFactors = FALSE, row.names = NULL)))
  det_bay <- merge(det_bay, bg, by = "plate_bay", all.x = TRUE)
  det_bay <- det_bay[order(det_bay$detectability_asis), ]

  lot <- ifelse(substr(sb, 1, 4) == "2026", "2026", "2025"); names(lot) <- names(sb)
  below <- marker_below_lod(npq, sb, lod, lot)
}

## --- write -------------------------------------------------------------------
w <- function(x, f) if (!is.null(x) && nrow(x)) {
  utils::write.csv(x, file.path(OUT_DIR, f), row.names = FALSE)
  cat(sprintf("  %-32s %6d rows\n", f, nrow(x))) }
cat("\nWrote:\n")
w(bg,      "nc_background_by_bay.csv")
w(lod,     "lod_comparison.csv")
w(loo,     "lod_leave_one_out.csv")
w(det_bay, "detectability_by_bay.csv")
w(below,   "marker_below_lod_by_lot.csv")

## --- summary -----------------------------------------------------------------
cat("\n================ SUMMARY ================\n")
cat(sprintf("\n-- §8: how fragile is a 4-blank LOD? --\n"))
cat(sprintf("  median leave-one-blank-out swing in LOD : %.3f log2\n",
            stats::median(loo$loo_lod_range, na.rm = TRUE)))
cat(sprintf("  90th percentile                          : %.3f log2\n",
            stats::quantile(loo$loo_lod_range, .90, na.rm = TRUE)))
cat(sprintf("  worst target x bay                       : %.3f log2\n",
            max(loo$loo_lod_range, na.rm = TRUE)))
cat(sprintf("  -> one blank well can move the threshold by this much, before any\n"))
cat(sprintf("     sample is considered. Stabilising keeps the per-bay MEAN and\n"))
cat(sprintf("     replaces the 4-point SD with the target's typical SD.\n"))

cat(sprintf("\n-- §2: bays whose BLANKS read high (the real 'low detectability' signal) --\n"))
hi <- bg[bg$flag_high_background %in% TRUE, ]
if (nrow(hi)) {
  for (i in order(-hi$nc_background_z))
    cat(sprintf("  %-26s background z %+5.2f\n", hi$plate_bay[i], hi$nc_background_z[i]))
} else cat("  none above +2 SD of the study norm\n")

if (!is.null(det_bay)) {
  cat(sprintf("\n-- detectability under each LOD (5 lowest bays) --\n"))
  cat(sprintf("  %-26s %8s %8s   %6s %6s  %s\n", "plate-bay", "as-is", "stabil.",
              "<0.90", "<0.90", "bg z"))
  for (i in seq_len(min(5, nrow(det_bay))))
    cat(sprintf("  %-26s %8.3f %8.3f   %5.1f%% %5.1f%%  %+5.2f\n",
                det_bay$plate_bay[i], det_bay$detectability_asis[i],
                det_bay$detectability_stabilised[i], det_bay$pct_below_090_asis[i],
                det_bay$pct_below_090_stabilised[i], det_bay$nc_background_z[i]))
  cat(sprintf("  correlation, bay detectability vs blank background: %.3f\n",
              stats::cor(det_bay$detectability_asis, det_bay$nc_background,
                         use = "complete.obs")))
}

if (!is.null(below)) {
  cat(sprintf("\n-- markers whose DETECTION CALLS move most between lots --\n"))
  wide <- reshape(below, idvar = "target", timevar = "group", direction = "wide")
  cn <- grep("^pct_below_lod", names(wide), value = TRUE)
  if (length(cn) == 2) {
    wide$shift <- wide[[cn[2]]] - wide[[cn[1]]]
    wide <- wide[order(-abs(wide$shift)), ]
    cat(sprintf("  %-16s %8s %8s %8s\n", "target", sub("pct_below_lod.", "", cn[1]),
                sub("pct_below_lod.", "", cn[2]), "shift"))
    for (i in seq_len(min(6, nrow(wide))))
      cat(sprintf("  %-16s %7.1f%% %7.1f%% %+7.1f\n", wide$target[i],
                  wide[[cn[1]]][i], wide[[cn[2]]][i], wide$shift[i]))
  }
}

cat("\n🚫 NEVER GATE OR DROP ON DETECTABILITY. A low-detectability bay means its\n")
cat("   BLANKS read high, not that its samples failed — an assay-background and\n")
cat("   analysis matter, not a wet-lab re-run. Use continuous NPQ, and never gate\n")
cat("   a marker on per-bay detection status without checking its blank trend.\n")
