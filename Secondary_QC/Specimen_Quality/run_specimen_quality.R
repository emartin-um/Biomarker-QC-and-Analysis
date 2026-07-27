#!/usr/bin/env Rscript
# =============================================================================
# run_specimen_quality.R  ·  Secondary QC — the per-specimen quality axis
# -----------------------------------------------------------------------------
# Standing panels for the three per-specimen properties the pipeline did not
# report: sequencing depth, cystatin C (renal), and whole-panel level.
#
# BASE R ONLY — no dplyr — so it runs anywhere the QC runs do.
# Run from the module directory:  Rscript run_specimen_quality.R
#
# Writes to output_files/:
#   specimen_axis_by_site.csv        per-site distribution of each axis
#   specimen_axis_by_ancestry.csv    same, by ancestry
#   specimen_axis_variance.csv       % of each axis between Site / Run / Bay / …
#   specimen_axis_cocluster.csv      do the axes single out the same sites?
#   specimen_axis_within_site.csv    does each axis survive within-site?
#   sign_balance_report.csv          exposure of every composite index we ship
#
# Implements: 2026-07-13 recommendations §2 (depth-by-site — never built until
# now) and the 2026-07-27 addendum §12 (cystatin C) and §6 (sign-balance guard).
# =============================================================================

HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
RUN_ROOT <- normalizePath(file.path(HERE, "..", ".."), mustWork = FALSE)
DATA_DIR <- file.path(RUN_ROOT, "datasets")
OUT_DIR  <- file.path(HERE, "output_files")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(HERE, "specimen_quality_helpers.R"))

## --- load the pooled cohort, exactly as the Extremes screen does -------------
u19o <- readRDS(file.path(DATA_DIR, "preprocessed_data_U19.rds"))
otho <- readRDS(file.path(DATA_DIR, "preprocessed_data_OtherProjects.rds"))
u19 <- as.data.frame(u19o$data); u19$dataset <- "U19"
oth <- as.data.frame(otho$data); oth$dataset <- "OtherProjects"
common <- intersect(names(u19), names(oth))
d <- rbind(u19[, common], oth[, common])
cat(sprintf("Pooled cohort: %d samples (U19 %d + Other %d)\n", nrow(d), nrow(u19), nrow(oth)))

## --- read depth: join Sample_QC.csv (the same join R/load_50plate.R uses) ----
SAMPLE_QC <- file.path(RUN_ROOT, "Primary_QC", "input_files", "Sample_QC.csv")
d$l2reads <- NA_real_
if (file.exists(SAMPLE_QC) && "SAMPLE_ALIQUOT" %in% names(d)) {
  q <- tryCatch(utils::read.csv(SAMPLE_QC, skip = 12, check.names = FALSE,
                                stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(q) && all(c("Sample Type", "Sample Name", "Reads") %in% names(q))) {
    q <- q[q[["Sample Type"]] == "Sample", ]
    reads <- .num(q[["Reads"]])[match(d$SAMPLE_ALIQUOT, q[["Sample Name"]])]
    d$reads <- reads; d$l2reads <- log2(reads)
    cat(sprintf("Read depth joined for %d of %d wells (median %.2fM reads)\n",
                sum(!is.na(reads)), nrow(d), stats::median(reads, na.rm = TRUE) / 1e6))
  }
}
if (all(is.na(d$l2reads)))
  cat("NOTE: no read depth available (Sample_QC.csv absent or unexpected layout) —\n",
      "      the depth panels will be skipped.\n")

## --- whole-panel level: reuse the Extremes column if the screen has run ------
intensity <- NULL
EX <- file.path(RUN_ROOT, "Secondary_QC", "Extremes", "output_files",
                "extreme_sample_master.csv")
if (file.exists(EX)) {
  xm <- utils::read.csv(EX, stringsAsFactors = FALSE, check.names = FALSE)
  if ("mean_INT" %in% names(xm)) {
    intensity <- xm$mean_INT[match(d$SAMPLE, xm$SAMPLE)]
    cat(sprintf("Whole-panel level taken from the Extremes screen (%d of %d matched)\n",
                sum(!is.na(intensity)), nrow(d)))
  } else {
    cat("NOTE: extreme_sample_master.csv has no mean_INT column — re-run the\n",
        "      Extremes screen (or its backfill) to include the panel-level axis.\n")
  }
} else {
  cat("NOTE: Extremes screen not run yet; the panel-level axis will be skipped.\n")
}

## --- build and write ---------------------------------------------------------
ax <- specimen_axes(d, intensity)
cat(sprintf("\nAxes available: %s\n", paste(names(ax), collapse = ", ")))

w <- function(x, f) if (!is.null(x) && nrow(x)) {
  utils::write.csv(x, file.path(OUT_DIR, f), row.names = FALSE)
  cat(sprintf("  %-34s %4d rows\n", f, nrow(x)))
}
cat("\nWrote:\n")
by_site <- axis_by_group(ax, d$Site, "Site");           w(by_site, "specimen_axis_by_site.csv")
by_anc  <- if ("Ancestry" %in% names(d))
             axis_by_group(ax, d$Ancestry, "Ancestry"); w(by_anc,  "specimen_axis_by_ancestry.csv")
vp      <- axis_variance(ax, d);                        w(vp,      "specimen_axis_variance.csv")
cc      <- axis_cocluster(ax, d$Site);                  w(cc,      "specimen_axis_cocluster.csv")
wi      <- axis_within_group(ax, d$Site);               w(wi,      "specimen_axis_within_site.csv")

## --- sign-balance guard: every composite index the pipeline ships ------------
PANEL_W <- list(
  ALS_MND            = c(pTDP43_409 = 1, TARDBP = 1),
  Synucleinopathy    = c(Oligo_SNCA = 1, pSNCA_129 = 1, DDC = 1, PARK7 = 1),
  AD_tau             = c(pTau_217 = 1, pTau_181 = 1, pTau_231 = 1, BD_pTau_217 = 1, MAPT = 1),
  Amyloid_low        = c(A42_A40_ratio = -1, A_42 = -1),
  BD_pTau_inflation  = c(BD_pTau_181 = 1, BD_pTau_217 = 1, BD_pTau_231 = 1, BD_MAPT = 1),
  Neurodegeneration  = c(GFAP = 1, NEFL = 1, NEFH = 1, GDF15 = 1),
  hemolysis_index    = c(HBA1 = 1, PGK1 = 1, MDH1 = 1, SOD1 = 1, ENO2 = 1))
sb <- sign_balance_report(PANEL_W); w(sb, "sign_balance_report.csv")

## --- console summary ---------------------------------------------------------
cat("\n================ SUMMARY ================\n")
cat("\n-- how much of each axis lies BETWEEN sites --\n")
vs <- vp[vp$variable == "Site", ]
for (i in order(-vs$pct_between)) cat(sprintf("  %-18s %5.1f%%\n", vs$axis[i], vs$pct_between[i]))

if (!is.null(cc)) {
  cat("\n-- do the axes single out the SAME sites? (r across site medians) --\n")
  for (i in seq_len(nrow(cc)))
    cat(sprintf("  %-18s vs %-18s %+6.2f\n", cc$axis_a[i], cc$axis_b[i], cc$r_across_groups[i]))
  cat("  -> a near-zero r means these are DIFFERENT sites, not one site effect.\n")
}

if (!is.null(wi)) {
  cat("\n-- does each axis survive comparing only WITHIN a site? --\n")
  for (i in seq_len(nrow(wi)))
    cat(sprintf("  %-18s %5.1f%% overall -> %5.1f%% within-site  (%.0f%% retained)\n",
                wi$axis[i], wi$pct_overall[i], wi$pct_within[i], wi$pct_retained[i]))
  cat("  -> retention near 100% means per-specimen, NOT a site artifact.\n")
}

cat("\n-- sign balance: which composite indices track the sample's overall level --\n")
for (i in order(-abs(sb$sign_balance)))
  cat(sprintf("  %-18s %2d markers  balance %+5.2f   %s\n", sb$score[i], sb$n_markers[i],
              sb$sign_balance[i], sb$verdict[i]))

cat("\n🚫 DESCRIPTIVE ONLY. No drop rule here, and cystatin C must never be added\n")
cat("   as a covariate to remove this — it is itself raised in clinical AD, so\n")
cat("   conditioning on it deletes disease signal. See ANALYSIS_PLAN §3a.\n")
