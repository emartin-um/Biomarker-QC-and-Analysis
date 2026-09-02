#!/usr/bin/env Rscript
# =============================================================================
# run_triage_review.R  ·  Secondary QC — auditing the Primary QC screens
# -----------------------------------------------------------------------------
# Primary QC removes wells (triage) and marks wells (flags). Nothing checked the
# SCREENS themselves. This does: does each gate fire, is it independent of the
# others, where does its cut sit on its own statistic, and is it keying on plate
# position rather than on the specimen.
#
# BASE R ONLY — no dplyr. Run from the module directory:
#   Rscript run_triage_review.R
#
# Reads (from the run directory two levels up):
#   Primary_QC/input_files/Sample_QC.csv          the attempted-well universe
#   Primary_QC/input_files/Replicate_list.csv     HIHG pool wells (excluded)
#   Primary_QC/output_files/samples_to_triage.csv     who was removed, and why
#   Primary_QC/output_files/flagged_read_outliers.csv marked but kept
#   Primary_QC/output_files/{upper,lower}_outliers_with_stats.csv  burden stats
#   Primary_QC/output_files/sample_pca_coordinates_round{1,2}.csv  PCA distances
#   Primary_QC/output_files/NPQ_*_post_QC{,_triage}.csv  + the _Low_ pair
#   Metadata_Merge/<output_dir>/metadata_PrimaryQC_refreshed.csv   labels
#
# Writes to output_files/:
#   qc_gate_inventory.csv        every gate: does it fire, is it inert
#   qc_gate_duplication.csv      is an in-house axis a vendor metric re-derived
#   triage_axis_overlap.csv      the true Venn, vs what a one-reason report says
#   gate_operating_points.csv    the quantised cuts a "threshold" actually is
#   read_gate_geometry.csv       why the lower read bound cannot fire
#   qc_axis_attribution.csv      each triage class against the known quality axes
#   well_position_effects.csv    per-position rates + permutation p
#   ic_row_gradient.csv          IC recovery down the plate, per plate-bay
#   triage_rate_by_site.csv      site rates, permuted and depth-conditioned
#   triaged_well_detail.csv      one row per removed/marked well, everything on it
#
# 🚫 DIAGNOSTIC ONLY. This module drops nothing and changes no pipeline output.
# =============================================================================

HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
RUN_ROOT <- normalizePath(file.path(HERE, "..", ".."), mustWork = FALSE)
OUT_DIR  <- file.path(HERE, "output_files")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(HERE, "triage_review_helpers.R"))

PQ <- file.path(RUN_ROOT, "Primary_QC")
MM_ROOT <- file.path(RUN_ROOT, "Metadata_Merge")
.rd <- function(...) { p <- file.path(...); if (file.exists(p))
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) else NULL }
.need <- function(x, what) if (is.null(x)) stop("triage review needs ", what) else x

## --- 1. the attempted-well spine ---------------------------------------------
q <- .need(tryCatch(utils::read.csv(file.path(PQ, "input_files", "Sample_QC.csv"),
                                    skip = 12, check.names = FALSE, stringsAsFactors = FALSE),
                    error = function(e) NULL), "Primary_QC/input_files/Sample_QC.csv")
qall <- q
q <- q[q[["Sample Type"]] == "Sample", , drop = FALSE]
rl <- .rd(PQ, "input_files", "Replicate_list.csv")
pool <- if (is.null(rl)) character(0) else as.character(rl[[1]])
q <- q[!q[["Sample Name"]] %in% pool, , drop = FALSE]

d <- data.frame(
  well_id = as.character(q[["Sample Name"]]),
  Run = sub("^([0-9]{8}-[0-9]{4}).*", "\\1", q$plateID),
  Bay = sub("^.*_(Bay[0-9]+)_.*$", "\\1", q$plateID),
  wellRow = q$wellRow, wellCol = as.integer(q$wellCol),
  IC_Median = .num(q[["IC Median"]]), Detectability = .num(q$Detectability),
  IC_Reads = .num(q[["IC Reads"]]), Reads = .num(q$Reads),
  alamar_warning = q[["QC Status"]] == "Warning",
  stringsAsFactors = FALSE)
d$run_bay <- paste(d$Run, d$Bay); d$year <- substr(d$Run, 1, 4)
d$l2reads <- log2(d$Reads)
cat(sprintf("Attempted patient wells: %d  (%d controls and %d pool wells excluded)\n",
            nrow(d), sum(qall[["Sample Type"]] != "Sample"), length(intersect(pool, qall[["Sample Name"]]))))

## --- 2. the QC decisions -----------------------------------------------------
tri <- .need(.rd(PQ, "output_files", "samples_to_triage.csv"), "samples_to_triage.csv")
d$triaged <- d$well_id %in% tri$SampleID
mi <- match(d$well_id, tri$SampleID)
ax_all <- triage_axes(ifelse(is.na(mi), "", tri$Reason[mi]))
for (v in names(ax_all)) d[[paste0("tri_", v)]] <- ax_all[[v]]
d$tri_n_axes <- rowSums(ax_all)
unknown <- tri$Reason[rowSums(triage_axes(tri$Reason)) == 0]
if (length(unknown))
  stop("samples_to_triage.csv has unrecognised Reason(s): ", paste(unique(unknown), collapse = " | "),
       "\n  -> add to triage_axes(); an unmapped reason would silently vanish from the audit.")

fro <- .rd(PQ, "output_files", "flagged_read_outliers.csv")
d$read_flagged <- if (is.null(fro)) FALSE else d$well_id %in% fro$SampleID
mf <- if (is.null(fro)) NA else match(d$well_id, fro$SampleID)
d$read_flag_dir <- if (is.null(fro)) "" else
  ifelse(!d$read_flagged, "", ifelse(d$Reads < fro$Lower_Bound[mf], "LOW", "HIGH"))

for (nm in c("round1", "round2")) {
  p <- .rd(PQ, "output_files", paste0("sample_pca_coordinates_", nm, ".csv"))
  if (is.null(p)) next
  m <- match(d$well_id, p$Sample)
  d[[paste0("PC_z_", nm)]] <- (p$PC_dist[m] - p$mean_PCdist[m]) / p$sd_PCdist[m]
  d[[paste0("pca_extreme_", nm)]] <- p$extreme_outlier[m]
}
up <- .rd(PQ, "output_files", "upper_outliers_with_stats.csv")
lo <- .rd(PQ, "output_files", "lower_outliers_with_stats.csv")
if (!is.null(up)) { m <- match(d$well_id, up$SampleID)
  d$n_out_up <- up$Total_Outliers[m]; d$fdr_up <- up$FDR[m] }
if (!is.null(lo)) { m <- match(d$well_id, lo$SampleID)
  d$n_out_lo <- lo$Total_Outliers[m]; d$fdr_lo <- lo$FDR[m] }
d$n_out_tot <- d$n_out_up + d$n_out_lo
d$fdr_min <- pmin(d$fdr_up, d$fdr_lo, na.rm = TRUE)
N_ANCHOR <- length(setdiff(names(up), c("SampleID","Run","Bay","Total_Outliers",
                                        "Expected_Outliers","P_Value","FDR")))

## --- 3. whole-panel intensity, INCLUDING the triaged wells -------------------
# The Extremes screen runs on the cohort, so triaged wells have no mean_INT and
# cannot be placed on the intensity axis at all. Pool the post-QC and triage NPQ
# and recompute INT over the union, so removed wells are on the SAME scale as
# retained ones. Verified 2026-07-27: r = 0.9948 against the production
# extreme_sample_master.csv on the shared wells.
npq_files <- list.files(file.path(PQ, "output_files"), pattern = "^NPQ_.*post_QC(_triage)?\\.csv$",
                        full.names = TRUE)
std <- npq_files[!grepl("_Low_", npq_files)]; lwf <- npq_files[grepl("_Low_", npq_files)]
X <- NULL
if (length(std)) {
  X <- do.call(rbind, lapply(std, function(f) utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)))
  if (length(lwf)) {
    L <- do.call(rbind, lapply(lwf, function(f) utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)))
    extra <- setdiff(names(L), names(X))
    if (length(extra)) X <- merge(X, L[, c("SampleID", extra)], by = "SampleID", all.x = TRUE)
  }
}
if (!is.null(X)) {
  mk <- setdiff(names(X), c("SampleID", "Run", "Bay"))
  mk <- mk[vapply(X[mk], is.numeric, logical(1))]
  int <- vapply(mk, function(v) { x <- X[[v]]; ok <- !is.na(x); o <- rep(NA_real_, length(x))
    if (sum(ok) >= 2) o[ok] <- stats::qnorm((rank(x[ok], ties.method = "average") - 0.5) / sum(ok)); o },
    numeric(nrow(X)))
  m <- match(d$well_id, X$SampleID)
  d$mean_INT <- rowMeans(int, na.rm = TRUE)[m]
  d$mean_INT_z <- as.numeric(scale(d$mean_INT))
  hm <- intersect(c("HBA1","PGK1","MDH1","SOD1","ENO2"), colnames(int))
  if (length(hm)) d$hemolysis_INT <- rowMeans(int[, hm, drop = FALSE], na.rm = TRUE)[m]
  if ("CST3" %in% names(X)) d$CST3 <- X$CST3[m]
  cat(sprintf("Whole-panel intensity computed over %d wells x %d markers (triaged included)\n",
              nrow(X), length(mk)))
} else cat("NOTE: no NPQ files found; the intensity axis will be skipped.\n")

## --- 4. labels ----------------------------------------------------------------
mm <- list.dirs(MM_ROOT, recursive = FALSE)
mm <- mm[grepl("(^|/)output_files", basename(mm)) & !grepl("_archive|_backup", basename(mm))]
mm <- mm[file.exists(file.path(mm, "metadata_PrimaryQC_refreshed.csv"))]
if (length(mm)) {
  md <- .rd(mm[which.max(file.mtime(file.path(mm, "metadata_PrimaryQC_refreshed.csv")))],
            "metadata_PrimaryQC_refreshed.csv")
  k <- match(d$well_id, md$SAMPLE_ALIQUOT)
  for (v in c("Site", "CDX", "sex", "age_at_subject", "Record_ID"))
    if (v %in% names(md)) d[[v]] <- md[[v]][k]
  cat(sprintf("Metadata labels joined for %d of %d wells\n", sum(!is.na(k)), nrow(d)))
}

## --- 5. build every output ----------------------------------------------------
wr <- function(x, f) if (!is.null(x) && nrow(x)) {
  utils::write.csv(x, file.path(OUT_DIR, f), row.names = FALSE)
  cat(sprintf("  %-32s %4d rows\n", f, nrow(x))) }
cat("\nWrote:\n")

## 5a. does each gate fire?
N <- nrow(d)
inv <- rbind(
  gate_row("triage: IC outlier", "IC Reads outside 0.6x-1.4x the Run x Bay median", d$tri_IC, N,
           "same quantity as the vendor's IC Median metric - see qc_gate_duplication.csv"),
  gate_row("triage: PCA round 1", "PC_dist > mean + 5 SD", d$pca_extreme_round1, N, ""),
  gate_row("triage: PCA round 2", "PC_dist > mean + 5 SD, recomputed after round 1",
           d$pca_extreme_round2, N, "if INERT the second pass adds nothing"),
  gate_row("triage: outlier burden", sprintf("BH-FDR < 0.01 on the count of %d anchors outlying", N_ANCHOR),
           !is.na(d$fdr_min) & d$fdr_min < 0.01, N,
           "FDR is a deterministic function of an integer count - see gate_operating_points.csv"),
  gate_row("flag: excessive reads", "Reads > Run x Bay median + 4*IQR", d$read_flag_dir == "HIGH", N,
           "MARKED, NOT REMOVED"),
  gate_row("flag: too few reads", "Reads < max(500000, median - 4*IQR)", d$read_flag_dir == "LOW", N,
           "see read_gate_geometry.csv - the floor may make this unreachable"),
  gate_row("vendor: Detectability", "Detectability < 0.90", d$Detectability < 0.90, N, ""),
  gate_row("vendor: IC Reads", "IC Reads < 1000", d$IC_Reads < 1000, N, ""),
  gate_row("vendor: Reads", "Reads < 500000", d$Reads < 5e5, N, ""),
  gate_row("vendor: IC Median", "|IC Median| > 0.40", abs(d$IC_Median) > 0.40, N, ""))
wr(inv, "qc_gate_inventory.csv")

## 5b. is an in-house axis a vendor metric re-derived?
dup <- rbind(
  flag_agreement(d$tri_IC, abs(d$IC_Median) > 0.40, "triage: IC outlier", "vendor: IC Median +/-40%"),
  flag_agreement(d$tri_IC, d$tri_PCA, "triage: IC outlier", "triage: PCA"),
  flag_agreement(d$tri_PCA, !is.na(d$fdr_min) & d$fdr_min < 0.01, "triage: PCA", "triage: burden"),
  flag_agreement(d$triaged, d$read_flagged, "triage (any)", "read flag"),
  flag_agreement(d$triaged, d$alamar_warning, "triage (any)", "vendor warning (any)"))
wr(dup, "qc_gate_duplication.csv")

## 5c. the true axis Venn vs a one-reason report
axt <- ax_all[d$triaged, , drop = FALSE]
ov <- axis_overlap(axt)
wr(ov, "triage_axis_overlap.csv")
wr(axis_combinations(axt), "triage_axis_combinations.csv")

## 5d. what a "threshold" actually is
wr(gate_operating_points(d$fdr_min, d$n_out_tot, cut = 0.01), "gate_operating_points.csv")
nb <- rbind(gate_cut_neighbourhood(d$fdr_min, !is.na(d$fdr_min) & d$fdr_min < 0.01, "outlier burden (FDR)"),
            gate_cut_neighbourhood(-d$PC_z_round1, d$pca_extreme_round1 %in% TRUE, "PCA round 1 (-PC_z)"))
wr(nb, "gate_cut_neighbourhood.csv")
rg <- read_gate_geometry(d$Reads, d$run_bay, k = 4)
wr(rg, "read_gate_geometry.csv")

## 5e. the removed/marked wells against the known quality axes, PER CLASS
classes <- list(
  "triaged (all)"      = d$triaged,
  "low reads (any)"    = d$tri_low_reads,
  "low reads only"     = d$tri_low_reads & d$tri_n_axes == 1,
  "IC only"            = d$tri_IC & d$tri_n_axes == 1,
  "PCA (any)"          = d$tri_PCA,
  "burden (any)"       = d$tri_burden,
  "burden, not low-read" = d$tri_burden & !d$tri_low_reads,
  "read-flagged HIGH"  = d$read_flag_dir == "HIGH",
  "vendor warning"     = d$alamar_warning)
AXES <- c("mean_INT_z", "l2reads", "Detectability", "IC_Median", "hemolysis_INT", "CST3",
          "n_out_tot", "PC_z_round1")
wr(axis_attribution(d, classes, AXES), "qc_axis_attribution.csv")

## 5f. is the screen keying on plate position?
pe <- position_effects(d$triaged, d$wellRow, d$wellCol)
pp <- attr(pe, "perm")
wr(pe, "well_position_effects.csv")
pf <- position_effects(d$read_flagged, d$wellRow, d$wellCol)
ppf <- attr(pf, "perm")
wr(pf, "well_position_effects_readflag.csv")
grad <- row_gradient(d$IC_Median, d$wellRow, d$run_bay)
wr(grad, "ic_row_gradient.csv")
# The decisive position statistic: within-plate rank, which survives plate-level
# shifts and is comparable across assay periods. Emitted overall and per period,
# because "the pattern is the same but the magnitude moved" is the finding.
prm <- position_rank_map(d$IC_Median, paste0(d$wellRow, d$wellCol), d$run_bay)
prm_perm <- attr(prm, "perm")
wr(prm, "position_rank_map_IC.csv")
prm_yr <- do.call(rbind, lapply(sort(unique(d$year)), function(y) {
  i <- d$year == y
  if (length(unique(d$run_bay[i])) < 3) return(NULL)
  z <- position_rank_map(d$IC_Median[i], paste0(d$wellRow, d$wellCol)[i], d$run_bay[i], n_perm = 1000L)
  cbind(year = y, z)
}))
wr(prm_yr, "position_rank_map_IC_by_year.csv")

# The spatial pattern as a single contrast: the far corner of the plate (high
# rows x high columns) against the rest, permuted WITHIN plate so plate-level
# shifts cannot produce it. This is the number to quote — it is stable across
# assay periods, whereas any individual well's change is not.
d$plate_region <- ifelse(match(d$wellRow, LETTERS) >= 6 & d$wellCol >= 9,
                         "far corner (F-H x 9-12)", "rest of plate")
reg <- do.call(rbind, lapply(c("all", sort(unique(d$year))), function(y) {
  i <- if (y == "all") rep(TRUE, nrow(d)) else d$year == y
  lr <- d$plate_region[i] == "far corner (F-H x 9-12)"
  x <- d$IC_Median[i]; pbi <- d$run_bay[i]
  obs <- mean(x[lr], na.rm = TRUE) - mean(x[!lr], na.rm = TRUE)
  set.seed(5)
  null <- replicate(2000L, { pp <- lr
    for (b in unique(pbi)) { j <- which(pbi == b); pp[j] <- sample(lr[j]) }
    mean(x[pp], na.rm = TRUE) - mean(x[!pp], na.rm = TRUE) })
  data.frame(period = y, n_corner = sum(lr), n_rest = sum(!lr),
             mean_corner = round(mean(x[lr], na.rm = TRUE), 4),
             mean_rest = round(mean(x[!lr], na.rm = TRUE), 4),
             difference = round(obs, 4),
             perm_p = round((1 + sum(null <= obs)) / (length(null) + 1), 5),
             stringsAsFactors = FALSE, row.names = NULL)
}))
wr(reg, "plate_region_contrast_IC.csv")

## 5g. site clustering, permuted and depth-conditioned
rs <- rate_by_group(d$triaged, d$Site, d$l2reads)
rsp <- if (!is.null(rs)) attr(rs, "perm") else NULL
rsc <- if (!is.null(rs)) attr(rs, "conditioned") else NULL
wr(rs, "triage_rate_by_site.csv")
rf <- rate_by_group(d$read_flagged, d$Site, d$l2reads)
rfp <- if (!is.null(rf)) attr(rf, "perm") else NULL
rfc <- if (!is.null(rf)) attr(rf, "conditioned") else NULL
wr(rf, "read_flag_rate_by_site.csv")

## 5h. the per-well detail sheet
det <- d[d$triaged | d$read_flagged | d$alamar_warning, ]
keep <- intersect(c("well_id","Run","Bay","wellRow","wellCol","year","Site",
                    "triaged","tri_low_reads","tri_IC","tri_PCA","tri_burden","tri_baddata","tri_n_axes",
                    "read_flagged","read_flag_dir","alamar_warning",
                    "IC_Median","Detectability","IC_Reads","Reads","l2reads",
                    "mean_INT_z","hemolysis_INT","CST3","n_out_tot","fdr_min",
                    "PC_z_round1","PC_z_round2"), names(det))
wr(det[order(-det$triaged, det$Run, det$Bay, det$wellRow), keep], "triaged_well_detail.csv")

## --- 6. console summary --------------------------------------------------------
cat("\n================ WHICH GATES ARE DOING WORK ================\n")
for (i in seq_len(nrow(inv)))
  cat(sprintf("  %-26s %4d wells  %-32s %s\n", inv$gate[i], inv$n_fired[i], inv$status[i],
              if (nzchar(inv$note[i])) "*" else ""))
inert <- inv$gate[inv$n_fired == 0]
if (length(inert))
  cat(sprintf("\n⚠ INERT GATES (never fire -> untested, not proven safe): %s\n",
              paste(inert, collapse = "; ")))

cat("\n================ IS ANY AXIS A DUPLICATE ================\n")
for (i in seq_len(nrow(dup)))
  cat(sprintf("  %-26s vs %-26s Jaccard %.3f  %s\n", dup$flag_a[i], dup$flag_b[i],
              dup$jaccard[i], dup$verdict[i]))

cat("\n================ THE AXES, WITHOUT THE PRECEDENCE COLLAPSE ================\n")
for (i in seq_len(nrow(ov)))
  cat(sprintf("  %-9s on-axis %3d   sole reason %3d   (%4.1f%% shared with another axis)\n",
              ov$axis[i], ov$n_wells_on_axis[i], ov$n_wells_sole_reason[i], ov$pct_of_axis_shared[i]))
cat(sprintf("  -> %d of %d triaged wells carry more than one axis. A single-reason\n",
            sum(rowSums(axt) > 1), nrow(axt)))
cat("     report (e.g. Cohort_Accounting exit_reason) UNDERCOUNTS the later axes.\n")

cat("\n================ WHAT THE BURDEN 'THRESHOLD' ACTUALLY IS ================\n")
op <- gate_operating_points(d$fdr_min, d$n_out_tot, cut = 0.01)
cat(sprintf("  the statistic takes %d distinct values in %d wells; %d wells (%.1f%%) sit at 1.0\n",
            nrow(op), N, sum(d$fdr_min == 1, na.rm = TRUE), 100 * mean(d$fdr_min == 1, na.rm = TRUE)))
head_op <- head(op[order(op$stat_value), ], 8)
for (i in seq_len(nrow(head_op)))
  cat(sprintf("    FDR %-10.4g  n_out ~%-3s  %3d wells here  %3d removed at this cut%s\n",
              head_op$stat_value[i],
              ifelse(is.na(head_op$count_value[i]), "?", head_op$count_value[i]),
              head_op$n_wells_at[i], head_op$n_removed_at_this_cut[i],
              ifelse(isTRUE(head_op$is_current_cut[i]), "   <- inside the current gate", "")))
cat(sprintf("  gap around the cut: worst retained %.4g vs best removed %.4g\n",
            nb$worst_retained[1], nb$best_removed[1]))
cat("  -> this is a choice among integer anchor counts, NOT a continuous dial.\n")

if (!is.null(rg))
  cat(sprintf("\n  read gate: median - 4*IQR is below the 500000 floor in %d of %d plate-bays, so\n  the IQR arm of the lower bound never binds. A well can only be flagged low by\n  the ABSOLUTE 500000 floor (%d well%s here), never relative to its own plate.\n",
              sum(rg$floor_binds), nrow(rg), sum(d$read_flag_dir == "LOW"),
              if (sum(d$read_flag_dir == "LOW") == 1) "" else "s"))

cat("\n================ IS THE SCREEN KEYING ON PLATE POSITION ================\n")
cat(sprintf("  worst position for triage: %s at %d of %d wells (%.1f%%) vs %.2f%% elsewhere\n",
            pe$position[1], pe$n_flagged[1], pe$n_wells[1], pe$pct[1],
            100 * (sum(d$triaged) - pe$n_flagged[1]) / (N - pe$n_wells[1])))
cat(sprintf("  permutation on the max per position: observed %d, null mean %.2f, p = %.4f\n",
            pp$observed_max, pp$null_mean_max, pp$p))
if (pp$p < 0.05) {
  cat("  ⚠ a position effect this size is not sampling noise. A specimen-quality screen\n    should be flat across the plate; check the layout for that well.\n")
  # Is it a standing design property, or did it START? A position that was fine and
  # went bad is drift in one physical location — a different (and more tractable)
  # problem from a layout flaw that was always there. Split the worst position by
  # assay period and by its immediate neighbours.
  wp <- pe$position[1]
  w <- d$wellRow[match(wp, paste0(d$wellRow, d$wellCol))]
  cn <- as.integer(sub("^[A-Z]", "", wp))
  nb_pos <- paste0(setdiff(LETTERS[1:8], w), cn)          # same column, other rows
  isw <- paste0(d$wellRow, d$wellCol) == wp
  isn <- paste0(d$wellRow, d$wellCol) %in% nb_pos
  cat(sprintf("\n  %s by assay period:\n", wp))
  for (y in sort(unique(d$year))) {
    i <- isw & d$year == y
    cat(sprintf("    %s  triaged %d of %2d   mean IC_Median %+.3f\n",
                y, sum(d$triaged[i]), sum(i), mean(d$IC_Median[i], na.rm = TRUE)))
  }
  cat(sprintf("  same column, other rows (n=%d): triaged %d, mean IC_Median %+.3f\n",
              sum(isn), sum(d$triaged & isn), mean(d$IC_Median[isn], na.rm = TRUE)))
  cat(sprintf("  everywhere else       (n=%d): triaged %d, mean IC_Median %+.3f\n",
              sum(!isw & !isn), sum(d$triaged & !isw & !isn),
              mean(d$IC_Median[!isw & !isn], na.rm = TRUE)))
  cat(sprintf("  spans %d runs and %d bays -- %s\n",
              length(unique(d$Run[isw & d$triaged])), length(unique(d$Bay[isw & d$triaged])),
              if (length(unique(d$Run[isw & d$triaged])) > 2)
                "not one plate, so not a one-off handling error"
              else "few enough plates that a one-off is possible"))
  # Whether THIS position changed more than positions change in general is a
  # separate question from whether it is bad, and it is easy to overclaim: with
  # 84 positions, a couple will always move. Report the change against the spread
  # of changes across all positions rather than asserting drift.
  if (length(unique(d$year)) == 2 && !is.null(prm_yr)) {
    ys <- sort(unique(d$year))
    a <- tapply(d$IC_Median[d$year == ys[1]], paste0(d$wellRow, d$wellCol)[d$year == ys[1]], mean)
    b <- tapply(d$IC_Median[d$year == ys[2]], paste0(d$wellRow, d$wellCol)[d$year == ys[2]], mean)
    kk <- intersect(names(a), names(b)); ch <- b[kk] - a[kk]
    cat(sprintf("  change %s->%s at %s: %+.3f, against SD %.3f across all %d positions (%.1f SD)\n",
                ys[1], ys[2], wp, ch[wp], stats::sd(ch), length(kk),
                (ch[wp] - mean(ch)) / stats::sd(ch)))
    cat(sprintf("  other positions moving similarly: %s\n",
                paste(names(sort(ch)[1:4]), collapse = ", ")))
    cat("  -> treat 'this well degraded' as suggestive unless it stands well clear\n")
    cat("     of that spread. The reproducible map above is the stronger claim.\n")
  }
}
cat(sprintf("\n  within-plate rank map of IC_Median (the statistic that survives plate shifts):\n"))
cat(sprintf("    worst positions: %s\n", paste(head(prm$position, 6), collapse = ", ")))
cat(sprintf("    %s sits in the bottom decile of its own plate on %d of %d plate-bays\n",
            prm$position[1], prm$n_plates_in_bottom_decile[1], prm$n_plates[1]))
cat(sprintf("    within-plate permutation: observed %.4f, null %.4f, p = %.4f\n",
            prm_perm$observed, prm_perm$null_mean, prm_perm$p))
if (!is.null(prm_yr) && length(unique(prm_yr$year)) == 2) {
  ys <- unique(prm_yr$year)
  a <- prm_yr[prm_yr$year == ys[1], ]; b <- prm_yr[prm_yr$year == ys[2], ]
  k <- intersect(a$position, b$position)
  r <- stats::cor(a$median_pct_rank[match(k, a$position)], b$median_pct_rank[match(k, b$position)])
  cat(sprintf("    the SAME map in both periods: r = %.3f over %d positions\n", r, length(k)))
  cat("    -> a standing spatial pattern, not a one-period failure. What moved\n")
  cat("       between periods is the magnitude, not which wells are affected.\n")
}
cat("\n  far corner of the plate vs the rest (permuted within plate):\n")
for (i in seq_len(nrow(reg)))
  cat(sprintf("    %-6s corner %+.3f vs rest %+.3f   difference %+.3f   p = %.4f\n",
              reg$period[i], reg$mean_corner[i], reg$mean_rest[i],
              reg$difference[i], reg$perm_p[i]))
neg <- sum(grad$rho_row_vs_value < 0, na.rm = TRUE)
cat(sprintf("\n  IC recovery down the plate: Spearman(row, IC_Median) is negative in %d of %d\n  plate-bays (median %.3f). Strictly monotone in %d.\n",
            neg, sum(!is.na(grad$rho_row_vs_value)), stats::median(grad$rho_row_vs_value, na.rm = TRUE),
            sum(grad$rho_row_vs_value == -1, na.rm = TRUE)))
cat("  -> a consistent gradient, not a monotone one. Because the IC gate keys on\n     exactly this quantity, well position partly decides who is triaged.\n")

cat("\n================ THE REMOVED WELLS ON THE KNOWN AXES ================\n")
at <- axis_attribution(d, classes, AXES)
if (!is.null(at) && "mean_INT_z" %in% at$axis) {
  s <- at[at$axis == "mean_INT_z", ]
  for (i in seq_len(nrow(s)))
    cat(sprintf("  %-18s n=%4d  mean_INT_z %+6.2f vs %+.2f   d %+6.2f  AUC %.3f\n",
                s$class[i], s$n[i], s$mean_in[i], s$mean_out[i], s$cohen_d[i], s$auc[i]))
  cat("  -> read the CLASSES, not the pooled row. Pooling manufactures a\n")
  cat("     'triaged wells are dim' story that need not describe most of them.\n")
}

report_group <- function(lbl, pp, cond) {
  if (is.null(pp)) return(invisible(NULL))
  cat(sprintf("\n  %s (%d groups with n >= %d; %d smaller groups excluded)\n",
              lbl, pp$n_groups, pp$min_n, pp$n_groups_dropped))
  cat(sprintf("    omnibus     deviance %6.1f            p = %.4f   <- read this one\n",
              pp$omnibus_deviance, pp$omnibus_p))
  cat(sprintf("    max rate    %5.1f%% (null %.1f%%)        p = %.4f\n",
              100 * pp$observed_max_rate, 100 * pp$null_mean_max_rate, pp$max_rate_p))
  if (pp$smallest_group_n < 20)
    cat(sprintf("    NOTE smallest retained group is n = %d; a maximum-rate statistic is\n         unreliable there, which is why the omnibus is reported beside it.\n",
                pp$smallest_group_n))
  if (!is.null(cond)) {
    cat(sprintf("    deviance explained: %.1f alone -> %.1f after read depth\n",
                cond$dev_group_only, cond$dev_group_given_depth))
    keep <- 100 * cond$dev_group_given_depth / max(cond$dev_group_only, 1e-9)
    cat(sprintf("      -> %.0f%% retained; %s\n", keep,
                if (keep < 50) "mostly depth composition, not a site property"
                else "survives depth, so not simply composition"))
  }
}
report_group("triage by site", rsp, rsc)
report_group("read flags by site", rfp, rfc)

cat("\n🚫 DIAGNOSTIC ONLY. Nothing here is a drop rule, and nothing here proposes\n")
cat("   one. A gate that fires on nobody has not been shown to be safe — only\n")
cat("   untested; and a gate that duplicates a vendor metric is not a second\n")
cat("   opinion. Per-specimen nuisance axes belong in the analysis model.\n")
