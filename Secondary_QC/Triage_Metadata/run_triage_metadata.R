#!/usr/bin/env Rscript
# =============================================================================
# run_triage_metadata.R  ·  Secondary QC — do the QC decisions know who the
#                            specimen came from?
# -----------------------------------------------------------------------------
# Primary QC removes wells (triage) and marks wells (flags) using assay
# statistics only — it is agnostic to phenotype by design, and the collection-site
# label does not even exist until the metadata merge. So nothing in the pipeline
# has ever asked whether those decisions land evenly across the people the
# specimens came from.
#
# Two different stakes:
#   CLINICAL COLLECTION SITE — a site whose specimens fail QC more often is a
#     collection, handling or shipping problem. Actionable, and a different owner
#     from a plate problem or a well-position problem.
#   DIAGNOSIS — differential triage by case/control status biases every
#     downstream effect estimate. This one can corrupt the analysis silently.
#
# BASE R ONLY — no dplyr. Run from the module directory:
#   Rscript run_triage_metadata.R
#
# Reads (from the run directory two levels up):
#   Primary_QC/output_files/well_level_QC.csv      preferred: one row per well
#     ...or, for runs made before that existed, reconstructs the same spine from
#   Primary_QC/input_files/Sample_QC.csv + output_files/samples_to_triage.csv
#                                        + output_files/flagged_read_outliers.csv
#   Metadata_Merge/<output_dir>/metadata_PrimaryQC_refreshed.csv   the labels
#
# 🚫 DIAGNOSTIC ONLY. Drops nothing, flags nothing, changes no pipeline output.
# =============================================================================

HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
RUN_ROOT <- normalizePath(file.path(HERE, "..", ".."), mustWork = FALSE)
OUT_DIR  <- file.path(HERE, "output_files")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(HERE, "triage_metadata_helpers.R"))
source(file.path(HERE, "triage_metadata_report.R"))

PQ      <- file.path(RUN_ROOT, "Primary_QC")
MM_ROOT <- file.path(RUN_ROOT, "Metadata_Merge")
.rd <- function(...) { p <- file.path(...); if (file.exists(p))
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) else NULL }

# ---------------------------------------------------------------------------
# Settings. Nothing below this block is specific to any one study.
# ---------------------------------------------------------------------------
MIN_N  <- 20L            # group size floor. Below it a rate is noise -- but the
                         # groups below the floor are still REPORTED, never dropped.
N_PERM <- 10000L         # permutation replicates
AGE_BREAKS <- c(60, 70, 80)   # age-band cutpoints; widen for a younger cohort

# Metadata columns to test, as name = column-in-the-merge. Add or remove freely;
# a column that is absent, constant, or all-NA is skipped with a note rather than
# erroring, so this list does not have to match every study's merge exactly.
GROUP_VARS <- c(Site = "Site", Diagnosis = "Case_Control", sex = "sex",
                Ancestry = "<derived>", age_band = "<derived>")

# Continuous covariates for the balance table.
BALANCE_CONTINUOUS <- c("age", "BMI", "l2reads", "mean_INT_z")

cat("=========================================================\n")
cat(" Triage/flag decisions vs specimen metadata\n")
cat("=========================================================\n\n")

## --- 1. the well spine, with every QC decision on it -------------------------
w <- .rd(PQ, "output_files", "well_level_QC.csv")
if (!is.null(w)) {
  cat("Spine: Primary_QC/output_files/well_level_QC.csv\n")
  d <- data.frame(well_id = as.character(w$SampleID), Run = w$Run, Bay = w$Bay,
                  Well = w$Well, Reads = .num(w$Reads),
                  triaged = as.logical(w$triaged), read_flagged = as.logical(w$read_flagged),
                  tri_IC = as.logical(w$tri_IC), tri_PCA = as.logical(w$tri_PCA),
                  tri_burden = as.logical(w$tri_burden),
                  mean_INT_z = .num(w$mean_INT_z), stringsAsFactors = FALSE)
} else {
  # Older runs: rebuild the same columns so this module still works on them.
  cat("Spine: reconstructed from Sample_QC.csv + samples_to_triage.csv\n")
  cat("  (no well_level_QC.csv -- this run predates the 2026-07-30 Primary QC)\n")
  p  <- file.path(PQ, "input_files", "Sample_QC.csv")
  pk <- readLines(p, n = 40, warn = FALSE)
  sk <- which(grepl('(^|,)\\s*"?Sample Name"?\\s*(,|$)', pk))[1] - 1L
  if (is.na(sk)) sk <- 12L
  q <- utils::read.csv(p, skip = sk, check.names = FALSE, stringsAsFactors = FALSE)
  q <- q[q[["Sample Type"]] == "Sample", , drop = FALSE]
  rl <- .rd(PQ, "input_files", "Replicate_list.csv")
  if (!is.null(rl)) q <- q[!q[["Sample Name"]] %in% as.character(rl[[1]]), , drop = FALSE]
  tri <- .rd(PQ, "output_files", "samples_to_triage.csv")
  fro <- .rd(PQ, "output_files", "flagged_read_outliers.csv")
  if (is.null(tri)) stop("need Primary_QC/output_files/samples_to_triage.csv")
  rs <- ifelse(is.na(match(q[["Sample Name"]], tri$SampleID)), "",
               tri$Reason[match(q[["Sample Name"]], tri$SampleID)])
  d <- data.frame(
    well_id = as.character(q[["Sample Name"]]),
    Run = sub("^([0-9]{8}-[0-9]{4}).*", "\\1", q$plateID),
    Bay = sub("^.*_(Bay[0-9]+)_.*$", "\\1", q$plateID),
    Well = if ("AUTO_WELLPOSITION" %in% names(q)) as.character(q$AUTO_WELLPOSITION)
           else paste0(q$wellRow, q$wellCol),
    Reads = .num(q$Reads),
    triaged = q[["Sample Name"]] %in% tri$SampleID,
    read_flagged = if (is.null(fro)) FALSE else q[["Sample Name"]] %in% fro$SampleID,
    tri_IC = grepl("IC outlier", rs), tri_PCA = grepl("PCA", rs),
    tri_burden = grepl("burden", rs), mean_INT_z = NA_real_,
    stringsAsFactors = FALSE)
}
d$plate   <- paste(d$Run, d$Bay)
d$l2reads <- log2(pmax(d$Reads, 1))
cat(sprintf("  %d wells, %d plate-bays, %d triaged, %d read-flagged\n\n",
            nrow(d), length(unique(d$plate)), sum(d$triaged), sum(d$read_flagged)))

## --- 2. the labels ------------------------------------------------------------
mm <- list.dirs(MM_ROOT, recursive = FALSE)
mm <- mm[grepl("(^|/)output_files", basename(mm)) & !grepl("_archive|_backup", basename(mm))]
mm <- mm[file.exists(file.path(mm, "metadata_PrimaryQC_refreshed.csv"))]
if (!length(mm)) stop("no Metadata_Merge/*/metadata_PrimaryQC_refreshed.csv found -- ",
                      "this module is entirely about the metadata; there is nothing to do without it.")
MM_DIR <- mm[which.max(file.mtime(file.path(mm, "metadata_PrimaryQC_refreshed.csv")))]
md <- .rd(MM_DIR, "metadata_PrimaryQC_refreshed.csv")
k <- match(d$well_id, md$SAMPLE_ALIQUOT)

.col <- function(nm) if (nm %in% names(md)) md[[nm]][k] else rep(NA, length(k))

# Diagnosis uses the merge's own curated CDX table when present. `CDX` is the
# column that table is keyed on; `Case_Control` is the fallback.
cdxrev <- .rd(dirname(MM_DIR), "review", "CDX_review.csv")
if (is.null(cdxrev)) cdxrev <- .rd(MM_ROOT, "review", "CDX_review.csv")
.dx_src <- if (!is.null(cdxrev) && "CDX" %in% names(md)) "CDX" else GROUP_VARS[["Diagnosis"]]
d$Diagnosis <- collapse_diagnosis(.col(.dx_src),
                                  cdx_review = if (.dx_src == "CDX") cdxrev else NULL)
cat(sprintf("Diagnosis collapsed from `%s` using %s.\n", .dx_src,
            if (.dx_src == "CDX") "review/CDX_review.csv" else "the fallback matcher"))

d$Site      <- .blank_to_na(.col(GROUP_VARS[["Site"]]))
d$sex       <- .blank_to_na(.col(GROUP_VARS[["sex"]]))
d$age       <- .num(.col("age_at_subject"))
d$age_band  <- age_band(d$age, AGE_BREAKS)
d$BMI       <- .num(.col("BMI"))
d$has_meta  <- !is.na(k)

# Ancestry comes from the merge's own curated rule table, not from a second
# derivation invented here. Without it the column is left NA and every ancestry
# test is skipped — better than a column of "Unknown" that looks like an answer.
ancrev <- .rd(dirname(MM_DIR), "review", "Ancestry_review.csv")
if (is.null(ancrev)) ancrev <- .rd(MM_ROOT, "review", "Ancestry_review.csv")
if (is.null(ancrev)) {
  d$Ancestry <- NA_character_
  cat("NOTE: Metadata_Merge/review/Ancestry_review.csv not found - ancestry tests skipped.\n")
} else {
  mdn <- data.frame(Group = if ("Group" %in% names(md)) md$Group[k] else NA,
                    Race = if ("Race" %in% names(md)) md$Race[k] else NA,
                    Ethnicity = if ("Ethnicity" %in% names(md)) md$Ethnicity[k] else NA,
                    stringsAsFactors = FALSE)
  d$Ancestry <- assign_ancestry_from_review(mdn, ancrev)
  d$Ancestry[!d$has_meta] <- NA_character_
}

# Wells with no metadata row are their own class, reported rather than dropped.
# They are not neutral: a well can be absent from the metadata for reasons that
# correlate with how it was collected.
cat(sprintf("Metadata joined for %d of %d wells (%.1f%%); %d have no metadata row.\n",
            sum(d$has_meta), nrow(d), 100 * mean(d$has_meta), sum(!d$has_meta)))
unl <- data.frame(
  class = c("has metadata", "no metadata row"),
  n = c(sum(d$has_meta), sum(!d$has_meta)),
  n_triaged = c(sum(d$triaged & d$has_meta), sum(d$triaged & !d$has_meta)),
  pct_triaged = round(100 * c(mean(d$triaged[d$has_meta]),
                              if (any(!d$has_meta)) mean(d$triaged[!d$has_meta]) else NA), 2),
  n_read_flagged = c(sum(d$read_flagged & d$has_meta), sum(d$read_flagged & !d$has_meta)),
  stringsAsFactors = FALSE)
utils::write.csv(unl, file.path(OUT_DIR, "unlabelled_wells.csv"), row.names = FALSE)
print(unl, row.names = FALSE)
cat("\n")

## --- 3. can the confounded questions even be asked? ---------------------------
# Site ships in batches, so "this site fails more" and "this site's samples
# landed on worse plates" can be the same data. Check separability BEFORE
# reporting any site effect -- if site were nested in plate, no amount of
# modelling would separate them and the honest answer would be "unanswerable".
cat("--- design separability -------------------------------------------------\n")
conf <- rbind(
  confounding_summary(d$Site, d$plate, "Site", "plate"),
  confounding_summary(d$Site, d$Diagnosis, "Site", "Diagnosis"),
  confounding_summary(d$Diagnosis, d$plate, "Diagnosis", "plate"))
utils::write.csv(conf, file.path(OUT_DIR, "design_confounding.csv"), row.names = FALSE)
for (i in seq_len(nrow(conf)))
  cat(sprintf("  %-10s x %-10s  V = %.3f   %s\n", conf$factor_a[i], conf$factor_b[i],
              conf$cramers_v[i], conf$verdict[i]))
cat("\n")

## --- 4. the axes, never pooled -----------------------------------------------
# The triage axes are near-disjoint and point in opposite directions on assay
# quality (Triage_Review S1). A metadata association tested on "triaged" as a
# lump would average three different screens into one uninterpretable number.
AXES <- list(
  "any triage"  = d$triaged,
  "IC"          = d$tri_IC,
  "PCA"         = d$tri_PCA,
  "burden"      = d$tri_burden,
  "read flag"   = d$read_flagged)
GROUPS <- list(Site = d$Site, Diagnosis = d$Diagnosis, sex = d$sex,
               Ancestry = d$Ancestry, age_band = d$age_band)
# A variable that is absent, all-NA, or constant cannot be tested. Drop it here
# with a note, so a study whose merge lacks one of these still runs.
.usable <- vapply(GROUPS, function(v)
  sum(!is.na(v)) >= MIN_N && length(unique(stats::na.omit(v))) >= 2, logical(1))
if (any(!.usable))
  cat(sprintf("Not testable on this run (absent, constant, or all missing): %s\n",
              paste(names(GROUPS)[!.usable], collapse = ", ")))
GROUPS <- GROUPS[.usable]

## 4a. per-group rates, with CIs and the below-floor groups kept visible
cat("--- rates by group ------------------------------------------------------\n")
rate_rows <- list(); floor_rows <- list()
for (gn in names(GROUPS)) for (an in names(AXES)) {
  r <- rates_by_group(AXES[[an]], GROUPS[[gn]], min_n = MIN_N)
  if (is.null(r) || !nrow(r)) next
  bf <- attr(r, "below_floor")
  r$variable <- gn; r$axis <- an; r$above_floor <- TRUE
  rate_rows[[length(rate_rows) + 1]] <- r
  if (!is.null(bf) && nrow(bf)) {
    bf$variable <- gn; bf$axis <- an; bf$above_floor <- FALSE
    floor_rows[[length(floor_rows) + 1]] <- bf
  }
}
rates <- do.call(rbind, c(rate_rows, floor_rows))
rates <- rates[, c("variable", "axis", "group", "n", "n_flagged", "pct",
                   "ci_lo_pct", "ci_hi_pct", "above_floor")]
utils::write.csv(rates[order(rates$variable, rates$axis, -rates$pct), ],
                 file.path(OUT_DIR, "triage_rates_by_metadata.csv"), row.names = FALSE)

site_any <- rates[rates$variable == "Site" & rates$axis == "any triage", ]
cat(sprintf("  Site x any triage: %d sites at or above n=%d, %d below the floor.\n",
            sum(site_any$above_floor), MIN_N, sum(!site_any$above_floor)))
cat("  Highest site rates (all sites, floor status shown):\n")
sa <- site_any[order(-site_any$pct), ]
for (i in seq_len(min(6, nrow(sa))))
  cat(sprintf("    %-6s %3d wells  %5.2f%%  [%5.2f, %5.2f]  %s\n", sa$group[i], sa$n[i],
              sa$pct[i], sa$ci_lo_pct[i], sa$ci_hi_pct[i],
              ifelse(sa$above_floor[i], "", "<-- BELOW FLOOR, not tested")))
cat("\n")

## 4b. the tests: two statistics, then conditioned
cat("--- association tests ---------------------------------------------------\n")
test_rows <- list()
for (gn in names(GROUPS)) for (an in names(AXES)) {
  fl <- AXES[[an]]; gv <- GROUPS[[gn]]
  if (sum(fl, na.rm = TRUE) < 2) next
  tt <- perm_group_test(fl, gv, stratum = NULL, min_n = MIN_N, n_perm = N_PERM)
  if (is.na(tt$observed_max_rate[1])) next
  # Stratified by plate: does the group still matter once each plate is compared
  # only with itself?
  ts <- perm_group_test(fl, gv, stratum = d$plate, min_n = MIN_N, n_perm = N_PERM)
  cd <- conditioned_deviance(fl, gv, plate = d$plate, depth = d$l2reads, min_n = MIN_N)

  base <- mean(fl, na.rm = TRUE)
  nn   <- rates_by_group(fl, gv, min_n = MIN_N)
  mde  <- if (!is.null(nn) && nrow(nn))
    min_detectable_rr(stats::median(nn$n), base, sum(nn$n)) else NA_real_

  test_rows[[length(test_rows) + 1]] <- data.frame(
    variable = gn, axis = an, n_groups_tested = tt$n_groups, n_flagged = tt$n_flagged,
    observed_max_rate = tt$observed_max_rate, p_max_rate = tt$p_max_rate,
    observed_deviance = tt$observed_deviance, p_deviance = tt$p_deviance,
    p_max_rate_within_plate = ts$p_max_rate, p_deviance_within_plate = ts$p_deviance,
    dev_group_only = if (is.null(cd)) NA_real_ else cd$dev_group_only,
    dev_given_plate = if (is.null(cd) || is.null(cd$dev_group_given_plate)) NA_real_
                      else cd$dev_group_given_plate,
    dev_given_depth = if (is.null(cd) || is.null(cd$dev_group_given_depth)) NA_real_
                      else cd$dev_group_given_depth,
    min_detectable_rate_ratio = mde,
    stringsAsFactors = FALSE)
}
tests <- do.call(rbind, test_rows)

# A null is only meaningful next to what it could have detected.
tests$reading <- ifelse(
  is.na(tests$p_deviance), "not testable",
  ifelse(tests$p_deviance < 0.05 | tests$p_max_rate < 0.05,
    ifelse(tests$p_deviance < 0.05,
           "association: omnibus, evidence spread across groups",
           "association: max-rate only -- ONE group is carrying it, check its n"),
    sprintf("no association detected; could only have seen a %.1fx elevation",
            tests$min_detectable_rate_ratio)))
utils::write.csv(tests, file.path(OUT_DIR, "metadata_association_tests.csv"), row.names = FALSE)

for (gn in names(GROUPS)) {
  cat(sprintf("  %s\n", gn))
  tg <- tests[tests$variable == gn, ]
  for (i in seq_len(nrow(tg)))
    cat(sprintf("    %-11s p_max %.4f  p_dev %.4f  | within-plate p_dev %.4f | dev %.1f -> %.1f (plate) / %.1f (depth)\n",
                tg$axis[i], tg$p_max_rate[i], tg$p_deviance[i], tg$p_deviance_within_plate[i],
                tg$dev_group_only[i], tg$dev_given_plate[i], tg$dev_given_depth[i]))
}
cat("\n")

## --- 5. the headline: collection site ----------------------------------------
cat("--- collection site: the headline ---------------------------------------\n")
ts <- tests[tests$variable == "Site" & tests$axis == "any triage", ]
if (nrow(ts)) {
  cat(sprintf("  %d sites at or above n=%d.  max rate %.2f%% (p = %.4f), omnibus deviance %.2f (p = %.4f)\n",
              ts$n_groups_tested, MIN_N, ts$observed_max_rate, ts$p_max_rate,
              ts$observed_deviance, ts$p_deviance))
  cat(sprintf("  within plate: p_dev = %.4f    deviance %.1f -> %.1f given plate, %.1f given depth\n",
              ts$p_deviance_within_plate, ts$dev_group_only, ts$dev_given_plate, ts$dev_given_depth))
  cat(sprintf("  -> %s\n", ts$reading))
  if (!is.na(ts$p_max_rate) && !is.na(ts$p_deviance) && ts$p_max_rate < 0.05 && ts$p_deviance >= 0.05)
    cat("  NOTE: the maximum rate is significant while the omnibus is not. One site is\n",
        "       carrying the whole result -- read its n and CI before acting on it.\n", sep = "")
}
cat("\n")

## --- 5b. which correlated label is carrying it? -------------------------------
# Site, ancestry, diagnosis and age are not independent — a site recruits a
# particular population, whose diagnosis and age composition follow. So the
# one-at-a-time table above will light up several of them off one underlying
# cause. Fit them together and report each given the others.
cat("--- mutual adjustment: each label given the others -----------------------\n")
adj_rows <- list()
for (an in names(AXES)) {
  if (sum(AXES[[an]], na.rm = TRUE) < 5) next
  ma <- mutual_adjustment(AXES[[an]], GROUPS, min_n = MIN_N, depth = d$l2reads)
  if (is.null(ma)) next
  ma$axis <- an
  adj_rows[[length(adj_rows) + 1]] <- ma
}
adj <- do.call(rbind, adj_rows)
if (!is.null(adj)) {
  adj <- adj[, c("axis", "variable", "df", "n_wells", "n_flagged", "dev_alone",
                 "dev_given_others", "p_given_others", "retained_share")]
  utils::write.csv(adj, file.path(OUT_DIR, "mutual_adjustment.csv"), row.names = FALSE)
  for (an in unique(adj$axis)) {
    cat(sprintf("  %s\n", an))
    aa <- adj[adj$axis == an, ]
    for (i in seq_len(nrow(aa)))
      cat(sprintf("    %-12s df %2d  dev alone %6.1f -> given others %6.1f  (keeps %s)  p = %.4g\n",
                  aa$variable[i], aa$df[i], aa$dev_alone[i], aa$dev_given_others[i],
                  ifelse(is.na(aa$retained_share[i]), "  NA",
                         sprintf("%3.0f%%", 100 * aa$retained_share[i])),
                  aa$p_given_others[i]))
  }
  cat("  A label that keeps its deviance is contributing something the others do not.\n")
  cat("  A label that collapses was a proxy for one of them.\n")
}
cat("\n")

## --- 6. differential triage: the bias question -------------------------------
# If triage removes cases at a different rate from controls, every downstream
# effect estimate is biased -- and unlike a site effect, nothing later in the
# pipeline would reveal it.
cat("--- differential triage by diagnosis ------------------------------------\n")
dx <- rates[rates$variable == "Diagnosis" & rates$axis == "any triage" & rates$above_floor, ]
dx <- dx[order(-dx$pct), ]
for (i in seq_len(nrow(dx)))
  cat(sprintf("    %-15s %4d wells  %5.2f%%  [%5.2f, %5.2f]\n",
              dx$group[i], dx$n[i], dx$pct[i], dx$ci_lo_pct[i], dx$ci_hi_pct[i]))
td <- tests[tests$variable == "Diagnosis" & tests$axis == "any triage", ]
if (nrow(td)) cat(sprintf("  -> %s\n", td$reading))

# Case vs Control specifically, as a rate ratio with a CI: the number that goes
# into a methods section.
cc <- d[d$Diagnosis %in% c("Case (AD)", "Control") & !is.na(d$Diagnosis), ]
if (nrow(cc) > 20) {
  a <- cc$Diagnosis == "Case (AD)"
  k1 <- sum(cc$triaged & a); n1 <- sum(a); k0 <- sum(cc$triaged & !a); n0 <- sum(!a)
  if (k1 > 0 && k0 > 0) {
    lrr <- log((k1 / n1) / (k0 / n0))
    se  <- sqrt(1 / k1 - 1 / n1 + 1 / k0 - 1 / n0)
    cat(sprintf("  Case vs Control triage rate ratio: %.2f  (95%% CI %.2f to %.2f)  [%d/%d vs %d/%d]\n",
                exp(lrr), exp(lrr - 1.96 * se), exp(lrr + 1.96 * se), k1, n1, k0, n0))
    cat("  A CI spanning 1 does not establish equivalence -- read it with the MDE above.\n")
  }
}
cat("\n")

## --- 7. who got removed vs who stayed ----------------------------------------
.bal_cont <- intersect(BALANCE_CONTINUOUS, names(d))
.bal_cat  <- setdiff(names(GROUPS), "Site")   # Site has its own per-site table
bal <- balance_table(d, d$triaged, continuous = .bal_cont, categorical = .bal_cat)
utils::write.csv(bal, file.path(OUT_DIR, "balance_triaged_vs_retained.csv"), row.names = FALSE)
balf <- balance_table(d, d$read_flagged, continuous = .bal_cont, categorical = .bal_cat)
utils::write.csv(balf, file.path(OUT_DIR, "balance_flagged_vs_unflagged.csv"), row.names = FALSE)

cat("--- balance: triaged vs retained (|std diff| >= 0.10 shown) --------------\n")
bb <- bal[!is.na(bal$std_diff) & abs(bal$std_diff) >= 0.10, ]
if (nrow(bb)) {
  for (i in seq_len(nrow(bb)))
    cat(sprintf("    %-12s %-16s  triaged %8.2f  retained %8.2f   std diff %+.3f\n",
                bb$variable[i], bb$level[i], bb$val_in[i], bb$val_out[i], bb$std_diff[i]))
} else cat("    nothing above 0.10.\n")
cat("\n")

## --- 8. the diagnosis collapse, written out so it can be argued with ---------
cm <- collapse_map(.col(.dx_src), cdx_review = if (.dx_src == "CDX") cdxrev else NULL)
utils::write.csv(cm, file.path(OUT_DIR, "diagnosis_collapse_map.csv"), row.names = FALSE)
cat(sprintf("--- diagnosis collapse: %d raw values -> %d levels via %s (see diagnosis_collapse_map.csv)\n\n",
            nrow(cm), length(unique(cm$collapsed)), cm$source[1]))

## --- 9. per-site detail for follow-up ----------------------------------------
# One row per site with every axis side by side, so a site conversation has the
# whole picture rather than one number.
sites <- sort(unique(stats::na.omit(d$Site)))
sd_rows <- do.call(rbind, lapply(sites, function(s) {
  i <- !is.na(d$Site) & d$Site == s
  ci <- wilson_ci(sum(d$triaged & i), sum(i))
  data.frame(Site = s, n_wells = sum(i), n_plates = length(unique(d$plate[i])),
             n_triaged = sum(d$triaged & i),
             pct_triaged = round(100 * mean(d$triaged[i]), 2),
             ci_lo_pct = round(100 * ci[1], 2), ci_hi_pct = round(100 * ci[2], 2),
             n_IC = sum(d$tri_IC & i), n_PCA = sum(d$tri_PCA & i),
             n_burden = sum(d$tri_burden & i),
             n_read_flagged = sum(d$read_flagged & i),
             median_l2reads = round(stats::median(d$l2reads[i], na.rm = TRUE), 2),
             median_mean_INT_z = round(stats::median(d$mean_INT_z[i], na.rm = TRUE), 2),
             above_floor = sum(i) >= MIN_N,
             stringsAsFactors = FALSE)
}))
utils::write.csv(sd_rows[order(-sd_rows$pct_triaged), ],
                 file.path(OUT_DIR, "site_detail.csv"), row.names = FALSE)

## --- 10. synthesis ------------------------------------------------------------
# Built from the fitted numbers, not written in advance, so it cannot drift away
# from what the module actually found on this run.
cat("=========================================================\n")
cat(" WHAT SURVIVED ADJUSTMENT\n")
cat("=========================================================\n")
if (!is.null(adj)) {
  syn <- adj[adj$variable != "log2(reads)" & !is.na(adj$p_given_others) &
               adj$p_given_others < 0.05, ]
  if (nrow(syn)) {
    cat("Associations that hold with the other labels and read depth in the model:\n")
    for (i in seq_len(nrow(syn)))
      cat(sprintf("  %-11s x %-10s  p = %.4g   (keeps %.0f%% of its unadjusted deviance)\n",
                  syn$axis[i], syn$variable[i], syn$p_given_others[i],
                  100 * syn$retained_share[i]))
  } else cat("  None. Every metadata association collapsed once the others were included.\n")

  coll <- adj[adj$variable != "log2(reads)" & !is.na(adj$retained_share) &
                adj$retained_share < 0.5, ]
  if (nrow(coll)) {
    cat("\nAssociations that COLLAPSED under adjustment - these were proxies:\n")
    for (i in seq_len(nrow(coll)))
      cat(sprintf("  %-11s x %-10s  deviance %.1f -> %.1f (keeps %.0f%%)\n",
                  coll$axis[i], coll$variable[i], coll$dev_alone[i],
                  coll$dev_given_others[i], 100 * coll$retained_share[i]))
  }
  dep <- adj[adj$variable == "log2(reads)", ]
  if (nrow(dep)) {
    cat("\nRead depth, for comparison (it is in the model above):\n")
    for (i in seq_len(nrow(dep)))
      cat(sprintf("  %-11s  deviance given the labels %.1f   p = %.3g\n",
                  dep$axis[i], dep$dev_given_others[i], dep$p_given_others[i]))
  }
  cat("\nA 'keeps' above 100% is not an error: with correlated labels a term can\n")
  cat("explain more once the others absorb variation it was competing with.\n")
}
cat("\n")

## --- 11. the run's findings as a standalone HTML ------------------------------
# Findings belong in a per-run report, not in the module's README — a README that
# reports results is out of date the moment it meets a different dataset.
run_label <- basename(RUN_ROOT)
rpt <- file.path(OUT_DIR, sprintf("triage_metadata_report_%s.html", run_label))
write_triage_metadata_report(
  path = rpt, run_label = run_label, d = d, unl = unl, conf = conf, rates = rates,
  tests = tests, adj = adj, bal = bal, cm = cm, sd_rows = sd_rows,
  min_n = MIN_N, n_perm = N_PERM, dx_source = cm$source[1],
  generated = format(Sys.time(), "%Y-%m-%d %H:%M"))
cat(sprintf("Report: %s\n\n", rpt))

cat("Wrote to output_files/:\n")
for (f in sort(list.files(OUT_DIR, pattern = "\\.csv$")))
  cat(sprintf("  %-38s %4d rows\n", f,
              nrow(utils::read.csv(file.path(OUT_DIR, f), stringsAsFactors = FALSE))))
cat("\nDIAGNOSTIC ONLY -- nothing was dropped, flagged or changed.\n")
