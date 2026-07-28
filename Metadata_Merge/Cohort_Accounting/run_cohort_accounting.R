#!/usr/bin/env Rscript
# =============================================================================
# run_cohort_accounting.R  ·  per-well exit reasons + the standing cohort-loss
#                             accounting
# -----------------------------------------------------------------------------
# Implements QC_PIPELINE_RECOMMENDATIONS_2026-07-13.md §1 and §4 — the last two
# open items from the 2026 batch investigation.
#
# BASE R ONLY — no dplyr — so it runs anywhere the QC runs do.
# Run from the module directory:  Rscript run_cohort_accounting.R
#
# Reads (from the run directory two levels up):
#   Primary_QC/input_files/Sample_QC.csv         the attempted-well universe
#   Primary_QC/input_files/Replicate_list.csv    the HIHG pool wells
#   Primary_QC/output_files/samples_to_triage.csv        who QC removed, and why
#   Primary_QC/output_files/flagged_read_outliers.csv    flagged, NOT removed
#   Primary_QC/output_files/per_plate_QC_summary.csv     for the side-by-side check
#   Metadata_Merge/<output_dir>/merged_combined_post_QC.csv
#   Metadata_Merge/<output_dir>/sample_exclusion_report.csv
#   Metadata_Merge/<output_dir>/filtered/filtered_combined_post_QC.csv
#   Metadata_Merge/review/{CDX,Ancestry}_review.csv      the drop configs
#
# Writes to output_files/:
#   well_exit_reasons.csv             one row per attempted well, one reason each
#   cohort_loss_summary.csv           attempted -> in-cohort, by reason, with owner
#   cohort_loss_by_year.csv           the §4 differential cuts
#   cohort_loss_by_site.csv
#   cohort_loss_by_ancestry.csv
#   per_plate_QC_and_cohort_summary.csv   the §1 side-by-side table
#
# 🚫 NO DROP RULE HERE. This module removes nothing and changes no output. It
# only says, for each well that is already absent, WHICH STAGE removed it.
# =============================================================================

HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
MM_DIR   <- normalizePath(file.path(HERE, ".."), mustWork = FALSE)
RUN_ROOT <- normalizePath(file.path(HERE, "..", ".."), mustWork = FALSE)
OUT_DIR  <- file.path(HERE, "output_files")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(HERE, "cohort_accounting_helpers.R"))

.rd <- function(...) {
  p <- file.path(...)
  if (!file.exists(p)) return(NULL)
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
}
.need <- function(x, what) if (is.null(x)) stop("cohort accounting needs ", what, " and it is missing") else x

PQ <- file.path(RUN_ROOT, "Primary_QC")

## --- locate the Metadata_Merge output directory ------------------------------
# The scaffold names it output_files_<run tag>, so find it rather than hard-code.
mm_out <- list.dirs(MM_DIR, recursive = FALSE)
mm_out <- mm_out[grepl("(^|/)output_files", basename(mm_out)) &
                 !grepl("_archive|_backup", basename(mm_out))]
mm_out <- mm_out[file.exists(file.path(mm_out, "merged_combined_post_QC.csv"))]
if (!length(mm_out)) stop("no Metadata_Merge output dir with merged_combined_post_QC.csv under ", MM_DIR)
MM <- mm_out[which.max(file.mtime(file.path(mm_out, "merged_combined_post_QC.csv")))]
cat("Metadata_Merge outputs: ", basename(MM), "\n", sep = "")

## --- 1. the attempted-well universe ------------------------------------------
# Sample_QC.csv is Alamar's own per-well table and the ONLY file that still lists
# a well after QC removed it, so it is the authoritative denominator. Its first
# 12 lines are a metric glossary, hence skip = 12.
q <- .need(tryCatch(utils::read.csv(file.path(PQ, "input_files", "Sample_QC.csv"),
                                    skip = 12, check.names = FALSE, stringsAsFactors = FALSE),
                    error = function(e) NULL),
           "Primary_QC/input_files/Sample_QC.csv")
w <- q[q[["Sample Type"]] == "Sample", , drop = FALSE]
w$well_id <- as.character(w[["Sample Name"]])
w$Run <- sub("^([0-9]{8}-[0-9]{4}).*", "\\1", w$plateID)
w$Bay <- sub("^.*_(Bay[0-9]+)_.*$", "\\1", w$plateID)
w$year <- substr(w$Run, 1, 4)
w$well <- w$AUTO_WELLPOSITION
cat(sprintf("Attempted wells: %d of %d rows in Sample_QC.csv (%d are NC/SC/IPC controls)\n",
            nrow(w), nrow(q), nrow(q) - nrow(w)))
if (anyDuplicated(w$well_id)) stop("Sample_QC.csv has duplicate Sample Names -- cannot key on them")

## --- 2. what each stage removed ----------------------------------------------
rl  <- .rd(PQ, "input_files", "Replicate_list.csv")
pool_ids <- if (is.null(rl)) character(0) else as.character(rl[[1]])
pool_ids <- intersect(pool_ids, w$well_id)

tri <- .need(.rd(PQ, "output_files", "samples_to_triage.csv"), "samples_to_triage.csv")
tflags <- split_triage_reasons(tri$Reason)
tri$primary_reason <- triage_primary_reason(tflags)
if (any(is.na(tri$primary_reason)))
  stop("samples_to_triage.csv has a Reason this module does not recognise: ",
       paste(unique(tri$Reason[is.na(tri$primary_reason)]), collapse = " | "),
       "\n  -> add it to split_triage_reasons() and the taxonomy; do not let it fall through.")

fro <- .rd(PQ, "output_files", "flagged_read_outliers.csv")
read_outlier_ids <- if (is.null(fro)) character(0) else as.character(fro$SampleID)

mc <- .need(.rd(MM, "merged_combined_post_QC.csv"), "merged_combined_post_QC.csv")
ex <- .need(.rd(MM, "sample_exclusion_report.csv"), "sample_exclusion_report.csv")
fc <- .need(.rd(MM, "filtered", "filtered_combined_post_QC.csv"), "filtered_combined_post_QC.csv")

# Everything that reached Primary QC's post-QC NPQ but has no metadata row.
post_qc_ids <- setdiff(w$well_id, c(pool_ids, as.character(tri$SampleID)))
not_in_meta <- setdiff(post_qc_ids, as.character(mc$SAMPLE_ALIQUOT))

# The exclusion report: three channels, distinguished by the reason text the
# pipeline writes. Anything excluded that is neither of the two automatic rules
# is a human decision -- an analyst hand-editing `exclude` and re-rendering with
# filter_only = TRUE, which is a documented workflow.
exr <- ex[ex$exclude %in% c(TRUE, "TRUE", "true"), , drop = FALSE]
is_dup  <- grepl("^Duplicate SAMPLE", exr$exclude_reason)
is_rng  <- grepl("^Implausible ",     exr$exclude_reason) & !is_dup
dup_ids <- as.character(exr$SAMPLE_ALIQUOT[is_dup])
rng_ids <- as.character(exr$SAMPLE_ALIQUOT[is_rng])
rev_ids <- as.character(exr$SAMPLE_ALIQUOT[!is_dup & !is_rng])

## --- 3. re-derive the covariate drops ----------------------------------------
d <- mc[!mc$SAMPLE_ALIQUOT %in% c(dup_ids, rng_ids, rev_ids), , drop = FALSE]
d <- normalise_ancestry_cols(d)

cdxrev <- .rd(MM_DIR, "review", "CDX_review.csv")
ancrev <- .rd(MM_DIR, "review", "Ancestry_review.csv")
if (is.null(cdxrev) || is.null(ancrev))
  stop("review/CDX_review.csv and review/Ancestry_review.csv are required.\n",
       "  Without them the pipeline falls back to the built-in config in\n",
       "  Metadata_Merge_Pipeline.Rmd Step 6b, which this module cannot read.\n",
       "  Re-render Metadata_Merge to emit them, or point at the run that has them.")

cdx_map <- build_cdx_map(cdxrev)
d$CDX_collapsed <- unname(cdx_map[as.character(d$CDX)])
cdx_ids <- as.character(d$SAMPLE_ALIQUOT[is.na(d$CDX_collapsed)])
unaccounted <- setdiff(unique(as.character(d$CDX[is.na(d$CDX_collapsed)])),
                       c(as.character(cdxrev$CDX), NA))
d <- d[!is.na(d$CDX_collapsed), , drop = FALSE]

anc_rules <- build_ancestry_rules(ancrev)
d$Ancestry <- assign_ancestry(d, anc_rules)
anc_ids <- as.character(d$SAMPLE_ALIQUOT[d$Ancestry == "Other"])
d <- d[d$Ancestry != "Other", , drop = FALSE]

n_na_rid <- sum(is.na(d$Record_ID) | trimws(as.character(d$Record_ID)) == "")
keep_v <- mark_most_recent_visit(d$Record_ID, d$age_at_subject)
rep_ids <- as.character(d$SAMPLE_ALIQUOT[!keep_v])

## --- 4. classify ---------------------------------------------------------------
# Order IS the pipeline order. First match wins.
exits <- list(
  pool_replicate          = pool_ids,
  qc_triage_baddata       = tri$SampleID[tri$primary_reason == "qc_triage_baddata"],
  qc_triage_IC            = tri$SampleID[tri$primary_reason == "qc_triage_IC"],
  qc_triage_PCA           = tri$SampleID[tri$primary_reason == "qc_triage_PCA"],
  qc_triage_outlier       = tri$SampleID[tri$primary_reason == "qc_triage_outlier"],
  not_in_metadata         = not_in_meta,
  metadata_review         = rev_ids,
  duplicate_aliquot_rerun = dup_ids,
  covariate_out_of_range  = rng_ids,
  missing_diagnosis       = cdx_ids,
  missing_ancestry        = anc_ids,
  repeat_visit_collapsed  = rep_ids)

cl <- classify_exits(w$well_id, exits)
w$exit_reason <- cl$exit_reason
tax <- exit_reason_taxonomy()
w$exit_stage <- tax$exit_stage[match(w$exit_reason, tax$exit_reason)]
w$exit_class <- tax$exit_class[match(w$exit_reason, tax$exit_reason)]

## --- 5. the independent axes, never folded into the reason -------------------
mi <- match(w$well_id, as.character(tri$SampleID))
for (v in names(tflags)) w[[paste0("flag_", v)]] <- ifelse(is.na(mi), FALSE, tflags[[v]][mi])
w$flag_read_outlier   <- w$well_id %in% read_outlier_ids
w$flag_alamar_warning <- w[["QC Status"]] == "Warning"
w$n_stages_qualified  <- cl$n_stages_qualified

## --- 6. labels for the cuts ---------------------------------------------------
# Deliberately taken from the METADATA, not from the surviving cohort, so a well
# that exited at QC triage still carries its Site and ancestry. That is what
# makes a differential-attrition cut possible at all: 74 of the 77 triaged wells
# on the 50-plate run have a metadata row.
md <- .rd(MM, "metadata_PrimaryQC_refreshed.csv")
if (is.null(md)) md <- mc
mm <- match(w$well_id, as.character(md$SAMPLE_ALIQUOT))
w$Site            <- if ("Site" %in% names(md)) md$Site[mm] else NA_character_
w$metadata_source <- if ("metadata_source" %in% names(md)) md$metadata_source[mm] else NA_character_
mdn <- normalise_ancestry_cols(md[, intersect(c("Group","Race","Ethnicity"), names(md)), drop = FALSE])
w$Ancestry <- if (ncol(mdn) == 3) assign_ancestry(mdn, anc_rules)[mm] else NA_character_
w$CDX_collapsed <- if ("CDX" %in% names(md)) unname(cdx_map[as.character(md$CDX)])[mm] else NA_character_

## --- 7. write -----------------------------------------------------------------
keepc <- c("well_id", "Run", "Bay", "well", "year", "plateID",
           "exit_reason", "exit_stage", "exit_class",
           "Site", "Ancestry", "CDX_collapsed", "metadata_source",
           "flag_qc_triage_baddata", "flag_qc_triage_IC", "flag_qc_triage_PCA",
           "flag_qc_triage_outlier", "flag_read_outlier", "flag_alamar_warning",
           "n_stages_qualified")
wells <- w[order(w$Run, w$Bay, w$well), intersect(keepc, names(w))]

wr <- function(x, f) if (!is.null(x) && nrow(x)) {
  utils::write.csv(x, file.path(OUT_DIR, f), row.names = FALSE)
  cat(sprintf("  %-38s %5d rows\n", f, nrow(x)))
}
cat("\nWrote:\n")
wr(wells, "well_exit_reasons.csv")
summ <- cohort_loss_table(w$exit_reason);                  wr(summ, "cohort_loss_summary.csv")

# A cut can only see wells that carry its label, and the label comes from the
# metadata -- so a well that exited BECAUSE it has no metadata row is invisible
# to the by-site and by-ancestry cuts by construction. Say how many rather than
# dropping them quietly: an unreported denominator is how a differential gets
# missed in the first place.
cut_of <- function(lab, name, min_n = 1) {
  n_unlab <- sum(is.na(lab) & in_denominator(w$exit_reason))
  if (n_unlab)
    cat(sprintf("    (%s: %d attempted wells have no %s and are outside this cut -- %s)\n",
                name, n_unlab, name,
                paste(sprintf("%s %d", names(table(w$exit_reason[is.na(lab) &
                        in_denominator(w$exit_reason)])),
                      table(w$exit_reason[is.na(lab) & in_denominator(w$exit_reason)])),
                      collapse = ", ")))
  cohort_loss_by(w$exit_reason, lab, name, min_n = min_n)
}
wr(cut_of(w$year,     "year"),            "cohort_loss_by_year.csv")
wr(cut_of(w$Site,     "Site",     5),     "cohort_loss_by_site.csv")
wr(cut_of(w$Ancestry, "Ancestry", 5),     "cohort_loss_by_ancestry.csv")
pp <- per_plate_accounting(w$exit_reason, w$Run, w$Bay, w$flag_read_outlier, w$flag_alamar_warning)
wr(pp, "per_plate_QC_and_cohort_summary.csv")

## --- 8. reconcile against the pipeline's own outputs -------------------------
# The whole module is worthless if it disagrees with the files it describes, so
# check rather than assert. These are the same equalities the acceptance test
# re-checks; failing them here surfaces the problem at build time.
cat("\n================ RECONCILIATION ================\n")
mine <- w$well_id[w$exit_reason == "in_cohort"]
chk <- function(lbl, cond, detail = "")
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(cond)) "ok" else "MISMATCH", lbl,
              if (nzchar(detail)) paste0("  -- ", detail) else ""))
chk("in_cohort set == filtered_combined_post_QC.csv",
    setequal(mine, as.character(fc$SAMPLE_ALIQUOT)),
    sprintf("%d vs %d", length(mine), nrow(fc)))
chk("post-QC count == NPQ post_QC row count",
    length(post_qc_ids) == nrow(mc) + length(not_in_meta),
    sprintf("%d = %d + %d", length(post_qc_ids), nrow(mc), length(not_in_meta)))
chk("every attempted well has exactly one reason",
    !any(is.na(w$exit_reason)) && nrow(w) == sum(summ$n))

opp <- .rd(PQ, "output_files", "per_plate_QC_summary.csv")
if (!is.null(opp)) {
  m2 <- match(paste(opp$Run, opp$Bay), paste(pp$Run, pp$Bay))
  chk("per-plate triage counts match Primary QC's own summary",
      all(opp$Triaged == pp$qc_triaged[m2], na.rm = TRUE))
}

## --- 9. console summary --------------------------------------------------------
den <- sum(in_denominator(w$exit_reason))
qc  <- sum(w$exit_class == "qc_triage")
coh <- sum(w$exit_class == "cohort_construction")
inc <- sum(w$exit_reason == "in_cohort")

cat("\n================ COHORT-LOSS ACCOUNTING ================\n")
cat(sprintf("\nAttempted patient wells: %d   (plus %d HIHG pool wells, outside the denominator)\n",
            den, sum(w$exit_reason == "pool_replicate")))
cat(sprintf("In cohort:               %d  (%.1f%%)\n", inc, 100 * inc / den))
cat(sprintf("Lost:                    %d  (%.1f%%)\n\n", den - inc, 100 * (den - inc) / den))
cat(sprintf("  QC triage           %5d   %5.1f%% of the loss   (wet lab)\n",
            qc, 100 * qc / (den - inc)))
cat(sprintf("  Cohort construction %5d   %5.1f%% of the loss   (data entry / study design)\n",
            coh, 100 * coh / (den - inc)))
cat("\n-- by reason --\n")
s <- summ[summ$n > 0 & summ$exit_reason != "in_cohort", ]
for (i in order(-s$n))
  cat(sprintf("  %-24s %5d  %5.2f%%   %s\n", s$exit_reason[i], s$n[i],
              s$pct_of_attempted[i], s$owner[i]))

if (length(unaccounted))
  cat(sprintf("\n⚠ CDX value(s) absent from CDX_review.csv entirely, dropped with no warning\n  in the rendered report: %s\n",
              paste(unaccounted, collapse = ", ")))
if (n_na_rid > 0)
  cat(sprintf("\n⚠ %d rows reached the visit collapse with a missing Record_ID; they pool into\n  one group and only one survives. Check before trusting the repeat-visit line.\n",
              n_na_rid))
ov <- sum(w$n_stages_qualified > 1)
if (ov > 0)
  cat(sprintf("\nNote: %d wells qualified for more than one exit; each is attributed to the\n  first stage that removed it. See n_stages_qualified in well_exit_reasons.csv.\n", ov))

cat("\n-- by year (QC and cohort kept apart) --\n")
by <- cohort_loss_by(w$exit_reason, w$year, "year")
if (!is.null(by)) for (i in seq_len(nrow(by)))
  cat(sprintf("  %-6s attempted %4d -> cohort %4d   QC %4.1f%%   cohort-construction %4.1f%%\n",
              by$group[i], by$n_attempted[i], by$n_in_cohort[i],
              by$pct_qc_triage[i], by$pct_cohort_construction[i]))

cat("\n🚫 NOTHING HERE IS A DROP RULE. Every well above was already in or out; this\n")
cat("   only records WHICH STAGE removed it. Read the QC and cohort columns\n")
cat("   separately -- a bay that is 40% missing on diagnoses is a data-entry\n")
cat("   problem, one that is 40% missing on IC failures is a wet-lab problem.\n")
