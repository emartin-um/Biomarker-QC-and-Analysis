#!/usr/bin/env Rscript
# =============================================================================
# test_cohort_accounting.R  ·  2026-07-27
# -----------------------------------------------------------------------------
# Acceptance test for the per-well exit reasons and the cohort-loss accounting.
# Base R only. Exits non-zero on the first failure.
#
#   Rscript test_cohort_accounting.R
#
# WHY A TEST AND NOT JUST A DIFF: output_files/ is gitignored, so checking out
# another branch does NOT revert the CSVs. Branch-switching cannot tell you
# whether a change was safe -- only a content check can.
#
# Two halves:
#   PART A  unit checks on the helpers, against hand-built data where the right
#           answer is known by inspection. These cover the edge cases that are
#           easy to get wrong and impossible to see in an aggregate count.
#   PART B  the real run. The load-bearing check is that re-deriving the cohort
#           from merged_combined_post_QC.csv reproduces
#           filtered_combined_post_QC.csv EXACTLY. If Metadata_Merge changes a
#           drop rule and this module does not, that equality breaks -- which is
#           the whole point, since a silently drifted accounting is worse than
#           none.
# =============================================================================

HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
source(file.path(HERE, "cohort_accounting_helpers.R"))
OUT     <- file.path(HERE, "output_files")
MM_DIR  <- normalizePath(file.path(HERE, ".."), mustWork = FALSE)
PQ      <- normalizePath(file.path(HERE, "..", "..", "Primary_QC"), mustWork = FALSE)

fails <- 0L
ok <- function(label, cond, detail = "") {
  pass <- isTRUE(cond)
  if (!pass) fails <<- fails + 1L
  cat(sprintf("  [%s] %s%s\n", if (pass) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0("  -- ", detail) else ""))
}
.rd <- function(...) { p <- file.path(...); if (file.exists(p))
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) else NULL }

# =============================================================================
cat("\n=== PART A: helper unit checks ===\n")
# =============================================================================

cat("\n-- taxonomy --\n")
tax <- exit_reason_taxonomy()
ok("every reason is unique", !anyDuplicated(tax$exit_reason))
ok("every reason has a stage, class and owner",
   all(nzchar(tax$exit_stage)) && all(nzchar(tax$exit_class)) && all(nzchar(tax$owner)))
ok("exit_class is one of the four known values",
   all(tax$exit_class %in% c("qc_triage", "cohort_construction", "by_design", "retained")))
ok("qc_read_outlier is NOT an exit reason (read outliers are flagged, not dropped)",
   !("qc_read_outlier" %in% tax$exit_reason))
ok("pool_replicate is the only reason outside the denominator",
   identical(tax$exit_reason[!in_denominator(tax$exit_reason)], "pool_replicate"))
ok("QC reasons precede cohort reasons in the taxonomy order",
   max(which(tax$exit_class == "qc_triage")) < min(which(tax$exit_class == "cohort_construction")))

cat("\n-- classify_exits: first match wins, and overlaps stay visible --\n")
cl <- classify_exits(c("a","b","c","d"), list(
  qc_triage_IC    = c("a","b"),
  not_in_metadata = c("b","c")))     # b qualifies for both; IC is earlier
ok("a well in two stages is attributed to the EARLIER one",
   cl$exit_reason[cl$well_id == "b"] == "qc_triage_IC")
ok("the overlap is still counted, not hidden",
   cl$n_stages_qualified[cl$well_id == "b"] == 2)
ok("an unlisted well is in_cohort", cl$exit_reason[cl$well_id == "d"] == "in_cohort")
ok("a well listed once qualifies once", cl$n_stages_qualified[cl$well_id == "c"] == 1)
ok("duplicate ids are rejected",
   inherits(try(classify_exits(c("a","a"), list()), silent = TRUE), "try-error"))
ok("an unknown reason name is rejected",
   inherits(try(classify_exits("a", list(not_a_reason = "a")), silent = TRUE), "try-error"))

cat("\n-- split_triage_reasons / triage_primary_reason --\n")
tr <- split_triage_reasons(c("IC outlier",
                             "PCA outlier; High outlier burden (FDR); IC outlier",
                             "High outlier burden (FDR)",
                             "Bad Data"))
ok("a multi-reason string sets every axis it names",
   all(unlist(tr[2, c("qc_triage_IC","qc_triage_PCA","qc_triage_outlier")])))
ok("'High outlier burden (FDR)' is not mistaken for a PCA outlier", !tr$qc_triage_PCA[3])
ok("precedence: bad data > IC > PCA > burden",
   identical(triage_primary_reason(tr),
             c("qc_triage_IC", "qc_triage_IC", "qc_triage_outlier", "qc_triage_baddata")))

cat("\n-- mark_most_recent_visit: the three inherited edge cases --\n")
ok("a single-visit subject is kept", all(mark_most_recent_visit(c("s1"), c(70))))
ok("the highest age wins",
   identical(mark_most_recent_visit(c("s1","s1"), c(70, 75)), c(FALSE, TRUE)))
ok("an NA age ranks LAST and is therefore the row KEPT (inherited from rank())",
   identical(mark_most_recent_visit(c("s1","s1"), c(70, NA)), c(FALSE, TRUE)))
ok("on tied ages the LATER row wins (ties.method='first')",
   identical(mark_most_recent_visit(c("s1","s1"), c(70, 70)), c(FALSE, TRUE)))
ok("exactly one row survives per subject",
   sum(mark_most_recent_visit(c("a","a","a","b"), c(1, 2, 3, 9))) == 2)

cat("\n-- build_cdx_map / build_ancestry_rules --\n")
cm <- build_cdx_map(data.frame(CDX = c("AD","NCI","Other"),
                               FINAL_CDX_collapsed = c("AD","NCI",""),
                               stringsAsFactors = FALSE))
ok("a blank FINAL_CDX_collapsed excludes that CDX", !("Other" %in% names(cm)))
# Single-bracket indexing is what the runner uses, and is the reason an
# unaccounted CDX silently becomes NA rather than erroring. `[[` would throw.
ok("a CDX absent from the table maps to NA (unaccounted -> dropped)",
   is.na(unname(cm["not_listed"])))
ok("an NA CDX maps to NA too", is.na(unname(cm[NA_character_])))
ok("a mapped CDX collapses correctly", unname(cm["NCI"]) == "NCI")

ar <- build_ancestry_rules(data.frame(
  Group = c("HI", "AA"), Race = c("BL", NA), Ethnicity = c("HI", "NH"),
  FINAL_Ancestry = c("HI_BL", "AA"), stringsAsFactors = FALSE))
ok("an NA component is DROPPED from the condition, widening the rule",
   !("Race" %in% names(ar[[2]]$cond)) && length(ar[[2]]$cond) == 2)
d <- data.frame(Group = c("HI","AA","AA","ZZ"), Race = c("BL","WH","BL","WH"),
                Ethnicity = c("HI","NH","NH","NH"), stringsAsFactors = FALSE)
a <- assign_ancestry(d, ar)
ok("the widened rule matches any Race", a[2] == "AA" && a[3] == "AA")
ok("an unmatched combination becomes 'Other' (and is what missing_ancestry counts)",
   a[4] == "Other")
ok("the FIRST rule in the table wins on overlap",
   assign_ancestry(data.frame(Group = "HI", Race = "BL", Ethnicity = "HI",
                              stringsAsFactors = FALSE),
     list(list(label = "first", cond = list(Group = "HI")),
          list(label = "second", cond = list(Race = "BL")))) == "first")

cat("\n-- normalise_ancestry_cols --\n")
nz <- normalise_ancestry_cols(data.frame(Group = c(NA, " AA "), Race = c("BL ", "WH"),
                                         Ethnicity = c("NH", "HI"), stringsAsFactors = FALSE))
ok("a missing Group becomes the STRING 'NA' (how OtherProjects rows classify)",
   nz$Group[1] == "NA")
ok("trailing whitespace is trimmed (Excel round-trips break equality otherwise)",
   nz$Race[1] == "BL" && nz$Group[2] == "AA")

# =============================================================================
cat("\n\n=== PART B: the real run ===\n")
# =============================================================================

# In the source repo there is no output_files/ at all -- code and docs only, no
# data. That is not a failure; there is simply nothing to check yet.
if (!dir.exists(OUT)) {
  cat("\nNo output_files/ here -- this is the source repo, not a QC run.\n")
  cat("Part A passed. Run this inside a QC run directory after the accounting.\n")
  cat(sprintf("\n%s  (%d failure%s)\n", if (fails == 0L) "ALL CHECKS PASSED" else "FAILED",
              fails, if (fails == 1L) "" else "s"))
  quit(status = if (fails == 0L) 0L else 1L)
}

w <- .rd(OUT, "well_exit_reasons.csv")
if (is.null(w)) { ok("well_exit_reasons.csv exists", FALSE)
  cat(sprintf("\nFAILED  (%d failures)\n", fails)); quit(status = 1L) }
summ <- .rd(OUT, "cohort_loss_summary.csv")
pp   <- .rd(OUT, "per_plate_QC_and_cohort_summary.csv")

cat("\n-- the partition --\n")
ok("one row per attempted well, no duplicate ids", !anyDuplicated(w$well_id))
ok("every well carries exactly one reason, and it is in the taxonomy",
   all(!is.na(w$exit_reason)) && all(w$exit_reason %in% tax$exit_reason))
ok("the summary counts sum to the well count", sum(summ$n) == nrow(w),
   sprintf("%d vs %d", sum(summ$n), nrow(w)))
den <- sum(in_denominator(w$exit_reason))
ok("attempted = in_cohort + QC triage + cohort construction",
   den == sum(w$exit_reason == "in_cohort") + sum(w$exit_class == "qc_triage") +
          sum(w$exit_class == "cohort_construction"))
ok("exit_stage and exit_class agree with the taxonomy for every row",
   all(w$exit_stage == tax$exit_stage[match(w$exit_reason, tax$exit_reason)]) &&
   all(w$exit_class == tax$exit_class[match(w$exit_reason, tax$exit_reason)]))

cat("\n-- the load-bearing check: does the accounting match the pipeline? --\n")
mm_out <- list.dirs(MM_DIR, recursive = FALSE)
mm_out <- mm_out[grepl("(^|/)output_files", basename(mm_out)) &
                 !grepl("_archive|_backup", basename(mm_out))]
mm_out <- mm_out[file.exists(file.path(mm_out, "merged_combined_post_QC.csv"))]
MM <- if (length(mm_out)) mm_out[which.max(file.mtime(file.path(mm_out, "merged_combined_post_QC.csv")))] else NA
fc <- if (!is.na(MM)) .rd(MM, "filtered", "filtered_combined_post_QC.csv") else NULL
if (!is.null(fc)) {
  mine <- w$well_id[w$exit_reason == "in_cohort"]
  ok("in_cohort is EXACTLY filtered_combined_post_QC.csv",
     setequal(mine, as.character(fc$SAMPLE_ALIQUOT)),
     sprintf("%d vs %d; %d only-mine, %d only-pipeline", length(mine), nrow(fc),
             length(setdiff(mine, fc$SAMPLE_ALIQUOT)),
             length(setdiff(fc$SAMPLE_ALIQUOT, mine))))
  # the standard/low files are re-synced BEFORE the visit collapse, so they should
  # equal the cohort plus exactly the wells attributed to repeat_visit_collapsed
  fs <- .rd(MM, "filtered", "filtered_standard_post_QC.csv")
  if (!is.null(fs))
    ok("in_cohort + repeat_visit_collapsed == filtered_standard_post_QC.csv",
       setequal(c(mine, w$well_id[w$exit_reason == "repeat_visit_collapsed"]),
                as.character(fs$SAMPLE_ALIQUOT)),
       sprintf("%d vs %d", length(mine) + sum(w$exit_reason == "repeat_visit_collapsed"), nrow(fs)))
  mc <- .rd(MM, "merged_combined_post_QC.csv")
  if (!is.null(mc))
    ok("everything after not_in_metadata sums to merged_combined",
       sum(!w$exit_reason %in% c("pool_replicate", "not_in_metadata") &
             w$exit_class != "qc_triage") == nrow(mc),
       sprintf("%d vs %d", sum(!w$exit_reason %in% c("pool_replicate","not_in_metadata") &
                                 w$exit_class != "qc_triage"), nrow(mc)))
} else ok("Metadata_Merge outputs found for comparison", FALSE)

cat("\n-- QC side against Primary QC's own files --\n")
tri <- .rd(PQ, "output_files", "samples_to_triage.csv")
if (!is.null(tri)) {
  qc_ids <- w$well_id[w$exit_class == "qc_triage"]
  ok("the QC-triaged set is exactly samples_to_triage.csv",
     setequal(qc_ids, as.character(tri$SampleID)),
     sprintf("%d vs %d", length(qc_ids), nrow(tri)))
  ok("every triage axis flag is set only on triaged wells",
     all(!(w$flag_qc_triage_IC | w$flag_qc_triage_PCA | w$flag_qc_triage_outlier |
             w$flag_qc_triage_baddata) | w$exit_class == "qc_triage"))
  ok("every triaged well has at least one axis flag set",
     all((w$flag_qc_triage_IC | w$flag_qc_triage_PCA | w$flag_qc_triage_outlier |
            w$flag_qc_triage_baddata)[w$exit_class == "qc_triage"]))
}

fro <- .rd(PQ, "output_files", "flagged_read_outliers.csv")
if (!is.null(fro)) {
  ok("the read-outlier flag is exactly flagged_read_outliers.csv",
     setequal(w$well_id[w$flag_read_outlier], as.character(fro$SampleID)),
     sprintf("%d vs %d", sum(w$flag_read_outlier), nrow(fro)))
  # The point of demoting it from an exit reason to an annotation: most flagged
  # read outliers are IN the analysis set. If that ever stops being true, the
  # pipeline started dropping on reads and the taxonomy needs a new exit.
  inc <- sum(w$flag_read_outlier & w$exit_reason == "in_cohort")
  ok("read outliers are kept, not dropped (they are an annotation, not an exit)",
     inc > 0.5 * sum(w$flag_read_outlier),
     sprintf("%d of %d flagged read outliers are in the cohort", inc, sum(w$flag_read_outlier)))
}

opp <- .rd(PQ, "output_files", "per_plate_QC_summary.csv")
if (!is.null(opp) && !is.null(pp)) {
  m <- match(paste(opp$Run, opp$Bay), paste(pp$Run, pp$Bay))
  ok("every plate-bay in Primary QC's summary is in the accounting", !anyNA(m))
  ok("per-plate triage counts equal Primary QC's own Triaged column",
     all(opp$Triaged == pp$qc_triaged[m], na.rm = TRUE))
  ok("per-plate Retained_postQC == attempted - triaged",
     all(opp$Retained_postQC == pp$n_attempted[m] - pp$qc_triaged[m], na.rm = TRUE))
}

cat("\n-- the §1 requirement: QC and cohort never merged into one number --\n")
if (!is.null(pp)) {
  ok("the per-plate table carries QC and cohort as separate columns",
     all(c("qc_triaged", "qc_pct", "cohort_removed", "cohort_pct") %in% names(pp)))
  ok("no column merges the two stages",
     !any(grepl("^(pct_removed|total_removed|n_removed|pct_missing)$", names(pp))))
  ok("per plate-bay: attempted = in_cohort + qc_triaged + cohort_removed",
     all(pp$n_attempted == pp$n_in_cohort + pp$qc_triaged + pp$cohort_removed))
  ok("the two stages genuinely differ per bay (they are not the same number)",
     any(pp$qc_triaged != pp$cohort_removed),
     sprintf("QC total %d vs cohort total %d", sum(pp$qc_triaged), sum(pp$cohort_removed)))
}

cat("\n-- the §4 differential cuts --\n")
for (f in c("cohort_loss_by_year.csv", "cohort_loss_by_site.csv", "cohort_loss_by_ancestry.csv")) {
  b <- .rd(OUT, f)
  if (is.null(b)) { ok(paste("exists:", f), FALSE); next }
  ok(paste0(f, ": attempted = in_cohort + QC + cohort per group"),
     all(b$n_attempted == b$n_in_cohort + b$n_qc_triage + b$n_cohort_construction),
     sprintf("%d groups", nrow(b)))
  ok(paste0(f, ": QC and cohort percentages are separate columns"),
     all(c("pct_qc_triage", "pct_cohort_construction") %in% names(b)))
}
# A cut is only meaningful if wells that exited EARLY still carry the label.
# 74 of the 77 triaged wells on the 50-plate run have a metadata row; if that
# coverage collapses, the site/ancestry cuts silently become cohort-only.
qcw <- w[w$exit_class == "qc_triage", ]
if (nrow(qcw))
  ok("QC-triaged wells still carry a Site label, so the cuts see both stages",
     sum(!is.na(qcw$Site)) > 0.8 * nrow(qcw),
     sprintf("%d of %d triaged wells have a Site", sum(!is.na(qcw$Site)), nrow(qcw)))

cat("\n-- sanity on the headline --\n")
qc <- sum(w$exit_class == "qc_triage"); coh <- sum(w$exit_class == "cohort_construction")
ok("cohort construction is the larger share of the loss, as the investigation found",
   coh > qc, sprintf("QC %d (%.1f%%) vs cohort %d (%.1f%%) of %d lost",
                     qc, 100*qc/(qc+coh), coh, 100*coh/(qc+coh), qc + coh))
ok("missing_diagnosis is a real category, not a rounding error",
   sum(w$exit_reason == "missing_diagnosis") > 0,
   sprintf("%d wells", sum(w$exit_reason == "missing_diagnosis")))

cat(sprintf("\n%s  (%d failure%s)\n",
            if (fails == 0L) "ALL CHECKS PASSED" else "FAILED",
            fails, if (fails == 1L) "" else "s"))
quit(status = if (fails == 0L) 0L else 1L)
