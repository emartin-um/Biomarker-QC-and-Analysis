# =============================================================================
# cohort_accounting_helpers.R  ·  added 2026-07-27
# -----------------------------------------------------------------------------
# Per-well exit reasons and the standing cohort-loss accounting, from
#   Batch_Effects/2026_Batch_Investigation/QC_PIPELINE_RECOMMENDATIONS_2026-07-13.md
# §1 (separate QC-triaged from cohort-filtered) and §4 (a standing accounting).
# Base R only, I/O-free — so the headless runner, the test and any future Rmd all
# get the SAME definition of every drop rule. There must be exactly one.
#
# THE PROBLEM THIS SOLVES. A well missing from the analysis set was repeatedly
# read as a well that "failed QC". They are different things with different
# owners, and on the 50-plate run they are not even the same order of magnitude:
# of 648 attempted patient wells that never reach the cohort, 77 (11.9%) are QC
# triage and 571 (88.1%) are cohort/metadata construction — 361 of those a
# missing diagnosis. A bay that is "40% missing" from missing diagnoses is a
# data-entry issue; one that is "40% missing" from IC failures is a wet-lab
# issue. Same number, opposite owner. Publishing one merged "% removed" per bay
# is what set off the entire 2026 batch investigation.
#
# WHAT THE TAXONOMY IS NOT. It is not the list in the 2026-07-13 recommendation.
# That list was written from the outside; this one was read off the code, and
# three things differ (all verified against the 50-plate frozen run):
#
#   * `qc_read_outlier` DOES NOT EXIST as an exit. Read outliers are FLAGGED and
#     kept — 94 of them, of which 89 are in the final data and 5 exit on an
#     unrelated triage reason. It is carried here as the annotation
#     `flag_read_outlier`, never as an exit reason.
#   * `qc_triage_baddata` had to be ADDED. Primary QC's triage table has a
#     fourth reason, "Bad Data" (non-numeric rows), that the recommendation did
#     not know about. It was 0 on this run, which is exactly why it is easy to
#     miss.
#   * `missing_covariate` SPLITS IN TWO, because the pipeline has no single
#     covariate filter. It drops on an unmatched ancestry combination
#     (`missing_ancestry`) and on an out-of-range age/BMI
#     (`covariate_out_of_range`). It never drops on a missing age, sex, APOE or
#     education — those reach the final dataset and are only reported.
#
# And `repeat_visit_collapsed` is split from `duplicate_aliquot_rerun`: the
# pipeline collapses a re-assayed specimen (most-recent RUN) at a different step
# from a genuine second visit (highest age per Record_ID), and they have
# different owners. This is the recommendation's own §5 — flag every axis
# independently — applied to its own §1 list.
# =============================================================================

## --- the taxonomy ------------------------------------------------------------

#' The exit-reason taxonomy, in the order the pipeline applies it.
#'
#' Every attempted patient well gets EXACTLY ONE of these. Order is load-bearing:
#' a well that would qualify for several exits is attributed to the first one it
#' reaches, because that is the one that actually removed it. Independent axes
#' (a triaged well that is also a read outlier, say) are carried as separate
#' flag columns, never folded into the reason.
#'
#' `owner` is the point of the whole module: it says who fixes it.
exit_reason_taxonomy <- function() {
  data.frame(
    exit_reason = c(
      "pool_replicate",
      "qc_triage_baddata", "qc_triage_IC", "qc_triage_PCA", "qc_triage_outlier",
      "not_in_metadata", "metadata_review", "duplicate_aliquot_rerun",
      "covariate_out_of_range", "missing_diagnosis", "missing_ancestry",
      "repeat_visit_collapsed",
      "in_cohort"),
    exit_stage = c(
      "Primary_QC",
      rep("Primary_QC", 4),
      rep("Metadata_Merge", 7),
      "-"),
    exit_class = c(
      "by_design",
      rep("qc_triage", 4),
      rep("cohort_construction", 7),
      "retained"),
    owner = c(
      "by design - HIHG pool well, never a patient",
      "wet-lab - non-numeric assay output",
      "wet-lab - internal-control failure",
      "wet-lab - multivariate profile outlier",
      "wet-lab - per-marker outlier burden",
      "data-entry - no metadata row for this aliquot",
      "data-entry - manual exclusion on review",
      "by design - specimen assayed more than once",
      "data-entry - age or BMI outside a plausible range",
      "data-entry - diagnosis missing or not codeable",
      "data-entry - race/ethnicity combination has no ancestry rule",
      "by design - subject seen more than once, collapsed to one visit",
      "-"),
    stringsAsFactors = FALSE)
}

#' Wells counted in the loss denominator.
#'
#' `pool_replicate` is excluded: the HIHG pool is a fixed control specimen on
#' every bay and was never a candidate for the cohort, so counting it as "lost"
#' would dilute the QC-vs-cohort ratio this module exists to make legible. It is
#' still tagged and still reported on its own line — it is out of the ratio, not
#' out of the accounting.
in_denominator <- function(exit_reason) exit_reason != "pool_replicate"

## --- the classifier ----------------------------------------------------------

#' Attribute each well to the first stage that removed it.
#'
#' @param ids   character vector of attempted well ids (one row per well).
#' @param exits named list, IN PIPELINE ORDER: each element is the vector of ids
#'   that exit at that reason. Names must be exit reasons from the taxonomy.
#'   Membership may overlap between stages — that is expected, and first match
#'   wins, because the earlier stage is the one that actually removed the well.
#' @return a data.frame with the id, its single `exit_reason`, and
#'   `n_stages_qualified` — how many stages the well would have qualified for.
#'   That last column is the honesty check: it is how you see an overlap instead
#'   of having the taxonomy hide it.
classify_exits <- function(ids, exits) {
  ids <- as.character(ids)
  if (anyDuplicated(ids))
    stop("classify_exits(): duplicate well ids -- the accounting must be one row per well")
  known <- exit_reason_taxonomy()$exit_reason
  bad <- setdiff(names(exits), known)
  if (length(bad))
    stop("classify_exits(): unknown exit reason(s): ", paste(bad, collapse = ", "))

  reason <- rep(NA_character_, length(ids))
  n_qual <- integer(length(ids))
  for (r in names(exits)) {
    hit <- ids %in% as.character(exits[[r]])
    n_qual <- n_qual + hit
    reason[hit & is.na(reason)] <- r
  }
  reason[is.na(reason)] <- "in_cohort"
  data.frame(well_id = ids, exit_reason = reason, n_stages_qualified = n_qual,
             stringsAsFactors = FALSE)
}

#' Split Primary QC's collapsed triage string into the four independent axes.
#'
#' `samples_to_triage.csv` carries one row per triaged well with the reasons
#' pasted together ("PCA outlier; High outlier burden (FDR); IC outlier"), so a
#' well can trip several at once — 16 of the 77 do on the 50-plate run. The
#' single exit reason has to pick one; these booleans keep all of them.
#'
#' @param reason the `Reason` column of samples_to_triage.csv.
#' @return data.frame of logical columns, one per triage axis.
split_triage_reasons <- function(reason) {
  r <- as.character(reason); r[is.na(r)] <- ""
  data.frame(
    qc_triage_baddata = grepl("Bad Data",                  r, fixed = TRUE),
    qc_triage_IC      = grepl("IC outlier",                r, fixed = TRUE),
    qc_triage_PCA     = grepl("PCA outlier",               r, fixed = TRUE),
    qc_triage_outlier = grepl("High outlier burden (FDR)", r, fixed = TRUE),
    stringsAsFactors = FALSE)
}

#' Which single triage reason a well is attributed to, when it trips several.
#'
#' Precedence runs most-fundamental first: bad data (the assay output could not
#' be read at all) before an internal-control failure, before a multivariate
#' profile outlier, before a per-marker burden count. The full set stays in the
#' four boolean columns; this only decides which line of the accounting the well
#' lands on.
triage_primary_reason <- function(tri) {
  ifelse(tri$qc_triage_baddata, "qc_triage_baddata",
  ifelse(tri$qc_triage_IC,      "qc_triage_IC",
  ifelse(tri$qc_triage_PCA,     "qc_triage_PCA",
  ifelse(tri$qc_triage_outlier, "qc_triage_outlier", NA_character_))))
}

## --- the Metadata_Merge drop rules, re-derived -------------------------------
# These reimplement, in base R, exactly what Metadata_Merge_Pipeline.Rmd does via
# covariate_explorer.R's apply_category_config() / apply_grouping_config() and
# the Step-8 longitudinal collapse. Verified 2026-07-27 on the 50-plate frozen
# run: applying them to merged_combined_post_QC.csv reproduces
# filtered_combined_post_QC.csv exactly -- same 3552 SAMPLE_ALIQUOTs, same
# CDX_collapsed, same Ancestry. Re-run test_cohort_accounting.R after any change
# to the pipeline; that equality IS the test.

#' Normalise the ancestry rule columns exactly as the pipeline's step6-normalize.
#'
#' Trailing whitespace from Excel round-trips silently breaks equality on
#' multi-word values like "BL WH", and a missing Group is matched as the STRING
#' "NA" (that is how OtherProjects rows, which have no Group, get classified on
#' Race + Ethnicity alone).
normalise_ancestry_cols <- function(d) {
  for (v in c("Group", "Race", "Ethnicity"))
    if (v %in% names(d)) d[[v]] <- trimws(as.character(d[[v]]))
  if ("Group" %in% names(d)) d$Group[is.na(d$Group)] <- "NA"
  d
}

#' CDX -> CDX_collapsed map from the review table. Blank FINAL = exclude.
#'
#' A CDX value that is absent from the table entirely is "unaccounted" and is
#' dropped too, by the same is.na() test. The pipeline warns about those, but the
#' chunk sets `warning: false`, so in the rendered report a brand-new CDX
#' spelling drops its samples with no visible signal at all. Here they are simply
#' `missing_diagnosis` like any other unmapped value, and named in the report.
build_cdx_map <- function(cdx_review) {
  f <- trimws(as.character(cdx_review$FINAL_CDX_collapsed))
  f[f == ""] <- NA
  keep <- !is.na(f) & !is.na(cdx_review$CDX)
  stats::setNames(f[keep], as.character(cdx_review$CDX)[keep])
}

#' Ancestry rules from the review table, in the pipeline's own construction.
#'
#' One rule per row with a non-blank FINAL_Ancestry. A component that is NA or
#' blank is DROPPED from the condition rather than matched — which makes that
#' rule broader, not narrower (a row with Ethnicity blank matches any
#' Ethnicity). Easy to misread, so it is spelled out here.
build_ancestry_rules <- function(anc_review) {
  f <- trimws(as.character(anc_review$FINAL_Ancestry))
  f[f == ""] <- NA
  keep <- which(!is.na(f))
  lapply(keep, function(i) {
    cond <- list(Group     = anc_review$Group[i],
                 Race      = anc_review$Race[i],
                 Ethnicity = anc_review$Ethnicity[i])
    cond <- cond[!vapply(cond, function(v) is.na(v) || v == "", logical(1))]
    list(label = f[i], cond = lapply(cond, as.character))
  })
}

#' Apply the ancestry rules. Unmatched rows get `other_label` and are dropped.
#'
#' Rules are applied in REVERSE order so the FIRST rule in the table wins on an
#' overlap — the pipeline's convention, preserved.
assign_ancestry <- function(d, rules, other_label = "Other") {
  out <- rep(other_label, nrow(d))
  for (i in rev(seq_along(rules))) {
    rl <- rules[[i]]
    mask <- rep(TRUE, nrow(d))
    for (cn in names(rl$cond)) {
      cm <- d[[cn]] == rl$cond[[cn]]
      cm[is.na(cm)] <- FALSE
      mask <- mask & cm
    }
    out[mask] <- rl$label
  }
  out
}

#' The Step-8 longitudinal collapse: one row per subject, the most recent visit.
#'
#' "Most recent" is the highest `age_at_subject` within a Record_ID. Three edge
#' cases are inherited deliberately, because the accounting has to match what the
#' pipeline does rather than what it ought to do:
#'   * ties.method = "first" -- among equal ages the LATER row in data order wins;
#'   * rank()'s na.last = TRUE -- a visit with NA age ranks last and is therefore
#'     the one KEPT, so a subject with one real age and one NA age keeps the NA;
#'   * a missing Record_ID pools every such row into a single group, of which one
#'     row survives. The runner counts those separately and says so.
#'
#' @return logical vector, TRUE for the row kept.
mark_most_recent_visit <- function(record_id, age) {
  keep <- logical(length(record_id))
  for (ix in split(seq_along(record_id), record_id)) {
    r <- rank(age[ix], ties.method = "first")
    keep[ix[r == length(ix)]] <- TRUE
  }
  keep
}

## --- the accounting tables ---------------------------------------------------

#' attempted -> in-cohort, one row per exit reason, with the owner attached.
cohort_loss_table <- function(exit_reason) {
  tax <- exit_reason_taxonomy()
  n <- vapply(tax$exit_reason, function(r) sum(exit_reason == r), integer(1))
  out <- data.frame(tax, n = unname(n), stringsAsFactors = FALSE)
  denom <- sum(exit_reason[in_denominator(exit_reason)] != "")
  out$pct_of_attempted <- round(100 * out$n / denom, 2)
  out$pct_of_attempted[!in_denominator(out$exit_reason)] <- NA_real_
  out[order(match(out$exit_reason, tax$exit_reason)), ]
}

#' The same accounting cut by a grouping variable (year, site, ancestry, bay).
#'
#' Reports the QC and the cohort loss as SEPARATE columns and never their sum.
#' A single "% removed" per group is the exact thing §1 forbids.
cohort_loss_by <- function(exit_reason, group, group_var = "group", min_n = 1) {
  tax <- exit_reason_taxonomy()
  cls <- tax$exit_class[match(exit_reason, tax$exit_reason)]
  keep <- in_denominator(exit_reason) & !is.na(group)
  g <- as.character(group)[keep]; cl <- cls[keep]; er <- exit_reason[keep]
  if (!length(g)) return(NULL)

  lv <- names(which(table(g) >= min_n))
  do.call(rbind, lapply(lv, function(k) {
    i <- g == k
    att <- sum(i)
    qc  <- sum(i & cl == "qc_triage")
    coh <- sum(i & cl == "cohort_construction")
    data.frame(
      group_var = group_var, group = k,
      n_attempted = att,
      n_in_cohort = sum(i & er == "in_cohort"),
      n_qc_triage = qc,
      n_cohort_construction = coh,
      pct_qc_triage = round(100 * qc / att, 2),
      pct_cohort_construction = round(100 * coh / att, 2),
      n_missing_diagnosis = sum(i & er == "missing_diagnosis"),
      n_not_in_metadata   = sum(i & er == "not_in_metadata"),
      stringsAsFactors = FALSE, row.names = NULL)
  }))
}

#' Per plate-bay, the QC columns and the cohort columns SIDE BY SIDE.
#'
#' This is the §1 deliverable. `per_plate_QC_summary.csv` from Primary QC counts
#' only triage and flags; it stops at the QC stage, so reading it as "what
#' happened to this bay" silently attributes the cohort loss to the wet lab.
#' Here both stages appear as their own columns and are never summed.
per_plate_accounting <- function(exit_reason, run, bay, flag_read_outlier = NULL,
                                 flag_alamar_warning = NULL) {
  tax <- exit_reason_taxonomy()
  cls <- tax$exit_class[match(exit_reason, tax$exit_reason)]
  keep <- in_denominator(exit_reason)
  key <- paste(run, bay, sep = "\r")[keep]
  cl <- cls[keep]; er <- exit_reason[keep]
  fro <- if (is.null(flag_read_outlier)) rep(FALSE, length(exit_reason)) else flag_read_outlier
  faw <- if (is.null(flag_alamar_warning)) rep(FALSE, length(exit_reason)) else flag_alamar_warning
  fro <- fro[keep]; faw <- faw[keep]

  lv <- sort(unique(key))
  out <- do.call(rbind, lapply(lv, function(k) {
    i <- key == k
    att <- sum(i)
    data.frame(
      Run = sub("\r.*$", "", k), Bay = sub("^.*\r", "", k),
      n_attempted            = att,
      # --- QC stage (Primary_QC) ---
      qc_triaged             = sum(i & cl == "qc_triage"),
      qc_pct                 = round(100 * sum(i & cl == "qc_triage") / att, 2),
      qc_flag_read_outlier   = sum(i & fro),
      qc_flag_alamar_warning = sum(i & faw),
      # --- cohort stage (Metadata_Merge) ---
      cohort_removed         = sum(i & cl == "cohort_construction"),
      cohort_pct             = round(100 * sum(i & cl == "cohort_construction") / att, 2),
      cohort_missing_diagnosis = sum(i & er == "missing_diagnosis"),
      cohort_not_in_metadata   = sum(i & er == "not_in_metadata"),
      cohort_repeat_or_rerun   = sum(i & er %in% c("repeat_visit_collapsed",
                                                   "duplicate_aliquot_rerun")),
      # --- the end of the line ---
      n_in_cohort            = sum(i & er == "in_cohort"),
      stringsAsFactors = FALSE, row.names = NULL)
  }))
  out[order(out$Run, out$Bay), ]
}
