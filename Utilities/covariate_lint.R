# =============================================================================
# covariate_lint.R  ·  added 2026-07-27
# -----------------------------------------------------------------------------
# Implements 2026-07-13 recommendation §7, guard 1: a covariate-name collision
# lint. Base R, no dependencies.
#
# THE BUG THIS PREVENTS, which we actually hit. A derived covariate named
# `APOE4_count` collided with the **APOE4 protein** biomarker column. Everything
# ran. No error, no warning — APOE4's true AD signal simply stayed masked until
# someone noticed by eye. That is the worst class of bug in this pipeline: it
# does not fail, it answers.
#
# THREE COLLISION KINDS, and the mechanism for each:
#
#   EXACT      covariate name == biomarker name.
#              Whichever is written second wins. Unambiguous bug.
#
#   CASE-ONLY  `cst3` vs `CST3`. R is case-sensitive, humans are not; this is a
#              bug waiting to be written.
#
#   PREFIX     `APOE4_count` begins with the biomarker `APOE4`. THIS IS THE ONE
#              THAT BIT US, and the mechanism is R's partial matching: on a data
#              frame, `d$APOE4` silently returns `APOE4_count` whenever the exact
#              column is absent from the frame being subscripted. It is also what
#              makes `grep("APOE4", names(d))` return both columns.
#              Verified: `data.frame(APOE4_count = 1:3)$APOE4` returns 1 2 3.
#
# SEVERITY. EXACT and CASE-ONLY are errors — there is no benign version. PREFIX
# is a warning by default, because a well-populated frame usually resolves the
# exact name first and this cohort legitimately carries seven such pairs
# (APOE.geno, APOE_WGS, APOE4_count, …). Pass strict = TRUE to escalate, or
# allow = to whitelist the ones you have checked.
#
# USAGE
#   source("Utilities/covariate_lint.R")
#   assert_no_covariate_collision(d, biomarker_cols)                  # before a fit
#   assert_formula_safe(pTau_217 ~ age + sex + APOE4_count, biomarker_cols)
#   covariate_lint_report(d, biomarker_cols)                          # in a report
# =============================================================================

#' Find every way a covariate name could shadow a biomarker column.
#'
#' @param covariates  covariate / derived-column names. If NULL and `data` is
#'   given, defaults to every column that is not a biomarker — which is what
#'   catches derived columns you forgot you made.
#' @param biomarker_cols  the assay's marker column names.
#' @param data  optional data frame.
#' @return data frame with covariate, biomarker, kind (EXACT / CASE-ONLY /
#'   PREFIX); zero rows if clean.
find_covariate_collisions <- function(covariates = NULL, biomarker_cols, data = NULL) {
  if (is.null(covariates)) {
    if (is.null(data)) stop("find_covariate_collisions(): give covariates or data")
    covariates <- setdiff(names(data), biomarker_cols)
  }
  covariates <- unique(covariates); biomarker_cols <- unique(biomarker_cols)
  out <- list()

  exact <- intersect(covariates, biomarker_cols)
  if (length(exact))
    out$exact <- data.frame(covariate = exact, biomarker = exact, kind = "EXACT",
                            stringsAsFactors = FALSE)

  lc <- tolower(covariates); lb <- tolower(biomarker_cols)
  ci <- setdiff(covariates[lc %in% lb], exact)
  if (length(ci))
    out$case <- data.frame(covariate = ci,
                           biomarker = biomarker_cols[match(tolower(ci), lb)],
                           kind = "CASE-ONLY", stringsAsFactors = FALSE)

  # PREFIX: covariate begins with a biomarker name followed by a separator.
  # Requiring the separator avoids flagging genuinely distinct markers that merely
  # share an opening substring (e.g. a hypothetical APOE40 against APOE4).
  pre <- do.call(rbind, lapply(biomarker_cols, function(b) {
    hit <- covariates[startsWith(covariates, b) & covariates != b &
                        grepl(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", b),
                                     "[._]"), covariates)]
    if (!length(hit)) return(NULL)
    data.frame(covariate = hit, biomarker = b, kind = "PREFIX", stringsAsFactors = FALSE)
  }))
  if (!is.null(pre)) out$prefix <- pre[!pre$covariate %in% c(exact, ci), , drop = FALSE]

  res <- do.call(rbind, out)
  if (is.null(res)) res <- data.frame(covariate = character(), biomarker = character(),
                                      kind = character(), stringsAsFactors = FALSE)
  rownames(res) <- NULL
  res
}

#' Fail loudly if any covariate shadows a biomarker column.
#'
#' Deliberately an error, not a warning: a warning in a long knit scrolls past
#' and the run completes with wrong numbers.
#'
#' @param strict  treat PREFIX collisions as errors too (default: warn).
#' @param allow   covariate names to exempt.
assert_no_covariate_collision <- function(data = NULL, biomarker_cols,
                                          covariates = NULL, allow = character(0),
                                          strict = FALSE) {
  cl <- find_covariate_collisions(covariates, biomarker_cols, data)
  cl <- cl[!cl$covariate %in% allow, , drop = FALSE]
  if (!nrow(cl)) return(invisible(TRUE))

  fmt <- function(x) paste(sprintf("  %-26s <- %-22s (%s)", x$covariate, x$biomarker, x$kind),
                           collapse = "\n")
  hard <- cl[cl$kind %in% c("EXACT", "CASE-ONLY"), , drop = FALSE]
  soft <- cl[cl$kind == "PREFIX", , drop = FALSE]
  if (strict) { hard <- cl; soft <- soft[0, , drop = FALSE] }

  if (nrow(soft))
    warning("covariate-name PREFIX collision — R's partial matching means `d$",
            soft$biomarker[1], "` can silently return `", soft$covariate[1],
            "` when the exact column is absent. This is the APOE4_count bug.\n",
            fmt(soft), "\nPass allow= once checked, or strict=TRUE to make this an error.",
            call. = FALSE)

  if (nrow(hard))
    stop("covariate-name collision — a derived column shares a name with a ",
         "biomarker column.\nThe model would silently fit the wrong variable.\n",
         fmt(hard), "\n\nRename the covariate, or pass allow= if intended.",
         call. = FALSE)

  invisible(TRUE)
}

#' Same check, driven by a model formula. The response is expected to BE a
#' biomarker, so only the right-hand side is checked.
assert_formula_safe <- function(formula, biomarker_cols, allow = character(0),
                                strict = FALSE) {
  rhs <- if (length(formula) == 3L) all.vars(formula[[3]]) else all.vars(formula)
  assert_no_covariate_collision(biomarker_cols = biomarker_cols, covariates = rhs,
                                allow = allow, strict = strict)
}

#' List rather than stop — for a QC report section that should inform, not halt.
covariate_lint_report <- function(data, biomarker_cols) {
  cl <- find_covariate_collisions(NULL, biomarker_cols, data)
  if (!nrow(cl)) {
    cat("Covariate-name lint: PASS — no derived column shadows a biomarker.\n")
  } else {
    n_hard <- sum(cl$kind %in% c("EXACT", "CASE-ONLY"))
    cat(sprintf("Covariate-name lint: %d collision(s) — %d blocking, %d prefix\n",
                nrow(cl), n_hard, sum(cl$kind == "PREFIX")))
    for (i in seq_len(nrow(cl)))
      cat(sprintf("  %-26s <- %-22s (%s)\n", cl$covariate[i], cl$biomarker[i], cl$kind[i]))
    if (sum(cl$kind == "PREFIX"))
      cat("  PREFIX: `d$<biomarker>` can partial-match these when the exact column\n",
          "         is absent from the frame being subscripted. Prefer d[[\"name\"]].\n")
  }
  invisible(cl)
}
