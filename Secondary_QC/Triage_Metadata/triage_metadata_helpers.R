# =============================================================================
# triage_metadata_helpers.R  ·  Secondary QC — triage/flag decisions vs metadata
# -----------------------------------------------------------------------------
# BASE R ONLY. Every inferential choice here exists because the obvious version
# of it is wrong on this data:
#
#   * chi-square is wrong  — 77 triaged wells over 23 sites leaves most cells at
#     0 or 1, where the asymptotic test is not conservative but anti-conservative.
#     Everything is permuted.
#
#   * the maximum rate is not enough — on the 50-plate run the largest site rate
#     is 40.0%, and it is 2 wells out of 5. A statistic that keys on the extreme
#     answers "is any one group unusual", which a tiny group can satisfy on
#     noise. An omnibus that pools evidence across groups is reported beside it,
#     and the two disagreeing is itself the finding.
#
#   * a raw group difference is not evidence of a group problem — site is
#     confounded with plate and with read depth. Deviance is refit conditioned on
#     each, and a signal that collapses was composition.
#
#   * "not significant" is not "no effect" — every null carries the smallest
#     effect it could have detected.
# =============================================================================

.num <- function(x) if (is.numeric(x)) x else
  suppressWarnings(as.numeric(gsub("[%,]", "", as.character(x))))

.blank_to_na <- function(x) { x <- as.character(x)
  x[is.na(x) | trimws(x) == "" | tolower(trimws(x)) %in% c("na", "n/a")] <- NA_character_; x }

# -----------------------------------------------------------------------------
# 1. Proportions with honest uncertainty
# -----------------------------------------------------------------------------

#' Wilson score interval. Used instead of a Wald interval because the rates here
#' are small (1-2%) and the group sizes range from 5 to 654 — Wald gives
#' impossible negative lower bounds and absurdly narrow intervals at n = 5,
#' which is exactly where the eye-catching rates live.
wilson_ci <- function(k, n, conf = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n; d <- 1 + z^2 / n
  c(max(0, (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d),
    min(1, (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d))
}

#' Smallest detectable rate RATIO at the given power, for a group of size n
#' against a baseline rate p0 over n0 wells.
#'
#' The point: a report that says "site is not associated with triage" is worth
#' nothing without this number. At a 1.8% base rate and n = 50, nothing short of
#' a ~4x elevation is detectable, so a null is close to uninformative — and it
#' must say so rather than reading as reassurance.
min_detectable_rr <- function(n, p0, n0, power = 0.80, alpha = 0.05) {
  if (!is.finite(n) || !is.finite(p0) || n < 1 || p0 <= 0 || p0 >= 1) return(NA_real_)
  za <- stats::qnorm(1 - alpha / 2); zb <- stats::qnorm(power)
  # Solve on the arcsine-root scale, which is variance-stabilising for a
  # proportion and avoids the normal approximation collapsing at small p.
  h <- (za + zb) * sqrt(1 / n + 1 / max(n0, 1))
  p1 <- sin(asin(sqrt(p0)) + h / 2)^2
  round(p1 / p0, 2)
}

# -----------------------------------------------------------------------------
# 2. Association of a binary flag with a grouping label
# -----------------------------------------------------------------------------

#' Per-group rates, with a size floor that REPORTS what it excluded.
#'
#' `min_n` groups are not silently dropped — they come back in the `below_floor`
#' attribute, because "the signal was in the sites we did not show you" is a
#' failure mode this module is specifically meant to prevent.
rates_by_group <- function(flag, group, min_n = 20L, conf = 0.95) {
  ok <- !is.na(group) & !is.na(flag)
  g <- as.character(group)[ok]; f <- as.logical(flag)[ok]
  if (!length(g)) return(NULL)
  tab <- table(g)
  mk <- function(lv) do.call(rbind, lapply(lv, function(v) {
    i <- g == v; k <- sum(f & i); n <- sum(i); ci <- wilson_ci(k, n, conf)
    data.frame(group = v, n = n, n_flagged = k, pct = round(100 * k / n, 2),
               ci_lo_pct = round(100 * ci[1], 2), ci_hi_pct = round(100 * ci[2], 2),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  keep <- names(tab)[tab >= min_n]; drop <- names(tab)[tab < min_n]
  out <- mk(keep)
  if (!is.null(out)) out <- out[order(-out$pct), ]
  attr(out, "below_floor") <- if (length(drop)) mk(drop) else NULL
  attr(out, "min_n") <- min_n
  out
}

#' Permutation test of flag-vs-group, with TWO statistics.
#'
#' @param stratum optional; if given, the permutation shuffles the flag WITHIN
#'   each stratum, so a stratum-level effect (e.g. a bad plate) cannot manufacture
#'   a group signal. This is the difference between "site matters" and "site
#'   matters for a reason other than which plates its samples landed on".
#'
#' Returns both `max_rate` (sensitive to one bad group, hostage to small ones)
#' and `deviance` (omnibus, pools evidence). They are reported side by side on
#' purpose: max_rate significant with deviance null means one small group is
#' carrying the whole result.
perm_group_test <- function(flag, group, stratum = NULL, min_n = 20L,
                            n_perm = 10000L, seed = 1L) {
  ok <- !is.na(group) & !is.na(flag)
  g <- as.character(group)[ok]; f <- as.logical(flag)[ok]
  st <- if (is.null(stratum)) rep("all", length(f)) else as.character(stratum)[ok]
  tab <- table(g); lv <- names(tab)[tab >= min_n]
  if (length(lv) < 2 || sum(f) < 2)
    return(data.frame(n_groups = length(lv), n_flagged = sum(f),
                      observed_max_rate = NA_real_, p_max_rate = NA_real_,
                      observed_deviance = NA_real_, p_deviance = NA_real_,
                      stringsAsFactors = FALSE))
  keep <- g %in% lv
  g <- g[keep]; f <- f[keep]; st <- st[keep]
  gi <- match(g, lv); ng <- tabulate(gi, nbins = length(lv))

  stat_max <- function(fl) max(tabulate(gi[fl], nbins = length(lv)) / ng)
  # G^2 against the null of a common rate. Pools evidence instead of keying on
  # the extreme; a small group can move it only as far as its own weight allows.
  stat_dev <- function(fl) {
    k <- tabulate(gi[fl], nbins = length(lv)); p <- sum(k) / sum(ng)
    if (p <= 0 || p >= 1) return(0)
    e1 <- ng * p; e0 <- ng * (1 - p); o1 <- k; o0 <- ng - k
    2 * (sum(ifelse(o1 > 0, o1 * log(o1 / e1), 0)) + sum(ifelse(o0 > 0, o0 * log(o0 / e0), 0)))
  }
  obs_max <- stat_max(f); obs_dev <- stat_dev(f)

  set.seed(seed)
  idx_by_st <- split(seq_along(f), st)
  k_by_st   <- vapply(split(f, st), sum, integer(1))
  null <- vapply(seq_len(n_perm), function(i) {
    fl <- logical(length(f))
    for (j in seq_along(idx_by_st)) {
      ii <- idx_by_st[[j]]; kk <- k_by_st[[j]]
      if (kk > 0) fl[ii[sample.int(length(ii), kk)]] <- TRUE
    }
    c(stat_max(fl), stat_dev(fl))
  }, numeric(2))

  data.frame(
    n_groups = length(lv), n_flagged = sum(f),
    observed_max_rate = round(100 * obs_max, 2),
    null_mean_max_rate = round(100 * mean(null[1, ]), 2),
    p_max_rate = (1 + sum(null[1, ] >= obs_max)) / (n_perm + 1),
    observed_deviance = round(obs_dev, 2),
    null_mean_deviance = round(mean(null[2, ]), 2),
    p_deviance = (1 + sum(null[2, ] >= obs_dev)) / (n_perm + 1),
    stratified_by = if (is.null(stratum)) "" else "yes",
    stringsAsFactors = FALSE)
}

#' Group deviance before and after conditioning on a confounder.
#'
#' Two confounders matter here and they are different in kind:
#'   plate  — a factor. Sites ship in batches, so "this site fails more" and
#'            "this site's samples landed on worse plates" are the same data.
#'   depth  — continuous. Read depth is partly a site property, so a site rate
#'            difference can be a depth composition difference.
#' If `dev_group_given_X` collapses toward 0, the association was composition.
conditioned_deviance <- function(flag, group, plate = NULL, depth = NULL, min_n = 20L) {
  ok <- !is.na(group) & !is.na(flag)
  g <- as.character(group)[ok]; f <- as.logical(flag)[ok]
  tab <- table(g); lv <- names(tab)[tab >= min_n]
  keep <- g %in% lv
  if (sum(keep) < 10 || length(lv) < 2 || sum(f[keep]) < 2) return(NULL)
  g <- factor(g[keep]); f <- f[keep]
  pl <- if (is.null(plate)) NULL else factor(as.character(plate)[ok][keep])
  dp <- if (is.null(depth)) NULL else .num(depth)[ok][keep]

  dev <- function(fml, data) suppressWarnings(
    stats::deviance(stats::glm(fml, data = data, family = stats::binomial())))
  d <- data.frame(f = f, g = g)
  out <- data.frame(dev_group_only = round(dev(f ~ 1, d) - dev(f ~ g, d), 2),
                    df_group = length(lv) - 1L, stringsAsFactors = FALSE)
  if (!is.null(pl) && nlevels(pl) > 1) {
    d$pl <- pl
    out$dev_group_given_plate <- round(dev(f ~ pl, d) - dev(f ~ pl + g, d), 2)
  }
  if (!is.null(dp) && sum(is.finite(dp)) > 10) {
    d$dp <- dp
    out$dev_group_given_depth <- round(dev(f ~ dp, d) - dev(f ~ dp + g, d), 2)
  }
  out
}

#' Mutual adjustment: which correlated label is actually carrying the signal?
#'
#' Site, ancestry and diagnosis are not independent — sites recruit particular
#' populations, and diagnosis composition follows. Testing each one alone will
#' therefore light up several of them off a single underlying cause, and the
#' one-at-a-time table cannot say which. This fits the full model and reports
#' each variable's deviance GIVEN the others.
#'
#' Read it as: a variable whose deviance survives adjustment is contributing
#' something the others do not. A variable that collapses was a proxy.
#'
#' Deliberately NOT a model-selection procedure — no stepwise, no p-value
#' shopping. The full model is fitted once and every term reported, so a term
#' that drops out is visible rather than absent.
mutual_adjustment <- function(flag, labels, min_n = 20L, depth = NULL) {
  vars <- names(labels)
  d <- data.frame(f = as.logical(flag), stringsAsFactors = FALSE)
  used <- character(0)
  for (v in vars) {
    x <- as.character(labels[[v]])
    tab <- table(x)
    x[!x %in% names(tab)[tab >= min_n]] <- NA   # rare levels out, as elsewhere
    d[[v]] <- factor(x)
    if (nlevels(d[[v]]) >= 2) used <- c(used, v)
  }
  if (!is.null(depth)) { d$.depth <- .num(depth); used <- c(used, ".depth") }
  d <- d[stats::complete.cases(d[, c("f", used), drop = FALSE]), , drop = FALSE]
  if (nrow(d) < 30 || sum(d$f) < 3 || !length(used)) return(NULL)

  dv <- function(terms) suppressWarnings(stats::deviance(stats::glm(
    stats::as.formula(paste("f ~", if (length(terms)) paste(terms, collapse = " + ") else "1")),
    data = d, family = stats::binomial())))
  full <- dv(used)
  do.call(rbind, lapply(used, function(v) {
    df_v <- if (v == ".depth") 1L else nlevels(d[[v]]) - 1L
    dev_alone <- dv(character(0)) - dv(v)
    dev_adj   <- dv(setdiff(used, v)) - full
    data.frame(variable = ifelse(v == ".depth", "log2(reads)", v),
               df = df_v, n_wells = nrow(d), n_flagged = sum(d$f),
               dev_alone = round(dev_alone, 2),
               dev_given_others = round(dev_adj, 2),
               p_given_others = signif(stats::pchisq(dev_adj, df_v, lower.tail = FALSE), 3),
               retained_share = if (dev_alone > 0) round(dev_adj / dev_alone, 2) else NA_real_,
               stringsAsFactors = FALSE, row.names = NULL)
  }))
}

# -----------------------------------------------------------------------------
# 3. Metadata tidying — explicit, and written out
# -----------------------------------------------------------------------------

#' Collapse the free-text diagnosis columns to analysis levels.
#'
#' `Case_Control` arrives with 15 distinct values on the 50-plate run, several of
#' which are administrative rather than clinical ("Missing", "Insufficient
#' Data"). Collapsing is unavoidable; doing it silently is not. The mapping is
#' returned so it can be written to disk and argued with.
collapse_diagnosis <- function(x) {
  v <- .blank_to_na(x); lo <- tolower(trimws(v))
  out <- rep(NA_character_, length(v))
  out[lo %in% c("control", "nci", "0")]                        <- "Control"
  out[lo %in% c("case", "ad")]                                 <- "Case (AD)"
  out[lo %in% c("mci", "cind", "cognitively impaired")]        <- "MCI/CIND"
  out[grepl("dementia|lewy|vascular", lo)]                     <- "Other dementia"
  out[lo %in% c("other", "other medical diagnosis")]           <- "Other"
  # Administrative non-answers are their OWN level, never folded into a clinical
  # one and never dropped: on this run they carry the highest triage rate seen,
  # which is a finding about data capture, not about dementia.
  out[lo %in% c("missing", "insufficient data")]               <- "Not codeable"
  out[is.na(v)]                                                <- "Not codeable"
  out[is.na(out)]                                              <- "Other"
  out
}

collapse_map <- function(x) {
  v <- .blank_to_na(x)
  m <- data.frame(raw = ifelse(is.na(v), "<NA>", v), collapsed = collapse_diagnosis(v),
                  stringsAsFactors = FALSE)
  a <- aggregate(list(n_wells = rep(1L, nrow(m))), by = list(raw = m$raw, collapsed = m$collapsed), FUN = sum)
  a[order(a$collapsed, -a$n_wells), ]
}

#' Ancestry, using the pipeline's OWN curated rule table.
#'
#' `Metadata_Merge/review/Ancestry_review.csv` is a Group x Race x Ethnicity
#' lookup that the merge already applies. Re-deriving ancestry from the raw
#' `Race`/`Ethnicity` codes here would create a second definition in the same
#' repo, and the first attempt at exactly that silently returned "Unknown" for
#' all 4087 wells — the codes are two-letter (BL/WH/MU/AI/HP, HI/NH), not the
#' spelled-out words a regex would expect. So: read the table, and if it is
#' missing say so in the label rather than guessing.
#'
#' Mirrors build_ancestry_rules()/assign_ancestry() in
#' Metadata_Merge/Cohort_Accounting/cohort_accounting_helpers.R, including the
#' two easily-misread conventions:
#'   * a blank condition column is a WILDCARD, making that rule broader
#'   * rules are applied in REVERSE order, so the FIRST matching row wins
#' If that file's conventions change, this must change with it.
assign_ancestry_from_review <- function(d, anc_review, other_label = "Other") {
  if (is.null(anc_review) || !"FINAL_Ancestry" %in% names(anc_review))
    return(rep(NA_character_, nrow(d)))
  f <- trimws(as.character(anc_review$FINAL_Ancestry)); f[f == ""] <- NA
  keep <- which(!is.na(f))
  out <- rep(other_label, nrow(d))
  for (i in rev(keep)) {
    cond <- list(Group = anc_review$Group[i], Race = anc_review$Race[i],
                 Ethnicity = anc_review$Ethnicity[i])
    cond <- cond[!vapply(cond, function(v) is.na(v) || trimws(as.character(v)) == "",
                         logical(1))]
    mask <- rep(TRUE, nrow(d))
    for (cn in names(cond)) {
      if (!cn %in% names(d)) { mask <- rep(FALSE, nrow(d)); break }
      cm <- as.character(d[[cn]]) == as.character(cond[[cn]])
      cm[is.na(cm)] <- FALSE
      mask <- mask & cm
    }
    out[mask] <- f[i]
  }
  out
}

age_band <- function(a) {
  x <- .num(a)
  ifelse(is.na(x), NA_character_,
  ifelse(x < 60, "<60", ifelse(x < 70, "60-69", ifelse(x < 80, "70-79", "80+"))))
}

# -----------------------------------------------------------------------------
# 4. Balance table: who got removed vs who stayed
# -----------------------------------------------------------------------------

#' Standardised difference for a continuous or binary variable. Reported instead
#' of a p-value because at n = 4200 vs 77 a p-value mostly measures the sample
#' size; the standardised difference measures the imbalance itself.
std_diff <- function(x, grp) {
  x <- .num(x); a <- x[grp & !is.na(x)]; b <- x[!grp & !is.na(x)]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  s <- sqrt((stats::var(a) + stats::var(b)) / 2)
  if (!is.finite(s) || s == 0) return(NA_real_)
  round((mean(a) - mean(b)) / s, 3)
}

balance_table <- function(d, grp, continuous, categorical) {
  rows <- list()
  for (v in continuous) if (v %in% names(d)) {
    x <- .num(d[[v]])
    rows[[length(rows) + 1]] <- data.frame(
      variable = v, level = "(mean)",
      n_in = sum(grp & !is.na(x)), n_out = sum(!grp & !is.na(x)),
      val_in = round(mean(x[grp], na.rm = TRUE), 2),
      val_out = round(mean(x[!grp], na.rm = TRUE), 2),
      std_diff = std_diff(x, grp), stringsAsFactors = FALSE)
  }
  for (v in categorical) if (v %in% names(d)) {
    lv <- sort(unique(stats::na.omit(as.character(d[[v]]))))
    for (l in lv) {
      ind <- !is.na(d[[v]]) & d[[v]] == l
      rows[[length(rows) + 1]] <- data.frame(
        variable = v, level = l,
        n_in = sum(grp & ind), n_out = sum(!grp & ind),
        val_in = round(100 * sum(grp & ind) / max(1, sum(grp & !is.na(d[[v]]))), 2),
        val_out = round(100 * sum(!grp & ind) / max(1, sum(!grp & !is.na(d[[v]]))), 2),
        std_diff = std_diff(as.numeric(ind), grp), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------
# 5. How separable are two design factors?
# -----------------------------------------------------------------------------

#' Cramer's V plus the spread that actually decides whether a conditioned test is
#' estimable. V alone can look fine while every group sits on one plate; the
#' `*_per_*` and concentration columns are what say whether the question can be
#' asked at all.
confounding_summary <- function(a, b, name_a = "A", name_b = "B") {
  ok <- !is.na(a) & !is.na(b)
  tb <- table(as.character(a)[ok], as.character(b)[ok])
  if (any(dim(tb) < 2)) return(NULL)
  chi <- suppressWarnings(stats::chisq.test(tb)$statistic)
  data.frame(
    factor_a = name_a, factor_b = name_b,
    n_a = nrow(tb), n_b = ncol(tb),
    cramers_v = round(sqrt(as.numeric(chi) / (sum(tb) * (min(dim(tb)) - 1))), 3),
    b_per_a_min = min(apply(tb, 1, function(r) sum(r > 0))),
    b_per_a_median = stats::median(apply(tb, 1, function(r) sum(r > 0))),
    b_per_a_max = max(apply(tb, 1, function(r) sum(r > 0))),
    max_share_of_a_on_one_b = round(max(apply(tb, 1, function(r) max(r) / sum(r))), 3),
    median_share_of_a_on_one_b = round(stats::median(apply(tb, 1, function(r) max(r) / sum(r))), 3),
    verdict = if (stats::median(apply(tb, 1, function(r) max(r) / sum(r))) > 0.8)
      paste0(name_a, " is effectively nested in ", name_b, " - a conditioned test is not estimable")
      else paste0(name_a, " spreads across ", name_b, " - a conditioned test is estimable"),
    stringsAsFactors = FALSE)
}
