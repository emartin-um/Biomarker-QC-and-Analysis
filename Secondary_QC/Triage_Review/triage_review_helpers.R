# =============================================================================
# triage_review_helpers.R  ·  added 2026-07-27
# -----------------------------------------------------------------------------
# What the Primary QC screens actually did, and whether their thresholds sit
# anywhere useful. Base R only, I/O-free, so the runner and the test share one
# definition of every rule.
#
# WHY THIS EXISTS. Primary QC removes wells (triage) and marks wells (flags), and
# until now nothing checked the screens themselves: whether an axis fires at all,
# whether two axes are the same test twice, where a cut sits relative to its own
# statistic, or whether a screen is keying on something other than sample quality.
# On the 50-plate run all four questions had a surprising answer.
#
# THE FOUR THINGS IT CHECKS
#
# 1. IS THE AXIS DOING WORK? A gate that never fires is not a safe gate, it is an
#    unmeasured one. Verified 2026-07-27: PCA round 2 flagged 0 of 4184 wells
#    (max PC_z 4.824 against a cut of 5.0); Alamar's IC_Reads < 1000 rule flagged
#    0 of 4200; its Reads < 500000 rule flagged 1, already caught by
#    detectability. The lower read bound cannot fire BY CONSTRUCTION -- see
#    read_gate_geometry().
#
# 2. IS THE AXIS INDEPENDENT? The in-house "IC outlier" test is the vendor's own
#    "IC Median within +/-40% of plate median" recomputed on the same mCherry
#    quantity: 57 of 58 wells agree. It supplies 58 of the 77 triage removals, so
#    three quarters of triage is a re-derivation of a metric the vendor already
#    reports. That is not wrong, but it must not be presented as independent
#    in-house QC.
#
# 3. WHERE DOES THE CUT SIT? The burden gate's FDR is a deterministic function of
#    an integer count over 18 anchors, so "FDR < 0.01" IS "5 or more anchors
#    out". It is not a dial. 91.5% of wells sit at fdr = 1.0 exactly and the
#    statistic takes only 14 distinct values in 4200 wells.
#
# 4. IS THE AXIS KEYING ON POSITION RATHER THAN SPECIMEN? Well F7 is triaged in
#    9 of its 50 plate-bays (18.0%) against 1.64% for every other position
#    (permutation p < 1e-4), all 9 through the IC axis. And IC recovery declines
#    from plate row A to row H on 50 of 50 plate-bays (median Spearman -0.738).
#    A screen that partly encodes well position is removing specimens for a
#    reason that has nothing to do with the specimen.
#
# 🚫 NO DROP RULE HERE, AND NO PROPOSAL TO ADD ONE. Everything is diagnostic. The
# study's standing position is that per-specimen nuisance axes belong in the
# analysis model, not on a drop list.
# =============================================================================

## --- small helpers -----------------------------------------------------------
.num <- function(v) suppressWarnings(as.numeric(gsub("[%,]", "", as.character(v))))

#' Standardised mean difference (Cohen's d) with a pooled SD.
cohen_d <- function(x, g) {
  a <- x[g & !is.na(x)]; b <- x[!g & !is.na(x)]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  s <- sqrt(((length(a) - 1) * stats::var(a) + (length(b) - 1) * stats::var(b)) /
              (length(a) + length(b) - 2))
  if (!is.finite(s) || s == 0) return(NA_real_)
  (mean(a) - mean(b)) / s
}

#' AUC of x separating g from !g, via the rank identity (no pROC dependency).
auc <- function(x, g) {
  ok <- !is.na(x); x <- x[ok]; g <- g[ok]
  na <- sum(g); nb <- sum(!g)
  if (na < 1 || nb < 1) return(NA_real_)
  (sum(rank(x)[g]) - na * (na + 1) / 2) / (na * nb)
}

#' Permutation p for any statistic of a logical label, holding the label count
#' fixed and reshuffling which wells carry it.
#'
#' Used wherever a cell count is too small for a chi-square to be valid -- 77
#' triaged wells over 84 plate positions or 23 sites means most cells are 0 or 1,
#' and an asymptotic test there is not merely conservative, it is wrong.
perm_test <- function(label, stat_fun, n_perm = 10000L, seed = 1L) {
  set.seed(seed)
  obs <- stat_fun(label)
  n <- length(label); k <- sum(label)
  null <- replicate(n_perm, {
    l <- logical(n); l[sample.int(n, k)] <- TRUE
    stat_fun(l)
  })
  list(observed = obs, null_mean = mean(null),
       p = (1 + sum(null >= obs)) / (n_perm + 1))
}

## --- 1. the triage axes, without the precedence collapse ---------------------

#' Split the collapsed triage Reason string into its independent axes.
#' Must stay in step with Primary QC's `triage_df` construction.
triage_axes <- function(reason) {
  r <- as.character(reason); r[is.na(r)] <- ""
  data.frame(
    IC      = grepl("IC outlier",                r, fixed = TRUE),
    PCA     = grepl("PCA outlier",               r, fixed = TRUE),
    burden  = grepl("High outlier burden (FDR)", r, fixed = TRUE),
    baddata = grepl("Bad Data",                  r, fixed = TRUE),
    stringsAsFactors = FALSE)
}

#' The true axis tabulation, and what a single-reason report loses.
#'
#' Cohort_Accounting (and any other one-reason-per-well report) resolves a
#' multi-axis well by a fixed precedence, so its per-axis counts are NOT axis
#' totals. On the 50-plate run the precedence is IC > PCA > burden and it
#' undercounts the PCA axis by 37.5% (10 reported vs 16 real) and the burden axis
#' by 64.0% (9 vs 25), while leaving IC exact. Both numbers are correct answers to
#' different questions; this function reports both so neither is mistaken for the
#' other.
#'
#' @param ax data.frame of logical axis columns, one row per TRIAGED well.
axis_overlap <- function(ax) {
  nm <- names(ax)
  n_axes <- rowSums(ax)
  combo <- apply(ax, 1, function(r) {
    s <- nm[as.logical(r)]
    if (!length(s)) "none" else paste(s, collapse = "+")
  })
  tot <- vapply(ax, sum, integer(1))
  sole <- vapply(nm, function(v) sum(ax[[v]] & n_axes == 1), integer(1))
  data.frame(
    axis = nm,
    n_wells_on_axis = unname(tot),
    n_wells_sole_reason = unname(sole),
    pct_of_axis_shared = round(100 * (1 - sole / pmax(tot, 1)), 1),
    stringsAsFactors = FALSE, row.names = NULL)
}

#' The full combination table (the Venn, flattened).
axis_combinations <- function(ax) {
  nm <- names(ax)
  combo <- apply(ax, 1, function(r) {
    s <- nm[as.logical(r)]; if (!length(s)) "none" else paste(s, collapse = "+") })
  t <- table(combo)
  data.frame(combination = names(t), n_wells = as.integer(t),
             n_axes = vapply(names(t), function(s)
               if (s == "none") 0L else length(strsplit(s, "\\+")[[1]]), integer(1)),
             stringsAsFactors = FALSE, row.names = NULL)
}

## --- 2. is an in-house axis a duplicate of a vendor metric? ------------------

#' Jaccard agreement between two logical flags, with the discordant counts.
#'
#' Used to ask whether an in-house screen is independent information. An axis at
#' Jaccard ~1 against a vendor metric is the vendor metric, and reporting it as a
#' separate QC axis double-counts the evidence.
flag_agreement <- function(a, b, name_a = "a", name_b = "b") {
  a <- !is.na(a) & a; b <- !is.na(b) & b
  both <- sum(a & b); only_a <- sum(a & !b); only_b <- sum(b & !a)
  data.frame(
    flag_a = name_a, flag_b = name_b,
    n_a = sum(a), n_b = sum(b), n_both = both,
    n_only_a = only_a, n_only_b = only_b,
    jaccard = round(both / max(both + only_a + only_b, 1), 4),
    verdict = if (both + only_a + only_b == 0) "neither fires"
              else if (both / (both + only_a + only_b) >= 0.9)
                "DUPLICATE - not independent information"
              else if (both / (both + only_a + only_b) >= 0.5) "substantial overlap"
              else "largely independent",
    stringsAsFactors = FALSE, row.names = NULL)
}

## --- 3. does a gate fire, and where does its cut sit? -----------------------

#' Inventory one gate: how many wells it fires on, and whether it is inert.
#'
#' `n_fired == 0` is the headline. A gate that never fires has not been shown to
#' be safe; it has been shown to be untested. Report it every run, because the
#' run where it starts firing is the run you want to know about.
gate_row <- function(gate, rule, fired, n_total, note = "") {
  n <- sum(!is.na(fired) & fired)
  data.frame(gate = gate, rule = rule, n_fired = n, n_wells = n_total,
             pct = round(100 * n / max(n_total, 1), 3),
             status = if (n == 0) "INERT - never fires on this run"
                      else if (n / n_total > 0.25) "fires on >25% of wells"
                      else "active",
             note = note, stringsAsFactors = FALSE, row.names = NULL)
}

#' Operating points of a quantised gate.
#'
#' The burden gate's FDR is a deterministic function of an integer anchor count,
#' so the "threshold" is a choice among a handful of integer cuts, not a
#' continuous dial. This enumerates them: for each distinct value of the
#' statistic, how many wells sit there and how many would be removed at that cut.
#' That table is what a threshold discussion should actually be about.
#'
#' @param stat  the gate statistic (lower = more extreme, e.g. an FDR).
#' @param count the integer count behind it, if any (else NULL).
#' @param cut   the cut currently in force.
gate_operating_points <- function(stat, count = NULL, cut = NULL, lower_is_worse = TRUE) {
  ok <- !is.na(stat)
  u <- sort(unique(stat[ok]), decreasing = !lower_is_worse)
  out <- do.call(rbind, lapply(u, function(v) {
    at <- stat == v & ok
    rm <- if (lower_is_worse) stat <= v & ok else stat >= v & ok
    data.frame(stat_value = v,
               count_value = if (is.null(count)) NA_integer_
                             else as.integer(stats::median(count[at], na.rm = TRUE)),
               n_wells_at = sum(at),
               n_removed_at_this_cut = sum(rm),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  out$is_current_cut <- if (is.null(cut)) NA
                        else if (lower_is_worse) out$stat_value < cut else out$stat_value > cut
  out
}

#' Is the gate's cut near anything, or is it isolated in a gap?
#'
#' A cut sitting in an empty gap is insensitive to its exact value (good), but it
#' also means the "threshold" is really the gap edge. Reports the nearest well on
#' each side.
gate_cut_neighbourhood <- function(stat, fired, label = "gate") {
  ok <- !is.na(stat)
  inside <- stat[ok & fired]; outside <- stat[ok & !fired]
  data.frame(gate = label,
             n_fired = length(inside),
             worst_retained = if (length(outside)) min(outside) else NA_real_,
             best_removed  = if (length(inside)) max(inside) else NA_real_,
             gap = if (length(inside) && length(outside)) min(outside) - max(inside) else NA_real_,
             stringsAsFactors = FALSE, row.names = NULL)
}

#' Does the IQR arm of the lower read bound ever bind?
#'
#' Primary QC sets Lower_Bound = max(500000, median - k*IQR) per Run x Bay. When
#' k*IQR exceeds the median the subtraction goes negative and the 500000 floor
#' always wins. Verified 2026-07-27: negative in 50 of 50 plate-bays (min
#' -31.7M against a bay median of 5.67M reads). So the lower bound is never
#' PLATE-RELATIVE -- a well can only be flagged low by the ABSOLUTE 500000 floor,
#' which caught exactly 1 well on this run. The upper bound has no floor and is
#' fully plate-relative, which is why the flag is 93 HIGH to 1 LOW.
#'
#' The consequence worth stating: "read outlier" sounds symmetric and is not. A
#' chronically low-input site cannot be flagged by this gate however far it sits
#' below its plate, unless it also drops under half a million reads.
read_gate_geometry <- function(reads, run_bay, k = 4) {
  sp <- split(.num(reads), run_bay)
  do.call(rbind, lapply(names(sp), function(b) {
    v <- sp[[b]][!is.na(sp[[b]])]
    if (length(v) < 4) return(NULL)
    q <- stats::quantile(v, c(.25, .75)); iqr <- q[2] - q[1]
    med <- stats::median(v)
    raw_lower <- med - k * iqr
    data.frame(run_bay = b, n = length(v), median_reads = med, iqr = unname(iqr),
               raw_lower_bound = unname(raw_lower),
               floor_binds = raw_lower < 5e5,
               min_reads = min(v),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
}

## --- 4. is the screen keying on plate position? ------------------------------

#' Per-well-position rates, with a permutation p on the maximum.
#'
#' 84 patient positions x 50 plate-bays. A specimen-quality screen should be flat
#' across positions; a plate-handling artifact is not. The test is on the MAXIMUM
#' per-position count, which is the right statistic for "is any one position
#' anomalous" and needs no multiple-testing correction of its own.
position_effects <- function(flag, well_row, well_col, n_perm = 10000L, seed = 1L) {
  pos <- paste0(well_row, well_col)
  lv <- unique(pos)
  tab <- table(factor(pos[flag], levels = lv))
  n_at <- table(factor(pos, levels = lv))
  set.seed(seed)
  k <- sum(flag)
  null_max <- replicate(n_perm, {
    l <- logical(length(pos)); l[sample.int(length(pos), k)] <- TRUE
    max(table(factor(pos[l], levels = lv)))
  })
  out <- data.frame(
    position = lv,
    n_wells = as.integer(n_at[lv]),
    n_flagged = as.integer(tab[lv]),
    pct = round(100 * as.integer(tab[lv]) / pmax(as.integer(n_at[lv]), 1), 2),
    stringsAsFactors = FALSE, row.names = NULL)
  out <- out[order(-out$n_flagged), ]
  attr(out, "perm") <- list(observed_max = max(tab), null_mean_max = mean(null_max),
                            p = (1 + sum(null_max >= max(tab))) / (n_perm + 1))
  out
}

#' Row-position gradient in a per-well quantity, per plate-bay.
#'
#' The IC-median gate keys on a quantity that declines down the plate: Spearman
#' of (row order, mean IC_Median) is negative on 50 of 50 plate-bays, median
#' -0.738. NB the decline is a consistent TENDENCY, not a monotone one -- 0 of 50
#' bays are strictly monotone -- so report the distribution of rho, not "it
#' decreases".
row_gradient <- function(x, well_row, run_bay) {
  sp <- split(seq_along(x), run_bay)
  rho <- vapply(sp, function(ix) {
    m <- tapply(x[ix], well_row[ix], mean, na.rm = TRUE)
    m <- m[order(names(m))]
    if (length(m) < 3 || all(is.na(m))) return(NA_real_)
    suppressWarnings(stats::cor(seq_along(m), m, method = "spearman", use = "complete.obs"))
  }, numeric(1))
  data.frame(run_bay = names(rho), rho_row_vs_value = round(unname(rho), 3),
             stringsAsFactors = FALSE, row.names = NULL)
}

## --- 5. what do the removed and marked wells look like on the known axes? ----

#' Effect size and AUC of each quality axis for a given well class.
#'
#' Reported PER TRIAGE CLASS, not for "triaged" as a whole. That matters: on the
#' 50-plate run triaged wells average mean_INT_z -1.04, but the AUC is only 0.408
#' and the deficit is concentrated in the 16 PCA-triaged wells (mean -4.56) --
#' the 50 IC-only wells sit at -0.21 and are near-indistinguishable from retained
#' wells on every axis except the one that removed them. Pooling the classes
#' manufactures a "triaged wells are dim" story that does not describe 65% of them.
axis_attribution <- function(d, classes, axes) {
  axes <- intersect(axes, names(d))
  do.call(rbind, lapply(names(classes), function(cl) {
    g <- classes[[cl]]
    if (sum(g) < 3) return(NULL)
    do.call(rbind, lapply(axes, function(a) {
      x <- .num(d[[a]])
      data.frame(class = cl, n = sum(g), axis = a,
                 mean_in = round(mean(x[g], na.rm = TRUE), 3),
                 mean_out = round(mean(x[!g], na.rm = TRUE), 3),
                 cohen_d = round(cohen_d(x, g), 3),
                 auc = round(auc(x, g), 3),
                 stringsAsFactors = FALSE, row.names = NULL)
    }))
  }))
}

#' Rate of a flag by group, with a permutation p and a depth-conditioned check.
#'
#' Site and depth are confounded (depth is ~37% a site property), so a raw
#' site-rate difference is not evidence of a site QC problem. The conditioned
#' column refits the rate against the group after regressing out depth; if the
#' association is composition it collapses.
rate_by_group <- function(flag, group, depth = NULL, min_n = 5, n_perm = 10000L, seed = 1L) {
  ok <- !is.na(group)
  g <- as.character(group)[ok]; f <- flag[ok]
  lv <- names(which(table(g) >= min_n))
  if (!length(lv)) return(NULL)
  set.seed(seed)
  k <- sum(f)
  out <- do.call(rbind, lapply(lv, function(v) {
    i <- g == v
    data.frame(group = v, n = sum(i), n_flagged = sum(f & i),
               pct = round(100 * sum(f & i) / sum(i), 2),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  # chi-square-like dispersion across groups, permuted (valid at tiny cell counts)
  stat <- function(lab) {
    e <- vapply(lv, function(v) sum(lab[g == v]), numeric(1))
    n <- vapply(lv, function(v) sum(g == v), numeric(1))
    max(e / n)
  }
  null <- replicate(n_perm, {
    l <- logical(length(f)); l[sample.int(length(f), k)] <- TRUE; stat(l) })
  attr(out, "perm") <- list(observed_max_rate = stat(f),
                            null_mean_max_rate = mean(null),
                            p = (1 + sum(null >= stat(f))) / (n_perm + 1))
  if (!is.null(depth)) {
    dd <- .num(depth)[ok]
    m0 <- stats::glm(f ~ 1, family = stats::binomial())
    m1 <- suppressWarnings(stats::glm(f ~ factor(g), family = stats::binomial()))
    m2 <- suppressWarnings(stats::glm(f ~ factor(g) + dd, family = stats::binomial()))
    attr(out, "conditioned") <- list(
      dev_group_only = round(stats::deviance(m0) - stats::deviance(m1), 2),
      dev_group_given_depth = round(stats::deviance(stats::glm(f ~ dd, family = stats::binomial())) -
                                      stats::deviance(m2), 2))
  }
  out[order(-out$pct), ]
}
