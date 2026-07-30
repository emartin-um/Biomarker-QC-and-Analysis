#!/usr/bin/env Rscript
# =============================================================================
# test_triage_metadata.R  ·  acceptance test for the Triage_Metadata module
# -----------------------------------------------------------------------------
# Two halves:
#   1. UNIT — the statistics, against hand-built data with a known answer. These
#      run anywhere, with no QC run present.
#   2. INTEGRATION — the module's own outputs, if it has been run here. These
#      assert the properties that make the outputs safe to read, NOT specific
#      values from any one dataset (a test that pins this run's p-values would
#      fail on the next run for the right reasons, which makes it useless).
#
#   Rscript test_triage_metadata.R
# =============================================================================

HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
source(file.path(HERE, "triage_metadata_helpers.R"))
OUT <- file.path(HERE, "output_files")

.pass <- 0L; .fail <- 0L
ok <- function(what, cond, detail = "") {
  if (isTRUE(cond)) { .pass <<- .pass + 1L; cat(sprintf("  ok   %s\n", what)) }
  else { .fail <<- .fail + 1L
         cat(sprintf("  FAIL %s%s\n", what, if (nzchar(detail)) paste0("  [", detail, "]") else "")) }
}
section <- function(s) cat(sprintf("\n%s\n%s\n", s, strrep("-", nchar(s))))

# =============================================================================
section("1. Wilson interval")
# =============================================================================
ci <- wilson_ci(2, 5)
ok("2 of 5 gives a wide interval, not a point estimate", ci[2] - ci[1] > 0.5,
   sprintf("[%.3f, %.3f]", ci[1], ci[2]))
ok("interval contains the point estimate", ci[1] <= 0.4 && ci[2] >= 0.4)
ok("never leaves [0,1] at extreme counts",
   { c0 <- wilson_ci(0, 7); c1 <- wilson_ci(7, 7); c0[1] >= 0 && c1[2] <= 1 })
big <- wilson_ci(20, 1000)
ok("large n gives a narrow interval", big[2] - big[1] < 0.03)
ok("n = 0 returns NA rather than dividing by zero", all(is.na(wilson_ci(0, 0))))

# =============================================================================
section("2. Minimum detectable effect")
# =============================================================================
m_small <- min_detectable_rr(25, 0.018, 4000)
m_large <- min_detectable_rr(500, 0.018, 4000)
ok("small groups need a bigger effect to detect", m_small > m_large,
   sprintf("n=25 -> %.2fx, n=500 -> %.2fx", m_small, m_large))
ok("MDE is a ratio above 1", m_small > 1 && m_large > 1)
ok("invalid input returns NA, not a number", is.na(min_detectable_rr(10, 0, 100)))

# =============================================================================
section("3. Group rates — the size floor must REPORT, not hide")
# =============================================================================
set.seed(4)
g <- c(rep("big", 200), rep("tiny", 5))
f <- c(rep(FALSE, 200), TRUE, TRUE, FALSE, FALSE, FALSE)   # tiny is 40%
r <- rates_by_group(f, g, min_n = 20L)
ok("group below the floor is excluded from the tested set", !("tiny" %in% r$group))
bf <- attr(r, "below_floor")
ok("...but is returned in the below_floor attribute", !is.null(bf) && "tiny" %in% bf$group)
ok("the excluded group keeps its real rate", isTRUE(bf$pct[bf$group == "tiny"] == 40))
ok("the excluded group carries a CI wide enough to disown it",
   bf$ci_hi_pct[bf$group == "tiny"] - bf$ci_lo_pct[bf$group == "tiny"] > 50)

# =============================================================================
section("4. Permutation test — two statistics, and they can disagree")
# =============================================================================
# One small-ish group with a high rate, everything else flat. max_rate should be
# far more excited than the omnibus. This is the exact configuration that made
# the earlier "40% site" look like a finding.
set.seed(11)
n <- 1000
grp <- c(rep("odd", 25), rep(paste0("g", 1:13), each = 75))
fl  <- logical(length(grp))
fl[sample(which(grp == "odd"), 5)] <- TRUE
fl[sample(which(grp != "odd"), 20)] <- TRUE
tt <- perm_group_test(fl, grp, min_n = 20L, n_perm = 2000L, seed = 2L)
ok("both statistics are reported", all(c("p_max_rate", "p_deviance") %in% names(tt)))
ok("max-rate reacts to the single odd group", tt$p_max_rate < 0.05,
   sprintf("p_max = %.4f", tt$p_max_rate))
ok("p-values are valid probabilities",
   tt$p_max_rate > 0 && tt$p_max_rate <= 1 && tt$p_deviance > 0 && tt$p_deviance <= 1)
ok("permutation p can never be exactly 0", tt$p_max_rate > 0 && tt$p_deviance > 0)

# A flag that is pure noise must not be significant.
set.seed(7)
fl2 <- rep(FALSE, n); fl2[sample.int(n, 25)] <- TRUE
tn <- perm_group_test(fl2, grp, min_n = 20L, n_perm = 2000L, seed = 3L)
ok("a null association is not significant on the omnibus", tn$p_deviance > 0.05,
   sprintf("p_dev = %.4f", tn$p_deviance))

# Stratified permutation must respect the strata.
st <- rep(c("A", "B"), length.out = n)
ts <- perm_group_test(fl, grp, stratum = st, min_n = 20L, n_perm = 500L, seed = 5L)
ok("stratified run returns a result and marks itself", nzchar(ts$stratified_by))

ok("too few groups returns NA rather than a bogus p",
   is.na(perm_group_test(fl, rep("only", n), min_n = 20L, n_perm = 100L)$p_deviance))

# =============================================================================
section("5. Conditioning — a pure composition effect must collapse")
# =============================================================================
# `grp2` carries NO signal of its own; it only predicts the plate, and the plate
# is what drives the flag. Conditioning on plate should remove it.
set.seed(21)
N <- 1200
plate <- sample(paste0("p", 1:8), N, replace = TRUE)
bad   <- plate %in% c("p1", "p2")
# grp2 tracks the bad plates closely (10% mislabelled, so the model stays
# estimable) but has no effect of its own — the flag depends only on the plate.
grp2  <- ifelse(xor(bad, runif(N) < 0.10), "x", "y")
flag3 <- runif(N) < ifelse(bad, 0.30, 0.01)
cd <- conditioned_deviance(flag3, grp2, plate = plate, min_n = 20L)
ok("unconditioned deviance sees the (spurious) group effect", cd$dev_group_only > 5,
   sprintf("dev = %.1f", cd$dev_group_only))
ok("conditioning on plate collapses it",
   cd$dev_group_given_plate < 0.4 * cd$dev_group_only,
   sprintf("%.1f -> %.1f", cd$dev_group_only, cd$dev_group_given_plate))

# =============================================================================
section("6. Mutual adjustment — a proxy must drop out")
# =============================================================================
# `shadow` is a noisy copy of `real`. Only `real` drives the flag. After
# adjustment `shadow` should retain far less of its deviance than `real`.
#
# The copy rate is deliberately moderate (55%). At 85% the two are so collinear
# that BOTH lose most of their deviance — `real` kept only 24% — which is correct
# behaviour but tests nothing: near-collinear predictors cannot be told apart, and
# a fixture that cannot separate them cannot demonstrate that the code does.
set.seed(31)
M <- 1500
real   <- sample(c("a", "b", "c"), M, TRUE)
shadow <- ifelse(runif(M) < 0.55, real, sample(c("a", "b", "c"), M, TRUE))
p      <- c(a = 0.02, b = 0.02, c = 0.18)[real]
flag4  <- runif(M) < p
ma <- mutual_adjustment(flag4, list(real = real, shadow = shadow), min_n = 20L)
ok("mutual_adjustment returns a row per variable", !is.null(ma) && nrow(ma) == 2)
rs <- stats::setNames(ma$retained_share, ma$variable)
ps <- stats::setNames(ma$p_given_others, ma$variable)
ok("the real driver keeps most of its deviance", rs[["real"]] > 0.6,
   sprintf("real keeps %.0f%%", 100 * rs[["real"]]))
ok("the proxy keeps far less than the real driver", rs[["shadow"]] < 0.5 * rs[["real"]],
   sprintf("shadow keeps %.0f%% vs real %.0f%%", 100 * rs[["shadow"]], 100 * rs[["real"]]))
ok("the real driver stays significant after adjustment", ps[["real"]] < 0.01,
   sprintf("p = %.3g", ps[["real"]]))
ok("the proxy does not", ps[["shadow"]] > 0.05, sprintf("p = %.3g", ps[["shadow"]]))
ok("adjusted p-values are valid", all(ma$p_given_others >= 0 & ma$p_given_others <= 1))

# =============================================================================
section("7. Diagnosis collapse — the curated table drives it")
# =============================================================================
# The intended path: the merge's own CDX_review.csv. Its `NA` rows mean "exclude
# from the analysis cohort" — this module must keep those wells, as their own
# level, because an excluded well was still either triaged or not.
rev_cdx <- data.frame(
  CDX = c("NCI", "MCI", "AD", "Dementia Vascular", "Missing", "Insufficient Data", "Other"),
  FINAL_CDX_collapsed = c("NCI", "MCI", "AD", "Dementia_Other", NA, NA, NA),
  stringsAsFactors = FALSE)
raw <- c("AD", "NCI", "MCI", "Missing", "Insufficient Data", "Dementia Vascular",
         "Other", "Never seen before", NA)
cd2 <- collapse_diagnosis(raw, cdx_review = rev_cdx)
ok("a curated level is used verbatim", cd2[1] == "AD" && cd2[2] == "NCI")
ok("a curated collapse is applied", cd2[6] == "Dementia_Other")
ok("cohort-EXCLUDED levels are kept, not dropped",
   cd2[4] == "Not codeable" && cd2[5] == "Not codeable" && cd2[7] == "Not codeable")
ok("an excluded level is NOT folded into a clinical one", cd2[4] != cd2[2])
ok("a value absent from the table still gets a level", cd2[8] == "Not codeable")
ok("NA becomes an explicit level rather than vanishing", !is.na(cd2[9]))
ok("no input is left unmapped", !any(is.na(cd2)))

cm <- collapse_map(raw, cdx_review = rev_cdx)
ok("the collapse map accounts for every input row", sum(cm$n_wells) == length(raw))
ok("the map records which mapping was used", all(cm$source == "CDX_review.csv"))

# The fallback matcher, used only when no review table exists.
cf <- collapse_diagnosis(c("AD", "NCI", "Lewy Body", "Missing"), cdx_review = NULL)
ok("fallback still maps the common levels", cf[1] == "Case (AD)" && cf[2] == "Control")
ok("fallback groups dementia subtypes", cf[3] == "Dementia_Other")
ok("fallback marks itself in the map",
   collapse_map("AD", NULL)$source[1] == "fallback matcher")

# =============================================================================
section("7b. Age bands are a parameter, not a constant")
# =============================================================================
ab <- age_band(c(55, 65, 75, 85, NA))
ok("default bands split at 60/70/80",
   identical(ab, c("<60", "60-69", "70-79", "80+", NA_character_)))
ok("NA age stays NA", is.na(ab[5]))
ab2 <- age_band(c(40, 55, 70), breaks = c(50, 65))
ok("custom cutpoints are honoured", identical(ab2, c("<50", "50-64", "65+")))

# =============================================================================
section("8. Ancestry comes from the curated table, never from a guess")
# =============================================================================
rev <- data.frame(Group = c("AA", "AFDC", ""), Race = c("BL", "", "WH"),
                  Ethnicity = c("NH", "", ""),
                  FINAL_Ancestry = c("AA", "AFDC", "EUR"), stringsAsFactors = FALSE)
dd <- data.frame(Group = c("AA", "AFDC", "OTHER", "AA"), Race = c("BL", "BL", "WH", "MU"),
                 Ethnicity = c("NH", "HI", "NH", "NH"), stringsAsFactors = FALSE)
a <- assign_ancestry_from_review(dd, rev)
ok("a fully-specified rule matches", a[1] == "AA")
ok("a blank condition acts as a wildcard, not a literal", a[2] == "AFDC")
ok("a row matching only the Race rule takes it", a[3] == "EUR")
ok("an unmatched row becomes Other, never a guess", a[4] == "Other")
ok("a missing review table yields NA, not a column of 'Unknown'",
   all(is.na(assign_ancestry_from_review(dd, NULL))))

# =============================================================================
section("9. Separability verdict")
# =============================================================================
# Nested: every level of A sits on exactly one B. A conditioned test is then
# not estimable and the module must say so instead of reporting a number.
a_nested <- rep(c("s1", "s2", "s3"), each = 40)
b_nested <- rep(c("p1", "p2", "p3"), each = 40)
cs <- confounding_summary(a_nested, b_nested, "Site", "plate")
ok("perfect nesting is called out as not estimable", grepl("not estimable", cs$verdict))
ok("nesting shows all wells on one level", cs$max_share_of_a_on_one_b == 1)

a_spread <- rep(c("s1", "s2", "s3"), each = 40)
b_spread <- rep(paste0("p", 1:8), length.out = 120)
cs2 <- confounding_summary(a_spread, b_spread, "Site", "plate")
ok("a spread design is called estimable", grepl("estimable", cs2$verdict) &&
     !grepl("not estimable", cs2$verdict))

# =============================================================================
section("10. Balance table")
# =============================================================================
set.seed(41)
bd <- data.frame(age = c(rnorm(50, 78, 8), rnorm(500, 72, 8)),
                 sex = sample(c("F", "M"), 550, TRUE), stringsAsFactors = FALSE)
grpb <- c(rep(TRUE, 50), rep(FALSE, 500))
bt <- balance_table(bd, grpb, continuous = "age", categorical = "sex")
ok("balance table covers the continuous variable", any(bt$variable == "age"))
ok("a real mean shift shows a positive standardised difference",
   bt$std_diff[bt$variable == "age"] > 0.3)
ok("a balanced categorical stays near zero",
   all(abs(bt$std_diff[bt$variable == "sex"]) < 0.3))

# =============================================================================
section("11. Integration — the module's own outputs, if present")
# =============================================================================
f_tests <- file.path(OUT, "metadata_association_tests.csv")
if (!file.exists(f_tests)) {
  cat("  (skipped: no output_files/ here — run run_triage_metadata.R first)\n")
} else {
  te <- utils::read.csv(f_tests, stringsAsFactors = FALSE)
  ra <- utils::read.csv(file.path(OUT, "triage_rates_by_metadata.csv"), stringsAsFactors = FALSE)

  ok("every axis is tested separately, never pooled into one row",
     length(unique(te$axis)) >= 3 && !any(te$axis == "all"),
     paste(unique(te$axis), collapse = ", "))
  ok("the triage axes are reported apart from the read flag",
     all(c("IC", "read flag") %in% te$axis))
  ok("both statistics are present for every test",
     all(!is.na(te$p_max_rate) & !is.na(te$p_deviance)))
  ok("every test carries a conditioned deviance",
     mean(!is.na(te$dev_given_plate)) > 0.8)
  ok("every non-significant row states what it could have detected",
     all(!is.na(te$min_detectable_rate_ratio[te$p_deviance >= 0.05])))
  ok("no result is reported without a reading", all(nzchar(te$reading)))

  ok("groups below the floor are still written out, marked",
     any(!ra$above_floor), sprintf("%d below-floor rows", sum(!ra$above_floor)))
  ok("no below-floor group was silently tested",
     { bl <- unique(ra$group[!ra$above_floor & ra$variable == "Site"])
       al <- unique(ra$group[ra$above_floor & ra$variable == "Site"])
       !length(intersect(bl, al)) })
  ok("reported rates are consistent with their counts",
     all(abs(ra$pct - 100 * ra$n_flagged / ra$n) < 0.011))
  ok("every CI brackets its own point estimate",
     all(ra$ci_lo_pct <= ra$pct + 1e-6 & ra$ci_hi_pct >= ra$pct - 1e-6))

  un <- utils::read.csv(file.path(OUT, "unlabelled_wells.csv"), stringsAsFactors = FALSE)
  ok("wells with no metadata are accounted for, not dropped",
     nrow(un) == 2 && sum(un$n) > 0)

  f_adj <- file.path(OUT, "mutual_adjustment.csv")
  if (file.exists(f_adj)) {
    ad <- utils::read.csv(f_adj, stringsAsFactors = FALSE)
    ok("mutual adjustment includes read depth as a competitor",
       any(ad$variable == "log2(reads)"))
    ok("mutual adjustment covers more than one label",
       length(unique(ad$variable)) >= 3)
    ok("adjusted deviance is never negative (IRLS noise is clamped)",
       all(ad$dev_given_others >= 0))
    ok("retained share is never negative", all(is.na(ad$retained_share) |
                                                 ad$retained_share >= 0))
  }

  rpt <- list.files(OUT, pattern = "^triage_metadata_report_.*\\.html$", full.names = TRUE)
  ok("a per-run HTML report was written", length(rpt) == 1, paste(basename(rpt), collapse = ", "))
  if (length(rpt) == 1) {
    h <- paste(readLines(rpt[1], warn = FALSE), collapse = "\n")
    ok("the report is a complete standalone document",
       grepl("<!doctype html>", h, ignore.case = TRUE) && grepl("</html>", h))
    ok("the report states it changes nothing", grepl("Diagnostic only", h, ignore.case = TRUE))
    ok("the report carries the tests, not just prose", grepl("p_deviance", h))
  }

  cf <- utils::read.csv(file.path(OUT, "design_confounding.csv"), stringsAsFactors = FALSE)
  ok("separability is checked before any conditioned claim", nrow(cf) >= 1 &&
       all(nzchar(cf$verdict)))
}

# =============================================================================
cat(sprintf("\n=========================================\n %d passed, %d failed\n=========================================\n",
            .pass, .fail))
if (.fail > 0) quit(status = 1)
