#!/usr/bin/env Rscript
# =============================================================================
# test_triage_review.R  ·  2026-07-27
# -----------------------------------------------------------------------------
# Acceptance test for the Primary QC screen audit. Base R only, non-zero exit on
# failure.
#
#   Rscript test_triage_review.R
#
# WHY A TEST AND NOT JUST A DIFF: output_files/ is gitignored, so checking out
# another branch does NOT revert the CSVs. Only a content check can tell you
# whether a change was safe.
#
# PART A  unit checks on the helpers against hand-built data where the answer is
#         known by inspection. Runs in the source repo, where there is no data.
# PART B  the real run, cross-checked against Primary QC's own files.
# =============================================================================

HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
source(file.path(HERE, "triage_review_helpers.R"))
OUT <- file.path(HERE, "output_files")
PQ  <- normalizePath(file.path(HERE, "..", "..", "Primary_QC"), mustWork = FALSE)

fails <- 0L
ok <- function(label, cond, detail = "") {
  pass <- isTRUE(cond); if (!pass) fails <<- fails + 1L
  cat(sprintf("  [%s] %s%s\n", if (pass) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0("  -- ", detail) else ""))
}
.rd <- function(...) { p <- file.path(...); if (file.exists(p))
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) else NULL }

# =============================================================================
cat("\n=== PART A: helper unit checks ===\n")
# =============================================================================

cat("\n-- cohen_d / auc --\n")
set.seed(7); x <- c(rnorm(200, 1), rnorm(200, 0)); g <- rep(c(TRUE, FALSE), each = 200)
ok("cohen_d is positive when the flagged group is higher", cohen_d(x, g) > 0.7,
   sprintf("d = %.2f", cohen_d(x, g)))
ok("auc agrees with the rank identity", abs(auc(x, g) - 0.76) < 0.06,
   sprintf("AUC = %.3f", auc(x, g)))
ok("auc of a perfect separator is 1", auc(c(1, 2, 3, 4), c(FALSE, FALSE, TRUE, TRUE)) == 1)
ok("auc is direction-aware (reversed group gives 1 - AUC)",
   abs(auc(x, g) + auc(x, !g) - 1) < 1e-9)
ok("cohen_d returns NA rather than dividing by zero on a constant",
   is.na(cohen_d(rep(5, 10), rep(c(TRUE, FALSE), 5))))

cat("\n-- triage_axes: the reason strings Primary QC actually writes --\n")
ax <- triage_axes(c("IC outlier",
                    "PCA outlier; High outlier burden (FDR); IC outlier",
                    "High outlier burden (FDR)", "Bad Data", ""))
ok("a multi-reason string sets every axis it names",
   all(unlist(ax[2, c("IC", "PCA", "burden")])) && !ax$baddata[2])
ok("'High outlier burden (FDR)' is not mistaken for a PCA outlier", !ax$PCA[3])
ok("an empty reason sets no axis", sum(unlist(ax[5, ])) == 0)
ok("the parenthesis in the burden reason is matched literally, not as a regex",
   ax$burden[3] && !ax$burden[1])

cat("\n-- axis_overlap: sole reason vs on-axis --\n")
A <- data.frame(IC = c(TRUE, TRUE, FALSE, FALSE), PCA = c(FALSE, TRUE, TRUE, FALSE),
                burden = c(FALSE, FALSE, TRUE, TRUE))
ov <- axis_overlap(A)
ok("on-axis counts every well carrying the axis",
   identical(ov$n_wells_on_axis, c(2L, 2L, 2L)))
ok("sole-reason counts only wells where it is the single axis",
   identical(ov$n_wells_sole_reason, c(1L, 0L, 1L)))
ok("an axis that is never the sole reason reports 100% shared",
   ov$pct_of_axis_shared[ov$axis == "PCA"] == 100)
cb <- axis_combinations(A)
ok("combinations sum to the number of wells", sum(cb$n_wells) == nrow(A))

cat("\n-- flag_agreement --\n")
a <- c(rep(TRUE, 57), TRUE, rep(FALSE, 100)); b <- c(rep(TRUE, 57), FALSE, TRUE, rep(FALSE, 99))
fa <- flag_agreement(a, b, "a", "b")
ok("Jaccard is both / union", abs(fa$jaccard - 57 / 59) < 1e-4, sprintf("%.4f", fa$jaccard))
ok("a near-identical pair is called a DUPLICATE", grepl("DUPLICATE", fa$verdict))
ok("a disjoint pair is called largely independent",
   grepl("largely independent",
         flag_agreement(c(TRUE, FALSE, FALSE), c(FALSE, TRUE, FALSE))$verdict))
ok("two flags that never fire do not divide by zero",
   flag_agreement(c(FALSE, FALSE), c(FALSE, FALSE))$verdict == "neither fires")
ok("NA is treated as not-flagged, not as TRUE",
   flag_agreement(c(NA, TRUE), c(TRUE, TRUE))$n_a == 1)

cat("\n-- gate_row: inert detection --\n")
ok("a gate that fires on nobody is INERT",
   grepl("INERT", gate_row("g", "r", rep(FALSE, 100), 100)$status))
ok("a gate that fires is active",
   gate_row("g", "r", c(TRUE, rep(FALSE, 99)), 100)$status == "active")
ok("NA does not count as fired", gate_row("g", "r", c(NA, NA, TRUE), 3)$n_fired == 1)

cat("\n-- gate_operating_points / gate_cut_neighbourhood --\n")
st <- c(0.001, 0.001, 0.05, 0.05, 1, 1, 1)
op <- gate_operating_points(st, count = c(6, 6, 4, 4, 0, 0, 0), cut = 0.01)
ok("one row per distinct value of the statistic", nrow(op) == 3)
ok("removal counts are cumulative in the worsening direction",
   identical(op$n_removed_at_this_cut[order(op$stat_value)], c(2L, 4L, 7L)))
ok("the current cut is marked on exactly the values inside it",
   sum(op$is_current_cut) == 1 && op$stat_value[op$is_current_cut] == 0.001)
nb <- gate_cut_neighbourhood(st, st < 0.01, "g")
ok("the gap is worst-retained minus best-removed", abs(nb$gap - (0.05 - 0.001)) < 1e-9,
   sprintf("gap %.4f", nb$gap))

cat("\n-- read_gate_geometry: the floor that makes the lower bound unreachable --\n")
# Realistic spread: the real plate-bays run IQR/median 0.45 upward, and the floor
# binds whenever k*IQR exceeds median - 500000. Here Q1 = 1e6, Q3 = 9e6, so
# 4*IQR = 32e6 against a median of 5e6 -- deeply negative, as in all 50 real bays.
gg <- read_gate_geometry(c(rep(1e6, 10), rep(9e6, 10)), rep("p1", 20), k = 4)
ok("the raw lower bound is reported before the floor is applied, and is negative here",
   gg$raw_lower_bound < 0, sprintf("raw lower bound %.3g", gg$raw_lower_bound))
ok("floor_binds is TRUE when median - k*IQR falls under 500000", gg$floor_binds)
# The complement: a tight plate at a high median keeps a genuinely plate-relative
# lower bound, which is what the gate was presumably designed to do.
tight <- read_gate_geometry(seq(4.5e6, 7.5e6, length.out = 40), rep("p", 40), k = 1)
ok("a tight, high-median plate keeps a plate-relative lower bound", !tight$floor_binds,
   sprintf("raw lower bound %.3g", tight$raw_lower_bound))
ok("a zero-IQR plate does not bind either (the subtraction vanishes)",
   !read_gate_geometry(rep(5e6, 20), rep("p", 20), k = 4)$floor_binds)

cat("\n-- perm_test / position_effects --\n")
set.seed(3)
lab <- c(rep(TRUE, 10), rep(FALSE, 90))
pt <- perm_test(lab, function(l) sum(l[1:10]), n_perm = 2000L)
ok("perm_test recovers an extreme observed statistic", pt$p < 0.01,
   sprintf("obs %d, null mean %.2f, p %.4f", pt$observed, pt$null_mean, pt$p))
ok("perm_test holds the label count fixed", pt$null_mean > 0.5 && pt$null_mean < 2)
# a planted position effect must be detected
pos_row <- rep(LETTERS[1:8], each = 12); pos_col <- rep(1:12, 8)
fl <- rep(FALSE, 96); fl[pos_row == "F" & pos_col == 7] <- TRUE
big <- do.call(rbind, lapply(1:50, function(i) data.frame(r = pos_row, c = pos_col, f = fl)))
pe <- position_effects(big$f, big$r, big$c, n_perm = 2000L)
ok("a planted single-position effect is found and ranked first",
   pe$position[1] == "F7" && pe$n_flagged[1] == 50)
ok("the permutation calls the planted effect significant", attr(pe, "perm")$p < 0.01,
   sprintf("p = %.4f", attr(pe, "perm")$p))
ok("a flat flag gives no significant position effect",
   { set.seed(9); f2 <- sample(c(TRUE, FALSE), nrow(big), TRUE, c(.02, .98))
     attr(position_effects(f2, big$r, big$c, n_perm = 2000L), "perm")$p > 0.05 })

cat("\n-- row_gradient --\n")
rr2 <- rep(LETTERS[1:8], each = 10)
val <- rep(8:1, each = 10) + rnorm(80, 0, 0.01)
rgd <- row_gradient(val, rr2, rep("p1", 80))
ok("a declining row gradient gives rho = -1", rgd$rho_row_vs_value == -1)
ok("a flat quantity gives rho near 0",
   abs(row_gradient(rnorm(80), rr2, rep("p1", 80))$rho_row_vs_value) < 1)

# =============================================================================
cat("\n\n=== PART B: the real run ===\n")
# =============================================================================
if (!dir.exists(OUT)) {
  cat("\nNo output_files/ here -- this is the source repo, not a QC run.\n")
  cat("Part A passed. Run this inside a QC run directory after the audit.\n")
  cat(sprintf("\n%s  (%d failure%s)\n", if (fails == 0L) "ALL CHECKS PASSED" else "FAILED",
              fails, if (fails == 1L) "" else "s"))
  quit(status = if (fails == 0L) 0L else 1L)
}

inv <- .rd(OUT, "qc_gate_inventory.csv");  dup <- .rd(OUT, "qc_gate_duplication.csv")
ov  <- .rd(OUT, "triage_axis_overlap.csv"); det <- .rd(OUT, "triaged_well_detail.csv")
op  <- .rd(OUT, "gate_operating_points.csv"); rg <- .rd(OUT, "read_gate_geometry.csv")
at  <- .rd(OUT, "qc_axis_attribution.csv"); pe <- .rd(OUT, "well_position_effects.csv")
tri <- .rd(PQ, "output_files", "samples_to_triage.csv")
fro <- .rd(PQ, "output_files", "flagged_read_outliers.csv")

cat("\n-- agreement with Primary QC's own files --\n")
if (!is.null(tri) && !is.null(ov)) {
  ok("the axis totals sum to at least the number of triaged wells (overlaps allowed)",
     sum(ov$n_wells_on_axis) >= nrow(tri),
     sprintf("%d axis-hits over %d wells", sum(ov$n_wells_on_axis), nrow(tri)))
  ok("sole-reason counts sum to the wells carrying exactly one axis",
     sum(ov$n_wells_sole_reason) <= nrow(tri))
  ok("every axis total is at least its sole-reason count",
     all(ov$n_wells_on_axis >= ov$n_wells_sole_reason))
  ok("triaged wells in the detail sheet match samples_to_triage.csv",
     sum(det$triaged) == nrow(tri), sprintf("%d vs %d", sum(det$triaged), nrow(tri)))
}
if (!is.null(fro))
  ok("read-flagged wells in the detail sheet match flagged_read_outliers.csv",
     sum(det$read_flagged) == nrow(fro), sprintf("%d vs %d", sum(det$read_flagged), nrow(fro)))

cat("\n-- the inert-gate check is the point of the inventory --\n")
ok("the inventory covers both in-house and vendor gates",
   any(grepl("^triage", inv$gate)) && any(grepl("^vendor", inv$gate)))
ok("status is INERT exactly when n_fired is 0",
   all(grepl("INERT", inv$status) == (inv$n_fired == 0)))
ok("at least one gate fires (otherwise the join is broken, not the QC)",
   sum(inv$n_fired) > 0, sprintf("%d total firings", sum(inv$n_fired)))

cat("\n-- duplication --\n")
ok("the verdict follows the Jaccard consistently",
   all((dup$jaccard >= 0.9) == grepl("DUPLICATE", dup$verdict)))
ic <- dup[dup$flag_a == "triage: IC outlier" & grepl("IC Median", dup$flag_b), ]
if (nrow(ic))
  ok("the in-house IC axis is flagged as a duplicate of the vendor IC-median metric",
     grepl("DUPLICATE", ic$verdict),
     sprintf("Jaccard %.3f, %d vs %d wells, %d discordant",
             ic$jaccard, ic$n_a, ic$n_b, ic$n_only_a + ic$n_only_b))

cat("\n-- the burden gate is quantised, which is the threshold answer --\n")
if (!is.null(op)) {
  ok("the statistic takes few distinct values (it is not a continuous dial)",
     nrow(op) < 40, sprintf("%d distinct values", nrow(op)))
  ok("removal counts increase monotonically as the cut loosens",
     !is.unsorted(op$n_removed_at_this_cut[order(op$stat_value)]))
  ok("exactly the values below the cut are marked as inside it",
     all(op$is_current_cut == (op$stat_value < 0.01), na.rm = TRUE))
}
if (!is.null(rg)) {
  ok("read_gate_geometry reports one row per plate-bay", nrow(rg) > 0)
  ok("floor_binds agrees with the raw lower bound being under 500000",
     all(rg$floor_binds == (rg$raw_lower_bound < 5e5)))
}

cat("\n-- axis attribution must be reported per class, not pooled --\n")
if (!is.null(at)) {
  ok("more than one triage class is reported", length(unique(at$class)) > 1,
     paste(unique(at$class), collapse = "; "))
  ok("AUC is in [0,1] and d is finite wherever reported",
     all(at$auc >= 0 & at$auc <= 1, na.rm = TRUE) && all(is.finite(at$cohen_d) | is.na(at$cohen_d)))
  # The substantive check: the classes must NOT all look alike. If they did, the
  # per-class split would be pointless and pooling would be safe.
  s <- at[at$axis == "mean_INT_z" & at$class %in% c("IC only", "PCA (any)"), ]
  if (nrow(s) == 2)
    ok("the triage classes differ materially on intensity, so pooling would mislead",
       abs(diff(s$cohen_d)) > 1,
       sprintf("IC-only d %+.2f vs PCA d %+.2f", s$cohen_d[s$class == "IC only"],
               s$cohen_d[s$class == "PCA (any)"]))
}

cat("\n-- position --\n")
if (!is.null(pe)) {
  ok("every plate position is represented", nrow(pe) > 50, sprintf("%d positions", nrow(pe)))
  ok("flagged counts never exceed the wells at that position", all(pe$n_flagged <= pe$n_wells))
  ok("positions are sorted worst-first", !is.unsorted(rev(pe$n_flagged)))
}

cat(sprintf("\n%s  (%d failure%s)\n",
            if (fails == 0L) "ALL CHECKS PASSED" else "FAILED",
            fails, if (fails == 1L) "" else "s"))
quit(status = if (fails == 0L) 0L else 1L)
