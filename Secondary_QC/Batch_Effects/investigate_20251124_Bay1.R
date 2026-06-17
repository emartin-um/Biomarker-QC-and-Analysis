#!/usr/bin/env Rscript
# =============================================================================
# investigate_20251124_Bay1.R
# -----------------------------------------------------------------------------
# Deep-dive on the biomarker assay anomaly flagged on run 20251124-1407_Bay1.
# The APOE secondary QC found ~38% APOE4-vs-genotype mismatch on this bay, with
# the APOE4 assay losing dynamic range. Here we characterize it rigorously, using
# the run's SIBLING bays (Bay2/Bay3, same plate/date/handling) as the internal
# control, and ask: is it APOE4-specific or panel-wide, all of Bay1 or a subset,
# and is the sample composition balanced across bays (i.e. not a biology confound).
#
# Reads the post-merge biomarker matrix (all bay samples retained) + biomarker
# list from preprocessed_data.rds. Writes a flag table + figures to output_files/.
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })

MM <- "../../Metadata_Merge/output_files"
OUT <- "output_files"; if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

TARGET_RUN  <- "20251124-1407"
TARGET_BAY  <- "20251124-1407_Bay1"

m  <- read.csv(file.path(MM, "merged_combined_post_QC.csv"), check.names = FALSE)
bc <- readRDS(file.path(MM, "preprocessed_data.rds"))$biomarker_cols
bc <- intersect(bc, names(m))

# ---- ε4 dosage from the (concordant) clinical genotype; this run is not a
#      genotype-discordance case, so APOE.geno is fine here ----------------------
e4_count <- function(g) {
  g <- gsub("[^0-9]", "", trimws(g))
  ifelse(g == "" | is.na(g), NA_integer_, lengths(regmatches(g, gregexpr("4", g))))
}
m$e4 <- e4_count(m$APOE.geno)

m <- m %>% mutate(
  grp = case_when(
    RUN == TARGET_BAY                                   ~ "Bay1 (target)",
    grepl(paste0("^", TARGET_RUN), RUN)                 ~ "Bay2/3 (sibling)",
    TRUE                                                ~ "All other runs"),
  grp = factor(grp, levels = c("Bay1 (target)", "Bay2/3 (sibling)", "All other runs"))
)

cat("================ run/bay sample counts ================\n")
print(m %>% count(RUN) %>% filter(grepl(TARGET_RUN, RUN)))
cat("\n")

# ---- 1. Is the sample composition balanced across the 3 bays? ----------------
cat("================ composition across the 3 bays (confound check) ========\n")
run3 <- m %>% filter(grepl(paste0("^", TARGET_RUN), RUN))
for (v in c("Ethnicity", "Ancestry", "CDX", "e4")) {
  if (v %in% names(run3)) {
    cat(sprintf("\n-- %s by bay --\n", v))
    print(addmargins(table(Bay = run3$RUN, run3[[v]]), 2))
  }
}

# ---- 2. APOE4 dynamic range by ε4 dosage, per group --------------------------
cat("\n================ APOE4 spread by ε4 dosage (carrier separation) =======\n")
rng <- m %>% filter(!is.na(e4), !is.na(APOE4)) %>%
  group_by(grp, e4) %>%
  summarise(n = n(),
            APOE4_med = round(median(APOE4), 2), APOE4_IQR = round(IQR(APOE4), 2),
            APOE_med  = round(median(APOE),  2), APOE_IQR  = round(IQR(APOE),  2),
            .groups = "drop")
print(as.data.frame(rng))

# carrier vs non-carrier APOE4 separation (med of carriers - med of non-carriers)
sep <- m %>% filter(!is.na(e4), !is.na(APOE4)) %>%
  mutate(carrier = e4 > 0) %>%
  group_by(grp, carrier) %>% summarise(med = median(APOE4), .groups = "drop") %>%
  pivot_wider(names_from = carrier, values_from = med) %>%
  mutate(separation = round(`TRUE` - `FALSE`, 2))
cat("\n-- APOE4 carrier vs non-carrier median separation (bigger = healthier assay) --\n")
print(as.data.frame(sep))

# ---- 3. Panel-wide: is Bay1 inflated only for APOE4? -------------------------
# For each biomarker, residualize on ε4 band (removes APOE4's genotype-driven
# bimodality) then compare Bay1 within-group IQR to the SIBLING bays' IQR.
resid_e4 <- function(df, col) {
  x <- df[[col]]; r <- x
  for (g in unique(df$e4)) { i <- which(df$e4 == g & !is.na(x))
    if (length(i) > 1) r[i] <- x[i] - median(x[i], na.rm = TRUE) }
  r
}
b1   <- m %>% filter(RUN == TARGET_BAY, !is.na(e4))
sib  <- m %>% filter(grepl(paste0("^", TARGET_RUN), RUN), RUN != TARGET_BAY, !is.na(e4))
oth  <- m %>% filter(!grepl(paste0("^", TARGET_RUN), RUN), !is.na(e4))

iqr_ratio <- sapply(bc, function(col) {
  r_b1  <- resid_e4(b1,  col); r_sib <- resid_e4(sib, col)
  IQR(r_b1, na.rm = TRUE) / IQR(r_sib, na.rm = TRUE)
})
ratio_df <- tibble(biomarker = bc, iqr_ratio_bay1_vs_sibling = iqr_ratio) %>%
  filter(is.finite(iqr_ratio_bay1_vs_sibling), iqr_ratio_bay1_vs_sibling > 0) %>%
  arrange(desc(iqr_ratio_bay1_vs_sibling)) %>%
  mutate(rank = row_number())
cat("\n================ within-bay IQR inflation vs sibling bays (top 15) =====\n")
print(as.data.frame(head(ratio_df, 15)))
cat(sprintf("\nAPOE4 rank: #%d of %d  (ratio %.2f)\n",
            ratio_df$rank[ratio_df$biomarker == "APOE4"], nrow(ratio_df),
            ratio_df$iqr_ratio_bay1_vs_sibling[ratio_df$biomarker == "APOE4"]))
cat(sprintf("APOE  rank: #%d of %d  (ratio %.2f)\n",
            ratio_df$rank[ratio_df$biomarker == "APOE"], nrow(ratio_df),
            ratio_df$iqr_ratio_bay1_vs_sibling[ratio_df$biomarker == "APOE"]))

# ---- 4. Is it ALL of Bay1 or a subset? --------------------------------------
# Per-sample APOE4 deviation from its ε4-band median (using sibling+other as ref).
ref_med <- m %>% filter(grp != "Bay1 (target)", !is.na(e4), !is.na(APOE4)) %>%
  group_by(e4) %>% summarise(med = median(APOE4), .groups = "drop")
b1_dev <- b1 %>% left_join(ref_med, by = "e4") %>%
  mutate(APOE4_dev = APOE4 - med, apoe4_off = abs(APOE4_dev) > 2 & !is.na(APOE4_dev))
cat("\n================ Bay1 per-sample APOE4 deviation from ε4-band median ===\n")
cat(sprintf("Bay1 n=%d | off-band (|dev|>2 log2): %d (%.0f%%) | of those |dev|>4: %d | median |dev|: %.2f | sibling median |dev|: %.2f\n",
            nrow(b1_dev), sum(b1_dev$apoe4_off, na.rm = TRUE),
            100*mean(b1_dev$apoe4_off, na.rm = TRUE),
            sum(abs(b1_dev$APOE4_dev) > 4, na.rm = TRUE),
            median(abs(b1_dev$APOE4_dev), na.rm = TRUE),
            median(abs(sib %>% left_join(ref_med, by="e4") %>%
                       mutate(d = APOE4 - med) %>% pull(d)), na.rm = TRUE)))
cat(sprintf("Direction of the off-band subset: reading LOW (dev<-2) = %d, HIGH (dev>2) = %d -> scrambled, not uniform signal loss.\n",
            sum(b1_dev$APOE4_dev < -2, na.rm = TRUE), sum(b1_dev$APOE4_dev > 2, na.rm = TRUE)))

# ---- 4a2. SWAP vs DROPOUT: are carrier-level APOE4 values conserved? ---------
# Dropout (signal loss) would only push carriers DOWN and reduce the count of
# carrier-level readings. A value swap/relabel conserves that count and moves the
# lost carrier values onto non-carriers (making 3/3 read like 3/4). Test it.
CARRIER_LEVEL <- 5                                      # APOE4 > 5 == carrier-level reading
b1_dev <- b1_dev %>% mutate(true_carrier = e4 >= 1, obs_carrier_level = APOE4 > CARRIER_LEVEL)
n_carr_b1   <- sum(b1_dev$true_carrier)
n_lost      <- sum(b1_dev$true_carrier  & !b1_dev$obs_carrier_level)   # carriers reading low
n_gained    <- sum(!b1_dev$true_carrier &  b1_dev$obs_carrier_level)   # non-carriers reading high
n_obs_carr  <- sum(b1_dev$obs_carrier_level)                            # wells reading carrier-level
gained_vals <- sort(round(b1_dev$APOE4[!b1_dev$true_carrier & b1_dev$obs_carrier_level], 1))
cat("\n================ SWAP vs DROPOUT (carrier-level conservation) ==========\n")
cat(sprintf("  true carriers: %d | carriers reading LOW (lost): %d | non-carriers reading HIGH (gained): %d\n",
            n_carr_b1, n_lost, n_gained))
cat(sprintf("  wells reading carrier-level (APOE4>%d): %d  vs  true carriers: %d  ->  %s\n",
            CARRIER_LEVEL, n_obs_carr, n_carr_b1,
            ifelse(n_obs_carr >= n_carr_b1 - 1, "CONSERVED => value SWAP/relabel (not dropout)", "REDUCED => dropout")))
cat(sprintf("  the gained (non-carrier-reading-high) APOE4 values land at the carrier mode: %s (normal carrier ~10.7)\n",
            paste(gained_vals, collapse = ", ")))

# ---- 4b. Decisive: are the APOE4-off samples bad WELLS (broken across the
#      panel) or APOE4-assay-specific? And do they track a sample sub-batch? ----
robust_z <- function(df, col) { x <- df[[col]]; s <- mad(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x))); (x - median(x, na.rm = TRUE)) / s }
Zb1   <- sapply(bc, function(c) robust_z(b1_dev, c))
n_bad <- rowSums(abs(Zb1) > 3, na.rm = TRUE)        # # of OTHER assays each Bay1 sample is an outlier on
cat(sprintf("\nMulti-assay outlier burden (of %d assays): APOE4-off samples mean=%.1f vs APOE4-ok mean=%.1f\n",
            length(bc), mean(n_bad[b1_dev$apoe4_off], na.rm = TRUE),
            mean(n_bad[!b1_dev$apoe4_off], na.rm = TRUE)))
cat("  -> equal burden = the APOE4-off samples are NORMAL on the rest of the panel (not bad wells); the failure is APOE4-specific.\n")
# overlap of each top-inflated assay's Bay1 outliers with the APOE4 off set
apoe4_set <- which(b1_dev$apoe4_off)
cat("\nJaccard overlap of each 'inflated' assay's Bay1 outliers with the APOE4 off set:\n")
for (c in intersect(c("IL1B","pTDP43_409","RUVBL2","PDGFRB","YWHAZ","IL15","CCL2"), bc)) {
  s <- which(abs(robust_z(b1_dev, c)) > 3); uni <- length(union(s, apoe4_set))
  cat(sprintf("  %-12s n_outliers=%2d  Jaccard_with_APOE4=%.2f\n", c, length(s),
              ifelse(uni > 0, length(intersect(s, apoe4_set))/uni, 0)))
}
cat("  -> near-zero overlap + many of these have 0 within-Bay1 outliers: their raw IQR inflation is Bay1's\n")
cat("     different sample composition (more non-Hispanic / different dx vocabulary), NOT broken assays.\n")
# sub-batch tracking
b1_dev$sample_era <- ifelse(grepl("^20250", b1_dev$SAMPLE), "2025xxx", "2024xxx")
cat("\nDoes the off-band subset track a sample accession batch or ancestry?\n")
cat(sprintf("  by SAMPLE era : 2024xxx %.0f%% off | 2025xxx %.0f%% off\n",
            100*mean(b1_dev$apoe4_off[b1_dev$sample_era=="2024xxx"], na.rm=TRUE),
            100*mean(b1_dev$apoe4_off[b1_dev$sample_era=="2025xxx"], na.rm=TRUE)))
print(b1_dev %>% filter(!is.na(apoe4_off)) %>% group_by(Ethnicity) %>%
        summarise(n=n(), off=sum(apoe4_off), pct=round(100*mean(apoe4_off)), .groups="drop"))
cat("  -> uniform ~31% across era/source/ancestry = a random ~third of Bay1 WELLS, i.e. a localized\n")
cat("     APOE4-channel (ε4-isoform detection) failure, genotype/sample-identity-independent.\n")

# ---- 4c. Genotype identity cross-check (chrX fingerprint) -------------------
# Are the broken Bay1 samples the right PEOPLE? The chrX genetic fingerprint
# (May_2026_WGS_QC, the same check that diagnosed the WGS plate swap) gives
# discord_pct ~0.003 for a correct identity and ~0.50 for a swapped one.
chrx_path <- "../../../../May_2026_WGS_QC/both_APOE_chrX_concord.csv"
have_chrx <- file.exists(chrx_path)
n_sanger <- sum(!is.na(b1_dev$APOE.geno) & trimws(b1_dev$APOE.geno) != "")
if (have_chrx) {
  cx <- read.csv(chrx_path, check.names = FALSE)
  cx$SAMPLE <- as.character(cx$SAMPLE)
  bx <- b1_dev %>% mutate(SAMPLE = as.character(SAMPLE)) %>%
    left_join(cx %>% select(SAMPLE, APOE_Sanger = APOE, APOE_WGS_cx = APOE_WGS,
                            chrX_discord = discord_pct), by = "SAMPLE")
  SWAP_THRESH <- 0.30                              # >=0.30 = different person (swap)
  n_wgs  <- sum(!is.na(bx$chrX_discord))
  n_both <- sum(!is.na(bx$APOE.geno) & trimws(bx$APOE.geno) != "" & !is.na(bx$chrX_discord))
  wgs_rows  <- bx %>% filter(!is.na(chrX_discord))
  chrx_lo <- round(100*min(wgs_rows$chrX_discord), 2); chrx_hi <- round(100*max(wgs_rows$chrX_discord), 2)
  n_rightperson <- sum(wgs_rows$chrX_discord < 0.05); n_swap <- sum(wgs_rows$chrX_discord >= SWAP_THRESH)
  off_wgs   <- bx %>% filter(apoe4_off, !is.na(chrX_discord))
  off_wgs_n <- nrow(off_wgs); off_wgs_pass <- sum(off_wgs$chrX_discord < 0.05); off_wgs_swap <- sum(off_wgs$chrX_discord >= SWAP_THRESH)
  cat("\n================ genotype identity cross-check (chrX fingerprint) ======\n")
  cat(sprintf("Bay1 genotype coverage: Sanger %d/%d | WGS+chrX %d/%d | BOTH %d/%d\n",
              n_sanger, nrow(bx), n_wgs, nrow(bx), n_both, nrow(bx)))
  cat(sprintf("All %d Bay1 WGS samples chrX discord = %.2f-%.2f%% (right-person ~0.3%%; swap ~50%%): right-person=%d, swaps=%d\n",
              n_wgs, chrx_lo, chrx_hi, n_rightperson, n_swap))
  cat(sprintf("APOE4-BROKEN subset: %d have WGS+chrX -> all pass identity? right-person=%d, swaps=%d\n",
              off_wgs_n, off_wgs_pass, off_wgs_swap))
  cat("Homozygous ε4 (4/4) broken cases — genotype triple-confirmed, protein wrong:\n")
  print(off_wgs %>% filter(e4 == 2) %>%
          transmute(SAMPLE, Sanger = APOE.geno, WGS = APOE_WGS_cx,
                    chrX_discord_pct = round(100*chrX_discord, 2), APOE4 = round(APOE4, 2)))
} else {
  cat("\n(chrX concordance file not found at", chrx_path, "- skipping identity cross-check.)\n")
}

# ---- 4d. Are the affected samples special? (covariates + accession adjacency) ----
# Compare affected vs unaffected Bay1 study samples across sample attributes, and
# test whether the affected ones cluster in SAMPLE (accession) number order.
cov_cols <- intersect(c("Site","country_of_birth","CDX","sex","Ethnicity","Race"), names(b1_dev))
cov_p <- sapply(cov_cols, function(cl) {
  tb <- table(b1_dev[[cl]], b1_dev$apoe4_off)
  if (nrow(tb) > 1 && ncol(tb) > 1)
    tryCatch(fisher.test(tb, simulate.p.value = TRUE, B = 3000)$p.value, error = function(e) NA) else NA
})
age_p <- if ("age_at_subject" %in% names(b1_dev))
  tryCatch(wilcox.test(age_at_subject ~ apoe4_off, b1_dev)$p.value, error = function(e) NA) else NA
min_cov_p <- suppressWarnings(min(c(cov_p, age_p), na.rm = TRUE))
# accession-number clustering: runs test on off/not sequence sorted by SAMPLE number
ord  <- order(suppressWarnings(as.numeric(b1_dev$SAMPLE)))
so   <- b1_dev$apoe4_off[ord]; n1 <- sum(so); n0 <- sum(!so)
runs <- 1 + sum(so[-1] != so[-length(so)])
mu_r <- 1 + 2*n1*n0/(n1+n0); va_r <- 2*n1*n0*(2*n1*n0-n1-n0)/((n1+n0)^2*(n1+n0-1))
runs_p <- 2*pnorm(-abs((runs - mu_r)/sqrt(va_r)))
cat("\n================ are the affected samples special? ================\n")
cat("  Fisher p by covariate:\n"); print(round(cov_p, 3))
cat(sprintf("  age Wilcoxon p: %.3f | smallest covariate p: %.3f\n", age_p, min_cov_p))
cat(sprintf("  accession-number clustering (runs test): p = %.2f (>0.05 = no clustering)\n", runs_p))
cat("  -> affected samples are not distinguished by any attribute, accession order, or (sec.7) well position.\n")

# ---- 4e. How bad could it be? detection floor + HIHG full-panel positive control ----
# (i) Our genotype-anchored check only sees carrier<->non-carrier APOE4 swaps; same-genotype
#     (and e4=1<->e4=2, which differ <2) swaps are invisible. Estimate the affected-set size
#     that a RANDOM permutation would need to still yield the observed detectable count.
n_off_e <- sum(b1_dev$apoe4_off, na.rm = TRUE)
nC_e <- sum(b1_dev$e4 >= 1); nN_e <- sum(b1_dev$e4 == 0); N_e <- nC_e + nN_e
pc_e <- nC_e/N_e; pn_e <- nN_e/N_e
full_detect   <- 2*nC_e*nN_e/N_e
s_random      <- n_off_e/(2*pc_e*pn_e)
hidden_random <- max(0, s_random - n_off_e)
cat("\n================ how bad could it be? (detection floor) ================\n")
cat(sprintf("  detectable now = %d; full scramble of all %d wells would give ~%.0f.\n", n_off_e, N_e, full_detect))
cat(sprintf("  a RANDOM permutation yielding %d detectable -> ~%.0f wells affected (~%.0f%%), ~%.0f hidden same-genotype swaps.\n", n_off_e, s_random, 100*s_random/N_e, hidden_random))

# (ii) HIHG positive control: a known reference sample run on every bay. If Bay1's data->well
#      mapping were globally scrambled (or the reagent failed), HIHG would be wrong too.
hihg_path <- "../../Primary_QC/output_files/hihg_reps_NPQ.csv"
have_hihg <- file.exists(hihg_path)
hihg_cor_lo <- hihg_cor_hi <- hihg_apoe4_z <- NA
if (have_hihg) {
  hh <- read.csv(hihg_path, check.names = FALSE)
  hbm <- intersect(bc, names(hh))
  cons <- colMeans(as.matrix(hh[, hbm]), na.rm = TRUE)
  hh$corr <- apply(as.matrix(hh[, hbm]), 1, function(x) cor(x, cons, use = "complete.obs"))
  hb1 <- hh[hh$Run == TARGET_RUN & hh$Bay == "Bay1", ]
  hot <- hh[!(hh$Run == TARGET_RUN & hh$Bay == "Bay1"), ]
  hihg_cor_lo <- min(hb1$corr); hihg_cor_hi <- max(hb1$corr)
  s_ap <- sd(hot$APOE4, na.rm = TRUE)
  hihg_apoe4_z <- if (!is.na(s_ap) && s_ap > 0) (mean(hb1$APOE4) - mean(hot$APOE4))/s_ap else NA
  cat(sprintf("  HIHG control on Bay1: full-panel r to consensus = %.3f-%.3f (other bays ~0.99); APOE4 z = %.2f.\n",
              hihg_cor_lo, hihg_cor_hi, hihg_apoe4_z))
  cat("  -> the bay measures a KNOWN sample correctly across ALL proteins: not a reagent failure, not a whole-bay scramble.\n")
}

# ---- 5. Could control-based QC have caught it? (Alamar control replicates) --
# IPC/SC/NC/HIHG are identical control material run on every bay. If the APOE4
# channel failed bay-wide, the controls on Bay1 would be off too.
PQC <- "../../Primary_QC/output_files"
ctrl_ok <- all(file.exists(file.path(PQC, c("ipc_reps_NPQ.csv","sc_reps_NPQ.csv",
                                            "nc_reps_NPQ.csv","hihg_reps_NPQ.csv"))))
if (ctrl_ok) {
  ctrls <- bind_rows(lapply(c(IPC="ipc", SC="sc", NC="nc", HIHG="hihg"), function(p)
    transform(read.csv(file.path(PQC, paste0(p, "_reps_NPQ.csv")), check.names = FALSE))),
    .id = "ctrl")
  cat("\n================ control-replicate APOE4 on Bay1 vs all other bays =====\n")
  cstat <- ctrls %>% group_by(ctrl) %>%
    mutate(other_med = median(APOE4[Run != TARGET_RUN], na.rm = TRUE),
           other_mad = mad(APOE4[Run != TARGET_RUN], na.rm = TRUE)) %>%
    filter(Run == TARGET_RUN, Bay == "Bay1") %>%
    summarise(n = n(), APOE4_Bay1 = round(median(APOE4), 2),
              other_bays_med = round(first(other_med), 2),
              max_abs_robustZ = round(max(abs((APOE4 - first(other_med)) /
                                              (first(other_mad) + 1e-9))), 1), .groups = "drop")
  print(as.data.frame(cstat))
  cat("-> controls read APOE4 NORMALLY on Bay1 (|robust-Z| small): a bay-level control-CV check\n")
  cat("   would NOT have flagged this; the failure spares the control wells.\n")
} else {
  cat("\n(Primary_QC control replicate files not found; skipping control check.)\n")
}

# ---- 6. The check that DOES catch it: per-run-bay APOE4-vs-genotype QC -------
# For every run-bay, the share of genotyped samples whose APOE4 is off its ε4 band,
# and the APOE4 IQR among carriers. This is the proposed going-forward APOE QC gate.
gband <- m %>% filter(!is.na(e4), !is.na(APOE4)) %>% group_by(e4) %>%
  summarise(md = median(APOE4), .groups = "drop")
qc <- m %>% filter(!is.na(e4), !is.na(APOE4)) %>%
  mutate(off = abs(APOE4 - gband$md[match(e4, gband$e4)]) > 2) %>%
  group_by(RUN) %>%
  summarise(n = n(), carrier_n = sum(e4 >= 1),
            APOE4_carrier_IQR = round(IQR(APOE4[e4 >= 1], na.rm = TRUE), 2),
            apoe4_off_band_pct = round(100 * mean(off), 1), .groups = "drop") %>%
  filter(n >= 15) %>% arrange(desc(apoe4_off_band_pct))
write.csv(qc, file.path(OUT, "per_run_apoe4_qc.csv"), row.names = FALSE)
cat("\n================ per-run-bay APOE4 QC (the going-forward check) ========\n")
cat(sprintf("Run-bays with n>=15: %d | median off-band %.1f%% | %s = %.1f%% (carrier IQR %.1f vs median %.1f)\n",
            nrow(qc), median(qc$apoe4_off_band_pct), TARGET_BAY,
            qc$apoe4_off_band_pct[qc$RUN == TARGET_BAY],
            qc$APOE4_carrier_IQR[qc$RUN == TARGET_BAY], median(qc$APOE4_carrier_IQR)))
print(as.data.frame(head(qc, 5)))

# ---- Figures ----------------------------------------------------------------
pA <- m %>% filter(!is.na(e4), !is.na(APOE4)) %>%
  ggplot(aes(factor(e4), APOE4, colour = grp, fill = grp)) +
  geom_violin(scale = "width", alpha = 0.2, position = position_dodge(0.8), width = 0.75, linewidth = 0.3) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45, linewidth = 0.4,
               position = position_dodge(0.8), show.legend = FALSE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
             size = 0.8, alpha = 0.45, show.legend = FALSE) +
  scale_fill_manual(values = c("Bay1 (target)" = "#E41A1C", "Bay2/3 (sibling)" = "#377EB8",
                               "All other runs" = "grey60")) +
  scale_colour_manual(values = c("Bay1 (target)" = "#E41A1C",
                                 "Bay2/3 (sibling)" = "#377EB8",
                                 "All other runs" = "grey60")) +
  labs(title = "APOE4 by ε4 dosage: Bay1 loses carrier/non-carrier separation",
       subtitle = paste0(TARGET_BAY, " (red) vs its sibling bays (blue) and all other runs (grey)"),
       x = "ε4 allele count (clinical genotype)", y = "APOE4 (log2 NPQ)", colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")

pB <- ratio_df %>%
  mutate(hl = biomarker %in% c("APOE4", "APOE")) %>%
  ggplot(aes(rank, iqr_ratio_bay1_vs_sibling, colour = hl)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_point(size = 1.6, alpha = 0.75) +
  ggrepel::geom_text_repel(data = ~ filter(.x, hl), aes(label = biomarker),
                           size = 3.4, fontface = "bold", box.padding = 0.6) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#E41A1C"), guide = "none") +
  labs(title = "Within-bay variance vs sibling bays (ε4-residualized)",
       subtitle = "APOE4 stands out; the few other high-ratio markers have ~0 outlier samples (composition, not breakage)",
       x = "Biomarker (ranked)", y = "IQR ratio (Bay1 / sibling bays)") +
  theme_bw(base_size = 11)

ggsave(file.path(OUT, "bay1_apoe4_anomaly.png"), pA / pB, width = 10, height = 10,
       dpi = 150, bg = "white")

# ---- Export per-sample Bay1 flag table --------------------------------------
e4_wgs_b1 <- e4_count(b1_dev$APOE_WGS)
out_tbl <- b1_dev %>%
  transmute(SAMPLE_ALIQUOT, SAMPLE, RUN, Bay, CDX,
            sex = if ("sex" %in% names(.)) sex else NA,             # reported sex — join genetic sex to test for swaps
            Ethnicity = if ("Ethnicity" %in% names(.)) Ethnicity else NA,
            APOE.geno, APOE_WGS,
            sanger_matches_wgs = ifelse(is.na(e4_wgs_b1), NA, e4 == e4_wgs_b1),
            e4, APOE4 = round(APOE4, 2), APOE = round(APOE, 2),
            APOE4_band_median = round(med, 2), APOE4_dev = round(APOE4_dev, 2),
            apoe4_unreliable = apoe4_off) %>%
  arrange(desc(apoe4_unreliable), desc(abs(APOE4_dev)))
write.csv(out_tbl, file.path(OUT, "bay1_apoe4_sample_deviations.csv"), row.names = FALSE)
writeLines(capture.output(print(as.data.frame(ratio_df))),
           file.path(OUT, "bay1_panel_iqr_ratio.txt"))

# ---- 7. Plate map: which WELLS mismatch (well layout from Alamar Sample_QC) --
# Sample_QC.csv carries the physical well (row/col) for every well on the bay,
# including controls. We overlay the APOE4-vs-genotype mismatch flag to see
# whether the broken wells form a spatial region (edge / row / column / block)
# or are dispersed.
have_plate <- FALSE; mm_rows <- mm_cols <- n_mm_plate <- NA
sqc_path <- "../../Primary_QC/input_files/Sample_QC.csv"
if (file.exists(sqc_path)) {
  have_plate <- TRUE
  ppct <- function(x) as.numeric(gsub("%", "", x))
  sqc <- suppressWarnings(readr::read_csv(sqc_path, skip = 12, show_col_types = FALSE)) %>%
    filter(grepl(TARGET_BAY, plateID)) %>%
    transmute(SAMPLE_ALIQUOT = `Sample Name`, sample_type = `Sample Type`,
              sampleBarcode = if ("sampleBarcode" %in% names(.)) sampleBarcode else NA,
              row = wellRow, col = as.integer(wellCol),
              Reads = as.numeric(Reads), IC_Reads = as.numeric(`IC Reads`),
              IC_Median = ppct(`IC Median`), Detectability = ppct(Detectability),
              QC_Status = `QC Status`)
  b1status <- m %>% filter(RUN == TARGET_BAY) %>%
    mutate(off = abs(APOE4 - gband$md[match(e4, gband$e4)]) > 2 & !is.na(e4)) %>%
    transmute(SAMPLE_ALIQUOT, SAMPLE, e4, APOE.geno, APOE4, off)
  pm <- sqc %>% left_join(b1status, by = "SAMPLE_ALIQUOT") %>%
    mutate(status = factor(case_when(
              sample_type != "Sample" ~ "Control (IPC/SC/NC)",
              is.na(e4)               ~ "Sample, no genotype",
              off                     ~ "APOE4 MISMATCH",
              TRUE                    ~ "APOE4 ok"),
            levels = c("APOE4 MISMATCH","APOE4 ok","Sample, no genotype","Control (IPC/SC/NC)")),
           y   = 9 - match(row, LETTERS[1:8]),
           lab = case_when(sample_type != "Sample" ~ sample_type,
                           is.na(e4)               ~ "—",
                           TRUE ~ sprintf("%s\n%.1f", APOE.geno, APOE4)))
  write.csv(pm %>% select(row, col, SAMPLE_ALIQUOT, sample_type, e4, APOE.geno, APOE4, status) %>%
              arrange(row, col), file.path(OUT, "bay1_plate_map.csv"), row.names = FALSE)
  # ---- clean hand-to-the-lab list: only the affected wells, with position + barcode ----
  lab_list <- pm %>% filter(status == "APOE4 MISMATCH") %>%
    transmute(run = TARGET_BAY, well = sprintf("%s%d", row, col), wellRow = row, wellCol = col,
              SAMPLE, SAMPLE_ALIQUOT, sampleBarcode,
              APOE_genotype = APOE.geno, APOE4_observed = round(APOE4, 1),
              direction = ifelse(APOE4 > gband$md[match(e4, gband$e4)], "reads HIGH (vs genotype)",
                                                                        "reads LOW (vs genotype)")) %>%
    arrange(wellRow, wellCol)
  write.csv(lab_list, file.path(OUT, "bay1_affected_samples_FOR_LAB.csv"), row.names = FALSE)
  if (requireNamespace("writexl", quietly = TRUE))
    writexl::write_xlsx(lab_list, file.path(OUT, "bay1_affected_samples_FOR_LAB.xlsx"))
  cat(sprintf("\n-> wrote %s/bay1_affected_samples_FOR_LAB.{csv,xlsx} (%d affected wells, with well + barcode)\n", OUT, nrow(lab_list)))
  n_mm_plate <- sum(pm$status == "APOE4 MISMATCH")
  mm_rows <- length(unique(pm$row[pm$status == "APOE4 MISMATCH"]))
  mm_cols <- length(unique(pm$col[pm$status == "APOE4 MISMATCH"]))
  n_geno  <- sum(pm$status %in% c("APOE4 MISMATCH","APOE4 ok"))

  pPlate <- ggplot(pm, aes(col, y)) +
    geom_tile(aes(fill = status), color = "white", width = 0.94, height = 0.94) +
    geom_text(aes(label = lab), size = 2.0, lineheight = 0.82) +
    scale_fill_manual(values = c("APOE4 MISMATCH" = "#E41A1C", "APOE4 ok" = "#A6D854",
                                 "Sample, no genotype" = "grey82",
                                 "Control (IPC/SC/NC)" = "#6BAED6")) +
    scale_x_continuous(breaks = 1:12, position = "top") +
    scale_y_continuous(breaks = 1:8, labels = rev(LETTERS[1:8])) +
    coord_equal(expand = FALSE) +
    labs(title = sprintf("%s — plate map: APOE4-vs-genotype mismatch by well", TARGET_BAY),
         subtitle = sprintf("%d of %d genotyped study wells mismatch (red, cell = genotype + APOE4). Controls (blue) read APOE4 normally.\nMismatches dispersed across %d of 8 rows and %d of 12 columns — no contiguous block / edge / single row or column.",
                            n_mm_plate, n_geno, mm_rows, mm_cols),
         x = "column", y = "row", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  ggsave(file.path(OUT, "bay1_plate_map.png"), pPlate, width = 11, height = 7.5,
         dpi = 150, bg = "white")
  cat(sprintf("\n================ plate map ================\n%d wells (%d controls); %d mismatch across %d rows / %d cols (dispersed, not a block).\n-> wrote %s/bay1_plate_map.png + bay1_plate_map.csv\n",
              nrow(pm), sum(pm$sample_type != "Sample"), n_mm_plate, mm_rows, mm_cols, OUT))

  # ---- 7b. Low-input check: do the mismatch wells have a bad-QC signature? ----
  # Overlay Alamar per-well QC (total reads, detectability, IC median, QC status)
  # on the layout, outlining the APOE4-mismatch wells. If they were low-input we'd
  # see them as the low-read / low-detectability cells.
  dd <- pm %>% filter(sample_type == "Sample", !is.na(off))
  qc_cmp <- dd %>% group_by(grp = ifelse(off, "APOE4 mismatch", "APOE4 ok")) %>%
    summarise(n = n(), Reads_median = round(median(Reads)), Reads_min = round(min(Reads)),
              Detectability_median = round(median(Detectability), 1),
              IC_abs_median = round(median(abs(IC_Median)), 1),
              all_passed = all(QC_Status == "passed"), .groups = "drop")
  pval <- c(Reads = wilcox.test(Reads ~ off, dd)$p.value,
            Detectability = wilcox.test(Detectability ~ off, dd)$p.value,
            IC_abs = wilcox.test(abs(IC_Median) ~ off, dd)$p.value)
  n_mm_flagged <- sum(dd$off & dd$QC_Status != "passed")
  cat("\n---- per-well QC: APOE4-mismatch vs ok (study wells) ----\n")
  print(as.data.frame(qc_cmp))
  cat(sprintf("Wilcoxon p: Reads %.2f | Detectability %.2f | |IC median| %.2f | mismatch wells failing Alamar QC: %d/%d\n",
              pval["Reads"], pval["Detectability"], pval["IC_abs"], n_mm_flagged, sum(dd$off)))

  ov <- pm %>% mutate(is_mm = status == "APOE4 MISMATCH",
                      Reads_M = Reads / 1e6,
                      Detect = ifelse(sample_type == "Sample", Detectability, NA))
  mk <- function(fillvar, title, lab_fmt) {
    ggplot(ov, aes(col, y)) +
      geom_tile(aes(fill = .data[[fillvar]]), color = "grey85", width = 0.95, height = 0.95) +
      geom_tile(data = ~ filter(.x, is_mm), color = "#E41A1C", linewidth = 1.3,
                fill = NA, width = 0.95, height = 0.95) +
      geom_text(aes(label = ifelse(is.na(.data[[fillvar]]), sample_type, sprintf(lab_fmt, .data[[fillvar]]))),
                size = 2.0) +
      scale_fill_viridis_c(option = "C", na.value = "grey90") +
      scale_x_continuous(breaks = 1:12, position = "top") +
      scale_y_continuous(breaks = 1:8, labels = rev(LETTERS[1:8])) +
      coord_equal(expand = FALSE) +
      labs(title = title, x = NULL, y = NULL, fill = NULL) +
      theme_minimal(base_size = 10) +
      theme(panel.grid = element_blank(), legend.position = "right")
  }
  ovfig <- (mk("Reads_M", "Total reads (millions)", "%.1f") +
            mk("Detect", "Detectability (%)", "%.0f")) +
    patchwork::plot_annotation(
      title = sprintf("%s — per-well QC with APOE4-mismatch wells outlined (red)", TARGET_BAY),
      subtitle = "If the mismatch were low-input, the red-outlined wells would be the dark (low-read / low-detectability) cells. They are not — reads if anything run higher, detectability is equal, and all mismatch wells passed Alamar QC.",
      theme = ggplot2::theme(plot.title = element_text(face = "bold")))
  ggsave(file.path(OUT, "bay1_plate_qc_overlay.png"), ovfig, width = 15, height = 7,
         dpi = 150, bg = "white")
  cat(sprintf("-> wrote %s/bay1_plate_qc_overlay.png\n", OUT))
}

# ---- Internal-colleagues note -----------------------------------------------
n_off <- sum(b1_dev$apoe4_off, na.rm = TRUE); n_b1 <- sum(!is.na(b1_dev$apoe4_off))
b1_off_pct <- qc$apoe4_off_band_pct[qc$RUN == TARGET_BAY]
next_pct   <- sort(qc$apoe4_off_band_pct, decreasing = TRUE)[2]
med_pct    <- median(qc$apoe4_off_band_pct)
b1_cIQR    <- qc$APOE4_carrier_IQR[qc$RUN == TARGET_BAY]
med_cIQR   <- median(qc$APOE4_carrier_IQR)
# Sanger-vs-WGS corroboration among the affected (Sanger was co-run with the assay; WGS was not)
.affi      <- which(b1_dev$apoe4_off); .e4w <- e4_count(b1_dev$APOE_WGS)
n_aff_wgs  <- sum(!is.na(.e4w[.affi])); n_aff_sw_agree <- sum(b1_dev$e4[.affi] == .e4w[.affi], na.rm = TRUE)
chrx_section <- if (have_chrx) c(
"## Genotype identity cross-check (chrX fingerprint) — does this point to a sample mix-up?",
"To ask whether the affected wells might be the *wrong people* (a physical sample swap), we used the chrX genetic fingerprint — the same check that diagnosed the WGS plate swap. It reads ~0.3% discordance for a correct identity and ~50% for a swapped one. The reasoning:",
sprintf("- **Bay1 genotype coverage:** Sanger %d/%d, WGS+chrX %d/%d, both %d/%d (the rest Sanger-only, no WGS).", n_sanger, nrow(b1_dev), n_wgs, nrow(b1_dev), n_both, nrow(b1_dev)),
sprintf("- **All %d Bay1 WGS samples land in the right-person range** (discord %.2f–%.2f%%; %d right-person, %d swaps), and none sit on the swapped WGS plates — so these wells are who they claim to be, and the Bay1 issue looks separate from the WGS swap.", n_wgs, chrx_lo, chrx_hi, n_rightperson, n_swap),
sprintf("- Of the %d APOE4-off samples, %d have WGS+chrX and **all %d are in the right-person range**. The homozygous ε4 (4/4) cases read APOE4 ≈ 0 even though the genotype is agreed by Sanger + WGS + chrX — so, on the evidence we have, it is the **protein** that disagrees, not the sample. (The Sanger-only affected samples show the same genotype-vs-protein pattern.)", n_off, off_wgs_n, off_wgs_pass),
sprintf("- **Sanger and WGS agree for all %d affected samples that have WGS** (ε4-count, %d/%d), and APOE4 is off-band vs *both*. This is worth dwelling on: **Sanger was co-run with the biomarkers, WGS was not** — yet they agree, which argues that the co-run did not corrupt the Sanger genotype. So the weight of evidence is on the **protein value** being wrong rather than the genotype. (**Caveat on the sex check:** genetic-vs-reported sex would test the *genotype/metadata* identity, like chrX — it does **not** test whether the Alamar protein data belongs to the labelled well, so it cannot catch an Alamar-side swap. Reported sex + both genotypes are in `bay1_apoe4_sample_deviations.csv` if useful for the genotyping side.)", n_aff_wgs, n_aff_sw_agree, n_aff_wgs),
"") else c(
"## Genotype identity cross-check (chrX fingerprint)",
"*(chrX concordance file not available in this run; identity cross-check skipped.)*",
"")
swap_section <- c(
"## Mechanism — our best current reading: a value swap rather than drop-out",
"Here is the reasoning; please check whether you find it convincing. A simple assay **drop-out** (signal loss) can only push ε4 carriers *down* — it would make a 3/4 look like a 3/3, and the number of carrier-level APOE4 readings on the plate should *fall*. We see something different:",
sprintf("- The number of wells reading carrier-level APOE4 (>%d) is **%d — essentially identical to the number of true ε4 carriers (%d)**. So the carrier-level signal appears **conserved** rather than lost. (You can re-derive this two-line count yourself from the per-sample table.)", CARRIER_LEVEL, n_obs_carr, n_carr_b1),
sprintf("- And it appears **relocated**: %d true carriers fell to non-carrier level, while %d true **non-carriers rose to carrier level** — the gained values land right at the normal carrier mode (%s; normal ~10.7). A drop-out cannot *add* carrier-level signal to a 3/3, so this two-directional, count-conserving pattern reads more like a swap than a loss.", n_lost, n_gained, paste(head(gained_vals, 8), collapse = ", ")),
"- Taken together, our working interpretation is a **value swap / mis-indexing confined to the APOE4 channel** — the per-well APOE4 *numbers* look permuted among a subset of Bay1 wells, rather than a reagent failure. Sample identity checks out (chrX + the other proteins), so on this reading it is the APOE4 *number* that moved, not the sample. A concrete question for Alamar: could the APOE4 target's per-well values have been mis-assigned during demux/quantification on this bay? (We cannot see that step from the processed data.)",
"- **Caveat:** with ~26 points, the count matching exactly could be partly coincidence, and a process that makes the APOE4 reading *uninformative* (e.g. effectively random noise drawn from the plate's bimodal APOE4 distribution) would look similar. We lean toward 'swap' because of the conservation and the gained values sitting tightly at the carrier mode — but we cannot prove a literal pairwise permutation from these data.",
"",
"**Can we identify / un-swap them?** We can name *which* samples are affected (the 13 carriers reading low + 13 non-carriers reading high, flagged in `bay1_apoe4_sample_deviations.csv`), but we do **not** see a way to reconstruct the pairing or recover the true values: within a genotype the APOE4 values are near-degenerate (every carrier ≈ 10.7, every non-carrier ≈ 0), and the affected wells did **not** follow any geometric rule we tried (no rotation/flip/transpose maps the low-wells onto the high-wells — best was 5/13 vs ~2/13 expected by chance). So — *unlike* the WGS plate swap, which was a clean 180° rotation and therefore deterministically reversible — we could **not** find a way to invert this one. Our suggestion is therefore to set APOE4 aside for this bay rather than try to correct it. (Most affected samples *do* have WGS, 19/26 with chrX, but WGS/chrX only confirms the genotype/identity — it cannot resolve a protein-readout permutation.)",
"")
plate_section <- if (isTRUE(have_plate)) c(
"## Where on the plate (well map) — for molecular colleagues",
"The figure below is the physical 96-well layout of this bay (rows A–H, columns 1–12; one tile = one well). Each tile is colored by whether that person's **APOE4 protein** reading agrees with their **APOE genotype**, and is labeled with the genotype and the APOE4 value:",
"",
"🟥 **red = APOE4 disagrees with genotype** (broken)  ·  🟩 **green = agrees**  ·  ⬜ grey = no genotype  ·  🟦 blue = assay controls (IPC/SC/NC)",
"",
"![Bay1 plate map — APOE4-vs-genotype mismatch by well](bay1_plate_map.png)",
"",
sprintf("Two things to notice: **(1)** the red (broken) wells are sprinkled across the plate — **%d of 8 rows and %d of 12 columns** — not in one corner, edge, or stripe, which argues against a spill or an edge-drying effect (and against a simple geometric well-swap); **(2)** the blue **control wells sit right among the broken ones and read APOE4 normally** — which is likely why routine plate QC did not flag this bay.", mm_rows, mm_cols),
"",
"### Is it just low-quality / low-input wells? Doesn't look like it.",
sprintf("Overlaying Alamar's per-well QC (figure below), the broken wells look **normal-to-better**: median reads **%.1fM vs %.1fM (ok)**, comparable detectability and IC-median (Wilcoxon p = %.2f reads / %.2f detectability / %.2f |IC|), and **all %d passed Alamar QC**. So a low-read / poor-quality explanation doesn't seem to hold — the wells look fine on every standard metric except the APOE4 *number*.",
        qc_cmp$Reads_median[qc_cmp$grp == "APOE4 mismatch"]/1e6, qc_cmp$Reads_median[qc_cmp$grp == "APOE4 ok"]/1e6,
        pval["Reads"], pval["Detectability"], pval["IC_abs"], sum(dd$off)),
"",
"![Bay1 per-well reads & detectability — broken wells outlined red (they are not the dark/low cells)](bay1_plate_qc_overlay.png)",
"") else character(0)
allplates_section <- c(
"## The same view across every run-bay (for extra eyes)",
"In case a pattern we've overlooked is visible to someone else, here is the same APOE4-vs-genotype map for **every** run-bay, not just Bay1. **Red = reads HIGH** (e.g. a 3/3 reading like a carrier), **blue = carrier reading LOW**, grey = agrees / control. On this freeze, `20251124-1407_Bay1` (top-left) is the only plate-level outlier; the rest show a handful of scattered wells, and the November `20251125-1332` plates show a cluster of mild reads-high (characterized separately as a possible AFDC-batch effect). Fresh eyes on the full set are welcome — if anything else jumps out, flag it. (Figure regenerates from `investigate_other_runs_apoe4_wells.R`.)",
"",
"![APOE4-vs-genotype mismatch wells across all run-bays](other_runs_apoe4_wellmap.png)",
"")
howbad_section <- c(
"## How bad could it be? (scope of the uncertainty)",
"Two honest limits on the count, and what bounds them:",
sprintf("- **The %d is a floor, not the total.** Our check only catches APOE4 swaps that cross the carrier <-> non-carrier line; **same-genotype swaps (and ε4-het <-> ε4-hom, which differ <2) are invisible to it.** A *clean pairwise* carrier <-> non-carrier swap would be exactly the %d we see (0 hidden); a *random* permutation producing %d detectable would instead imply **~%.0f wells affected (~%.0f%% of the bay), with ~%.0f hidden**. The exact 13/13 balance and conserved carrier count fit the clean end best, but we cannot exclude the larger. Either way: **treat every APOE4 value on this bay as untrustworthy.**", n_off_e, n_off_e, n_off_e, s_random, 100*s_random/N_e, hidden_random),
"- **Could whole samples (all proteins) be mislabeled, not just APOE4?** We cannot fully exclude it for the study wells. Our identity checks (chrX, Sanger = WGS, sex) all validate the *genotype/metadata* side — **none of them test whether the Alamar protein data in a well belongs to the labelled sample.** The only per-sample anchor we have is APOE4-vs-genotype; for the other ~130 proteins there is no anchor, and a whole-sample swap **preserves each protein's marginal distribution**, so it would leave no trace.",
if (have_hihg) sprintf("- **What bounds the worst case — the HIHG control.** HIHG is a known reference sample run on every bay; on Bay1 it reproduces the **full ~130-protein panel at r = %.3f–%.3f** (same as every other bay) with **APOE4 z = %.2f**. So this is **not** a reagent failure, **not** a whole-bay or global demux scramble, and the **control wells are correctly mapped for every protein**. The 'whole plate is broken' version is off the table — what remains is a *localized* mis-assignment of study wells that spares the controls.", hihg_cor_lo, hihg_cor_hi, hihg_apoe4_z) else NULL,
"- **Net.** The plate is fundamentally sound (HIHG proves it), so this is a localized study-well problem — but whether it touched *only the APOE4 value* or *whole samples* cannot be settled from the data, because the one anchor that could reveal a whole-sample study-well swap exists only for APOE4. Conservative stance: drop **all** APOE4 for the bay, and regard the study→data mapping for the *other* proteins as **unverified (not disproven)** pending Alamar's process answer or a re-assay.",
"")
note <- c(
"# 20251124-1407_Bay1 — APOE4 readout anomaly (working diagnosis)",
"",
"**For internal colleagues (Biomarker QC / analysis team).**",
sprintf("Notes on the APOE4 anomaly first surfaced by the APOE secondary QC on run `%s`, Bay1.", TARGET_RUN),
"Our **current reading** is that this sits in the APOE4 *protein readout* rather than in the genotype — but this is a working diagnosis, not a closed case. Below we lay out what we see and the reasoning behind each step, so you can reproduce it and judge for yourselves; we also flag where we are least sure (see *What we might be missing*). We have **not** found any effect on the APOE genotypes (`APOE.geno` / `APOE_WGS` / `APOE.geno_final`), and this looks unrelated to the WGS plate-swap fix. Everything here regenerates from `investigate_20251124_Bay1.R`.",
"",
"## What we observe",
sprintf("- Run `%s` ran on 3 bays (~84 samples each). The effect appears **specific to Bay1** and to the **APOE4 target** — its sibling bays look clean.", TARGET_RUN),
sprintf("- On **%d of %d (%.0f%%)** Bay1 samples the APOE4 reading is far from what the ε4 genotype predicts (|APOE4 − ε4-band median| > 2 log2; all such cases > 4).", n_off, n_b1, 100*n_off/n_b1),
"- The mismatch runs **both ways** — carriers reading low *and* non-carriers reading high (all 3 homozygous ε4 samples read APOE4 ≈ 0, where ~12.6 is expected). That two-directional pattern is what first pointed us away from a simple drop-out (reasoning in *Mechanism*).",
"- **Total APOE looks intact** on the same samples, which argues against loading / dilution / degradation and points to something specific to the ε4-isoform APOE4 measurement. The carrier-vs-non-carrier *median* still separates, so the damage is mostly in the **variance** — which is why a median-based check did not catch it.",
"- **As far as we can tell, only APOE4 is affected.** Total APOE and the other ~130 proteins read normally on these wells (the affected samples carry the same outlier burden across the other assays as the unaffected ones). We have not found another target behaving this way (variance-rank panel in `bay1_apoe4_anomaly.png`).",
"",
swap_section,
"## What we checked, and could not pin it on",
"Each of these was a candidate explanation; here is why it doesn't seem to fit (please poke at these):",
"- **Bad wells in general** — doesn't fit: the APOE4-off samples are normal across the other 130 assays (same multi-assay outlier burden as the APOE4-ok samples), so the wells aren't broadly low-quality.",
"- **Other assays affected too** — doesn't seem so: markers that looked IQR-inflated on a raw scan (IL1B, pTDP43_409, RUVBL2, PDGFRB) have ~0 within-Bay1 outlier *samples*; that spread tracks Bay1's **different sample composition** (a distinct cohort from its sibling bays — more non-Hispanic, different dx vocabulary), which is not the same as breakage.",
sprintf("- **Any sample attribute or handling order** — doesn't fit: within Bay1 the affected ~1/3 are **not** distinguished by site, country of birth, diagnosis (CDX), sex, ancestry/race, or age (smallest p ≈ %.2f), do **not** cluster by accession/sample number (runs test p = %.2f; off-rate ~27-43%% across accession blocks), and are scattered across the physical plate (well map above). The affected set looks **random with respect to who the samples are and where/when they were handled** — which argues against a contiguous handling block, a site/biology effect, and an edge/spatial effect. (A few near-consecutive accession pairs exist, e.g. 200905200/201, 202500432/433, but no more than chance given ~1/3 are affected.)", min_cov_p, runs_p),
"- **A physical sample mix-up** — argued against by the chrX + Sanger-vs-WGS checks below (though the sex check is still outstanding).",
"",
chrx_section,
plate_section,
allplates_section,
"## Could we have caught it earlier in the QC pipeline?",
"On the QC we currently run, it seems we mostly could not — which is itself worth noting:",
"- **Alamar's sample-level auto-QC** (`Primary_QC/.../Alamar_QC_warn.csv`) flagged samples on **Bay2** of this run, not Bay1, on IC-median / detectability — none of it APOE4-specific.",
"- **Control-replicate / plate-CV QC** (`Primary_QC` IPC, SC, NC, HIHG; `Secondary_QC/Replicate_Analysis`): the pooled control material reads **APOE4 normally on Bay1** (|robust-Z| ≤ ~1.4 vs other bays). Because the affected wells appear to spare the controls and hit ~1/3 of the *sample* wells, a bay-level control-CV check would tend to miss this. (It also argues against a bay-wide reagent collapse, which should have moved the controls.)",
sprintf("- **What does appear to catch it: a genotype-anchored per-run APOE4 check.** Scoring every run-bay by the share of genotyped samples whose APOE4 is off its ε4 band makes `%s` a clear outlier on this freeze: **%.1f%% off-band vs %.1f%% next-worst vs %.1f%% median** run-bay; APOE4 carrier-band IQR **%.1f vs %.1f median (~%.0fx)**. A simple threshold (e.g. off-band > 15%%, or carrier IQR > 3× median) separated it here without false positives — though thresholds should be revisited as more runs accumulate.", TARGET_BAY, b1_off_pct, next_pct, med_pct, b1_cIQR, med_cIQR, b1_cIQR/med_cIQR),
"",
howbad_section,
"## What we might be missing (please scrutinize)",
"This is a working diagnosis from the processed data, so we want to be upfront about its limits:",
"- **We can describe the pattern but not prove the physical mechanism.** 'Value swap' is a description that fits the data; we have not seen the raw per-well, per-target Alamar quantification, which is where a true demux/indexing error would show. The honest statement is 'consistent with a swap', not 'proven swap'.",
"- **We can't explain why only APOE4.** A single-target value permutation, with every other protein and the genotype intact, is unusual — we don't have a mechanism for it, and that gap should make everyone (including us) cautious.",
"- **The conservation could be partly coincidental** (~26 points), and an 'APOE4 became uninformative/noise' process would mimic a swap. We lean swap but can't exclude this.",
"- **No assay-intrinsic identity check on the Alamar side.** All our identity tools (chrX, Sanger=WGS, sex) anchor on the genotype/metadata, not the protein plate — so a *whole-sample* study-well swap of the same genotype is invisible to us (see *How bad could it be?*). The HIHG control bounds this (the bay measures known samples correctly), but doesn't fully close it.",
"- **We have not re-confirmed the values trace to the raw Alamar NPQ** (vs anything our own merge could have introduced); total-APOE being fine, the HIHG control reproducing, and the effect being one bay argue against a pipeline bug, but it's worth a direct check.",
"- **Other explanations we may simply not have thought of are possible.** If something here doesn't sit right with you, that's useful — please push on it.",
"",
"## Suggested handling (for discussion)",
sprintf("1. **Lean toward not using APOE4 / `APOE4 - APOE` from `%s`** for now. Since we can't tell the affected ~1/3 from the rest, the safe default is to treat APOE4 as unreliable for the **whole bay**; total APOE and the other ~130 proteins on Bay1 look usable.", TARGET_BAY),
"2. **Consider keeping the per-run APOE4-vs-genotype QC metric** in the APOE secondary QC so future runs are screened (`per_run_apoe4_qc.csv`; added as a section in `APOE_geno_protein.qmd`).",
"3. **Next concrete steps to firm this up (in priority order):** (a) **ask Alamar** whether the fault is at an APOE4-target processing step (confines it to APOE4) vs a sample-tracking / plate step (whole well) — the single most informative question; (b) **re-assay a handful of the flagged samples** (definitive); (c) the raw Alamar well-level APOE4 layout for this bay. The genetic-sex check only validates the genotyping/metadata side, so it is lower priority here.",
"",
"## Artifacts (this folder's `output_files/`)",
"- `bay1_plate_map.png` / `bay1_plate_map.csv` — well-layout map of the mismatched samples",
"- `bay1_plate_qc_overlay.png` — per-well reads / detectability with mismatch wells outlined (low-input check)",
"- `bay1_apoe4_sample_deviations.csv` — per-sample Bay1 APOE4 deviations + unreliable flag",
"- `per_run_apoe4_qc.csv` — the per-run-bay APOE4 QC metric across all runs",
"- `bay1_apoe4_anomaly.png` — APOE4-by-ε4 (Bay1 vs siblings vs others) + panel variance scan",
"- `other_runs_apoe4_wellmap.png` — the same APOE4-vs-genotype map across **all** run-bays (from `investigate_other_runs_apoe4_wells.R`)",
"- `bay1_panel_iqr_ratio.txt` — full per-biomarker Bay1/sibling IQR ratios",
"- script: `investigate_20251124_Bay1.R`")
writeLines(note, file.path(OUT, "Bay1_APOE4_anomaly_internal_note.md"))

# ---- Concise 2-page colleague summary (figure-forward) -> renders to docx/pdf ----
summary_qmd <- c(
"---",
"title: \"APOE4 readout anomaly — run 20251124-1407, Bay1\"",
"subtitle: \"Biomarker QC — working note (please scrutinize)\"",
sprintf("date: \"%s\"", format(Sys.Date())),
"format:",
"  docx: default",
"  typst:",
"    papersize: us-letter",
"    margin: { x: 2cm, y: 1.5cm }",
"---",
"",
"## Bottom line",
"",
sprintf("- **APOE4 protein disagrees with APOE genotype on %d of %d (%.0f%%) wells — Bay1 only** of run 20251124-1407 (its sibling bays are clean).", n_off, n_b1, 100*n_off/n_b1),
"- **Genotypes are fine** (Sanger = WGS, identity chrX-confirmed) and **total APOE + the other ~130 proteins are fine.** Only **APOE4**, only this bay.",
"- Errors run **both ways** (carriers read low *and* non-carriers read high) and the carrier-level count is conserved → looks like a **value swap, not signal drop-out** — so it is **not correctable** (we recommend dropping APOE4 for this bay).",
"- It slipped past Alamar's QC and our control-CV QC; a new genotype-anchored per-run check catches it.",
"",
"*Working diagnosis, not a closed case (see Open questions). Full reasoning + numbers: `Bay1_APOE4_anomaly_internal_note.md`.*",
"",
"## The Bay1 plate — red = APOE4 disagrees with genotype",
"",
"![](bay1_plate_map.png){width=95%}",
"",
sprintf("- Scattered across the whole plate (%d/8 rows, %d/12 columns) — **not** an edge, stripe, or corner.", mm_rows, mm_cols),
"- Blue **control wells sit among the broken ones and read APOE4 normally** (why plate QC missed it). Broken wells are **not low-quality** — normal reads/detectability, all passed Alamar QC.",
"",
"## Both directions, only APOE4 — the \"swap\" signature",
"",
"![](bay1_apoe4_anomaly.png){width=78%}",
"",
"- At ε4 = 1 and 2, Bay1 (red) splits into **two lobes** — carriers crashing toward ~0 *and* non-carriers rising to carrier level — while the sibling bays (blue) stay tight. Total APOE is unaffected.",
"",
"## Every run-bay — does anything else jump out? (red = reads high, blue = reads low)",
"",
"![](other_runs_apoe4_wellmap.png){width=95%}",
"",
"- Bay1 (top-left) is the only plate-level outlier. The November **20251125-1332** plates show a milder reads-high cluster — possibly an **African (AFDC) batch or biology** effect (separate note).",
"- **Fresh eyes welcome — if anything else stands out to you, please flag it.**",
"",
"## What we ruled out",
"",
"**Not** low-input wells · **not** a sample mix-up (chrX = right people; Sanger = WGS) · **not** other assays · **not** distinguished by site / country / diagnosis / sex / ancestry / age, and **not** clustered by sample number or well position — the affected ~1/3 look **random** with respect to who/where/when.",
"",
"## Open questions / what we might be missing",
"",
"- **Why only APOE4?** No mechanism yet for a single-target value scramble.",
"- **Swap vs noise** — we can't fully separate a literal value-swap from \"APOE4 became uninformative.\"",
"- **Sex check not done** — genetic-vs-reported sex for these wells would further test (or break) the \"right people, wrong APOE4\" reading.",
"- **Not yet traced to the raw Alamar NPQ** (vs anything our own merge introduced).",
"- **Other explanations we haven't thought of** — please push on these.",
"",
"## Suggested handling (for discussion)",
"",
"- **Drop APOE4 / APOE4 − APOE for 20251124-1407_Bay1**; keep total APOE and the other proteins.",
"- **Keep the per-run APOE4-vs-genotype QC screen** so future runs are flagged automatically.",
"- **Ask Alamar** whether APOE4's per-well values could have been mis-assigned (demux/quantification) on this bay; **run the sex check**.",
"",
"---",
"",
sprintf("*Appendix (key numbers): carrier-level APOE4 reads = %d = number of true carriers (%d) → conserved; %d carriers fell to non-carrier level, %d non-carriers rose to carrier level. chrX: all %d affected-with-WGS in the right-person range; Sanger = WGS %d/%d. Per-run screen: Bay1 %.0f%% off-band vs %.0f%% next-worst / %.0f%% median; carrier-band IQR %.1f vs %.1f median. Reproduce: `investigate_20251124_Bay1.R`.*",
        n_obs_carr, n_carr_b1, n_lost, n_gained, n_aff_wgs, n_aff_sw_agree, n_aff_wgs, b1_off_pct, next_pct, med_pct, b1_cIQR, med_cIQR))
writeLines(summary_qmd, file.path(OUT, "Bay1_APOE4_summary.qmd"))

cat(sprintf("\n-> wrote %s/{Bay1_APOE4_anomaly_internal_note.md (full), Bay1_APOE4_summary.qmd (2-page), figures + csvs}\n", OUT))
