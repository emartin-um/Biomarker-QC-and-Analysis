# biomarker_batch_run_investigation.R
# Characterize the biomarker-assay anomaly on run 20251124-1407_Bay1.
# The biomarker-vs-genotype "mismatch" is computed from APOE4 - APOE, so we ask
# WHICH protein breaks. Finding: APOE4's variance explodes on this bay (carrier
# IQR ~9x), collapsing the carrier/non-carrier separation the e4 inference needs;
# APOE (total) and the rest of the 132-plex panel are essentially normal -> a
# target-specific APOE4 readout failure, not a sample mix-up or plate-wide drift.

suppressMessages({library(tidyverse); library(patchwork)})

fc <- read.csv("../../Metadata_Merge/output_files/filtered/filtered_combined_post_QC.csv",
               check.names = FALSE)
norm_g <- function(g) ifelse(is.na(g) | g == "", NA, gsub("[/ ]", "", trimws(g)))
e4c    <- function(g) ifelse(is.na(g) | g == "", NA_integer_,
                             lengths(regmatches(g, gregexpr("4", g))))
TARGET <- "20251124-1407_Bay1"

d <- fc %>% mutate(
  sanger = norm_g(APOE.geno), wgs = norm_g(APOE_WGS),
  e4 = e4c(sanger), dm = APOE4 - APOE,
  concord = !is.na(sanger) & !is.na(wgs) & sanger == wgs,
  target = RUN == TARGET)

# ── biomarker-implied e4 (nearest reference band median, from concordant) ─────
ref <- d %>% filter(concord, !is.na(dm), !is.na(e4)) %>%
  group_by(k = e4) %>% summarise(med = median(dm), .groups = "drop")
med <- setNames(ref$med, ref$k)
bio_e4 <- function(x) {
  dmat <- sapply(names(med), function(k) abs(x - med[[k]]))
  as.integer(names(med)[max.col(-dmat, ties.method = "first")])
}
conc <- d %>% filter(concord, !is.na(dm), !is.na(e4)) %>%
  mutate(e4_bio = bio_e4(dm), mm = e4_bio != e4)

# ════════ PLOT A: mismatch rate by RUN, target highlighted ════════
run_mm <- conc %>% group_by(RUN) %>%
  summarise(n = n(), mismatch = sum(mm), pct = 100 * mean(mm), .groups = "drop") %>%
  filter(n >= 15) %>% arrange(pct) %>%
  mutate(RUN = factor(RUN, levels = RUN), hl = RUN == TARGET)

pA <- ggplot(run_mm, aes(RUN, pct, fill = hl)) +
  geom_col(width = 0.75) +
  geom_text(data = subset(run_mm, hl),
            aes(label = sprintf("%.0f%% (%d/%d)", pct, mismatch, n)),
            hjust = -0.1, size = 3.3, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#E41A1C"), guide = "none") +
  expand_limits(y = max(run_mm$pct) * 1.15) +
  labs(title = "Biomarker–genotype mismatch by run (concordant samples)",
       subtitle = paste0(TARGET, " stands far above every other run (n >= 15)."),
       x = NULL, y = "Mismatch rate (%)") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 6.5))

# ════════ PLOT B: the mild Bay1>Bay2>Bay3 gradient (pooled) ════════
bay_mm <- conc %>% filter(!is.na(Bay)) %>% group_by(Bay) %>%
  summarise(n = n(), mismatch = sum(mm), pct = 100 * mean(mm), .groups = "drop")
pB <- ggplot(bay_mm, aes(Bay, pct, fill = Bay)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct, mismatch, n)),
            vjust = -0.2, size = 3.2) +
  expand_limits(y = max(bay_mm$pct) * 1.25) +
  labs(title = "Mild Bay gradient (all runs pooled)",
       x = NULL, y = "Mismatch rate (%)") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

# ════════ PLOT C: APOE4 variance blow-up (the mechanism) ════════
dc <- d %>% filter(!is.na(e4)) %>%
  mutate(grp = ifelse(target, paste0(TARGET), "All other runs"),
         e4_lab = paste0(e4, " ε4"))
pC <- ggplot(dc, aes(factor(e4), APOE4, colour = grp)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), width = 0.65) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
             size = 1, alpha = 0.5) +
  scale_colour_manual(values = setNames(c("#E41A1C", "grey55"),
                                        c(TARGET, "All other runs"))) +
  labs(title = "APOE4 by ε4 dosage: the assay loses dynamic range on this run",
       subtitle = "APOE4 scatters wildly on the target run (carrier IQR ~9x), erasing the carrier/non-carrier separation.",
       x = "ε4 allele count (genotype)", y = "APOE4 (log2 NPQ)", colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")

# ════════ PLOT D: APOE4 is the outlier — variance ratio across the panel ════════
# Residualize every biomarker on e4 band first, so APOE4's genotype-driven
# bimodality (identical on all runs) is removed and only measurement scatter
# remains. For non-genotype proteins, band centering is harmless.
cand <- names(fc)[28:166]
bm <- setdiff(cand[sapply(fc[cand], is.numeric)],
              c("APOE4_count", "APOE4_carrier", "APOE_WGS", "APOE_WGS_norm",
                "APOE_WGS_changed", "APOE.geno_final"))
dr <- d %>% filter(!is.na(e4))
resid_band <- function(col) {
  x <- dr[[col]]; out <- x
  for (g in unique(dr$e4)) { idx <- dr$e4 == g; out[idx] <- x[idx] - median(x[idx], na.rm = TRUE) }
  out
}
vr <- sapply(bm, function(c) {
  r <- resid_band(c)
  IQR(r[dr$target], na.rm = TRUE) / IQR(r[!dr$target], na.rm = TRUE)
})
vrdf <- tibble(bm = bm, ratio = vr) %>% filter(is.finite(ratio)) %>%
  mutate(hl = bm %in% c("APOE4", "APOE")) %>% arrange(desc(ratio))
labdf <- vrdf %>% mutate(rank = row_number()) %>% filter(hl)
pD <- ggplot(vrdf %>% mutate(rank = row_number()),
             aes(rank, ratio, colour = hl)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_point(size = 1.6, alpha = 0.7) +
  ggrepel::geom_text_repel(data = labdf, aes(label = bm),
                           size = 3.4, fontface = "bold", box.padding = 0.6) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#E41A1C"), guide = "none") +
  labs(title = "Within-run spread vs the rest of the panel (132 proteins)",
       subtitle = "Ratio of target-run IQR to other-runs IQR. APOE4 has by far the most inflated variance.",
       x = "Biomarker (ranked by variance ratio)", y = "IQR ratio (target / others)") +
  theme_bw(base_size = 11)

ggsave("biomarker_batch_run_anomaly.png", (pA | (pB / pC)) + plot_layout(widths = c(1, 1.1)),
       width = 14, height = 8, dpi = 150)
ggsave("apoe4_variance_panel.png", pD, width = 9, height = 5, dpi = 150)
ggsave("apoe4_mechanism.png", pC, width = 9, height = 5, dpi = 150)

cat("APOE4 variance ratio (target/others):",
    round(vrdf$ratio[vrdf$bm == "APOE4"], 1), "x | rank",
    which(vrdf$bm == "APOE4"), "of", nrow(vrdf), "\n")
cat("APOE  variance ratio:", round(vrdf$ratio[vrdf$bm == "APOE"], 1), "x | rank",
    which(vrdf$bm == "APOE"), "\n")
cat("Saved 3 PNGs.\nDONE\n")
