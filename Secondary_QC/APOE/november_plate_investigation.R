# november_plate_investigation.R
# QUESTION (user): the November plates are "all Africans" -- maybe the
# biomarkers legitimately plot LOWER for Africans rather than the genotypes
# being wrong.
#
# CLARIFICATION from the data:
#   * 20251125-1332 IS the all-African plate (246 AFDC: Tanzania/Nigeria/Uganda),
#     but it is SANGER-ONLY (no WGS) -> it has NO genotype discordance at all.
#   * 20251124-1407 is mostly HISPANIC (HI_MU); only its Bay1 carries WGS, and
#     that Bay1 is the SEPARATE biomarker-batch outlier (~38% bio-vs-geno
#     mismatch among concordant samples), unrelated to African ancestry.
#
# So the genotype-discordance story does NOT live on the African plate. The
# only thing left to test there is the biological question: do APOE proteins
# plot lower for Africans? We answer with the TRUSTWORTHY Sanger genotype
# (no WGS required), so the full African plate is included.

suppressMessages({library(tidyverse); library(patchwork)})

fc <- read.csv("../../Metadata_Merge/output_files/filtered/filtered_combined_post_QC.csv")

e4_count <- function(g) ifelse(is.na(g) | g == "", NA_integer_,
                               lengths(regmatches(g, gregexpr("4", g))))
norm_g <- function(g) ifelse(is.na(g) | g == "", NA, gsub("[/ ]", "", trimws(g)))

d <- fc %>%
  mutate(
    sanger = norm_g(APOE.geno),
    e4 = e4_count(sanger),
    dm = APOE4 - APOE,
    nov = grepl("^20251124|^20251125", RUN),
    all_afr_plate = grepl("^20251125", RUN),
    african = Ancestry %in% c("AFDC", "AA")
  )

s <- d %>% filter(!is.na(sanger), !is.na(dm), !is.na(e4)) %>%
  mutate(grp = case_when(
    all_afr_plate              ~ "African plate (20251125)",
    african & !nov             ~ "African, other plates",
    !african & !nov            ~ "Non-African, other plates",
    TRUE                       ~ "Other November"),
    e4_lab = paste0(e4, " ε4"))

# medians for the answer table
cat("=== APOE4-APOE medians by e4 band x group (Sanger genotype) ===\n")
print(s %>% group_by(e4, grp) %>%
        summarise(n = n(), median_dm = round(median(dm), 2), .groups = "drop") %>%
        arrange(e4, grp) %>% as.data.frame(), row.names = FALSE)

# ── PLOT 1: the requested "gradient" — APOE4-APOE by e4 dosage, by group ──────
band_meds <- s %>% group_by(e4_lab) %>% summarise(med = median(dm), .groups = "drop")

p1 <- ggplot(s, aes(x = grp, y = dm, colour = grp)) +
  geom_hline(data = band_meds, aes(yintercept = med),
             linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6) +
  geom_jitter(width = 0.18, size = 0.9, alpha = 0.45) +
  facet_wrap(~ e4_lab, scales = "free_y") +
  labs(title = "APOE4 − APOE by ε4 dosage: is the African plate lower?",
       subtitle = paste("Trustworthy Sanger genotype (no WGS required).",
                         "Dashed line = overall band median.",
                         "African plate is NOT lower — if anything slightly higher in carrier bands."),
       x = NULL, y = "APOE4 − APOE (log2)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1))

# ── PLOT 2: African ancestry across ALL plates (does ancestry track lower?) ───
anc <- s %>% mutate(anc = ifelse(african, "African (AFDC/AA)", "Non-African"))
p2 <- ggplot(anc, aes(x = factor(e4), y = dm, colour = anc)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 0.8, alpha = 0.4) +
  labs(title = "Pooled across all plates: African vs non-African by ε4 band",
       subtitle = "ε4-dosage separation is preserved in both ancestry groups; medians nearly overlap.",
       x = "ε4 allele count (Sanger genotype)", y = "APOE4 − APOE (log2)",
       colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")

ggsave("november_african_vs_plate.png", p1, width = 11, height = 4.8, dpi = 150)
ggsave("african_ancestry_pooled.png", p2, width = 8.5, height = 5, dpi = 150)
cat("\nSaved november_african_vs_plate.png and african_ancestry_pooled.png\n")
cat("\nDONE\n")
