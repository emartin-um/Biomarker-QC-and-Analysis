#!/usr/bin/env Rscript
# =============================================================================
# investigate_20251125_African_plate.R
# -----------------------------------------------------------------------------
# Characterizes the "reads-high" APOE4 cluster on run 20251125-1332 (2nd-worst on
# the per-run APOE4 QC screen, ~6-11% off-band, all "reads HIGH"). 20251125-1332 is
# the all-African plate (AFDC: Tanzania/Nigeria/Uganda). Question: is the higher
# carrier APOE4 a RUN/batch effect or genuine African (AFDC) biology?
#
# IMPORTANT CONFOUND: every AFDC (African-born) sample is on this one run, so AFDC
# ancestry and the run cannot be separated with these data. The only "African" on
# other runs is African-American (AA) — a different, admixed population. So we can
# only bound the question, not resolve it. Output says so honestly.
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })

MM  <- "../../Metadata_Merge/output_files"
OUT <- "output_files"; if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
PLATE <- "20251125-1332"

m  <- read.csv(file.path(MM, "merged_combined_post_QC.csv"), check.names = FALSE)
fc <- read.csv(file.path(MM, "filtered/filtered_combined_post_QC.csv"), check.names = FALSE)
m  <- m %>% left_join(fc %>% select(SAMPLE_ALIQUOT, Ancestry) %>% distinct(), by = "SAMPLE_ALIQUOT")
e4 <- function(g){ g <- gsub("[^0-9]", "", trimws(g)); ifelse(g == "" | is.na(g), NA, lengths(regmatches(g, gregexpr("4", g)))) }
m <- m %>% mutate(e4 = e4(APOE.geno), plate = sub("_Bay.*", "", RUN),
                  grp = factor(case_when(
                    Ancestry == "AFDC"                 ~ "AFDC, African-born (ALL on run 20251125-1332)",
                    Ancestry == "AA"                   ~ "AA, African-American (other runs)",
                    TRUE                               ~ "Non-African (other runs)"),
                    levels = c("AFDC, African-born (ALL on run 20251125-1332)",
                               "AA, African-American (other runs)", "Non-African (other runs)")))

# ---- the confound, made explicit ---------------------------------------------
cat("=== AFDC sample location (is AFDC confounded with the run?) ===\n")
print(m %>% filter(Ancestry == "AFDC") %>% mutate(loc = ifelse(plate == PLATE, "on 20251125-1332", "other runs")) %>% count(loc))
afdc_other <- m %>% filter(Ancestry == "AFDC", plate != PLATE) %>% nrow()
cat(sprintf("  AFDC on OTHER runs: %d  => AFDC ancestry is %s with run %s.\n\n",
            afdc_other, ifelse(afdc_other == 0, "PERFECTLY CONFOUNDED", "partly separable"), PLATE))

tab <- m %>% filter(!is.na(e4), !is.na(APOE4)) %>%
  group_by(grp, e4) %>%
  summarise(n = n(), APOE4_med = round(median(APOE4), 2), APOE4_IQR = round(IQR(APOE4), 2),
            APOE_med = round(median(APOE), 2), .groups = "drop")
cat("=== APOE4 by ε4 dosage x group ===\n"); print(as.data.frame(tab))
g1 <- function(grp) tab$APOE4_med[grp == levels(m$grp)[grp] & FALSE]  # placeholder
med1 <- function(i) tab$APOE4_med[tab$grp == levels(m$grp)[i] & tab$e4 == 1]
cat(sprintf("\nHeterozygote (e4=1) APOE4 median: AFDC/1125 %.2f | AA-other %.2f | non-African-other %.2f\n",
            med1(1), med1(2), med1(3)))

# uniform across the run's bays?
bay <- m %>% filter(plate == PLATE, e4 == 1, !is.na(APOE4)) %>%
  group_by(RUN) %>% summarise(n = n(), APOE4_med = round(median(APOE4), 2), .groups = "drop")
cat("\n=== 1125 e4=1 carriers by bay ===\n"); print(as.data.frame(bay))

# ---- honest note -------------------------------------------------------------
writeLines(c(
"# 20251125-1332 (all-African / AFDC plate) — higher ε4-carrier APOE4: batch effect OR true biology (CONFOUNDED)",
"",
"**Summary.** Run `20251125-1332` is the all-African plate (248 samples, AFDC African-born: Tanzania 119, Nigeria 79, Uganda 41). On the per-run APOE4 QC screen it is 2nd-worst (~6-11% off-band), and unlike `20251124-1407_Bay1` its off-band wells all **read HIGH** (carriers above their ε4 band), not crashed.",
"",
"## The signal",
sprintf("- ε4-heterozygous (e4=1) APOE4 median **%.1f on AFDC/1125** vs **%.1f for non-Africans** — a mild **~+%.1f log2** carrier elevation; e4=2 similar. Non-carriers (~0) and **total APOE are normal**, so it is specific to the APOE4 (ε4-isoform) channel.", med1(1), med1(3), med1(1) - med1(3)),
"- The shift is roughly **uniform across all 3 bays** of the run (~11.2 / 11.5 / 11.8 for e4=1).",
"",
"## Can we tell batch from biology? **No — they are confounded.**",
sprintf("- **Every AFDC (African-born) sample sits on this one run** (AFDC on other runs: %d). So AFDC ancestry and the 20251125-1332 run are **perfectly confounded** — this dataset cannot separate them.", afdc_other),
sprintf("- The only \"African\" on other runs is **African-American (AA)**, a *different, admixed* population, which reads like non-Africans (e4=1 ≈ %.1f vs %.1f). That argues against an *admixed*-African effect but says little about **African-born (AFDC)** ε4 biology, which is untested off this plate.", med1(2), med1(3)),
"",
"## Interpretation (noted, not resolved)",
"- **Could be true AFDC biology:** AFDC is a distinct, less-admixed population; a real ε4-related APOE4 difference is plausible, and the uniform-across-bays, total-APOE-normal pattern is *consistent* with a population property (not only with a batch artifact).",
"- **Could be a run/batch offset:** it is a single run; run-level calibration offsets are common; AA (related ancestry) elsewhere does not show it.",
"- **To resolve, we would need** AFDC samples on a second run, or an orthogonal APOE4 measurement on these samples (e.g. re-assay a subset on another plate / a different platform).",
"",
"## Practical guidance",
"- **Mild** (~+1 log2 in carriers); genotype is **not** decoupled (carriers clearly carrier-level, non-carriers ~0), so carrier/non-carrier APOE4 calls remain valid — nothing like the Bay1 break.",
"- For quantitative ε4-dosage / `APOE4 - APOE` analyses, **flag run `20251125-1332` (= the AFDC samples) as a potential batch/ancestry offset** and handle with a RUN/ancestry covariate; do not treat the absolute APOE4 there as directly comparable to other runs, and do not assume it is purely biological **or** purely artifact."),
  file.path(OUT, "African_plate_20251125_batch_note.md"))

# ---- figure (violins) --------------------------------------------------------
dd <- m %>% filter(!is.na(e4), !is.na(APOE4))
pA <- ggplot(dd, aes(factor(e4), APOE4, fill = grp, colour = grp)) +
  geom_violin(scale = "width", alpha = 0.25, position = position_dodge(0.85), width = 0.8, linewidth = 0.3) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.4,
               position = position_dodge(0.85), show.legend = FALSE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.85),
             size = 0.6, alpha = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "grey55")) +
  scale_colour_manual(values = c("#E41A1C", "#4DAF4A", "grey55")) +
  labs(title = "APOE4 by ε4 dosage — AFDC/1125 carriers run high, but AFDC is confounded with the run",
       subtitle = "AA (admixed African, green) matches non-Africans (grey). AFDC (red) is higher — but ALL AFDC are on this one run, so batch vs biology can't be separated.",
       x = "ε4 allele count (clinical genotype)", y = "APOE4 (log2 NPQ)", fill = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")

pB <- m %>% filter(plate == PLATE, e4 == 1, !is.na(APOE4)) %>%
  ggplot(aes(RUN, APOE4)) +
  geom_violin(scale = "width", fill = "#FBB4AE", alpha = 0.6, width = 0.7, linewidth = 0.3) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, linewidth = 0.4, colour = "grey20") +
  geom_jitter(width = 0.1, size = 1, alpha = 0.6) +
  geom_hline(yintercept = med1(3), linetype = "dashed", colour = "grey40") +
  annotate("text", x = 0.7, y = med1(3) - 0.6, hjust = 0, label = "non-African baseline", size = 3, colour = "grey40") +
  labs(title = "Uniform across all 3 bays of the run (e4=1 carriers)",
       subtitle = "Run-level offset — consistent with either a batch calibration or an AFDC population property",
       x = NULL, y = "APOE4 (log2 NPQ)") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(OUT, "African_plate_20251125_batch.png"), pA / pB, width = 11, height = 10, dpi = 150, bg = "white")
cat(sprintf("\n-> wrote %s/{African_plate_20251125_batch.png, African_plate_20251125_batch_note.md}\n", OUT))
