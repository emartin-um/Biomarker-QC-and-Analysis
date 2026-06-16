#!/usr/bin/env Rscript
# =============================================================================
# investigate_other_runs_apoe4_wells.R
# -----------------------------------------------------------------------------
# Companion to investigate_20251124_Bay1.R. Maps the APOE4-vs-genotype mismatch
# ("off-band") wells on EVERY run-bay (well positions from Alamar Sample_QC),
# to check whether any OTHER plate shows a Bay1-like spatial/well pattern or
# just sparse, scattered background outliers. Splits the mismatches by direction
# (carrier crashed LOW vs read HIGH).
#
# Writes: other_runs_apoe4_mismatch_wells.csv  (every off-band well, all runs)
#         other_runs_apoe4_wellmap.png         (faceted plate maps, off>=2 runs)
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

MM  <- "../../Metadata_Merge/output_files"
SQC <- "../../Primary_QC/input_files/Sample_QC.csv"
OUT <- "output_files"; if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

m <- read.csv(file.path(MM, "merged_combined_post_QC.csv"), check.names = FALSE)
e4 <- function(g){ g <- gsub("[^0-9]", "", trimws(g)); ifelse(g == "" | is.na(g), NA, lengths(regmatches(g, gregexpr("4", g)))) }
m$e4 <- e4(m$APOE.geno)
band <- m %>% filter(!is.na(e4), !is.na(APOE4)) %>% group_by(e4) %>% summarise(md = median(APOE4), .groups = "drop")
m <- m %>% mutate(dev = APOE4 - band$md[match(e4, band$e4)],
                  off = abs(dev) > 2 & !is.na(e4),
                  direction = ifelse(dev < 0, "carrier crashed LOW", "reads HIGH"))

sq <- suppressWarnings(readr::read_csv(SQC, skip = 12, show_col_types = FALSE)) %>%
  transmute(SAMPLE_ALIQUOT = `Sample Name`, sample_type = `Sample Type`,
            row = wellRow, col = as.integer(wellCol))

dat <- m %>% filter(!is.na(e4)) %>%
  transmute(RUN, SAMPLE_ALIQUOT, e4, APOE.geno, APOE4 = round(APOE4, 1), dev = round(dev, 1), off, direction) %>%
  left_join(sq, by = "SAMPLE_ALIQUOT")

# ---- table of every off-band well across all runs ----
offtab <- dat %>% filter(off) %>%
  transmute(RUN, well = paste0(row, col), row, col, APOE.geno, APOE4, dev, direction) %>%
  arrange(RUN, row, col)
write.csv(offtab, file.path(OUT, "other_runs_apoe4_mismatch_wells.csv"), row.names = FALSE)

# ---- per run-bay summary + spatial dispersion (no clustering => max-in-row/col small) ----
summ <- dat %>% group_by(RUN) %>%
  summarise(n = n(), n_off = sum(off),
            off_pct = round(100 * mean(off), 1),
            crashed = sum(off & direction == "carrier crashed LOW"),
            high    = sum(off & direction == "reads HIGH"),
            rows_hit = length(unique(row[off])), cols_hit = length(unique(col[off])),
            max_in_a_row = ifelse(any(off), max(table(row[off])), 0),
            max_in_a_col = ifelse(any(off), max(table(col[off])), 0),
            .groups = "drop") %>%
  filter(n >= 15) %>% arrange(desc(n_off))
cat("=== off-band APOE4 wells per run-bay (study samples, n>=15) ===\n")
print(as.data.frame(summ), row.names = FALSE)
cat(sprintf("\nBay1 = %d off (%.0f%%); every other run-bay <= %d off; none cluster (max-in-row/col <= %d).\n",
            summ$n_off[summ$RUN == "20251124-1407_Bay1"], summ$off_pct[summ$RUN == "20251124-1407_Bay1"],
            max(summ$n_off[summ$RUN != "20251124-1407_Bay1"]),
            max(summ$max_in_a_row[summ$RUN != "20251124-1407_Bay1"],
                summ$max_in_a_col[summ$RUN != "20251124-1407_Bay1"])))

# ---- faceted plate maps for run-bays with >=2 off-band wells ----
runs_show <- summ %>% filter(n_off >= 2) %>%
  mutate(facet = sprintf("%s  (%d off, %.0f%%)", RUN, n_off, off_pct))
grid_all <- dat %>% filter(RUN %in% runs_show$RUN) %>%
  left_join(runs_show %>% select(RUN, facet, n_off), by = "RUN") %>%
  mutate(facet = fct_reorder(facet, -n_off),
         y = 9 - match(row, LETTERS[1:8]),
         cls = case_when(!off ~ "ok / control", direction == "carrier crashed LOW" ~ "crashed LOW", TRUE ~ "reads HIGH"),
         cls = factor(cls, levels = c("crashed LOW", "reads HIGH", "ok / control")))

p <- ggplot(grid_all, aes(col, y)) +
  geom_tile(aes(fill = cls), color = "grey90", width = 0.95, height = 0.95) +
  facet_wrap(~ facet, ncol = 4) +
  scale_fill_manual(values = c("crashed LOW" = "#377EB8", "reads HIGH" = "#E41A1C",
                               "ok / control" = "grey88")) +
  scale_x_continuous(breaks = c(1, 6, 12)) +
  scale_y_continuous(breaks = c(1, 8), labels = c("H", "A")) +
  coord_equal(expand = FALSE) +
  labs(title = "APOE4-vs-genotype mismatch wells across run-bays (Alamar well layout)",
       subtitle = "Bay1 (top-left) is the only plate-level failure (32%); elsewhere just 2-9 scattered wells.\nblue = carrier reads low,  red = reads high,  grey = agrees / control.",
       x = "column", y = "row", fill = "APOE4 mismatch") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        strip.text = element_text(size = 7))
ggsave(file.path(OUT, "other_runs_apoe4_wellmap.png"), p, width = 12, height = 9, dpi = 150, bg = "white")
cat(sprintf("\n-> wrote %s/other_runs_apoe4_mismatch_wells.csv (%d wells) + other_runs_apoe4_wellmap.png (%d run-bays)\n",
            OUT, nrow(offtab), nrow(runs_show)))
