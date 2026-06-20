#!/usr/bin/env Rscript
# OPTIONAL tool. By default APOE_perbatch_QC.Rmd §A bands from the batch's OWN clean bays (a two-pass
# that excludes scrambled bays) — no cross-batch reference is used, because APOE4-APOE is not perfectly
# comparable batch-to-batch (e.g. the 8-plate carriers run ~+1 log2 vs ALZ123, which would over-flag).
# Use THIS only for the rare single-scrambled-plate case with no clean bay to anchor on, and build it
# from a SAME-NORMALIZATION clean batch; then pass the output via params$ref_bands.
# Output: reference_apoe4_bands.csv (keyed by APOEgeno, the cleaned 2-digit genotype).
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
SRC  <- if (length(args) >= 1) args[1] else
  "/Users/emartin1/Library/CloudStorage/Box-Box/AD/BiomarkerProject2024/NewAlamarDataAugust2025/Alamar_QC_Protocol/QC_Runs/QC_ALZ123_repro_2026Jun/Metadata_Merge/output_files_ALZ123_repro_remap_2026Jun/filtered/filtered_combined_post_QC.csv"
# bays known to be APOE4-scrambled — EXCLUDE from the reference (would inflate the bands)
EXCLUDE_RUN <- c("20251124-1407_Bay1")
MIN_N <- 20
out  <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "reference_apoe4_bands.csv")
if (is.na(out) || !nzchar(out)) out <- "reference_apoe4_bands.csv"

d <- read.csv(SRC, check.names = FALSE)
.clean <- function(x){ y <- gsub("\\s|/", "", trimws(as.character(x))); y[y==""] <- NA; y }
d$APOEgeno         <- .clean(d$APOE.geno)
d$APOE4_minus_APOE <- d$APOE4 - d$APOE
ref <- d %>% filter(!is.na(APOEgeno), !is.na(APOE4_minus_APOE),
                    !("RUN" %in% names(d) & RUN %in% EXCLUDE_RUN))
tab <- ref %>% group_by(APOEgeno) %>%
  summarise(n_ref = n(),
            Q1 = quantile(APOE4_minus_APOE, .25), Q3 = quantile(APOE4_minus_APOE, .75),
            IQR = Q3 - Q1, lower = Q1 - 3*IQR, upper = Q3 + 3*IQR,
            genotype_median = median(APOE4_minus_APOE), .groups = "drop") %>%
  filter(n_ref >= MIN_N) %>%
  mutate(across(c(Q1,Q3,IQR,lower,upper,genotype_median), ~round(.,3)),
         source = paste0("ALZ123_clean(", basename(SRC), ") excl ", paste(EXCLUDE_RUN, collapse=","),
                         "; built ", format(Sys.Date())))
write.csv(tab[, c("APOEgeno","n_ref","lower","upper","genotype_median","source")], out, row.names = FALSE)
cat("reference rows (genotypes with n>=", MIN_N, "):\n", sep=""); print(tab[,c("APOEgeno","n_ref","lower","upper","genotype_median")])
cat("\nwrote ", out, "\n", sep="")
