# investigate_discordance.R
# What — beyond biomarker RUN date — distinguishes Sanger-vs-WGS discordant samples?
# RUN date is the biomarker assay date and cannot CAUSE genotype-vs-genotype
# discordance, so we look for sample-level attributes that travel with the
# discordant set. We also separately test the genuine biomarker-batch signal:
# biomarker-vs-genotype disagreement among CONCORDANT (trustworthy) samples.

suppressMessages({library(tidyverse)})

fc_path  <- "../../Metadata_Merge/output_files_n2851_2026May/filtered/filtered_combined_post_QC.csv"
old_path <- "U19_Alamar_metadata_2026Jan16_with_WGSAPOE.csv"

fc  <- read_csv(fc_path,  show_col_types = FALSE)
old <- read_csv(old_path, show_col_types = FALSE)

e4_count <- function(g) ifelse(is.na(g) | g == "", NA_integer_,
                               lengths(regmatches(g, gregexpr("4", g))))

# old lookup: which SAMPLE_ALIQUOT had a WGS call in the Jan-2026 502-sample file
old_wgs_col <- intersect(c("APOE_WGS","APOE_WGS_norm","WGS_APOE"), names(old))[1]
old_wgs_ids <- old %>% filter(!is.na(.data[[old_wgs_col]])) %>% pull(SAMPLE_ALIQUOT) %>% unique()
cat("Old (Jan2026) WGS-typed samples:", length(old_wgs_ids), "\n\n")

sw <- fc %>%
  mutate(
    sanger   = na_if(gsub("[/ ]", "", trimws(APOE.geno)), ""),
    wgs      = na_if(gsub("[/ ]", "", trimws(APOE_WGS)),  ""),
    samp_num = suppressWarnings(as.numeric(sub("-.*", "", SAMPLE_ALIQUOT))),
    e4_sanger= e4_count(sanger),
    e4_wgs   = e4_count(wgs),
    call_category = case_when(
      !is.na(sanger) & !is.na(wgs) & sanger == wgs ~ "Both agree",
      !is.na(sanger) & !is.na(wgs) & sanger != wgs ~ "Both disagree",
      !is.na(sanger) &  is.na(wgs)                 ~ "Sanger only",
       is.na(sanger) & !is.na(wgs)                 ~ "WGS only",
      TRUE                                          ~ "Neither"),
    wgs_newly_added = !(SAMPLE_ALIQUOT %in% old_wgs_ids) & !is.na(wgs),
    APOE4_minus_APOE = APOE4 - APOE
  )

dual <- sw %>% filter(call_category %in% c("Both agree","Both disagree")) %>%
  mutate(disc = call_category == "Both disagree")

cat("Dual-typed samples:", nrow(dual),
    " | discordant:", sum(dual$disc),
    " | concordant:", sum(!dual$disc), "\n\n")

# ── Helper: compare a categorical attribute discordant vs concordant ─────────
compare_cat <- function(df, col) {
  if (!col %in% names(df)) { cat("  (no column", col, ")\n"); return(invisible()) }
  tab <- df %>%
    mutate(grp = ifelse(disc, "discordant", "concordant")) %>%
    count(.data[[col]], grp) %>%
    pivot_wider(names_from = grp, values_from = n, values_fill = 0) %>%
    mutate(
      concordant = if (!"concordant" %in% names(.)) 0L else concordant,
      discordant = if (!"discordant" %in% names(.)) 0L else discordant,
      total      = concordant + discordant,
      pct_disc   = round(100 * discordant / total, 1)
    ) %>%
    arrange(desc(pct_disc))
  cat("\n## ", col, " (discordance rate by level)\n", sep = "")
  print(as.data.frame(tab), row.names = FALSE)
}

for (col in c("metatada_source","wgs_newly_added","Ancestry","Race",
              "Ethnicity","Site","Country/State","CDX","CDX_collapsed")) {
  compare_cat(dual, col)
}

# ── SAMPLE number range ───────────────────────────────────────────────────────
cat("\n\n## SAMPLE number (accession order) — discordant vs concordant\n")
print(dual %>% group_by(disc) %>%
        summarise(n=n(), min=min(samp_num,na.rm=TRUE),
                  q25=quantile(samp_num,.25,na.rm=TRUE),
                  median=median(samp_num,na.rm=TRUE),
                  q75=quantile(samp_num,.75,na.rm=TRUE),
                  max=max(samp_num,na.rm=TRUE), .groups="drop") %>%
        as.data.frame(), row.names = FALSE)

# ── Newly-added WGS is the key suspect: cross newly_added × discordance ───────
cat("\n## newly-added WGS  ×  discordance\n")
print(dual %>% count(wgs_newly_added, disc) %>%
        pivot_wider(names_from = disc, values_from = n, values_fill = 0) %>%
        as.data.frame(), row.names = FALSE)

# ── RUN: confound check. Is RUN associated with newly_added & metadata_source? ─
cat("\n\n## RUN × discordance (top 15 by discordant count) — the CONFOUND\n")
print(dual %>% count(RUN, disc) %>%
        pivot_wider(names_from = disc, values_from = n, values_fill = 0) %>%
        rename(concordant=`FALSE`, discordant=`TRUE`) %>%
        mutate(total=concordant+discordant, pct=round(100*discordant/total,1)) %>%
        arrange(desc(discordant)) %>% head(15) %>%
        as.data.frame(), row.names = FALSE)

cat("\n## Within the high-discordance RUNs, is it driven by metadata_source / newly_added?\n")
hi_runs <- dual %>% count(RUN, disc) %>% filter(disc) %>% arrange(desc(n)) %>% head(6) %>% pull(RUN)
print(dual %>% filter(RUN %in% hi_runs) %>%
        count(RUN, metatada_source, wgs_newly_added, disc) %>%
        arrange(RUN, desc(disc)) %>% as.data.frame(), row.names = FALSE)

# ── Is RUN still predictive AFTER conditioning on ancestry? ──────────────────
# If the Oct-2025 runs are discordant ONLY because they're Hispanic-heavy, then
# within Hispanic samples the Oct runs should look like other Hispanic runs.
cat("\n\n## Ancestry composition of the high-discordance Oct-2025 runs vs the rest\n")
print(dual %>%
        mutate(oct = grepl("^202510", RUN),
               hispanic = Ethnicity == "HI") %>%
        count(oct, hispanic) %>%
        pivot_wider(names_from = hispanic, values_from = n, values_fill = 0) %>%
        as.data.frame(), row.names = FALSE)

cat("\n## Discordance rate by RUN-era, restricted to HISPANIC samples only\n")
print(dual %>% filter(Ethnicity == "HI") %>%
        mutate(era = ifelse(grepl("^202510", RUN), "Oct-2025 runs", "other runs")) %>%
        group_by(era) %>%
        summarise(n=n(), discordant=sum(disc),
                  pct=round(100*mean(disc),1), .groups="drop") %>%
        as.data.frame(), row.names = FALSE)

cat("\n## Discordance rate by RUN-era, restricted to NON-Hispanic samples only\n")
print(dual %>% filter(Ethnicity == "NH") %>%
        mutate(era = ifelse(grepl("^202510", RUN), "Oct-2025 runs", "other runs")) %>%
        group_by(era) %>%
        summarise(n=n(), discordant=sum(disc),
                  pct=round(100*mean(disc),1), .groups="drop") %>%
        as.data.frame(), row.names = FALSE)

# ── SEPARATE TEST: genuine biomarker-batch effect ───────────────────────────
# Among CONCORDANT (trustworthy genotype) samples, does the biomarker-implied
# e4 (nearest reference median of APOE4_minus_APOE) disagree with the agreed
# genotype more often in certain RUNs/Bays? THAT is what a biomarker batch
# effect would actually look like.
cat("\n\n#### SEPARATE: biomarker-vs-genotype disagreement among CONCORDANT samples\n")
if ("APOE4_minus_APOE" %in% names(sw)) {
  ref <- sw %>% filter(call_category=="Both agree", !is.na(APOE4_minus_APOE)) %>%
    group_by(e4=e4_sanger) %>% summarise(med=median(APOE4_minus_APOE,na.rm=TRUE),.groups="drop")
  med <- setNames(ref$med, ref$e4)
  bio_e4 <- function(x){
    if (all(is.na(x))) return(rep(NA_integer_,length(x)))
    d <- sapply(names(med), function(k) abs(x - med[[k]]))
    as.integer(names(med)[max.col(-d, ties.method="first")])
  }
  conc <- sw %>% filter(call_category=="Both agree", !is.na(APOE4_minus_APOE)) %>%
    mutate(e4_bio = bio_e4(APOE4_minus_APOE),
           bio_mismatch = e4_bio != e4_sanger)
  cat("Concordant samples with biomarker data:", nrow(conc),
      " | biomarker-vs-genotype mismatches:", sum(conc$bio_mismatch,na.rm=TRUE),
      sprintf(" (%.1f%%)\n", 100*mean(conc$bio_mismatch,na.rm=TRUE)))
  cat("\nBiomarker-mismatch rate by RUN (top 15 by mismatch count):\n")
  print(conc %>% count(RUN, bio_mismatch) %>%
          pivot_wider(names_from=bio_mismatch, values_from=n, values_fill=0) %>%
          rename(match=`FALSE`, mismatch=`TRUE`) %>%
          mutate(total=match+mismatch, pct=round(100*mismatch/total,1)) %>%
          arrange(desc(mismatch)) %>% head(15) %>% as.data.frame(), row.names=FALSE)
  if ("Bay" %in% names(conc)) {
    cat("\nBiomarker-mismatch rate by Bay:\n")
    print(conc %>% count(Bay, bio_mismatch) %>%
            pivot_wider(names_from=bio_mismatch, values_from=n, values_fill=0) %>%
            rename(match=`FALSE`, mismatch=`TRUE`) %>%
            mutate(total=match+mismatch, pct=round(100*mismatch/total,1)) %>%
            arrange(desc(pct)) %>% as.data.frame(), row.names=FALSE)
  }
}
cat("\nDONE\n")
