# Hemolysis Findings — Summary for Collaborator Discussion

**Prepared:** 2026-05-05
**Data version:** v4 (n=2621 samples, post-QC; full source `U19_otherprojects_metadata_combined_n2871_May2026.csv`)
**Source analysis:** `Alamar_Biomarker_QC_Repo/Secondary_QC/Hemolysis_Check/Hemolysis_by_Site.html`

---

## TL;DR

We found **systematic, large-magnitude differences in hemolysis levels across collection sites** that are confounded with our ancestry/country comparisons. Fifteen biomarkers (HBA1 + 14 others) are demonstrably affected, and ten of them appear in the headline "ancestry-significant" Paper 1 Table 1 results. **Some apparent ancestry differences are partly or wholly pre-analytical artifacts.**

We need protocol-level decisions before publication:
1. How to align collection/processing between sites going forward.
2. How to flag or correct affected biomarkers in the current dataset.
3. Whether to add a post-hoc QC criterion that excludes severely hemolyzed samples.

---

## 1. The smoking gun — Hemolysis Index by Site

The Hemolysis Index (mean of HBA1, PGK1, MDH1, SOD1, ENO2 NPQ values) varies by **~4 log-units across sites**, far beyond technical variance.

| Site | n | Country/Region | Mean Hemo Index | Mean HBA1 |
|---|---|---|---|---|
| **MULH** | 41 | Uganda | **14.2** | 12.0 |
| MAKU | 3 | Uganda | 12.9 | 8.8 |
| DUK | 45 | US (Duke) | 12.9 | 10.2 |
| NC | 286 | US | 12.7 | 9.6 |
| CU | 507 | US (Columbia) | 12.7 | 9.2 |
| IHG | 247 | US (Miami) | 12.6 | 8.6 |
| CMUL | 33 | Nigeria (Lagos) | 12.5 | 8.2 |
| CUMC | 209 | US (Columbia) | 12.4 | 8.7 |
| CWR | 160 | US (Case Western) | 12.4 | 9.2 |
| ABTH | 46 | Nigeria (Abuja) | 12.3 | 10.3 |
| IHG HI | 428 | US Hispanic | 12.2 | 7.9 |
| AAU | 5 | Ethiopia | 11.6 | 6.3 |
| IHG AA | 263 | US AA | 11.5 | 7.2 |
| **PER** | 230 | Peru | **10.7** | 7.7 |
| **MNPH** | 117 | Tanzania | **10.3** | 6.9 |

The TZA-vs-UGA gap (10.3 → 14.2) is the largest gradient in the data and corresponds **exactly** to the within-AFDC heterogeneity we observed in Paper 1 results (TZA samples consistently lowest, UGA highest, for all hemolysis-correlated biomarkers).

---

## 2. Affected biomarkers (15 total, including HBA1 reference)

These are biomarkers whose plasma levels correlate strongly with HBA1 *and* whose sample-level variance exceeds technical-replicate variance — i.e., the variation is real but driven by hemolysis.

| Biomarker | r(HBA1) | Tech CV | Sample CV | Excess CV | Currently in Table 1 sig? |
|---|---|---|---|---|---|
| HBA1 | 1.00 | 1.45 | 29.9 | 28.5 | reference |
| PRDX6 | 0.62 | 2.04 | 11.0 | 8.9 |  |
| MDH1 | 0.50 | 2.05 | 8.6 | 6.5 | ✓ |
| **GOT1** | **0.48** | 2.05 | 4.7 | 2.6 | (high r, low excess CV — borderline) |
| SOD1 | 0.48 | 2.03 | 9.7 | 7.7 |  |
| ARSA | 0.47 | 1.69 | 13.2 | 11.5 | ✓ (5th largest F) |
| SNCA | 0.45 | 1.12 | 6.8 | 5.7 | ✓ |
| PGK1 | 0.43 | 5.00 | 19.2 | 14.2 | ✓ |
| S100A12 | 0.39 | 1.03 | 16.1 | 15.1 | ✓ |
| IL16 | 0.38 | 1.02 | 9.6 | 8.6 | ✓ |
| CXCL8 | 0.35 | 1.16 | 16.8 | 15.7 |  |
| S100B | 0.34 | 1.18 | 8.8 | 7.6 | ✓ |
| ENO2 | 0.34 | 1.98 | 8.8 | 6.8 | ✓ |
| PSEN1 | 0.34 | 2.13 | 9.4 | 7.2 | ✓ |
| HTT | 0.32 | 1.02 | 9.4 | 8.4 | ✓ |
| PARK7 | 0.31 | 1.27 | 12.9 | 11.6 |  |

**Compared to prior list (Feb 2026 hemolysis run):** PSEN1, HTT, PARK7 are newly flagged in v4. GOT1 borderline.

**Of the 15 flagged: 10 appear in Paper 1 Table 1 v4** (the headline 64 ancestry-significant biomarkers): ARSA, ENO2, SNCA, HTT, S100A12, MDH1, PGK1, PSEN1, S100B, IL16.

---

## 3. What this likely means

### Pre-analytical (sample-handling) artifact — most likely
Site-level differences in hemolysis with **no consistent sex effect** (overall M=12.2, F=12.2 — see §4) and high inter-marker correlations point to **collection or processing differences** rather than biology. Candidates:
- Time from blood draw to centrifugation
- Centrifuge speed/temperature
- Needle gauge / draw technique
- Number of freeze-thaw cycles
- Storage temperature/duration before shipment to HIHG

### Biological hemolysis — partial contribution possible
Some sex-by-site variability suggests possible G6PD-deficiency contribution at specific sites (Nigeria, parts of Africa), which could reflect *true in-vivo* hemolysis. However, the overall Sex × Hemolysis effect is small and inconsistent. Worth a deeper look at G6PD genotype in WGS data.

### Why this matters for Paper 1
- The "Africa > Americas" finding for some markers is **partly explained by hemolysis** rather than biology. ARSA, the 5th-largest-F-stat finding in Table 1, is an example.
- The within-AFDC heterogeneity we documented (TZA vs NGA vs UGA differing on 84/131 biomarkers) is **largely a within-Africa hemolysis-by-site gradient**, not African genetic heterogeneity.
- We cannot publish the current Table 1 as-is without addressing this, per scientific honesty.

---

## 4. Sex effect (briefly checked, no smoking gun)

Overall: F=12.2, M=12.2, U=12.8 (n=2 only). At the site level, max M-vs-F difference is 0.7 log-units (ABTH; M lower) — not the consistent males-high-females-low pattern that G6PD deficiency would predict.

This **mostly rules out G6PD as the dominant driver** at the population level — though it doesn't rule out G6PD contributions in specific subsets, and we have not yet checked G6PD genotype directly.

---

## 5. Recommended discussion items for collaborators

### Protocol decisions (going forward)
1. **Standardize blood-collection protocol across sites.** Document exact draw-to-centrifuge time, centrifuge parameters, aliquot/freeze procedure, and shipping conditions. Prioritize sites with highest hemolysis (MULH/Uganda).
2. **Visit-day hemolysis assessment.** Add a quick visual hemolysis score at the collection site (free hemoglobin in serum is straightforward — clear/yellow/pink/red index). This data should travel with the sample.
3. **Re-collect a subset?** Especially MULH samples — their Hemo Index of 14.2 is so high that some assays may be unreliable.

### Analytic decisions (current dataset)
1. **Three options for handling 15 flagged biomarkers:**
   - **(a) Drop them entirely** from Paper 1 — clean but loses some real signal (e.g., GFAP-related biology).
   - **(b) Adjust for HBA1** as a covariate in the linear models — the standard approach. Risks over-correcting and removing real biology if hemolysis itself differs by ancestry for biological reasons (G6PD).
   - **(c) Two-track reporting** — present Table 1 with and without flagged biomarkers as parallel result sets and let the reader compare. Most transparent.
2. **Add Site as a covariate** to all Paper 1 models. This will deflate within-AFDC contribution but is more honest.
3. **Drop-one-country sensitivity** (drop TZA, NGA, UGA in turn) — biomarkers that survive in all three drops are the robust ancestry findings.
4. **Sample-level outlier exclusion** — flag samples in upper tail of Hemolysis Index (e.g., >mean+3SD) and re-run with them excluded. Numbers TBD; depends on threshold.

### G6PD follow-up (if biological hemolysis is suspected)
- Pull G6PD variants from WGS and check carrier frequency by site.
- Specifically: G6PD A− (rs1050828) — common in West African ancestry; G6PD Med (rs5030868) — Mediterranean; G6PD Cosenza — small numbers.
- If carriers cluster by site, that supports a partial biological component to the site-level hemolysis pattern.

### What I (Eden) recommend going in
- Go with **(c) two-track reporting** in Paper 1 (most transparent).
- Add **Site as covariate** in primary models, not as sensitivity. (Note: HBA1-as-covariate did *not* work — see §6.2.)
- **Drop-one-country sensitivity is reassuring** (93/98 robust, only 5 borderline) — see §6.1.
- Tighten **future collection protocol** with one harmonized SOP across all collection sites.
- **Sample doubling expected ~early June 2026** — borderline findings should be flagged as "may shift on full data" rather than treated as final.

---

## 6. Sensitivity analyses (added 2026-05-05)

Two follow-up analyses, both run on v4 data with the Round-3 ANOVA model
`biomarker ~ age + Ancestry + CDX_collapsed`. Output CSVs in
`AAIC2026_Repo/output_files/`.

### 6.1 Drop-one-country (TZA / NGA / UGA in turn)

Goal: are the within-AFDC sites driving the headline ancestry findings?

| | Count |
|---|---|
| Sig in FULL run | 98 |
| Robust in all 3 country drops | **93** (95%) |
| Lost in ≥1 drop | 5 |

The 5 borderline biomarkers — **TREM2, CXCL10, REST, MSLN, FCN2** — were all
weak FDR (0.0007–0.048 in FULL) and **none are hemolysis-flagged**. The 14
hemolysis-flagged biomarkers (HBA1 excluded) are all robustly significant in
every drop, confirming that the hemolysis pattern is consistent *across* the
three African countries — dropping one doesn't fix the problem.

**Conclusion:** Within-AFDC heterogeneity does not undermine the broad
ancestry findings. Country-driven effects are minimal.

Output: `Paper1_DropCountry_Sensitivity.csv`.

### 6.2 HBA1-as-covariate

Goal: can adjusting for HBA1 alone "correct" the hemolysis-confounded
biomarkers and still leave the underlying ancestry signal interpretable?

| | Count |
|---|---|
| Sig unadjusted | 97 |
| Sig HBA1-adjusted | **95** |
| Lost after adjustment | **2** (REST, BD_pTau_231) — neither hemolysis-flagged |

**All 14 hemolysis-flagged biomarkers stayed significant. Their F-stats
*increased* (more significant) after HBA1 adjustment — typically by 25–46%.**

This is a **negative result for HBA1-as-covariate as a fix.** Interpretation:

1. HBA1 captures only part of the hemolysis variance. Removing it concentrates
   the residual hemolysis signal in the other markers (PGK1, MDH1, SOD1...).
2. There's likely a coordinated multi-marker hemolysis pattern (or other
   site-handling effects) beyond what HBA1 alone reflects.
3. **HBA1-as-covariate is the wrong tool.** Better candidates:
   - **Site as a fixed effect** (captures country-of-collection directly)
   - **Drop hemolysis-flagged biomarkers from headline** results (option (a) from §5)
   - **Two-track reporting** (option (c) from §5)

Output: `Paper1_HBA1_Covariate_Sensitivity.csv`.

### 6.3 Updated recommendation
Going with **Site-as-covariate + drop-or-flag the 15 hemolysis biomarkers from
headline** is more honest than HBA1 adjustment. The non-hemolysis top hits
(TARDBP, NRGN, CD40LG, FGF2, CD63, CCL26, MME, IL15, TIMP3, ANXA5, CCL4,
CCL11, DDC) are stable across both sensitivities — these are the real ancestry
findings to feature in the manuscript.

---

## 7. Appendix: Full hemolysis assessment file

`output_files/biomarker_hemolysis_assessment.csv` — per-biomarker r(HBA1), tech CV, sample CV, excess CV, FDR, Hemolysis_Affected flag for all 118 biomarkers.

`output_files/hemo_data_with_outliers.csv` — sample-level Hemolysis Index with site-specific outlier flags. Useful for the §5 sample-level exclusion option.

---
