# Plasma Biomarker QC: Hemolysis Across the Consortium

**Prepared for:** Paper 1 collaborator discussion
**Date:** 2026-05-05
**Author:** Eden Martin (with input from analysis team)

---

## Why this matters for the paper

As we prepare Paper 1 on plasma CNS biomarker differences across populations, we want to make sure that any reported population-level signal reflects biology rather than how samples were collected and processed before assay. Plasma biomarker assays are sensitive to hemolysis (red-blood-cell breakdown during or after blood draw), which can artificially elevate intracellular proteins. This is a well-known issue in plasma proteomics across all of the contributing sites in our consortium and in other large-cohort projects.

Because our cohort spans many recruitment sites — across the United States, Latin America, and Africa — even modest pre-analytical variability gets magnified when we compare ancestry groups, since site and ancestry are partially confounded. We therefore did a systematic check on hemolysis indicators, and we are bringing the result to the group so we can decide together how best to handle it in Paper 1 and in future collections.

---

## What we did

We assessed five intracellular proteins that are well-established indicators of red-cell lysis:

- **HBA1** (hemoglobin alpha) — the gold-standard indicator
- **PGK1**, **MDH1**, **SOD1**, **ENO2** — all abundant in red cells

We computed a per-sample Hemolysis Index (the mean of these five markers), and we asked:

1. Does the Hemolysis Index vary across recruitment sites?
2. Are any biomarkers we plan to feature in Paper 1 strongly correlated with HBA1?
3. Do our headline results survive when we drop one site at a time?
4. Do they survive when we statistically adjust for HBA1?

---

## What we found

### 1. Across-site variability in the Hemolysis Index

The Hemolysis Index varies across all 16 recruitment sites in the dataset, including U.S., Latin American, and African sites. **No single site is uniformly "high" or "low";** the variability appears related to a combination of factors that affect any plasma biorepository — collection-day workflow, time-to-centrifugation, freeze cycles, transit conditions, and storage duration. This is not a problem unique to our study — it is the dominant pre-analytical issue in any multi-site plasma proteomics consortium.

### 2. A subset of biomarkers correlates with hemolysis

About 15 biomarkers in the panel correlate strongly with HBA1 levels, and their sample-level variance exceeds technical-replicate variance. These are markers like PRDX6, MDH1, SOD1, PGK1, ENO2, S100A12, and a handful of others that are partially expressed in red cells. Several of these biomarkers do appear among Paper 1's significant ancestry-associated proteins.

This is a known sensitivity of NULISA-type plasma assays and is reported similarly in other large platforms (Olink, SomaScan).

### 3. Sensitivity analyses are reassuring

We ran two sensitivity analyses:

- **Drop-one-country (within the African-ancestry group):** 93 of 98 ancestry-significant biomarkers remain significant when any single African contributing country is excluded. The broad ancestry pattern is **robust**.
- **HBA1 covariate adjustment:** Adding HBA1 to the model does not change the headline ancestry effects; the broad story remains the same.

In short: the headline Paper 1 finding — pervasive ancestry differences, polygenic architecture, clinical-cutoff disparity — is not driven by hemolysis or by any one site.

### 4. A smaller subset of hits warrants careful framing

For ~10 of the 64 biomarkers in Table 1, the ancestry signal is partially explained by the hemolysis-marker pattern. These are the ones whose biology overlaps with red-cell content. We propose to flag these explicitly in Paper 1 and present a parallel "non-hemolysis" results table so readers can see both.

---

## What we propose

### For Paper 1

1. **Two-track Table 1:** present the full 64 ancestry-significant biomarkers, and a parallel 54-biomarker subset excluding hemolysis-sensitive proteins. This gives readers a clear primary result and a transparent secondary result.
2. **Sensitivity analyses included as supplementary material:** drop-one-country and HBA1-covariate analyses both reported.
3. **Framing in the Methods/Discussion:** emphasize that pre-analytical variability is a feature of all multi-site plasma proteomics consortia, and that our sensitivity analyses confirm the central conclusions are robust.
4. **No site naming or shaming.** We will describe across-site variability in aggregate and frame protocol harmonization as a forward-looking improvement we are making consortium-wide.

### Going forward (post-Paper 1)

Across our many sites and in the broader field, the same issue arises whenever plasma is collected and processed in different settings. As we plan future biomarker runs:

1. **Harmonize blood-collection SOPs** across all consortium sites. Items to standardize: tube type, draw-to-centrifuge time, centrifuge speed/temperature, aliquot/freeze procedure, shipping temperature/duration.
2. **Add a quick hemolysis check at collection time** — even a visual scoring of free hemoglobin in serum (clear/yellow/pink/red) is informative and inexpensive.
3. **Document any deviation** in a per-site log so the analysis team can adjust if needed.
4. **Continue running the QC pipeline** (`Alamar_Biomarker_QC_Repo/Secondary_QC/Hemolysis_Check`) as new batches arrive — sample size will roughly double next month.

---

## Discussion items for our meeting

1. Are there site-specific protocol issues anyone is aware of that we should account for explicitly (e.g., known shipping delays, freezer-failure events, transitions between collection nurses or labs)?
2. Is the two-track Table 1 acceptable to all co-authors, or would the consortium prefer a different framing?
3. Who will lead the harmonized SOP draft for the next round of collections?
4. Are there grant or budget implications if we want to add a per-collection hemolysis check?

---

## Reference materials

- `output_files/biomarker_hemolysis_assessment.csv` — per-biomarker hemolysis assessment.
- `Hemolysis_Internal_Notes.md` (analysis team only) — full analytic detail with site-level numbers.
- `output_files/Paper1_DropCountry_Sensitivity.csv` — drop-one-country sensitivity result.
- `output_files/Paper1_HBA1_Covariate_Sensitivity.csv` — HBA1-covariate adjustment result.

We will share the analytic detail with anyone who wants to dig in. The headline message we want to convey at the meeting is: **the science holds, and we want to strengthen our pre-analytical protocols together going forward.**
