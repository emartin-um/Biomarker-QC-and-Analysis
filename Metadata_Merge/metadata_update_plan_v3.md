# Metadata Update Plan: n2200 → n2850
**Repo:** `Alamar_Biomarker_QC_Repo`  
**Date:** 2026-03

---

## Step 1: Protect current state (git)

Check for untracked files before doing anything:

```bash
git status
git stash list
```

If there are untracked files you want to keep, stash or commit them first:

```bash
# Option A: commit them
git add .
git commit -m "WIP: saving untracked files before metadata update"

# Option B: stash
git stash push -u -m "untracked before metadata-update-n2850"
```

Then tag and branch:

```bash
git checkout main
git pull
git tag v1.0-n2200-meta-2026Jan16
git push origin v1.0-n2200-meta-2026Jan16
git checkout -b metadata-update-n2850-OtherProjects-03052026
```

---

## Step 2: Merge metadata files (R)

**Input files:**
- Original: `Datasets/U19_Alamar_metadata_2026Jan16.csv`
- Update: `Datasets/OtherProjects_metadata_03052026.xlsx`

```r
library(tidyverse)
library(readxl)

meta_orig <- read_csv(
  "Datasets/U19_Alamar_metadata_2026Jan16.csv"
) |> mutate(metadata_source = "U19_2026Jan16")

meta_new <- read_xlsx(
  "Datasets/OtherProjects_metadata_03052026.xlsx"
) |> mutate(metadata_source = "OtherProjects_03052026")

meta_combined <- bind_rows(meta_orig, meta_new)

write_csv(
  meta_combined,
  "Datasets/U19_OtherProjects_metadata_combined_n2850_2026Mar.csv"
)
```

**Do NOT overwrite** `U19_Alamar_metadata_2026Jan16.csv` or `OtherProjects_metadata_03052026.xlsx`.  
Duplicate checking and SAMPLE_ALIQUOT validation handled by the Metadata_Merge pipeline.

---

## Step 3: Re-run Metadata_Merge pipeline

- Point pipeline input to:  
  `Datasets/U19_OtherProjects_metadata_combined_n2850_2026Mar.csv`
- Output to new Box directory:  
  `QC_output_n2850_2026Mar/`
- Original outputs stay untouched in:  
  `QC_output_n2200_2026Jan16/`

---

## Step 4: Validate

- Re-run n2200 subset and confirm outputs match prior run exactly
- Check sample counts at each pipeline step sum to ~2850
- Check missingness in the n650 new samples
- Confirm `metadata_source` column tracks correctly through outputs

---

## Step 5: Merge to main

```bash
git add .
git commit -m "Metadata update: OtherProjects_03052026 added, combined n=2850"
git checkout main
git merge metadata-update-n2850-OtherProjects-03052026
git tag v2.0-n2850-meta-2026Mar
git push origin main --tags
```

---

## Step 6: Downstream analyses

- Update downstream repos/scripts to point to `QC_output_n2850_2026Mar/`
- Note N change in analysis headers/READMEs
- Keep `QC_output_n2200_2026Jan16/` on Box — do not delete

---

## Files to preserve (do not overwrite)

| File | Notes |
|---|---|
| `Datasets/U19_Alamar_metadata_2026Jan16.csv` | original metadata |
| `Datasets/OtherProjects_metadata_03052026.xlsx` | new metadata, read-only input |
| `QC_output_n2200_2026Jan16/` | Box, prior QC outputs |
| git tag `v1.0-n2200-meta-2026Jan16` | repo snapshot before update |
