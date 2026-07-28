# Cohort_Accounting — per-well exit reasons and the standing cohort-loss accounting

Implements `QC_PIPELINE_RECOMMENDATIONS_2026-07-13.md` **§1** (separate QC-triaged
from cohort-filtered everywhere) and **§4** (a standing cohort-loss accounting) —
the last two open items from the 2026 batch investigation.

**Nothing here is a drop rule.** The module removes nothing and changes no
pipeline output. It says, for each well that is *already* absent, **which stage
removed it**.

```bash
Rscript run_cohort_accounting.R     # from this directory, inside a QC run
Rscript test_cohort_accounting.R    # acceptance test
```

Base R only — no dplyr — so it runs anywhere the QC runs do.

---

## Why

A well missing from the analysis set was repeatedly read as a well that "failed
QC". They are different things with different owners, and they are not even the
same order of magnitude. On the 50-plate frozen run:

| | wells | share of the loss |
|---|---:|---:|
| **QC triage** (Primary QC) | 77 | **11.9%** |
| **Cohort construction** (Metadata_Merge) | 571 | **88.1%** |
| lost overall, of 4200 attempted | 648 | 15.4% |

A bay that is "40% missing" from missing diagnoses is a data-entry issue; one
that is "40% missing" from IC failures is a wet-lab issue. **Same number,
opposite owner.** Publishing a single merged "% removed" per bay is what set off
the entire 2026 batch investigation.

The clearest single case on this run is `20250804-1712 Bay2`: **3.6% QC loss and
52.4% cohort loss.** Merged, it reads as a 56%-dead bay and goes on a re-run
list. Split, it is a plate that ran fine and whose patients mostly lack a
codeable diagnosis.

---

## The chain, verified end to end

Every number below was re-derived on 2026-07-27 from the 50-plate frozen run and
is re-checked by `test_cohort_accounting.R` on every run.

```
4800   rows in Primary_QC/input_files/Sample_QC.csv
-500   NC / SC / IPC control wells      (Sample Type != "Sample")
=4300  PLASMA wells
-100   HIHG pool replicates             -> pool_replicate
=4200  ATTEMPTED PATIENT WELLS  <- the denominator
- 77   Primary QC triage                -> qc_triage_IC (58) / _PCA (10) / _outlier (9)
=4123  == NPQ_*_post_QC.csv
-110   no metadata row                  -> not_in_metadata
=4013  == merged_combined_post_QC.csv
- 57   specimen assayed twice           -> duplicate_aliquot_rerun
-  5   age/BMI out of range             -> covariate_out_of_range
=3951
-361   diagnosis missing/not codeable   -> missing_diagnosis
=3590
- 19   no ancestry rule matches         -> missing_ancestry
=3571  == filtered_standard_post_QC.csv
- 19   second visit collapsed away      -> repeat_visit_collapsed
=3552  == filtered_combined_post_QC.csv   IN COHORT
```

The three `==` lines are independent confirmations, not assumptions: the module
re-derives the cohort from `merged_combined_post_QC.csv` and lands on the same
3552 `SAMPLE_ALIQUOT`s the pipeline did, with the same `CDX_collapsed` and the
same `Ancestry`.

---

## The taxonomy is read off the code, not off the recommendation

The 2026-07-13 list was written from the outside. This one was read from what
`Metadata_Merge_Pipeline.Rmd` and `covariate_explorer.R` actually do, and three
things differ. All three matter.

**`qc_read_outlier` does not exist as an exit.** Read outliers are **flagged and
kept** — 94 of them, of which 65 are in the final cohort and the rest leave for
unrelated reasons. Carrying it as an exit reason would have created a
permanently-empty category *and* implied that low-read wells are dropped, which
they are not. It is here as the annotation **`flag_read_outlier`**.

**`qc_triage_baddata` had to be added.** Primary QC's triage table has a fourth
reason, `"Bad Data"` (non-numeric rows), that the recommendation did not list. It
was 0 on this run — which is exactly why it was easy to miss. The runner
**stops** if it meets a triage reason it does not recognise rather than letting
it fall through into `in_cohort`.

**`missing_covariate` splits in two,** because the pipeline has no single
covariate filter. It drops on an unmatched race/ethnicity combination
(`missing_ancestry`) and on an out-of-range age or BMI
(`covariate_out_of_range`). It **never drops on a missing covariate**: rows with
a missing age, sex, APOE or education reach the final dataset and are only
reported in `missing_data/`.

And `repeat_visit_collapsed` is split from `duplicate_aliquot_rerun` — the
pipeline collapses a re-assayed specimen (most-recent `RUN`) at a different step
from a genuine second visit (highest age per `Record_ID`), and they have
different owners. That is the recommendation's own §5 (*flag every axis
independently*) applied to its own §1 list.

### Full taxonomy

| exit_reason | stage | class | owner |
|---|---|---|---|
| `pool_replicate` | Primary_QC | by_design | HIHG pool well, never a patient |
| `qc_triage_baddata` | Primary_QC | qc_triage | wet-lab — non-numeric assay output |
| `qc_triage_IC` | Primary_QC | qc_triage | wet-lab — internal-control failure |
| `qc_triage_PCA` | Primary_QC | qc_triage | wet-lab — multivariate profile outlier |
| `qc_triage_outlier` | Primary_QC | qc_triage | wet-lab — per-marker outlier burden |
| `not_in_metadata` | Metadata_Merge | cohort_construction | data-entry — no metadata row |
| `metadata_review` | Metadata_Merge | cohort_construction | data-entry — manual exclusion |
| `duplicate_aliquot_rerun` | Metadata_Merge | cohort_construction | by design — specimen assayed twice |
| `covariate_out_of_range` | Metadata_Merge | cohort_construction | data-entry — age/BMI implausible |
| `missing_diagnosis` | Metadata_Merge | cohort_construction | data-entry — diagnosis not codeable |
| `missing_ancestry` | Metadata_Merge | cohort_construction | data-entry — no ancestry rule matches |
| `repeat_visit_collapsed` | Metadata_Merge | cohort_construction | by design — collapsed to one visit |
| `in_cohort` | — | retained | — |

**Order is load-bearing.** A well that would qualify for several exits is
attributed to the **first** one, because that is the stage that actually removed
it — and the overlap stays visible rather than being absorbed. Two columns do
that. `n_stages_qualified` counts how many exits each well qualified for; on this
run it is **0 wells above 1**, i.e. the stages happen to be cleanly nested, which
is worth knowing precisely because it cannot be assumed. The four
`flag_qc_triage_*` booleans keep every triage axis a well tripped, and there **16
of the 77 triaged wells trip more than one** — so the single reason genuinely is
a choice, and the full set is still on the row.

---

## Outputs

| file | rows | what it is |
|---|---:|---|
| `well_exit_reasons.csv` | 4300 | one row per attempted well, exactly one reason, plus the independent flags and the cut labels |
| `cohort_loss_summary.csv` | 13 | attempted → in-cohort by reason, with the owner attached |
| `cohort_loss_by_year.csv` | 2 | the §4 differential cuts — QC and cohort as **separate** columns |
| `cohort_loss_by_site.csv` | 23 | |
| `cohort_loss_by_ancestry.csv` | 7 | |
| `per_plate_QC_and_cohort_summary.csv` | 50 | the §1 side-by-side table |

### `per_plate_QC_and_cohort_summary.csv` supersedes reading `per_plate_QC_summary.csv` alone

Primary QC's `per_plate_QC_summary.csv` is not wrong — it is *upstream*. It
counts triage and flags and stops there, because at that point in the pipeline
the cohort filter has not run. Reading it as "what happened to this bay" is what
silently attributes the cohort loss to the wet lab. This file carries the same QC
columns (`qc_triaged` is asserted equal to its `Triaged`) with the cohort columns
beside them, and **never sums the two**. The test asserts no merged
`pct_removed` / `total_removed` column exists.

### The §4 cuts, on this run

**By year — the differential the recommendation predicted, now standing output:**

| year | attempted | in cohort | QC triage | cohort construction |
|---|---:|---:|---:|---:|
| 2025 | 2940 | 2606 | 1.6% | **9.7%** |
| 2026 | 1260 | 946 | 2.3% | **22.6%** |

2026 loses about 2.3× the share 2025 does, and essentially all of the difference
is cohort construction, not QC. That is the pattern that quietly biases a
downstream analysis if nobody is looking.

**By site**, the highest cohort-construction losses are the large US sites —
IHG 17.8%, NC 16.9%, CUMC 15.1%, CU 14.1% — driven almost entirely by
`missing_diagnosis`, while the African teaching hospitals sit near 0%. This is
the opposite of the intuition the depth and detectability axes create, and it is
worth stating plainly: on this run the sites that lose the most wells lose them
to **data entry**, not to the assay.

---

## Limits, stated rather than discovered later

**A cut can only see wells that carry its label.** Site and ancestry come from
the metadata, so the 110 `not_in_metadata` wells are invisible to the by-site and
by-ancestry cuts *by construction* — a well with no metadata row has no site. The
runner prints how many wells each cut could not label and which reasons they
carry, rather than dropping them quietly. On this run: 113 wells (110
`not_in_metadata` + 3 QC-triaged with no metadata row).

**`Ancestry == "Other"` in the by-ancestry cut is not a population.** It is the
set of wells whose Group × Race × Ethnicity combination matches no rule — the
same condition that produces `missing_ancestry`. It reads as 96.9% lost, which is
close to tautological. It is kept because it *names* the 32 wells whose
demographics the study has never classified, which is the useful part.

**The cuts label QC-triaged wells too.** Labels are taken from the metadata, not
from the surviving cohort, so a well that exited at Primary QC still carries its
Site (74 of 77 do). Without that the cuts would silently become cohort-only. The
test asserts the coverage stays above 80%.

**Two upstream behaviours are inherited deliberately,** because the accounting
has to match what the pipeline does rather than what it ought to do:

- The visit collapse ranks `age_at_subject` with `rank()`'s defaults, so a visit
  with **NA age ranks last and is therefore the row kept**, and tied ages keep
  the later row in data order. If any rows reach the collapse with a **missing
  `Record_ID`**, they pool into one group and only one survives — the runner
  counts those and warns.
- A `CDX` value absent from `CDX_review.csv` entirely is "unaccounted" and gets
  dropped by the same test as a deliberately excluded one. The pipeline warns,
  but the chunk sets `warning: false`, so in the rendered report a brand-new CDX
  spelling drops its samples with **no visible signal at all**. The runner names
  any such value in its console summary.

**Requires the review tables.** `review/CDX_review.csv` and
`review/Ancestry_review.csv` are the pipeline's data-driven config and this
module reads the same files. If they are absent the pipeline falls back to a
built-in config inside `Metadata_Merge_Pipeline.Rmd` Step 6b, which this module
cannot read — so it stops with that explanation rather than guessing.

---

## The test

`test_cohort_accounting.R` — 59 checks on a full run, base R, non-zero exit on
failure.

`output_files/` is gitignored, so checking out another branch does **not** revert
the CSVs. Branch-switching cannot tell you whether a change was safe; only a
content check can.

**Part A** unit-checks the helpers against hand-built data where the answer is
known by inspection: first-match-wins and overlap counting in `classify_exits`,
triage-reason precedence, the three `rank()` edge cases in the visit collapse,
and the ancestry rule engine's reverse-order precedence and NA-component
widening. Part A runs in the source repo, where there is no data.

**Part B** is the load-bearing half: re-deriving the cohort from
`merged_combined_post_QC.csv` must reproduce `filtered_combined_post_QC.csv`
**exactly** — same 3552 ids, and the 3571 intermediate must equal
`filtered_standard_post_QC.csv`. If Metadata_Merge changes a drop rule and this
module does not, that equality breaks. A silently drifted accounting is worse
than none.

It caught one real defect on its first run: the CDX-map test used `[[`, which
throws on a missing name, where the runner correctly uses `[`, which yields the
`NA` that makes an unaccounted diagnosis a `missing_diagnosis`.

---

## Where this lives

Canonical: `Alamar_Biomarker_QC_Repo/Metadata_Merge/Cohort_Accounting/`.

It sits under `Metadata_Merge` rather than `Secondary_QC` because that is where
the cohort filter lives and where the drop rules it mirrors are defined —
`scaffold_qc_run.R` copies `Metadata_Merge` by default, so the module ships with
every new run without anyone naming it.

`scaffold_qc_run.R` copies code from the canonical repo into each run directory,
so **editing the copy under `QC_Runs/` is lost on the next scaffold.** Edit here,
then sync the run-directory copy.

Run it **after** `Metadata_Merge_Pipeline.Rmd`, since it reads that pipeline's
outputs.
