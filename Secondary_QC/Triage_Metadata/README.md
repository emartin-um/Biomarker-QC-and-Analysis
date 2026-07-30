# Triage_Metadata — do the QC decisions know who the specimen came from?

Primary QC triages and flags wells on assay statistics alone. It is agnostic to
phenotype **by design**, and the clinical collection-site label does not exist
until the metadata merge. So nothing in the pipeline asks whether those decisions
land evenly across the people the specimens came from.

This module asks. **Diagnostic only — it drops nothing, flags nothing, and
changes no pipeline output.**

```bash
Rscript run_triage_metadata.R     # from this directory, inside a QC run
Rscript test_triage_metadata.R    # acceptance test
```

Base R only, so it runs anywhere the QC runs do.

**Results for a given run go to
`output_files/triage_metadata_report_<run>.html`**, not into this file. A README
that reports findings is stale the moment it meets a different dataset.

---

## The two questions, and why they have different stakes

**Clinical collection site.** A site whose specimens fail QC more often has a
collection, handling or shipping problem. That is actionable, and it has a
different owner from a plate problem or a well-position problem.

**Diagnosis.** If triage removes cases at a different rate from controls, every
downstream effect estimate is biased. Unlike a site effect, nothing later in the
pipeline would reveal it.

Sex, ancestry and age are tested alongside them, both in their own right and
because they are the obvious confounders of the first two.

---

## How it works

### Per axis, never pooled

The triage axes are near-disjoint and point in opposite directions on assay
quality (see `Triage_Review`). Testing "triaged" as a lump averages several
different screens into one uninterpretable number. IC, PCA, outlier burden, the
read flag and their union are each tested separately.

This is not a stylistic preference — the axes genuinely disagree about metadata,
and a pooled test hides which screen is responsible.

### Two statistics, reported side by side

| statistic | strength | weakness |
|---|---|---|
| **maximum group rate** | reacts to a single unusual group | hostage to small groups — a group of 5 can top the table on noise |
| **omnibus deviance** | pools evidence across all groups | a genuinely bad single group gets diluted |

Both are permuted. **Them disagreeing is itself a finding:** max-rate significant
with the omnibus null means one small group is carrying the entire result, and
the report says so rather than reporting the smaller p-value.

### Permutation, not chi-square

Triage counts spread over ~20 sites leave most cells at 0 or 1, where the
asymptotic test is not conservative but anti-conservative. Tests are also run in
a within-plate-stratified form, which holds each plate's flag count fixed so a
bad plate cannot manufacture a group signal.

### A size floor that reports what it excluded

Groups below `MIN_N` are not tested — but they are written to
`triage_rates_by_metadata.csv` with `above_floor = FALSE`, and they appear in the
report. "The signal was in the groups we didn't show you" is a failure mode this
module exists to prevent. Every rate carries a **Wilson** interval, which stays
inside [0,1] and stays honestly wide at small n, unlike a Wald interval.

### Separability checked before any conditioned claim

Sites ship in batches, so "this site fails more" and "this site's samples landed
on worse plates" can be the same data. `design_confounding.csv` reports Cramér's
V *and* the spread that actually decides estimability — how many plates each site
spans, and what share of a site sits on its single largest plate. If a label is
effectively nested inside plate, the module reports the question as
**unanswerable** rather than returning a number.

### Mutual adjustment, not model selection

Site, ancestry, diagnosis and age are correlated: a site recruits a particular
population, whose age and diagnosis composition follow. Testing them one at a
time lights several of them up off a single underlying cause.

The full model — every label plus `log2(reads)` — is fitted **once**, and each
term is reported given the others. No stepwise, no p-value shopping, so a term
that drops out is visible rather than absent.

> A `retained_share` above 1 is not an error. With correlated terms a variable
> can explain *more* once the others absorb variation it was competing with.

### Every null is bounded

A non-significant result reports the smallest rate ratio it had 80% power to
detect. "No association" without that number is not a result — at a ~2% base rate
and small groups, a design can easily be unable to see anything short of a
several-fold elevation, and a bare null reads as reassurance it has not earned.

### Categorical labels come from the pipeline's own curated tables

* **Diagnosis** from `Metadata_Merge/review/CDX_review.csv`
* **Ancestry** from `Metadata_Merge/review/Ancestry_review.csv`

Re-deriving either here would put a second, drifting definition in the same repo,
and a hand-written matcher only knows the vocabulary of the cohort it was written
against. A fallback matcher exists for runs without the review tables, and
labels itself as such in `diagnosis_collapse_map.csv`.

> **One deliberate difference from the merge.** The review tables mark some
> levels for *exclusion from the analysis cohort*. This module keeps those wells
> as an explicit "Not codeable" level. A well excluded from the cohort was still
> either triaged or not, and how those wells behave is a finding about data
> capture.

---

## Settings

All at the top of `run_triage_metadata.R`; nothing below that block is specific
to any one study.

| setting | default | what it controls |
|---|---|---|
| `MIN_N` | 20 | group size floor for testing (below-floor groups are still reported) |
| `N_PERM` | 10000 | permutation replicates |
| `AGE_BREAKS` | 60, 70, 80 | age-band cutpoints — widen for a younger cohort |
| `GROUP_VARS` | Site, Diagnosis, sex, Ancestry, age_band | which metadata columns to test |
| `BALANCE_CONTINUOUS` | age, BMI, l2reads, mean_INT_z | covariates in the balance table |

A variable that is absent from the merge, constant, or entirely missing is
skipped with a note rather than erroring, so a study whose metadata differs still
runs.

## Inputs

**Preferred:** `Primary_QC/output_files/well_level_QC.csv` — one row per well
with every gate flag. This module is the reason that file exists.

**Older runs:** the same spine is reconstructed from `Sample_QC.csv` +
`samples_to_triage.csv` + `flagged_read_outliers.csv`, so the module still runs
on runs made before that file existed (without `mean_INT`).

**Labels:** `Metadata_Merge/<output_dir>/metadata_PrimaryQC_refreshed.csv`, plus
the two `review/` tables above.

## Outputs

| file | what it answers |
|---|---|
| `triage_metadata_report_<run>.html` | **the run's findings**, as a standalone report |
| `metadata_association_tests.csv` | per axis × variable: both permutation statistics, within-plate versions, conditioned deviances, MDE, and a plain-language reading |
| `mutual_adjustment.csv` | each label given the others and read depth — which is carrying an association, which was a proxy |
| `triage_rates_by_metadata.csv` | per-group rates with Wilson CIs, below-floor groups marked |
| `site_detail.csv` | one row per site, every axis side by side |
| `design_confounding.csv` | whether the conditioned question is estimable at all |
| `balance_triaged_vs_retained.csv` | who was removed vs who stayed, as standardised differences |
| `balance_flagged_vs_unflagged.csv` | the same for read flags |
| `diagnosis_collapse_map.csv` | the raw → collapsed mapping, and which table produced it |
| `unlabelled_wells.csv` | wells with no metadata row, and their triage rate |

Standardised differences rather than p-values in the balance tables: against
thousands of retained wells a p-value mostly measures the sample size, while the
standardised difference measures the imbalance itself.

---

## What this module does not claim

* **It cannot establish that triage is unbiased with respect to diagnosis.** It
  bounds the effect and reports the interval. That is weaker, and it is what the
  data support.
* **It does not explain any association it finds.** Whether a cause is
  collection, handling, shipping, storage time, or the specimens themselves is a
  wet-lab and logistics question this module cannot reach.
* **It proposes no drop rule, and no result here should become one.** Removing
  wells because of who they came from manufactures exactly the bias the
  diagnosis section checks for.
* **It says nothing about plate position** — the other meaning of "site". That is
  Primary QC §D.5.4 and `Triage_Review` §5.

## Where this lives

Canonical: `Alamar_Biomarker_QC_Repo/Secondary_QC/Triage_Metadata/`.

`scaffold_qc_run.R` copies code from the canonical repo into each run directory,
so **editing the copy under `QC_Runs/` is lost on the next scaffold.** Edit here,
then sync. `Secondary_QC` is not in the scaffold's default module list — pass it
explicitly:

```bash
Rscript scaffold_qc_run.R "QC_myrun" "Primary_QC,Metadata_Merge,Secondary_QC"
```

Run after Primary QC and Metadata_Merge, since it reads both.
