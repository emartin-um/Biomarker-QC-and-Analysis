# Triage_Metadata — do the QC decisions know who the specimen came from?

Primary QC removes wells (triage) and marks wells (flags) using assay statistics
only. It is agnostic to phenotype **by design**, and the clinical collection-site
label does not exist until the metadata merge. So nothing in the pipeline has
ever asked whether those decisions land evenly across the people the specimens
came from.

```bash
Rscript run_triage_metadata.R     # from this directory, inside a QC run
Rscript test_triage_metadata.R    # acceptance test (61 assertions)
```

Base R only. **Diagnostic only — drops nothing, flags nothing, changes no
pipeline output.**

Two questions, with different stakes:

* **Clinical collection site.** A site whose specimens fail QC more often is a
  collection, handling or shipping problem. Actionable, and a different owner
  from a plate problem or a well-position problem.
* **Diagnosis.** Differential triage by case/control status biases every
  downstream effect estimate. Unlike a site effect, nothing later in the
  pipeline would reveal it.

---

## What it found on the 50-plate frozen run

### 1. The "40% site" is not a finding. The real site signal is smaller and spread out.

`Triage_Review` reported triage rate by site reaching **40.0%** with permutation
p = 0.0040 on the maximum rate. That 40% is **2 wells out of 5** (site UJTH), and
its 95% Wilson interval is **[11.8%, 76.9%]** — an interval that spans from "less
than average" to "most of the plate".

This module puts a size floor at n = 20 and reports two statistics instead of
one. The result changes:

| statistic | value | reading |
|---|---|---|
| max site rate (n ≥ 20) | 5.32% (MCRH), **p = 0.40** | no single site stands out |
| omnibus deviance across 19 sites | 30.57, **p = 0.039** | a real but diffuse association |
| omnibus, permuted **within plate** | p = 0.063 | marginal once plates are held fixed |

So: **there is a modest site association, and it is spread across sites rather
than driven by one bad one.** The headline 40% was small-sample noise, and a
statistic that keys on the maximum will find that kind of noise every time.

Groups below the floor are still written to `triage_rates_by_metadata.csv` with
`above_floor = FALSE` — excluded from testing, never hidden.

### 2. The read flag's site and diagnosis associations are read depth, entirely

Tested alone, the read flag looks strongly associated with both site
(deviance 76.7, p = 0.0001) and diagnosis (20.8, p = 0.0067). Fitted together
with read depth:

| label | deviance alone | given the others | keeps |
|---|---:|---:|---:|
| Site | 76.7 | 10.3 | **13%** (p = 0.92) |
| Diagnosis | 20.8 | 4.7 | **23%** (p = 0.45) |
| Ancestry | 11.6 | 5.5 | 47% (p = 0.48) |
| `log2(reads)` | 566.6 | 491.6 | 87% (**p = 6e-109**) |

The read flag *is* a depth statistic, and sites differ in depth. Both apparent
associations are composition. **Do not treat the read flag as evidence that a
site has a specimen problem.**

### 3. The IC axis carries genuine, mutually independent associations

The IC gate supplies about three quarters of all triage. On that axis, three
labels survive adjustment for each other **and** for read depth:

| label | deviance alone | given the others | keeps | p |
|---|---:|---:|---:|---:|
| Site | 31.4 | 29.2 | 93% | **0.046** |
| Ancestry | 20.3 | 15.7 | 77% | **0.016** |
| age band | 10.1 | 8.9 | 89% | **0.030** |
| Diagnosis | 7.4 | 6.4 | 86% | 0.27 |

None of these is a proxy for the others. This is the substantive result of the
module, and it belongs in any analysis that stratifies on ancestry or adjusts
for age.

### 4. Triage is age-associated, and that is new

Triaged wells are **3.6 years older** on average (75.7 vs 72.1, standardised
difference **+0.35**). Age band survives adjustment on `any triage` (p = 0.015),
IC (p = 0.030) and burden (p = 0.044). It is not a site effect and not a depth
effect.

### 5. Diagnosis does *not* differentially drive triage — with a stated bound

| diagnosis | n | triage rate | 95% CI |
|---|---:|---:|---|
| Other dementia | 174 | 3.45% | [1.59, 7.32] |
| **Not codeable** | 432 | **3.24%** | [1.94, 5.37] |
| Case (AD) | 681 | 1.91% | [1.12, 3.24] |
| MCI/CIND | 980 | 1.73% | [1.09, 2.76] |
| Control | 1878 | 1.38% | [0.95, 2.02] |

Case vs Control rate ratio **1.38 (95% CI 0.71–2.67)**; diagnosis does not
survive adjustment on any axis (p = 0.18–0.50).

> **This is not proof of no effect.** At these group sizes and a 1.8% base rate
> the design could only have detected a **2.1×** elevation. A CI spanning 1 is
> not equivalence. Every null in `metadata_association_tests.csv` carries its own
> minimum detectable rate ratio for exactly this reason.

The interesting row is **Not codeable** — wells whose diagnosis is "Missing" or
"Insufficient Data" are triaged at 3.24%, more than double the Control rate, and
are over-represented among triaged wells (18.2% vs 10.1%, std diff +0.23). Poor
metadata capture and poor assay outcome co-occur. That is a finding about data
collection, not about dementia, which is why the collapse keeps administrative
non-answers as their own level instead of folding them into a clinical one.

### 6. Wells with no metadata are not neutral

113 of 4200 wells have no metadata row. They are triaged at **2.65%** against
1.81% for labelled wells — reported in `unlabelled_wells.csv` rather than
silently dropped from the denominator.

---

## Method notes

**Per axis, never pooled.** The triage axes are near-disjoint and point in
opposite directions on assay quality (`Triage_Review` §1). A metadata test on
"triaged" as a lump averages three different screens into one uninterpretable
number. IC, PCA, burden, the read flag, and their union are tested separately —
and §2 and §3 above are precisely why that matters: the read flag's site
association is fake and the IC axis's is real.

**Permutation, not chi-square.** 77 triaged wells over 23 sites leaves most cells
at 0 or 1, where the asymptotic test is not conservative but anti-conservative.

**Two statistics, reported side by side.** `max_rate` is sensitive to a single
bad group but hostage to small ones; the omnibus deviance pools evidence and
cannot be moved far by a group of 5. *Them disagreeing is itself the finding* —
max significant with omnibus null means one small group is carrying everything.

**Separability is checked before any conditioned claim.** Sites ship in batches,
so "this site fails more" and "this site's samples landed on worse plates" could
be the same data. `design_confounding.csv` reports it first: Site × plate
Cramér's V = 0.375, sites span a median of 15 plates, and a site's largest plate
holds a median of 19% of its wells. **The conditioned test is estimable here.**
Had site been nested in plate, the module would say the question is unanswerable
rather than report a number.

**Mutual adjustment, not model selection.** Site, ancestry, diagnosis and age are
correlated, so testing them one at a time lights several of them up off one cause.
The full model is fitted **once** and every term reported given the others — no
stepwise, no p-value shopping, so a term that drops out is visible rather than
absent.

> A `keeps` above 100% is not a bug. With correlated predictors a term can
> explain *more* once the others absorb variation it was competing with.

**Every null is bounded.** A non-significant result reports the smallest rate
ratio it had 80% power to detect. "No association" without that number is not a
result.

**Ancestry comes from the pipeline's own curated table**
(`Metadata_Merge/review/Ancestry_review.csv`), not from a second derivation
invented here — including its two easily-misread conventions (a blank condition
is a *wildcard*; rules apply in reverse so the first row wins). The first attempt
at re-deriving it from `Race`/`Ethnicity` returned "Unknown" for all 4087 wells,
because the codes are two-letter (`BL`/`WH`/`MU`/`AI`/`HP`, `HI`/`NH`) and not
the spelled-out words a regex expects. If the review table is absent the column
is left `NA` and the ancestry tests are skipped — better than a column of
"Unknown" that looks like an answer.

**The diagnosis collapse is written out.** 15 raw values → 6 levels, in
`diagnosis_collapse_map.csv`, so it can be argued with.

---

## Inputs

Preferred: `Primary_QC/output_files/well_level_QC.csv` — one row per well with
every gate flag, written by the 2026-07-30 Primary QC. This module is the reason
that file exists.

Older runs: the same spine is reconstructed from `Sample_QC.csv` +
`samples_to_triage.csv` + `flagged_read_outliers.csv`, so the module still runs
on them (without `mean_INT`).

Labels: `Metadata_Merge/<output_dir>/metadata_PrimaryQC_refreshed.csv` and
`Metadata_Merge/review/Ancestry_review.csv`.

## Outputs

| file | what it answers |
|---|---|
| `metadata_association_tests.csv` | every axis × variable: both permutation statistics, within-plate versions, conditioned deviances, MDE, and a plain-language reading |
| `mutual_adjustment.csv` | each label given the others and read depth — **which one is actually carrying it** |
| `triage_rates_by_metadata.csv` | per-group rates with Wilson CIs, including groups below the size floor, marked |
| `site_detail.csv` | one row per site, every axis side by side, for a site conversation |
| `design_confounding.csv` | is the conditioned question estimable at all |
| `balance_triaged_vs_retained.csv` | who got removed vs who stayed, as standardised differences |
| `balance_flagged_vs_unflagged.csv` | the same for read flags |
| `diagnosis_collapse_map.csv` | the 15 → 6 mapping |
| `unlabelled_wells.csv` | wells with no metadata row, and their triage rate |

---

## What this module does *not* claim

* **It does not establish that triage is unbiased with respect to diagnosis.** It
  bounds the effect (≤ 2.1× detectable) and reports the CI. That is weaker, and
  it is what the data support.
* **It does not explain the IC/site/ancestry/age associations.** It establishes
  that they are mutually independent and survive depth adjustment. Whether the
  cause is collection, handling, shipping, storage time, or something about the
  specimens themselves is a wet-lab and logistics question this module cannot
  reach.
* **It does not propose a drop rule, and no result here should become one.**
  Removing wells because of who they came from is how you manufacture the exact
  bias §5 is checking for.
* **It says nothing about plate position.** That is the *other* meaning of
  "site" — see `Primary_QC` §D.5.4 and `Triage_Review` §5.

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
