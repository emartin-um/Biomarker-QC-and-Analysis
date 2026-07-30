# Triage_Review — auditing the Primary QC screens themselves

Primary QC removes wells (triage) and marks wells (flags). Every other module
asks *what happened to the samples*. This one asks **whether the screens are any
good**: does each gate fire at all, is it independent of the others, where does
its cut sit on its own statistic, and is it keying on plate position rather than
on the specimen.

```bash
Rscript run_triage_review.R      # from this directory, inside a QC run
Rscript test_triage_review.R     # acceptance test
```

Base R only. **Diagnostic only — drops nothing, proposes no new drop rule.**

---

## The five answers, on the 50-plate frozen run

### 1. What do the triaged and flagged wells have in common? Almost nothing.

They are near-disjoint (Jaccard 0.030, 5 wells shared) and point in **opposite**
directions on the one axis that separates them:

| class | n | mean_INT_z | Cohen's d | AUC |
|---|---:|---:|---:|---:|
| triaged (all) | 77 | −1.04 | −1.07 | 0.408 |
| — IC only | 50 | **−0.21** | −0.21 | 0.462 |
| — PCA (any) | 16 | **−4.56** | −4.77 | 0.116 |
| — burden (any) | 25 | −2.92 | −3.01 | 0.262 |
| read-flagged HIGH | 93 | **+1.37** | +1.43 | 0.876 |

**Do not read the pooled row.** "Triaged wells are dim (−1.04)" is true of the
average and false of most of them: the deficit lives almost entirely in the 16
PCA wells, while the 50 IC-only wells — 65% of all triage — sit at −0.21 with
AUC 0.462, i.e. indistinguishable from retained wells on intensity. The module
reports per class for exactly this reason, and the test asserts the classes
differ (IC-only d −0.21 vs PCA d −4.77) so that pooling can never quietly return.

### 2. Do they relate to the problems already found? Yes, and one is new.

The PCA and burden axes are catching the **dim** tail of the per-sample intensity
axis (addendum §9); the read flag is catching the **bright** tail. So Primary QC
already touches the axis the Extremes screen was blind to — from both ends, by
accident, through statistics designed for something else.

What that does to the analysis set is smaller than it looks. Across triage the
moment skewness of `mean_INT_z` goes −0.503 → +0.069 and the range −7.16…+3.79 →
−4.23…+3.22, but **robust skewness does not move at all** (Bowley +0.020 →
+0.020). QC removes about a dozen catastrophically dim wells; it does not reshape
the body of the distribution. The nuisance axis is still there and still belongs
in the model.

The new thing is positional — see §5.

### 3. Should thresholds be tightened? Not the burden gate. It is not a dial.

`FDR` is a deterministic function of an integer count over 18 anchor markers, so
**"FDR < 0.01" *is* "5 or more anchors out"**. The statistic takes 14 distinct
values in 4200 wells and **91.5% sit at exactly 1.0**. Around the cut:

| FDR | ≈ anchors out | wells here | removed at this cut |
|---|---:|---:|---:|
| 0.007357 | 5 | 3 | **25 ← current gate** |
| 0.011938 | 4 | 1 | 26 |
| 0.054750 | 4 | 13 | 39 |

There is nothing between 0.0074 and 0.0119 to tune. "Tightening" means moving to
4-of-18, which admits 14 more wells in one step.

And it would not fix the blind spot, because **the burden count is a distance
statistic, not a level statistic**: it counts markers outside a per-plate IQR
fence, in either direction, so a sample that is moderately high on all 131
markers never makes any single marker extreme. That is structural. No FDR value
recovers it — the verified sweep found that loosening all the way to FDR < 0.20
raises bright-well coverage only from 10/94 to 14/94 while costing 16 cohort
wells, and FDR < 0.75 costs 208 cohort wells and still does not get there.

**The instrument that would work already exists and is not being used as one.**
`mean_INT` is a level statistic and separates the extremes the burden test
cannot. The recommendation stands as in addendum §9: report it, model it, do not
gate on it.

> **Loosening PCA is not the alternative, and this is worth stating because it
> looks like one.** PCA round 1's `PC_z` does separate extreme-intensity wells
> better than the burden count (AUC 0.933 vs 0.676) — but that gap is **131
> markers beating 18 anchors, not PCA beating counting**: applying the same
> counting logic to the full panel gives AUC 0.913, indistinguishable from
> PC_z's 0.933 (bootstrap difference 0.020, 95% CI −0.004 to +0.046). And
> loosening the 5 SD cut cannot reach the bright tail in any case: **the maximum
> `PC_z` among all 94 bright wells is 4.226**, below the cut, so lowering it
> admits dim wells before it admits a single bright one. The sensitive
> instrument is not merely turned down — it is pointed elsewhere.

Two gates should be **fixed rather than tuned**, because they cannot fire:

- **PCA round 2 flagged 0 of 4184 wells** (max `PC_z` 4.824 against a 5.0 cut).
  The second pass adds nothing on this run.
- **The lower read bound is never plate-relative.** `median − 4×IQR` is negative
  in **50 of 50** plate-bays, so `max(500000, …)` always resolves to the absolute
  floor. A chronically low-input site cannot be flagged however far below its
  plate it sits, unless it drops under half a million reads — which happened to
  exactly 1 well. "Read outlier" sounds symmetric; 93 of 94 flags are *high*.

### 4. Are the flagged wells fine? Yes, keep them — but "redundant" is not proven.

93 of 94 are high-read (`l2reads` d = 2.50, AUC 0.982). They look bright
(`mean_INT_z` d = 1.43) and hemolysed (d = 1.86), but their cystatin C is dead
normal (d = 0.041), so this is not renal-driven brightness. 65 are in the
analysis cohort.

On "is it just an indicator to consider depth" the answer is **probably, but the
data cannot settle it**. Adding `read_flagged` to a model already carrying a
flexible spline in `l2reads` plus plate-bay moves R² for `mean_INT_z` by 0.0009 —
which looks like nothing and is not. At 2.16% prevalence the *conditional
ceiling* on incremental R² is about 0.0099, so 0.0009 is the arithmetic signature
of a **≈0.3 SD** shift in 2% of wells, not of a null. The minimum detectable
effect here is ≈0.43 SD and the 95% CI on the residual shift is **[+0.05, +0.49]**
— a half-SD residual brightness effect is *not* excluded.

So: keep the 65 in-cohort flagged wells, model `log2(reads)` as already planned,
and do not claim the flag has been shown to be redundant. It has been shown to be
small. Deciding it is redundant needs either more flagged wells or a targeted
test, not this one.

Two further caveats. Because the bound is *plate-relative*, "high-read" is a
property of the well's neighbours as much as of the well — the same specimen on
another plate need not be flagged.

And depth does not act on NPQ as a uniform multiplier. NPQ is IC-normalised and
`Reads` is empirically decoupled from `IC_Reads` (r = −0.015), so a pure
sequencing-depth effect should **divide out and give slope 0**. It does not: six
of eight markers reject slope = 0 at Bonferroni, and the slopes are strongly
heterogeneous — PGK1 1.590, MDH1 0.839, SOD1 0.689, HBA1 0.509, ENO2 0.376
against CST3 0.019, NEFL 0.019, GFAP −0.014 (**34-fold in SD units**, robust to
leave-one-site-out, leave-one-run-out, within-plate-bay and cohort-only). So
`log2(reads)` is a real and marker-specific term, not a formality — and one
common covariate will not absorb it identically across the panel.

### 5. Do they cluster by site? Weakly, and not where the first pass said.

> ⚠️ **Corrected 2026-07-30.** An earlier version of this section reported
> "triage by site reaches 40.0%, permutation p = 0.0040". That was **2 wells out
> of 5** at one site (Wilson 95% CI 11.8%–76.9%), found by a maximum-rate
> statistic with a size floor of n ≥ 5 — and a maximum keys on exactly that kind
> of extreme. `rate_by_group()` now defaults to **n ≥ 20** and reports an
> **omnibus** statistic beside the maximum, so the maximum can no longer be read
> on its own. See also `Secondary_QC/Triage_Metadata/`, which separates site from
> depth, ancestry and age.

Over the 19 sites with n ≥ 20 (4 smaller sites excluded, and named in the output):

| | omnibus deviance | omnibus p | max rate | max-rate p | after depth |
|---|---:|---:|---:|---:|---:|
| triage | 30.6 | **0.039** | 5.3% | 0.40 | 34.5 (**113%** retained) |
| read flags | 78.1 | **0.0001** | 15.0% | 0.003 | 20.7 (**26%** retained) |

The two behave oppositely once depth is partialled out, and that is the useful part:

- **Triage** is a weak, diffuse association across sites rather than one bad site
  — and it is **not** depth composition. Conditioning on read depth leaves it
  slightly larger. The top site is now MCRH at 5.3% (5 of 94), with no single
  site standing out (max-rate p = 0.40).
- **Read flags** look strongly site-structured, but **26% of that survives depth**
  — so it is mostly which sites submit deep samples, not a site QC problem. The
  top site is MULH at 15.0% (12 of 80), which is the deepest site in the study
  (~14M reads). Since the flag is a depth threshold, that is the expected
  direction. The old maximum-rate test missed this association entirely
  (p = 0.165).

The statistic mattered in both directions: it manufactured a triage signal that
is not there and hid a read-flag signal that is. The positional finding below is
unaffected — it rests on 50 observations of one well position, not on a small site.

> **Well F7 is triaged in 9 of its 50 plate-bays — 18.0%, against 1.64% at every
> other position.** Permutation on the maximum per-position count: observed 9,
> null mean 3.87, **p = 0.0001**. All 9 exit through the IC axis.

**The robust finding is a standing spatial gradient, not one broken well.** Median
`IC_Median` by well position shows a low corner — rows F–H × columns 9–12 — that
is present *identically in both reagent lots*:

| | corner (F–H × 9–12) | rest of plate | difference | perm p |
|---|---:|---:|---:|---:|
| 2025 | −0.124 | +0.021 | **−0.145** | 0.0002 |
| 2026 | −0.130 | +0.030 | **−0.161** | 0.0002 |

The whole position map correlates **r = 0.932** between the two lots, and the
permutation shuffles positions *within each plate*, so plate-level effects cannot
produce it. Raw IC reads also fall from row A to row H on **50 of 50** plate-bays
(median Spearman −0.75) — a consistent tendency, though **0 of 50** bays are
strictly monotone.

**F7 is the extreme of that gradient**, in the bottom decile of its own plate on
**15 of 15** 2026 plate-bays (16 of 35 in 2025), across 7 runs and 3 bays. In 2026
its whole population lands on the ±0.40 line: all 15 wells span −0.34 to −0.47, so
6 cross the flag and the other 9 sit within 0.06 of it. The gate splits one
uniformly shifted population near-arbitrarily.

> **What this module does *not* claim.** An earlier draft called F7 "2026 drift at
> one coordinate". That does not survive: across 84 positions the 2025→2026 change
> has SD 0.20 and F7's −0.38 is 1.8 SD of it, with two other positions (A8 −0.31,
> A1 −0.30) moving comparably. With 84 positions, a couple always will. The runner
> now prints the change against that spread instead of asserting drift.

It is also **not** control adjacency — the control block is columns 6 and 8 in rows
D–H, so D7/E7/G7/H7 are flanked identically and sit at +0.019 — and **not** sample
composition: F7 draws from 12 of 23 sites, and the corner shows no significant
association with site, diagnosis or ancestry (p = 0.078 / 0.94 / 0.44).

**It reaches the reported values, not just the flag.** Because NPQ is normalised to
the IC, corner wells read **+0.116 vs −0.019** on whole-panel level (p = 0.002).
Small — 0.13 SD — but systematic, and keyed to position rather than specimen.

This is a wet-lab and vendor question rather than an analysis one, and it is the
single most actionable thing here.

### 6. And one that changes how another module should be read

**Three quarters of triage is a vendor metric recomputed.** The in-house "IC
outlier" test flags `IC Reads` outside 0.6–1.4× the plate median; Alamar's own
`IC Median` column is the same check on the same quantity. Agreement is **57 of
58, Jaccard 0.966**, with one discordant well each way sitting within 0.01 of the
0.40 bar. It supplies 58 of the 77 removals. It is not wrong, but it is not a
second opinion, and methods text should not present it as independent in-house QC.

Relatedly, **only two of Alamar's four QC metrics do any work**: `Detectability
< 0.90` (190 wells) and `|IC Median| > 0.40` (58). `IC Reads < 1000` fires on
zero wells and `Reads < 500000` on one, already caught by detectability.

---

## ⚠️ A defect this module found in `Metadata_Merge/Cohort_Accounting/`

`Cohort_Accounting` assigns one exit reason per well by a fixed precedence
(IC > PCA > burden). That is correct for a partition — the loss must sum — but its
per-axis numbers are **not axis totals**, and 16 of the 77 triaged wells carry
more than one axis:

| axis | on-axis (true) | `exit_reason` says | sole reason |
|---|---:|---:|---:|
| IC | 58 | 58 | 50 |
| PCA | **16** | 10 | 2 |
| burden | **25** | 9 | 9 |

So `exit_reason` undercounts the PCA axis by 37.5% and the burden axis by 64.0%.
Both tabulations are right answers to different questions; the failure mode is
reading one as the other. `triage_axis_overlap.csv` reports both side by side,
and `Cohort_Accounting/README.md` now points here.

The total of 77, and the 77/648 = 11.9% headline, are unaffected.

---

## Outputs

| file | what it answers |
|---|---|
| `qc_gate_inventory.csv` | does each gate fire, or is it **INERT** (untested, not proven safe) |
| `qc_gate_duplication.csv` | is an in-house axis a vendor metric re-derived (Jaccard + verdict) |
| `triage_axis_overlap.csv` | the true per-axis totals vs what a one-reason report says |
| `triage_axis_combinations.csv` | the full Venn, flattened |
| `gate_operating_points.csv` | what a "threshold" actually is — the quantised cuts available |
| `gate_cut_neighbourhood.csv` | nearest retained and nearest removed well on each gate |
| `read_gate_geometry.csv` | per plate-bay, whether the 500000 floor makes the lower bound unreachable |
| `qc_axis_attribution.csv` | each triage class against the known quality axes (d, AUC) |
| `well_position_effects.csv` | per-position rates + permutation p — this is what finds F7 |
| `well_position_effects_readflag.csv` | the same for read flags |
| `ic_row_gradient.csv` | IC recovery down the plate, per plate-bay |
| `triage_rate_by_site.csv` | site rates, permuted and depth-conditioned |
| `read_flag_rate_by_site.csv` | the same for read flags |
| `triaged_well_detail.csv` | one row per removed or marked well, everything on it |

---

## Method notes

**Triaged wells are placed on the intensity axis for the first time.** The
Extremes screen runs on the cohort, so removed wells have no `mean_INT` at all.
This module pools the post-QC and triage NPQ files and recomputes INT over the
union (4200 wells × 131 markers), so removed and retained wells are on the same
scale. Validated at **r = 0.9948** against the production
`extreme_sample_master.csv` on the 3608 shared wells.

**Permutation, not chi-square.** 77 triaged wells over 84 plate positions or 23
sites leaves most cells at 0 or 1, where an asymptotic test is not conservative
but wrong. Every clustering claim uses a 10,000-rep permutation holding the flag
count fixed, and the statistic is the **maximum** per-group rate, which needs no
further multiple-testing correction.

**Per class, never pooled.** See §1.

**`gate_operating_points.csv` caveat:** `count_value` is the median
`n_out_tot` among wells at that FDR. Because `fdr_min = min(fdr_up, fdr_lo)` and
the two tails carry different BH adjustments, the count is not a strict function
of the FDR across tails — it is a guide to the operating point, not an identity.

---

## Provenance

Findings came from a six-dimension multi-agent analysis on 2026-07-27 with
adversarial verification. **That run degraded**: 26 of 66 agents hit a session
limit, including both the synthesis agent and 24 of the verifiers. Of 58 raw
findings, 6 survived refutation, 28 were refuted, and **24 were never
adjudicated**. Everything asserted in this README is either (a) among the
verified survivors, or (b) re-verified directly — F7 and the row gradient were
recomputed by hand, and the row-gradient claim was **corrected**: the original
"monotone in 50 of 50 plate-bays" is false (0 of 50 are monotone); what holds is
a negative rank correlation in 50 of 50.

Refuted claims are deliberately not repeated here. The most common failure mode
among them was an interpretation outrunning arithmetic that reproduced exactly —
several were killed for circularity (explaining a PCA outlier with a panel-wide
statistic that PCA is computed from) rather than for wrong numbers.

---

## Where this lives

Canonical: `Alamar_Biomarker_QC_Repo/Secondary_QC/Triage_Review/`.

`scaffold_qc_run.R` copies code from the canonical repo into each run directory,
so **editing the copy under `QC_Runs/` is lost on the next scaffold.** Edit here,
then sync. Note `Secondary_QC` is not in the scaffold's default module list —
pass it explicitly:

```bash
Rscript scaffold_qc_run.R "QC_myrun" "Primary_QC,Metadata_Merge,Secondary_QC"
```

Run after Primary QC and Metadata_Merge, since it reads both.
