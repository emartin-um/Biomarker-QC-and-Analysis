# Detectability & LOD (Secondary QC)

**Detectability measures the blanks, not the samples.** This module makes that
visible and fixes the estimator behind it.

Implements the 2026-07-13 recommendations **§2** (standing detectability-by-bay
panel, with NC background and pooled-LOD columns) and **§8** (variance-stabilised
LOD). Base R only.

```bash
Rscript run_detectability_lod.R
```

## What detectability actually is

```
detectability  = (# targets with NPQ > targetLOD_NPQ[target, bay]) / n_targets
targetLOD_NPQ  = log2( mean(2^NC) + 3 · SD(2^NC) )
```

The NC term comes from just **four blank wells per plate-bay**. So a bay "fails"
detectability when **its blanks read high**, not when its samples fail. That makes
it an **assay-background** signal — an analysis and Alamar question — never a
sample-quality flag or a re-run trigger.

Verified on the 50-plate run: bay detectability correlates **−0.645** with blank
background. The metric is largely reporting the blanks.

## Two structural problems, both quantified here

**1. A four-point SD is fragile (§8).** Leave one blank out and recompute:

| leave-one-blank-out swing in a bay's LOD | log2 |
|---|---|
| median | **0.33** |
| 90th percentile | **1.36** |
| worst target × bay | **6.73** |

One well moves the threshold that far *before any sample is considered*.

The fix keeps the per-bay **mean** — the real background signal — and replaces the
four-point SD with that target's **typical across-bay SD** (median, so one
pathological bay cannot inflate everyone else's LOD). Effect on the worst bays:

| bay | detectability as-is | stabilised | wells under 0.90, as-is → stabilised | blank-background z |
|---|---|---|---|---|
| 20260520-1309_Bay3 | 0.890 | **0.903** | 67.4% → **41.9%** | +3.56 |
| 20260519-1245_Bay2 | 0.894 | **0.907** | 55.8% → **30.2%** | +3.30 |
| 20260514-1222_Bay3 | 0.895 | 0.906 | 55.8% → 41.9% | +1.85 |

**2. The gate is a knife-edge.** Detectability is quantised at k/n_targets and the
0.90 line falls between k = 116 (0.8992, FAIL) and k = 117 (0.9070, PASS) — it
bisects the mode, so a small background shift flips a large share of a bay's wells.
That is why the two flagged bays cross the gate under stabilisation without any
change to the samples.

## 🚫 Never gate or drop on detectability

And never gate a *marker* on per-bay detection status without first checking that
marker's blank trend across the batch. `marker_below_lod_by_lot.csv` is that check.
On this run the calls that move most between lots:

| target | 2025 | 2026 | shift |
|---|---|---|---|
| **SNAP25** | 5.8% | **91.5%** | −85.8 |
| **YWHAG** | 13.8% | 55.1% | −41.3 |
| **pTau-217** | 2.4% | 33.0% | −30.6 |
| BD-pTau-217 | 52.6% | 25.5% | +27.1 |

pTau-217 is the known case, but **SNAP25 and YWHAG are worse** — SNAP25's detection
calls are essentially dead in 2026 while its quantitation is unaffected. Use
continuous NPQ for any near-floor marker.

## Outputs

| file | what it answers |
|---|---|
| `nc_background_by_bay.csv` | how high does each bay's blanks read? z vs the study norm, flagged above +2 SD |
| `lod_comparison.csv` | per bay × target: LOD as shipped vs stabilised, and the gap |
| `lod_leave_one_out.csv` | how far one blank well can move each LOD |
| `detectability_by_bay.csv` | detectability under both LODs, with the blank background beside it |
| `marker_below_lod_by_lot.csv` | per marker, % of wells under their own bay's LOD, by lot |

## Provenance

The module reproduces Alamar's own `targetLOD_NPQ` exactly — **median |difference|
0.000000 across 6,450 plate × target pairs** — so the stabilised variant differs
from the shipped one only in the SD term, by construction.

Inputs: `Secondary_QC/Replicate_Analysis/output_files/nc_reps_NPQ.csv` (blanks),
`Primary_QC/input_files/Sample_QC.csv` (the sample → plate-bay map), and optionally
`Replicate_Analysis/input_files/NPQ_20260623.csv` and `Annotation_Targets.csv`. It
exits cleanly with a note if the required two are absent.
