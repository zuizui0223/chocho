# Post-result source-history stratum and cluster sensitivity audit (2026-11-08)

**Status:** exploratory **post hoc analysis** of existing, fixed, public GBIF-index model predictions. No original raw Parquet recovery, fitted coefficient change, biological colonisation claim or manuscript alteration.

## Exact input (not a new raw source)
- GitHub Actions [successful source-locked run 42458801598](https://github.com/zuizui0223/chocho/actions/runs/42458801598), artifact ID `13609335251`, name `ncomms-public17-source-access-heldout-all225`.
- Downloaded ZIP SHA-256: `4518bc5502567fd0a1d57f25a54e206284ef5b6e1a7a0e9e1a8741777c7732c5`.
- All six 225-species×17-island heldout prediction CSVs from that artifact were read, and their B/D/C/DC Brier/log loss was recomputed against the source's `source_prediction_receipt.json`. Agreement within `1e-10`.
- The original **private** `work/barriers/global_original_species_island_year.parquet` is **NOT** in the artifact or this PR; original SHA-256 `8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f` is a provenance reference and was **not reverified** today.

## Source-history stratification of one-island-out *B→D* improvement

B: existing GBIF effort, prior-source counts and baseline geography. D: B plus geodesic nearest already-indexed conspecific source distance. Percentages are relative **Brier prediction-error reductions**, always evaluated on exactly the same held-out cells. Positive = D is better.

| Era | Previously indexed same species on a different island | Eligible cells | New GBIF first-index cases | B→D Brier gain |
| --- | --- | ---: | ---: | ---: |
| 2011–2017 | **Yes** | 864 | 57 | **+30.2407%** |
| 2011–2017 | No indexed source | 2,720 | 23 | +0.7225% |
| 2018–2022 | **Yes** | 1,039 | 17 | **+6.4619%** |
| 2018–2022 | No indexed source | 2,465 | 1 | −0.1623% |

The overall one-island-out B→D source-receipt effects are **+21.3784%** (2011–2017, 80/3,584 cases) and **+6.0910%** (2018–2022, 18/3,504). This stratification is a new diagnostic of *where* that already established retrospective score gain comes from, not a newly discovered independent colonisation outcome.

## Recomputed uncertainty and influence of original island groups

Re-sampled full held-out **island ID** clusters 5,000 times with replacement (seed `20261108`) and **never refitted models**. Intervals are on *absolute* Brier gain `mean((y−p_B)^2−(y−p_D)^2)`, not relative percentages. Because out-of-fold predictions remain fixed, these are conditional cluster uncertainty diagnostics, not source-survey or model-refitting uncertainty.

| Era and contrast | Brier percent gain | Absolute gain 95% island resample interval |
| --- | ---: | ---: |
| 2011–2017 B→D, all cells | **+21.3784%** | **[+0.002148,+0.007684]** |
| 2011–2017 B→D, prior source exists | **+30.2407%** | **[+0.009028,+0.031143]** |
| 2018–2022 B→D, all cells | +6.0910% | [−0.000018,+0.000697] |
| 2018–2022 B→D, prior source exists | +6.4619% | [−0.000057,+0.002372] |
| 2011–2017 **D→DC** ocean circuit beyond distance | **−1.8846%** | [−0.000857,+0.000105] |
| 2018–2022 **D→DC** ocean circuit beyond distance | +0.0833% | [−0.000088,+0.000086] |

Deleting each of the 17 original islands **without refitting** leaves all-period B→D Brier gains between **+18.055% and +23.677%** (early) and **+4.115% and +7.600%** (late). This influence diagnostic is narrower than uncertainty from geographic transfer: the original 2011–2017 **component-held-out** B→D still **worsens** Brier by 75.79%, and log-loss worsens when adding the ocean term D→DC in all six original holdout comparisons. Those unfavorable findings remain unchanged.

## Scientific interpretation and limitations
The earlier-period within-source-history proximity association is quantitatively stronger than the no-prior-source slice, consistent with indexed geographic source availability being informative in this selected 17-island observation panel. Neither the presence of an earlier indexed conspecific island nor geographic proximity proves actual movement, reproductive source activity, maritime travel, or physical shipborne host relocation. First-index records are not biological arrival or colonisation events; no-index rows are **not ecological absences**. The 17-island set was selected from a geographically restricted near-neighbor study and is not the original full 51-island or global population. The latest 18 events offer little independent statistical power.

The global **+54.9% potential host-resource expansion**, island/mainland **+118.9%/+40.9%**, and prior 155-species positive host-forecast results remain intact at their original respective population definitions. This diagnostic does not test a historical species-specific host × ocean-access interaction. A dated authentic host panel and source-controlled all-51-island observations are still needed before that ecological mechanism can be claimed.

### Reproduce from exact Actions artifact
Download the named artifact from run `42458801598`; for each `{era}_{holdout}_heldout_predictions.csv` compute squared error and Bernoulli log loss for `pred_B`, `pred_D`, `pred_C`, `pred_DC`. Check against the corresponding `source_prediction_receipt.json` score to tolerance `1e-10`. For island holdout, stratify `prior_source_islands>0`, resample entire `landmass_id` clusters with seed `20261108` for 5,000 draws and calculate the mean per-cell Brier gain each draw; take 2.5 and 97.5 percentiles. The original source and model coefficients do not change.

**Scope:** follow-up inference audit only; independent GEB paper, JBI Draft #151, and NC prior-source main comparison remain unchanged.
