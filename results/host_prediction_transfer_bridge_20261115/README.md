# Butterfly host geography: observation-based validation of potential resource patterns

15 November 2026 | Source-locked reproducibility audit of the original 155-species GBIF forecasting analysis.

## Main scientific result

An existing larval host-resource map provides independent **predictive information** about later and geographically held-out **recording distributions** for butterfly species represented in model training. This is a distinct observational target from the 239-species resource reconstruction. It does **not** directly test plant introductions, resource consumption, butterfly colonization, shipping, or actual changes in island population sizes.

The original model lock predates training (29 October 2026). The response is the probability of an indexed butterfly record standardized to ten target-group recording days. Models compare background/previous recording, climate and land cover (M3) with host-resource indices (M4) and added host-family composition (M5). Primary test dates are 2019–2022; the host catalogue is **static contemporary** and was obtained retrospectively.

| Model evaluation | Species | Test species×cell pairs | Baseline M3 MSE | Host M5 MSE | M5 vs M3 relative MSE reduction | Species improved | Descriptive relative 95% interval |
|---|---:|---:|---:|---:|---:|---:|---|
| **Known species / later period** | 155 | 208,482 | 0.03206388 | 0.03074279 | **+4.120%** | 119 | +1.66 to +7.05% |
| **Known species / held spatial folds** | 155 | 208,482 | 0.03232691 | 0.03135865 | **+2.995%** | 115 | +0.56 to +5.47% |
| **Same held spatial folds, ≥100 km from training** | 155 | 102,728 | 0.03346781 | 0.03255225 | **+2.736%** | 111 | +0.26 to +5.23% |
| **Unseen butterfly species** | 155 | 208,482 | 0.03263301 | 0.03270106 | **−0.209%** | 95 | −2.95 to +2.63% |

Intervals above are algebraic conversions of the original 2,000-draw species × spatial-block fixed-prediction bootstrap **difference** intervals, dividing by M3 MSE. They do not include model-refitting, data sourcing, phylogenetic or plant-introduction uncertainty. The 100-km comparison is a **post hoc nested subset**, not an independent experiment.

## Does family composition matter beyond host counts?

For spatial held-out folds, M4 MSE was **0.03181539** and M5 MSE was **0.03135865**, a further **1.44% reduction relative to M4**. The original fixed-prediction interval for M5−M4 absolute MSE was **−0.0009219 to −0.0000403**. The same additional composition difference was **not clearly supported** in the temporal window (95% includes zero) or in untrained species.

This supports location-conditioned information in known host identities/family composition beyond their aggregate abundance indicators. It is not a causal explanation of island-specific butterfly establishment, nor proof that port-associated *introduced* hosts specifically caused the improvement.

## Previously unrecorded in the immediate preceding period

The model lock pre-specified 2015–2018 focal-species recording probability equal to zero. Those cells need not have been truly absent before that time. **154 species and 133,557 species×cell evaluations** qualify.

| Transfer mode | M3 MSE | M5 MSE | Relative reduction | Relative MAE reduction |
|---|---:|---:|---:|---:|
| Later held period | 0.02382452 | 0.02257809 | **+5.232%** | +2.69% |
| Held spatial block | 0.02453172 | 0.02323541 | **+5.284%** | +2.11% |
| Unseen species | 0.02431415 | 0.02449929 | **−0.761%** | −0.81% |

The positive held-period and held-spatial MSE differences had original descriptive bootstrap intervals wholly below zero. The unseen-species confidence interval crossed zero. Improvement is **not confined to cells where the butterfly was recorded in the preceding period**, but does not demonstrate actual new immigration, breeding or previously absent butterfly populations.

## Boundaries and robustness

Four original 155-species observation-quality sensitivities all retained positive MSE reductions: 20-day standardization **+2.98%** (159,084 pairs); 250-km candidate domain **+3.76%** (180,654); largest record provider excluded **+5.55%** (191,004); coordinate uncertainty at most 100 m **+5.34%** (135,532). These reuse the same archive and are **not independent replications**. The simple persistence model sometimes had lower **MAE**, so avoid universal claims of superiority on every loss.

The main 239-species host-resource reconstruction (+54.9% global, +118.9% strict islands vs +40.9% mainland) is also an **opportunity** calculation. The 155-species validation offers a different, observational measure of information in resource geography, **not a joint island×port×marine causal pathway**. The 51-island retrospective marine Circuitscape target and its strict geographic transfer failures remain separate.

## Reproduction and source provenance

- [Successful exact-source audit, CI, and frozen-model cohort checks](https://github.com/zuizui0223/chocho/actions/runs/43670497879)
- Audit implementation: scripts/audit_host_prediction_transfer_20261115.py. Original inputs are the frozen forecast validation all-metrics, paired comparisons, initially-unrecorded comparisons, model-lock and original metric validator in provenance/reviewer_defenses/results/host_forecast_validation_20261029/.
- The archived CI artifact contains the complete source-digest-bound summary and transfer_contrasts.csv.
- Original butterfly occurrence archive DOI: 10.15468/dl.umqmsf. The model inputs use contemporary host geography, not a plant map frozen before 2019.

**Editorial recommendation:** Promote the observation-based validation to the integrated GEB main Results, explicitly retain the non-transfer result for unknown butterfly species, and keep earlier port host indexing and ocean access separate from claims of new island establishment.
