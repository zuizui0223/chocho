# Original 51 islands: candidate-host identity beyond richness

## Result (exploratory sensitivity, 2026-11-07)

**Candidate-host Jaccard retains held-island prediction information after adjusting for both island-level candidate-host taxon counts and the native AIS/wind marine circuit.** Likewise, the circuit terms retain an increment after host identity and marginal richness are represented. Both improvements occur for the later butterfly community target, 2011–2022.

| Later (2011–2022) held-island model comparison | MSE reference | MSE augmented | MSE reduction | Islands better/worse | MAE reduction |
|---|---:|---:|---:|---:|---:|
| Prior C → HC: host Jaccard beyond circuits | 0.00600023 | 0.00589635 | +1.7312% | 31/20 | +1.4752% |
| C → CR: candidate-host taxon counts beyond circuits | 0.00600023 | 0.00584999 | +2.5039% | 31/20 | +1.4660% |
| **CR → HCR: candidate-host composition beyond richness + circuits** | **0.00584999** | **0.00571656** | **+2.2809%** | **34/17** | **+1.5127%** |
| Prior H → HC: circuits beyond host Jaccard | 0.00630545 | 0.00589635 | +6.4880% | 39/12 | +4.5286% |
| **HR → HCR: circuits beyond host Jaccard + richness** | **0.00595497** | **0.00571656** | **+4.0036%** | **36/15** | **+2.5097%** |
| HC → HCR: richness beyond the prior joint model | 0.00589635 | 0.00571656 | +3.0492% | 35/16 | +1.5035% |

The prespecified primary sensitivity is CR → HCR. Its positive increment is *not* a direct increase from the prior 1.7312%: those percentages have different reference models and denominators. The same signs for MAE provide a second error-metric check rather than an independent population confirmation. Candidate-host taxon counts also have predictive information; they are not rendered irrelevant by a positive composition increment.

### Earlier (2001–2010) temporal context

| Comparison | MSE reduction | Islands better/worse | MAE reduction |
|---|---:|---:|---:|
| CR → HCR: host identity beyond richness and circuits | +4.7484% | 33/18 | +2.6192% |
| HR → HCR: circuits beyond host identity and richness | −0.9678% | 24/27 | +0.1332% |
| C → CR: richness beyond circuits | +4.8006% | 31/20 | +3.9289% |

The earlier circuit increment is adverse on MSE and slightly positive on MAE. The different responses are not a formal test of changing migration mechanisms, and the present-day AIS/wind surfaces are not historical routes.

## What exactly was tested

Use the [frozen protocol](PROTOCOL.md) and [read-only runner](../../scripts/run_island51_richness_control.py), stacked on [island atlas PR #129](https://github.com/zuizui0223/chocho/pull/129). For each candidate host list, take log(1+number of introduced larval-host candidate taxa). Two response-free pair controls are its endpoint sum and absolute difference. They enter C, H and HC as CR, HR and HCR, respectively. C is the parent geography + uniform-ocean + native AIS/wind model; H includes host Jaccard instead of AIS/wind; HC includes both. Parent geographic controls already include log-distance terms, original port-count sum/difference, and the earlier butterfly similarity only for the later target.

This reuses **51 islands**, **1,275 unordered pairs**, **194 original biological ports** and the same prior 182 marine-qualified terminal ports. Each fold trains on 1,225 other-island pairs and tests the 50 pairs involving the held-out island. Every model and target retains two endpoint-specific predictions for each unique source pair. There are **612 model folds**, **30,600 prediction records** and **14 fixed comparison summaries**.

The source files are verified by the previous parent input lock; the parent fold table is pinned to the prior atlas lock. All **306 C/H/HC parent model-fold MSE and MAE values reproduced**, with maximum MSE difference 7.12e-17 and MAE difference 2.09e-16. Only two host-richness controls are added. The prior observations, original 225-butterfly panel, candidate-host lists, circuits and previous model files remain unchanged.

Run (repository root):

    python -m pip install -e ".[test,analysis]"
    python -m pytest -q tests/test_run_island51_richness_control.py
    python scripts/run_island51_richness_control.py --out-dir /tmp/chocho-richness
    python results/island51_richness_control_20261107/verify_predictions.py --results-dir /tmp/chocho-richness

The dedicated [GitHub workflow](../../.github/workflows/chocho-island51-richness-control.yml) archives all 612 model folds, 30,600 individual predictions, 14 comparisons, receipt hashes, and an independent source-code-separated loss reconstruction. Artifact retention follows GitHub Actions settings. See [Draft PR #130](https://github.com/zuizui0223/chocho/pull/130) for the run and its review checks.

## Interpretation and limits

The stronger claim supported here is **incremental out-of-island predictive information for candidate-host species identity beyond the chosen two marginal richness terms**, conditional on the specified geography and circuit predictors. This does *not* fully control plant sampling effort, climate, source-region histories, island species-area effects, or host quality. Recorded candidate hosts are not demonstrated larval feeding resources. Host lists pool 2001–2024 and can postdate the later 2011–2022 butterfly response. First-indexing does not demonstrate actual arrival, survival, reproduction, or colonization. The AIS density proxy does not identify actual bilateral ship routes and the wind circuit is not directional dispersal. All folds share source pairs and overlapping training data; island improvement counts, fold percentages and error metrics are descriptive, not independent p-values or causal evidence.

The 51-island richness sensitivity belongs to the island–port paper track. It does not alter the GEB global +54.9% resource reconstruction or claim that the same reused input is a distinct global discovery.
