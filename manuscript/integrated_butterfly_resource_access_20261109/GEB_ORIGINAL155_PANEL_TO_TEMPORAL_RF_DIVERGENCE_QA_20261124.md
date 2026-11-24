# Original155 rebuilt panel → frozen temporal RF: numerical failure and bounded diagnostic (24 November 2026)

**FAIL as an end-to-end reproduction, despite complete model execution.** The [first full derived-source-to-155-model CI #45320236263](https://github.com/zuizui0223/chocho/actions/runs/45320236263) executed the unmodified original panel builder from the 13 source inputs (old 10 and newly recovered 3), then the unmodified original temporal RF model source. All **842,964 original keyed rows and 96 columns** numerically matched the frozen panel under the declared panel comparison tolerance. The resulting temporal RF predictions **did not** match the originally archived 208,482 keyed observations.

| Diagnostic | Archive | Fresh model from regenerated panel |
|---|---:|---:|
| M5 minus M3 MSE reduction | **+4.1201756079%** | **+3.9356441905%** |
| 208,482 keyed target rows | Same | Same |
| Full panel keyed values (96 columns) | Matched within tolerance | Matched within tolerance |
| M0–M5 prediction comparisons above tolerance | – | **208,469 / 208,472 / 208,478 / 208,480 / 208,478 / 208,478** |
| Persistence model predicted values | Original | **0 mismatches** |
| Numerical comparison groups | – | **66/175 failed** |

The tool saved an aggregate-only receipt as [artifact #14887466734](https://github.com/zuizui0223/chocho/actions/runs/45320236263/artifacts/14887466734). The original frozen model results and scientific paper were **not** overwritten. Positive direction of the M5 relative MSE gain alone is **not** an exact reproduction.

**Competing explanations:** the all-column panel check first sorted by species×cell×period keys and applied a numerical tolerance, so it did **not** verify original row ordering or exact float32 training-feature identity. Fixed-seed random-forest bootstrap samples are indexed by training row position; if the newly recovered original candidate file changes species order, an identical keyed table with different physical ordering can yield many different fitted predictions. This remains an explicitly **untested hypothesis**, not a demonstrated explanation from this failed run. A second, predeclared [row-order/float32 diagnostic #45321468031](https://github.com/zuizui0223/chocho/actions/runs/45321468031) compares ordered keys and float32 values, then (only as a counterfactual diagnostic) trains the newly rebuilt panel ordered by the original saved panel's **keys only**, never its predictor values. This reference-order intervention must not be represented as a fully independent historical builder-order reproduction.

**Science and data-rights gates remain:** original155 temporal and spatial refits **from the archived original panel** were already successful, and all 96 predictor columns were reproduced from the now recovered sources; this additional chain test fails. Unseen-species fold 3 and two sensitivity numerical discrepancies remain unresolved. The three new source CSVs' historical creation chain, raw original GBIF-to-recording-day preprocessing, whole SI, source rights and stable double-anonymous reviewer access remain outside the evidence. PRs #179 and #189 are Draft, original effect values unchanged.
