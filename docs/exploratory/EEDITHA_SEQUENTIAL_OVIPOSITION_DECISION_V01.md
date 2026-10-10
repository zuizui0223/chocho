# Short-term repeated host choice: exact-original sequential prediction test

**Date:** 2026-10-10  
**Status:** SOURCE VERIFIED / PRIMARY PREDICTIVE GATE FAILED / NO EVIDENCE FOR CAUSAL OR EVOLUTIONARY HYSTERESIS  
**Scope:** exploratory only; the GEB manuscript and PR #38 are untouched.

## Independently sourced evidence

- Haan, Bowers & Bakker (2021), *Scientific Reports* 11:992, DOI https://doi.org/10.1038/s41598-020-80413-y.
- Original deposited oviposition data: https://zenodo.org/records/4318182/files/oviposition.csv?download=1 ; DOI https://doi.org/10.5061/dryad.612jm642h ; **verified MD5 `6d9770cce9cacfccd98f7d19c67d0ec5`**.
- Exact run: https://github.com/zuizui0223/chocho/actions/runs/38012350290 ; result artifact: https://github.com/zuizui0223/chocho/actions/runs/38012350290/artifacts/11654840600.
- Protocol was committed before actual source download or analysis results: `EEDITHA_SEQUENTIAL_OVIPOSITION_PROTOCOL_V01.json` (initial `d900156`, corrected source-semantic gate `5c19e17`); final analysis script `scripts/analyze_eeditha_sequential_oviposition.py`.

**Original experimental design:** 29 mated females, repeatedly assigned three-way choice pots containing each of two ancestral *Castilleja* species and adopted exotic *Plantago lanceolata*. There were 145 source trial rows, with 5,417 eggs recorded on these three plants and 438 eggs on other surfaces. The original paper already demonstrated that host preferences changed over time. Source columns `CAHI.cl`, `CALE.cl`, `PLLA.cl` count **egg clutches**, not how many host plants were offered. They were never used as same-trial predictors.

## Source-checked prospective predictive comparison

Within each female's trial sequence, require eggs laid on the currently observed trial and on an earlier egg-producing trial. Label the current choice Plantago-dominant if Plantago eggs **strictly exceed both ancestral hosts combined**. The previous trial's Plantago share is the incremental predictor. The baseline uses calendar date and trial number. Models used fixed C=1 logistic regularization, standard scaling and **leave-one-female-out** (matriline-blocked) predictions; independent replication unit is the female, not an egg.

| Quantity | Exact result |
|---|---:|
| Source trial rows | 145 |
| Unique original females | 29 |
| Eligible repeat trials | 93 |
| Females contributing repeat trials | 27 |
| Plantago-dominant repeat trials | 33 |
| Baseline heldout log loss | 0.6700996262 |
| Add preceding Plantago egg share: heldout log loss | 0.6656197616 |
| Baseline minus added-history log loss (positive = better) | **+0.0044798645** |
| Female-cluster bootstrap 95% interval (2,000 resamples, seed 20261010) | **[−0.0181745, +0.0256644]** |
| Brier improvement | +0.00128857 |

The improvement is small and uncertain. The interval includes zero. **The preregistered predictive gate FAILED.**

A directional negative control uses future rather than past choice, never represented as a valid prospective model. On the same 66 trial subset, adding *past* choice changes baseline log loss by **−0.002083** (worse) and adding *future* choice by **−0.005135** (also worse). The past model is less poor than the future one on this subset but **neither improves on baseline**. It does not support behavioral memory.

### Limits and decision

1. This is **one butterfly population**. The 29 females are not 29 independent populations.
2. Each trial includes all three species but pot identity, plant architecture, recent oviposition, female age and plant physiology can vary. Without randomized history exposure, even a large autoregressive effect could reflect stable individual preference or nonrandom resource experience.
3. These are **egg-placement** outcomes, not individually tracked offspring, fitness, local demographic resilience, or evidence of actual post-removal fallback.
4. Original 2021 publication already showed temporal shift in host preference; there is no claim to discover that fact here.
5. We do **not** adjust the model, select favorable dates/mothers, redefine prior choice, optimize regularization, re-bin the outcome or retest alternative windows after seeing the failed gate.
6. The exploratory short-term behavioral carryover route is closed. This result neither falsifies nor proves the stronger cross-generational *history × controlled-removal* ecological hysteresis hypothesis; that requires an independently designed experiment with randomized resource-history treatment where feasible, local plant-stage controls and egg-to-flight-capable-adult response.

**Interpretation:** Global plant redistribution mechanically creates mapped resource opportunities, but this experiment shows why mechanistic behavioral inheritance cannot be safely inferred from host identity or repeated egg counts. The scientifically testable next step is not more species-level geography or an ad hoc lag model: it is to cross replicated experimental resource histories with a later standardized ancestral-host-only environment, measuring viable offspring and switching behavior.

**GEB manuscript stays unchanged.**
