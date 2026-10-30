# Island × port × host: exact recorded-host overlap and source-distance sensitivity

**2026-10-30 | Exploratory continuation of draft PR #102.** Primary recorded-source and original GEB manuscript results are unchanged. This is a new calculation on the **same** selected observational records, not new experimental or biological replication.

**Result:** early documented catalogue-compatible host recording overlaps later butterfly recording slightly more than the fixed-margin uniform table law predicts, **20 versus 17.0602 cells**, with exact inclusive upper-tail **6,487,980/147,117,600 = 0.04410064**. This numerical threshold is fragile to removing any informative geographic block. Under a post-result proximity-weighted model that reproduces the strong already observed nearest-recorded-source distance tendency, the expected overlap rises to **18.0493** and the descriptive upper-tail becomes **0.130546**. The host mechanism is therefore **not independently identified**.

## 1. Fixed units and what has actually been measured

The unchanged original candidate set contains 2,385 species×target-island cells, 66 cells with focal butterfly records in 2018–2025, and 93 early-recorded butterfly species in a selected 96-island source panel (28 frozen targets). A cell has documented prior known-host recording (`R=1`) if at least one conservatively taxonomically linked source-catalogue host species was found in the **indexed 2010–2017 plant records inside that same island polygon**. All other cells are `R=0`, **not evidence of true host absence**; for every 72 unmatched candidate cells in the currently variable support, some catalogue identities remain unresolved.

Observed recorded-host overlap: `H_obs = sum(Y_late * R_early) = 20` (9 forced-positive and 11 variable-positive cells). Only **17 recorded-host labels occur in the 89 cells with positive but non-unit occupancy** under the existing fixed geographic margins. All 88 fully resolved-catalogue no-match cells lie outside those variable placements and offer no true host-absence comparison.

The [uniform reference was fixed before calculating its H distribution](ISLAND_RECORDED_HOST_EXACT_PLAN_20261030.md). It preserves all structural zeros, every species×longitude×latitude-zone positive total, and each island's later-positive total. The 12 old latitude-zone blocks contain **147,117,600** admissible binary tables. Exact row-choice dynamic programming gives every integer-H frequency, and separate direct table enumeration independently reproduces the 52,542 / 10 / 28 / 10 nontrivial block table counts and H moments. No new data or simulations were used.

| Indicator | Exact result |
| --- | ---: |
| Observed H | **20** |
| Expected H under uniform frozen margins | **17.06018423** |
| Observed minus null mean | **+2.93981577 cells** |
| Permissible H support | **12–22** |
| Tables with H >= 20 | **6,487,980 of 147,117,600** |
| Exploratory exact uniform inclusive upper-tail | **0.044100637857** |

These figures quantify a *recording co-location* pattern relative to this selected reference. They do not establish a plant resource mechanism, true larval utilization, wild versus planted host status, island colonization, individual abundance or community growth. The >0.05 principle should not be used mechanically here: this is a post-pilot and multi-hypothesis exploratory continuation, with pronounced geographic concentration and incomplete taxonomy.

## 2. Why the marginal signal is not geographically robust

Only **three** of the 12 blocks permit H to vary; other blocks' H contribution is invariant even if some of their tables vary in butterfly–island placement. We computed all three leave-one-informative-block-out results as **diagnostics, not three new independent experiments**.

| Omitted informative block | Observed remaining H | Expected remaining H | Exact upper-tail in the remaining table law |
| --- | ---: | ---: | ---: |
| `longitude_1:tropical` (90°W–0° tropics) | 14 | 12.10 | **0.06428571** |
| `longitude_3:north_extratropical` (90°E–180° north) | 16 | 14.5602 | **0.21338739** |
| `longitude_3:tropical` (90°E–180° tropics) | 19 | 16.4602 | **0.06644276** |

The dependence on a few informative blocks argues against a universal "island hosts facilitate colonization" inference. Geographic resampling here does not create independent islands or remove observation bias.

## 3. Direct check of the proximity alternative

The preexisting nearest-*recorded*-island distance test found later butterfly positives preferentially close to earlier indexed butterfly occurrences. In the same **uniform finite fiber**, recorded-host H and the sum D of `log(1 + source coast gap km)` have correlation **−0.14884779** (an exact finite-fiber covariance/moment calculation). Some of the H pattern is therefore statistically aligned with the distance gradient.

A separately [fixed post-result distance-tilt sensitivity](ISLAND_RECORDED_HOST_DISTANCE_SENSITIVITY_PLAN_20261030.md) assigns each of the **same** admissible tables T relative weight `exp(−lambda * D(T))`, lambda >= 0. Enumeration of every complete block table allows exact weighted expectations at any **fixed** lambda. A single illustrative lambda is fitted **using D only** so the model reproduces the observed full D score, `490.55057679` (uniform expectation `507.19277559`). The fitted lambda is **0.76110801**. No H observation is used to choose lambda, but D and H still come from the same recorded outcomes.

| Source-distance preference | Expected H | Weighted probability of H >= 20 |
| --- | ---: | ---: |
| lambda = 0, uniform | 17.0602 | 0.044101 |
| lambda = lambda*/2 = 0.380554 | 17.5587 | 0.080799 |
| **lambda = lambda* = 0.761108** | **18.0493** | **0.130546** |
| lambda = 2lambda* = 1.522216 | 18.3981 | 0.171672 |

These are *model-dependent numerical sensitivity probabilities*, not additional exact p-values testing a composite null with fitted nuisance parameters. This is not a fitted biological dispersal law: the source-island panel omits many populations and recording histories, and near-source distances could be related to vegetation and observer coverage. However, **the originally borderline H-overlap is compatible with source-distance organization**, so it is not defensible to present it as demonstrated independent host facilitation.

## 4. Why the port hypothesis cannot be rescued by this test

Every target-island positive margin `c_i = sum_species Y_si` is **fixed** in the current conditional law. Thus, for any island-only port score `P_i`, the proposed additive statistic is

`sum_species,island Y_si * P_i = sum_island c_i * P_i`

and therefore is **mathematically invariant** across the reference tables. This is an identification limitation of the design, in addition to the earlier insufficient island-level recording support (23 supported late-period islands against the original 30-island requirement; near-continent/far-WPI class only one against three). The current WPI2025 geographic proximity does not itself quantify 2010–2025 shipping traffic. No port coefficient or claim about human-assisted introduction is warranted.

## 5. Scientific conclusion and next evidence needed

The nontrivial achievement is a reproducible transition from an undated host catalogue and broad geographical patterns to actual prior indexed host–butterfly co-recording and an exact within-panel comparison. The new evidence is *suggestive*, but **not yet sufficient to distinguish resource facilitation from previously occupied-source proximity, plant/butterfly record effort, and taxonomy coverage**. Actual butterfly breeding and plant introduction origin remain outside the measured endpoints.

To make the biological test possible, prioritize (i) better local plant records and defensible identity decisions for the 72 unresolved no-match cells in the variable support, especially the mixed species×zone rows, (ii) independent occupancy/detection or larval-use evidence where it exists, and (iii) actual historical port traffic plus more independently supported islands for the arrival pathway. Do **not** infer host absence from unmatched catalogue entries or drop poorly sampled islands to achieve a preferred result.

**Reproduce:** `python scripts/analyze_island_recorded_host_exact_20261030.py --output data/exploratory/island_recorded_host_exact_20261030/recorded_host_exact_full.json`; unit tests are in `tests/test_island_recorded_host_exact_20261030.py`. The source script reproduces both exact integer frequency and the distance-tilt numerical companion. The manuscript/anonymous-review bundle is explicitly excluded from this exploratory continuation.
