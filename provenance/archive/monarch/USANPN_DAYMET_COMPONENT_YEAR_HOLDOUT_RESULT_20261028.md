# Temperature + frost forecast monitored milkweed leaves beyond seasonal calendar alone

**2026-10-28, source-complete original-data predictive analysis, not future plants or causality.** The present study is a posthoc refinement of a completed positive weather result on PR #93. No original GEB main/manuscript edits.

## Exact-source study
USA-NPN 2018–2025 six Asclepias taxa, 27,057 direct Leaves YES/NO plant×date records; all eight NPN yearly source hashes matched previous source receipts. Daymet 1-km original Tmax/Tmin/precipitation 2018–2025: all **72/72** original effort-selected monitoring sites re-retrieved with identical combined Daymet SHA256. **20,770** full-weather observed plant visits; no site selection or original sample counts changed. Original station IDs, locations and raw source observations were never committed.

## Unseen-site comparison: original five whole-station folds

| Covariate increment over species+season+geography+year | Site-equal held-out log-loss reduction | Positive folds |
|---|---:|---:|
| 30-day mean temperature + frost days | +0.024562 | 5/5 |
| 30-day precipitation only | -0.001734 | 1/5 |
| Temperature + frost + precipitation | +0.025406 | 5/5 |

**Temperature/frost improved 5/5 held-out station folds; precipitation-only improved 1/5**, with a negative source-mean difference. The full three-weather predictor model reproduces the original PR #93 5/5 improvement (+0.025406). Positive predictive performance is not causal evidence for temperature or zero precipitation effects.

## Whole-calendar-year extrapolation (2018–23 training → 2024–25 testing)

The held-out 2024–2025 subset contains **5,198** plant×date visits at **49** original stations, five taxa; 2,521 recorded YES and 2,677 recorded NO. Earlier period training contains 15,572 visits. Time out-of-sample is NOT site out-of-sample (the separate spatial folds test that).

| Model | Held-out year site-equal log-loss | Improvement vs baseline |
|---|---:|---:|
| Season+species+geography+year | 0.304754 | +0.000000 |
| Plus temperature/frost | 0.297329 | +0.007425 |
| Plus precipitation only | 0.305226 | -0.000472 |
| Plus full weather | 0.300247 | +0.004507 |

**Temperature plus frost passes the frozen 5-fold and whole-year holdout gates** (2024–2025 log-loss 0.304754 → 0.297329; +0.007425). This supports real interannual weather predictability of monitored live-leaf status, not a fitted future response function.

## Required forecast transport caveat
If advancing to conditional CMIP6 warming scenarios, use real CMIP6 future-minus-1970–2000 monthly temperature deltas added to the actual 2018–25 daily Daymet values, with daily precipitation fixed and frost days recomputed from daily Tmin. Report the original observation dates as held fixed: never interpret predicted YES fractions as a population-wide leafing calendar, use counts of monarch generations, or project future leaf senescence without outcome validation. Monthly GCM climatologies have no daily frost-days data. Audit climate extrapolation and source locality inclusion.

## Reproducibility and citations
- [Original validated action](https://github.com/zuizui0223/chocho/actions/runs/40620259854), SUCCESS, three source-independent tests.
- Frozen original source result at `provenance/archive/monarch/results/USANPN_DAYMET_COMPONENT_YEAR_HOLDOUT_20261028.json` SHA256: `4a935e60b4a88d81520acf32d7deaec66fee792b3ebb0a169a3f6cd9add76bc3`.
- Frozen pre-outcome follow-up contract at `provenance/archive/monarch/USANPN_DAYMET_COMPONENT_AND_YEAR_HOLDOUT_20261028.md`. No raw station/observer/plant/day-data exports.
- USA-NPN credit: Data were provided by the USA National Phenology Network and the many participants who contribute to its Nature’s Notebook program. Daymet single-pixel DOI:10.3334/ORNLDAAC/2361.
