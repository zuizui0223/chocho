# Do two host plants buffer local butterfly persistence after loss of one host?

**Date:** 2026-10-10  
**Scope:** independent ecological exploration; no changes to submission PR #38  
**Decision:** SOURCE GATE PASSED; OBSERVATIONAL INTERACTION PREDICTION FAILED; STOP ADDITIONAL POST-HOC MODEL SEARCH.

## Recovered independently archived original field data

The source-blocked DiLeo et al. (2024) 2004–2013 file is **not** used in this analysis. Instead, a separate official [Zenodo deposit 4987060](https://zenodo.org/records/4987060), based on [Schulz, Vanhatalo & Saastamoinen 2020](https://doi.org/10.1111/ecog.04799), was downloaded through GitHub Actions with source MD5 **69122a1d82b1fb970fb6638b02da3db4** and SHA256 **c70f43544f61e7273afebfa42ed1d488c3b945efe4d5c68e6557030087ee7fd9** verified.

Source: `data/survey_data.tsv` in `ECOG-04799.zip`: **74,065 unique patch × survey-year records**, **4,902 patches**, calendar years **1999–2018** (20 calendar years in file, despite paper's description as a 19-year study). 775 survey rows have missing scores for each focal host; missing scores were **not** treated as host absence. Independent source audit workflow: [38014613652](https://github.com/zuizui0223/chocho/actions/runs/38014613652); frozen population/event-support audit: [38014769330](https://github.com/zuizui0223/chocho/actions/runs/38014769330).

**The empirical target is `Melitaea cinxia` local patch non-detection one year later, not GBIF range shift, exotic-plant introduction, fitness, or evolutionary hysteresis.**

## Genuine three-year risk set

For each patch, link exact consecutive autumn surveys `t−1,t,t+1`, require larval nests detected in autumn `t`, and classify `Veronica spicata` decline between `t−1` and `t`. The candidate alternative host `Plantago lanceolata` is classified high when its original 0–3 abundance code is 2 or 3, low otherwise. The next-autumn response is presence of **zero observed nests**, not proven local extinction.

Of **10,061** eligible annual triples (from **2,043 occupied-risk-set patches**, **18 index years**), the decline cells contain:

| Decline of Veronica, current Plantago | Patch-years | Different patches | Future non-detections | Raw fraction |
| --- | ---: | ---: | ---: | ---: |
| Decline, backup-high | 730 | 429 | 262 | 35.89% |
| Decline, backup-low | 163 | 115 | 68 | 41.72% |
| No decline, backup-high | 8,418 | 1,911 | 3,642 | 43.26% |
| No decline, backup-low | 750 | 485 | 409 | 54.53% |

The ~5.8 percentage-point decline-stratum raw contrast is **not an estimated functional-buffering effect**: it may reflect existing host abundance, population size, site area, spatial location, management, season/drought, detection and strong imbalance between support groups. The screen required 40 events and 30 distinct patches in both decline strata, 3 years and ≥20 future zero outcomes; all gates passed.

## Adjusted fixed model: two replication-bound validation regimes

Model features set *after source/observed event-count audits, before adjusted fit*: current log nest count, current/previous host abundances, source decline and backup indicators, year terms, patch area/coordinates, grazing, host desiccation, and training-only missing-value imputation. Incremental parameter: **Veronica decline × high Plantago**, with main effects already present. Logistic L2 penalty C=1, not outcome-tuned. Raw data and exact source are retrieved afresh by both verified workflows.

### 1. Leave-one-index-year-out retrospective interpolation — NOT a forward prediction

[Run 38014972470](https://github.com/zuizui0223/chocho/actions/runs/38014972470) used all 18 index years as heldout folds, but **its training folds included future years**, so it is only a calendar-year-interpolation diagnostic and must not be called prospective or an independent temporal hindcast.

- Baseline heldout logloss: **0.6525925**.
- Augmented heldout logloss: **0.6528456**.
- Added-term heldout logloss benefit: **−0.0002531** (worse, rather than better).
- Patch cluster bootstrap 95% interval: **[−0.0005289,+0.0000014]**; year cluster **[−0.0006628,+0.0001256]**.
- Benefit positive in **8/18** heldout index years. Frozen model-comparison gate failed.

### 2. Strict rolling-origin forward validation — PRIMARY decision

A source-independent methodological correction was [committed before the first adjusted model result was inspected](https://github.com/zuizui0223/chocho/blob/exploration/realized-host-window-v01/docs/exploratory/MELITAEA_19YEAR_STRICT_FORWARD_VALIDATION_V01.json). For each test year 2010–2017, fit only earlier index-year triples whose `t+1` outcomes had already been observed by the test year. No training year follows or coincides with its test year.

[Exact forward workflow 38015094308](https://github.com/zuizui0223/chocho/actions/runs/38015094308), source matched, three forward chronology tests and model production successful:

- Forward test **4,543** records from **1,522** patches over **8** calendar years.
- Baseline forward logloss **0.7086670** vs augmented **0.7089529**.
- Incremental forward logloss benefit **−0.0002859** (negative means augmented model worse). Brier improvement **−0.0000914**.
- Patch-cluster bootstrap 95% interval **[−0.0009382,+0.0003398]**; year-cluster **[−0.0016323,+0.0006923]**.
- **5/8** forward years individually positive, below frozen minimum of 6/8; neither CI excludes zero.
- **Primary forward gate FAILED.**

Notably test year 2017 predicts 2018, the year of a documented major climatic population crash in this same Åland metapopulation (Opedal et al. 2020). This is a substantial climate/detection sensitivity, not permission to omit 2018 post hoc. The recorded future outcome is nest non-detection, subject to imperfect field detection.

## Ecology and novelty assessment

1. **Already established:** raw host amount increases butterfly occupancy and can reduce local extinction; Schulz et al. (2020) found a positive simultaneous-two-host association after environmental controls, and Opedal et al. (2020), DOI https://doi.org/10.1002/ecy.3186, explicitly studied host-specific abundance, connectivity, colonization and extinction, including parasitoid metapopulation processes.
2. **Actually tested here:** whether a recently declining *Veronica* host modifies the predictive value of high *Plantago* beyond static and lagged plant counts, nest abundance and local environmental metrics.
3. **Not supported:** additional *generalizable, forward predictive* host-decline × backup-host signal under the frozen response and model, notwithstanding the large independent sample.
4. **Not tested:** experimentally observed host switching, larval survival after a causal removal, host-specific enemy exposure, changed maternal preference, local genetic history or introduced-vs-native host adaptation. Static host-option geography is not equivalent to these endpoints.
5. **Future evidence needed for a new ecological mechanism:** paired host-manipulation experiments controlling edible biomass and phenological stage, replicated butterfly families/populations, parasitoid exposure, individual egg-to-flight-capable-adult outcomes, and randomized current loss following randomized history regimes where feasible. The already defined 2×3 factorial design remains the distinct experimental route, not a result.

**STOP RULE ACTIVATED:** Do not choose different future years, alter host score thresholds, tune penalties, remove drought years or choose a more favorable lag after seeing these null results. A new hypothesis must introduce genuinely independent biological observations or an experiment, not another rearrangement of the same occupancy/host-cover associations.

## Code, audit and manuscript separation

- Exact source inventory: `scripts/audit_melitaea_zenodo_alt_source.py`.
- Actual matched-source support counts: `scripts/audit_melitaea_19year_fallback_support.py`.
- Nonforward sensitivity (not to be called prospective): `scripts/analyze_melitaea_19year_host_fallback_prediction.py`.
- Primary forward analysis: `scripts/analyze_melitaea_19year_strict_forward.py`.
- Frozen original protocol and fit specification: `MELITAEA_19YEAR_DYNAMIC_FALLBACK_PROTOCOL_V01.json`, `MELITAEA_19YEAR_FALLBACK_FIT_ADDENDUM_V01.json`, `MELITAEA_19YEAR_STRICT_FORWARD_VALIDATION_V01.json`.
- This branch must not be merged into GEB PR #38 absent separate scientific justification. GEB concerns potential resource geography and remains unchanged.
