# Historical port connectedness and butterfly community prediction: 17 original islands and an influential mixed-country island group

**Completed 2026-11-08. Exploratory observational source-locked test, NOT an independent causal shipping experiment.** No existing global GEB manuscript or original host/port/previous AIS finding changed. Draft stacked on PortWatch original source-geography [PR #148](https://github.com/zuizui0223/chocho/pull/148), which contains 24 independently WPI geo-verified harbors across 17 original island/landmass groups.

## Research question

Does **historical (2006–2010) UNCTAD PLSCIQ PORT-LEVEL LINER SHIPPING CONNECTEDNESS INDEX** have predictive information for the recorded similarity of **later (2011–2022) first-indexed butterfly communities**, beyond original geography, historical WPI botanical port counts, 2001–2010 general-plant GBIF buffer observation effort, and earlier 2001–2010 butterfly community similarity?

This is a legitimate temporal-ordering contrast (indexed shipping connectedness before the later butterfly record window), **not** a voyage, cargo, confirmed larval host arrival or biological colonization test. Earlier 2001–2010 butterfly community similarity is shown only as a **contemporaneous index comparison** because it overlaps 2006–2010 PLSCIQ.

## Source quality and fixed population

The original recovered UNCTAD PLSCIQ old provider edition was reacquired in [PR #146](https://github.com/zuizui0223/chocho/pull/146) with **all 3 original response SHA-256 pages exactly matching** the earlier 2026-10-30 archived original database. It covers 1,300 port time series, 81,900 original port×quarter records and 57,312 nonmissing observations from 2006-Q1 through 2021-Q3 on old index base (Maximum Q1 2006=100). It is **not a direct port-to-port voyage or transported seed count**.

From the original independently PortWatch-verified 24 WPI ports, **22 have at least 16 of the 20 2006–2010 quarterly PLSCIQ source values** and span **17 original island/landmass groups**. The two remaining verified ports with insufficient historical index data are retained in the [original 24-source-port quality ledger](original24_geoverified_ports_historical_PLSCIQ_quality.csv) and *never called traffic-free*. The 17 original island groups generate **136 unique unordered island pairs** from the original independently indexed 1,275-pair butterfly population. The full [original 51-island index eligibility ledger](original51_historical_PLSCIQ_island_eligibility.csv) preserves original group IDs and all missing source-index groups (not biologically unconnected islands).

For each of 17 source-eligible island groups, the **prespecified primary source metric is the median of the 2006–2010 port means** for index-qualified original WPI harbors, with a prespecified maximum-port alternative for hub sensitivity. Important: original island group **9** is the geographical landmass containing **four historically well-indexed source WPI ports in both Indonesia and Malaysia** (Samarinda, Banjarmasin, Pontianak, Kuching). That original source group was retained as one 51-island study unit; excluding it was a source-frozen sensitivity, not a redefinition of the island population after seeing outcomes. All other historical-index-eligible original islands remain fixed.

The original source response is pairwise Jaccard of butterflies **FIRST INDEXED** in the source periods, not abundances, butterfly breeding, actual arrival or individual dispersal. The historical port index may proxy human trade/global connectedness, vegetation introductions, observation effort or geography, not identify any one mechanism.

## Strict temporal and geographic holdout

For **each** test island pair (A,B), its source butterfly response is predicted from a regression trained ONLY on pairs with **neither A nor B present in either training-pair endpoint**. Across all 17 eligible islands there are **136 separate original pair tests, each fitted on the other C(15,2)=105 island pairs**. For a source-frozen sensitivity **excluding island group 9**, 16 islands yield **120 separate original pair tests, each trained on the other C(14,2)=91 pairs**. Feature standardization uses training rows only, with OLS intercept and 0–1-clipped Jaccard predictions.

Primary prefit BASE predictors (same original source inputs for both models): original great-circle port distance and its squared log, sum/difference of log original WPI port counts, representative island latitude separation, early 2001–2010 general-plant WPI buffer RECORD-count sampling effort sum/difference. The later 2011–2022 butterfly outcome additionally adjusts source earlier 2001–2010 butterfly Jaccard. **Original 2001–2024 retrospective host network Jaccard, modern AIS2015–2021 and GLSN2018 modelled co-service network are explicitly excluded** from these time-oriented primary regressions, so no future candidate-host record is used.

Four precommitted models: BASE, MED_SUM (+ median historical port index endpoint sum), MED_SUM_DIFF (+ endpoint absolute index difference), MAX_SUM_DIFF (+ maximum source port connectivity sum/difference). Both earlier and later butterfly responses; 17-island full source and prespecified 16-island no-group9 source. Exactly **2,048 unique heldout model predictions**, **256 distinct two-island-unseen training/test folds** and **20 original heldout model contrasts**.

## Main result: index predictive gain is not stable across source-island composition

| PLSCIQ historical feature added beyond BASE | **Later butterfly** MSE gain, all 17 islands/136 pairs | **Later butterfly** MSE gain, excluding group9/120 pairs |
|---|---:|---:|
| Median 2006–2010 index endpoint **sum** | **−2.9233%** | **+8.1461%** |
| Median index endpoint **sum + absolute difference** | −4.2844% | +6.5834% |
| Max(harbor) index endpoint **sum + absolute difference** | −3.0938% | +6.6920% |

A positive value means lower heldout MSE. The main original full17 analysis **does not support a broad positive increment** from historical port connectedness to later butterfly similarity under these controls. A **prespecified** sensitivity that omits geographically composite group9 DOES show a positive out-of-sample predictive increment, including **+8.1461% MSE / +5.4636% MAE** for the simpler historical median sum model. This is heterogeneity/instability, not proof of an ecological port-transport effect.

**Earlier 2001–2010 butterflies:** median historic index sum adds **+3.2351% MSE but −0.1984% MAE** across all17; omitting group9 and refitting gives **−0.5944% MSE and −1.8773% MAE**. The earlier target and port index partly overlap in time. Do not call this a 2001–2010 temporal forecast, or use it as clean proof of a real later historical network shift.

All original **20** contrasts, including the negative ones, MSE, MAE and per-pair better/worse counts, are saved in [comparison_scope_summary.csv](comparison_scope_summary.csv), with original full response and all group9-excluded responses preserved. These are correlated island pairs, NOT 136 independent biogeographic experiments or pair-independent p-values.

## Post-result diagnostic: the training membership is important, not merely test inclusion

The group9 sensitivity removed group9 from BOTH the training and the testing cohort. This could change the test distribution or fitted model. Therefore a separate, clearly **POST-RESULT** diagnostic ([protocol](POSTRESULT_GROUP9_DIAGNOSTIC_PROTOCOL.md)) reused the **same saved 2,048 predictions without fitting a single additional model**. It partitioned group9-containing test pairs and compared forecasts for the **same 120 original group9-free test pairs** while keeping original full17 versus no-group9 trained models separate:

| Later butterfly MSE gain, median historical PLSCIQ endpoint sum | Heldout test pairs | Model fitting population |
|---|---:|---|
| −2.9233% | Full original 136 pairs | Original 17 islands |
| **−1.8943%** | Same group9-free 120 test pairs | Original **17-island trained** models (group9 may appear in training) |
| **+8.1461%** | Exactly the SAME 120 test pairs | Models refit using **16 islands only** (group9 excluded from every training fold) |
| −8.2212% | The 16 test pairs involving group9 | Original 17-island trained models |

**The sign reversal cannot be explained solely by omitting group9's 16 pairs from the test set.** It is associated with *re-estimating models on an island-group-specific training population*. The full17-trained model still worsens prediction when evaluated on the no9-only 120 test pairs (−1.89%); using a training population that also excludes group9 gives +8.15%. This demonstrates sensitivity of the historical port index predictor to one composite source group under strict LO2 rather than proving that any country/island group causally suppresses maritime dispersal. Changes to the baseline fit also contribute: after excluding group9, BASE absolute MSE changes from **0.014181** on the original full17-trained no9 tests to **0.013170** for the refitted no9 models; thus both predictors and baseline are training-population dependent.

On the original 16 group9-containing test pairs, the historical index addition worsens MSE **−8.2212%**; 16 dependent pairs are too few for independent causal claims. The original group9 geography spans multiple countries and index-measured ports but that alone cannot be singled out as the verified mechanism. Historical PLSCIQ is an island-hub opportunity index not a measured realized movement or larval establishment driver.

All **40** diagnostic scores and MAE as well as MSE are preserved in [group9_training_vs_testing_loss_decomposition_40.csv](group9_training_vs_testing_loss_decomposition_40.csv), and the original 20 prior source-score rows were reproduced exactly. See [source-only diagnostic receipt](group9_diagnostic_receipt.json) and [independent group9 diagnostic audit](group9_independent_diagnostic_audit.json).

## Reproducible files, exact chronology and independent verification

- [Original 24 PortWatch geography-verified source ports and historical PLSCIQ qualification](original24_geoverified_ports_historical_PLSCIQ_quality.csv)
- [All 51 original landmass groups and historical port-index sampling availability](original51_historical_PLSCIQ_island_eligibility.csv)
- [Original 17 island pairs and both historic/later Jaccard responses, all 136 test identities](original17_136_source_pair_manifest.csv)
- [Every one of 256 source-only strict double-island-unseen train/test fold identities](double_unseen_source_folds.csv)
- [All 2,048 individual target/model/fold heldout predictions (compressed)](LO2_heldout_predictions.csv.gz)
- [The 20 prefit source contrasts including both later and early targets](comparison_scope_summary.csv)
- [Source/feature/hash execution receipt](execution_receipt.json)
- [Source-only independent audit of original 22 historical index ports, all 51 island IDs, 136/120 folds, 2,048 saved predictions and 20 MSE/MAE contrasts](independent_source_and_prediction_audit.json)
- [Post-result exact same-test-set decomposition: 40 source contrasts](group9_training_vs_testing_loss_decomposition_40.csv)
- [Independent post-result source forecast audit](group9_independent_diagnostic_audit.json)

The [source model](../../scripts/run_island17_historical_PLSCIQ_butterfly_LO2.py) and [separately implemented source audit](verify_LO2_source.py) reproduce original source port qualification, target years, every original island-pair test and all prediction loss values. Training data always omit both islands of the test pair, with exactly 105 or 91 other-island pairs per refit. [Independent source-run and 2048-prediction audit succeeded](https://github.com/zuizui0223/chocho/actions/runs/42437746577); the initial scientific pipeline encountered a syntax error in the independently written validator, not a failure of the fit, and the repaired validator independently passed without altering source observations or model coefficients. A separate [post-result source-only diagnostic and independent audit](https://github.com/zuizui0223/chocho/actions/runs/42438000209) succeeded on the original saved prediction rows; 40 diagnostic comparisons and prior 20 locked comparisons are identical when recomputed.

## What can be claimed and what cannot

**Data-backed finding:** chronological port-level connectedness information is not a robust universal source of later butterfly assemblage similarity across the original data-supported 17-island cohort. A prespecified sensitivity involving a highly composite original island grouping gives an appreciable *conditional* predictive improvement, but the full cohort does not; a post-result same-120-test decomposition shows that **island composition in model TRAINING is important**. This is useful in identifying which dimensions of geography and network source opportunity require better measurement, rather than selling a misleading uniformly positive maritime mechanism.

**Not supported:** proof of historical 2006–2010 bilateral port voyages or their 2011–2022 increase; verified shipborne seed/larval import; true breeding/colonization events; "Borneo blocks shipping-driven expansion"; independent island-pair p-values; global generalization from selected 17 source-available islands; a universal negative about maritime accessibility. The selection of source-qualified PortWatch/WPI ports may concentrate certain countries and settlement histories, and original GBIF butterfly/plant observation processes are incomplete.

The previous GEB global *potential* larval-host resource expansion (+54.9%), other original island/port candidate-host associations, earlier historical plant source, original later AIS spatial predictive gain and the strict original AIS local-750km node-shuffle specificity failure are all preserved in their separate original populations. **No GEB manuscript change or journal submission.**
