# Candidate-host composition and marine circuits across pairs of two unseen islands

## Results (2026-11-07, exploratory double-unseen geographic holdout)

**All 1,275 original unordered island pairs were predicted once with BOTH endpoints excluded from their respective training data.** The fixed 51 islands, original candidate host Jaccard, native AIS and wind, marginal host richness and pre-2011 all-plant recording-effort proxy remain unchanged. Longitude-only grouping used five blocks of 11,10,10,10,10 islands, without inspecting response scores. The original prior analyses and positive results remain unchanged; their leave-one-island-out estimand differs from this new double-unseen task.

### Later butterfly first-index similarity (2011–2022)

| Added term in the fixed three-model comparison | All 1,275 pairs: MSE reduction | Within-block 235 pairs | Between-block 1,040 pairs | All-pair MAE reduction |
|---|---:|---:|---:|---:|
| **Native AIS + wind circuit fields** (HRE → HCRE) | **+4.4744%** | +4.0654% | **+4.7442%** | +2.6192% |
| Candidate-host Jaccard (CRE → HCRE) | +0.8808% | +1.9729% | +0.1417% | +0.2978% |

The circuit increment is positive under the more demanding evaluation, including pairs from different original longitude blocks. The candidate-host identity increment remains weakly positive overall but is only +0.1417% MSE for between-block pairs; **that between-block host MAE increment is negative (−0.3255%)**. This means that the regional transfer of candidate host identity is not consistent across error metrics. The common joint model's later MSE over all 1,275 held-out pairs is 0.007107995; compare with 0.007440932 for the circuit-absent model and 0.007171157 for the host-absent model.

Earlier first-indexed butterflies (2001–2010) are a secondary, contemporaneous contrast: added host Jaccard gives **+2.9550%** MSE, while AIS+wind gives **-10.2742%** MSE and −4.3391% MAE. The earlier adverse circuit result is prominent; it must not be deleted or treated as direct evidence that shipping and wind mechanisms changed in time.

### What was held out

The source-backed [island block map](longitude_island_blocks.csv) retains every representative longitude, latitude, landmass ID and original example port. Islands were sorted by original (pre-existing) representative longitude, breaking ties by numeric island ID, then divided into 5 fixed blocks. The [protocol](PROTOCOL.md) freezes this rule before the current fit. These are **longitude strata**, not true biogeographic regions. Some blocks span large longitudinal distances because the 51-island input population is uneven.

For each unordered pair of longitude blocks (15 folds: five same-block and ten cross-block), the training set contained only pairs among islands NOT in either held-out block. This yielded training sizes of **435–820 pairs**, test sizes 45–110 pairs, and collectively **235 within-block and 1,040 between-block predictions per model and target**. All source pairs were evaluated once, unlike the prior 51-island leave-one-island-out protocol that generated two endpoint-specific predictions for each pair and trained with the *other* endpoint represented. Therefore the two percentages must not be interpreted as changes on a fixed denominator/training set.

The model family is identical to the previous PR131:
- CRE contains geographic and original port-count controls, uniform ocean resistance, native AIS and wind resistances, candidate-host marginal richness and pre-2011 general botanical port-buffer record effort.
- HRE has the same controls and original candidate-host Jaccard in place of native AIS/wind.
- HCRE includes both host Jaccard and native AIS/wind with all the same controls.

The later response also includes early (2001–2010) butterfly Jaccard as a baseline predictor, while the early response does not include its own values as a predictor. OLS uses train-only standardization, intercept and [0,1] clipped predictions. The [saved scope summaries](contrast_scope_summary.csv) list all 12 target / comparison / geographic-scope results, including the MSE/MAE signs and per-pair improvement counts. Model-fold and contrast-fold scores, full 7,650 predictions, geographic/source check receipts and independent audit are retained in the [GitHub Actions artifact](https://github.com/zuizui0223/chocho/actions/runs/42325176725) (30-day retention).

### Independent verification and interpretation limits

An independently implemented standard-library [audit](verify_double_unseen.py) reconstructs all **1,275 test pairs, 7,650 predictions, 90 model-fold scores, 60 comparison-fold scores and 12 scope summaries** from original source values. Every excluded training block is reconstituted from the source; no test-pair endpoint enters the training population for its fold. The audit rejects changed butterfly observations, source pair IDs, fold assignment, train/test counts and numerical errors. Dedicated CI, manuscript paper CI, anonymous snapshot and blinded document build passed for the initial validated results. The final source/published-table verification is additionally enforced by CI.

**Main limitation:** This reuses the same retrospectively selected 51 islands, so it is a stronger holdout *task* rather than an independent new population. Host candidates are pooled across 2001–2024, potentially later than the 2011–2022 butterfly observations; an actual pre-2011 host-species network is still not reconstructed because the raw historical species×port files are absent from this GitHub repo. Adding historical botanical record-count controls does not fix this data-date mismatch. AIS and undirected wind circuits are connectivity proxies, not direct shipping or directional dispersal tracks. First-indexed butterfly similarity is not arrival, successful larval feeding, population growth or colonization. The longitude-block folds and test pairs are not independent ecological replicates; no p-values or naive intervals are estimated. Previous global GEB +54.9% host opportunity expansion and PR128–#131 island results retain their original population/estimands.
