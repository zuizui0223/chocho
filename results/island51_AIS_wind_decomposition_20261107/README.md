# Separating AIS vessel-density and wind proxies on 1,275 two-unseen-island comparisons

## Result (exploratory, 2026-11-07)

The fixed **51 islands, 1,275 unique island pairs, and 15 target-independent longitude-block folds** from [PR132](https://github.com/zuizui0223/chocho/pull/132) were reused unchanged. Every pair was predicted once with **both endpoints absent from training**, with 235 same-longitude-block and 1,040 different-longitude-block pairs. Original candidate-host Jaccard, 2001–2010 all-plant recording effort, marginal host richness, geographic predictors, 194 original ports, and Circuitscape source matrices were retained unchanged.

**Both marine resistance fields contribute conditional prediction information for the later 2011–2022 butterfly-record communities**, including after conditioning on each other. Source-derived AIS density is not a ship route, and wind-speed resistance is not a directional migration path.

### Later butterflies: full fixed 2×2 source-factor comparisons

| Added field | All 1,275 pairs MSE gain | Within-block 235 pairs | Between-block 1,040 pairs | All-pair MAE gain |
|---|---:|---:|---:|---:|
| AIS added without wind | +2.2226% | +0.3529% | +3.4559% | +1.1166% |
| Wind added without AIS | +1.3593% | +2.7155% | +0.4646% | +1.1073% |
| Wind added beyond AIS | +2.3030% | +3.7256% | +1.3344% | +1.5196% |
| AIS added beyond wind | +3.1581% | +1.3876% | +4.2996% | +1.5289% |
| Both added | +4.4744% | +4.0654% | +4.7442% | +2.6192% |

The largest separate between-block conditional gain is **AIS beyond wind: +4.2996% MSE**, and wind beyond AIS remains **+1.3344% MSE**. This supports **distinct additional predictive information** from these two represented physical layers under the frozen additive model. It does not prove either source transports butterflies or host plants. Some MAE contrasts differ in sign: AIS added alone to within-longitude-block pairs worsens MAE by 1.2271%, and AIS conditional on wind in the same within-block class worsens MAE by 0.3412%. All values are retained in the [complete 30-row table](contrast_scope_summary.csv).

These gains have different reference-model denominators. AIS-only gain +2.2226% and wind added conditional on AIS +2.3030% do **not** sum directly to the combined +4.4744% because each percentage has a different baseline.

### Earlier butterflies: contradictory temporal context must remain visible

The 2001–2010 first-indexed butterfly comparison is explicitly reported although it is unfavorable for the recent marine circuit fields. The same model comparisons, using the same modern AIS/wind sources, give:

| Added field | All-pair MSE change (positive means improvement) | All-pair MAE change |
|---|---:|---:|
| AIS added without wind | -6.8798% | -3.1618% |
| Wind added without AIS | -3.4356% | -1.3310% |
| Wind added beyond AIS | -3.1759% | -1.1412% |
| AIS added beyond wind | -6.6115% | -2.9686% |
| Both added | -10.2742% | -4.3391% |

All five MSE improvements are **negative** for earlier butterflies. Modern commercial AIS (2015–2021) and the associated wind period do not reconstruct actual 2001–2010 marine movement. These descriptive differences across response windows **are not** evidence that real historical vessel use, dispersal or biological mechanisms changed over time.

## Original source, fit and sample

The original 51 island IDs and all source butterfly Jaccard observations are exactly those of the earlier [observed circuit analysis](../island51_observed_ocean_circuit_20261106/README.md) and [joint host/circuit analysis](../island51_host_circuit_joint_20261106/README.md). The already-audited original longitude blocks and train/test partitions are inherited from [PR132](../island51_double_unseen_blocks_20261107/README.md); no source-driven grouping was selected using butterfly model outcomes.

All four models retain the original **HRE baseline**: geographic distance and original port-count variables, uniform ocean resistance, candidate-host Jaccard, candidate-host marginal richness terms, 2001–2010 general botanical port-buffer observation-count controls, and—for later butterflies only—earlier butterfly similarity. Extra features:

- HRE: native AIS absent, wind absent
- HRE_AIS: add native AIS resistance only
- HRE_WIND: add wind resistance only
- HCRE: add both

Two original Circuitscape conductivity fields are the previously fixed models, not new AIS/wind calibrations. All fits use existing train-only-standardized OLS with intercept and predictions clipped to [0,1]. The source is **native AIS activity density and mean wind-speed fields**, not route-specific or directional transport. Source labels are defined explicitly in the [frozen protocol](PROTOCOL.md). The prior HRE→HCRE result was independently reproduced for **both responses and all three scope categories** before accepting the five new contrasts.

Full reproducibility: [read-only script](../../scripts/run_island51_AIS_wind_decomposition.py), [independent stdlib verifier](verify_AIS_wind_decomposition.py), [complete summaries](contrast_scope_summary.csv), and [full 10,200 source-pair predictions / 120 model-fold / 150 comparison-fold artifacts](https://github.com/zuizui0223/chocho/actions/runs/42325732617) (30-day artifact retention). The source-based verifier reconstructs original butterfly observations, 15 geography folds, no overlap into any training island, every numeric error, and every saved contrast.

The existing 51-island LOIO and two-unseen-endpoint results are exploratory retrospectively inspected analyses; this new 2×2 study was specified **after** the combined +4.47% result was known. It has no new external validation population, no naive fold/pair p-values or CIs, and no causally identified transport mechanism. **2001–2024 pooled host Jaccard can postdate butterflies**, and historical plant sampling controls do not reconstruct an actual pre-2011 host-identity panel. Actual plant introduction, larval feeding, local reproduction, ship passage and butterfly establishment remain unobserved. Earlier positive host associations and GEB +54.9% resource-opportunity estimates remain distinct.
