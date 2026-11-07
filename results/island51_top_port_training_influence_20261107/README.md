# Training influence of the two port-heaviest original island groups

## Main new result

**In the fixed original 51-island double-unseen longitude-block test, AIS-derived prediction improvements remain positive after removing source pairs associated with either or both of the two largest original port groups from the training data.** The original latitude-adjusted combination of AIS plus wind remains adverse in every training scenario.

This completed *post-hoc exploratory* test is stacked on [PR134](https://github.com/zuizui0223/chocho/pull/134). Its frozen [protocol](PROTOCOL.md) specified the four training-population perturbations before the current new refits, but all source and earlier predictive results had already been seen. No independent external island population or causal effect is claimed.

The original source-only port-count roster identifies landmass ID **15 with 39 original botanical ports** and ID **9 with 21**, together accounting for 60 of the 194 original biological ports. Port grouping labels are numeric landmass IDs, not authoritative named islands or measured real shipping routes. [Original complete 51-island ranked roster](original_port_roster_rank.csv) retains the basis of the selection.

### Later butterflies: 2011–2022, same 1,275 held-out island pairs in all four scenarios

| Training source population | AIS alone MSE improvement | AIS alone MAE improvement | AIS beyond wind MSE improvement | Wind alone MSE improvement | AIS+wind combined MSE improvement |
|---|---:|---:|---:|---:|---:|
| Full original training | +1.0020% | +0.8485% | +0.5446% | -4.6618% | -4.0918% |
| Exclude landmass 15 (39 ports) | +1.3953% | +1.0754% | +1.1465% | -4.5053% | -3.3071% |
| Exclude landmass 9 (21 ports) | +0.8609% | +0.6576% | +0.3366% | -4.7081% | -4.3556% |
| Exclude both port-heavy landmasses (60 ports) | +1.2959% | +0.8564% | +0.9848% | -4.5755% | -3.5456% |

All four AIS-alone MSE increments are positive, **+0.8609% to +1.3953%**. All four AIS-alone MAE increments are also positive. They describe conditional out-of-sample predictive information, **not the causal effect of ports**. Wind-alone and joint AIS+wind MSE changes are adverse throughout the same fixed conditions. AIS added after wind is positive but has a different reference-model error denominator from AIS alone.

### The difference between conditional improvement and model prediction accuracy

Reducing the number of training islands changes absolute forecast difficulty. The prediction error in this table is a model-specific MSE, not a percentage comparison:

| Training source population | Latitude-controlled baseline HRE MSE | Baseline + AIS MSE | Baseline + AIS + wind MSE |
|---|---:|---:|---:|
| Full original training | 0.007271385 | 0.007198528 | 0.007568916 |
| Exclude landmass 15 (39 ports) | 0.007377632 | 0.007274692 | 0.007621617 |
| Exclude landmass 9 (21 ports) | 0.007382636 | 0.007319081 | 0.007704195 |
| Exclude both port-heavy landmasses (60 ports) | 0.007506698 | 0.007409418 | 0.007772854 |

Removing source observations can worsen the absolute MSE even as AIS retains a positive conditional improvement. Therefore these results should NOT be narrated as "excluding the large port groups makes predictions better." The fixed training exclusion diagnoses the influence of source data; it does not estimate a local biological mechanism.

### Geographic scope of the AIS-only comparison

| Training source population | All 1,275 pairs | Within longitude block (235 pairs) | Between longitude blocks (1,040 pairs) |
|---|---:|---:|---:|
| Full original training | +1.0020% | +0.3812% | +1.3347% |
| Exclude landmass 15 (39 ports) | +1.3953% | +0.4355% | +1.8973% |
| Exclude landmass 9 (21 ports) | +0.8609% | +0.3446% | +1.1345% |
| Exclude both port-heavy landmasses (60 ports) | +1.2959% | +0.3738% | +1.7705% |

These are target-independent *longitude-defined* strata, not independent biogeographic regions. Full-sample ratios weight the 1,275 individual squared errors and are not simple averages of the within-/between-block percentages.

### Earlier butterfly records: 2001–2010 context (unfavorable findings retained)

| Training source population | AIS alone MSE improvement | AIS + wind combined MSE improvement |
|---|---:|---:|
| Full original training | -8.1571% | -7.7788% |
| Exclude landmass 15 (39 ports) | -7.7547% | -7.7755% |
| Exclude landmass 9 (21 ports) | -7.0273% | -6.2864% |
| Exclude both port-heavy landmasses (60 ports) | -6.6261% | -6.2483% |

These earlier outcomes were observed before or during different marine/environmental conditions than the source AIS activity density (2015–2021). An adverse predictive score on older indexed communities is not a verified negative effect of actual ship traffic or evidence that biological transport changed in time.

## Frozen design and verification

The exact original **51 eligible landmasses, 194 botanical ports, 1,275 unique unordered island pairs**, original 225-butterfly host-candidate panel, source-derived 15 longitude-block pair folds, 182 eligible marine ports and three frozen uniform/native-AIS/wind resistance matrices remain unchanged. All models retain original great-circle distance, original port-count controls, candidate host similarity pooled 2001–2024, candidate-host marginal richness, early 2001–2010 *all-plant recording effort* and the four original source-latitude/hemisphere controls. Later models also control 2001–2010 butterfly composition.

For each of the 15 folds, every source pair involving either held-out longitude block stays excluded from training. The four fixed scenarios additionally remove training pairs involving ID15, ID9 or both. **The same 1,275 source pairs are evaluated once per scenario and per model, regardless of source exclusions.** Training never contains a test pair, either held-out longitude-block endpoint or a scenario-excluded source island, and can contain fewer original pairs when the two groups are not already held out.

There are precisely **40,800 saved model-specific source-pair predictions**, **480 model folds**, **600 contrast/fold losses**, **120 source-factor/scope summaries**, plus 60 source-specified scenario/fold training ledgers. The original 30 latitude-adjusted PR134 MSE/MAE comparison rows reproduce within **4.790e-13** absolute error against the publicly pinned parent CSV.

An independently implemented Python-standard-library auditor reconstructed **40,800 individual prediction errors**, all 480 model-fold losses, 600 fold contrasts and 120 comparison summaries. It also independently checks original port ranks, all 60 train/test support scenarios, the absence of every excluded island from training, original butterfly source-pair identities and complete identical test-pair populations. Source and prior results remain unchanged. The [full CI reproduction](https://github.com/zuizui0223/chocho/actions/runs/42332615939) archives the full 17-digit numerical outputs, source manifests and independent audit (GitHub retention 30 days).

Permanent sources: [fixed protocol](PROTOCOL.md), [read-only runner](../../scripts/run_island51_top_port_training_influence.py), [standard-library auditor](verify_port_influence.py), [120-row comparison summary](contrast_scope_summary.csv), [source port roster](original_port_roster_rank.csv).

## Inference limits

The strongest supported inference is **source-training-population robustness of AIS's modest incremental predictive information after coarse latitude control**. Training exclusions do not prove movement via ports, an environmental or resource-mediated causal effect, independence from all regional composition, or that candidate hosts in 2001–2024 existed before butterfly observations in 2011–2022. AIS is a vessel-density proxy, not actual directional shipping routes; wind is an undirected mean-speed circuit, not an air-current dispersal trajectory. Local butterfly reproductive success, larval feeding and colonization remain unobserved. Source 2001–2010 botanical record-count controls are not 2001–2010 species-identity data. All four scenarios evaluate the SAME original 51 islands, with correlated pairs and overlapping training folds: these are NOT independent region-level p-value replicates. Prior positive GEB global candidate-host resource expansion and original island–port findings retain their own populations and estimands.
