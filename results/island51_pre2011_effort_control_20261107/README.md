# Archived pre-2011 botanical recording effort and fixed island51 prediction

**Completed exploratory analysis (2026-11-07), stacked on PR #130.** All 51 original islands, 194 original port search buffers, 1,275 island pairs and both response windows are preserved. The source original host panel, 2001–2024 pooled candidate-host Jaccard, 182 marine-qualified ports and circuit graph are unchanged. The [locked protocol](PROTOCOL.md), [code](../../scripts/run_island51_pre2011_effort_control.py), [source table](../../provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_port_early_inventory_status.csv), [original historic audit](../../provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_port_early_inventory_audit.json) and [full reproducible CI](https://github.com/zuizui0223/chocho/actions/runs/42323799520) provide traceable evidence.

## Main later-butterfly result (first indexed 2011–2022)

| Contrasts on exactly 51 leave-one-island-out folds | MSE reference | MSE augmented | MSE improvement | Improved/worse | MAE improvement |
|---|---:|---:|---:|---:|---:|
| Prior candidate-host Jaccard beyond circuits + richness, CR → HCR | 0.00584999 | 0.00571656 | +2.2809% | 34/17 | +1.5127% |
| **Candidate-host Jaccard beyond circuits + richness + early botanical effort, CRE → HCRE** | **0.00578137** | **0.00566452** | **+2.0211%** | **32/19** | **+1.3619%** |
| Prior circuits beyond host identity + richness, HR → HCR | 0.00595497 | 0.00571656 | +4.0036% | 36/15 | +2.5097% |
| **Native AIS + wind circuits beyond host identity + richness + early botanical effort, HRE → HCRE** | **0.00584922** | **0.00566452** | **+3.1576%** | **33/18** | **+1.3841%** |
| Early botanical effort beyond existing joint model, HCR → HCRE | 0.00571656 | 0.00566452 | +0.9103% | 37/14 | +0.2814% |

**Positive predictive increments remain after the fixed historical botanical-recording proxy is included.** They are not independent tests of colonization or historical host availability. Reference errors differ, so changes in percentage points across rows are not a decomposition and percentages cannot simply be added.

The earlier-butterfly response, 2001–2010, is only contemporaneous context: with the same early botanical effort controls, host composition CRE → HCRE gives 4.9380% (33/51); MAE 2.7541% and circuit HRE → HCRE gives -0.4585% (29/51); MAE 0.5392% on MSE and MAE respectively. The latter MSE remains adverse, and both response windows should be reported without choosing the favorable one.

[All 14 comparison summaries](comparison_summary.csv) are published with 12-decimal log precision. Full precision scores, every fold, all 30,600 endpoint-specific predictions, hashes and the independent audit are in the [GitHub Actions artifact](https://github.com/zuizui0223/chocho/actions/runs/42323799520) (30-day retention).

## What the source actually contains

The original archived year-restricted GBIF plant query successfully covered **194/194 port buffers, 51/51 islands**. There are **141 port buffers with at least one indexed 2001–2010 plant record, 53 with zero indexed records**. By island group, **44/51 have at least one earlier plant record** and seven IDs (**12, 19, 21, 47, 134, 218, 242**) do not. These are observation gaps, **not botanical absence**. All original island assignments and source counts were checked against their archived Git blobs and original port totals.

Summing per-port buffer counts across all islands gives **155,002 counted plant-record occurrences**, but this is **not a unique GBIF record total** because nearby port buffers may overlap. The largest source group, landmass **12**, has 0 of these (0.00%) across 5 original ports. The new pairwise controls use log1p-transformed island sums (sum and absolute difference per island pair) rather than raw counts. This substantial effort imbalance is retained and disclosed, not removed through outcome-based island deletion.

The [per-island historical source coverage](historical_botanical_coverage_by_island.csv) has all 51 island groups and original port-buffer counts.

**Crucial temporal limitation:** the completed 2001–2010 *plant species-by-port* facet files, source-specific host/accepted-name mapping, and historic plant-identity network are in original local archived paths and are not in the public GitHub tree. The historical plant counts verified here are survey-exposure proxies. The predictor candidate-host Jaccard STILL pools **2001–2024**, and may postdate the 2011–2022 butterfly response. Consequently this is **NOT** a reconstruction of a pre-2011 host network or validation that resource identity existed before butterflies. The [separate archived earlier-host analysis](../../provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_dated_host_match_results_ja.md) found early compatible host records for 5/1,163 shared later butterfly pairs (different unit); the [independent early-recurrence study](../../provenance/reviewer_defenses/results/original230_dated_host_join_20261031/original230_prior_host_recurrence_results_ja.md) found 34/113 earlier compatible-host records among previously unselected butterfly first-index events, with a structured-null mean 25.2% versus observed 30.1%, but no added 5-year butterfly re-recording benefit. These earlier populations and estimands remain distinct and neither is superseded.

## Reproduction / audit

This addition is exploratory and specified after inspecting prior PR130 results, with new test definitions frozen before scoring. The pipeline reused 51 folds × 2 targets × 6 models, checked 306 previous C/H/HC original held-island MSE/MAE values and the previous six CR/HR/HCR global MSE means before accepting the new analysis.

The parent maximum absolute MSE discrepancy was **7.11e-17**, MAE discrepancy **2.08e-16**. The independent [standalone verifier](verify_predictions.py) checks **30,600 individual errors, 612 model folds, 14 comparisons, and two endpoint-specific predictions per source pair/model/target**. All five exact-commit workflow checks succeeded, including paper CI, anonymous review bundle, blinded GEB DOCX, and push/PR scientific reproductions.

Run from repository root:
```bash
python -m pip install -e ".[test,analysis]"
python -m pytest -q tests/test_run_island51_pre2011_effort_control.py
python scripts/run_island51_pre2011_effort_control.py --out-dir /tmp/island51_early_effort
python results/island51_pre2011_effort_control_20261107/verify_predictions.py --results-dir /tmp/island51_early_effort
```

No original scientific artifacts or the global GEB submission model were modified. These results are conditional predictive differences, not p-values or proof of larval feeding, movement routes, butterfly establishment, host quality, historical regional alien status or ecological causality.
