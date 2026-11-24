# Original155 time-split RF: isolate source-panel distance precision (24 November 2026)

**Main new result: causal numeric discrepancy isolated, not yet independently source-bit reproducible.** The original panel builder reconstructed **842,964 rows and all 96 columns** within its declared numeric tolerance. First temporal RF fit from that generated Parquet produced M5-over-M3 MSE improvement **+3.9356441905%** rather than the archived **+4.1201756079%** and failed 66 of 175 prediction/score groups ([CI #45320236263](https://github.com/zuizui0223/chocho/actions/runs/45320236263)).

The second source-bound comparison found all **842,964 row keys** in the exact same physical row positions, including species-first-appearance order. Exactly one numeric source-panel column differed in *exact* comparison: `distance_km` (60,639 of 842,964 rows; 16,083 differences remain after float32 casting). The 95 other panel columns have no exact numeric differences. So the original row-order hypothesis is **rejected** ([CI #45321468031](https://github.com/zuizui0223/chocho/actions/runs/45321468031)).

## Predeclared one-column counterfactual

The [successful source-bound CI #45322533520](https://github.com/zuizui0223/chocho/actions/runs/45322533520) ran unchanged original scripts, the same 13 recovered builder inputs, and the same numerical comparison. Before fitting, it **copied only the original archived `distance_km` column** into the newly rebuilt 842,964-row panel; the other **95 columns remained from the reconstructed sources**. This is a diagnostic intervention, not an independent feature reconstruction.

| Controlled check | Result |
|---|---|
| Original model test support | **208,482 keyed rows** |
| Original-vs-counterfactual numeric comparison groups | **175/175 PASS** |
| Per-model predictions above fixed tolerance | **0 for all 6 RF models + persistence** |
| M5-vs-M3 MSE reduction from counterfactual refit | **+4.120175607931824%** |
| Original archived reduction | **+4.120175607931824%** |
| Numerical comparison | **RTOL 1e−6, ATOL 1e−9** |
| Archived predictor columns reused | **only `distance_km`** |

**Inference:** The tiny differences in `distance_km` introduced during the source-conditional panel reconstruction are the cause of the RF numerical prediction discrepancy; the mechanism by which the geographic distance values acquired those tiny differences is not yet established. The unmodified builder derives distances from original past occurrence cell centres using `pyproj.Transformer.from_crs(6933, 4326)` and `sklearn.neighbors.BallTree` with haversine distance. A software environment/geometry numeric effect is plausible but **untested for the historical original**. A small difference can affect a fixed-seed tree-fitting algorithm, but this result does not isolate the exact library responsible.

The counterfactual **does not** establish that the 13 originally available upstream sources can reproduce the historically serialized distance values bit-for-bit; it explicitly reused one original source variable. The raw occurrence partitions were not processed through the original recording-day aggregation, the three 24-Nov recovered CSVs' historical creation lineage remains unverified, and source-licence/anonymous-data repository gates remain open. The original unseen-species fold 3 and two sensitivity numerical mismatches are independent questions.

This report was added to **Draft PR #189**, leaving paper **Draft PR #179** and its frozen scientific estimates unchanged. No private rows or models were uploaded; the audit artifact is aggregated only, not a reviewer deposit. [Machine-readable receipt](../../../results/geb_v05_private_source_recovery_20261120/ORIGINAL155_DISTANCE_ONLY_COUNTERFACTUAL_PASS_20261124.json).
