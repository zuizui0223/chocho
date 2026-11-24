# Original155 geographic distance precision — controlled pyproj test (24 November 2026)

**Completed experiment, negative version comparison:** [CI #45323737296](https://github.com/zuizui0223/chocho/actions/runs/45323737296), independent private original-data source inputs, **two pinned pyproj releases**, fixed sklearn 1.7.2 and NumPy 2.3.5, SHA-verified original ZIP and original frozen panel oracle. Each controlled arm independently reconstructed the original 155-species geographical domain of **280,988 species×cells**, without copying archived original distance values into its calculations. The reference was the period-0 distance values in the verified archived panel. Each species-cell distance is repeated in three model-panel periods (842,964 rows).

| Comparison against original `distance_km` | pyproj **3.7.1** | pyproj **3.7.2** |
|---|---:|---:|
| Internal PROJ version | **9.5.1** | **9.5.1** |
| Matching species×cell key set / physical ordering | Yes / Yes | Yes / Yes |
| Exact float64 differences | **20,213 / 280,988** | **20,213 / 280,988** |
| Raw float32 differences | **5,361 / 280,988** | **5,361 / 280,988** |
| **Actual RF `log1p(distance_km)` float32 differences** | **5,361 / 280,988** | **5,361 / 280,988** |
| Maximum absolute distance discrepancy | **1.4495071809506044e−12 km** | **1.4495071809506044e−12 km** |

These counts multiplied over all three periods correspond to **60,639** exact distance discrepancies and **16,083** post-float32 discrepancies observed in the preceding full-panel audit. The significant learning sensitivity was independently demonstrated in the [single-distance-column intervention #45322533520](https://github.com/zuizui0223/chocho/actions/runs/45322533520): it restored the original **+4.120175607931824%** MSE improvement, **175/175** prediction/score groups and **zero** out-of-tolerance M0–M5 predictions.

**Interpretation:** Changing the *tested* pyproj minor versions from 3.7.2 to 3.7.1 while both have PROJ 9.5.1 did **not** eliminate the tiny geographic distance precision differences. It is not justified to claim which historical library or CPU routine produced the archived values; the historical environment was not fully pinned. The all-96-column numeric reproduction, temporary one-column original-distance counterfactual and source-independent exact numeric reconstruction have **different evidentiary scopes**.

This audit does not examine alternate PROJ builds, original source recording-day construction from GBIF raw records, the creation lineage of three late-recovered CSVs, unseen-species fold 3 or two remaining sensitivity refit differences. Do not modify frozen results to force a fit, call the end-to-end chain independently reproduced or claim reviewer redistribution rights. Both PRs remain Draft; no source rows exported.

[Aggregate source receipt](../../../results/geb_v05_private_source_recovery_20261120/ORIGINAL155_GEOMETRIC_PYPROJ_VERSIONS_20261124.json) · [CI artifact](https://github.com/zuizui0223/chocho/actions/runs/45323737296/artifacts/14888288922).
