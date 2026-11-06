# Single native-AIS source-quality sensitivity

This is a **post hoc source-quality sensitivity**, chosen after the primary overview-based biological results had been seen. The separate [plan](../NATIVE_SOURCE_SENSITIVITY_PLAN.md) and [lock](../native_source_sensitivity_plan_lock.json) were frozen at **2026-11-06 08:54:25 UTC**, before fitting any native-AIS outcome model. The original protocol, source manifest, matrices, scores and six official numerical checks remain unchanged.

The source audit found nearest-pixel subsampling in the first overview. This one check replaces only the AIS 0.5° input with an area-overlap mean of the original 0.005° native pixels. There was no conductivity tuning, model sweep or further geographic selection. The 51 islands, 1,275 pairs, 182 marine terminals, original port-count controls, uniform/wind matrices, transformations and folds were retained.

## Results alongside the primary analysis

Positive values below mean lower held-out MSE relative to the identical base-plus-uniform model. The two columns differ only in AIS source averaging.

| Comparison | Target | Primary overview | Native source | Change, percentage points |
|---|---|---:|---:|---:|
| Additionally AIS | Candidate host similarity | −0.1522% | −0.1056% | +0.0466 |
| Additionally AIS | Later butterfly similarity | +0.5582% | +0.6595% | +0.1013 |
| Additionally AIS | Earlier butterfly similarity | −0.0395% | −0.1248% | −0.0853 |
| Additionally AIS and wind | Candidate host similarity | +0.3211% | +0.3517% | +0.0305 |
| Additionally AIS and wind | Later butterfly similarity | +6.9802% | +7.0065% | +0.0263 |
| Additionally AIS and wind | Earlier butterfly similarity | −1.1770% | −1.2536% | −0.0766 |

The combined later-butterfly comparison improves **40/51** island holdouts in both versions. Its MSE changes from **0.0060019218** to **0.0060002266**; its MAE becomes slightly worse, from **0.0541204888** to **0.0541318508**. The earlier-butterfly AIS-containing results also deteriorate, as shown above. The full side-by-side table retains both metrics and all targets.

The uniform and wind-only comparisons are **exactly unchanged**, including every individual fold loss and identity: the maximum absolute numeric difference is **0**. In particular, the later-butterfly wind-only MSE reduction remains **5.4961%**. Six fixed control comparisons across the two unaffected model comparisons and three outcomes were checked recursively, rather than comparing only headline scores.

The source-quality change does not reverse the declared predictive pattern in this single check: modern AIS adds little for later butterflies and does not improve host MSE on its own, while the wind-containing models retain their larger later-butterfly improvement. This is not an independent validation set, a demonstration of physical transport, or a reason to replace the preserved primary analysis with whichever score is larger.

## Source and numerical gates

Before fitting, the runner independently verified:

- The frozen sensitivity plan and lock, plus the original source manifest and protocol.
- Native grid SHA256 **`cbb3c963bd4f466453937f5792f2bbd47786fb7f380d38a4b128ab86295a475f`**, its derivation receipt and the unchanged 39-cell source audit.
- Original ZIP identity, the native TIFF member's exact CRC and byte length, native 0.005° geometry, scale 1, offset 0, and absence of detected overview factors.
- **All 39** independently calculated native pixel-overlap means; maximum relative error **5.10 × 10⁻⁸**, below the fixed `3e-7` float32 tolerance.
- Original51 source recovery and all original network export hashes, by rerunning the primary `prepare_inputs` in a separate temporary directory.
- Byte-identical uniform and wind network exports and terminal weights. Every graph cell remains observed; no missingness-driven domain change is allowed.

The new AIS matrix has a maximum linear-system relative residual of **1.24 × 10⁻¹³**. Across all 1,275 pairs, its correlation with the primary raw AIS resistance is **0.9998381**, its median native/primary resistance ratio is **0.9851161**, and its largest absolute relative pair change is **5.5529%**. These are numerical comparisons of the two represented networks, not independent source calibration.

Actual **Julia 1.13.3 / Circuitscape 6.1.0** then checked both fixed pairs on the new native-AIS network. The maximum relative resistance difference from SciPy is **3.97 × 10⁻¹¹**, and the maximum all-node Kirchhoff relative residual is **4.37 × 10⁻¹⁰**. The unchanged portable six-case CLI also repeated four uniform/wind controls. This new [receipt](native_AIS_official_validation.json) separates the **two new native-AIS checks** from those four controls and from the original six-case evidence.

## Reproduce the one sensitivity

From the repository root, with the declared Python analysis dependencies installed:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python scripts/run_island51_native_AIS_source_sensitivity.py \
  --input-root data/exploratory/island51_observed_ocean_20261106 \
  --primary-results-dir results/island51_observed_ocean_circuit_20261106 \
  --native-grid data/exploratory/island51_observed_ocean_20261106/ais_native_quality_check/commercial_AIS_native_mean_0p5deg_2015_2021.npz \
  --native-receipt data/exploratory/island51_observed_ocean_20261106/ais_native_quality_check/commercial_AIS_native_mean_0p5deg_2015_2021_receipt.json \
  --native-audit results/island51_observed_ocean_circuit_20261106/ais_source_quality/published_AIS_source_quality_receipt.json \
  --out-dir /tmp/chocho-native-source-sensitivity \
  --network-dir /tmp/chocho-native-quality-networks
```

The frozen native NPZ, receipt and independent audit are sufficient; the approximately 9.8 GB native TIFF is only needed to rebuild the source derivative. `--plan` and `--plan-lock` default to the preserved files beside this result directory. Primary outputs can also come from a fresh run in a different directory, provided its verified source identities and fixed design match.

To repeat the separate official checks using the pinned environment described in [official validation instructions](../official_validation/README.md):

```bash
"$JULIA_BINARY" --project="$CS_CHECKOUT" \
  scripts/verification/validate_island51_observed_ocean_circuits.jl \
  /tmp/chocho-native-quality-networks \
  /tmp/chocho-native-source-sensitivity \
  /tmp/chocho-native-source-official-checks
```

The sensitivity-specific tests are `tests/test_run_island51_native_AIS_source_sensitivity.py`: six tests cover altered plans, source checksums/CRC/overview use, all39 independent means, control fold values and control fold identities. The repository suite after this addition reports **189 passed, 1 skipped**; the existing optional upstream wind-builder test requires `xarray`.

## Outputs and unchanged interpretation limits

| File | Content |
|---|---|
| [native_source_preflight.json](native_source_preflight.json) | Pre-fit gate timestamp, exact plan/source hashes, all39 native overlap comparisons and unchanged geometry/network evidence |
| `native_AIS_distributed_port_resistance.npz` | New51×51 distributed-port AIS matrix |
| [island51_observed_circuit_pairs.csv](island51_observed_circuit_pairs.csv) | All1,275 pairs, native AIS predictor, unchanged uniform/wind predictors and original AIS resistance |
| [model_comparison_summary.csv](model_comparison_summary.csv) | All12 target/comparison summaries |
| [all_island_holdout_losses.csv](all_island_holdout_losses.csv) | All612 fold records, including adverse results |
| [primary_vs_native_comparison.csv](primary_vs_native_comparison.csv) | Direct primary/native MSE, MAE, reduction and island-count comparison |
| [native_source_sensitivity.json](native_source_sensitivity.json) | Complete comparison definitions, invariants, diagnostics and protected-primary-file hashes |
| [native_AIS_official_validation.json](native_AIS_official_validation.json) | Actual new native-network verification scope and errors |
| `official_validation/` | Small regenerated source/ground/configuration files, numerical checks and successful log |

The circuit remains an undirected distributed-port energy measure. AIS is a modern published observation-density index, not measured bilateral cargo or historical plant carriage. Wind is an undirected monthly10m proxy; source quality sensitivity does not validate a dispersal mechanism. Pairs and folds are dependent, observed first records do not establish colonization dates, and the plants can have been recorded after the butterfly period. No additional sensitivity or parameter search was performed after this check.
