# GEB original155: row-order hypothesis rejected, distance_km float32 differences localized (24 November 2026)

**Directly verified:** [Ordered-panel and model rerun #45321468031](https://github.com/zuizui0223/chocho/actions/runs/45321468031) found **842,964/842,964 exact source→rebuild row-position key matches**. The original saved panel and the 96-column rebuilt panel also have identical species first-appearance order. Reordering the rebuilt panel to the original row-key order did not change the failed temporal RF reproduction: **+3.935644%** fresh versus **+4.120176%** archived M5-vs-M3 MSE reduction, with 66 of 175 checked score/prediction groups outside fixed numerical tolerance. The initial row-order causal hypothesis is therefore **rejected**.

Sorting the original and regenerated tables by independent `speciesKey,cell,period` keys and comparing **exact** floating-point values identified exactly one discrepant column: `distance_km`.

| Differential | Rows (of 842,964) |
|---|---:|
| Original vs rebuilt physical row order differs | **0** |
| Numeric `distance_km` values differ exactly (float64 representation) | **60,639** |
| `distance_km` values still differ after model's float32 cast | **16,083** |
| All other original model-panel numeric columns with any exact/float32 mismatch | **0** |

The previous all-96-column PASS used a declared panel numeric tolerance (`rtol=1e−6`, `atol=1e−8`). That value-consistency result was not wrong, but did not imply **exact float32 inputs** to the random forest. Different `distance_km` float32 values are a concrete, sufficient *candidate* for the extensive forest prediction differences. **Causal attribution is not yet confirmed**: a [single-column counterfactual #45322533520](https://github.com/zuizui0223/chocho/actions/runs/45322533520) is now fitting the original RF from the rebuilt panel with only the original `distance_km` column substituted. It does not reuse the original target or the remaining 95 predictors. The counterfactual's use of an original source predictor prevents labeling it independently source-reproducible even if all model predictions match.

- Distinct initial numeric mismatch receipt: [run #45320236263](https://github.com/zuizui0223/chocho/actions/runs/45320236263).
- Confirmed all-column reconstruction: [run #45318106188](https://github.com/zuizui0223/chocho/actions/runs/45318106188).
- Original 13-input sources remain checksum-verified and private; full raw GBIF occurrence→recording-day source preprocessing and acquisition chain of the newly recovered CSVs remain open.
- All original scientific paper results are unchanged. Original unseen-species fold-3 and two sensitivity numerical discrepancies are distinct open issues. Neither submission nor data release is approved.

