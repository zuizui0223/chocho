# GEB v0.5 — original 155-species source recovery and refit QA

**22 November 2026. Temporal and spatial numerical reproduction passed; the unseen-species numerical comparison failed after successful model execution. Stage A remains open.**

## Source identity and scope

The existing private draft release supplied the original analysis archive; no replacement cohort was substituted. Four split files reconstructed a 5,552,304,553-byte ZIP with SHA-256 `9e39b0068c3617741f9b275a1edfddd574c93d620853e956bd706ad067a53f9f`. The checksum-verified manifest contained 3,032 source records and the ZIP 3,033 members. Selected required members were individually verified before execution; the receipt does not claim an independent extraction or rerun of every archive member.

The model panel has **842,964 rows and 155 species**, SHA-256 `5ae5b48b752fad05067537de5fb95d4ee74f17f3ba1114443d337bd7fd3695a1`. The model lock has SHA-256 `0c7fc34d7858d6f62fe31a914188aee63e56ab1d41872abe40c314f6c5ed8b1e`. The fit script at original commit `42f8fcc3a388e96456c76180d7e56d75dd2ed48c` has SHA-256 `a34ad4cd1e4f40460c4f6bf6ed0b22ba3885b17d20ac17f74bd6564d0618ce6d`. The recording-day input retains `datasetKey`.

The complete transitions contain 193,703 evaluable training rows and 208,482 test rows before mode-specific fold exclusion. All refits retain 160 trees, depth 12, minimum leaf 40, feature fraction 0.7, seed 20261029 and four jobs. Relative tolerance `1e-6` and absolute tolerance `1e-9` were declared before comparison and have not been relaxed.

## Actual refit results

[Run 44970180129](https://github.com/zuizui0223/chocho/actions/runs/44970180129), at source-audit code head `1483144a2b8d89b992635dfe6507a088a0ab7d42`, completed fresh execution in all three jobs. Its overall failure records the species numerical-comparison failure; it is not a failed archive recovery or an unexecuted fit.

| Evaluation | Original M5-versus-M3 MSE reduction | Fresh result | Numerical comparison |
|---|---:|---:|---|
| Temporal | +4.120175607932% | +4.120175607932% | PASS: 144/144 comparison groups; maximum RF prediction difference `5.551115123125783e-16` |
| Spatial block | +2.995205951699% | +2.995205951699% | PASS: 174/174 groups including the ≥100-km checks; maximum RF prediction difference `4.440892098500626e-16` |
| Spatial ≥100 km | +2.735657911003% | +2.735657911003% | PASS: identical membership, 102,728 rows / 1,788 cells / 155 species; same spatial fits, not a fourth model |
| Unseen species | −0.208520153993% | −0.208933309709% | FAIL: 68/144 groups pass and 76 fail; adverse direction retained |

Each comparison includes keyed predictions, seven aggregate score groups and species-level metrics. For each mode, none of the six original per-species CSVs was found at its expected output path: reference species metrics were recomputed from the checksum-verified original prediction rows and compared with the fresh six RF species CSVs. This distinction is preserved in the receipts.

For the species holdout, prediction comparisons outside tolerance total **137,322 model–row comparisons across six models**. This is not a count of unique test observations. M0–M5 counts are 28,575; 25,278; 31,544; 26,716; 24,147; and **1,062**, respectively. M5's maximum absolute difference is `0.005421749577093515`. All compared values are finite. Rounding both aggregate reductions to −0.209% does not establish numerical reproduction.

The replay environment was Python 3.12.16, pandas 2.2.3, NumPy 2.3.5, scikit-learn 1.7.2, PyArrow 23.0.1, joblib 1.5.3 and SciPy 1.18.1. Static historical scikit-learn 1.7.2 metadata comes from the archived **temporal M0 model only**; no pickle was executed, and this hint does not bind the historical spatial/species environments or supply a full runtime lock.

## Archived prediction and receipt diagnosis

[No-fit diagnostic run 44973229972](https://github.com/zuizui0223/chocho/actions/runs/44973229972), at `4e5f42e6e934bf73f074e00e8d46d751754b54f5`, verified the following independently of new model fitting:

- All **49/49 original aggregate score comparisons** passed when recalculated from the original species predictions; maximum absolute difference `6.938893903907228e-18`.
- All 208,482 reconstructed test keys matched the original predictions in exactly the same order. Targets and baselines matched exactly.
- The native original score JSON matched its manifest SHA-256. Normalizing its 82 CRLF line endings to LF reproduced the frozen public JSON SHA-256 exactly.
- Both the panel and species-prediction Parquet files record pandas **2.3.3** and Arrow writer **24.0.0**. These are concrete writer-version differences from replay; they do not prove the complete historical fitting environment or establish a cause of the numerical discrepancy.
- The canonical `scripts/temporal/fit_forecast.py` and allowlisted environment files were absent at the inspected archive paths. This bounded lookup does not exclude a differently named file elsewhere. Selected ZIP timestamps were all 1980-01-01 and cannot establish execution order.
- The native score receipt lacks `sensitivity`, `input_panel` and `rarefaction_n`, which the frozen final fit script writes unconditionally. This establishes a difference in receipt-generation history, not a demonstrated change in fitting logic. The older-schema temporal receipt nevertheless reproduced successfully.

This rules out a mismatch between original predictions and public scores, and a mismatch in the reconstructed test keys, targets or baselines. The historical fit-script revision and complete runtime remain unbound. The controlled M0 probe below tested the observed reader-version differences; no all-mode reproduction claim follows from the no-fit check.

## Controlled M0 reader/runtime probe

[Run 44974105799](https://github.com/zuizui0223/chocho/actions/runs/44974105799), at `e85e41624857d031caf32cb376e02839c4af874c`, completed the predeclared four-arm comparison on one Linux x86-64 runner. pandas 2.2.3 versus 2.3.3 and PyArrow 23.0.1 versus 24.0.0 were varied independently and together; Python 3.12.16, the other pinned packages, thread controls, source, panel and RF settings were held fixed.

All four arms produced **identical complete preflight fold records**: ordered training/test keys, float32 feature arrays including layout, targets, raw inverse-frequency weights, normalization scalars and normalized weights. The decoded, canonically keyed archived prediction oracle was also identical. Only the baseline five-fold M0 was actually fitted; the other three fits were skipped because their measured estimator inputs were identical. Skipped arms establish input equality, not additional fitted-model reproductions.

| M0 species fold | Test rows | Predictions outside tolerance | Maximum absolute difference from archive |
|---|---:|---:|---:|
| 0 | 50,404 | 0 | `4.440892098500626e-16` |
| 1 | 54,697 | 0 | `4.440892098500626e-16` |
| 2 | 41,190 | 0 | `4.440892098500626e-16` |
| 3 | 33,431 | **28,575** | **`0.007789458622820877`** |
| 4 | 28,760 | 0 | `3.3306690738754696e-16` |

The overall count and maximum error exactly repeat the earlier M0 comparison summary. The prior fresh prediction vector was unavailable for a direct row-by-row comparison between replays. The result localizes M0's discrepancy to fold 3 and finds no change in its measured inputs across these reader versions. It does not establish the historical fitting cause, the historical code revision or complete runtime, or validate the other five RF models. The original temporal/spatial PASS and species FAIL remain the governing outcomes.

`SPECIES_M0_CONTROLLED_RUNTIME_20261122.json` has SHA-256 `daa90f340ac98debf45ccfe228ab33b9379f7632e933741f7e6930a0e92d7fb2`. Its aggregate-only artifact `14734179419` is 4,851 bytes, ZIP SHA-256 `d4717919819b660a7d8bc9df83061d4a51e05433425b075deb13d98e2f7154f6`, expiring 20 February 2027. Private inputs, temporary predictions and source logs were removed by the workflow cleanup.

## Durable aggregate evidence

The nine original restoration/execution/comparison JSONs are preserved byte-for-byte under `results/geb_v05_private_source_recovery_20261120/`. `ORIGINAL155_PRIMARY_REFIT_STATUS_20261122.json` records their SHA-256 values, each mode's actual outcome and the separately verified no-fit diagnostic. `SPECIES_ARCHIVED_LINEAGE_20261122.json` has SHA-256 `0c31c1a7ea00460fc33b0278391b5cd8c66d8d1c3b6f40ca8935f0200bd3703c`.

| Temporary Actions artifact | Bytes | ZIP SHA-256 |
|---|---:|---|
| Temporal `14733188847` | 5,381 | `dcfa9860bd5a2a6b549695dc78e1599b070954ceb6c6c637d6e61beef7ba5022` |
| Spatial `14733942476` | 5,872 | `b26ddc17f4fb7e27e1a40a629a1b0e4a0e07508ace093984e6179dc6f657b492` |
| Species `14734715640` | 6,062 | `a78abfd7d20fe394e9e8d1ecff17668b17361d6df8650dd2496f63ce6295fbc4` |
| No-fit lineage `14734109643` | 2,414 | `9f6da386cf6e2735c177ef9c464b35afb11e0dcbef8edc0319056dec4dccff7b` |

Artifacts expire on 20 February 2027; they are not permanent reviewer deposits. Private observations, panels, predictions and serialized models are absent from these public aggregate receipts. The source release remains Draft.

## Remaining gates and manuscript decision

Recovery found 12 of 13 previously sought paths. `work/remote_sensing/host_region_counts.csv` was absent at its recorded path. Recovering its five columns from the final panel would be circular evidence for raw-to-panel reconstruction; no such reconstruction is claimed. The four additional original sensitivity fits (`effort20`, `radius250`, `exclude_top_provider`, `precision100`) are outside this run.

The [journal export at paper head `03e9450`](https://github.com/zuizui0223/chocho/blob/03e9450c832bb769c42ca983bffdc327ecc8ccea/manuscript/integrated_butterfly_resource_access_20261109/GEB_REPRODUCIBILITY_UPDATE_QA_20261122.md) now identifies the two verified primary refits and explicitly records the adverse, numerically unmatched species result. Scientific results and figures remain unchanged. The 39-file internal Fig1/2/3/4 package is a separate limited replay: its Fig4 component replays saved scores and does not execute the private 155/224 models. Full-SI replay, source rights, indirect anonymity, stable permitted reviewer access and author approvals remain separate open gates. No coauthor message, review request, Ready state, merge, release publication or submission is authorized by these checks.
