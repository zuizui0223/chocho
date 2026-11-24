# GEB original155 source model serialization audit — 24 November 2026

**PASS within audited scope.** A recovery worker read the existing original private draft-release archive and SHA-verified each of 30 saved `joblib` model-object files against the registered 3,032-entry source manifest. The tool extracted only bounded **static `_sklearn_version` string metadata** and did **not** call `pickle.load`, `joblib.load`, run any saved model, export model bytes, export occurrence/panel rows, or publish the original release.

| Original saved model group | Hash-verified original model objects | Static scikit-learn version string |
|---|---:|---|
| Primary temporal | 6 | **1.7.2** for all six |
| Effort threshold 20 provider days | 6 | **1.7.2** for all six |
| Radius ≤250 km | 6 | **1.7.2** for all six |
| Exclude largest provider | 6 | **1.7.2** for all six |
| Coordinate precision ≤100 m | 6 | **1.7.2** for all six |
| **Total** | **30/30** | **all readable; all 1.7.2** |

**Interpretation and negative diagnostic:** Recent 155-species source-panel refits also ran under scikit-learn 1.7.2, including the two numerically discordant temporal sensitivity results. A mismatch of the **recorded scikit-learn version label** therefore does not explain those sensitivity prediction differences. This does **not** prove identical fit code revision, BLAS/joblib/numpy runtime, training sequence or linkage between the serialized models and particular archived prediction rows. The unseen-species folds have **no corresponding models among these 30 temporal saved objects**; species fold 3 remains unresolved and this static result does not bind its historical library version.

The first workflow attempt failed at a split archive part size/SHA check, **before reading any model**. The single retry passed all original archive and individual-model SHA checks. This is a source-transfer failure rather than a biological or random-forest predictive result. Do not ignore or relabel the first failure as an explained network outage.

- Reproducible scanner: `scripts/audit_GEB_v05_original155_saved_model_versions_20261124.py` (synthetic static-parser tests in `tests/test_GEB_v05_original155_saved_model_versions_20261124.py`).
- [Successful attempt 2 CI #45311652838](https://github.com/zuizui0223/chocho/actions/runs/45311652838).
- Original private archive SHA-256: `9e39b0068c3617741f9b275a1edfddd574c93d620853e956bd706ad067a53f9f`.
- Aggregate-only Actions artifact `14883324615` has ZIP SHA-256 `73c9962f7d1ebdc9b0391c453c8bc02881a15381d09c1f28c2894399a4e038b0`, and contains all 30 source SHA/size and version entries for independent inspection. It is **not** a long-term reviewer data deposit.

Scientific original results remain untouched. Remaining model reproducibility issues are source-to-training-runtime lineage for the unseen-species fold 3 and two temporal sensitivity discrepancies, raw-to-panel reconstruction, complete SI coverage and permitted reviewer access. This report does not authorize coauthor review, an archive DOI, branch merge or journal submission.
