# Surviving Windows runtime reproduces original distances and temporal refit

The original local execution log records the panel builder running on 29 October 2026 at 08:01:39 UTC with C:/Program Files/Python310/python.exe. Its output reports 280,988 species-cell pairs and 842,964 period rows. The attached excerpt is restricted to three forecast commands and their outputs; the full chat history is excluded.

A single distance-only replay on the surviving local Windows runtime compared all 280,988 keys. There were zero differing float64 distances, float32 distances, or float32 log1p(distance) values. Maximum absolute distance difference was zero.

The fixed original panel was copied to an isolated folder, its distance column replaced with those freshly computed values, and all six primary temporal RF models were fitted once using the unchanged original fitting script and lock. The existing audit comparator passed all 175 groups across 208,482 test rows and 155 species. Maximum prediction difference was 5.551115123125783e-16, within the pre-existing rtol=1e-6 and atol=1e-9. Host improvement was 4.120175607931823%.

The surviving runtime is Windows / Python 3.10.11 / NumPy 2.2.6 / SciPy 1.15.3 / scikit-learn 1.7.2 / pandas 2.3.3 / pyproj 3.7.1 / PROJ 9.5.1 / pyarrow 24.0.0. Current package directory dates predate the original execution. Dates establish context, not proof that all installed bytes remained unchanged. The build metadata was collected in a separate inspection process; its reported thread counts are not a training-process trace.

This resolves the local distance-reproduction question. It does not isolate which Windows/Linux, NumPy, library-build or other numerical difference caused the previous mismatch. The old log excerpts include preliminary fits and do not certify the revision of every historical run. This check retains 95 original panel columns and reconstructs the distance only; it is not a fresh raw-to-all-96-column replay. Unseen-species and sensitivity refits were not rerun in this check.

Source hashes:
- original panel: 5ae5b48b752fad05067537de5fb95d4ee74f17f3ba1114443d337bd7fd3695a1
- fitting script: a34ad4cd1e4f40460c4f6bf6ed0b22ba3885b17d20ac17f74bd6564d0618ce6d
- distance/panel-builder script: 99af0755252f2bb81c53f1b07a504e34d15092ccbd8a10a0bc16dde0e0b57ec6

Run replay_distance.py from the recovered analysis root. compare_local_refit.py expects the same isolated refit layout documented by its source. The audit comparator was retrieved from audit PR #189 at 8f56dc9feff6dd82f9237c673243b61403041229.
