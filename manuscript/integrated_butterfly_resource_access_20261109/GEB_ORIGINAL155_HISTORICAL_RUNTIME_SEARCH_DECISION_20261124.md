# Original155 — historical runtime retrieval decision (24 November 2026)

**Decision: search the original October 2026 execution workspace or its backup, not another substitute ecological dataset.** The canonical original 842,964-row × 96-column panel, original frozen builder/model scripts, 155-species cohort, original 208,482-row temporal predictions, and original RF model settings have already been recovered and SHA-verified. Re-fitting the original canonical panel successfully reproduced the published temporal/scientific scores. Rebuilding all panel columns from upstream source files gives complete numeric agreement within predeclared tolerance, but exactly one geographic predictor (`distance_km`) differs at floating-point bits and drives an RF refit discrepancy; replacing that sole column in a **diagnostic counterfactual** restores all 175 temporal checks. Thus the unresolved question is original geographic numeric runtime/precision, not lack of original model output.

## Archive-only search now completed

The 3,032-entry, SHA-256-verified 21-Nov *original* 5.55GB Draft-release source manifest was scanned without downloading raw source files or exporting private paths. [CI #45325608482](https://github.com/zuizui0223/chocho/actions/runs/45325608482) **PASS** found **zero candidates** in each of the following conservative filename classes:

- Environment files (`requirements*.txt`, `environment.y*ml`, `conda-lock`, `uv.lock`, `Pipfile.lock`, `poetry.lock`, etc.)
- pip/conda package snapshots
- Original forecast/panel/distance `*.log` execution logs
- Separate `*distance*`/`*nearest*` CSV/Parquet/Numpy cache files
- Explicit pyproj/PROJ/sklearn environment-version metadata files

This **bounded filename inventory** does not prove a log cannot exist under another name or on an earlier local machine. The original October 29 PR #101 stores inspection scripts and source-model summaries, not an independently verified full historical geographic runtime. Later source manifests should not be misdescribed as 2026-10-29 runtime evidence.

## The one source worth recovering now

An original **2026-10-29 period** workspace backup or a copy of its Python virtual environment/conda metadata, ideally including:

1. Exact Python/environment receipts: `pip freeze`, `conda list --explicit`, `environment.yml`, `requirements*.txt`, `pyproject.toml`, `conda-meta/`, `.venv/pyvenv.cfg`, shell/notebook logs, or session metadata with actual `pyproj`, bundled **PROJ**, `numpy`, `scikit-learn` and `pandas` versions.
2. The original geographic distance construction or saved intermediate `work/temporal/forecast/species_cell_period_panel.parquet` (already recovered) and any historical pre-panel nearest-record/centre-coordinate tables and generation logs.
3. Exact build-script source or Git commit from the time of execution (published script was recovered, but historical provenance of the execution environment is missing).

**Privacy:** retain original inputs/venv/logs privately; do not upload source occurrence files, provider identifiers, local usernames, credentials, or raw archives to public PR #189. Source logs can be locally inspected and only release a sanitized package/version receipt.

## Bounded validation and stop condition

If historical environment evidence is recovered, run **one** source-locked original geographic reconstruction on the same known 280,988 species×cells, compare exact float64 and actual model `log1p(distance_km)` float32 values against the archived panel, and then refit the original 155-species temporal model (208,482 heldout test rows) without any predictor substitutions.

**If the environment/backups do not exist**, do not block the ecological manuscript indefinitely on byte-level distance reproduction. Use the SHA-verified canonical original panel as the explicitly named analysis input, report the successful archived-panel model refit, maintain separate numerical provenance for the 13-input upstream reconstruction, disclose the model's small numeric sensitivity (M5-vs-M3 +4.120% archived vs +3.936% from independent derived-panel rebuild), and preserve the positive/adverse sensitivity and unseen-species results. This does **not** clear independent GBIF-source/raw-to-panel rebuild, license, reviewer data/code anonymity, or author permissions; those remain Stage A gates.

### Scope cautions

- [PR189 distance intervention #45322533520](https://github.com/zuizui0223/chocho/actions/runs/45322533520) is an explanation of a model-input numeric difference, **not a biological mechanism** or independent distance reconstruction.
- [Controlled pyproj 3.7.1 vs 3.7.2 check #45323737296](https://github.com/zuizui0223/chocho/actions/runs/45323737296) had identical PROJ 9.5.1 and identical remaining exact/numeric distance differences, so repeatedly varying these two identical-PROJ builds is not a useful next step.
- Previously archived RF models expose scikit-learn 1.7.2; this alone does not identify the historic Python/numpy/PROJ runtime.
- Original scientific paper PR #179 unchanged; source-audit PR #189 stays Draft. No coauthor contact, public original data release, or journal submission authorized.
