# GEB original 155-species recovery and refit implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Recover the existing original 155-species forecast inputs, verify their source identity, and independently rerun the original fixed model before changing the manuscript's reproducibility claim.

**Architecture:** Read the existing private draft release through an authenticated, read-only Actions job. Verify immutable asset and source hashes before extracting the required inputs. Run the pinned original model in a separate directory and publish only aggregate verification receipts.

**Tech Stack:** GitHub Actions, Python, Parquet, original scikit-learn random forests.

**Spec:** `manuscript/integrated_butterfly_resource_access_20261109/GEB_V05_PRECOAUTHOR_PREPARATION_GATE_20261121.md`, especially A4, and `results/geb_v05_private_source_recovery_20261120/ORIGINAL155_FORECAST_SOURCE_AND_PANEL_GATE_20261121.json`.

## Global constraints

- Stage A preparation only: keep PRs Draft; no merge, submission, coauthor message or review request.
- Original source commit: `42f8fcc3a388e96456c76180d7e56d75dd2ed48c`.
- Original model panel: 842,964 rows, 155 species; 193,703 training and 208,482 temporal test rows.
- Original random forest: 160 trees, depth 12, minimum leaf 40, feature fraction 0.7, seed 20261029, 4 jobs.
- Original 155 and 224 analyses share GBIF DOI `10.15468/dl.umqmsf` but have distinct cohorts and estimands.
- Draft release 453667492 is a source recovery location, not an anonymous reviewer deposit or permission to redistribute raw inputs.
- Existing original results and input archive remain separate from newly fitted outputs.
- Do not execute recovery-provided scripts or load pickled models to inspect their metadata.

## Review focus

- A download returning JSON/HTML instead of asset bytes must fail its pinned SHA-256 check.
- Split-file order, archive digest, ZIP traversal and individual source hashes must be verified before use.
- Published receipt fields must exclude machine paths, original observations, provider identifiers and credentials.
- Matching stored scores alone must never set `model_refit_completed` to true.
- Runtime version differences or actual numerical discrepancies must be reported rather than hidden by changed tolerances.

## Task 1: Recover authenticated source inventory

**Files:** `.github/workflows/audit-GEB-v05-original155-recovered-archive-20261122.yml`; later a focused recovery utility under `scripts/`.

**Interfaces:** Read draft release 453667492 assets using the job token; emit only a sanitized inventory to an Actions artifact. A small first run establishes the actual receipt schema before archive extraction is implemented.

- [ ] Verify the receipt asset (775476029) against SHA-256 `b8eff912a43373424e5f31ee8714c89d672f6d44c4011bbd81514e8a9afc3aa2` and SHA256SUMS asset (775476002) against `3cf15b18d68f27425e3b78630769b13e8a4e0c7232c6dd7016f7deaf857bd140`.
- [ ] Read the receipt structure without exposing source-machine metadata. Record only safe relative paths, file hashes, counts and recognized scientific input summaries.
- [ ] Download and verify the four original split assets, reconstruct the ZIP and extract required original inputs into ephemeral runner storage.
- [ ] Verify panel/lock identity and exact cohort/transition sizes; capture source and current environment evidence.

## Task 2: Independently refit and compare

**Files:** `scripts/verify_GEB_v05_original155_refit_20261122.py`, focused comparison/recovery tests, and the Actions workflow.

**Interfaces:** Consume the verified panel, model lock and original `fit_forecast.py`; produce a sanitized source-versus-refit numeric comparison.

- [ ] Add meaningful failing checks for a changed cohort, changed settings, corrupted input and altered aggregate/per-species predictions before implementation.
- [ ] Install the byte-verified fit script at `<run>/scripts/temporal/fit_forecast.py`; keep the source archive outside that directory.
- [ ] Execute original temporal and spatial fits; derive the 100-km evaluation from the new spatial predictions, not a new buffer-fit model.
- [ ] Compare all seven model score groups, original per-species metrics and prediction identities where available, with predeclared tolerance; record actual library versions.
- [ ] Run the independent reviewer over the retrieval and comparison code and inspect fresh CI evidence.

## Task 3: Record the verified current state

**Files:** Source-audit receipt, pre-coauthor gate and current manuscript export source-status text where warranted.

- [ ] Update the missing-panel gate only after real source input checks; update refit claims only after actual successful model execution and comparison.
- [ ] Preserve separate reviewer access, rights and anonymity gates.
- [ ] Keep a concise handoff of verified work and the exact remaining blocker.

## Execution record

- 2026-11-22: Existing isolated source-audit checkout restored at original PR189 head `be3bf40a59023b2336da3876d50fd2c9e8ca76ae`.
- Ruling: Proceed without a design-approval pause under the user's explicit action-first instruction; all changes stay in the existing Draft preparation branch.
- Ruling: Use one small authenticated inventory run before implementing extraction because the recovery receipt schema and original runtime versions are not yet readable through the connector. This prevents guessing source identities.
- Original paper PR179 head `9a718997bce21f5b99329a91844bd31c375c6d16` has 12/12 successful workflow runs, including the previously pending 104-island audit.
- Task 1 integration check: run `44969368087` stopped at reading draft release metadata with `contents: read`; no private asset was exported. GitHub documents that only credentials with push access can see draft releases (`https://docs.github.com/en/rest/releases/releases`). The next run changes only the required repository-content permission to `write`, keeps every API operation GET-only, and adds a sanitized HTTP status for diagnosis. This does not publish or modify the draft release.
- Task 1 inventory verified by successful CI `44969486664`, artifact `14733300476`: 3,032 source records, 59 occurrence-part paths, original panel SHA256 `5ae5b48b752fad05067537de5fb95d4ee74f17f3ba1114443d337bd7fd3695a1`, original model-lock hash matching the frozen public source. Full archive digest `9e39b0068c3617741f9b275a1edfddd574c93d620853e956bd706ad067a53f9f` is now independently read from the byte-verified checksum file. Actual panel restoration and fitting remain pending this step.
- Whole recovery/refit code received an independent fresh-context review. Important logging finding fixed with a dedicated approved-diagnostic exception and a failing-then-passing private-marker regression. The model-version step title was made conditional-neutral because a source metadata hint is not the complete historical runtime lock. The reviewer found no critical archive or numerical-integrity blocker; actual archive contents, runtime headroom and fit agreement are integration gates.
