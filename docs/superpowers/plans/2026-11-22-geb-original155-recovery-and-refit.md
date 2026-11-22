# GEB original 155-species recovery and refit implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Recover the existing original 155-species forecast inputs, verify their source identity, and independently rerun the original fixed model before changing the manuscript's reproducibility claim.

**Architecture:** Read the existing private draft release through an authenticated, read-only Actions job. Verify immutable asset and source hashes before extracting the required inputs. Run the pinned original model in a separate directory and publish only aggregate verification receipts.

**Tech Stack:** GitHub Actions, Python, Parquet, original scikit-learn random forests.

**Spec:** `manuscript/integrated_butterfly_resource_access_20261109/GEB_V05_PRECOAUTHOR_PREPARATION_GATE_20261121.md`, especially A4, and `results/geb_v05_private_source_recovery_20261120/ORIGINAL155_FORECAST_SOURCE_AND_PANEL_GATE_20261121.json`.

## Current state — 22 November 2026

1. **Completed:** original archive/panel recovered and hash-verified; temporal/spatial fits and the spatial ≥100-km evaluation reproduced. The complete species fit ran and its numerical comparison failed. All outcomes are retained.
2. **Completed:** original species predictions reproduced all 49 frozen score fields; ordered test keys, targets and baselines are exact. The 39-file four-component internal replay passed. Corrected manuscript/SI exports passed 32-page visual/source QA and are saved.
3. **Completed:** controlled M0 runtime diagnosis, run `44974105799`. All four arms had identical measured estimator inputs; only baseline M0 was fitted. Its 28,575 discrepancies were confined to species fold 3. The historical fitting cause remains unresolved.
4. **Recorded outcome:** completed diagnostics, source receipts and current gates are tracked on Draft PR189 with the bounded reproduction status. Paper head `03e9450` passed all 12 workflow runs / 13 check jobs. The next source-reproduction blocker is the unbound historical species fit revision/runtime; broader Stage A rights, raw-to-panel, full-SI, stable reviewer-access and author-approval gates stay separate.

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

- [x] Verify the receipt asset (775476029) against SHA-256 `b8eff912a43373424e5f31ee8714c89d672f6d44c4011bbd81514e8a9afc3aa2` and SHA256SUMS asset (775476002) against `3cf15b18d68f27425e3b78630769b13e8a4e0c7232c6dd7016f7deaf857bd140`.
- [x] Read the receipt structure without exposing source-machine metadata. Record only safe relative paths, file hashes, counts and recognized scientific input summaries.
- [x] Download and verify the four original split assets, reconstruct the ZIP and extract required original inputs into ephemeral runner storage.
- [x] Verify panel/lock identity and exact cohort/transition sizes; capture source and current environment evidence.

## Task 2: Independently refit and compare

**Files:** `scripts/verify_GEB_v05_original155_refit_20261122.py`, focused comparison/recovery tests, and the Actions workflow.

**Interfaces:** Consume the verified panel, model lock and original `fit_forecast.py`; produce a sanitized source-versus-refit numeric comparison.

- [x] Add meaningful failing checks for a changed cohort, changed settings, corrupted input and altered aggregate/per-species predictions before implementation.
- [x] Install the byte-verified fit script at `<run>/scripts/temporal/fit_forecast.py`; keep the source archive outside that directory.
- [x] Execute original temporal and spatial fits; derive the 100-km evaluation from the new spatial predictions, not a new buffer-fit model.
- [x] Compare all seven model score groups, original per-species metrics and prediction identities where available, with predeclared tolerance; record actual library versions.
- [x] Run the independent reviewer over the retrieval and comparison code and inspect fresh CI evidence.

## Task 3: Record the verified current state

**Files:** Source-audit receipt, pre-coauthor gate and current manuscript export source-status text where warranted.

- [x] Update the missing-panel gate only after real source input checks; update refit claims only after actual successful model execution and comparison.
- [x] Preserve separate reviewer access, rights and anonymity gates.
- [x] Keep a concise handoff of verified work and the exact remaining blocker.

## Execution record

- 2026-11-22: Existing isolated source-audit checkout restored at original PR189 head `be3bf40a59023b2336da3876d50fd2c9e8ca76ae`.
- Ruling: Proceed without a design-approval pause under the user's explicit action-first instruction; all changes stay in the existing Draft preparation branch.
- Ruling: Use one small authenticated inventory run before implementing extraction because the recovery receipt schema and original runtime versions are not yet readable through the connector. This prevents guessing source identities.
- Original paper PR179 head `9a718997bce21f5b99329a91844bd31c375c6d16` has 12/12 successful workflow runs, including the previously pending 104-island audit.
- Task 1 integration check: run `44969368087` stopped at reading draft release metadata with `contents: read`; no private asset was exported. GitHub documents that only credentials with push access can see draft releases (`https://docs.github.com/en/rest/releases/releases`). The next run changes only the required repository-content permission to `write`, keeps every API operation GET-only, and adds a sanitized HTTP status for diagnosis. This does not publish or modify the draft release.
- Task 1 inventory verified by successful CI `44969486664`, artifact `14733300476`: 3,032 source records, 59 occurrence-part paths, original panel SHA256 `5ae5b48b752fad05067537de5fb95d4ee74f17f3ba1114443d337bd7fd3695a1`, original model-lock hash matching the frozen public source. Full archive digest `9e39b0068c3617741f9b275a1edfddd574c93d620853e956bd706ad067a53f9f` is now independently read from the byte-verified checksum file. Actual panel restoration and fitting remain pending this step.
- Whole recovery/refit code received an independent fresh-context review. Important logging finding fixed with a dedicated approved-diagnostic exception and a failing-then-passing private-marker regression. The model-version step title was made conditional-neutral because a source metadata hint is not the complete historical runtime lock. The reviewer found no critical archive or numerical-integrity blocker; actual archive contents, runtime headroom and fit agreement are integration gates.

- Task 1 full restoration passed in all three original155 jobs of run `44970180129`: all four part hashes and full ZIP hash match; the 842,964-row/155-species panel and the model lock match the archived source; recording-day datasetKey is retained. Twelve of 13 old sought paths are present; `work/remote_sensing/host_region_counts.csv` is absent at that path. This is not raw-to-panel reconstruction.
- Temporal fresh execution and comparison completed successfully: 193,703 training / 208,482 test rows, all seven score groups and six RF species summaries; 144 comparison groups, zero outside fixed tolerance, maximum prediction absolute difference 5.551115123125783e-16. Source static sklearn metadata and replay version are both 1.7.2. Spatial/species fits remain active until their own completion receipts pass.
- In parallel, the four-component, 39-file internal Fig1/2/3/4 limited replay passed local independent extraction and CI `44970859551` at source head `81b869e37a9ca2ec711c99d775c6d1d34d9caf5e`; inner ZIP SHA256 `40369ba5820f0e40ed56758f64165f951c82c6d3d64e5613175533d3c2055c95`. It remains separate from private RF refits and whole-SI/reviewer access.
- Bounded upstream follow-up: the original consumer reads only `species`, `region`, `native_host_count`, `introduced_host_count`, `introduced_share` from `work/remote_sensing/host_region_counts.csv`. Original CSV bytes/writer/hash were not located in the checked Git snapshots, small October download-receipt artifacts or the nine-member barrier/port recovery ZIP. Adjacent frozen source code locates `source_augmented_model` and `work/remote_sensing/host_source/extracted`, with original accepted HOSTS edges (SHA256 `0a084fb5273e4780b03e015205c87e4545c69f0f8e517feb1a2a3879889508e9`), native distribution (SHA256 `c731906315f7452f83ce302c1d36c75ad94afc24ce7e24d461360312244ac558`), contemporary distribution, and BCE additional edge/native/contemporary files. Exact supplement hashes and serialization remain unverified. Recovering these five columns from the model panel would be circular evidence for a raw-to-panel claim; no such reconstruction is claimed. A future bounded alias check can inspect the already verified 569,496-byte recovery receipt asset for `host_region_counts`/`host_source`, without downloading the large archive again.

- Actual run `44970180129` completed with temporal PASS (144/144), spatial PASS (174/174 including ≥100-km) and species numerical comparison FAIL (68/144 passing) after successful fresh execution. Species model-row comparisons exceeded unchanged tolerance in 137,322 cases across six models; M5 affected 1,062 rows. The adverse aggregate direction remained −0.208520153993% original vs −0.208933309709% refit. The failure is retained, not converted to a successful refit verification.
- A new no-fit diagnostic at `4e5f42e6e934bf73f074e00e8d46d751754b54f5` checks byte-bound original species predictions against frozen scores, Parquet writer/version metadata and selected code/receipt hashes. Twenty focused tests and independent source/privacy review passed. No original source, tolerance, model setting or three-fit workflow was modified. Source receipts reveal a different historical receipt writer: the native CRLF species/temporal receipts lack three fields always written by the archived final script. This proves receipt-generation history differs, not that fitting logic changed.

- No-fit lineage run `44973229972` completed successfully: all 49 original-score fields match archived predictions, all 208,482 test keys match in order and targets/baselines are exact. Both archived Parquet writers report pandas 2.3.3 / Arrow 24.0.0, compared with replay 2.2.3 / 23.0.1. Writer metadata does not prove the historical fitting runtime. Canonical script/environment paths were absent and selected ZIP timestamps were fixed at 1980-01-01. No causal conclusion was made.
- The controlled M0 probe passed six new fixtures and independent review, including synthetic execution of the five actual source boundaries with estimator fitting prohibited. Code `e85e41624857d031caf32cb376e02839c4af874c` (reviewed tree `4c09e6629a6193d5cc5ab3831e91631ce9f534a1`) is executing in run `44974105799`: same runner, four pandas/Arrow combinations, other versions/settings held fixed, original prediction oracle checked across readers, actual X/y/weight/order hashes, baseline plus only distinct-input M0 fits. Skipped arms cannot be reported as fitted.
- Paper source-status update pushed to existing Draft PR179 as `03e9450c832bb769c42ca983bffdc327ecc8ccea` after actual 32-page local visual validation. All scientific numerical tokens, source references and six figure images remain unchanged. Final manuscript wording refers to the archived script, because the precise historical executable is not yet bound. The four document deliverables and expanded internal 39-file replay ZIP are saved.

- Controlled M0 run `44974105799` completed: all four full preflight records and decoded archived prediction oracles are identical; only baseline five-fold M0 was fitted. Fold 3 contains all 28,575 discrepancies; the other four folds agree within tolerance. Prior outside-count/max-error summaries repeat exactly, but the prior fresh prediction vector was unavailable for direct vector equality. Historical fit-code/runtime causation remains unresolved. Source-safe receipt SHA256 `daa90f340ac98debf45ccfe228ab33b9379f7632e933741f7e6930a0e92d7fb2`.
- Final paper head `03e9450c832bb769c42ca983bffdc327ecc8ccea` passed 12/12 workflow runs and 13/13 check jobs. Final local PDFs are 21 main + 11 SI pages and remain distinct from the CI render. Source-current handoff is `results/geb_v05_private_source_recovery_20261120/GEB_V05_STAGE_A_VERIFIED_HANDOFF_20261122.json`.
