# Joint host and marine-circuit comparison Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Execute the fixed conditional host/circuit prediction comparison for the original51 islands, preserving prior evidence and exposing all predictions.

**Architecture:** A separate runner validates a source/protocol lock, reuses the source loader, circuit transformation and OLS prediction function, and fits the five declared feature sets. Prior G/C predictions are a reproduction gate. Existing source/circuit code and results are not modified.

**Tech Stack:** Existing Python3.12+, NumPy, pandas, SciPy and optional matplotlib analysis environment; no new external datasets or dependencies.

**Spec:** [PROTOCOL.md](PROTOCOL.md).

## Global constraints

- Exactly51 eligible islands,1,275 unordered pairs,51 folds,1,225 training and50 test pairs per fold.
- Five models G,H,C,HC,I; two butterfly targets only; five comparisons per target.
- Native AIS matrix and original uniform/wind matrices; original194-port host inventory and port-count controls; original182 circuit terminals.
- All source hashes and protocol must pass before fitting; prior G/C fold values must match within1e-12 absolute and relative tolerance.
- No model/transform selection by new outcomes; no independent-pair or independent-fold inferential claims.

## Review focus

1. Matrix labels can be permuted while dimensions remain correct: reject ID-order mismatch rather than silently pair wrong islands.
2. Host similarity can be zero, while an earlier response must not predict itself: preserve host zeros and restrict response-specific earlier control.
3. Held-out observations can accidentally enter training: changing a withheld island's responses must leave its own predictions unchanged.
4. A changed frozen source can produce apparently reasonable numbers: fail the input gate before output fitting.
5. Summary percentages can hide a denominator or duplicated-prediction error: reconstruct all summaries from saved prediction rows and reproduce the prior G/C folds.

## Task1: Fixed joint runner and scientific regression tests

**Files:** create `scripts/run_island51_host_circuit_joint.py`; create `tests/test_run_island51_host_circuit_joint.py`.

**Interfaces:** consume `source()` and `train_predict()` from the existing geometric-holdout script, `make_base_features()` from the observed runner and `resistance_feature()` from the circuit module. Produce fixed-model feature dictionaries, model/fold and comparison summaries, individual predictions and a source/validation receipt. CLI takes `--out-dir` and optional `--input-lock` with repository-relative defaults; output is separate from source and prior result directories.

- [ ] Write failing tests for required model membership, secondary interaction construction, zero-host preservation and the earlier-control rule.
- [ ] Write failing tests for reordered matrix IDs and tampered input/protocol bytes.
- [ ] Write a failing leakage test using complete synthetic51-island pairs: alter only the chosen island's held-out responses and require its predictions to remain equal.
- [ ] Implement the minimal source gate, feature sets, shared train-only prediction loop and summary aggregation; reuse existing mathematics.
- [ ] Run the targeted tests, then the complete repository suite.
- [ ] Run the real fixed comparison only after `input_lock.json` is frozen; require all prior G/C fold values and summaries to reproduce.

## Task2: Independent result reconstruction and review

**Files:** add independent review/verification receipts under this results directory. Inputs are the runner's complete prediction and fold outputs, protocol and source lock.

- [ ] Independently group25,500 saved prediction rows into510 model-fold losses and10 comparison summaries.
- [ ] Verify test/training support,51-island coverage, both representations of each unique pair, and no omitted adverse folds.
- [ ] Review the runner/spec for material methodological defects; resolve actual defects without outcome-driven model changes.

## Task3: Report, figure, CI and persistence

**Files:** add `README.md`, an exact-data summary figure and reproducible renderer in the new result scope; add a scoped GitHub workflow for this continuation.

- [ ] Explain both conditional comparisons, source scales and retrospective limits; include both favorable and adverse results and secondary temporal/interaction comparisons.
- [ ] Render and inspect a scientific figure showing model errors and the two conditional increments, with no independent-fold uncertainty bars.
- [ ] Run the documented command from frozen inputs and the full repository tests; record actual environment and results.
- [ ] Preserve prior source/output hashes, save the continuation commit and create a draft PR targeting the existing circuit branch; verify remote CI.

The user's action-first continuation instructions authorize implementation after recording this plan. Routine analysis and implementation choices are made by the assistant; no additional approval pause is introduced.
