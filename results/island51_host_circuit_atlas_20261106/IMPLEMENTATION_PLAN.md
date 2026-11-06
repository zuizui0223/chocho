# Island host/circuit gain atlas Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to execute the scoped tasks. Frozen checklist bytes are retained; completion is recorded separately.

**Goal:** Turn the fixed 51-island joint predictions into a complete geographic account of conditional gains, adverse contributions and their concentration.

**Architecture:** One standalone summary script lazily invokes the existing independent parent auditor, validates a new source/specification lock, joins original geographic metadata, and calculates descriptive identities from saved losses. Separate rendering and independent review use only the resulting tables. No fitting code is invoked.

**Tech Stack:** Existing Python environment; standard-library CSV/JSON/math for calculation and synthetic tests; NumPy/matplotlib for figures. No new datasets or dependencies.

**Spec:** [PROTOCOL.md](PROTOCOL.md).

## Global Constraints

- Parent commit17a4768cbb6cddf249449bd67e3ea1124e0db2f5; unchanged51 islands,1275 pairs,194 biological ports and182 marine terminals.
- Retain all510 comparison-fold records and both responses; primary geographic categories use only C→HC and H→HC.
- No refitting, model selection, new covariate regressions, political/island-name inference or independent-fold inference.
- Global contributions use the common reference-error sum; gross positive and adverse masses remain separate. Concentration k values are1,3,5.
- Preserve all parent files. Imports and unit tests must work without results/ artifacts; production auditing loads them lazily.

## Review Focus

1. Saved delta has the opposite sign from improvement: independently verify reference minus augmented and global additive identities.
2. Zero/tiny reference errors and cancellation can create misleading percentages: preserve undefined values and use gross mass for concentration.
3. Ports/names/countries can be joined to the wrong island: strict OBJECTID mapping and canonical support are required.
4. Port-group points can cross the dateline or omit marine-excluded ports: reproduce spherical means from all original ports.
5. A geographic map can imply local ecology or independent evidence: report the50-pair estimand, shared HC prediction and existing temporal/effort limitations.

## Task1: Summary runner and focused regression tests

**Create:** `scripts/summarize_island51_host_circuit_atlas.py`, `tests/test_summarize_island51_host_circuit_atlas.py`.

**CLI:** required `--out-dir`; optional `--input-lock` defaults to this directory's lock. Production source paths remain fixed by that lock. All output goes to the new directory; reject writes to parent result/source locations.

**Outputs:** the five CSV tables named in the protocol, `parent_prediction_audit.json`, `atlas_receipt.json`. All tables join using canonical string `island_id`; port provenance retains original `landmass_id`/`port_id` too. Decomposition rows retain parent keys/models/role/train/test counts and source loss values. Provide explicit signed MSE/MAE gains, island-relative percentages, global mean-loss contributions, global percentage-point contributions and strict signs. Geography carries original/marine counts, representative coordinates, all source port IDs/names/country codes and a smallest-ID example port. Conditional-status rows carry both contrasts and the five-class MSE category.

- [ ] Observe failing synthetic tests for sign reversal and reference-weighted percentage: references[0.01,0.09], gains[0.005,0.009] give global14%, not30%.
- [ ] Observe failing tests for zero denominators, all-zero/all-adverse masses and cancellation[+0.6,+0.4,−0.9,0]: gross positive1, adverse0.9, net0.1 and top-positive share60%.
- [ ] Test rank ties by numeric island ID, all sign combinations and MSE/MAE disagreement, missing/duplicate/misaligned joins, and dateline-safe spherical means including excluded marine ports.
- [ ] Implement the minimal lazy parent-audit gate, fixed joins, accounting, outputs and checksums. Synthetic tests use no actual results/ fixtures and remain active in the anonymous snapshot.
- [ ] Run focused tests and actual atlas only after the lock is frozen; preserve the parent hashes.

## Task2: Independent numerical and source review

**Create:** compact independent audit receipt under this new result directory. Do not import the new summary script to calculate expected identities.

- [ ] Independently reconstruct all510 gains and10 concentration summaries from the parent fold CSV, and verify all102 conditional-status rows and category counts.
- [ ] Verify geography coverage, exact original194-port membership,182 marine status, truthful labels and map coordinates.
- [ ] Review implementation/spec compliance; resolve material defects without changing the outcome-derived model set or concentration rules.

## Task3: Scientific report, figures, CI and persistence

**Create:** `README.md`, `render_atlas.py`, source-bound PNG/SVG figures and receipts under this result directory; `.github/workflows/chocho-island51-host-circuit-atlas.yml`.

- [ ] Report where both gains occur, complete positive/adverse accounting, fixed concentration summaries and selected source-port examples, with all51 records retained.
- [ ] Render and inspect both-target maps and all-island primary-contrast ranking figure. Do not invent island names or causal port effects.
- [ ] Run the complete repository suite, a real anonymous-stage compatibility check and a separate full atlas reproduction; verify numerical outputs and all parent files unchanged.
- [ ] Save a continuation commit/draft PR against the PR128 branch and verify remote CI. Record exact completion evidence separately from this frozen plan.

The user's instruction to continue without unnecessary confirmation authorizes this execution. These method/display choices are analyst decisions, not claimed user approval of individual statistical details.
