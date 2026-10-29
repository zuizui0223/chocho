# Baseline support, recorded sources and island host co-recording — plan

**Goal:** resolve the demonstrated far-WPI baseline bottleneck, examine a species-specific recorded-source mechanism in the already observed panel, and acquire actual point-to-island host co-recording evidence.

**Starting state:** repository `757d7568b4f57ac074c240f78441cfb20eb3b09f`, draft PR #102. The original 96-island richness/port models remain not fitted under their immutable plan. This continuation is user-authorized. It does not alter the original GEB manuscript or the previous eligibility decisions.

## 1. New baseline-only audit: exactly 120 unseen deficient-cell candidates

Before any new outcome, select every original source island satisfying: area100–20,000 km²; representative latitude within±55°; no Antarctic ambiguity; continental coast-gap upper bound<100 km; WPI-distance lower bound>100 km; and absent from the already observed96 roster. The geography census yields117 medium and3 large islands. Freeze identities, exact geometry hashes, original source hashes and this plan hash before querying.

Acquire **2010–2017 only**, using the same COL XR Papilionoidea concept, exact unsimplified polygons and complete species/year/dataset/basis facets. Retain the separate first-page check, original-byte receipts, complete-zero versus failed-query distinction and no island replacement. All120 remain in outputs. This is a baseline recording feasibility audit, not a survey or biological absence test.

Apply the previous baseline recording screen: >=20 species-assigned records, >=90% species assignment, >=3 recorded years. If fewer than3 candidates pass, stop this new port-cohort route without acquiring its late-period outcomes. If >=3 pass, inspect their fixed geographic components and counterpart availability before fixing any subsequent cohort; passing this cell-count threshold alone is not sufficient for the full model. No late outcome informs this new eligibility decision. Request starts remain throttled; the baseline and point-acquisition clients run sequentially to avoid competing for the same public source.

## 2. Separate exploratory test: proximity to prior recorded sources

The96 early/late outcomes have already been seen. The following is explicitly **post-pilot exploration**, not an untouched temporal validation or a replacement for the failed richness/port model.

Primary taxa are all230 conservatively bridged original names, restricted to those recorded on at least one selected island in2010–2017. A candidate species–target-island cell has no early focal-species record, an early record on another selected island, >=3 total early recorded species at the target, and >=4 total late recorded species. The last condition supplies at least3 OTHER late species whether the focal label is0 or1, avoiding a permutation-dependent eligibility mask. Targets and controls are not called surveyed absences. The outcome is at least1 late indexed record of the focal species.

Measure minimum exact source-coast gap from each target to the selected islands with an early focal-species record, excluding the target. Use the original <=1 km source-coast sampling and retain the discretization interval. The sampled96 source universe omits mainland and many island populations; these are **recorded source candidates**, not verified breeding populations or demonstrated sources of dispersal. Do not substitute representative-point distance.

The primary statistic is mean `log(1 + nearest recorded-source coast-gap km)` over positive later-record cells. The directional question is whether later records are closer to earlier recorded sources than expected while preserving both (a) each species' number of later records **within each fixed90° longitude block** and (b) each target island's number of later records. These margins control species/block recording gains and target recording propensity. They do not control all climate, habitat, observer or source-sampling bias. Island-constant port/area/mainland variables cannot be estimated from this fixed-island-margin test.

Use a uniform fixed-margin binary-matrix null on the allowed risk mask. A reversible proposal chooses distinct existing positive edges and cyclically permutes their columns only if all proposed cells are allowed and currently zero; include every possible cycle length, so structural-zero fibers are not limited to2×2 swaps. Sample on **attempted steps**, including rejections, not only accepted moves. Verify detailed balance and connectivity on a small enumerated fixture. Four fixed seeds,100,000 attempted burn-in steps each, then4,000 samples per chain every50 attempts. Report moves/acceptance, independent-chain stability, split-Rhat and effective sample size. If nondegenerate sampling is not supported, split-Rhat>=1.05 or pooled ESS<1,000, do not report a qualified tail probability. The MCMC tail area is an approximate conditional randomization result, not an exact test or a causal effect.

Primary inference uses the conservative original panel and the region-conditioned margins above. Two declared robustness descriptions use source-island support>=3 early focal records and late focal outcome>=2 records, respectively; keep all other rules/margins and source labels explicit. Do not search alternative specifications for a favorable result. This diagnostic can justify further source/host acquisition, not a colonization or establishment claim.

## 3. Fresh point acquisition of the original host-linked cases

Keep the previously source-locked26 butterfly–country cases,44 plant–country query definitions,46 candidate edges,16 butterfly taxa and24 plant taxa. Use their original numeric GBIF keys in the explicitly declared old Backbone checklist, with2010–2025 filters. This is a new retrieval; source count changes are recorded, never silently equated with the earlier snapshot.

Use count probes and paginated point reads, bounded at5,000 records per query,12,000 rows total and200 actual HTTP attempts including retries. Any cap, source inconsistency or failure is incomplete, not zero. Each query keeps request identity, source-body hashes, counts and time/dataset diagnostics. Raw exact coordinates and record bodies remain **only in transient process memory**: none in repository, cache, logs or output artifacts. This new transient-local pipeline leaves the old Actions-only script/workflow unchanged and preserves the original protocol's no-raw-persistence boundary.

Primary spatial output is exact membership in the fixed96 source polygons. Separately inventory membership among all8,260 source islands>=5 km², using the same acquired points and no extra source queries. That inventory describes the original positive-country case cohort; it does not replace the96 roster or support representative global inference.

Distinguish coordinate-point membership from a separate conservative containment screen for the full **reported** uncertainty envelope. Unknown uncertainty stays unknown. Preserve source polygon holes, boundaries and dateline behavior; do not claim that reported uncertainty or source shorelines are perfectly accurate. Report same-island butterfly–candidate-plant records by early/late window and an absolute year-difference<=3 rule. Temporal cross-window pairs remain explicitly separate. Co-recording is not feeding, successful reproduction, first arrival or host-mediated expansion.

## Implementation and review

1. Root implements/runs `audit_unseen_island_baseline_20261029.py` with source/selection/eligibility tests; outputs `data/exploratory/island_baseline_20261029/`.
2. Recorded-source worker implements exact distances, immutable risk-set accounting and a verified fixed-margin null in `analyze_island_recorded_sources_20261029.py`; outputs `data/exploratory/island_recorded_sources_20261029/`.
3. Host worker implements/runs `acquire_island_host_corecords_20261029.py`, numerical/spatial/cap tests and derived-only outputs in `data/exploratory/island_host_corecords_20261029/`.
4. Independent review checks source identity, missingness, sampling semantics, geometry and interpretation. Root verifies actual outputs, updates the current report/PR, excludes these scripts/tests from the original reviewer bundle, and runs tests/CI.

The broader original-host roster remains a known source dependency: the live, SHA-pinned author artifact can supply catalogued host edges through the established runner-only derivation route. Local summaries cover only7 of the105 island-recorded original taxa. That broader recovery is not silently substituted with inferred botanical associations in this continuation.
