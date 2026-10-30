# Exact conditional source-distance and original host-roster implementation plan

> **For agentic workers:** use the existing parallel implementation and independent review assignments. Follow each test-first step before new production code. User-authorized continuation; no additional design approval is required.

**Goal:** independently resolve numerical uncertainty in the frozen recorded-source comparison and recover the original host-edge roster for all230 conservatively bridged butterfly taxa.

**Architecture:** keep the complete20261029 analysis, its inputs, sampler and failed qualification result unchanged. Factor the same fixed-margin binary-table space into its four longitude blocks, count completions with exact integers and draw independent uniform tables by random-rank unranking. Separately derive a bounded230-taxon host roster from the SHA-pinned author archive inside GitHub Actions; export derived scientific fields only.

**Tech stack:** Python3.12, standard-library integer arithmetic/gzip/hashlib, NumPy/SciPy for checks and binomial intervals; existing GitHub Actions source-restoration route.

**Spec:** the accepted next steps in `ISLAND_BASELINE_HOST_RESULTS_20261029.md`, and the unchanged scientific estimand in `ISLAND_BASELINE_HOST_PLAN_20261029.md`. Starting repository commit: `417f6f733e7412634664b155a0f2f1fed0ea87c0`. This new plan is saved before scored exact-DP results or new random draws. Prior MCMC contrasts and unscored mask/margin feasibility are already known; this is **post-pilot numerical validation**, not a new untouched ecological validation.

## Global constraints and source identities

- Preserve every prior scientific input, plan, first-run receipt and sampled outcome byte. Do not extend the old chains, relax their gates, modify the risk set or acquire late outcomes for the stopped120-island roster.
- Use the primary and two previously declared scenarios with exactly the same species, islands, structural zeros, species-by-longitude-block margins, island margins, source-record thresholds and upper coastline-distance statistic. Mainland and unsampled sources remain omitted. Early nondetection is not absence; later recording is not arrival or establishment.
- Original recorded-source diagnostic SHA256: `372fec4a1035132c17eec67c5392c5eab0727d3c596edb2d4f2b45b88b6ecebb`; risk cells `a8c7905f1ef6e187b2647d5ed39f1a5a67f901b94a5aab3a86b6358969f72fb4`; fixed margins `e472a2dcb0977a02dda1f313bc66ab32d0c096a2b74d505c764791b1b15ecbd0`; coast pairs `d83803dd3f976069533da40ff1c20f6602386b699de3dd744f37096e63f03545`; first-run lock `a0124d968516c1d9a5b5cffb1ba6d6dd7f5826cedccd81c188c9d2e295735619`.
- Existing conservative230-name bridge SHA256: `cd8b2cb8e0279a5b940e46267c67f0bb9843d2e1c54f38bc405ec30a95fe6eb2`. Expanded96 species facets SHA256: `061f9ce3cf21e4d64572456e0ffbbb88022fbccfe9aa924d5f02062027c97db6`.
- The original checkout contains a truncated uncommitted facet CSV; preserve it. All work uses an isolated checkout of the verified starting commit, where the complete CSV matches its published hash. Never pass the truncated file to analysis or silently include it in a commit.
- No new biological-coordinate acquisition or persistence in this continuation. Do not export the raw HOSTS/WCVP archive, tables or unrelated taxon data. No source archive is downloaded into this local checkout.

## Review focus

1. Exact integer completion counts and unranking must be a bijection over **labeled** feasible tables, including structural zeros and forced cells; do not treat row/column identities as exchangeable without multiplicities.
2. Independent draws are uniform **with replacement** from the whole finite table space. Do not deduplicate ranks, use floating branch weights, condition proposals on validity, or stop when a desired tail count appears.
3. Zero tail hits have a nonzero uncertainty bound. Preserve the old MCMC failure and explain why a valid binomial bound is available for the new independent draws.
4. Original230 host selection uses the conservative taxonomy bridge, with all9 excluded original concepts accounted for. Baseline93 coverage uses2010–2017 only; neither105 ever-recorded taxa nor102 late taxa replaces that denominator.
5. Source identity, complete count reconciliation and safe derived-only outputs must gate host export; missing/conflicting concepts cannot silently become zero, and failure produces an explicit incomplete status.

## Task1 — exact finite table counts, independent sampling and numerical bounds

**Files:** create `scripts/analyze_island_source_exact_null_20261030.py`, `tests/test_island_source_exact_null_20261030.py`, and `data/exploratory/island_source_exact_null_20261030/`. The exact-null worker owns these until handoff.

**Interfaces:** `ExactBlockFiber` consumes a labeled binary mask and fixed row/column margins and exposes an integer completion count and deterministic `unrank(rank)`; the analysis consumes the pinned old risk/margin/diagnostic files and emits block counts, cell occupancy counts, conditional moments, every sampled rank/statistic and qualification receipts.

- [ ] Test exhaustive tiny feasible fibers, every rank's unique table, empty/forced/zero-margin cases, block-product rank bijection, infeasible margins, resource caps, score-independent table law, toy moments and fixed-seed replay. Watch failures before implementation.
- [ ] Implement suffix dynamic programming `C(i, residual_columns)`, summing exact integer child counts over lexicographic allowed row subsets of the fixed row sum. Keep column labels; order active rows by number of allowed subsets, allowed degree and original index. Zero-margin deletion is a reversible reduction. No-negative, residual-sum and suffix-capacity pruning may remove only impossible completions.
- [ ] Cap at100,000 memoized states and2,000,000 attempted row-pattern transitions **per block**; any exceeded cap fails explicitly. Unscored feasibility found at most6,401 states and67,074 transitions, so the caps are conservative but fixed before scored execution.
- [ ] The global fiber is the Cartesian product of four block fibers. Draw one unbiased arbitrary-integer rank uniformly from `[0, product(block_counts))`; split into block ranks with integer division and unrank each block. Use **100,000 independent full-table draws per scenario**, with replacement. Fixed scenario order/seeds: primary2026103001, source_records_ge3=2026103002, late_records_ge2=2026103003. Do not adapt budgets after seeing a result.
- [ ] Compute conditional mean/variance by stable centered-mixture recursion with exact integer completion counts and floating score arithmetic. Independently obtain cell occupancy counts from forward prefix multiplicities times suffix counts; reconcile expected row/column sums and the mean statistic. Report forced and variable cells and descriptive block contributions; no additional ecological hypothesis test is introduced.
- [ ] Preserve the original mean log(1+upper coast-gap) statistic and inclusive closer-source tail, with the existing1e-12 numerical tie tolerance. Check the observed recomputation against the old receipt. Validate every sampled table's mask and margins. Store every rank and statistic in deterministic compressed CSV, plus input/script/plan hashes and actual run details.
- [ ] Report the empirical tail count and a **one-sided95% exact binomial Clopper–Pearson upper bound** for its Monte Carlo uncertainty; for zero hits use `-expm1(log(.05)/N)`, never report probability zero. Also predeclare a simultaneous three-scenario95% upper-bound family using per-scenario alpha=.05/3. Bounds concern uniform conditional-table tail probabilities, not ecological effect sizes. No old20-tail-hit gate is transplanted to these independent binomial draws; justify the changed numerical inference in the report.
- [ ] Qualify only after all exact-count/bijection/source/margin checks and all fixed draws complete. Incomplete execution produces no qualified bound. Sample mean/variance versus DP moments are diagnostics, not a reason to choose a different specification or extend the budget.

## Task2 — recover original230 host edges inside the established source runner

**Files:** create `scripts/derive_original230_host_roster_20261030.py`, `tests/test_original230_host_roster_20261030.py`, `.github/workflows/island-original230-host-roster.yml`; retrieve derived results into `data/exploratory/island_host_roster_20261030/`. The host worker owns the new script/test/workflow; root publishes and retrieves the derived artifact.

**Inputs:** successful source run39006554609, head `d53e04fb361701265139ba77b8bca56689da6532`, artifact12090584135. Outer ZIP digest `a784770f5a215bffb32393c127b66e7320e8ed7ac9b36f921bb9c056f3497a36`; inner `butterfly_resource_novelty_anonymous_OFFLINE_source_v01.zip`, size3,255,228 bytes and SHA256 `408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323`. Original `inputs/original/insect_host_accepted.csv` SHA256 `0a084fb5273e4780b03e015205c87e4545c69f0f8e517feb1a2a3879889508e9`; original239 metrics SHA256 `7b2a20387d656dbfcd7b3c38ec6a7fe2505474ad51f783dd1f201eb6b7eabd30`.

- [ ] Test source/count drift, accepted-ID duplicate and name ambiguity handling, taxonomy collisions/exclusions, early-only coverage, archive tampering/unsafe paths and the exact derived-output allowlist before implementing.
- [ ] Restore the pinned artifact in GitHub Actions using the existing authorized download route; verify source run/head/artifact identity, inner SHA/size and source-manifest member hashes. Read the host table from the verified ZIP without extracting raw tables. Use the original accepted HOSTS table only, without the BCE-augmented old-country edge subset.
- [ ] Reconcile all239 original per-species unique accepted-host-ID counts against pinned metrics. Expected aggregate counts from existing summaries are2,596 all239 and2,524 conservative230; treat them as checks, not substitute recovered edges. Resolve230 exactly from the bridge and explicitly account for9 exclusions.
- [ ] Export only `original230_host_edges.csv`, `original239_host_count_alignment.csv`, `plant_taxon_query_candidates.csv`, `baseline93_host_coverage.csv`, and `original230_host_roster_manifest.json`. Fields are scientific butterfly identities, current concept keys/status, accepted WCVP plant IDs/names/families, count/provenance/ambiguity flags and hash/validation metadata. Do not copy raw source tables or excluded-taxon edges. Preserve conflicting names or IDs as explicit unresolved records, never inferred botanical links.
- [ ] Derive baseline93 exclusively from early expanded96 facets; compute recovered catalog-edge coverage. Plant candidates are a query inventory, not GBIF-resolved records or proof of plant presence. No new occurrence-source queries in this task.
- [ ] Use a narrow push trigger on the current PR102 branch and this workflow/script/test/plan paths. Read-only repository/actions permissions; tests before source restoration; upload only the named derived outputs. Fail closed on source/schema/count errors.

## Task3 — integrate, verify and publish the scientific advance

- [ ] Root and a separate reviewer reconcile exact table counts, selected rank-to-table reconstructions, tail bounds, cell/block support and the distinction between conditional recording structure and dispersal. Independently verify host output counts and source hashes.
- [ ] Save a report and useful numerical figure, explaining that the same conditional comparison is now evaluated without MCMC mixing uncertainty. Preserve the original failed run as historical evidence. State remaining climate, recording and source-universe limitations.
- [ ] Exclude all new exploratory scripts/tests from the original reviewer bundle; run the complete project test suite and verify the original222-file bundle identity. Publish to draftPR102 and check CI on the actual new commit; do not merge the manuscript branch.

**Scientific stop rule:** if exact numerical/source checks fail, document the specific failure without claiming the distance hypothesis is confirmed. If host recovery fails, keep the concrete workflow and failure evidence, complete unaffected exact-null work and identify the unresolved source dependency. Do not replace either task with literature-only review or a newly selected favorable model.
