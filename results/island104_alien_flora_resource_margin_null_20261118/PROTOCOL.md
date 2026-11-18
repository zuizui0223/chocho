# Calibrate source-overlapping alien-flora → butterfly-resource prediction against exact gain-margin null

**Fixed exploratory protocol, 2026-11-18, after the +26.363% (E) / +33.768% (S) direct bridge was seen.** Do not call this preregistered or independent. Source WCVP native/current plants and fixed HOSTS 239 larval-host assignment are reused, so a large observed geographical predictive association is partly built into source definitions.

## Immutable sources and outcome

Reopen original PR #118 commit `bb4d8f00fbdb97fa2d90f606a8c9e7e0620646f2` exact 239×369 bit matrix, SHA256 `ec6a5530d65ac0b04f230e0c51fd68bf543555641a17b94795996bf11ce257d0`, native=26,530, added=14,553, current=41,083 resource cells. Frozen original 104 GSHHG-connected WCVP island regions and 5,355 plant pairs from CI #44010307551/#44010673580. Exactly 2 island pairs have undefined current **butterfly potential resource Jaccard**, leaving 5,353 paired source rows. Both source WCVP whole-alien-plant Jaccard and baseline native resource Jaccard and margins are held fixed. Exact five preexisting SHA256 island-hash folds produce 1,088 double-unseen island test pairs. Never optimize their membership, alpha, exposed overlap or source subset after seeing null.

## Counterfactual structure (no real biology generated)

Fix all original native butterfly resources; preserve exact added-resource count of **each butterfly×Level1×island/mainland/mixed class**, exact added count in each botanical region, exact added count of each butterfly globally, and forbidden overlap with native bits. Generate alternative potential-resource assignments using original source `swap_new_resources` 2×2 transitions (same code as original global restricted null), with `72,765` swap proposals per retained draw. Class roster is fixed main `results/island-mainland-resource-q1/region_metrics.csv`. The 104 island subset is never permuted independently; changes in other 369-region nodes are allowed as dictated by original exact-margin network. This is a null for **placement of added opportunities**, not a null for the WCVP plant composition or observed butterfly movement.

Three independent PRNG chains with fixed seeds `20261119`, `20261120`, `20261121`, each `5×72,765` initial burn-in proposals plus **75 retained states** separated by `72,765` proposals, total **225 correlated draw scores**. Report actual accepted transitions and between-chain range; source nodes/dyads are not iid. If insufficient legal swaps or source counts drift, fail closed rather than silently loosening controls.

## Original prediction model fixed from prior bridge

- Primary score **S**: 5 source-hashfold equal-weight retrospective MSE reduction % from adding WCVP **whole introduced flora Jaccard** to distance/Level1/latitude/area, modern port index, native WCVP flora Jaccard, native butterfly resource Jaccard, original native and alien flora richness, native butterfly resource richness, and **source-current butterfly resource richness margins**. This last feature is derived from potential-resource response and fixed by exact source null; it is **not available in a prospective new-region forecast**.
- Secondary score **E**: same source geography/plant/port/native controls, *without* current butterfly resource richness.
- Fixed source ridge alpha20, train-only standardisation, predictions clipped to [0,1], same 1,088 source-held-both-endpoints test pairs, same original plant/survey universe and original 104 island coastal subset.
- For each null draw compute 239-bit current potential-resource Jaccard for exactly the original source dyads and the E/S incremental fold gains. Compute 225-draw median, central 95% **descriptive** quantiles, max/min and empirical high-tail fraction vs observed. **Retain original positive, negative and zero differences**, no direction-dependent selection.
- Also report source-overlapping raw tied-rank Spearman of alien-plant Jaccard against current-resource Jaccard as a descriptive calibration, not the primary test.

## Decision and interpretation

If observed S gain is large but **typical of the restricted exact-margin null**, the +33.8% is predominantly a construction/geographic common-source signal, not evidence that real butterfly communities track traded plants. If observed S is above this null, it only shows that the exact placement of resources in the existing HOSTS×WCVP network aligns more with whole-alien-flora composition than random placements with the same species/region/Level1×island margins. This remains partly **structural** because the response and predictor share WCVP taxa and the null cannot preserve every host-family, phylogenetic or climatic structure. Do not call it causal, independent, or verified historical plant and butterfly migration.

The external observed GBIF 224 species control has only **+0.157%** MSE from host composition at the 169-region global panel and uncertain intercontinental transport. The island-endemic *Pieris cheiranthi* larval photographs are direct local use for a single species, not a colonization effect. Any source result must be reported alongside these limits, and cannot erase the original +54.9% global/+118.9% island resource-opportunity positives.
