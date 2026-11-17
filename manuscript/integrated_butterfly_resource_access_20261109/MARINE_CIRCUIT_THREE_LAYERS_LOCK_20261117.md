# ONE GEB manuscript — marine network results and source/claim gate (2026-11-17)

**Single coauthor review branch:** `analysis/butterfly-trait-public51-20261109` ([PR #157](https://github.com/zuizui0223/chocho/pull/157)). **Source branches are archives, not competing manuscript authorities.** Do not bulk merge/close other exploratory PRs. No journal submission or causal marine transport claim.

## Main result: circuit resistance contributes differently to contemporary plant floras

Native versus introduced WCVP accepted-species Jaccard use the **same 104 TDWG3 botanical regions and same 5,355 unordered pairs**, independently coastline-anchored to original GSHHG 2.3.7 marine geometry. All test outcomes are contemporary indexed regional plant assemblages, **not seed arrivals or abundance**. The same original uniform marine network `0.0625°` water mask and `0.5°` finite-volume graph (189,550 giant-component nodes) produced **undirected** effective resistance, using 1–3 verified *nonport* coastal terminals per botanical region. This is neither directional ocean-current velocity nor shipping OD.

Plant models use source-locked ridge alpha 20, standardization from training pairs only, fixed base predictors (great-circle log distance + square, TDWG Level1 equality, latitude gap, original WCVP marginal richness min/max for the *outcome-specific* flora, island area min/max), and append one feature at a time. Each island is omitted, including all its incident pairs, from training; fold errors average equally across **104 held islands**. The held island's *contemporary flora richness is known* in the base model, so this is not prospective floristic prediction. The outcome was already used in exploratory choice of physical predictor comparisons.

| Increment beyond matching baseline | Native WCVP flora MSE gain | Introduced WCVP flora MSE gain |
|---|---:|---:|
| **Uniform-ocean effective resistance (Circuitscape-type)** | **−0.1096%**, descriptive island bootstrap95 **[−0.3906,+0.1652]** | **+4.2390%**, bootstrap95 **[+2.5114,+6.0294]** |
| Shortest-water path excess | +2.2343% | +1.7751% |
| Contemporary WPI port-count product proxy | −0.0771% | +0.0429% |
| Full uniform circuit **after** port proxy | −0.1060% | +3.9723% |
| Port-count proxy **after** uniform circuit | −0.0734% | −0.2356% |

The native circuit negative value means **no additional improvement under this predictor and baseline**, not absence of historical marine dispersal. Shortest-water distance improves both floras, so do NOT claim the native flora is unrelated to marine geography. Port numbers are NOT an observed vessel-origin–destination network. Only 7/7,875 source WCVP pairs have geoverified 2020 AIS OD from the restricted WPI coverage.

### New same-island paired effect-specificity check

Using the original **all 104 paired per-island loss records** from successful [GSHHG circuit held-island CI #44013825017](https://github.com/zuizui0223/chocho/actions/runs/44013825017), source artifact `14298086879` (ZIP SHA256 `16a64ea83d4c0ad9fe6f9961cdc303afc6fce79c65ef3c3efb02f1db44224252`), the introduced-minus-native **relative MSE-gain contrast is +4.3486 percentage points**. Exactly the same held botanical region indices were resampled **20,000 paired draws**, fixed RNG seed `20261117`; descriptive paired bootstrap 95% **[+2.7107,+6.0873] pp**. This explicitly compares normalized improvements on the same regional set, *not* the absolute Jaccard or raw MSE scales across two floras. It is **not a confirmatory between-outcome interaction p-value**: the predictors were explored on the same underlying WCVP outcomes, fold errors share training pairs, and the resampling is not independent regional replication.

The original null of 104 island geometry circuit was executed source-verifiably [#44013658955](https://github.com/zuizui0223/chocho/actions/runs/44013658955); original shortest-water source and WCVP taxon source are [#44010673580](https://github.com/zuizui0223/chocho/actions/runs/44010673580), [#44010011288](https://github.com/zuizui0223/chocho/actions/runs/44010011288). The exact original model is [#44013825017](https://github.com/zuizui0223/chocho/actions/runs/44013825017). This GEB branch contains a copied full source script and direct source workflow; the original individual-fold result ZIP remains a source artifact, not fabricated in the repository.

## Butterfly layer: the same *uniform-ocean circuit method*, but another sample and baseline

The original 51 selected landmass groups / 1,275 unordered pairs target the **2011–2022 GBIF first-indexed butterfly assemblage Jaccard**, not the WCVP floras. The predictor comes from the **same original GSHHG finite-volume marine source**, but has **182 source-validated WPI port terminal weights** rather than up to three botanical coastal terminals per TDWG region, so it is not the numerically identical 104-island matrix. The baseline includes geographic distance and its square, port-count sums/differences, and **earlier butterfly (2001–2010) Jaccard**; OLS not plant-model ridge. One held island's 50 incident pairs is removed from model training; 33/51 held islands improve.

Uniform-ocean circuit MSE gains **+5.081886845%** over the separately specified butterfly baseline (MSE **0.00679776077 → 0.00645230626**). The archived original 51-fold loss CSV ([source](../../results/island51_observed_ocean_circuit_20261106/all_island_holdout_losses.csv), blob SHA `98177bf46a2274bea68015b87a1ccb7b6fc7f32c`) was read and its 51 butterfly base→uniform rows independently aggregated; 20,000 fixed-seed (xorshift32 20261117) descriptive held-island bootstrap draws gave **95% approximately [+1.36,+8.30]%**. This interval is a supplementary visualization only and not a hypothesis-confirming independent 51-fold confidence interval.

**Cross-layer synthesis:** The two models show that undirected marine network resistance adds predictive information to recorded introduced flora and butterfly composition, whereas it does not improve the native-flora model in the matched botanical analysis. They **do not** estimate a plant→butterfly mediation effect; the 51/104 island populations, response definitions, terminal allocations, baseline controls and estimator families differ. Do not claim that identical circuit *values* were fitted to identical plant and butterfly pairs or that current direction transports the organisms.

## Figure 3 and editorial interpretation

- Main candidate: [Fig3_circuit_three_layers_v01.svg](../../figures/integrated_island_port_butterfly_v03/Fig3_circuit_three_layers_v01.svg), Panel A shows three separate relative prediction changes with descriptive bootstrap intervals; Panel B holds the *same 104 botanical islands* constant and distinguishes circuit, shortest-water, and contemporary port-number models.
- Correct headline: **"Introduced floras and later butterfly records share additional predictive structure associated with marine network connectivity"**; no "currents caused homogenization" or proven shipping/cargo route.
- The NOAA old-OSCAR 2005 directional-particle source (12 model coastal-approach directions and only 4 other-region coarse land-mask contacts) belongs to exploratory Supplement only. It does not validate viable seed landing. The EuDiS thalassochory 22-island composition difference is independent morphology evidence, not natural transport history.
- Original resource and butterfly positive findings remain unchanged: potential resource geography +54.9% global; 239 butterflies and +118.9% island relative gain; 155 known-butterfly GBIF record forecast MSE +4.12%; 51-island host-vs-butterfly similarity rho +0.300. No causal arrival or successful breeding inference from GBIF first indices.

## Reviewer-ready source and priority

**Canonical source choice:** this file plus the [single manuscript](ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md), [figure captions](FIGURE_CAPTIONS_V03.md), the [machine results ledger](../../results/port_weed_current_review_20261117/ONE_PAPER_CIRCUIT_RESULTS_LEDGER_20261117.json) and the copied scripts on this SAME branch. Keep original exploratory PRs (including #175) as source archives, not as independent papers or replaceable latest numbers.

**Remaining blocking claim check:** 104-island matched circuit contrasts and 51-island butterfly source have separate populations. A stricter botanical **both-island-endpoints-unseen** 5-fold test exists for *shortest-water geometry* but has not been verified for the source-equivalent uniform circuit; DO NOT mislabel the original +4.239% as a strict both-unseen circuit result. Likewise, future true current/cargo historical tests require better directed data and independently dated floras. Do not claim PR merge, submission, coauthor approval, or independent three-community causal mechanism.
