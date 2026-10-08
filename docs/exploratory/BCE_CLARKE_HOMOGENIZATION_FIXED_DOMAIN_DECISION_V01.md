# Source-completeness sensitivity of butterfly resource homogenization

**Date:** 2026-10-08. **Status:** COMPLETE / VERIFIED PAIRWISE ORIGINAL-vs-BCE-AUGMENTED FIXED-MARGIN ANALYSIS / POST-HOC STRUCTURAL SENSITIVITY, NOT DEMOGRAPHIC OR CAUSAL INFERENCE. Original GEB submission manuscript and PR #38 unchanged.

## Why this checks the primary GEB biological reconstruction

The original headline that anthropogenic host redistribution homogenizes **potential** butterfly larval-resource assemblages rests on:
1. a large and partly arithmetic raw increase in regional Jaccard similarity, native `0.2769818977` → contemporary `0.4620834005`;
2. a **small** positive excess above a row-and-column-preserving added-resource fixed-margin null (`~0.00860` absolute Jaccard, Monte Carlo one-sided `p=0.002`);
3. potential shared-resource exposure, **not** observed competition or population performance.

An additional independently compiled evidence-filtered source (the Butterfly Conservation Europe, BCE, pages *derived from* Clarke 2024, DOI `10.1002/ece3.10834`) identifies 1,027 candidate accepted-WCVP-ID butterfly–host links absent from the original fixed HOSTS relation among 81 exact European species. This source is geographically selected and can share primary references with HOSTS. The preceding partial one-direction botanical expansion sensitivity changed the absolute added region units `14,553→14,967` while the proportional native-to-added expansion decreased `54.85%→50.14%`. See `BCE_CLARKE_ONE_DIRECTION_GEO_SENSITIVITY_DECISION_V01.md`. Whether the **nonrandom homogenization excess** remained was still untested.

## Frozen comparable analysis

The protocol `BCE_CLARKE_HOMOGENIZATION_FIXED_DOMAIN_PROTOCOL_V01.json` was written **before the new homogenization outcome**. We reran the native/contemporary plant-region reconstruction from:
- original HOSTS commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`;
- WCVP commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`;
- exact accepted-ID candidate links from source WCVP reconciliation [run 37777881774](https://github.com/zuizui0223/chocho/actions/runs/37777881774);
- original frozen 239-species resource metric source and original homogenization receipt `provenance/reviewer_defenses/results/butterfly_resource_homogenization_v0.1.json`.

**Important comparability decision:** fix exactly the **original 355 WGSRPD3 regions containing any native resource**. The augmented native source model supports **356** such regions globally, but the one newly active region is not introduced into this paired comparison. We do not redefine the native-active domain and then call the difference a biological signal. Region overlap and butterfly-resource overlap were recalculated in each scenario on the same 355 regions.

For **each** scenario independently, the original `analyze_butterfly_resource_homogenization.py` fixed-margin null was rerun with 499 permutations, seed `20261007`. It preserves (1) the number of introduced-added resource regions for each butterfly, (2) added butterfly opportunities per WGSRPD3 region, and (3) native source incidences as forbidden added positions. These are **scenario-conditional** nulls with different input margins, not paired causal interventions.

**Successful [GitHub Actions run 37781612222](https://github.com/zuizui0223/chocho/actions/runs/37781612222)**; source code `scripts/analyze_bce_clarke_homogenization_fixed_domain.py`; reproducibility artifact `chocho-bce-clarke-fixed355-homogenization-null-v01`. The exact original four mean Jaccard values (native/current regional and butterfly-pair resource geography) and the original 239-species/14,553 added-unit totals were validated against frozen receipts before the augmented analysis was permitted.

## Numerical result

| Same 355-region, 239-butterfly comparison | Original frozen HOSTS | One-direction BCE evidence-filtered accepted-ID host augmentation |
| --- | ---: | ---: |
| Native regional mean pairwise Jaccard | **0.2769818977** | **0.2993793713** |
| Contemporary regional mean pairwise Jaccard | **0.4620834005** | **0.4850609467** |
| Fixed-margin null median contemporary regional Jaccard | **0.4534793774** | **0.4771550540** |
| **Observed-minus-null excess, regional Jaccard** | **+0.0086040231** | **+0.0079058927** |
| One-sided Monte Carlo `p` (`499` draws) | **0.002** | **0.002** |
| Observed-minus-null excess, butterfly resource-geography Jaccard | **+0.0064287919** | **+0.0063997856** |
| One-sided Monte Carlo butterfly resource-overlap `p` | **0.002** | **0.002** |

Within these same original 355 regions, the fixed HOSTS reconstruction contained 14,121 butterfly × region introduced-added cells; the BCE-augmented projection 14,444. These are **truncated-to-common-domain** values and differ from the all-regions headline 14,553/14,967; never substitute them as the new global totals.

The outcome is notable in a bounded methodological sense: the **small nonrandom resource-homogenization excess persists**, albeit modestly attenuated, despite the 1,027 accepted-ID candidate butterfly-host additions, widespread source-data disagreement, and a substantial change in the inferred native resource baseline. **It is not evidence of realized consumer homogenization, competition or direct demographic impacts**; within-network structural effects remain the measured quantity.

## Scientific interpretation and limits

1. **Structural signal retained.** The excess above each scenario's fixed-margin null is positive and Monte Carlo-resolved. Thus an original-only HOSTS lexicon is not strictly necessary for the sign of this **specific** regional potential-resource structural result.
2. **Most of the raw gain is expected.** The original `+0.1851` raw mean-Jaccard gain dwarfs the `+0.0086` excess. After augmentation the primary excess remains small (`+0.0079`). Do not claim the full raw homogenization comes from special host identity or independent interaction rewiring.
3. **Comparability requires the original 355-region domain.** BCE adds one additional native-active region, and that change is reported separately. On a changing region domain the Jaccard averages are not directly comparable.
4. **Source selection persists.** The supplement modifies the 83 BCE evidence-bearing European butterfly species (81 with additional links) and leaves other original focal butterflies untouched. BCE data use Clarke's evidence grades 1–3 but do not expose individual source/country/tissue; grade 3 egg/oviposition is not proof of adult survival.
5. **The null's estimand differs across scenarios.** Each randomization preserves the margins of its own scenario. Similar `p` and slightly different excess are **not** a paired significance test or proof of ecological effect-size equivalence.
6. **A more geographically conservative challenge remains.** The original manuscript also tested broad geography with WGSRPD Level1 constrained fixed-margin swaps. The source-augmented result has not yet been tested with that stricter continent-preserving null. Do not extend source-robustness claims to all possible null specifications or to observed species interactions.
7. **The GEB paper is still about butterflies.** This is evidence that the way globalized plants are linked to actual larval consumers matters for reconstructing potential geography. No novel local biological host-use/demographic mechanism has emerged from the secondary checklist.

## Decision

**Retain**, with explicit conditional scope, the current GEB claim that plant redistribution simultaneously expands and structurally homogenizes potential butterfly larval-resource geography. This central fixed-margin excess does **not** collapse under the tested one-direction host-source augmentation. **Do not inflate** the main manuscript claim, silently replace the original dataset, or rewrite the paper from the source-augmented numbers. Keep the source sensitivity in the exploratory branch as a reviewer-defense receipt and consider a carefully labeled appendix only after assessing the stricter continent-preserving null and taxonomy/provenance scope.

The submission PR #38 and its successful four workflows are deliberately untouched.
