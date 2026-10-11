# Corrected core resource ecology — v04

Date: 2026-10-11. Retrospective sensitivity, not preregistration, new field evidence or a complete submission replacement.

## Decision

The v03 hybrid-marker reconciliation does not overturn the core reconstructed results. On the fixed 239-butterfly panel, resource expansion is 55.2444%, regional resource-assemblage similarity increases by 67.2113%, and exact shared-host butterfly-pair × region exposure increases by 128.6444%. These are resource-geography quantities, not observed range expansion, community homogenization or competition.

## Corrected observed quantities

Native and contemporary butterfly-region totals are 26,533 and 41,191, giving 14,658 added cells; 206 species expand. On the same 355 original native-resource-active regions, mean regional Jaccard is 0.277019 to 0.463207. Mean butterfly resource-envelope Jaccard across all 28,441 pairs is 0.222409 to 0.328707. The full source grid has 367 region columns, but only 365 contain focal contemporary resources; ANT and HMD are zero-support columns. These domains must not be conflated.

Exact accepted-host-ID shared geography increases from 62,473 to 142,841 pair-region units, compared with 141,885 contemporary units before correction. There are 996 host-sharing pairs, 930 of which gain shared-resource regions (originally 986 and 920). Across-family pairs account for 42,051 of 80,368 novel shared pair-region units, 52.3231%. Recovered records are not newly evolved or newly observed interactions.

Concentration remains: 38 of 677 contributing taxonomic IDs account for 50.1703% of fractional added-resource credit. IDs may be infraspecific, so they are not necessarily 677 independent plant species. The single-ID-supported share is 58.2208%; grouping to canonical species labels only inside supported cells gives 58.7051%. The latter grouping imports neither unrecorded form feeding nor parent-species distribution elsewhere.

## Ecological nulls rerun, not inherited from the old result

Both global and continent-constrained butterfly-resource incidence nulls were rerun on original and corrected inputs. All retain native resource cells as structural exclusions and preserve regional addition counts. The continent null additionally preserves each butterfly's added-region count within continent.

| Variant and null | Observed regional Jaccard | Null mean | Excess |
| --- | ---: | ---: | ---: |
| Original, global | 0.462083 | 0.453470 | +0.008613 |
| Corrected, global | 0.463207 | 0.454622 | +0.008585 |
| Original, continent | 0.462083 | 0.458469 | +0.003615 |
| Corrected, continent | 0.463207 | 0.459618 | +0.003589 |
| Corrected, continent, independent seed and doubled mixing | 0.463207 | 0.459621 | +0.003586 |

Each run has 499 constrained-chain states. Upper-tail Monte Carlo p is 0.002 in all cases, the minimum attainable with the +1 correction. The corrected continent chain's central 95% null range is 0.459520–0.459712; this is not a causal confidence interval. The long-chain replicate gives a similar mean and lower lag-1 dependence. Margins and native exclusions were verified, but complete state-space mixing is not proven.

The absolute excess is small relative to the raw native-to-contemporary increase. Do not describe the 67.2% increase as an equally large excess over null expectation. These tests concern homogenization of butterfly-resource incidence, whereas the earlier climate/seasonality models concerned introduced plant alignment with native-host resource envelopes. The present result does not overturn the earlier bounded, very small saturation residual.

## Verification

An independent pandas source-link join and integer SciPy reconstruction reproduced the original and corrected totals, exact and species-grouped sharing, contributor credits and fixed-domain Jaccards. Every saved randomized matrix was unpacked and checked: 2,495 states and 4,990 independently recomputed Jaccard scores, including row/column/continent margins, structural exclusions and file checksums. The full local package passed 147 tests with no failures or skips. Three standalone PNG/SVG comparison figures and tables were regenerated and visually checked.

Two numerical bookkeeping safeguards were added during verification: distinguish the 367-column source grid from 365 focal occupied regions; and use exact rational expansion ratios for descriptive diet-breadth ranks so floating-point log subtraction does not split mathematical ties. Corrected exact-tie Spearman rho is 0.009613; the legacy scalar-log result is 0.009755. No zero-effect or equivalence claim is made.

## What remains before replacing the submitted numbers

BCE/Clarke source augmentation (S9), crop/Poaceae exclusions, occurrence validation, phylogenetic/equivalence tests and plant-level host-bias/climate/seasonality nulls have not been rerun with the corrected host set. The three figures are audit comparisons, not the complete final manuscript/SI/anonymous package. GEB PR #38 remains at 88209b8a215ac1bc9c4c5e67622fa104179c8895. Do not swap only 54.9% for 55.2% while leaving incompatible downstream material.

The result receipt is ../results/butterfly_corrected_core_ecology_v04_20261011.json. Full code, input sources, five raw traces and packed matrices, reproducibility commands, comparison tables and an incremental v03 patch are in the delivered v04 bundle. Receipt and note commits do not imply the entire local code package was merged into this remote branch.
