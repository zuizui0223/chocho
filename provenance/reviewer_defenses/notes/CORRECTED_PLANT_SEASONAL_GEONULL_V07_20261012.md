# v07: Crossed host-source plant-level seasonal and geographic null

Date: 2026-10-12 JST. **Exploratory, post-hoc; not a revision of the frozen GEB paper.**

## Question and controlled comparison

Can the observed overlap of introduced known larval host plants with butterflies' pre-existing native-host resource regions exceed chance when the same plants are permitted to move only between comparable regions? Four host-resource reconstructions use the same *botanical* null pool and the same sampled shuffled plant matrices:

1. Original literal HOSTS taxon names.
2. HOSTS plus unique WCVP hybrid-marker reconciliation.
3. Original HOSTS plus the fixed, Europe-biased BCE/Clarke host candidates.
4. Hybrid-reconciled HOSTS plus the same BCE/Clarke candidates.

The native plant ranges and host-region interaction lists per source remain fixed. The null conserves each plant's introduction count within family × climate × subcontinental blocks and each region's number of introductions in each family. BIO1/BIO12 climate tertiles, BIO7/BIO15 seasonal-proxy halves, WGSRPD2 are applied in common, using previously checked frozen climate and mapping data. Real seasonal tissue availability, larval development, adult dispersal, and demographic colonization are not measured.

## Fixed 191-butterfly primary quality-panel result

| Host data | Observed share of introduced host incidences inside native resource geography | Expected null share | Conditional excess |
|--|--:|--:|--:|
| Original HOSTS | 62.505% | 62.429% | +0.075 percentage points |
| Marker-corrected HOSTS | 62.294% | 62.216% | +0.077 percentage points |
| Original HOSTS + BCE | 67.266% | 67.167% | +0.099 percentage points |
| Marker-corrected + BCE | 67.055% | 66.956% | +0.099 percentage points |

Different seed, twice the burn-in and between-sample mixing: marker-corrected+BCE excess **+0.0997 percentage points** (primary seed +0.0994). Both have 199 sampled Markov-chain states, so the upper-tail p=0.005 is the possible minimum and cannot be interpreted as a large demographic/causal effect.

**Source-completeness contrast:** Adding BCE to the marker-corrected portfolio increases observed pre-existing-resource overlap share by **4.7615 percentage points**, but the same fixed botanical shuffles predict increases of **4.7393 and 4.7365 points** in the two seeds. The extra difference is just **+0.0223–0.0250 points**; paired upper-tail p=0.09, 0.075. This does not demonstrate a new independent sorting effect.

## Selection-dependent new-resource count

In all 239 species the source correction changes unique newly supported cells 14,658 (marker) → 15,072 (marker+BCE). In the original frozen 191-quality subset those cells instead change **13,634 → 13,302**, as new records reveal native resources in previously introduced-only cells. Neither is a measured biological range gain or loss. Europe-biased supplementary host records are not a new globally representative ecology sample.

## Validation and limitations

- **9,280** botanical accepted taxonomic IDs in common, **78,751** plant×region introduced edges, **36,439** in eligible swap blocks.
- Exact family mapping for the 689 WCVP IDs lacking complete local family metadata was retrieved directly from frozen WCVP commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`, with no guessed genera. [Source acquisition run](https://github.com/zuizui0223/chocho/actions/runs/38150410940).
- Two chains: 199 draws each; **1,237,297** and **2,474,312** accepted switches. Final focal edges moved only about **6.14%** and **6.21%** from their original positions because of tight constraints. This is a conditional comparison within a narrow swap space; multiple source variants are not independent replicates.
- Independent SciPy integer matrices reproduced all 4×3 scenario/panel observed scores. Null traces checked for **4,776** panel×draw incidence-accounting identities and **960** comparison-statistic values across seeds. Prior v06 suite 169 tests passed and v07 scoring tests 4 passed.
- The complete, reproducible v07 ZIP with the exact v06 source archive, WCVP family artifact, both traces, code, tests, audit and SHA256SUMS has SHA256 `2f00b7004304e348b94dd767d309528a6f97acf25c176ed5a1eb888b869a9626`.

**Scientific decision:** A tiny null-relative overlap excess persists after multiple source corrections, but the strong additional ecological mechanism is not established. Do not promote the effect in the frozen manuscript, over-interpret p-values, or treat the regional potential-host reconstruction as measured butterfly breeding. The next strongest discriminating evidence would be verified, temporally followed larval development on introduced known host plants in regions lacking recorded native hosts.

The frozen [GEB PR #38](https://github.com/zuizui0223/chocho/pull/38) was not modified. This receipt documents independently executed post-hoc source-condition sensitivity and does not mean v06/v07 cumulative local code has been merged into the analysis branch.
