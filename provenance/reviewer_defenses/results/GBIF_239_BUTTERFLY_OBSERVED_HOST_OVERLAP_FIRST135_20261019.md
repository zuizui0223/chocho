# Real butterfly observations versus redistributed larval host-resource geography: first global GBIF cohort

**2026-10-19. Status: SOURCE-VERIFIED LOCAL OFFLINE RECOMPUTATION, INITIAL 135-OF-239 COMPLETED TAXA, NOT AN INDEPENDENTLY COMPLETE GLOBAL EFFECT.**

## Ecological question
Do observed butterfly presences coincide with added larval HOSTS/BCE resource-region opportunities more than expected after broad geography, the climate of originally native hosts and an internal background observer-effort proxy are held comparable? This is macroecological **sampled occurrence association**, never local larval use, reproduction or evidence that introduced plants caused butterfly ranges to increase.

## What actually happened
- New original GBIF taxon-keyed public search, years **2010–2017** and **2018–2025** (systematically separated max 2×300 records per window, max 1,200 sampled records/species). Original workflow [#39233336804](https://github.com/zuizui0223/chocho/actions/runs/39233336804) completed taxon/data retrieval for **all 239 focal species attempted**, yielding **135 complete**, **96 PARTIAL_FETCH_NOT_ELIGIBLE**, **8 TAXON_UNRESOLVED**.
- Original source artifact [#12178748642](https://github.com/zuizui0223/chocho/actions/runs/39233336804/artifacts/12178748642) preserved the audit and **123,764** valid GBIF sampled occurrence records (from complete taxa only) across **6,411** species×WGSRPD3×period×basis aggregates and **268** mapped distinct regions. There were no raw coordinates in the artifact.
- First 239-run model stage failed on a **code interface error**: `source_augmented_model(...)` returns **4** objects but the first downstream script unpacked **5**. The four-object interface was repaired at commit `8e711c9c4dd1`, and independently exact-source smoke-tested via [workflow #39237023131](https://github.com/zuizui0223/chocho/actions/runs/39237023131), without requerying the 239 species.
- A fresh *local* independent aggregation of the archived COMPLETE 135-taxon GBIF source, using the pinned author-only WCVP/HOSTS/BCE/CHELSA source/manifest and the original predeclared 72 species×period×effort×climate null contrasts, now yields the following numeric results. This is a reproducible local calculation on public GBIF sampled records, **not** an independently successful cloud analysis at this point. Do not describe a GitHub workflow as green unless its derived evidence artifact exists.

## Fixed principal comparison, completed 135 species only

The originally defined effort gate was **>=3 other independently sampled focal butterflies** in the same region and time period; native HOST-PLANT CHELSA BIO1/BIO12 **5–95% climate analogue**; comparison among regions outside the focal butterfly's native host range, conditional on butterfly species × WGSRPD Level1. This is not global survey absence.

| Metric | Real 135-complete-species sampled data |
|---|---:|
| Completed species with adequate original native-host climate | **131** |
| Species with comparable effort-supported resource presence/absence | **100** |
| Informative species × Level1 strata | **276** |
| Observed butterfly-positive, effort-supported candidate species×region cells | **119** |
| Positives in **added** larval-resource host regions | **74/119 = 62.18%** |
| Region-conditioned hypergeometric expected coincidence | **69.618/119 = 58.50%** |
| Excess observed compared with conditional expectation | **+4.382 cells; +3.682 percentage points** |
| Model-dependent 9,999-draw right-tail fraction, **NOT population significance** | **0.1127** |
| Genus cluster bootstrap 95% percentile interval (2,999 draws) | **−0.935 to +8.966 percentage points** |

**Decision:** Positive point estimate, **NOT statistically/biologically distinguishable from the sampled geographical null** under the predeclared primary comparison. Do not claim observed butterfly range expansion or introduced-host demographic benefit based on this incomplete-GBIF result.

## Predeclared sensitivity directions

| Same sampled source; other focal butterfly effort >=3 | Observed | Null expectation | Excess (pp) | Serial structural conditional-null tail |
|---|---:|---:|---:|---:|
| Broad q=.05 and Level1 only (primary) | 74/119 | 69.618/119 | **+3.682** | .113 |
| Broad q=.05 with additional Level1-local BIO1/BIO12 2×2 matching | 39/67 | 37.918/67 | **+1.615** | .389 |
| Narrow q=.10 and Level1 only | 45/78 | 43.251/78 | **+2.242** | .335 |
| Narrow q=.10 plus coarse climate-bin matching | 26/49 | 26.585/49 | **−1.195** | .656 |
| 2010–2017, broad q=.05, Level1 | 63/98 | 58.390/98 | +4.704 | .094 |
| 2018–2025, broad q=.05, Level1 | 61/96 | 56.014/96 | +5.193 | .069 |

Additional climate-bin conditioning and narrow climate-window choices weaken the positive point effect. These are **not** independently randomized causal tests, and the two historical record windows **cannot be interpreted as butterfly expansion over time** because host plant introduction time and sampling effort changed.

## Taxon contribution and negative-control diagnosis (same original 135 source records)

A separate offline aggregation of the same predeclared primary species × Level1 strata (no new GBIF retrieval) confirms all **74 observed vs 69.618 null-expected cells** and adds a taxonomic-concentration diagnostic: among **100 informative species**, only **26** have at least one effort-supported positive GBIF region in a primary stratum (15 positive, 9 negative, 2 effectively balanced/nonzero; the other 74 contribute no observed positives). All 100 remain in the candidate exposure denominator, including zero-positive strata. Leave **one entire butterfly genus** out at a time (73 genera) changes the pooled excess from +3.682 pp to a range of **+3.076 to +4.783 pp**; its sign remains positive. Thus a single genus does not explain the *point estimate*, yet the genus-cluster bootstrap interval **[−0.935,+8.966]** spans zero; neither genus leave-one-out nor a raw positive ratio establishes a biologically independent presence response.

Major mismatch is **ascertainment precision rather than one spectacular outlier**: under the first 135 complete taxa, only 119 GBIF-positive candidate resource regions remain after matching climate, native-resource unsupported geography and other-butterfly effort. In 96 other taxa, HTTP200/curl28/18-second hard-body timeout prevented clean completion; this is differential source missingness and cannot be treated as butterfly absence. The rescue pass preserving 300-ordinal sample windows through 75-record transport chunks is separately staged.

## Critical source selection and quality gates
- **96/239 butterflies failed complete source sampling**, mostly due to the original GBIF `limit=300` page responding HTTP200 but its payload exceeding an 18-second hard network read timeout (curl exit 28), which the code initially misclassified as a nonretryable error. Some windows were successfully sampled but were deliberately excluded when their same-species other period failed. This is missingness, **not true absence**, potentially correlated with GBIF volume/region.
- **Eight** GBIF butterfly taxon strings lacked a safe accepted exact species match under the project's strict rule; not silently converted to arbitrary taxa or treated as ecological absence.
- Systematic 300-record ordinal GBIF pages for 2010–2025 are nonrandom, nonexhaustive and may emphasize heavily observed regions. The internal `E_-i` background proxy of other butterflies is not equal to all biodiversity observer effort.
- Full 239-species generality, region-area correction, climate niche physiology, measured adult/egg/larval demographic outcomes and causal introduced host resources **remain unavailable**.

## Next exact technical repair to justify another macro pass
Refetch **ONLY the 96 incomplete species**, retaining the exact original selected two 300-record ordinal windows per time period but fetching each in four transport-only 75-record chunks. Retry curl timeout 28 despite HTTP200. Aggregate successfully re-fetched species exactly once and join with the 135 already completed species; do not repeat all original 239 GBIF source searches or silently increase the sample budget. Then rerun these exact 72 preregistered combinations from a versioned merged original GBIF aggregate with species completeness and query timestamps.

**Paper positioning:** The positive resource-climate structure remains well supported at the PLANT-RESOURCE layer, but actual butterfly observations do **not yet** give strong independent support beyond geography/climate/effort from the partial 135-species cohort. No individual published butterfly experiment paper reanalysis was needed. Draft [PR #46](https://github.com/zuizui0223/chocho/pull/46) remains separate from frozen original GEB PR #38 and main.
