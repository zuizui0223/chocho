# Plant-network preserving null for butterfly climate-analogue opportunity
**Fixed design 2026-10-19, prior to qualified long-chain results.**
**Status: EXPLORATORY, original 239-species HOSTS/WCVP/CHELSA macro analysis, NOT review of individual publications.**

## Core scientific question
The source-verified 239-butterfly climate-analogue result found 9,451/14,481 added butterfly×region units inside native-host climatic analogues (65.265%), versus a within-continent uniformly allocated opportunity benchmark (56.984%). A stricter BUTTERFLY×REGION fixed-margin null preserved each butterfly's number of added regions per continent and each region's number of newly supported butterflies, giving medians 8,908–8,910 (additional +~3.74 pp). Its two 199-draw chains passed preset ESS/mixing checks in [run 39211273830](https://github.com/zuizui0223/chocho/actions/runs/39211273830).

**Critical unresolved inference:** A butterfly-resource climate alignment residual could be entirely inherited from the ecological climate niches of the already-naturalized HOST PLANTS, without any additional trophic-network-specific matching. A butterfly×region permutation does NOT preserve plant×region geography or host identities.

Here test whether the observed climate-aligned *fraction of added resource opportunity* differs from independent null distributions that preserve (i) fixed butterfly→accepted-host associations; (ii) every host's natural distribution; (iii) introduced plant×WGSRPD3 incidence degrees within geographical/climate strata; and (iv) region×plant-family introduced-host incidence degrees. These are plant-side ecological opportunity nulls, **not colonisation/fitness hypotheses**.

## Same pinned source and denominators
- Retrieve and verify the author-only 3.25MB source ZIP from run 39006554609 (SOURCE_MANIFEST + SHA256 member checks) through GitHub Actions. Do not redistribute raw WCVP/HOSTS derived data without rights review; upload only derived statistical aggregate.
- Rebuild the same 239×365 accepted-host native/introduced resource reconstruction, same CHELSA BIO1/BIO12 region vectors, and species-specific climate analogue masks defined using native host geography and 5–95%, with predeclared 10–90% sensitivity.
- Original all-species gains 14,553. Native climatic support >=10 eligible for 234 species and 14,481 added opportunity cells; observed aligned cells 9,451 (5–95%) / 6,948 (10–90%). Gate these exact totals; never select only butterflies displaying positive alignment.
- The plant swaps change the number of distinct butterfly×region union cells, because host occurrences overlap differently, so **primary response is the ratio `aligned_added/evaluable_added`**, not raw aligned count. Always report null denominator change (this is essential and cannot be hand-waved as plant supply fixed at the butterfly level).

## Prespecified complementary nulls (same observed response)
1. **PLANT_FAMILY×CONTINENT** — for each accepted plant, retain number of introduced WGSRPD3 cells per Level1, keep plant family and each region's count of introduced plants per family. Swap accepted-host-plant × region edges within plant family × Level1, rejecting placements in a host's own native range and duplicate introduced edges. All original butterfly–host associations remain unchanged. Uses original `make_groups` and `attempt_range_swaps` with the same native exclusions.
2. **PLANT_FAMILY×CONTINENT×CLIMATE_2x2** — same controls plus the host plant's continent-local BIO1 and BIO12 binned climate occupancy, using 2×2 quantile classes defined on all climate-covered regions within each Level1. Swaps only within original plant family, Level1 and same climate class, preserving introduction counts by coarse climate stratum. Use pre-existing `climate_strata` and `grouped_introduced_edges`. This conditions explicitly on broad plant climatic sorting (which the first null does not).
These plant nulls do **not** preserve each butterfly's gained-region count, or each region's butterfly opportunity degree; the prior fixed butterfly margin test targets a different null.

Each null: two fresh chains with independent seeds 20261019/20261020, 10×original introduced plant-edge burn attempts, 199 recorded samples, 4×original edge-count attempted swaps between draws. No reuse of chain-mutated edges between chains. Random generator fixed; algorithm deterministic for pinned source. Record both thresholds from the **same graph draw**. Full serialized draws must be retained in derived receipt (not just transcribed means).

## Required gates and interpretation hierarchy
- Frozen source SHA, exact observed original counts and original plant family/network relationship counts PASS.
- Every accepted graph swap preserves a plant's introduction cells per specified family/Level1(/climate stratum) and native plant exclusion. Verify restored counts per plant×stratum and region×plant family at burn/end and verify baseline reconstructed butterfly counts at start.
- MCMC numeric checks per threshold and mode: both chains approx ESS>=30; chain null medians differ by no more than 2×maximum within-chain null SD; burn moves >=10% original introduction edges, and >=50 accepted swaps per sampling spacing. If ANY check fails, status `NOT_EVALUABLE`, no null-tail claims.
- Show observed difference in *percentage points of aligned added opportunity*, 2.5/50/97.5% null percentiles, effect in pp, and separate chain-tail diagnostics. Serial MCMC nominal tails are NOT exact independent permutation p-values.
- If plant family×continent null explains bulk of butterfly×region margin excess, demote interpretation to **plant-side biogeographic sorting**.
- If additional climate_2x2 condition removes the remaining excess, interpret as plant climatic conservation, not a novel butterfly-specific interaction. If it survives both under good mixing, a residual is compatible with nonrandom **fixed host association × introduced host geography** alignment, but still is not a causal interaction, demographic benefit or current butterfly climate niche.
- Report sensitivity to 5–95% vs 10–90% and the changing null denominator. Do not select a positive version or remove apparently uninteresting climate restrictions post-hoc.
- No direct climate physiology measured, no environmental **time series**, no larval recruitment, no observed butterfly population trends, no local management prescription. Do not merge results into frozen GEB PR #38. Stop if comprehensive plant-network preservation eliminates the signal instead of inventing a new hypothesis after the fact.

## Main biological decision to enable
Compare whether climate-similar *potential larval-resource availability* is an emergent property of **plant distributions and their introduction filters** (host identity-free geographic mechanism) or a genuinely **interaction-network-specific sorting residual** after fixing plant family, continental movement, broad climatic propensity, and all original butterfly–host edges.
