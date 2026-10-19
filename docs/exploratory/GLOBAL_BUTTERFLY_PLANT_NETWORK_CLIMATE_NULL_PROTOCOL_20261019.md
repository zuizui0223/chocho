# Host-identity-preserving macroecological climate-analogue null

**Protocol:** exploratory posthoc added 19 October 2026, after positive geographic and butterfly×region margin tests were observed. **NOT pre-registered, NOT a new species observation, and NOT a randomized historical introduction experiment.** Fixed before the formal independently-seeded 2-chain source-verified null run. This uses only the same pinned original WCVP/HOSTS/CHELSA databases and their derived species resources; no more individual literature reanalysis.

## Central falsification
The previous climate-analogue excess (+8.28 pp beyond a uniform within-continent reference; ~+3.74 pp beyond a butterfly×region fixed-margin null) might be completely explained by **which plant families naturalize where**, rather than unusually well-matched consumer host opportunities. To assess this, hold the entire *butterfly–host network* fixed and permute plant introduction cells only among source-identical **plant-family × WGSRPD Level1 × CHELSA BIO1/BIO12 2×2 stratum** groups.

### Exact model inputs / values that must fail closed
- 239 butterfly species, HOSTS-recorded accepted host plants, 365 region units, 26,530 native potential butterfly-region cells, 14,553 introduced-added cells; source archive author-only, manifest checksum verified.
- 30,373 introduced **host plant × region** edges for the relevant source graph, irrespective of whether those host plants are used by one or many butterflies. More than 85% of these edges must belong to a swappable family–continent–climate stratum.
- 234 species with ≥10 native-host regions covering CHELSA BIO1/BIO12; these species originally contribute 14,481 added butterfly×region cells, including 9,451 in BIO1/BIO12 native-host 5–95% analogues, and 6,948 in 10–90% analogues.
- All climate data are 1981–2010 region summaries, and all sources are source-pinned. Native plant regions, insect host identities and the list of focal species NEVER change during randomization.

### Swap invariants
Within each plant-family × Level1 × discretized climate group, swap (plant A→region 1, plant B→region 2) into (A→2, B→1) only if:
1. Both plant IDs and both region IDs differ;
2. Neither added plant region is already a native or introduced cell of its recipient plant;
3. No duplicate host plant × region edge is created.

This preserves the **number of introduced regions of every individual plant in each broad climate bin and continent**, every regional count of introduced plants *within plant family and climate bin*, and all original HOSTS butterfly-to-plant edges. It does **not** preserve fine-grained climate placement within each bin, regional spatial autocorrelation, a given butterfly's count of new cells, or exact local biomass.

### Response and inference gate
Because randomizing plant geography changes overlapping host unions, the total newly supported butterfly cells A can change. Therefore **DO NOT compare raw analogue counts O as if A were fixed**. Compare **O/A**, the proportion of newly added butterfly×region resource opportunities climate-analogous to that butterfly's fixed *native-host* climate envelope. Compute O and A for all 234 evaluable butterflies jointly for each 5–95 and 10–90 threshold. Record A each draw, not just the ratio; no fixed-denominator claims.

Two independent MCMC chains with seeds 20261019/20261020, burn 10× number of host-plant introduction edges, spacing 4× edge count, 199 observations per chain. Run only on the frozen source graph; no ecological outcome resampling or source reclassification. Require:
- ≥85% source introduced plant edges climate-covered, ≥40% in viable swaps;
- Each chain ≥10% of originally introduced host plant edges moved after burn, ≥100 accepted burn swaps, ≥100 accepted per sample-spacing on average;
- Every chain >30 positive-sequence approximate effective samples per outcome (O/A), between-chain medians within 2× larger chain standard deviation;
- Source counts, total introduced host-incidence budget and basic host→plant group sizes auditable. A nonconverged chain cannot support a formal claim. If gates fail, report NOT EVALUABLE; do not choose shorter climate bins or alter the ecological response after seeing signs.

The 199 serial draws are NOT 199 independent permutations; add-one serial-tail frequency is *not* an exact randomization p. Compare observed-minus-median in percentage points and present both two-chain null medians, ESS and qualifications; interpret small residuals in terms of effect size, not nominal tail alone.

## Additional ecological nulls not addressed here
The model preserves known host-plant identities and coarse location/climate family margins, but not native host **plant** climatic tolerance within introduced climate-bin neighborhoods, actual butterfly larval performance, phenology, real butterfly geographic occurrences, historical introduction timing, native butterfly dispersal, or sampling effort. It could leave a biologically small geometric residual even without any selective increase in butterfly populations. Results from this null are conditional on HOSTS completeness; the Europe-focused BCE augmented candidate host catalogue has not been rerun for this exact test.

## Decision gates
- If plant-network-constrained null explains the original resource alignment, demote the "butterfly-specific" climatic association to plant/region geometry: no novel consumer-residual pattern.
- If positive residual persists at both predeclared thresholds and mixing gates pass, describe a **small host-network-conditional macroecological residual** in potential climate-analog resource geography, not new species' reproduction.
- Keep exploratory branch/PR Draft; do not modify main or the frozen original GEB v0.2 submission bundle. Source rights for Kew/TDWG-derived data remain subject to manual review; upload only derived numeric JSON.
