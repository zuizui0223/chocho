# New ecological finding: introduced larval hosts reinforce pre-existing resource geography

**Date:** 2026-10-11  
**Status:** post-hoc, model-based spatial reconstruction. Do not silently add to the frozen GEB submission.  
**Source paper head:** `88209b8a215ac1bc9c4c5e67622fa104179c8895`.  
**Result:** `provenance/reviewer_defenses/results/butterfly_resource_redundancy_geonull_decomposition_20261011.json`.  
**Successful exact workflow:** https://github.com/zuizui0223/chocho/actions/runs/38112574884 (199 permutations, seed 20261011).

## Main result

Plant introduction expands known larval-host resource geography in the frozen reconstruction (206/239 butterflies gain resource regions; aggregate +54.9%). Nevertheless, the **observed configuration** produces significantly *less* novel butterfly × WGSRPD3 resource geography than a constrained reshuffling of the same introduced plant–region incidences.

| Panel | Observed novel resource cells | Null median | Difference | Relative difference | Null one-sided Monte Carlo p |
| --- | ---: | ---: | ---: | ---: | ---: |
| All 239 butterflies | 14,553 | 15,490 | -937 | -6.05% | 0.005 |
| 191 quality-filtered, expanded butterflies | 13,529 | 14,445 | -916 | -6.34% | 0.005 |
| 80 butterflies with 10+ resolved hosts | 6,224 | 6,681 | -457 | -6.84% | 0.005 |

Fixed by the null: the set of butterfly-host links and plant taxonomic families; each plant's introduced WGSRPD3 count within each WGSRPD level1 region; the per-WGSRPD3 introduction count of plants in each botanical family; and all native host ranges as forbidden introduction cells. The observed native-only butterfly envelopes are unchanged.

## Exactly why the observed novel geography is smaller

Across the 239-butterfly panel, define each occurrence of an introduced *known* host for a given focal butterfly and region as one introduced host-region incidence. The 75,346 incidences split exactly into:

1. **Already covered by another native known host:** observed **46,387**, null median **45,097** (observed +1,290; upper-tail p=0.005).
2. **Extra introduced-host overlaps outside the native resource envelope:** observed **14,406**, null median **14,757** (observed -351; lower-tail p=0.005).
3. **Unique newly resource-covered butterfly × region units:** observed **14,553**, null median **15,490** (observed -937; lower-tail p=0.005).

The three observed counts sum exactly to 75,346. Because individual null medians are taken separately, medians need not sum identically.

### Interpretation

**The shortfall of novel resource geography comes from an excess of introduced-host records in regions where the butterfly already has another native known host.** It is *not* evidence that introduced hosts are unusually piled on top of one another within previously unsupported regions: the latter component is below null expectation.

This is a previously unresolved dimension of the original fixed-host redistribution counterfactual: introduced plant ranges can **reinforce known host availability at already resource-covered places** rather than creating proportional geographic expansion opportunities. The data do not distinguish introduction routes, climatic filtering, spatially autocorrelated recording effort or human activity as causes.

The former single-host finding also needs explicit downgrading: the pooled share of newly added resource cells with only one known introduced host is 58.9%, close to the null median of 59.2%; among butterflies with 10+ known hosts it is 45.36%, close to the null median of 45.38%. This is **not an extraordinary spatial deficit of resource redundancy**.

## Execution and checks

- Exact original frozen host-native-contemporary reconstruction reproduced: 26,530 native butterfly × region units, 14,553 newly added units, and 8,574 single-host-supported added units.
- 74,055 unique introduced plant × region edges across the fixed HOSTS-WCVP pool; 72,065 edges in blocks eligible for some swaps.
- 7,536,378 accepted double-edge swaps / 44,103,780 attempts (17.1% acceptance), and all row/column/native-exclusion invariants passed.
- Geographic-null simulation: 199 samples from **one constrained swap chain**. Monte Carlo p values have resolution 0.005 and are exploratory; additional independent seeds, convergence diagnostics, geographic/climatic constraints and source-completeness sensitivity are not yet assessed.
- The strict post-hoc nature prevents interpreting the result as an a priori confirmatory test or as demographic or adaptive evidence.

## Ecological next tests

1. **Climate-conditioned swap null.** Within each botanical family and continent, additionally constrain interchange to comparable temperature and moisture bins to determine whether abiotic niche sorting generates the above-native overlap.
2. **Anthropogenic exposure / land use.** Test whether human footprint and cultivation/trade history explain excess native-envelope reinforcement after climate matching.
3. **Observed butterfly range response.** Compare time-resolved occurrence gain in truly new resource cells to already resource-covered cells (with spatial observation-effort adjustment). The current static reconstruction cannot identify historical arrival order or colonization causality.
4. **Source completeness.** Test whether additional BCE/Clarke host links reproduce this overlap decomposition in the geographically biased European sensitivity panel.

**Publication position:** continue the existing GEB submission at its tested, frozen exact head. Treat the native-resource saturation result as a promising but post-hoc second-stage ecological analysis until the climate/recording sensitivity is independently tested.
