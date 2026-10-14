# v16 — The early-flowering Thai Acacia provenance remains distinctive after clone-count matching

**2026-10-14. Post-hoc plant-only analysis. Not a measured butterfly cohort, proven causal genotype effect or new host-contributed butterfly range expansion.** Frozen GEB PR #38 unchanged.

## Question

In the Iriomote *Acacia auriculiformis* experiment of Kato et al. (2012), a Thailand source-g group of four clones/six ramets opens the stand flowering period before all other seven sources. Could this extension be a mechanical consequence of group size or an exceptional single clone? Earlier phenological differentiation of source g was **already reported by Kato et al.**, and is **not** claimed as a novel discovery.

## Exact observed union

| Three April–March seasons | Source g only 10-day bins | Others only | Shared | Total union |
| --- | ---: | ---: | ---: | ---: |
| 2008–09 | 9 | 2 | 8 | 19 |
| 2009–10 | 6 | 4 | 7 | 17 |
| 2010–11 | 8 | 3 | 10 | 21 |
| **Totals (season-bin units)** | **23** | **9** | **25** | **57** |

A total of 23 flowering **10-day categories across three seasons** is **not 230 measured flowering days**, bud biomass, insect habitat suitability or continuous larval food supply.

## Exact size-and-clone-count counterfactual

The 23 clones consist of seven one-ramet clones and sixteen two-ramet clones. The true Thai group has **two of each**. Thus we exhaustively assessed **choose(7,2) × choose(16,2) = 2,520** hypothetical *four-clone/six-ramet* groups, each scored as the number of 10-day flowering bins in which the selected group, but no clone outside it, was flowering.

Across the 2,520 possible labels, the **mean is 1.0103**, median **0**, nearest-quantile 95% interval **0–5** seasonal-bin units, second-highest **10**, maximum **23**. The actual Thai set (146/147/148/149) is the **only** group attaining 23. Its 1/2520 rank/frequency is a *descriptive combinatorial rarity*, not a meaningful genetic-treatment significance level: actual provenance groups were not randomly assigned, the interesting group was chosen after published inspection, plant phenology is correlated across years, and labeling randomly mixes true source categories.

### Critical flowering-consistency control

All six Thai-source plants bloomed in each of three years. If group labels are additionally restricted to **six plants that flowered in all three years**, only **15** alternative clone-group configurations exist. The Thai group is still the maximum, but **1/15** is a much weaker descriptive frequency. This condition *necessarily includes Thai single-ramet clones 146 and 149* because no other one-ramet clone flowered in every season. Consequently **flowering probability/persistence and early onset cannot be independently identified** with this dataset.

### Leave-one-clone-out robustness

Physically excluding each of the four source-g clones leaves **22, 23, 21, 23** group-exclusive 10-day bins across the three years, respectively. Each variant retains earlier flowering in all three seasons. Kato's 2009 narrative (maximum duration eight units) conflicts with original Table 7 records of clone 146 at 13 and clone 149 at 10. A sensitivity clipping those two to eight units changes neither 2009's six group-exclusive bins nor the three-year total 23.

## Why it matters for butterfly ecology, but remains a hypothesis

*Nacaduba biocellata* larvae in Singapore feed on *A. auriculiformis* flower buds and flowers for roughly 9.5–11 days in rearing (ButterflyCircle 2010). The Iriomote experiment samples **only plant flowering**, not the butterfly, flower bud biomass, or daily flower-bud replacement. Nonetheless it defines a testable path: botanical presence alone need not imply flowering-tissue availability during an entire larval-development window. Early and annually consistent flowering sources could affect temporal resource coverage of florivores, but **no such demographic effect is quantified here**.

Source: [Kato et al. 2012 original article](https://doi.org/10.1515/sg-2012-0029), Tables 4/5/7. Reproducible repository code: [matched_clone_null.py](../../research/nacaduba_stage_v16/matched_clone_null.py) pinned to the two original Git blobs in v15. [GitHub Action #38506434982](https://github.com/zuizui0223/chocho/actions/runs/38506434982) passed, independently reproducing all key totals and exporting JSON. Complete local ZIP `chocho_v16_balanced_clone_phenology_verified_20261014.zip` contains all 2,520 combinations, source-pinned CSVs, nine test cases, independent NumPy reconstruction, diagnostic plots and a SHA256 manifest (ZIP SHA256 `f41d6ef21de01aacd343c9ae4dc9b9786d3bc09c0388e00b42b99903db1e45de`).

**Boundary:** This result strengthens the **plant-side time-window premise**, not a demonstrated novel butterfly ecological effect. The fixed HOSTS/WCVP GEB paper and all source-conditioned estimates remain separate.
