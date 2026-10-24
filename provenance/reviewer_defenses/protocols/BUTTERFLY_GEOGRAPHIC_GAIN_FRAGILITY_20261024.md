# Fixed exploratory analysis: geographic gain versus host-support fragility
**24 October 2026.** Frozen before evaluating this specific paired association. This is a NEW exploratory macroecological hypothesis, not an independent confirmation or a preregistered effect. Preserve main GEB manuscript/PR38.

## Core question
Across butterflies that gain any potential resource geography from introduced host ranges, do large proportional gains come with **low spatial redundancy** of documented introduced hosts? The ecological hypothesis predicts a positive association between relative resource expansion and single-host dependence within newly added butterfly × region resource cells, conditional on recorded host list breadth and pre-existing native resource geography.

We distinguish:
- Extent: `x = log(contemporary_resource_units / native_resource_units)`. Original native geography is not butterfly-native range.
- Fragility **proxy**: `y = 1 − fraction_added_units_with_multiple_contributing_hosts`, the fraction of newly available regions for which exactly one *recorded* introduced host supports that butterfly. This does **not** establish that a butterfly depends on that host in the field.
- Host portfolio concentration as secondary: `maximum_single_host_fractional_share` of contributed added-cell credit across plant taxa. It is a separate derived catalogue metric, not directly interchangeable with site-level redundancy.
- Both x and y are deterministic functions of the same host-resource reconstruction. Statistical association is descriptive **and may reflect set-union/mathematical coupling**; no randomization or independent realized-use response supports a causal process.

## Frozen inputs, unit, exclusions
Exact three unchanged Git blobs (same as prior trait work) on `analysis/gbif-macro-recovery-20261020`:
1. `data/frozen/figure_sources/anthropogenic_species_metrics.csv` = `9f83d2d9ac8e46f1cfd4ded572c7a023273503dd`
2. `data/frozen/figure_sources/host_contribution_metrics.csv` = `4cf4851499ca2233ff2adadf1c4ee7f0e863f23d`
3. `data/frozen/figure_sources/s1_resource_descriptors.csv` = `09541b5894529cc8629c7515ea9022ded9edead2`.

Primary sample all 239 butterflies with exact-match joins and **positive `introduced_added_units`** (expected 206), positive native resource denominator and finite redundancy fraction. Zero-added taxa are excluded from the main regression, because their proportion of "added" units has no denominator—not because they are zero-redundant. No host link addition, time slicing, actual occupancy or quality filtering. Each butterfly species weighted once. Species' scientific genus (first word) defines the resampling unit.

## Primary, frozen regression
Species-level OLS `y ~ 1 + log_growth + log(native_resource_units) + log1p(resolved_host_species) + host_family_count`, treating n=1 species as unit and no cell-level pseudoreplication. Main directional association is the coefficient for `log_growth`, **positive hypothesized**. Report raw Spearman rank correlation (descriptive), adjusted coefficient, the pairwise Pearson residual correlation between x and y after controls, and genus-cluster 5,000-resample 95% percentile CI (fixed xorshift32 seed `20261024`). Do not interpret bootstrap positive-tail fraction as an independent p value.

Diagnostics are fixed *before* regression:
- repeat same control set on host-taxonomy adequate subset (`host_taxonomy_lower_bound_adequate=True`)
- repeat on butterfly species without any Poaceae host only if a frozen per-species source Poaceae flag is present; otherwise mark **not evaluable**, not inferred from family-level diet class
- leave-one-genus-out coefficient range
- secondary outcome `maximum_single_host_fractional_share` with same predictor/controls and same bootstrap.
- Show descriptive high/low expansion quantile groups for interpretation **without treating chosen strata as inferential replications**.

## Interpretation stop
Claim "gain–fragility trade-off" only as a within-model association if positive and genus-cluster CI excludes zero in both primary and adequate-host subset. If either interval includes zero or the sign flips, stop the positive assertion. A positive within-catalogue association, however strong, is **not a biological resilience or conservation prescription**; actual host-use, host biomass/phenology, butterfly demography, local source dependence, interaction-network homogenization and host-list incompleteness remain unverified.
