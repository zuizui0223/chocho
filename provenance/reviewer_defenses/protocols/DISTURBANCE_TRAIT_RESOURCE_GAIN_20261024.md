# Fixed analysis plan: disturbance affinity and geographic butterfly-resource gain (2026-10-24)

**Status:** new exploratory hypothesis, fixed before evaluating its regression coefficients or group outcome summaries. **Scope:** potential larval resources, not adult abundance, realized host use or breeding. This branch is separate from the frozen GEB manuscript.

## Question and directional hypothesis
Do butterflies annotated as strongly disturbance-associated receive greater proportional potential larval-resource gains from documented host-plant redistribution than strongly disturbance-avoidant butterflies, after adjusting for original native host-resource geography and host catalogue breadth?

A positive result would indicate ecological trait sorting of reconstructed *opportunity*, **not** higher fitness, realized range expansion, competitive exclusion or causal facilitation of disturbance-associated butterflies.

## Source-locked design
Join exactly by full accepted butterfly species name:
- `data/frozen/figure_sources/anthropogenic_species_metrics.csv`
- `data/frozen/figure_sources/s1_resource_descriptors.csv`
- `data/frozen/figure_sources/host_contribution_metrics.csv`

Keep the original 239-species resource panel intact for QC. Primary inferential subset: exact descriptors `Disturbance-associated (strong)` versus `Disturbance-avoidant (strong)` only. The weak, mixed and variable labels are **excluded from the two-group inferential comparison**, never recoded as avoidant. The source species name first token defines the genus cluster. Do not silently introduce another host source or a new geography.

Primary response: `log(contemporary_resource_units/native_resource_units)` for positive native resource units; direction associated minus avoidant. Count regions once per butterfly, irrespective of how many host plants support a region. Secondary response (separately labelled): fraction of introduced-added units that have two or more documented contributing hosts among butterflies with strictly positive added units. It measures **catalogue-recorded support redundancy**, not demographic resilience.

## Primary estimation
1. Report complete-data counts, exact two groups, unadjusted per-group mean/median, raw difference and number with zero gain.
2. Fit equal-species-weight linear regression (OLS), `y ~ disturbance_strong + log(native_resource_units) + log1p(resolved_host_species) + host_family_count + wing_size_proxy + I(voltinism=B) + I(voltinism=U)`. M is the other explicit category. Drop only rows missing a required model covariate; report missingness. No added-resource predictor appears among controls.
3. Estimate 95% percentile intervals by **5000 genus-cluster bootstraps**, with deterministic seed `20261024`, sampling genera with replacement and retaining their exact species groups. Discard degenerate resamples with both groups absent or a singular design and report the count. Report the observed coefficient, CI and fraction of bootstrap coefficients greater than zero, but **do not interpret it as a confirmatory p-value**.
4. Two predeclared sensitivity estimators: reduced control set (same but no wing size or voltinism), and full-control model restricted to `host_taxonomy_lower_bound_adequate=True`. Report whether their signs and cluster-bootstrap uncertainty agree; if a subset is underpowered, say so.
5. Compute leave-one-genus-out ranges for the full model, and descriptive mean expansion in the two groups. A multi-family phylogenetic effect cannot be inferred with genus bootstrap.

## Secondary, not eligible to rescue a null primary
Repeat a fully labelled descriptive adjusted analysis for `fraction_added_units_with_multiple_contributing_hosts` using the same adjustment set **plus** `log1p(introduced_added_units)`, restricting to positive-added taxa. Genus-cluster bootstrap CI as above. This is exploratory dependence on plant records, not demonstrated local redundancy.

## Stop rules and inferential status
- Never mine weaker/mixed disturbance classes, edge/canopy/moisture categories, arbitrary thresholds or time periods to reverse an unconvincing primary sign.
- A positive point estimate without a sufficiently precise genus-cluster interval is `suggestive/NOT ESTABLISHED`.
- Association may reflect life-history trait coding, broad environment affinity, geographic opportunity, plant source sampling and phylogeny. Source data do not identify settlement, larval performance or actual competition.
- Any result here is a **new exploratory within-catalogue analysis**. Preserve main manuscript and prior GEB claims. If source joins or arithmetic QC fail, return `NOT EVALUABLE`.
