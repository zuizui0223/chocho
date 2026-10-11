# Post-hoc: taxonomic host breadth versus geographic redundancy of added butterfly resources

**Date:** 2026-10-11. **Status:** exploratory, not a preregistered result and not part of the currently frozen GEB submission package.

## Ecological question

Does being a taxonomic larval-host generalist protect a butterfly from relying on only one recorded host plant in newly available regions after human-mediated host redistribution?

**Counterfactual caution:** The species-level host family count is a global attribute; region-level redundancy is computed from *geographic union cells* within WGSRPD3. Neither represents an observed population-level realized diet, relative host biomass, larval performance or persistence.

## Frozen material

- Source head: `88209b8a215ac1bc9c4c5e67622fa104179c8895` (GEB PR #38 exact-head handoff, not modified).
- `data/frozen/figure_sources/anthropogenic_species_metrics.csv`
- `data/frozen/figure_sources/host_contribution_metrics.csv`
- Results receipt: `provenance/reviewer_defenses/results/butterfly_resource_redundancy_hostcount_20261011.json`
- Recompute: `python scripts/analyze_butterfly_resource_redundancy_hostcount.py --output-json /tmp/redundancy.json`.

The 206 butterflies with introduced-resource expansion reduce to **191** after applying the frozen `host_taxonomy_lower_bound_adequate=True` filter.

## Result A: the apparent taxonomic breadth effect is explained by recorded host species count

The response is the *fraction of added butterfly × region resource cells supported by more than one recorded introduced host species*.

1. Marginal Spearman correlation between host-family count and this resource-redundancy fraction: **rho = +0.3731** (191 species).
2. Average-tied ranks of both variables were separately demeaned within each **exact resolved-host-species-count** group; correlation of rank residuals: **r = -0.00041**.
3. The fixed-host-count conditional null permutes host-family ranks **within each exact host-count group** (9,999 permutations, deterministic seed 20261011). Two-sided **p = 0.9954**, null 2.5th–97.5th quantiles **[-0.1651, 0.1670]**.

Interpretation: A higher number of host families does not show an independently identifiable relationship with reconstructed *regional* resource redundancy after controlling the exact number of recorded host species. This is **not** proof of no effect, a causal demonstration, or a phylogenetically adjusted result. Host-family information is partly constrained by recorded host-species richness; only 158 species are in host-count strata with any family-breadth variation. A very similar adjustment already removed portfolio-*concentration* signals in the older `butterfly_architecture_hostcount_adjusted_v0.1.json`; the present test extends the check to the *fraction of newly added regions with a single recorded host*.

## Result B: high global host richness does not ensure redundant hosts in newly available regions

| Recorded host-species count | Butterfly species | Added butterfly×region units | Units with only one recorded host | Single-host share |
| --- | ---: | ---: | ---: | ---: |
| 1–2 | 32 | 1,874 | 1,653 | 88.2% |
| 3–4 | 37 | 2,035 | 1,426 | 70.1% |
| 5–9 | 42 | 3,396 | 1,831 | 53.9% |
| 10+ | 80 | 6,224 | 2,823 | **45.4%** |

For butterflies with 10+ known hosts, a butterfly-species bootstrap (10,000 replicates, resampling species and weighting by their added units) yielded a descriptive 95% interval of **40.8%–50.1%** for the single-host-unit share. It does not account for shared hosts, phylogenetic relatedness, or species-selection uncertainty.

Examples within the >=10-host subset with >=90% of *their added cells* supported by one recorded host include `Thorybes bathyllus` (12 recorded hosts, 56 added cells, 98.2% single-host), `Apodemia mormo` (17 hosts, 52 added cells, 100%), and `Graphium doson` (14 hosts, 5 host families, 30 added cells, 90%). These are **flags for data validation**, not claims about realized dependence.

## Relation to current literature

Gross, Kawahara & Daru (2026, *Nature Communications*, doi:10.1038/s41467-026-73236-4) explain global variation in larval host breadth using plant diversity, climate and island context. This exploratory analysis measures a *different variable*: geographic redundancy of fixed, already-recorded host links specifically in introduced-host-added WGSRPD3 regions. Earlier studies also distinguish continental and local host use (e.g. studies of local lycaenid specialization). The distinction itself is not unprecedented; the new quantitative result is the 191-butterfly global-scale redistribution reconstruction.

## Boundaries and next falsifiable test

The 45.4% is **not** the percentage of populations facing a one-host extinction bottleneck. WGSRPD3 is geographically coarse; WCVP range inclusion does not assure local occurrence, seasonally suitable tissue, oviposition, larval performance or biomass, and incompletely recorded hosts can bias apparent single-host support upward. The statistic is also not compared with a plant-level geographic null that controls each host's occupancy.

**Next test (not yet executed):** For each added butterfly × WGSRPD3 cell, randomize host identities with the same plant-family composition, native range size and contemporary introduced-area occupancy; compare the observed single-host share to expected overlap. Follow with empirical butterfly occurrences and phenology/performance checks in candidate species/regions. A positive geography-adjusted excess of single-host cells would support unusual spatial segregation of known hosts; a null result would mean most of the pattern is predicted by the geometry and completeness of host records.

The ecological claim suitable *now*: **Global diet breadth and geographic backup options should not be treated as equivalent.** Do not promote a conservation vulnerability or demographic mechanism claim from this reconstruction alone.

## Branch policy

These files are isolated on `analysis/resource-redundancy-20261011` and should **not** silently alter the publication-ready exact-head GEB PR #38 or its anonymous bundles. Publication requires deliberate author review and a clear post-hoc label.
