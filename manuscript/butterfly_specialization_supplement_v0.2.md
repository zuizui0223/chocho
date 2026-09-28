# Supplementary Information — butterfly resource geography v0.2

**Associated manuscript:** *Anthropogenic expansion of butterfly resource geography is concentrated in network-prominent host plants*

This file collects post-hoc robustness analyses that define the manuscript's inference boundaries. The main text reports the ecological results; this supplement provides the parameter and sensitivity details.

## Supplementary Table S1. Strict host-bias null parameter sensitivity

The strict null preserves butterfly host-species richness and plant-family composition, matches candidate hosts on native WGSRPD3 breadth, and weights candidates by the number of other Lepidoptera species recorded using each plant in HOSTS.

Observed values in the 207-species strict-null panel:
- total introduced-added units = **12,853**
- mean log expansion = **0.4390**
- median log expansion = **0.3330**
- mean other-Lepidoptera consumers per observed host = **25.28**

| Native-range bandwidth | Degree exponent | Mean sampled degree | Null median total added | p(total ≥ observed) | p(mean ≥ observed) | p(median ≥ observed) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 0.5 | 11.52 | 9,809 | 0.0033 | 0.0033 | 0.0033 |
| 0.1 | 1.0 | 19.37 | 12,783 | 0.433 | 0.023 | 0.023 |
| 0.1 | 2.0 | 33.00 | 16,911 | 1.000 | 1.000 | 1.000 |
| 0.2 | 0.5 | 11.69 | 9,765 | 0.0033 | 0.0033 | 0.0033 |
| 0.2 | 1.0 | 21.10 | 13,038 | 0.653 | 0.087 | 0.083 |
| 0.2 | 2.0 | 37.59 | 17,631 | 1.000 | 1.000 | 1.000 |
| 0.4 | 0.5 | 11.72 | 9,777 | 0.0033 | 0.0033 | 0.0033 |
| 0.4 | 1.0 | 22.41 | 13,349 | 0.810 | 0.373 | 0.270 |
| 0.4 | 2.0 | 42.08 | 18,153 | 1.000 | 1.000 | 1.000 |

Each cell used 299 deterministic randomizations. The sensitivity is strongly governed by how closely the null reproduces observed host-network prominence: exponent 0.5 under-matches prominence and retains an apparent host-identity excess, whereas exponent 2 over-matches prominence and generates null portfolios more expansion-prone than observed. A separate prominence-calibrated null therefore selects the degree-weight exponent using only the observed host-degree target, not expansion outcomes.

## Supplementary Table S2. Plant network prominence and anthropogenic geographic expansion

Across **8,909 host plants in 278 plant families**, plant network degree is the number of distinct Lepidoptera species in the fixed HOSTS-WCVP reconstruction recorded using the plant.

| Lepidoptera consumer degree | Plants | Fraction with introduced-range expansion | Median added WGSRPD3 units | Median log expansion |
|---|---:|---:|---:|---:|
| 1 consumer | 4,060 | 0.262 | 0 | 0 |
| 2 consumers | 1,565 | 0.357 | 0 | 0 |
| 3–5 consumers | 1,647 | 0.430 | 0 | 0 |
| 6+ consumers | 1,637 | 0.659 | 4 | 0.176 |

Association diagnostics:
- Spearman log degree vs log expansion: **rho = 0.313**
- partial rank association controlling log native breadth: **0.267**
- family-adjusted rank correlation: **0.341**
- partial rank association with absolute added units controlling native breadth: **0.272**
- within plant-family × native-breadth-quintile strata: **r = 0.295**
- stratified permutation: **0/4,999** correlations as extreme; Monte Carlo **p < 0.001**
- stratified null 95% interval: **-0.0222 to 0.0221**

HOSTS consumer degree is interpreted as a joint axis of ecological host prominence/commonness and study/recording intensity, not as a pure biological trait.

## Supplementary Table S3. Species-level occurrence robustness

The occurrence analysis reuses the 32-species panel originally stratified for the climate analysis; it was not prospectively designed as a validation panel. Twenty-three species had at least one contemporary butterfly occurrence outside the native host-resource envelope.

Overall:
- recovered units = **66/115 (57.4%)**
- mean species recovery fraction = **0.583**
- median species recovery fraction = **0.600**
- region-matched null median mean species recovery = **0.385**
- null 95% interval = **0.308–0.470**
- **0/99,999** null draws reached the observed mean species recovery
- species-cluster bootstrap pooled recovery 95% interval = **0.349–0.774**

High-leverage exclusions:
- excluding *Pyrgus communis*: **44/93 = 47.3%**, region-matched null median 30 recovered units (24–36), **0/99,999** exceedances
- excluding *Pieris brassicae*: **62/111 = 55.9%**, null median 46 (40–52)

Examples of incomplete or failed recovery:
- *Erynnis tristis*: **0/11**
- *Historis acheronta*: **5/13**

## Supplementary Table S4. Taxonomic and phylogenetic non-independence

Butterfly-Family random intercept, n = 239:
- standardized rank host-breadth coefficient = **0.021**
- 95% CI = **-0.111 to 0.154**
- Wald p = **0.752**
- Family random-intercept variance estimated at approximately zero

Kawahara et al. (2023) species-level phylogenetic GLS:
- exact tree matches = **124/239 (51.9%)**
- Brownian rank-PGLS beta = **0.140**, 95% CI **-0.036 to 0.317**, p = **0.119**
- Brownian raw-log PGLS beta = **-0.0024**, p = **0.976**
- estimated Pagel lambda = **0.033**
- Pagel rank-PGLS beta = **-0.051**, 95% CI **-0.234 to 0.131**, p = **0.581**

Unmatched panel species were excluded rather than phylogenetically imputed.

## Supplementary Table S5. Climate-distance sensitivity and effect-size precision

Twenty-four species were climate-informative.

| Comparison | Median filtering score | Species > 0.5 |
|---|---:|---:|
| Original cross-fit | 0.801 | 23/24 |
| Same WGSRPD level-1 region | 0.846 | 23/24 |
| Distance matched, 250 km | 0.714 | 21/24 |
| Distance matched, 500 km | 0.750 | 23/24 |
| Distance matched, 1,000 km | 0.750 | 22/24 |

Primary host-breadth release test:
- partial Spearman rho = **-0.166**
- one-sided permutation p = **0.2237**
- 30,000-species-bootstrap 95% interval = **-0.583 to 0.261**
- approximate magnitude required for 80% power: **|partial rho| ≈ 0.505**

This analysis is secondary in v0.2 and is shown as Supplementary Figure S1.

## Supplementary inference boundaries

1. Resource envelopes are reconstructed geographic opportunities, not realized local host use.
2. WCVP introduced status does not provide dates or causal introduction pathways.
3. Network degree cannot separate ecological host commonness from HOSTS recording intensity.
4. Occurrence recovery does not establish local larval use or host-caused colonization.
5. The occurrence panel was climate-stratified rather than designed as a validation sample.
6. Null-model, ceiling, phylogenetic and geographic sensitivity analyses were added post hoc after manuscript review.
