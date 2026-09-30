# Supplementary Information — butterfly resource geography v0.2

**Associated manuscript:** *Human redistribution of host plants expands butterfly resource geography across the specialization spectrum*

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

Each cell used 299 deterministic randomizations. The sensitivity is strongly governed by how closely the null reproduces observed host-network prominence: exponent 0.5 under-matches prominence and retains an apparent host-identity excess, whereas exponent 2 over-matches prominence and generates null portfolios more expansion-prone than observed.

## Supplementary Table S2. Exploratory plant network prominence and anthropogenic geographic expansion

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

## Supplementary Table S4. Diet-breadth precision, taxonomic and phylogenetic non-independence

Full-panel host-family breadth versus log resource expansion, n = 239:
- Spearman **rho = 0.008**
- 49,999-bootstrap 90% interval = **-0.092 to 0.108**
- 49,999-bootstrap 95% interval = **-0.111 to 0.128**
- post-hoc equivalence margin = **|rho| < 0.10**
- approximate Fisher-z TOST = **p = 0.078**

The equivalence diagnostic is deliberately strict and post hoc. It does not establish that the effect is exactly zero; instead, it shows why the title and conclusions use "largely independently" and "little relationship" rather than an exact-null claim.

Butterfly-Family random intercept, n = 239:
- standardized rank host-breadth coefficient = **0.021**
- 95% CI = **-0.111 to 0.154**
- Wald p = **0.752**
- Family random-intercept variance estimated at approximately zero

Kawahara et al. (2023) species-level phylogenetic GLS:
- exact tree matches = **124/239 (51.9%)**
- estimated Pagel lambda = **0.033**
- estimated-lambda rank-PGLS beta = **-0.051**, 95% CI **-0.234 to 0.131**, p = **0.581**
- fixed Brownian lambda = 1 rank-PGLS beta = **0.140**, 95% CI **-0.036 to 0.317**, p = **0.119**
- fixed Brownian raw-log PGLS beta = **-0.0024**, p = **0.976**

Because the fitted Pagel lambda is close to zero, the estimated-lambda model is treated as the data-adaptive phylogenetic correction; the fixed Brownian model is retained as a boundary sensitivity and demonstrates that the selective 124-species subset cannot exclude a moderate positive rank effect under lambda = 1.

Unmatched panel species were excluded rather than phylogenetically imputed.

### Poaceae host-guild sensitivity

Because seven of the ten largest plant contributors were Poaceae, we tested whether grass-feeding butterflies—especially one-family Poaceae specialists—were creating the near-zero diet-breadth slope. **Fifty-eight** of the 239 butterflies had at least one resolved Poaceae host. Removing all 58 left **rho = 0.026** (n = 181; 49,999-bootstrap 95% interval **-0.116 to 0.168**). Removing only the **34 one-family Poaceae specialists** left **rho = 0.013** (n = 205; 95% interval **-0.122 to 0.147**).

Poaceae users did contribute more absolute added geography (median **69.5** added WGSRPD3 units versus **54.0** among non-Poaceae users), but their median log proportional expansion was slightly lower (**0.310** versus **0.329**). Thus the redistributed-grass guild contributes substantially to absolute resource geography without explaining the near-zero relationship between family-level diet breadth and proportional gain.

### Butterfly-family distribution of resource expansion

Resource expansion was not confined to one major butterfly lineage. Among families represented by at least 10 focal species, the fraction gaining at least one introduced-host resource region ranged from **83.5% to 90.2%**.

| Butterfly family | Species | Expanded | Fraction expanded | Median host-family breadth | Median added WGSRPD3 units | Median log expansion |
|---|---:|---:|---:|---:|---:|---:|
| Nymphalidae | 103 | 86 | **83.5%** | 2 | 51 | 0.285 |
| Hesperiidae | 51 | 46 | **90.2%** | **1** | 56 | 0.359 |
| Pieridae | 43 | 37 | **86.0%** | **1** | 78 | 0.405 |
| Lycaenidae | 25 | 22 | **88.0%** | 3 | 73 | 0.367 |
| Papilionidae | 16 | 14 | **87.5%** | 2 | 49.5 | 0.388 |

Riodinidae was represented by a single species and is not interpreted as a family-level estimate. These summaries are descriptive; unequal family sample sizes preclude treating the table as a comparative test among butterfly families.

Tree coverage was not fully representative of expansion magnitude. Exact tree matches and unmatched species had the same median host-family breadth (2 vs 2; Wilcoxon p = 0.677), but matched species had greater median log resource expansion (0.376 vs 0.285; p = 0.0012) and a higher fraction with any expansion (0.952 vs 0.765; Fisher p < 0.001). Importantly for the focal slope, host-family breadth versus log expansion was near zero in both groups (rho = -0.055 among matched species and +0.053 among unmatched species). We therefore use PGLS only as a sensitivity for the breadth–expansion association, not as an estimator of panel-wide expansion magnitude.

## Supplementary Methods S1. Climate filtering within contemporary resource opportunity

This secondary analysis tested the pre-specified prediction that broader host-family diets weaken climatic filtering within contemporary host-resource opportunity. Species entered only when occurrence coverage, host-taxonomy consistency and within-envelope sampling met fixed quality criteria; **24 species** qualified.

Climate was represented by CHELSA BIO1, BIO7, BIO12 and BIO15 (Karger et al. 2017, 2021). Within each species, observed resource units were split deterministically into training and evaluation sets. Climate centre and scale were estimated from training occurrences. Held-out observed resource units were then compared with effort-supported contemporary host units in which the butterfly was not observed.

The climate-filtering score was the probability that a never-observed resource unit had greater climatic mismatch from the training niche than a held-out observed unit, with ties receiving half weight. A score of 0.5 is neutral. The primary host-breadth test used partial Spearman correlation between host-family breadth and filtering score while controlling contemporary resource breadth, with 9,999 residual permutations and a one-sided negative alternative.

To assess geographic confounding, sensitivities restricted comparisons to occupied WGSRPD level-1 regions and matched held-out observed and never-observed resource units by distance to the nearest training-observed region using 250-, 500- and 1,000-km calipers. A species bootstrap summarized effect-size precision. Because this analysis concerns realization of reconstructed opportunity rather than the primary resource-redistribution question, it is reported entirely in Supplementary Information.

## Supplementary Table S5. Climate-distance sensitivity and effect-size precision

Twenty-four species were climate-informative.

| Comparison | Median filtering score | Species > 0.5 |
|---|---:|---:|
| Primary cross-fit | 0.801 | 23/24 |
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

## Supplementary Table S6. Concentration of added opportunity across host plants and structural portfolio diagnostics

### Across-plant contribution concentration

For every butterfly × WGSRPD3 unit that entered the contemporary resource envelope only after introduced host distributions were included, one unit of credit was divided equally among all host species that contributed that unit. This preserves the aggregate total exactly: plant credits sum to **14,553 added butterfly × WGSRPD3 units**.

Across **670 contributing host species**:

| Ranked contributors | Cumulative share of all added units |
|---:|---:|
| Top 1 | **4.55%** |
| Top 5 | **15.70%** |
| Top 10 | **25.10%** |
| Top 20 | **37.27%** |
| Top 50 | **57.00%** |
| Top 100 | **73.64%** |

The smallest prefix reaching 50% contained **38 host species**. The plant-level Gini coefficient was **0.777**, HHI was **0.0108**, and the inverse-Herfindahl effective contributor number was **92.5 species**. Thus the distribution is strongly right-skewed but not reducible to only one or two dominant plants.

This is a contribution decomposition, not an independent plant-level effect analysis. A plant accumulates credit through the combination of being used by focal butterflies and adding introduced geographic units. It therefore answers **which resources generate the reconstructed opportunity**, not why those plants were redistributed.

The ten largest species-level contributors were:

| Host species | Family | Share of all added butterfly × region units |
|---|---|---:|
| *Medicago sativa* | Fabaceae | **4.55%** |
| *Poa pratensis* | Poaceae | **3.29%** |
| *Oxalis corniculata* | Oxalidaceae | **3.01%** |
| *Cynodon dactylon* | Poaceae | **2.60%** |
| *Dactylis glomerata* | Poaceae | **2.25%** |
| *Zea mays* | Poaceae | **2.03%** |
| *Avena sativa* | Poaceae | **2.03%** |
| *Urtica urens* | Urticaceae | **1.82%** |
| *Oryza sativa* | Poaceae | **1.80%** |
| *Digitaria sanguinalis* | Poaceae | **1.72%** |

Seven of the ten largest species-level contributors were Poaceae. This taxonomic composition is descriptive and is not used as a test of grass-associated redistribution.

### Genus-level sensitivity

To test whether species-level concentration could be an artifact of taxonomic splitting, contributing host species were collapsed to unique genera separately within every butterfly × WGSRPD3 added unit before fractional credit was assigned. The total again summed exactly to **14,553 units**.

Across **431 contributing genera**, the top 10 accounted for **28.69%**, the top 20 for **45.17%**, and the top 50 for **67.04%**. Only **25 genera** were required to reach half of all added opportunity. The genus-level Gini coefficient was **0.774**, HHI was **0.0143**, and the inverse-Herfindahl effective contributor number was **69.9 genera**. Taxonomic aggregation therefore retained, rather than erased, the strongly uneven contribution pattern.

| Genus | Family | Share of all added butterfly × region units |
|---|---|---:|
| *Medicago* | Fabaceae | **5.02%** |
| *Poa* | Poaceae | **3.57%** |
| *Senna* | Fabaceae | **3.46%** |
| *Oxalis* | Oxalidaceae | **3.01%** |
| *Cynodon* | Poaceae | **2.60%** |
| *Passiflora* | Passifloraceae | **2.49%** |
| *Avena* | Poaceae | **2.27%** |
| *Dactylis* | Poaceae | **2.26%** |
| *Zea* | Poaceae | **2.04%** |
| *Plantago* | Plantaginaceae | **1.99%** |

### Within-butterfly portfolio architecture is structurally constrained

Among 191 expanded butterflies in the host-taxonomy-adequate subset, the unadjusted association between host-family breadth and effective contributor number was **rho = 0.492**, while the association with maximum single-host share was **rho = -0.486**. After conditioning on exact resolved host-species richness, these fell to **0.059** and **-0.071**, with conditional permutation p-values of **0.502** and **0.415**. Joint adjustment for resolved host richness and butterfly Family gave **-0.039** and **0.015**.

These within-butterfly architecture metrics are therefore retained only as structural diagnostics. They are distinct from the across-plant concentration result above: the former asks how each butterfly's added units are partitioned among its hosts, whereas the latter asks which host plants account for the aggregate 14,553-unit expansion across the whole panel.

## Supplementary inference boundaries

1. Resource envelopes are reconstructed geographic opportunities, not realized local host use.
2. WCVP introduced status does not provide dates or causal introduction pathways.
3. Network degree cannot separate ecological host commonness from HOSTS recording intensity.
4. Occurrence recovery does not establish local larval use or host-caused colonization.
5. The occurrence panel was climate-stratified rather than designed as a validation sample.
6. Null-model, equivalence, host-contribution, ceiling, phylogenetic and geographic sensitivity analyses are post-hoc analyses and are interpreted as such.
