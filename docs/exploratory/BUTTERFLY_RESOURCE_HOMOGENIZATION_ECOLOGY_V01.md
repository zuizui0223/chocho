# Exploratory ecological reframing — resource expansion, homogenization and interaction exposure

Status: **EXPLORATORY / DO NOT MERGE INTO THE SUBMISSION MANUSCRIPT WITHOUT REVIEW**

## Frozen execution

- Workflow: `butterfly resource homogenization`
- Successful run: `37611759965`
- Head SHA: `dee4ee956e25e7e9a5d3d1e8fb7572040e475c51`
- Artifact: `11477849498`
- Focal butterflies: 239
- Added butterfly × WGSRPD3 resource units: 14,553

## 1. Individual opportunity expands while regional resource assemblages homogenize

Across 355 WGSRPD3 regions with native resource opportunity, mean pairwise Jaccard similarity of the butterfly resource assemblage increased from **0.27698** under native host distributions to **0.46208** under contemporary host distributions, a **66.8% relative increase**.

A fixed-margin swap null preserved both:

1. each butterfly's number of added resource regions, and
2. each region's number of added butterfly resource opportunities.

The null median contemporary similarity was **0.45349** (95% interval 0.45332–0.45364). The observed mean was **0.00860 higher** than this fixed-margin expectation (`p = 0.002`, 499 permutations).

Interpretation: most raw homogenization follows from the large amount of added opportunity, but the actual identity structure of host redistribution produces a small, detectable additional homogenizing effect beyond row and column totals.

## 2. Expansion magnitude is concentrated, but homogenization is not produced only by the top hosts

Thirty-eight of 670 contributing hosts provide 50.6% of fractional added opportunity.

Among these 38 hosts:

- 26 contribute to butterflies in at least two butterfly families;
- 8 contribute to at least three butterfly families;
- 2 contribute to four butterfly families;
- 71.5% of top-38 fractional credit comes from hosts contributing to at least two butterfly families.

Examples include `Medicago sativa` (17 focal butterfly species in four families) and `Robinia pseudoacacia` (seven species in four families).

However, homogenization is not solely a top-host phenomenon:

- native mean Jaccard: 0.27698
- adding only the top 38 introduced-host contributions: 0.40582 (Δ = +0.12884)
- adding only the remaining hosts: 0.40358 (Δ = +0.12660)
- all introduced-host contributions: 0.46208 (Δ = +0.18510)

Thus the **volume** of added opportunity is strongly concentrated, whereas the **homogenization** signal is distributed across both dominant and long-tail hosts.

## 3. Human redistribution greatly expands potential interspecific resource-sharing space

Define a potential resource-sharing unit as a butterfly-pair × WGSRPD3 combination in which both butterflies are documented to use at least one exact same host species present in that region.

- native shared-resource pair × region units: **62,473**
- contemporary shared-resource pair × region units: **141,885**
- novel shared-resource pair × region units created by introduced host ranges: **79,412**
- relative increase: **+127.1%**

Of 986 butterfly pairs that share resource geography, **920** gain at least one novel shared-resource region.

Novel shared-resource units are not confined within butterfly families:

- within-family: 37,869
- between-family: **41,543 (52.3%)**

The average number of other focal butterflies sharing at least one exact host in the same regional resource opportunity is:

- native resource units: mean **4.71**, median **3**
- introduced-added resource units: mean **5.71**, median **4**

The fraction with no other focal host-sharing butterfly falls from **18.3%** in native resource units to **9.9%** in introduced-added units.

### Boundary

These quantities measure **potential resource co-use / interaction exposure**, not realized competition. They do not contain abundance, density, host quality, phenology, parasitoid data or direct demographic responses.

They can motivate hypotheses about exploitative competition, plant-mediated interactions and apparent competition, but cannot establish any of those mechanisms.

## 4. Counterfactual host-removal stress test

Introduced-range contributions of host plants were removed in descending order of their fractional contribution to the 14,553 added butterfly × region units, while native host ranges were retained.

| Removed top hosts | Added units lost | Fraction lost | Butterflies affected | Lose ≥50% of added opportunity | Lose all added opportunity |
|---:|---:|---:|---:|---:|---:|
| 1 | 446 | 3.1% | 14 | 3 | 0 |
| 5 | 1,816 | 12.5% | 40 | 17 | 6 |
| 10 | 2,893 | 19.9% | 61 | 24 | 7 |
| 20 | 4,199 | 28.9% | 93 | 36 | 9 |
| 38 | **5,812** | **39.9%** | **123** | **56** | **16** |
| 50 | 6,615 | 45.5% | 133 | 66 | 20 |
| 100 | 8,781 | 60.3% | 164 | 101 | 33 |

For comparison, removing 38 randomly chosen contributing hosts loses a median of only **432** added units in 499 randomizations.

This is a **network stress test**, not a recommendation to remove or retain any plant. The top hosts were selected because they contribute strongly, so the targeted-versus-random contrast is descriptive of concentration/fragility rather than an independent causal test.

## 5. Ecological interpretation worth pursuing

A stronger ecological framing is an **expansion–homogenization trade-off**:

> Human redistribution of host plants expands geographic resource opportunity for individual butterflies while simultaneously homogenizing the resource landscape and increasing the spatial overlap of potential resource use among butterfly species.

This produces three distinct predictions:

1. **Range-expansion prediction.** Where climate and dispersal permit, butterflies should be more likely to colonize regions in which introduced known hosts remove a previous resource barrier.
2. **Interaction-overlap prediction.** Newly created resource regions should expose butterflies to more host-sharing species, increasing opportunities for direct resource competition, plant-mediated effects or shared-enemy apparent competition.
3. **Management-sensitivity prediction.** Butterflies whose contemporary added opportunity is concentrated on a small number of introduced hosts should be especially sensitive to host removal unless suitable native resources are restored locally.

The second prediction is consistent with ecological theory showing that host sharing can structure butterfly co-occurrence and that phytophagous insects may interact through both resource competition and shared natural enemies. It remains a prediction here because abundance and natural-enemy data are absent.

## 6. Conservation use

The current data support **risk screening**, not site-level prescriptions.

Useful conservation outputs from the next analysis are:

- butterfly species with high introduced-resource dependence and low host redundancy;
- WGSRPD3 regions with high novel host-sharing overlap;
- host plants whose loss removes unusually large amounts of resource opportunity;
- regions/species where introduced-host removal should trigger locality-level validation and native-host restoration planning.

The management principle is therefore:

> Do not ask only whether an introduced plant is present or harmful in general. Ask which native consumers currently depend on it, how redundant that resource is, and whether removal occurs before or after suitable native resources are restored.

## 7. Next analysis

The active branch additionally tests:

- normalized butterfly **host × region resource-niche Jaccard overlap** before vs after host redistribution;
- the fraction of added butterfly × region opportunities supported by only one introduced host;
- species-level low-redundancy dependence;
- targeted versus random host-loss robustness.

If normalized niche overlap also increases, the ecological claim can be sharpened from regional resource homogenization to **consumer resource-niche convergence**.


## 8. Normalized species-pair resource-niche overlap: a useful boundary

The latest exact run (`37612135885`, head `c1eaa7eacd12e766df52fea7b910fe34d7321084`) additionally represented each butterfly niche as its set of exact **host species × WGSRPD3** resource cells and calculated pairwise Jaccard overlap across all 28,441 butterfly pairs.

- mean native Jaccard: **0.004032**
- mean contemporary Jaccard: **0.004885**
- mean change: **+0.000853**
- pairs increasing: **616**
- pairs decreasing: **366**
- pairs unchanged: **27,459**
- pairs gaining at least one newly shared host × region cell: **920**
- of these, **450 (48.9%)** are cross-family pairs

This result prevents an overclaim. Anthropogenic host redistribution does **not** produce a wholesale convergence of the resource niches of all butterflies. Most butterfly pairs use different hosts and remain resource-disjoint.

The stronger, biologically defensible interpretation is:

> Host redistribution greatly expands the **geographic interaction opportunity among butterflies that already share hosts**, while regional resource assemblages become more similar overall.

Thus “resource-niche convergence” should not be used as the headline without qualification. “Regional resource homogenization” and “expanded interspecific resource-sharing geography” are better supported.

## 9. Added opportunity is often locally non-redundant in the reconstructed network

Across the 14,553 butterfly × region units added by introduced host ranges:

- **8,574 (58.9%)** are supported by exactly **one** introduced known host;
- median number of contributing introduced hosts per added unit = **1**.

Several butterflies have all of their reconstructed added geography supported by a single introduced host at each added region, including high-magnitude cases such as:

- `Pseudozizeeria maha`: 259/259 added regions single-host supported;
- `Ypthima asterope`: 133/133;
- `Lasiommata megera`: 126/126;
- `Oeneis norna`: 101/101;
- `Carcharodus alceae`: 91/91.

This is a **known-host network redundancy** result, not proof of demographic dependence. HOSTS incompleteness and unrecorded local hosts can only make apparent redundancy lower than reality.

Nevertheless, it motivates a practical field-validation priority:

> high introduced-resource gain + low reconstructed host redundancy = high priority for checking local host use and native replacement resources before invasive-host removal.

## 10. Competition simulation decision

A fully parameterized Lotka–Volterra or MacArthur consumer–resource model is **not yet justified** by these data because local host biomass, butterfly abundance, host quality, phenology and parasitoid/predator networks are absent.

The defensible sequence is:

1. empirical structural result: regional resource homogenization;
2. empirical structural result: >2× expansion of shared-resource geography among host-sharing butterfly pairs;
3. empirical robustness result: 58.9% of added units are single-host supported;
4. counterfactual network stress test: targeted host loss disproportionately contracts added resource opportunity;
5. **then**, only as a sensitivity analysis, use a transparent fixed-abundance/null allocation model to ask how resource expansion could redistribute potential crowding from within species to among host-sharing species.

Any competition simulation must be labeled a scenario model rather than a fitted demographic prediction.


## 11. Crop and grass sensitivity: homogenization is not a crop/Poaceae artifact

A completed sensitivity run (`37613267331`) reconstructed the resource matrices after excluding major anthropogenic host categories.

### Regional resource-assemblage similarity

| Reconstruction | Native mean Jaccard | Contemporary mean Jaccard | Relative increase |
|---|---:|---:|---:|
| Baseline | 0.2770 | 0.4621 | **+66.8%** |
| Exclude FAO crop binomials | 0.2704 | 0.4279 | **+58.2%** |
| Conservative crop-genus exclusion | 0.2646 | 0.4180 | **+58.0%** |
| Exclude all Poaceae hosts | 0.2600 | 0.4327 | **+66.5%** |

### Butterfly-pair resource-geography overlap

| Reconstruction | Native mean Jaccard | Contemporary mean Jaccard | Relative increase |
|---|---:|---:|---:|
| Baseline | 0.2224 | 0.3248 | **+46.0%** |
| Exclude FAO crop binomials | 0.2126 | 0.2927 | **+37.7%** |
| Conservative crop-genus exclusion | 0.1942 | 0.2646 | **+36.2%** |
| Exclude all Poaceae hosts | 0.1310 | 0.1906 | **+45.6%** |

Thus the resource-homogenization signal is not reducible to crop hosts or redistributed grasses. These exclusions are sensitivity reconstructions, not causal estimates of agriculture or Poaceae.

## 12. Minimal fixed-demand prediction for intraspecific crowding

Before fitting any competition dynamics, a parameter-light null prediction can be obtained from resource geography alone.

Assume each butterfly species has a fixed total larval demand of one unit and allocates that demand uniformly across all available resource regions. The probability that two independent units of conspecific demand fall in the same resource region is then proportional to (1/R_i), where (R_i) is the number of resource regions available to species (i).

Across 239 butterflies:

- mean (1/R_i), native resource geography: **0.02126**
- mean (1/R_i), contemporary resource geography: **0.01416**
- relative change in the mean concentration proxy: **−33.4%**
- median species-level change: **−28.0%**
- 206/239 species show lower concentration; 33 are unchanged.

This creates a clear field-testable prediction:

> If total abundance and resource quality were held constant, host redistribution should reduce within-species spatial concentration of larval demand while increasing the geography over which host-sharing species can potentially interact.

The first half is a null allocation prediction, not an observed demographic effect. The second half is supported structurally by the >2× increase in exact shared-resource pair × region opportunities.

## 13. Conservation prediction framework

The exploratory results motivate two independent screening axes:

1. **introduced-resource dependence** — how much added geography is supported by one or a few introduced known hosts;
2. **novel interaction exposure** — how strongly the new resource geography overlaps host-sharing butterflies.

This generates four useful field-validation classes:

- low dependence / low exposure: comparatively buffered;
- high dependence / low exposure: removal-sensitive resource opportunity;
- low dependence / high exposure: interaction-rich but resource-redundant;
- high dependence / high exposure: highest priority for local validation of host use, larval performance, competitors and shared natural enemies.

The framework does not produce site-level management decisions from WGSRPD3 data. It identifies where fieldwork would most efficiently distinguish a true life-raft effect from competition, apparent competition, ecological traps or unused potential resources.
