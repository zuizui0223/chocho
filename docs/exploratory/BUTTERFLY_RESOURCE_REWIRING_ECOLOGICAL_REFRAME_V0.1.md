# Ecological reframe: anthropogenic rewiring of butterfly resource landscapes

Status: exploratory synthesis; do not promote to the submission manuscript until the focused connectivity and dependency tests pass.

## Why the previous headline is insufficient

“Introduced hosts expand reconstructed resource geography” is partly structural: contemporary resource envelopes are native envelopes plus introduced host distributions. The ecological question must therefore concern **how the added resources are organized and what ecological consequences that organization creates**, not the direction of expansion itself.

## Candidate ecological question

**Does plant globalization merely add larval resources, or does it reorganize the spatial structure, sharing and robustness of butterfly resource niches?**

This separates four ecological consequences:

1. **resource homogenization** — do geographically distinct regions become similar in which butterflies have larval resources?
2. **niche overlap / interaction exposure** — do butterflies that share exact hosts gain more geographic overlap in which those resources are jointly available?
3. **resource-network robustness** — are new opportunities redundant across hosts or dependent on a small number of redistributed plants?
4. **range-frontier opportunity** — where do introduced hosts create climate-compatible, effort-supported resource regions that remain unobserved for the butterfly?

## Verified results

### 1. Regional larval-resource assemblages strongly converge

Across 355 WGSRPD3 regions with native resource opportunity, mean pairwise Jaccard similarity among regional butterfly-resource assemblages increased from **0.27698 to 0.46208**, a **66.8% relative increase**.

This increase is largely expected from the amount and placement of added opportunity, but not completely. A fixed-margin switch null preserving each butterfly’s number of added regions, each region’s number of added butterfly opportunities, and all native incidences gave a median contemporary Jaccard of **0.45347**; observed similarity was higher (**p = 0.002**).

A stricter null additionally preserving, for every butterfly, the number of added regions within each WGSRPD level-1 region gave a median of **0.45845** (95% null interval 0.45836–0.45855); observed remained higher by 0.00363 (**p = 0.002**).

Interpretation: most homogenization follows from the broad redistribution process itself, while the actual host-redistribution configuration is modestly but consistently more homogenizing than matched alternatives.

### 2. Butterfly resource niches converge as well

Mean pairwise overlap of butterfly resource geography increased from **0.22237 to 0.32788** (+47.4%) in the full reconstruction. Under the level-1-constrained fixed-margin null, expected contemporary overlap was **0.32516**; observed was again higher (**p = 0.002**).

This is potential resource-niche overlap, not realized competition.

### 3. Existing host-sharing relationships are geographically amplified

For butterfly pairs documented to use at least one exact same host species, shared-resource butterfly-pair × region units increased from **62,473 to 141,885**. Introduced host ranges therefore created **79,412 novel shared-resource pair × region units** (+127.1%).

The number of pair identities did not increase: the same **986** resource-sharing butterfly pairs existed in native and contemporary reconstructions. Instead, **920/986 pairs** gained at least one new shared-resource region. Of the novel pair × region units, **41,543 (52.3%)** linked butterflies from different families.

This is best described as **spatial amplification of existing trophic overlap**, not creation of novel trophic links.

### 4. New resource geography is more crowded by exact-host co-users

Among native butterfly × region resource opportunities, a focal butterfly shared at least one exact host with a mean **4.71** other focal butterflies (median 3). In introduced-added opportunities, the corresponding mean was **5.71** (median 4). Resource units with zero exact-host co-users declined from 18.3% of native opportunities to 9.9% of introduced-added opportunities.

This directly establishes greater **potential resource-sharing exposure** in added geography. It does not establish realized competition or apparent competition.

### 5. Resource gains are large but weakly redundant

Of the 14,553 added butterfly × region units, **8,574 (58.9%)** were supported by only one introduced host species in that butterfly × region combination; the median number of contributing hosts per added unit was one.

Removing introduced-range contributions of the top-ranked hosts in a global counterfactual stress test produced strongly concentrated losses:

- top 5 hosts: 1,816 units lost (12.5%); 40 butterflies affected; 6 lose all added opportunity;
- top 10: 2,893 units (19.9%); 61 affected; 7 lose all;
- top 38 (the set supplying half of fractional contribution): **5,812 units (39.9%)**, 123 butterflies affected, 56 lose at least half of their added opportunity and 16 lose all;
- random removal of 38 contributing hosts loses a median **433** units (95% 168.9–1123.5); targeted loss is above all 499 randomizations (**p = 0.002**).

Thus contribution concentration translates into **network fragility**, not merely a descriptive rank curve.

### 6. Opportunity and potential competition can move in opposite directions

A transparent sensitivity model down-weighted each butterfly × region opportunity by the number of other focal butterflies sharing an exact host there. The penalty parameter is hypothetical and is **not** estimated from ecological data.

Under reciprocal weighting, aggregate effective opportunity remains greater than native opportunity across the tested range, but species-level reversals accumulate: 11 species have negative effective gain at alpha=0.1, 26 at 0.25, 42 at 0.5 and 56 at 1. The aggregate reciprocal-model break-even occurs only at alpha = **9.55**. Under an exponential penalty, the aggregate break-even occurs at alpha = **1.51**, after which crowding can outweigh the nominal geographic gain.

Species-level break-even thresholds are heterogeneous. Under the reciprocal model, 98/206 expanded species have a finite threshold within the search range (median alpha = **0.666**), whereas 108 never reverse within the search. Under the exponential model, 157/206 have a finite threshold (median alpha = **0.502**). These thresholds are dimensionless and model-dependent; they quantify how strong crowding costs would need to be to erase the reconstructed gain, not measured competition coefficients.

This is a mechanistic **sensitivity boundary**, not evidence that competition actually cancels resource gain.


### 6b. The biggest resource gains also carry the biggest crowding increases

Across 239 butterflies, proportional resource gain was positively associated with the increase in mean exact-host co-user exposure (**Spearman rho = 0.317; 19,999-bootstrap 95% CI 0.187–0.438**). The relationship remained positive in rank-residual sensitivities controlling native resource breadth (~0.397), native resource breadth plus host-family breadth (~0.397), or native resource breadth, host-family breadth and native co-user exposure together (~0.372).

Thus the resource and interaction consequences of plant globalization are not independent. Butterflies receiving larger proportional expansions of host-resource geography also tend to enter geography in which more other focal butterflies share their exact hosts.

This gives a sharper ecological trade-off than either result alone:

> **the largest apparent resource winners also acquire the largest increases in potential interspecific resource-sharing exposure.**

The association is not evidence of demographic competition. It shows coupling between reconstructed resource release and a spatial proxy for potential competitive interaction.

### 7. A falsifiable range-frontier prediction set exists

Within the independently constructed 24-species climate panel, **806 species × region units** have:
- introduced-only reconstructed host opportunity,
- no butterfly occurrence in the existing panel,
- sufficient other-butterfly sampling effort,
- an estimable climate mismatch.

Of these, **159** have climate compatibility rank >=0.5 and **44** >=0.8 relative to held-out observed regions.

These should be frozen as **prospective detection/colonization-frontier candidates**, not described as guaranteed future colonizations. Examples among the highest-ranked, high-effort candidates include *Melanitis leda* in Brazil Southeast and *Biblis hyperia* in Florida.


## Candidate ecological abstract

**Aim:** Human transport redistributes plant resources across biogeographic barriers, but whether this merely enlarges herbivore resource ranges or reorganizes the spatial structure and robustness of consumer niches is unclear. We asked how globalization of known larval host plants changes butterfly resource landscapes, resource sharing among consumers and dependence on particular redistributed hosts.

**Methods:** For 239 butterflies, we held known larval-host identities fixed and reconstructed native and contemporary host-resource geography across WGSRPD3 regions. We quantified regional resource-assemblage homogenization, butterfly resource-niche overlap, geographic amplification of exact-host sharing, redundancy of newly added resource regions and counterfactual sensitivity to loss of high-contribution introduced hosts. Fixed-margin nulls preserved butterfly-specific, region-specific and broad continental amounts of added opportunity.

**Results:** Human host redistribution increased aggregate butterfly × region resource opportunity by 54.9%, but its stronger ecological effect was structural. Mean similarity among regional butterfly-resource assemblages rose by 66.8% and butterfly resource-geography overlap by 47.4%; both exceeded level-1-constrained fixed-margin expectations (p = 0.002). Existing exact-host-sharing relationships were geographically amplified: shared-resource butterfly-pair × region units increased from 62,473 to 141,885 (+127.1%), with 920/986 host-sharing pairs gaining new shared-resource regions. Added resource opportunities were more crowded by exact-host co-users than native opportunities (mean 5.71 versus 4.71 other butterflies). At the same time, 58.9% of added butterfly × region opportunities were supported by only one introduced host in that butterfly-region combination. Removing the 38 highest-contribution hosts eliminated 39.9% of all added opportunity, compared with a median 3.0% under random removal of 38 hosts (p = 0.002).

**Main conclusion:** Plant globalization does more than enlarge butterfly resource ranges. It spatially homogenizes larval-resource landscapes, amplifies pre-existing trophic overlap and creates geographically extensive but often weakly redundant resource dependencies. Thus anthropogenic redistribution can simultaneously relax resource limitation and increase potential resource-sharing exposure, creating a trade-off between opportunity, interaction exposure and robustness that is invisible from consumer diet breadth alone.

## What is genuinely non-trivial

The paper should not sell “resources increased.” The non-trivial pattern is:

> **Human plant redistribution simultaneously enlarges and homogenizes butterfly larval-resource landscapes, geographically amplifying resource sharing among consumers while concentrating much of the new opportunity in weakly redundant host dependencies.**

This creates a biologically interesting tension:

**resource release ↔ resource-sharing exposure ↔ dependency/fragility**

The data directly support the first, the geography of the second, and the structural fragility of the third. Realized competition, host quality, abundance, enemy-mediated apparent competition and demographic benefit remain unmeasured.

## Candidate main hypotheses

H1 — **Homogenization:** Actual host redistribution makes regional butterfly-resource assemblages more similar than matched redistribution with the same butterfly, region and broad-continent margins.

H2 — **Interaction exposure:** Added geography is more densely shared among exact-host co-users than native resource geography, spatially amplifying pre-existing consumer overlap.

H3 — **Fragility:** Because added opportunity is concentrated and weakly redundant, removing high-contribution introduced hosts eliminates far more resource opportunity than random host removal.

H4 — **Connectivity (technical rerun pending):** Introduced resource regions act as stepping stones that connect native resource components more than level-1-matched random additions. The first execution failed only because the workflow omitted the Shapely dependency; the workflow has been repaired and the ecological test itself is unchanged.

## Conservation interpretation

The management implication is not “retain invasive plants.” Instead, the analyses can identify **where an introduced plant may sit inside a low-redundancy consumer resource network**. Removal or restoration decisions should then prioritize local validation of larval use, host quality, abundance, alternative native hosts and butterfly demography.

At WGSRPD3 resolution these are screening hypotheses and prioritization tools, not site-level prescriptions.

## Literature position

Global introductions are already known to homogenize species composition and plant–frugivore interaction networks (e.g. Fricke & Svenning 2020), and recent work has examined exotic plants in local or multi-site plant–herbivore networks. Butterfly work also shows that host availability and host sharing can structure geographic co-occurrence. The distinct contribution here is the **resource-side spatial rewiring of a herbivore interaction network while consumer host identities are held fixed**: host globalization changes where pre-existing trophic relationships can potentially occur, amplifies their geographic overlap and exposes the resulting network to concentration and redundancy tests without requiring host switching.

## Candidate title direction

**Plant globalization homogenizes butterfly resource landscapes and amplifies shared-resource exposure**

Alternative, more conservative:

**Anthropogenic host redistribution homogenizes butterfly larval-resource geography**

Do not use “competition” in the title unless abundance or demographic competition is independently demonstrated.
