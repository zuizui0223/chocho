# Manuscript-ready ecological narrative — introduced larval-host geography

**Status (16 October 2026): ALTERNATIVE SINGLE-PAPER FRAMING / DISCUSSION ONLY.** This file proposes a revised scientific narrative for the existing 239-butterfly GEB paper; it is **not** a second manuscript, a change to the frozen submitted draft, a preregistered hypothesis or proof of population establishment. The earlier public manuscript and original sources remain unchanged.

## Biological question and novelty

**Question:** When humans move the known larval host plants of butterflies into new regions, how much of that geographic resource supply is genuinely *new to the butterfly*, and how much falls within the butterfly's already available native-host resource envelope?

The natural expectation “more host plant range produces more mapped opportunity” is essentially set-theoretic. The nontrivial test is whether the introduced host×consumer×region arrangement is **geographically more coincident with existing native resource access than a matched spatial rearrangement would produce**, and whether that remains true when plant taxonomic group, destination regional richness, broad climate and geographical proximity are constrained.

This explicitly concerns *resource geography*, not observed plant–larva edges in the introduced region or butterfly population demography. Before asserting novelty, distinguish from (i) plant–frugivore meta-network homogenization (Fricke & Svenning 2020, DOI 10.1038/s41586-020-2640-y), and (ii) local empirical plant–herbivore network topology under exotic plant prevalence (Araújo et al. 2024, DOI 10.1016/j.scitotenv.2024.173949). The distinguishing design here is keeping known butterfly–host link identities fixed and reconstructing **where potential resources overlap or expand**.

## Recommended English title and abstract

**Working title:** *Plant globalization often overlaps existing butterfly larval-resource geography rather than creating novel regional opportunity*

**Alternative (more neutral):** *Geographic novelty and redundancy of butterfly larval-host resources under plant globalization*

### Abstract — revised for ecological specificity

Human introductions redistribute host plants across biogeographic regions, but geographical expansion of hosts need not translate into proportionally novel resources for consumers. We reconstructed native and contemporary larval-host resource geography for **239 butterflies** across **365 WGSRPD3 regions**, fixing recorded butterfly–plant host associations. We distinguished introduced host–butterfly–region incidences that overlap already available native-host resources, duplicate introduced-host support in newly opened resource regions, and support distinct new butterfly–region opportunities. Of **75,346** candidate introduced-host incidences, **46,387 (61.6%)** were superimposed on regions where the butterfly already had native-host resource access. A further **14,406 (19.1%)** represented additional introduced-host supports within newly supplied regions, leaving **14,553 distinct butterfly–region additions (19.3% of the incidence budget)**. In post-hoc plant-family-, continent-, climate- and local-distance-constrained redistribution nulls, the observed placement remained modestly more coincident with the native resource envelope. Under the fully evaluated **1,500-km** spatial sensitivity, the observed native-overlap count exceeded the two null medians by **271–276 host incidences**, while distinct additions fell **264–267 butterfly–region units** below null medians. These deviations were small (overlap excess **0.36 percentage points** of candidate incidence volume), and butterfly host-family breadth was only weakly associated with proportional resource expansion (**Spearman rho = 0.008**). Introduced hosts therefore mostly supplied geographically overlapping opportunities rather than distinct new mapped regions within this focal resource network. The results separate geographic resource novelty from range expansion but do not establish local larval host use, butterfly colonization, fitness, competition or causality of the observed overlap.

**Word-choice caveat:** “Introduced hosts therefore mostly supplied geographically overlapping opportunities” describes the union-map reconstruction and incidence accounting. “More than expected” is *specific to the named post-hoc spatial randomization*. It must **not** imply that plants prefer butterfly-rich regions, or that individuals feed successfully on introduced plants.

## Main Results: proposed order

### 1. Widespread larval-resource increase masks considerable overlap

**Panel scope:** 239 butterflies; 365 geographically delineated regions; fixed HOSTS records and WCVP v13 native/introduced statuses; 670 plant species supply at least one distinct new resource unit.

Observed native butterfly×region resource cells **26,530**; contemporary **41,083**, an increase of **14,553 (+54.9%)**. Species gaining any extra region: **206/239 (86.2%)**. These are the scale and context, **not** the lead surprise.

To unpack the geographic novelty of candidate introduced host distribution:

\[
H = A + N,\qquad N = U + R
\]

where `H=75,346` introduced host×butterfly×region incidences, `A=46,387` incidences falling in regions already within the butterfly's native larval-resource envelope, `N=28,959` incidences falling into novel butterfly×region resource cells, `U=14,553` distinct novel butterfly×region cells and `R=14,406` extra supporting host incidences beyond one per new resource cell. **H=A+R+U is a bookkeeping identity with one incidence counted as the first support of each new cell, not a ratio of population growth or resource biomass.**

The last two terms are not independent: highly overlapping introduced plants can supply multiple supporting incidences within one new cell. None of these counts measures edible plant biomass, insect acceptance or reproductive success.

### 2. Introduced ranges are more likely to coincide with native-resource geography than the constrained structural null would predict — but only slightly

The initial null preserves planted-host degree and continental plant-family region margins. Of the fixed 75,346 candidate incidences, observed native-envelope overlap is **46,387**, versus **45,301/45,298** in two 499-draw null chains (**+1,086/+1,089**).

CHELSA temperature/precipitation-bin matching reduces the excess to **+552/+558**. A more stringent *locally tethered* comparison (1,500 km, original introduced host footprint, the same plant-family × continent × 2×2 climate bins, two independent longer chains) retains **+271/+276**, or **+0.36 percentage points** of the fixed candidate incidence budget. Independent-chain simulation diagnostic thresholds were satisfied: autocorrelation lag-1 ≈ 0.06–0.10; estimated ESS ≈ 75–85 for the primary A response; chain median difference 5 units. Distinct new resource additions remain **−264/−267** lower than the respective null medians (**14,817/14,820**).

These are **small conditional differences** in a coarse mapped resource system. Long-chain numerical qualification does not convert a post-hoc comparison into a confirmatory permutation test or into identification of climate, geography or introduction mechanisms.

The short-chain local null (199 serial draws per seed) had lag-1 around 0.78 at 1,500 km and ESS ≈ 9–10, and its nominal p-values were rejected as evidentially unreliable. The longer-chain sensitivity was run once with a finite frozen mixing gate; results were stable, not chosen by repeated threshold fishing.

### 3. Specialization does not capture geographic novelty well

Across focal butterflies, the association between host-family breadth and proportional introduced resource gain is **rho=0.008 (bootstrap 95% CI −0.111 to 0.128)**. This does **not** prove a mathematically exact zero effect; a post-hoc narrow ±0.10 equivalence margin was not met (approximate TOST p=0.078). Host-family breadth and native geographic resource breadth are only partially correlated (**rho=0.276** full sample). In the more host-taxonomy-adequate 215-species subset, **18 one-family host specialists** nevertheless occupy the top quartile of geographic resource opportunity.

For the introduced-added cells themselves, **8,574/14,553 (58.9%)** are supported by a single introduced host, but the proportion is not unusually high relative to a plant-range null with relevant margins. More host species typically allow more redundant added opportunities; **taxonomic host-family breadth by itself does not capture that spatial portfolio structure**.

### 4. Evidence that resource opportunity is biologically relevant is partial, and conservation claims are deliberately bounded

A **secondary, non-prospective** 32-species climate-stratified GBIF panel included **115** butterfly outside-native resource-region occurrence units from 23 informative species; **66** coincided with introduced host-resource gains versus a geography/region-matched null median **49**. Species sensitivity and uneven GBIF coverage warrant careful reporting; the result is **not** proof of introduced-host larval use or successful colonization.

Among the reconstructed regions, regional butterfly **potential resource assemblage** Jaccard similarity increased **0.277→0.462** (+66.8%), but a null fixing both butterfly-specific and region-specific additions predicted about **0.453**; the residual is only **+0.0086**, falling to **+0.0036** under continental margin controls. It is inaccurate to transfer the raw +66.8% to *observed butterfly community homogenization*, or claim the new shared-host exposure caused competition. Likewise, 920 of 986 host-sharing butterfly pairs acquire new shared-host regions, even though **at least 300 pairs** must have gained shared-host geography while their normalized host×region overlap declined; these are distinct quantities.

## Discussion — central implication

The primary ecological distinction is *access versus additional access*, not area versus population growth. Human redistribution creates many new potential plant occurrences, but the same introduced plant can (a) lie in a region with existing native host resources, (b) supply a new region alongside several other introduced plants, or (c) be the first introduced host to supply an otherwise unmatched butterfly×region resource cell. Conventional resource-envelope expansion reports only the union, making these different ecological situations look alike. The data show that most candidate supply maps to (a), with the observed arrangement producing a modest excess of (a) over plant-range alternatives that also preserve broad climate and source-geographic locality.

The remaining overlap could reflect spatial autocorrelation, climatic niches at finer scales, land use, introduction history, botanical recording, host-database incompleteness or correlated native plant richness. These are not statistically eliminated. Importantly, no population success or reduction in fitness can be inferred without local plant abundance, larval performance, resource phenology, and verified host use.

## Final figure story (proposed, *not* currently adopted in submission)

- **Figure 1:** Global analysis overview, exact 239/365/75,346 scope and the **A+R+U** resource novelty accounting. Native, added and contemporary species×region totals shown as scale, not as headline inference.
- **Figure 2:** Reconstructed observed-minus-null native-resource absorption and new-cell deficit under plant-family/continent, plus climate, plus qualified 1,500-km locality controls. Show each chain separately; source receipts and the method's ESS/gates in caption. This figure should make the decreasing residual visible rather than exaggerate a small effect.
- **Figure 3:** Butterfly heterogeneity and biological relevance: specialist versus geographic resource-breadth discordance, plus outside-native observed butterfly occurrence recovery against a region-matched structural null. State separate sample sizes and biases.
- **Supplement:** Native-only versus contemporary reconstructed potential resource-assemblage Jaccard, matched structural null and exact shared-host geography, with negative-result interpretation of fragility and competition claims.

Submission-ready title/abstract and plot storyboard are a *candidate author decision*. Do not change the GEB submission bundle or silently merge PR #43 into PR #38.

## Source hierarchy and checks

- All reconstruction and null inputs remain pinned to the exact commits and source artifacts documented in `docs/REPRODUCIBILITY.md`.
- Base geographic null: `provenance/reviewer_defenses/results/butterfly_plant_range_structural_null_two_chains_v0.2.json`.
- Climate sensitivity: `provenance/reviewer_defenses/results/butterfly_climate_matched_absorption_v0.1.json`.
- Numerically qualified 1,500-km local null: `provenance/reviewer_defenses/results/butterfly_spatial_null_longchain_qualified_v0.2.json`.
- Figure generator: `scripts/render_butterfly_resource_novelty_figure.py`.
- Plant–frugivore global interaction homogenization: Fricke & Svenning (2020), *Nature*, DOI 10.1038/s41586-020-2640-y.
- Exotic plants and empirical plant–herbivore network topology: Araújo et al. (2024), *Science of the Total Environment*, DOI 10.1016/j.scitotenv.2024.173949.

**Final stop rule:** There is no need to keep searching new predictors until a larger post-hoc result appears. The supported conclusion is source-conditional, effect-size-qualified resource geography. If the paper requires a mechanism, it needs *independent plant/life-stage or temporal evidence*, not more re-shuffled versions of the same host map.
