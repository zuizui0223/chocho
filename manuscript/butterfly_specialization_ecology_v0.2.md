# Anthropogenic host redistribution expands butterfly resource geography beyond host-richness and plant-family expectations

**Working manuscript v0.2 — 2026-09-28**

## Abstract

**Aim:** To test whether human redistribution of larval host plants creates non-random geographic resource opportunity for butterflies, whether this opportunity is realized in contemporary butterfly records, and whether taxonomic generalism predicts its magnitude or release from climatic filtering.

**Location:** Global.

**Time period:** Contemporary host and butterfly distributions; butterfly occurrences 2010–2026; CHELSA climatology 1981–2010.

**Major taxa studied:** Butterflies and their larval host plants.

**Methods:** For 239 butterflies, we reconstructed native and contemporary larval-resource geography from LepTraits, HOSTS and WCVP. We compared observed host portfolios with 999 same-plant-family, same-host-count random portfolios per species, thereby separating host identity from the arithmetic effects of host richness and family composition. In an independent 32-species butterfly panel, we tested whether introduced-added host regions recovered modern occurrences lying outside native host envelopes, with an exact geographic-overlap null. Finally, we evaluated climate-associated filtering in 24 quality-qualified species and repeated the comparison within the same WGSRPD level-1 regions.

**Results:** Introduced host distributions expanded reconstructed opportunity for 206/239 species (86.2%) and increased aggregate resource geography by 54.9%. In the conservative 215-species subset, observed host portfolios exceeded the matched host null by a median 0.140 log units (156/215 species above the null median); this excess was almost unrelated to host-family breadth (Spearman rho = 0.041). Even one-family butterflies showed positive excess expansion (median 0.111; 63/97 above the null median). In the independent occurrence panel, introduced-added host regions recovered 66 of 115 butterfly species × region observations outside native host envelopes, versus 38.0 expected under the geographic-overlap null (p = 1.1 × 10^-11). Climate-associated filtering remained strong within the same broad geographic regions (median score = 0.798; 21/23 species >0.5), but broader family-level diets did not show weaker filtering.

**Main conclusions:** Anthropogenic plant redistribution changes butterfly resource geography through the identities of hosts that butterflies use, not simply through how many host species or families they use. Documented butterfly host portfolios are disproportionately composed of plants whose introduced ranges add geographic opportunity, and those added regions are enriched for contemporary butterfly occurrences. Taxonomic generalism does not confer a clear proportional advantage, and host availability remains further filtered by climate-associated constraints.

**Keywords:** butterflies; biological invasions; ecological specialization; global change; host plants; introduced plants; resource geography


---

## 1. Introduction

Specialization is one of the most familiar axes used to compare herbivorous insects. Butterflies are routinely described as specialists or generalists according to the number or phylogenetic breadth of larval host plants they use, and diet breadth has been linked to geographic range size, range dynamics, diversification, and environmental gradients (Slove & Janz 2011; Lancaster 2020; Gross et al. 2026). Yet a single taxonomic measure of diet breadth conflates several ecological quantities that need not vary together.

A butterfly that feeds on a single plant family may use one host species or dozens. Those host species may themselves be geographically restricted or nearly cosmopolitan. Consequently, taxonomic host breadth need not equal the geographic breadth of larval resources. Under global change, this distinction becomes dynamic because host plants themselves are redistributed. Regional work has long shown that introduced hosts can support butterfly range expansion, persistence, or phenological change (Graves & Shapiro 2003), while recent global work documents extensive human-mediated redistribution of Lepidoptera themselves (Couto et al. 2026). A complementary resource-side problem remains unresolved: across butterfly species, how much does anthropogenic host redistribution expand potential larval-resource geography, does that proportional gain favor broad generalists, and can similar aggregate gains be assembled through different host portfolios?

The geographic consequence of plant redistribution depends on host identity as well as host number. A butterfly using one exceptionally widespread introduced plant could acquire more new resource geography than a nominal generalist whose many hosts remain largely within native ranges. Conversely, simply drawing more host species mechanically increases the chance of including a widely redistributed plant. Separating those possibilities requires a null model that fixes host richness and plant-family composition while changing the identities of host species.

A further ecological layer separates potential resource opportunity from realized species geography. Even when suitable larval hosts are present, butterflies may fail to occupy parts of that opportunity because of climate, dispersal, phenology, habitat structure, biotic interactions, or incomplete colonization. Recent work shows both extensive climate-associated butterfly range shifts and important host constraints on thermal niche margins, while global analyses increasingly separate climate effects from functional traits and resource specialization (Rashid et al. 2026; Chowdhury et al. 2026). These findings motivate an explicit decomposition: hosts define one component of geographic opportunity, whereas climate and other processes filter which portions become realized butterfly distributions.

Here we reconstruct these layers within a single butterfly ecology framework. Our central question is not merely whether introduced hosts enlarge a potential resource envelope, but whether the particular host species used by butterflies generate more anthropogenic geographic opportunity than alternative host sets with the same richness and plant-family composition. This distinction turns a structural comparison into an ecological one: a butterfly can benefit strongly from plant redistribution only if its documented hosts are themselves among the plants redistributed into new regions.

We ask four linked questions. First, how much does retaining introduced host ranges enlarge butterfly resource geography, and does proportional expansion depend on taxonomic host breadth? Second, after fixing plant-family identity and the number of host species, do observed butterfly host portfolios produce more expansion than random documented host plants from the same families? Third, is the introduced-added part of the reconstructed host envelope geographically relevant to contemporary butterflies—specifically, are modern butterfly records outside native host envelopes enriched within those introduced-added host regions? Fourth, after resource opportunity is defined, how strongly is realized butterfly geography filtered by climate, and is that filtering weaker in family-level generalists?

The large resource reconstruction is descriptive and the matched-host null is a post-hoc reviewer-defense analysis motivated by the structural bounds of the original concentration metrics. The occurrence validation uses a separately acquired butterfly panel that was not used to construct the 239-species resource result. The host-breadth/climate prediction was generated from an earlier pilot and then evaluated on an independent panel with pre-specified species selection, quality thresholds and inference.

We expected introduced host ranges to enlarge resource opportunity across the specialization spectrum. The matched-host null asks a sharper question without assuming its answer: whether observed host identity carries information beyond host count and plant-family composition. For the independent occurrence validation, we expected introduced-added host regions to recover more modern butterfly records than expected from their share of geographic headroom alone. The independent climate prediction remained directional: broader family-level diets were predicted to weaken climatic filtering within contemporary host-resource opportunity.


---

## 2. Methods

### 2.1 Study overview

The analysis combines five nested ecological dimensions:

1. **Taxonomic host breadth**: LepTraits number of larval host-plant families.
2. **Species-level host portfolio richness**: number of resolved larval host species.
3. **Geographic resource breadth**: the union of WGSRPD level-3 units occupied by resolved host plants.
4. **Anthropogenic host-identity effect**: deviation of observed resource expansion from same-plant-family, same-host-count random portfolios.
5. **Realized geographic validation and climate filtering**: enrichment of contemporary butterfly records in introduced-added host regions, followed by climatic mismatch within contemporary host opportunity.

Analysis scale followed a nested funnel. The descriptor set contained 339 species; 239 met resource-eligibility criteria and 215 met the conservative host-taxonomy lower-bound criterion used for the matched-host null. Host-contribution structure was retained as a secondary diagnostic in the 191 conservative species whose envelopes expanded. Geographic validation used a separately acquired 32-species butterfly panel; 24 of those passed the pre-climate quality gate and were climate-informative.


### 2.2 Butterfly trait panel

We used the butterfly reconstruction derived from LepTraits 1.0 (Shirey et al. 2022). The resource descriptor table contained 339 species with major trait axes, including host-family breadth, wing-size proxies, voltinism, and habitat-affinity traits. Resource analyses were restricted to species with positive host-family breadth and at least one reconstructed native host-resource WGSRPD3 unit, yielding 239 species.

Because host-interaction databases are incomplete, we also defined a conservative host-taxonomy lower-bound subset in which the number of resolved HOSTS-WCVP host species was at least as large as LepTraits host-family count. This condition is necessary but not sufficient for interaction completeness. It yielded 215 species.

### 2.3 Host interaction and plant distribution reconstruction

Larval host records came from a fixed snapshot of the HOSTS database of lepidopteran host plants (Robinson et al. 2023). Plant taxonomy and distribution came from WCVP v13 (Govaerts et al. 2021). Exact source identities and repository commits are provided in the reproducibility materials rather than repeated in the ecological Methods.

Species-level host names were resolved to accepted WCVP plant identifiers. The **native resource envelope** for each butterfly was the union of WGSRPD3 units in which any resolved host was recorded as native, extant, and non-doubtful. The **contemporary resource envelope** retained both native and introduced records while continuing to exclude extinct and location-doubtful records.

Geographic units followed level 3 of the World Geographical Scheme for Recording Plant Distributions (WGSRPD; Brummitt et al. 2001). WGSRPD3 is intentionally coarse. We therefore interpret these envelopes as regional resource opportunity, not as local host occupancy or confirmed butterfly habitat.

### 2.4 Taxonomic versus geographic specialization

For each of the 239 resource-eligible butterflies, we compared LepTraits host-family count with native host-resource WGSRPD3 breadth. We calculated Spearman rank correlations for the full resource-eligible panel and the 215-species host-taxonomy lower-bound subset.

To quantify discordance between taxonomic and geographic specialization, both variables were converted to percentile ranks. Large discordance was defined a priori in the reconstruction as an absolute rank difference >= 0.5. We also identified conservative "specialist-wide" species: one-family butterflies whose geographic resource breadth was in the upper quartile of the host-taxonomy-adequate panel.

### 2.5 Anthropogenic expansion of resource opportunity

For each butterfly we calculated:

- native resource units;
- contemporary resource units;
- introduced-added units = contemporary - native;
- contemporary/native ratio;
- log resource expansion = log(1 + contemporary) - log(1 + native);
- introduced share of contemporary geography.

We summarized these metrics across four host-family strata: one family, two families, three to five families, and six or more families. We tested descriptive rank associations between host-family breadth and proportional resource expansion. We repeated the association in the host-taxonomy lower-bound subset.

These analyses quantify potential geographic resource opportunity generated by host redistribution. They do not show that butterflies currently use every introduced host population or that butterfly ranges expanded because hosts were introduced.

### 2.6 Same-family, same-host-count host-identity null

The original host-contribution metrics are mathematically bounded by portfolio size, so we added a matched randomization to isolate ecological host identity from those structural constraints. For each resource-eligible butterfly, we preserved the exact number of resolved accepted host species contributed by each plant family. Within each family, we then sampled without replacement from all WCVP-resolved plant species recorded anywhere in the fixed HOSTS database. Thus every randomized portfolio had the same butterfly identity, plant-family composition and family-specific host counts as the observed portfolio, but different host species.

For each randomized portfolio we rebuilt native and contemporary WGSRPD3 resource unions and recalculated introduced-added units, log proportional expansion and host-contribution metrics. We used 999 deterministic draws per butterfly. Primary reviewer-defense quantities were the observed-minus-null-median residual in log expansion and the fraction of species exceeding their null median or 97.5th percentile. The candidate pool contains documented HOSTS plants rather than all vascular plants; the null therefore asks whether butterfly hosts are unusual relative to alternative documented host plants in the same families, not whether host choice is random in nature.

### 2.7 Host-contribution architecture

For butterflies whose resource envelopes expanded, we decomposed introduced-added WGSRPD3 units among contributing host species.

We summarized:

- number of contributing host species;
- **maximum single-host fractional share**, the fraction of added units attributable to the most dominant contributing host;
- **top-two-host fractional share**;
- **effective contributor number**, an inverse-concentration measure of how broadly added resource opportunity was distributed among host species;
- fraction of added resource units supported by multiple contributing host species.

We examined associations of family-level host breadth and resolved host-species richness with these architecture metrics, focusing on the conservative host-taxonomy-adequate subset.

### 2.8 Hierarchical specialization within host-family strata

To test whether species-level portfolio structure remains informative after holding family-level breadth approximately fixed, we repeated associations within each host-family stratum.

The most interpretable comparison used 82 expanded species that all had exactly one host family. Within this fixed family-breadth category, resolved host-species richness ranged from 1 to 37. We related host-species richness to effective contributor number, maximum single-host share, and total introduced-added units. For visualization, one-family species were also divided into low, middle, and high tertiles of resolved host-species richness.

This follow-up is exploratory and was motivated by the preceding host-contribution result.

### 2.9 Independent climate-filtering panel

An exploratory 10-species pilot was used only to generate the later host-breadth/climate prediction. Those ten species were excluded from the independent panel.

Before independent occurrence or climate responses were examined, we specified a 32-species panel balanced across four host-family strata (8 species each: 1, 2, 3-5, and 6+ families). Within strata, species were selected deterministically to span native resource breadth, resolved host-species richness, and wing size.

The pre-specified hypothesis was:

> Within larval host-resource opportunity, climatic filtering of realized butterfly distribution weakens as taxonomic host breadth increases.

### 2.10 Butterfly occurrence acquisition, quality gates and geographic validation

Butterfly occurrences were acquired through the GBIF Occurrence API (GBIF.org 2026) for 2010-2026, requiring coordinates, no flagged geospatial issue and occurrence status PRESENT. Each species was queried with the same deterministic sampling rule. Acquisition was initially incomplete for 14 species; the identical query rules were reapplied to those species without replacement of taxa or changes to ecological thresholds. After completion, 31/32 species had complete occurrence acquisition; one species remained rejected at taxon resolution.

A species passed the pre-climate quality gate only when all of the following held:

- occurrence transport was complete;
- resolved host-species count was at least host-family count;
- no more than 10% of observed WGSRPD3 units lay outside the contemporary host envelope;
- at least five contemporary host units had >=10 GBIF records from other independent-panel butterflies;
- among those effort-supported units, at least two contained the target butterfly and at least two did not.

Twenty-four species passed after uniform technical completion, exceeding the pre-specified minimum of 12.

We additionally used the full 32-species panel as an independent geographic validation of the reconstructed resource expansion. For every species, we counted butterfly-observed WGSRPD3 units outside its native host envelope and asked how many were recovered only when introduced host distributions were admitted. To diagnose whether overlap was expected simply because introduced host ranges occupy geographic space, we conditioned on each species' number of observations outside the native envelope. Under the diagnostic null, those observations were distributed uniformly without replacement over all WGSRPD3 units outside the native envelope; overlap with introduced-added host units is therefore hypergeometric. We convolved the species-specific distributions exactly to obtain a combined upper-tail probability. We repeated the summary for the 24 quality-qualified species. This is a geographic enrichment test, not a dispersal model, and an occurrence inside an introduced-added envelope does not by itself confirm local use of the introduced host.

### 2.11 Climate cross-fitting and within-region sensitivity

Climate was represented by CHELSA v2.1 1981-2010 variables BIO1, BIO7, BIO12, and BIO15 (Karger et al. 2017, 2021).

Within each quality-qualified species, observed WGSRPD3 units inside the contemporary host-resource envelope were deterministically split 50:50 into training and evaluation sets using a SHA256 species-by-unit rule.

The species climate center and scale were estimated only from occurrence-level climate values in training-observed units. Climate mismatch for each WGSRPD3 resource unit was calculated as the root mean squared standardized deviation from that training climate.

Evaluation compared effort-supported held-out observed units with effort-supported contemporary host units in which the butterfly was never observed. The species-level **climate-filtering score** was

`P(mismatch_never-observed > mismatch_held-out-observed) + 0.5 * P(tie)`.

A value of 0.5 is neutral; larger values indicate that butterfly-unobserved parts of the reconstructed resource envelope are climatically farther from the species' training niche than held-out observed portions.

Because climate mismatch could covary with broad geographic separation, we repeated the pairwise comparison using only held-out observed and never-observed resource units belonging to the same WGSRPD level-1 region. This coarse within-region restriction reduces intercontinental comparisons while preserving the original training niche and effort support. It is an accessibility sensitivity rather than a mechanistic dispersal model.

Species required at least 30 training occurrence records with climate, at least two effort-supported held-out observed units, and at least two effort-supported never-observed units. All 24 pre-climate-qualified species were climate-informative.

### 2.12 Independent test of host-breadth climate release

The primary response was the species climate-filtering score. The predictor was LepTraits host-family count, and contemporary host-resource WGSRPD3 breadth was the pre-specified control.

The observed statistic was partial Spearman correlation on ranks. The null distribution used 9,999 deterministic, species-identity-invariant Freedman-Lane-style reduced-response residual permutations. The one-sided alternative predicted a negative association: broader host-family breadth should reduce climate filtering. Alpha was 0.05. Because only 24 species were climate-informative, we also report an approximate Fisher-z 95% interval for the observed partial correlation and an approximate effect magnitude required for 80% power at one-sided alpha = 0.05; these are precision diagnostics rather than replacements for the permutation test.

We used a row-order-invariant reduced-response permutation that preserved predictor-control structure. The inferential specification was finalized before the independent ecological responses were examined; implementation history is documented in the reproducibility supplement.

---

## 3. Results

### 3.1 Introduced host ranges expand resource opportunity across the specialization spectrum

When introduced host distributions were retained, 206 of 239 species (86.2%) gained reconstructed geographic resource opportunity. Total species × region units increased from 26,530 under native-only host distributions to 41,083 under contemporary distributions, an addition of 14,553 units (+54.9%).

Proportional resource expansion was almost unrelated to host-family breadth (Spearman rho = 0.008). The conservative 215-species subset gave the same result (rho = 0.015), with 191 species (88.8%) expanding. Median contemporary/native ratios were 1.390, 1.442, 1.411 and 1.361 for the one-, two-, three-to-five- and six-plus-family strata, respectively.

The near-zero proportional association was not created by the finite number of WGSRPD3 units. After controlling native resource breadth, partial Spearman rho was 0.012; after additionally removing butterfly-family mean differences it was 0.022. Absolute added units showed only a weak positive association with family breadth (rho = 0.127). Thus broad diets had larger resource geographies on average but did not receive a clear proportional advantage from host redistribution.

### 3.2 Observed host identity produces more expansion than matched alternative host sets

The same-family, same-host-count null showed that host identity matters beyond host richness and plant-family composition. In the conservative 215-species subset, observed log expansion exceeded the species-specific null median by a median 0.140 log units, equivalent to a 15.1% higher expansion ratio on the exponentiated scale. Overall, 156/215 species (72.6%) exceeded their null median, 33 exceeded the null 97.5th percentile, and only two fell below the 2.5th percentile.

The effect occurred throughout the specialization spectrum. Among 97 one-family butterflies, the median observed-minus-null log expansion was 0.111 (11.7% on the exponentiated scale), with 63/97 above the null median and nine above the 97.5th percentile. Corresponding median residuals were 0.143 for two-family species, 0.134 for three-to-five-family species and 0.146 for six-plus-family species.

Crucially, this null-adjusted excess expansion was itself almost unrelated to host-family breadth (Spearman rho = 0.041). Human redistribution therefore did not preferentially amplify broad family-level generalists. Instead, butterflies across the specialization spectrum disproportionately used host identities whose introduced distributions generated more geographic opportunity than alternative documented plants from the same families.

### 3.3 Introduced-added host geography is enriched for contemporary butterfly occurrences

The independent occurrence panel provided a direct geographic validation of this reconstructed opportunity. Across 32 species, 115 butterfly species × WGSRPD3 observations lay outside native host envelopes. Sixty-six of these (57.4%) fell inside regions added only by introduced host distributions. Given the amount of non-native geographic headroom occupied by those introduced host ranges, the random-overlap expectation was 38.0 units; the observed 66 was strongly enriched (exact combined hypergeometric upper-tail p = 1.1 × 10^-11).

The result strengthened in the 24 species that passed the independent quality gate. Of 61 observed units outside native host envelopes, 49 (80.3%) were recovered by introduced-added host regions, compared with 30.6 expected under the same geographic null (p = 5.5 × 10^-9). Only 12 of 929 observed species × region units in this qualified panel lay outside the contemporary host envelope.

These results do not demonstrate which host was used at a locality, but they show that the resource geography created by introduced hosts coincides with contemporary butterfly geography far more often than expected from its geographic extent alone.

### 3.4 Taxonomic specialization only partly tracks geographic resource specialization

Across 239 resource-eligible species, host-family breadth and native geographic resource breadth were positively but weakly associated (Spearman rho = 0.276). In the 215-species host-taxonomy lower-bound subset, the association increased to rho = 0.398 but remained far from one-to-one.

Within the conservative subset, 33 species (15.3%) differed by at least 0.5 between their percentile ranks for taxonomic and geographic specialization. Eighteen one-family specialists nevertheless occupied the upper quartile of geographic resource breadth. Thus low family-level diet breadth did not imply a geographically narrow larval-resource base.

### 3.5 Portfolio concentration gradients are mostly structural, but observed portfolios are non-random

The original concentration metrics showed strong raw gradients: in 191 conservative expanded species, host-family breadth correlated positively with effective contributor number (rho = 0.492) and negatively with maximum single-host share (rho = -0.486). However, family breadth also strongly predicted the number of contributing host species (rho = 0.501), and contributor count nearly determined the concentration metrics. After controlling contributor count, the family-breadth associations fell to rho = 0.071 for effective contributor number and -0.093 for maximum single-host share. The corresponding butterfly-family-blocked values were 0.132 and -0.159.

We therefore do not interpret the raw specialist-generalist concentration gradient as an independent mechanism. The matched host-identity null nevertheless revealed non-random portfolio composition: among the 191 conservative expanded species, effective contributor number exceeded the matched null 97.5th percentile in 46 species and fell below it in only one; maximum single-host share fell below its null 2.5th percentile in 26 species and exceeded the 97.5th percentile in only one. Observed butterfly host sets therefore tended to distribute anthropogenic opportunity across more contributing hosts than alternative same-family, same-count host sets, even though the between-category concentration gradient itself was largely structural.

### 3.6 Climate-associated filtering persists within broad geographic regions

After uniform technical completion of occurrence acquisition, 24 independent-panel species passed the pre-climate gate, and all 24 were climate-informative. The median climate-filtering score was 0.801; 23/24 species exceeded the neutral value of 0.5.

Restricting comparisons to held-out observed and never-observed resource units within the same WGSRPD level-1 region left 23 informative species. The median filtering score remained 0.798, with 21/23 species above 0.5. Thus the strong climate-associated separation was not solely a consequence of comparing host opportunity on different continents.

### 3.7 Broader host-family diets did not detectably weaken climate filtering

The prospectively specified climate-release prediction was not supported. After controlling contemporary geographic resource breadth, the partial Spearman association between host-family breadth and climate-filtering score was -0.166 (one-sided residual-permutation p = 0.2237; n = 24).

The approximate 95% interval for this partial correlation was broad (-0.54 to +0.26), and an absolute partial correlation of roughly 0.51 would have been required for 80% power under a one-sided alpha = 0.05 approximation. The independent test therefore rules out neither moderate negative effects nor moderate positive effects.

The within-WGSRPD-level-1 sensitivity gave no indication that geographic restriction rescued the predicted direction: partial rho = +0.207, one-sided p = 0.8232 (n = 23). The data support widespread climate-associated filtering but not the claim that broader family-level diets consistently weaken it.


---

## 4. Discussion

Human redistribution of plants has altered butterfly resource geography in a way that cannot be reduced to host richness. Three results form the core of the study. First, retaining introduced host distributions increased reconstructed larval-resource geography for most butterflies. Second, the actual host species used by butterflies generated more expansion than alternative host sets matched for plant-family composition and host count. Third, the introduced-added portion of those host envelopes was strongly enriched for independent contemporary butterfly occurrences outside native host geography. Together, these results move the inference beyond the arithmetic statement that more hosts create more possible regions: **host identity determines which butterflies inherit geographic opportunity from plant redistribution.**

### 4.1 Host identity, not generalism alone, links plant redistribution to butterfly opportunity

The matched null changes the interpretation of the resource result. If the apparent expansion were only a consequence of drawing more host species, or of using plant families that happen to be widely distributed, observed portfolios should resemble random HOSTS plants after family-specific host counts are fixed. They did not. Nearly three quarters of the conservative species exceeded their own null median, and excess expansion occurred even among one-family butterflies.

At the same time, null-adjusted expansion remained almost unrelated to host-family breadth. This combination is informative. Anthropogenic host redistribution is not a uniform bonus that simply increases with taxonomic generalism. Instead, butterflies across the specialist-generalist spectrum differ in whether their particular hosts are among the plants that have acquired introduced distributions. A family-level specialist can therefore receive substantial anthropogenic opportunity when one or more of its accepted hosts have been widely redistributed.

This interpretation is consistent with regional observations that exotic plants can become larval resources for butterflies, but extends the problem to global resource geography. The global signal lies in the non-random identities of redistributed hosts, not merely in counting host taxa.

### 4.2 The reconstructed opportunity is geographically relevant to modern butterflies

Potential resource maps are easy to overinterpret, so the independent occurrence validation is central. More than half of all butterfly observations lying outside native host envelopes were recovered when introduced host ranges were admitted, and the overlap greatly exceeded a null based on the geographic extent of those added host regions. In the quality-qualified panel, four fifths of outside-native observations were recovered.

This does not prove local feeding on an introduced plant. A butterfly occurrence and an introduced host can coexist in the same WGSRPD3 unit without interacting, and both may track a third environmental gradient. Nevertheless, the enrichment rejects a weaker criticism: that introduced host ranges merely add large amounts of irrelevant geographic area. The added resource geography disproportionately coincides with where contemporary butterflies are actually recorded.

The residual observations outside contemporary host envelopes are also informative. They can arise from incomplete HOSTS records, incomplete plant-distribution data, mapping error, vagrancy or host use absent from the frozen reconstruction. They should not be interpreted as host-independent occupancy.

### 4.3 Taxonomic specialization and spatial resource specialization are distinct

Family-level host breadth only partly tracked geographic resource breadth. Some one-family butterflies possessed very broad native resource envelopes because their hosts themselves were widespread. Conversely, adding host families did not guarantee proportionally larger gains from contemporary plant redistribution.

This distinction matters because "specialist" is often used as if it simultaneously describes phylogenetic diet breadth, number of host species, spatial resource extent and sensitivity to environmental filtering. Our results show that these dimensions can separate. For macroecological prediction, the identity and geography of host species can be at least as important as a single taxonomic breadth count.

### 4.4 Portfolio concentration should not be mistaken for an independent breadth effect

The original analysis suggested that specialists gain through concentrated host contributions whereas generalists gain through more distributed portfolios. The raw pattern is real descriptively, but the reviewer-defense analyses show why it is not sufficient as a central mechanism: effective contributor number cannot exceed contributor count, and the minimum possible maximum share declines as contributors increase. Once contributor count was held constant, most of the association with family breadth disappeared.

This is a useful negative correction rather than a loss of the main ecological signal. The matched host-identity null addresses the more relevant question. Observed host sets often involved more effective contributors and lower single-host dominance than alternative same-family, same-count portfolios, indicating non-random host composition. What cannot be defended is the stronger claim that family-level generalism itself independently causes portfolio dispersion.

### 4.5 Climate filters opportunity after hosts define it

Contemporary host geography was not sufficient to reproduce realized butterfly geography. Never-observed resource units were generally more climatically mismatched than held-out observed units, and this pattern persisted when comparisons were restricted within the same broad WGSRPD level-1 regions.

The directional prediction that family-level generalists would show weaker climatic filtering was not supported. The small independent panel also provides limited precision, so the result should not be read as evidence that host breadth has exactly zero climatic effect. Rather, it shows that a simple monotonic release-from-climate model is not sufficient at the effect sizes this design could resolve.

A more useful hierarchy is therefore: host identity and plant redistribution determine a changing geographic opportunity set; climate, dispersal, habitat and other processes filter which parts of that opportunity are realized. Different axes of specialization can enter at different levels of that hierarchy.

### 4.6 Implications for global change biogeography

Biogeographic responses to global change are often framed around the movement of focal species. Our results highlight a complementary route: human transport moves resources, thereby changing the opportunity landscape before or independently of consumer range change. For butterflies, introduced plants can create potential larval-resource regions that were absent from native host geography, and those regions are already enriched for contemporary butterfly records.

This mechanism creates asymmetric consequences across consumers because host identity matters. Two butterflies with the same number of host species can inherit different geographic opportunities if one uses hosts that humans have redistributed widely and the other does not. Conversely, a taxonomic specialist can acquire substantial new opportunity without becoming a generalist.

Future tests should connect these macrogeographic patterns to local interaction data, establishment dates and temporal range change. The strongest next prediction is not that generalists always benefit more, but that butterfly expansion into newly available regions should be greatest where the species' realized larval hosts have undergone unusually large anthropogenic redistribution and where climatic and habitat filters are weak.


---

## 5. Limitations

First, HOSTS is incomplete and geographically biased. The conservative lower-bound filter removes obvious cases in which resolved host-species richness is smaller than the independent family-breadth count, but it cannot recover undocumented interactions. The same-family null is therefore conditional on plants represented in HOSTS and should not be interpreted as randomization over all physiologically possible hosts.

Second, WCVP introduced status describes distribution status rather than the date, pathway or local abundance of establishment. We reconstruct anthropogenic geographic resource opportunity, not a historical causal effect of plant introduction on butterfly range expansion.

Third, WGSRPD3 is coarse. Co-occurrence of a butterfly and an introduced host in a region is not proof of local interaction. The occurrence enrichment therefore validates geographic relevance of the added envelope, not realized larval host use at individual records.

Fourth, GBIF non-observation is not confirmed absence. The climate analysis used an independent-panel effort proxy and held-out occurrences, but it remains a macroecological filtering analysis rather than an occupancy model from standardized surveys.

Fifth, broad geographic accessibility is imperfectly represented. Restricting the climate comparison within WGSRPD level-1 regions preserved the filtering result, but level-1 regions are much coarser than species-specific dispersal domains. Explicit accessible-area models would provide a stronger test.

Sixth, species are phylogenetically non-independent. Removing butterfly-family mean differences produced nearly identical resource-expansion conclusions and did not restore the raw architecture effects, but family blocking is not a substitute for a resolved species-level phylogenetic covariance model. A PGLS or phylogenetic mixed model should be added if a defensible dated phylogeny covering the panel can be assembled without substantial taxonomic loss.

Seventh, the matched-host null, structural architecture diagnostic and occurrence-enrichment null were added after the original v0.1 manuscript in response to reviewer-style concerns. They are explicitly post-hoc robustness analyses. Their definitions were fixed before their corresponding outputs were opened, but they should not be described as prospectively confirmatory.

Finally, the independent climate panel is small. The approximate 95% interval around the host-breadth partial correlation is wide, and the design had useful power only for relatively large effects. The unsupported directional prediction is therefore a boundary on the present evidence, not evidence of exact absence.


---

## 6. Conclusions

Human redistribution of larval host plants has expanded reconstructed butterfly resource geography on a global scale, but the key signal is not that generalists mechanically have more hosts. After fixing plant-family composition and host-species counts, observed butterfly host portfolios still generated more anthropogenic expansion than alternative documented host sets, including among one-family butterflies.

That reconstructed opportunity is geographically consequential: introduced-added host regions recovered contemporary butterfly observations outside native host envelopes far more often than expected from their geographic extent alone. Host identity therefore links plant redistribution to consumer opportunity.

Taxonomic host breadth remains only one layer of specialization. It weakly predicts proportional anthropogenic gain after the matched null, the original concentration gradient is largely structural, and broader family-level diets did not detectably weaken climate-associated filtering. The emerging hierarchy is: **host identity determines which redistributed plants can create opportunity; plant biogeography determines where that opportunity occurs; and climate-associated and other filters determine how much is realized.**


---

## Figure legends

**Figure 1. Introduced host distributions expand reconstructed resource opportunity across the specialization spectrum.** **a**, Native versus contemporary host-resource breadth for 239 resource-eligible butterflies. The dashed line is the 1:1 expectation; points above it gain WGSRPD3 resource units when introduced host distributions are retained. Overall, 206/239 species expanded, and aggregate species × region units increased from 26,530 to 41,083 (+14,553; +54.9%). **b**, Log proportional resource expansion by host-family breadth class. Host-family breadth was nearly unrelated to proportional expansion (Spearman rho = 0.008), including after controlling native resource breadth (partial rho = 0.012).

**Figure 2. Observed host identities create more anthropogenic expansion than matched alternative host sets, and the added geography is realized disproportionately often.** **a**, Observed-minus-null-median log resource expansion for the 215 conservative species under 999 same-plant-family, same-host-count random portfolios per species; zero denotes the matched null median. **b**, Butterfly species × region observations outside native host envelopes that were recovered by introduced-added host ranges in the full 32-species independent panel and the 24-species quality-qualified subset, compared with exact random geographic-overlap expectations.

**Figure 3. Taxonomic host breadth only partly predicts geographic resource breadth.** Host-family count is plotted against native host-resource WGSRPD3 breadth for the 215 species meeting the conservative host-taxonomy lower-bound criterion. One-family butterflies can occupy the upper quartile of geographic resource breadth, illustrating that taxonomic and spatial specialization are not equivalent.

**Figure 4. Raw host-contribution concentration gradients are structurally bounded by contributor richness.** For 191 conservative expanded species, host-family breadth is related to the number of contributing hosts and therefore to effective contributor number and maximum single-host share. The figure contrasts raw architecture values with residual or conditional quantities after contributor count is accounted for; the raw specialist-generalist gradient is treated as descriptive rather than as an independent mechanism.

**Figure 5. Climate-associated filtering persists within broad regions, but broader diets do not show weaker filtering.** Climate-filtering scores are shown for the original 24 climate-informative species and the 23-species same-WGSRPD-level-1 sensitivity. The original primary test controlled contemporary resource breadth and did not support weaker climate filtering in broader host-family generalists (partial Spearman rho = -0.166; one-sided p = 0.2237). The within-region sensitivity retained strong filtering (median = 0.798; 21/23 >0.5) and did not restore the predicted negative association.


## Data and Code Availability

All analysis code and manuscript-facing result summaries are versioned in the public `chocho` repository. The reconstruction uses fixed public releases or repository snapshots of LepTraits, HOSTS, WCVP, WGSRPD, GBIF and CHELSA. Exact source identities, hashes, workflow records and historical analysis contracts are documented separately in the repository's reproducibility and provenance materials so that the main Methods can remain focused on ecological design. A permanent archival snapshot and DOI will be supplied before submission.


## References (working)

- Brummitt, R. K., Pando, F., Hollis, S. & Brummitt, N. A. 2001. *World Geographical Scheme for Recording Plant Distributions*, 2nd edn. Hunt Institute for Botanical Documentation, Carnegie Mellon University.
- Chowdhury, S. et al. 2026. Extensive climate-induced range shifts in butterflies across the globe. *Nature Ecology & Evolution*. https://doi.org/10.1038/s41559-026-03117-y
- Couto, H. et al. 2026. The Lepidopteran Hitchhiker's Guide to the Globe: The Spread and Dispersal of Non-Native Moths and Butterflies. *Global Ecology and Biogeography*. https://doi.org/10.1111/geb.70292
- GBIF.org. 2026. GBIF Occurrence API, version 1. Global Biodiversity Information Facility. Occurrence queries accessed September 2026.
- Govaerts, R., Nic Lughadha, E., Black, N., Turner, R. & Paton, A. 2021. The World Checklist of Vascular Plants, a continuously updated resource for exploring global plant diversity. *Scientific Data* 8: 215. https://doi.org/10.1038/s41597-021-00997-6
- Graves, S. D. & Shapiro, A. M. 2003. Exotics as host plants of the California butterfly fauna. *Biological Conservation* 110: 413–433. https://doi.org/10.1016/S0006-3207(02)00233-1
- Gross, C., Kawahara, A. & Daru, B. 2026. Climate and regional plant richness drive diet specialization in butterfly caterpillars. *Nature Communications*. https://doi.org/10.1038/s41467-026-73236-4
- Guo, F., McKirdy, S. J., Gao, L. & Gao, G. 2026. Climate and traits are differentially associated with range extent and range geometry in global butterflies. *Ecological Indicators* 189: 115231. https://doi.org/10.1016/j.ecolind.2026.115231
- Karger, D. N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R. W., Zimmermann, N. E., Linder, H. P. & Kessler, M. 2017. Climatologies at high resolution for the earth's land surface areas. *Scientific Data* 4: 170122. https://doi.org/10.1038/sdata.2017.122
- Karger, D. N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R. W., Zimmermann, N. E., Linder, H. P. & Kessler, M. 2021. Climatologies at high resolution for the earth's land surface areas. EnviDat. https://doi.org/10.16904/envidat.228
- Lancaster, L. T. 2020. Host use diversification during range shifts shapes global variation in Lepidopteran dietary breadth. *Nature Ecology & Evolution* 4: 963–969. https://doi.org/10.1038/s41559-020-1199-1
- Rashid, S., Wessely, J., Hausharter, J., Moser, D., Gattringer, A., Fiedler, K., Hülber, K. & Dullinger, S. 2026. Food Plant Availability Constrains Climatic Niches of Host-Specialized Europe-Centred Butterflies. *Diversity and Distributions* 32: e70245. https://doi.org/10.1111/ddi.70245
- Robinson, G. S., Ackery, P. R., Kitching, I., Beccaloni, G. W. & Hernández, L. M. 2023. HOSTS - a Database of the World's Lepidopteran Hostplants [Data set]. Natural History Museum. https://doi.org/10.5519/havt50xw
- Shirey, V., Larsen, E., Doherty, A. et al. 2022. LepTraits 1.0: A globally comprehensive dataset of butterfly traits. *Scientific Data* 9: 382. https://doi.org/10.1038/s41597-022-01473-5
- Slove, J. & Janz, N. 2011. The relationship between diet breadth and geographic range size in the butterfly subfamily Nymphalinae: a study of global scale. *PLoS ONE* 6: e16057. https://doi.org/10.1371/journal.pone.0016057

## Repository provenance

Primary result receipts:

- `benchmarks/exploratory/butterfly_specialization_dimensionality_result_v0.1.json`
- `benchmarks/exploratory/butterfly_anthropogenic_resource_expansion_result_v0.1.json`
- `benchmarks/exploratory/butterfly_resource_expansion_mechanism_result_v0.1.json`
- `benchmarks/exploratory/butterfly_host_specialization_hierarchy_result_v0.1.json`
- `benchmarks/exploratory/butterfly_climate_release_postgate_independent_result_v0.1.json`
- `benchmarks/exploratory/butterfly_reviewer_defense_result_v0.1.json`

Integrated synthesis:

- `benchmarks/exploratory/butterfly_specialization_ecology_synthesis_v0.1.json`
- `docs/exploratory/BUTTERFLY_SPECIALIZATION_ECOLOGY_SYNTHESIS_V0_1.md`
