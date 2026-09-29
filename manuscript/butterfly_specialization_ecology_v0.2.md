# Human redistribution of host plants expands butterfly resource geography largely independently of diet breadth

**Working manuscript v0.2 — 2026-09-28**

## Abstract

**Aim:** To quantify how human redistribution of larval host plants changes butterfly resource geography, whether proportional gains depend on taxonomic diet breadth, and whether introduced-host geography aligns with contemporary butterfly occurrence.

**Location:** Global.

**Time period:** Contemporary host distributions; butterfly occurrences 2010–2026.

**Major taxa studied:** Butterflies and larval host plants.

**Methods:** We reconstructed native and contemporary host-resource geography for 239 butterflies. We decomposed added butterfly × region opportunity across contributing host plants using unit-preserving fractional attribution, with genus-level aggregation as a taxonomic sensitivity. We tested host-family breadth against proportional resource expansion using finite-area, regional, butterfly-Family and species-level phylogenetic sensitivities. A 32-species panel originally stratified for climate analysis was reused secondarily to test whether introduced hosts recover occurrences outside native host-resource envelopes beyond structural-overlap expectations.

**Results:** Introduced hosts expanded reconstructed resource geography for 206/239 species (86.2%) and increased aggregate species × region coverage by 54.9%. The added geography was unevenly generated across 670 contributing host species: the top 10 accounted for 25.1% of all added butterfly × region units, the top 50 for 57.0%, and 38 species together accounted for half. Host-family breadth was nearly unrelated to proportional expansion (Spearman rho = 0.008; bootstrap 95% CI -0.111 to 0.128); a strict post-hoc ±0.10 equivalence test was narrowly inconclusive (TOST p = 0.078). Introduced hosts recovered 66/115 outside-native species × region observations (57.4%). Mean species-level recovery was 0.583 versus a region-matched null median of 0.385, with no exceedance in 99,999 draws. Excluding *Pyrgus communis*, the largest contributor, still left 44/93 units recovered versus a null median of 30.

**Main conclusions:** Human redistribution of host plants has substantially expanded butterfly resource opportunity across the specialization spectrum, while taxonomic diet breadth provides little information about proportional gain. The expansion is nevertheless strongly uneven on the resource side, with a minority of host plants accounting for much of the added geography. Secondary occurrence evidence shows that contemporary introduced-host geography captures distributions missed by native-only resource maps without establishing local host use or causal range expansion.

**Keywords:** biotic redistribution, butterflies, host plants, introduced species, resource geography, specialization
---

## 1. Introduction

Specialization is one of the most familiar axes used to compare herbivorous insects. Butterflies are routinely described as specialists or generalists according to the number or phylogenetic breadth of larval host plants they use, and diet breadth has been linked to geographic range size, range dynamics, diversification and environmental gradients (Slove & Janz 2011; Lancaster 2020; Gross et al. 2026). Yet taxonomic diet breadth and the geography of those resources are not the same ecological quantity.

A butterfly feeding within one plant family may use one host species or dozens, and those hosts may themselves range from local to nearly cosmopolitan. Human transport adds a second source of decoupling because plant geography is no longer determined only by native biogeographic history. Regional studies have shown that introduced plants can become larval hosts and can accompany butterfly persistence or range change (Graves & Shapiro 2003), while recent global work documents extensive human-mediated redistribution of Lepidoptera themselves (Couto et al. 2026). Thus the same taxonomic diet breadth can correspond to very different amounts of geographic resource opportunity and very different anthropogenic gains.

A resource-side reconstruction is biologically useful only if the added geography bears some relationship to contemporary butterfly distributions. A native-only resource envelope may classify an observed butterfly region as host-unavailable even when a known host has been introduced there. Conversely, the presence of a known host does not guarantee butterfly occurrence because climate, dispersal, habitat, phenology and biotic interactions can filter otherwise available resource geography. These processes motivate a layered view in which plant distributions define resource opportunity and additional filters shape its realization.

Here we ask three primary questions. First, how strongly do introduced host distributions expand butterfly resource geography, and how concentrated is that added opportunity among host plants? Second, is proportional expansion greater in butterflies with broader family-level diets? Third, does contemporary introduced-host geography recover butterfly occurrences that fall outside native host-resource envelopes more often than structural overlap predicts? Matched-host and plant-prominence analyses are used as post-hoc diagnostics of possible host-identity and data-structure effects rather than as headline mechanisms. We retain a smaller climate analysis as a secondary test of whether available resource geography is fully realized.

The resource reconstruction is descriptive and uses fixed LepTraits, HOSTS and WCVP inputs. The 32-species panel was assembled separately for the climate test, using eligibility and analysis rules fixed before its ecological responses were evaluated; the occurrence analysis is a secondary use of that panel. All later null-model, ceiling, phylogenetic and geographic sensitivities are identified as post-hoc robustness analyses.

We expected anthropogenic host redistribution to enlarge resource opportunity across much of the specialization spectrum, but did not assume that broader taxonomic diets must receive a larger proportional gain. We also expected the contemporary resource envelope to recover at least some butterfly occurrences missed by a native-only envelope. 

## 2. Methods

### 2.1 Study design

We separated four layers of butterfly specialization and biogeography: taxonomic host breadth, geographic host-resource breadth, anthropogenic expansion of that resource geography, and realization of contemporary resource opportunity by butterflies.

The LepTraits-derived descriptor panel contained 339 species. Resource analyses used 239 species with positive host-family breadth and at least one reconstructed native host-resource WGSRPD3 unit. A conservative subset of 215 species also required the number of resolved host species to equal or exceed LepTraits host-family count. A separate 32-species panel was used for occurrence validation and climate analysis; it had been selected to balance the subsequent climate test, not optimized for validation. Twenty-four species passed the climate quality criteria.

### 2.2 Host-resource reconstruction

Larval host records came from HOSTS and plant taxonomy and distributions from WCVP (Robinson et al. 2023; Govaerts et al. 2021). Species-level host names were resolved to accepted WCVP taxa.

For each butterfly, the **native resource envelope** was the union of WGSRPD level-3 units in which any resolved host was native, extant and non-doubtful. The **contemporary resource envelope** additionally retained introduced host distributions. WGSRPD3 is a coarse regional framework, so these envelopes represent geographic resource opportunity rather than local host occupancy.

We calculated native and contemporary resource breadth, introduced-added units, the contemporary/native ratio, log resource expansion, and the introduced share of contemporary geography. In a post-hoc descriptive decomposition, each butterfly × WGSRPD3 unit newly added by introduced host distributions contributed one unit of credit divided equally among the host species that added that unit. Credits therefore summed exactly to the aggregate number of added butterfly × region units. We summarized the concentration of those credits across host plants using cumulative top-k shares, the Gini coefficient and the inverse-Herfindahl effective contributor number. As a taxonomic sensitivity, contributing plant species were first collapsed to unique genera within each butterfly × region unit and the unit credit was then divided equally among those genera. This decomposition is explicitly contribution-weighted: a plant can accumulate credit by being used by more focal butterflies, by adding more geographic units, or both. We therefore use it to describe which resources generate the reconstructed opportunity, not as a sample of independent plant-level effects.

### 2.3 Robustness analyses for host identity, specialization and taxonomic structure

We used matched-host nulls to test whether apparent host-identity effects could arise from host richness, plant-family composition, starting native range and the prominence of candidate plants in HOSTS. Because conclusions about host-identity excess were sensitive to how strongly network prominence was weighted, we do not interpret these nulls as identifying a butterfly-specific host-choice mechanism; the full parameter grid is reported in Supplementary Information.

We tested whether the host-family-breadth association with proportional expansion was sensitive to the finite 369-unit WGSRPD3 support, geographic composition and taxonomic non-independence. Post-hoc analyses progressively excluded broad native resource envelopes, repeated associations within dominant WGSRPD level-1 regions, fitted butterfly Family as a random intercept, and used the Kawahara et al. (2023) time-calibrated butterfly phylogeny for exact species matches. For the phylogenetic sensitivity we report both a fixed Brownian model (lambda = 1) and a Pagel-lambda model in which residual phylogenetic covariance is estimated from the data. Because the estimated lambda was close to zero, we treat the estimated-lambda model as the data-adaptive phylogenetic correction and the fixed Brownian model as a boundary sensitivity rather than privileging lambda = 1.

An exploratory plant-level analysis related HOSTS consumer degree to plant geographic expansion across 8,909 HOSTS-WCVP plants. Because HOSTS includes moths as well as butterflies, this is a Lepidoptera-wide network measure rather than a butterfly-specific trait. We treat it only as a possible explanatory axis because degree may reflect ecological host prominence, study intensity, and reverse causation if introduced plants acquire additional consumer records after redistribution.

Within-butterfly portfolio-concentration metrics were retained only as a secondary structural diagnostic because effective contributor number and maximum single-host share are bounded by host number. This diagnostic is distinct from the across-plant concentration of total added opportunity described above.

### 2.4 Secondary occurrence validation

Butterfly occurrences for 2010–2026 were obtained from GBIF for the separately assembled 32-species panel and mapped to WGSRPD3 units. We counted occurrence units outside each species' native host-resource envelope and asked how many entered the contemporary envelope after introduced host ranges were retained.

Two post-hoc nulls tested whether recovery exceeded geometric overlap expected from larger envelopes. The first preserved native-envelope size, the number of introduced-added units and the number of outside-native occurrence units, but placed added units uniformly among non-native WGSRPD3 units. The second additionally preserved the number of added units within each WGSRPD level-1 region. Aggregate null distributions used 199,999 hypergeometric draws. We then treated butterfly species as the replication unit: 99,999 region-matched null draws were used for the mean and median species-level recovery fractions, 50,000 species-cluster bootstrap resamples quantified uncertainty, and every informative species was removed in turn to assess leave-one-species-out sensitivity.

### 2.5 Climate filtering within contemporary resource opportunity

The climate analysis tested the pre-specified prediction that broader host-family diets weaken climatic filtering within contemporary host-resource opportunity. Species entered this analysis only when occurrence coverage, host-taxonomy consistency and within-envelope sampling met fixed quality criteria; 24 species qualified.

Climate was represented by CHELSA BIO1, BIO7, BIO12 and BIO15. Within each species, observed resource units were split deterministically into training and evaluation sets. Climate center and scale were estimated from training occurrences. Held-out observed resource units were compared with effort-supported contemporary host units in which the butterfly was not observed.

The species climate-filtering score was the probability that a never-observed resource unit had greater climatic mismatch from the training niche than a held-out observed unit, with ties receiving half weight. A score of 0.5 is neutral.

The primary test used partial Spearman correlation between host-family breadth and filtering score while controlling contemporary resource breadth, with 9,999 residual permutations and a one-sided negative alternative.

To assess geographic confounding, post-hoc sensitivities restricted comparisons to occupied WGSRPD level-1 regions and matched held-out observed and never-observed resource units by distance to the nearest training-observed region using 250-, 500- and 1,000-km calipers. A species bootstrap summarized effect-size precision.

## 3. Results

### 3.1 Human redistribution broadly expands butterfly resource geography

Introduced host distributions expanded reconstructed resource opportunity for 206 of 239 species (86.2%). Aggregate species × WGSRPD3 coverage increased from 26,530 under native-only host distributions to 41,083 under contemporary distributions, an addition of 14,553 units (+54.9%).

Those 14,553 added butterfly × region units were generated unevenly across host plants. Fractional attribution identified 670 contributing host species. The top 10 species accounted for 25.1% of total credit, the top 20 for 37.3% and the top 50 for 57.0%; only 38 host species were needed to reach half of all added opportunity. The plant-level Gini coefficient was 0.777 and the inverse-Herfindahl effective contributor number was 92.5 species. Thus expansion had a long tail of contributors but a strongly concentrated upper tier rather than being distributed uniformly across known hosts.

A simple same-family matched null suggested that actual host identities were unusually expansion-prone, but this interpretation weakened once candidate plants were matched on native range and weighted by their use by other Lepidoptera. Under proportional degree weighting, observed total expansion lay inside the combined null; weaker weighting retained an excess. Because this attribution depended on the prominence-weighting strength, we do not infer an independent butterfly host-choice mechanism. The full 3 × 3 parameter grid and prominence-calibrated null are reported in Supplementary Information.

An exploratory Lepidoptera-wide plant analysis suggested that network-prominent plants are more often anthropogenically redistributed, but consumer degree can reflect both ecological prominence and database visibility, and redistribution itself may generate additional host-use records. We therefore treat network prominence as a possible explanation for the matched-null behavior rather than as the paper's primary result.

### 3.2 Expansion shows little relationship to broad family-level diet breadth

Host-family breadth was nearly unrelated to proportional expansion in the full panel (Spearman rho = 0.008; 49,999-bootstrap 95% CI -0.111 to 0.128) and conservative subset (rho = 0.015). A deliberately strict post-hoc equivalence diagnostic using |rho| < 0.10 as the equivalence region was narrowly inconclusive: the bootstrap 90% interval was -0.092 to 0.108 and the approximate Fisher-z TOST gave p = 0.078. We therefore interpret the full-panel association as small and imprecisely bounded rather than as proof of an exact zero effect. The stricter host-bias nulls did not reveal a hidden positive specialization gradient.

The near-zero relationship persisted after accounting for finite map support and taxonomic structure. Controlling native resource breadth gave a rank association of 0.012. A butterfly-Family random-intercept model estimated a standardized rank effect of 0.021 (95% CI -0.111 to 0.154; p = 0.752). The Kawahara et al. (2023) tree contained exact matches for 124/239 panel species. The estimated-lambda model found very little residual phylogenetic covariance (Pagel lambda = 0.033) and a rank-PGLS effect of -0.051 (95% CI -0.234 to 0.131; p = 0.581). A fixed Brownian lambda = 1 sensitivity gave beta = 0.140 (95% CI -0.036 to 0.317; p = 0.119), leaving moderate positive effects unresolved in that selective tree-matched subset. Restricting analyses to progressively narrower native envelopes produced correlations of 0.028–0.067.

Regional stratification likewise showed that the global result was not solely driven by Northern America. Within dominant native-resource regions, rho was -0.005 in Northern America (n = 102), 0.013 in Africa (n = 38), 0.029 in Temperate Asia (n = 24), 0.051 in Tropical Asia (n = 9) and 0.089 in Southern America (n = 32). Europe showed a modest positive association (rho = 0.233; n = 34).

### 3.3 Introduced host geography aligns with separate-panel butterfly occurrence

Thirty-one of 32 species in the separately assembled climate-stratified panel yielded contemporary GBIF records. Across the panel, 115 butterfly species × WGSRPD3 observations occurred outside native host-resource envelopes. Adding introduced host distributions recovered 66 of these units (57.4%).

Twenty-three species had at least one outside-native occurrence, and 18 recovered at least one unit after introduced hosts were added. The median species-level recovery fraction was 0.60; eight species recovered all outside-native observed units. Seven of those eight complete-recovery cases involved only 1–4 outside-native units; the exception was *Pyrgus communis*, with 22/22 units recovered.

Recovery exceeded geometric expectation at both pooled and species levels. A uniform envelope-enlargement null had median recovery of 38 units (95% interval 30–46), while a null preserving each species' broad regional distribution of added host geography had median 49 (43–56); none of 199,999 draws reached the observed 66.

Treating species rather than species × region units as the replication level gave the same conclusion. Mean species-level recovery was 0.583 versus a region-matched null median of 0.385 (95% interval 0.308–0.470), with no exceedance in 99,999 draws. A species-cluster bootstrap placed the pooled recovery fraction between 0.349 and 0.774 (95% percentile interval).

*Pyrgus communis* contributed 22 of the 66 recovered units, but removing it still left 44/93 units recovered (47.3%) versus a region-matched null median of 30 (24–36), again with no exceedance in 99,999 draws. Removing *Pieris brassicae* left 62/111 (55.9%) versus a null median of 46 (40–52). Recovery was heterogeneous rather than universal: *Erynnis tristis* recovered 0/11 outside-native units and *Historis acheronta* recovered 5/13.

### 3.4 Secondary climate analysis

Climate-associated filtering remained evident in the 24 climate-informative species. Median filtering scores were 0.801 in the original analysis and 0.714–0.750 after 250–1,000-km distance matching; 21/24, 23/24 and 22/24 species remained above the neutral score of 0.5 at the three matching scales.

The pre-specified prediction that broader host-family diets weaken filtering was not supported (partial Spearman rho = -0.166, one-sided p = 0.2237) and was imprecisely estimated (bootstrap 95% interval -0.583 to +0.261). Full diagnostics are shown in Supplementary Figure S1.

## 4. Discussion

Human redistribution of larval host plants has altered butterfly resource geography at a scale that is visible across hundreds of species. Introduced host distributions increased aggregate reconstructed resource opportunity by 54.9% and expanded the envelope of more than 86% of resource-eligible butterflies. The proportional gain showed little relationship to family-level diet breadth, while the resource-side contribution was highly uneven: 38 of 670 contributing host species accounted for half of all added butterfly × region opportunity. A secondary occurrence analysis then showed that introduced hosts recovered 57.4% of contemporary species × region observations that fell outside native host-resource envelopes.

The central ecological result is therefore not a specialist–generalist contrast in portfolio concentration. It is a separation between **taxonomic diet breadth, resource geography and realized distribution** under global change.

### 4.1 Resource redistribution, not broad diet breadth, structures anthropogenic opportunity

The most stable result is that human redistribution of host plants has enlarged butterfly resource geography across the specialization spectrum. The full-panel host-breadth association was close to zero and remained weak under finite-area, regional and random-effect sensitivities, although the strict ±0.10 equivalence test was narrowly inconclusive and the fixed-Brownian tree subset retained substantial uncertainty. Taxonomic diet breadth therefore provides little information about proportional anthropogenic resource gain, without requiring the stronger claim that its effect is exactly zero.

The clearer asymmetry appears on the plant side. Hundreds of host species contributed some added opportunity, but a small upper tier contributed disproportionately: the top 10 plants supplied one quarter of all fractional credit and 38 species supplied one half. Human redistribution can therefore produce broad consumer-side opportunity through a comparatively concentrated set of highly redistributed resources. This concentration is descriptive rather than mechanistic, but it gives a concrete resource-side explanation for how specialists and generalists can both experience substantial geographic gains.

The matched-host analyses help delimit, rather than establish, a mechanism. Actual host portfolios looked unusually expansion-prone under a simple same-family null, but the excess changed substantially when candidate plants were weighted by their prominence in the wider HOSTS network. That sensitivity means the host-choice component cannot be cleanly separated from plant-level commonness, human association and database visibility.

Network prominence remains a plausible explanatory hypothesis, not a causal conclusion. In Supplementary analyses, plants recorded with more Lepidoptera consumers were more likely to have introduced-range expansion even after native-breadth and plant-family adjustment. However, HOSTS includes moths as well as butterflies, and the direction of causation is unresolved: common, cultivated or weedy plants may both attract study and be moved by humans, while introduced plants may also acquire new Lepidoptera associations after arrival. The main ecological conclusion therefore does not depend on this association.

### 4.2 Secondary occurrence validation shows that added resource geography is ecologically relevant

A resource-envelope analysis can otherwise remain purely potential. This secondary occurrence comparison provides a bridge to realized biogeography.

More than half of species × region observations lying outside native host-resource envelopes were brought inside the envelope by adding introduced host ranges. Crucially, the result persisted when butterfly species were treated as the replication unit and when the largest-contributing species was removed. The validation therefore cannot be reduced to the fact that contemporary envelopes are larger or to a single high-leverage butterfly. This does not prove that the recorded butterfly used the introduced host at that locality, nor that host introduction caused colonization, but it demonstrates that a native-only view of larval resources can miss contemporary geographic opportunities that coincide with butterfly presence.

That distinction matters for macroecological analyses that combine consumer distributions with resource distributions. When resources themselves have been redistributed, native resource maps can create apparent consumer–resource mismatches that are partly artifacts of treating present-day interaction opportunity as if plant geography were still native.

### 4.3 Resource opportunity is not realized geography

The occurrence validation shows that contemporary host geography matters, but it does not imply that all available resource geography is realized. In the smaller climate panel, climatically mismatched resource units were less often observed even after geographic matching. At the same time, broader host-family diets did not detectably weaken this filtering, and the effect estimate was too imprecise to exclude moderate associations.

We therefore treat climate as a secondary boundary on resource realization rather than a co-equal mechanism in this paper. Dispersal history, habitat, phenology, adult resources, biotic interactions and imperfect detection can also generate unrealized resource opportunity.

### 4.4 Implications for global change biogeography

Global change redistributes interacting species as well as climate. Consumer biogeography therefore cannot always be interpreted against static or native-only resource templates.

For butterflies, introduced plants can enlarge potential larval-resource geography across both specialists and generalists. Some of that added opportunity coincides with contemporary butterfly occurrence, while much remains unrealized and environmentally filtered. Similar logic should apply to other consumers whose resources are transported, cultivated, invaded or otherwise redistributed by humans.

A useful next step is to replace coarse opportunity envelopes with local realized host-use data and dated introductions. That would allow direct tests of whether introduced resources facilitate colonization, whether consumers switch among hosts after arrival, and whether resource redistribution changes persistence rather than merely geographic potential. 

## 5. Limitations

HOSTS is incomplete and geographically uneven, and the conservative filter cannot make host-interaction sampling globally uniform. Regional stratification reduces concern that the main breadth–expansion result is solely a North American artifact, but tropical strata remain small. Species-level PGLS was possible for only 124/239 butterflies directly represented in the Kawahara et al. (2023) tree. Tree-matched species had greater expansion than unmatched species, although host-family breadth was similar and the breadth–expansion correlation was near zero in both subsets; PGLS is therefore used only as a slope sensitivity, not to estimate panel-wide expansion magnitude.

WCVP introduced status identifies contemporary distribution status rather than the timing or pathway of introduction. We therefore infer changed **resource opportunity**, not a causal historical effect of host introduction on butterfly range expansion. The host-contribution decomposition is also not a plant-trait effect analysis: fractional credit combines a plant's use across focal butterflies with the number of geographic units its introduced distribution adds.

Occurrence validation is presence-side and opportunistic: the 32-species panel was constructed for a stratified climate test, not specifically for validating introduced-host opportunity. A butterfly observation in a region covered by an introduced known host does not prove larval use of that host population, and GBIF non-observation is not confirmed absence. The panel also includes species with documented non-native populations, including *Pieris brassicae* (Phillips et al. 2020); excluding that species did not change the occurrence-overlap conclusion.

WGSRPD3 is coarse and imposes a finite geographic ceiling. Sensitivity analyses indicate that this ceiling does not explain the near-zero specialization gradient, but finer spatial data would improve inference.

The climate test contains only 24 informative species. Its broad bootstrap interval means moderate effects remain plausible despite the unsupported directional prediction.

Finally, the matched-host nulls, occurrence-overlap nulls, equivalence diagnostic, host-contribution concentration analysis, ceiling, phylogenetic and spatial sensitivities are post-hoc analyses and are interpreted accordingly. Other-Lepidoptera consumer count is an imperfect joint proxy for ecological host prominence and HOSTS recording intensity. The fixed HOSTS reconstruction also does not identify whether a consumer record predates plant introduction or occurred within the plant's native range, so a pre-introduction/native-area network degree cannot be reconstructed here. Moreover, the 3 × 3 parameter grid showed that attenuation of the simple host-identity excess depends on the strength of degree weighting. We therefore use these nulls to identify prominence as a confounding axis, not to claim that it fully mediates or quantitatively explains the host-identity pattern.

## 6. Conclusions

Human redistribution of host plants has substantially expanded reconstructed butterfly resource geography: more than 86% of resource-eligible species gained opportunity and aggregate species × region coverage increased by 54.9%. Host-family breadth was a poor predictor of proportional gain, although a strict ±0.10 equivalence test was narrowly inconclusive. In contrast, the plant-side contribution was strongly uneven: 38 of 670 host species accounted for half of all added butterfly × region opportunity.

The secondary occurrence analysis shows that the added geography is not merely cartographic potential. Introduced host distributions recovered 57.4% of butterfly species × region observations that lay outside native host-resource envelopes, more than expected under species-level region-matched nulls and after removal of the largest contributor.

The defensible picture is therefore layered: **taxonomic diet breadth provides little information about proportional anthropogenic resource gain; a concentrated subset of redistributed host plants contributes much of the added opportunity; human redistribution changes where larval resources occur; and only part of that expanded opportunity is reflected in realized butterfly geography.** Plant network prominence may help explain which hosts are redistributed, but current data cannot separate ecological prominence, recording bias and reverse causation.

## References (working)

- Brummitt, R. K., Pando, F., Hollis, S. & Brummitt, N. A. 2001. *World Geographical Scheme for Recording Plant Distributions*, 2nd edn. Hunt Institute for Botanical Documentation, Carnegie Mellon University.
- Chowdhury, S. et al. 2026. Extensive climate-induced range shifts in butterflies across the globe. *Nature Ecology & Evolution*. https://doi.org/10.1038/s41559-026-03117-y
- Couto, H. et al. 2026. The Lepidopteran Hitchhiker's Guide to the Globe: The Spread and Dispersal of Non-Native Moths and Butterflies. *Global Ecology and Biogeography*. https://doi.org/10.1111/geb.70292
- GBIF.org. 2026. GBIF Occurrence API, version 1. Global Biodiversity Information Facility. Occurrence queries accessed September 2026.
- Govaerts, R., Nic Lughadha, E., Black, N., Turner, R. & Paton, A. 2021. The World Checklist of Vascular Plants, a continuously updated resource for exploring global plant diversity. *Scientific Data* 8: 215. https://doi.org/10.1038/s41597-021-00997-6
- Graves, S. D. & Shapiro, A. M. 2003. Exotics as host plants of the California butterfly fauna. *Biological Conservation* 110: 413–433. https://doi.org/10.1016/S0006-3207(02)00233-1
- Gross, C., Kawahara, A. & Daru, B. 2026. Climate and regional plant richness drive diet specialization in butterfly caterpillars. *Nature Communications*. https://doi.org/10.1038/s41467-026-73236-4
- Guo, F., McKirdy, S. J., Gao, L. & Gao, G. 2026. Climate and traits are differentially associated with range extent and range geometry in global butterflies. *Ecological Indicators* 189: 115231. https://doi.org/10.1016/j.ecolind.2026.115231
- Kawahara, A. Y., Storer, C., Carvalho, A. P. S. et al. 2023. A global phylogeny of butterflies reveals their evolutionary history, ancestral hosts and biogeographic origins. *Nature Ecology & Evolution* 7: 903–913. https://doi.org/10.1038/s41559-023-02041-9
- Karger, D. N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R. W., Zimmermann, N. E., Linder, H. P. & Kessler, M. 2017. Climatologies at high resolution for the earth's land surface areas. *Scientific Data* 4: 170122. https://doi.org/10.1038/sdata.2017.122
- Karger, D. N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R. W., Zimmermann, N. E., Linder, H. P. & Kessler, M. 2021. Climatologies at high resolution for the earth's land surface areas. EnviDat. https://doi.org/10.16904/envidat.228
- Lancaster, L. T. 2020. Host use diversification during range shifts shapes global variation in Lepidopteran dietary breadth. *Nature Ecology & Evolution* 4: 963–969. https://doi.org/10.1038/s41559-020-1199-1
- Phillips, C. B., Brown, K., Green, C., Toft, R., Walker, G. & Broome, K. 2020. Eradicating the large white butterfly from New Zealand eliminates a threat to endemic Brassicaceae. *PLoS ONE* 15: e0236791. https://doi.org/10.1371/journal.pone.0236791
- Rashid, S., Wessely, J., Hausharter, J., Moser, D., Gattringer, A., Fiedler, K., Hülber, K. & Dullinger, S. 2026. Food Plant Availability Constrains Climatic Niches of Host-Specialized Europe-Centred Butterflies. *Diversity and Distributions* 32: e70245. https://doi.org/10.1111/ddi.70245
- Robinson, G. S., Ackery, P. R., Kitching, I., Beccaloni, G. W. & Hernández, L. M. 2023. HOSTS - a Database of the World's Lepidopteran Hostplants [Data set]. Natural History Museum. https://doi.org/10.5519/havt50xw
- Shirey, V., Larsen, E., Doherty, A. et al. 2022. LepTraits 1.0: A globally comprehensive dataset of butterfly traits. *Scientific Data* 9: 382. https://doi.org/10.1038/s41597-022-01473-5
- Slove, J. & Janz, N. 2011. The relationship between diet breadth and geographic range size in the butterfly subfamily Nymphalinae: a study of global scale. *PLoS ONE* 6: e16057. https://doi.org/10.1371/journal.pone.0016057

## Data and Code Availability

Analysis code and the inputs required to reproduce the reported results and figures are versioned in the study repository. Public source datasets include LepTraits, HOSTS, WCVP, GBIF, CHELSA and WGSRPD. Exact source identities, retrieval rules, robustness-analysis specifications and execution records are documented in the supplementary reproducibility materials rather than the main text. A permanent archival snapshot and DOI will be supplied with the final submission.

## Figure legends

**Figure 1. Human redistribution of host plants broadly expands butterfly resource geography with little proportional gradient across family-level diet breadth.** **a**, Native versus contemporary host-resource breadth across 239 butterflies; 206 species expanded and aggregate coverage increased by 54.9%. **b**, Log proportional resource expansion across four host-family-breadth classes. The full-panel association was nearly zero (Spearman rho = 0.008). The matched-host and plant-prominence diagnostics are reported in Supplementary Information.

**Figure 2. Introduced host geography recovers butterfly occurrences beyond structural overlap expectations.** **a**, Fraction of outside-native occurrence units recovered for each of 23 informative species. **b**, Observed aggregate and species-level recovery versus structural nulls, with leave-one-species-out sensitivity. Removing *Pyrgus communis* left 44/93 units recovered versus a region-matched null median of 30; no region-matched null draw equalled or exceeded the observed value in 99,999 replicates.

**Figure 3. Diet-breadth associations remain small across finite-support and regional sensitivities.** **a**, Host-family breadth versus proportional expansion after progressively excluding broad native resource envelopes. **b**, Within-region associations for species grouped by dominant native host-resource region.

**Supplementary Figure S1. Climate-associated filtering persists after geographic controls, whereas the predicted host-breadth release is unsupported.** **a**, Filtering scores under the original analysis, within-region restriction and three distance-matching calipers. **b**, Primary partial host-family-breadth effect with bootstrap 95% interval and approximate 80%-power detectable-effect threshold.
