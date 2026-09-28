# Anthropogenic host redistribution expands butterfly resource geography beyond matched-host expectations

**Working manuscript v0.2 — 2026-09-28**

## Abstract

**Aim:** To test how human redistribution of larval host plants changes butterfly resource geography, whether the added geography aligns with contemporary butterfly occurrences, and whether climate-associated filtering persists within that opportunity.

**Location:** Global.

**Time period:** Contemporary distributions; butterfly occurrences 2010–2026; CHELSA climatology 1981–2010.

**Major taxa studied:** Butterflies and larval host plants.

**Methods:** We reconstructed native and contemporary host-resource geography for 239 butterflies. For 215 conservative-panel species, matched null portfolios preserved exact host-species richness and plant-family composition while randomizing host identity. A separate 32-species occurrence panel tested whether introduced hosts recover butterfly observations outside native host envelopes beyond structural overlap expectations. Twenty-four quality-qualified species entered an independently specified climate cross-fit.

**Results:** Introduced hosts expanded resource geography for 206/239 species and increased aggregate species × region coverage by 54.9%. Actual host identities generated 13,529 added units versus a matched-null median of 8,046; mean and median proportional expansion also exceeded null expectations (all p = 0.0005). This excess was not concentrated in broad family-level generalists (observed rho = 0.015; matched-null p = 0.7865). Introduced hosts recovered 66/115 independent butterfly observations outside native host envelopes, exceeding a region-matched overlap null (median 49; p = 5 × 10^-6). Climate-associated filtering persisted after geographic matching, whereas the predicted weakening of filtering with broader diets was unsupported.

**Main conclusions:** Human redistribution of particular host species has created butterfly resource opportunity beyond that expected from host richness or plant-family breadth alone. Taxonomic diet breadth, resource geography and realized butterfly distributions are distinct dimensions of global-change response.

**Keywords:** butterflies; biotic redistribution; climate filtering; host plants; introduced species; resource geography

---

## 1. Introduction

Specialization is one of the most familiar axes used to compare herbivorous insects. Butterflies are routinely described as specialists or generalists according to the number or phylogenetic breadth of larval host plants they use, and diet breadth has been linked to geographic range size, range dynamics, diversification and environmental gradients (Slove & Janz 2011; Lancaster 2020; Gross et al. 2026). Yet taxonomic diet breadth and the geography of those resources are not the same ecological quantity.

A butterfly feeding within one plant family may use one host species or dozens, and those hosts may themselves range from local to nearly cosmopolitan. Human transport adds a second source of decoupling because plant geography is no longer determined only by native biogeographic history. Regional studies have shown that introduced plants can become larval hosts and can accompany butterfly persistence or range change (Graves & Shapiro 2003), while recent global work documents extensive human-mediated redistribution of Lepidoptera themselves (Couto et al. 2026). What remains less clear at global scale is how strongly host redistribution has changed the geographic template of potential butterfly resources.

A resource-side reconstruction is biologically useful only if the added geography bears some relationship to contemporary butterfly distributions. A native-only resource envelope may classify an observed butterfly region as host-unavailable even when a known host has been introduced there. Conversely, the presence of a known host does not guarantee butterfly occurrence because climate, dispersal, habitat, phenology and biotic interactions can filter otherwise available resource geography. These processes motivate a layered view in which plant distributions define resource opportunity and additional filters shape its realization.

Here we ask four linked questions. First, how strongly do introduced host distributions expand butterfly resource geography, and is proportional expansion greater in taxonomic generalists? Second, does adding introduced host geography recover contemporary butterfly occurrences that fall outside native host-resource envelopes? Third, how closely does family-level diet breadth correspond to geographic resource breadth, and do apparent specialist–generalist differences in host-contribution architecture persist after the mechanically relevant number of resolved host species is held fixed? Fourth, within contemporary host-resource opportunity, is butterfly geography climatically filtered, and does broader host-family breadth weaken that filtering?

The resource reconstruction is descriptive and uses fixed LepTraits, HOSTS and WCVP inputs. The occurrence validation uses a separately assembled 32-species panel. The host-breadth/climate prediction was generated from an earlier exploratory pilot and then evaluated on a separately frozen panel with prospectively fixed quality gates and a response-blind statistical correction before independent ecological responses were opened. Post-hoc ceiling, taxonomic-family and geographic-distance sensitivities are reported explicitly as robustness analyses rather than as preregistered tests.

We expected anthropogenic host redistribution to enlarge resource opportunity across much of the specialization spectrum, but did not assume that a larger number of host families must produce a larger proportional gain. We further expected the contemporary envelope to recover at least some butterfly occurrences that a native-only envelope misses. For the independent climate test, we prospectively predicted that broader taxonomic diets would weaken climatic filtering within contemporary host-resource opportunity. 

## 2. Methods

### 2.1 Study design

We separated four layers of butterfly specialization and biogeography: taxonomic host breadth, geographic host-resource breadth, anthropogenic expansion of that resource geography, and realization of contemporary resource opportunity by butterflies.

The LepTraits-derived descriptor panel contained 339 species. Resource analyses used 239 species with positive host-family breadth and at least one reconstructed native host-resource WGSRPD3 unit. A conservative subset of 215 species also required the number of resolved host species to equal or exceed LepTraits host-family count. A separate 32-species panel was used for occurrence validation and climate analysis; 24 species passed the climate quality criteria.

### 2.2 Host-resource reconstruction

Larval host records came from HOSTS and plant taxonomy and distributions from WCVP (Robinson et al. 2023; Govaerts et al. 2021). Species-level host names were resolved to accepted WCVP taxa.

For each butterfly, the **native resource envelope** was the union of WGSRPD level-3 units in which any resolved host was native, extant and non-doubtful. The **contemporary resource envelope** additionally retained introduced host distributions. WGSRPD3 is a coarse regional framework, so these envelopes represent geographic resource opportunity rather than local host occupancy.

We calculated native and contemporary resource breadth, introduced-added units, the contemporary/native ratio, log resource expansion, and the introduced share of contemporary geography.

### 2.3 Matched-host null and structural sensitivities

To test whether observed expansion followed automatically from host richness or plant-family composition, we randomized host identity within the 215-species conservative subset. Each null portfolio preserved the focal butterfly's exact number of resolved host species and exact number of hosts in each WCVP plant family. Alternative HOSTS-WCVP species were sampled without replacement within those families, and native and contemporary resource envelopes were reconstructed. We used 1,999 deterministic randomizations.

Primary null endpoints were the number of expanding butterflies, mean and median log expansion, total introduced-added species × WGSRPD3 units, and the association between host-family breadth and log expansion.

We also tested whether the near-zero specialization–expansion relationship was sensitive to the finite 369-unit WGSRPD3 support, butterfly Family, or geographic composition. These post-hoc analyses progressively excluded species with broad native resource envelopes, adjusted rank associations by butterfly Family, and repeated associations within dominant WGSRPD level-1 resource regions.

Host-contribution concentration was retained only as a secondary structural diagnostic. Because effective contributor number and maximum single-host share are bounded by host number, their family-breadth associations were recalculated after holding exact resolved host-species richness fixed.

### 2.4 Independent occurrence validation

Butterfly occurrences for 2010–2026 were obtained from GBIF for the separately assembled 32-species panel and mapped to WGSRPD3 units. We counted occurrence units outside each species' native host-resource envelope and asked how many entered the contemporary envelope after introduced host ranges were retained.

Two post-hoc nulls tested whether recovery exceeded geometric overlap expected from larger envelopes. The first preserved native-envelope size, the number of introduced-added units and the number of outside-native occurrence units, but placed added units uniformly among non-native WGSRPD3 units. The second additionally preserved the number of added units within each WGSRPD level-1 region. Aggregate null distributions used 199,999 hypergeometric draws.

### 2.5 Climate filtering within contemporary resource opportunity

The climate analysis tested the pre-specified prediction that broader host-family diets weaken climatic filtering within contemporary host-resource opportunity. Species entered this analysis only when occurrence coverage, host-taxonomy consistency and within-envelope sampling met fixed quality criteria; 24 species qualified.

Climate was represented by CHELSA BIO1, BIO7, BIO12 and BIO15. Within each species, observed resource units were split deterministically into training and evaluation sets. Climate center and scale were estimated from training occurrences. Held-out observed resource units were compared with effort-supported contemporary host units in which the butterfly was not observed.

The species climate-filtering score was the probability that a never-observed resource unit had greater climatic mismatch from the training niche than a held-out observed unit, with ties receiving half weight. A score of 0.5 is neutral.

The primary test used partial Spearman correlation between host-family breadth and filtering score while controlling contemporary resource breadth, with 9,999 residual permutations and a one-sided negative alternative.

To assess geographic confounding, post-hoc sensitivities restricted comparisons to occupied WGSRPD level-1 regions and matched held-out observed and never-observed resource units by distance to the nearest training-observed region using 250-, 500- and 1,000-km calipers. A species bootstrap summarized effect-size precision.

## 3. Results

### 3.1 Actual host identities produce excess anthropogenic resource expansion

Introduced host distributions expanded reconstructed resource opportunity for 206 of 239 species (86.2%). Aggregate species × WGSRPD3 coverage increased from 26,530 under native-only host distributions to 41,083 under contemporary distributions, an addition of 14,553 units (+54.9%).

In the conservative 215-species subset, this magnitude was substantially greater than expected from host richness and plant-family composition alone. Observed host portfolios produced 13,529 introduced-added species × region units, compared with a matched-null median of 8,046 (95% interval 7,153–8,944; p = 0.0005). Mean log expansion was 0.444 versus a null median of 0.291, and median log expansion was 0.329 versus 0.195 (both p = 0.0005). Forty-eight species individually exceeded their matched null at p <= 0.05.

The number of species with any expansion was less exceptional: 191 species expanded versus a null median of 185 (177–192; p = 0.0705). The strongest non-random signal was therefore the **magnitude** of expansion.

### 3.2 Excess expansion is not a broad-generalist effect

Host-family breadth was nearly unrelated to proportional expansion in the full panel (Spearman rho = 0.008) and conservative subset (rho = 0.015). Under the matched-host null, the observed conservative-panel correlation was also unexceptional: null median rho = 0.004, 95% interval -0.102 to 0.110, two-sided p = 0.7865.

The near-zero relationship persisted after accounting for finite map support and taxonomic structure. Controlling native resource breadth gave a rank association of 0.012; adjusting by butterfly Family gave 0.021. Restricting analyses to species with progressively narrower native envelopes produced correlations of 0.028–0.067.

Regional stratification likewise showed that the global result was not solely driven by Northern America. Within dominant native-resource regions, rho was -0.005 in Northern America (n = 102), 0.013 in Africa (n = 38), 0.029 in Temperate Asia (n = 24), 0.051 in Tropical Asia (n = 9) and 0.089 in Southern America (n = 32). Europe showed a modest positive association (rho = 0.233; n = 34).

### 3.3 Introduced host geography aligns with independent butterfly occurrence

Thirty-one of 32 independent-panel species yielded contemporary GBIF records. Across the panel, 115 butterfly species × WGSRPD3 observations occurred outside native host-resource envelopes. Adding introduced host distributions recovered 66 of these units (57.4%).

Twenty-three species had at least one outside-native occurrence, and 18 recovered at least one unit after introduced hosts were added. The median species-level recovery fraction was 0.60; eight species recovered all outside-native observed units.

Recovery exceeded geometric expectation. A uniform envelope-enlargement null had median recovery of 38 units (95% interval 30–46), while a null preserving each species' broad regional distribution of added host geography had median 49 (43–56). The observed value of 66 exceeded both nulls (p = 5 × 10^-6).

### 3.4 The original portfolio-concentration gradient was largely structural

Before structural adjustment, host-family breadth correlated with effective contributor number (rho = 0.492) and maximum single-host share (rho = -0.486) among 191 expanded conservative-panel species. After exact resolved host-species richness was held fixed, these associations fell to 0.059 and -0.071; conditional permutation p-values were 0.502 and 0.415. Joint adjustment for host-species richness and butterfly Family gave -0.039 and 0.015.

Thus the striking marginal specialist–generalist architecture contrast largely reflects finer-scale host-species portfolio richness and is not treated as an independent mechanism.

### 3.5 Climate-associated filtering persists, but the host-breadth prediction was unsupported

The 24 climate-informative species had a median filtering score of 0.801, and 23/24 exceeded the neutral value of 0.5. Restricting comparisons to occupied WGSRPD level-1 regions gave a median score of 0.846. Distance matching reduced the median to 0.714, 0.750 and 0.750 at 250, 500 and 1,000 km, respectively, but the majority of species remained above 0.5 in each analysis.

The pre-specified prediction that broader host-family diets weaken climate filtering was not supported. The primary partial Spearman correlation was -0.166 (one-sided permutation p = 0.2237). Spatial sensitivities did not recover the predicted negative direction.

The effect estimate was imprecise: a 30,000-replicate species bootstrap gave a 95% interval of -0.583 to +0.261. An approximate effect of |partial rho| ~= 0.50 would be required for 80% power under the current design.

## 4. Discussion

Human redistribution of larval host plants has altered butterfly resource geography at a scale that is visible across hundreds of species. Introduced host distributions increased aggregate reconstructed resource opportunity by 54.9% and expanded the envelope of more than 86% of resource-eligible butterflies. This expansion was not concentrated in family-level generalists, survived finite-geography and butterfly-Family sensitivities, and—most importantly—was reflected in independent butterfly occurrences: introduced hosts recovered 57.4% of contemporary species × region observations that fell outside native host-resource envelopes.

The central ecological result is therefore not a specialist–generalist contrast in portfolio concentration. It is a separation between **taxonomic diet breadth, resource geography and realized distribution** under global change.

### 4.1 Plant redistribution changes the geographic meaning of specialization

A species can be taxonomically specialized yet geographically resource-rich when its few host lineages are widespread. Human transport amplifies this possibility by moving known host plants far beyond their native biogeographic ranges. Consequently, the number of host families does not specify how much geographic opportunity a butterfly gains when plants are redistributed.

The matched-host null reveals the more interesting non-random pattern. Butterflies do not use arbitrary plants drawn from the same host families: their actual host identities are disproportionately associated with species whose introduced distributions enlarge geographic resource opportunity. Matching exact host-species richness and exact plant-family composition reduced expected total added opportunity from the observed 13,529 units to a median of 8,046. Yet the near-zero host-family-breadth gradient was ordinary under the same null. The biologically informative signal is therefore **which host species butterflies use**, not whether butterflies cross many plant-family boundaries.

The near-zero relationship between family breadth and proportional expansion is robust to the obvious finite-area objection. Species with large native envelopes have less room to expand on a 369-unit map, but excluding progressively broader native envelopes did not reveal a hidden generalist advantage. Likewise, the small positive association with absolute added units disappeared after starting resource breadth was controlled. The result is best read as a property of **relative opportunity gain**: broad taxonomic diets do not automatically translate into disproportionate benefit from host globalization.

### 4.2 Independent occurrences show that added resource geography is ecologically relevant

A resource-envelope analysis can otherwise remain purely potential. The independent occurrence comparison provides a stronger bridge to realized biogeography.

More than half of species × region observations lying outside native host-resource envelopes were brought inside the envelope by adding introduced host ranges. Crucially, the aggregate overlap exceeded both an envelope-size null and a stricter null preserving the broad regional distribution of added resource units. The validation therefore cannot be reduced to the fact that contemporary envelopes are larger. This does not prove that the recorded butterfly used the introduced host at that locality, nor that host introduction caused colonization, but it demonstrates that a native-only view of larval resources systematically misses contemporary geographic opportunities that coincide with butterfly presence.

That distinction matters for macroecological analyses that combine consumer distributions with resource distributions. When resources themselves have been redistributed, native resource maps can create apparent consumer–resource mismatches that are partly artifacts of treating present-day interaction opportunity as if plant geography were still native.

### 4.3 Portfolio architecture should not be mistaken for an independent family-breadth mechanism

The original descriptive analysis suggested a striking contrast: one-family butterflies appeared to obtain added opportunity from more concentrated host portfolios than broad generalists. The structural diagnostics show why that pattern must be interpreted cautiously.

Effective contributor number cannot exceed the number of contributing hosts, and the minimum possible largest-host share falls as more hosts are available. Once exact resolved host-species richness was fixed, the apparent family-breadth gradient essentially vanished; jointly controlling butterfly Family did not restore it. Thus the marginal architecture result mostly expresses the finer-scale number of known hosts rather than an independent property of deep taxonomic generalism.

This is still biologically informative because it shows that "host breadth" contains nested levels. But the relevant conclusion is not that family specialists and generalists intrinsically assemble anthropogenic opportunity differently. Rather, species-level host richness determines much of the portfolio geometry that family-level categories only imperfectly summarize.

### 4.4 Resource opportunity is not realized geography

The contemporary host envelope was much broader than observed butterfly geography for most species. Climate-associated mismatch separated held-out observed from never-observed resource units in 23 of 24 species, and this pattern persisted after coarse regional restriction and explicit distance matching.

The accessibility sensitivity is important because climate and distance are correlated in global data. Restricting comparisons geographically reduced the median filtering score under strict distance matching, but did not erase it. The result therefore supports climate as one informative axis of unrealized resource opportunity without claiming that climate is the only filter or that the score is causal.

Dispersal history, habitat, phenology, adult resources, biotic interactions and imperfect detection can all generate residual non-realization. The appropriate interpretation is layered: introduced plants alter where larval resources can occur, while additional ecological filters help determine where butterflies are actually observed.

### 4.5 Taxonomic generalism does not provide a universal shortcut

The prospectively frozen climate prediction was not supported. Broad host-family diets did not detectably weaken climate filtering after resource breadth was controlled, and spatial sensitivities did not reveal the predicted negative relationship.

This negative result is useful but should not be overread. With 24 informative species, the test has limited precision for moderate effects. What the data reject is the idea that family-level generalism supplies an obvious, strong and portable predictor of release from climate-associated filtering in this panel.

Taken together with the resource analysis, the broader message is that family-level diet breadth is not a sufficient proxy for either the geography of resource opportunity or the filters acting on that opportunity.

### 4.6 Implications for global change biogeography

Global change redistributes interacting species as well as climate. Consumer biogeography therefore cannot always be interpreted against static or native-only resource templates.

For butterflies, introduced plants can enlarge potential larval-resource geography across both specialists and generalists. Some of that added opportunity coincides with contemporary butterfly occurrence, while much remains unrealized and environmentally filtered. Similar logic should apply to other consumers whose resources are transported, cultivated, invaded or otherwise redistributed by humans.

A useful next step is to replace coarse opportunity envelopes with local realized host-use data and dated introductions. That would allow direct tests of whether introduced resources facilitate colonization, whether consumers switch among hosts after arrival, and whether resource redistribution changes persistence rather than merely geographic potential. 

## 5. Limitations

HOSTS is incomplete and geographically uneven, and the conservative filter cannot make host-interaction sampling globally uniform. Regional stratification reduces concern that the main breadth–expansion result is solely a North American artifact, but tropical strata remain small. Butterfly Family is only a coarse correction for phylogenetic non-independence and is not equivalent to species-level PGLS.

WCVP introduced status identifies contemporary distribution status rather than the timing or pathway of introduction. We therefore infer changed **resource opportunity**, not a causal historical effect of host introduction on butterfly range expansion.

Occurrence validation is presence-side: a butterfly observation in a region covered by an introduced known host does not prove larval use of that host population. Likewise, GBIF non-observation is not confirmed absence.

WGSRPD3 is coarse and imposes a finite geographic ceiling. Sensitivity analyses indicate that this ceiling does not explain the near-zero specialization gradient, but finer spatial data would improve inference.

The climate test contains only 24 informative species. Its broad bootstrap interval means moderate effects remain plausible despite the unsupported directional prediction.

Finally, the matched-host null, occurrence-overlap nulls, ceiling, Family and spatial sensitivities were added after manuscript review and are therefore post-hoc robustness analyses. They delimit the supported claims but do not convert the study into a fully prospective test of host redistribution.

## 6. Conclusions

Human redistribution of host plants has substantially expanded reconstructed butterfly resource geography: more than 86% of resource-eligible species gained opportunity and aggregate species × region coverage increased by 54.9%. In the conservative panel, actual known host identities generated substantially more expansion than random host sets with the same host-species count and plant-family composition. That excess was not preferentially concentrated in broad family-level generalists.

Independent occurrences show that the added geography is not merely cartographic potential. Introduced host distributions recovered 57.4% of contemporary butterfly species × region observations that lay outside native host-resource envelopes, significantly more than expected after preserving both envelope size and broad regional placement. Yet resource availability alone did not reproduce realized geography: climate-associated filtering persisted after regional and distance controls.

The resulting picture is layered rather than one-dimensional: **host taxonomy describes interaction breadth; plant biogeography determines where resource opportunity exists; human redistribution changes that opportunity; and climate-associated plus other filters shape how much is realized.** 

## Figure legends

**Figure 1. Introduced host distributions expand butterfly resource geography beyond matched-host expectations.** **a**, Aggregate native-only and contemporary resource geography across 239 butterflies; 206 species expanded and aggregate coverage increased by 54.9%. **b**, Observed mean and median log expansion in the conservative 215-species subset compared with 1,999 null portfolios preserving exact host-species richness and plant-family composition (both p = 0.0005). **c**, Observed total introduced-added units (13,529) versus matched-null expectation (median 8,046; p = 0.0005); the host-family-breadth gradient itself was not different from null expectation (p = 0.7865).

**Figure 2. Introduced host geography recovers independent butterfly occurrences beyond structural overlap expectations.** **a**, Fraction of outside-native occurrence units recovered for each of 23 informative species. **b**, Observed aggregate recovery (66/115 units) versus uniform and region-matched structural nulls; both p = 5 × 10^-6.

**Figure 3. The absence of a broad-generalist advantage is robust to finite geographic support and regional composition.** **a**, Host-family breadth versus proportional expansion after progressively excluding broad native resource envelopes. **b**, Within-region associations for species grouped by dominant native host-resource region.

**Figure 4. Climate-associated filtering persists after geographic controls, whereas the predicted host-breadth release is unsupported.** **a**, Filtering scores under the original analysis, within-region restriction and three distance-matching calipers. **b**, Primary partial host-family-breadth effect with bootstrap 95% interval and approximate 80%-power detectable-effect threshold.

## Data and Code Availability

Analysis code and the inputs required to reproduce the reported results and figures are versioned in the study repository. Public source datasets include LepTraits, HOSTS, WCVP, GBIF, CHELSA and WGSRPD. Exact source identities, retrieval rules, robustness-analysis specifications and execution records are documented in the supplementary reproducibility materials rather than the main text. A permanent archival snapshot and DOI will be supplied with the final submission.

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
