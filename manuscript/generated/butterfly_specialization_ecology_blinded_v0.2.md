# Anthropogenic expansion of butterfly resource geography is concentrated in network-prominent host plants
**Running title:** Host-plant prominence and resource expansion


## Abstract

**Aim:** To test how human redistribution of larval host plants changes butterfly resource geography, what accounts for that expansion, and whether the added geography aligns with contemporary butterfly occurrences.

**Location:** Global.

**Time period:** Contemporary distributions; butterfly occurrences 2010–2026.

**Major taxa studied:** Butterflies and larval host plants.

**Methods:** We reconstructed native and contemporary host-resource geography for 239 butterflies. Matched null models preserved butterfly host richness and plant-family composition; stricter sensitivities additionally matched host native-range breadth and weighted candidate plants by use by other Lepidoptera in HOSTS. A 32-species panel originally stratified for the climate test was reused to test whether introduced-host geography recovers butterfly occurrences outside native host envelopes beyond structural overlap expectations.

**Results:** Introduced hosts expanded resource geography for 206/239 species and increased aggregate species × region coverage by 54.9%. A simple same-family null suggested excess expansion, but the excess disappeared after matching native range and weighting other-Lepidoptera use (12,853 observed added units vs null median 12,971; p = 0.584). Across 8,909 host plants, network degree remained positively associated with geographic expansion within plant-family and native-breadth strata (rank r = 0.295; p = 0.0002). Host-family breadth remained nearly unrelated to butterfly proportional expansion after phylogenetic correction. Introduced hosts recovered 66/115 butterfly observations outside native host envelopes; species-level recovery remained above region-matched expectation after resampling and removal of the largest contributor.

**Main conclusions:** Human redistribution of host plants has greatly expanded butterfly resource geography. Network-prominent host plants are disproportionately redistributed, whereas broad butterfly diets do not predict proportional gain. Taxonomic diet breadth is therefore a poor proxy for anthropogenic resource gain.

**Keywords:** butterflies; biotic redistribution; host plants; introduced species; resource geography; specialization

---

## 1. Introduction

Specialization is one of the most familiar axes used to compare herbivorous insects. Butterflies are routinely described as specialists or generalists according to the number or phylogenetic breadth of larval host plants they use, and diet breadth has been linked to geographic range size, range dynamics, diversification and environmental gradients (Slove & Janz 2011; Lancaster 2020; Gross et al. 2026). Yet taxonomic diet breadth and the geography of those resources are not the same ecological quantity.

A butterfly feeding within one plant family may use one host species or dozens, and those hosts may themselves range from local to nearly cosmopolitan. Human transport adds a second source of decoupling because plant geography is no longer determined only by native biogeographic history. Regional studies have shown that introduced plants can become larval hosts and can accompany butterfly persistence or range change (Graves & Shapiro 2003), while recent global work documents extensive human-mediated redistribution of Lepidoptera themselves (Couto et al. 2026). Host-plant use is also highly uneven across plant lineages: a relatively small set of plant genera can support a large fraction of Lepidoptera diversity, and butterfly–plant networks repeatedly concentrate interactions on particular host clades (Narango et al. 2020; Braga et al. 2018). This raises a biogeographic question that taxonomic diet breadth alone cannot answer: are plants that are prominent across Lepidoptera–host networks also disproportionately redistributed by humans?

A resource-side reconstruction is biologically useful only if the added geography bears some relationship to contemporary butterfly distributions. A native-only resource envelope may classify an observed butterfly region as host-unavailable even when a known host has been introduced there. Conversely, the presence of a known host does not guarantee butterfly occurrence because climate, dispersal, habitat, phenology and biotic interactions can filter otherwise available resource geography. These processes motivate a layered view in which plant distributions define resource opportunity and additional filters shape its realization.

Here we ask three primary questions. First, how strongly do introduced host distributions expand butterfly resource geography, and does any apparent host-identity excess remain after accounting for plant native-range breadth and network-wide host prominence? Second, is proportional expansion greater in butterflies with broader family-level diets? Third, does contemporary introduced-host geography recover butterfly occurrences that fall outside native host-resource envelopes more often than structural overlap predicts? We retain a smaller climate analysis as a secondary test of whether available resource geography is fully realized.

The resource reconstruction is descriptive and uses fixed LepTraits, HOSTS and WCVP inputs. The 32-species panel was assembled separately for the climate test, using eligibility and analysis rules fixed before its ecological responses were evaluated; the occurrence analysis is a secondary use of that panel. All later null-model, ceiling, phylogenetic and geographic sensitivities are identified as post-hoc robustness analyses.

We expected anthropogenic host redistribution to enlarge resource opportunity across much of the specialization spectrum, but did not assume that broader taxonomic diets must receive a larger proportional gain. We also expected the contemporary resource envelope to recover at least some butterfly occurrences missed by a native-only envelope. 

## 2. Methods

### 2.1 Study design

We separated four layers of butterfly specialization and biogeography: taxonomic host breadth, geographic host-resource breadth, anthropogenic expansion of that resource geography, and realization of contemporary resource opportunity by butterflies.

The LepTraits-derived descriptor panel contained 339 species. Resource analyses used 239 species with positive host-family breadth and at least one reconstructed native host-resource WGSRPD3 unit. A conservative subset of 215 species also required the number of resolved host species to equal or exceed LepTraits host-family count. A separate 32-species panel was used for occurrence validation and climate analysis; it had been selected to balance the subsequent climate test, not optimized for validation. Twenty-four species passed the climate quality criteria.

### 2.2 Host-resource reconstruction

Larval host records came from HOSTS and plant taxonomy and distributions from WCVP (Robinson et al. 2023; Govaerts et al. 2021). Species-level host names were resolved to accepted WCVP taxa.

For each butterfly, the **native resource envelope** was the union of WGSRPD level-3 units in which any resolved host was native, extant and non-doubtful. The **contemporary resource envelope** additionally retained introduced host distributions. WGSRPD3 is a coarse regional framework, so these envelopes represent geographic resource opportunity rather than local host occupancy.

We calculated native and contemporary resource breadth, introduced-added units, the contemporary/native ratio, log resource expansion, and the introduced share of contemporary geography.

### 2.3 Matched-host nulls and structural sensitivities

We first tested whether observed expansion followed automatically from host richness and plant-family composition. In the conservative panel, null portfolios preserved each butterfly's exact number of resolved host species and exact number of hosts in each WCVP plant family while replacing host identities with alternative HOSTS-WCVP species from the same families.

Because HOSTS is not an exchangeable plant pool, we then added stricter post-hoc nulls addressing two plant-level biases. Candidate hosts were either (1) matched probabilistically to each observed host's native WGSRPD3 breadth, (2) weighted by the number of other Lepidoptera species recorded using that plant in HOSTS, or (3) subjected to both constraints simultaneously. Observed hosts themselves were excluded from candidate pools. The combined null therefore asks whether actual butterfly hosts retain excess anthropogenic expansion after accounting for both starting geographic breadth and network-wide host prominence/recording intensity. This analysis retained 207 species with sufficient alternative same-family host pools and used 999 randomizations.

We also tested whether the specialization–expansion relationship was sensitive to the finite 369-unit WGSRPD3 support, geographic composition, and taxonomic non-independence. Post-hoc analyses progressively excluded broad native resource envelopes, repeated associations within dominant WGSRPD level-1 regions, fitted butterfly Family as a random intercept, and used the Kawahara et al. (2023) time-calibrated butterfly phylogeny for exact species matches. Finally, across all 8,909 HOSTS-WCVP plants with valid native and contemporary distributions, we related the number of distinct Lepidoptera consumers to plant geographic expansion, controlling native resource breadth and testing the association within plant-family × native-breadth-quintile strata by permutation.

Host-contribution concentration was retained only as a secondary structural diagnostic because effective contributor number and maximum single-host share are bounded by host number.

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

### 3.1 Anthropogenic resource expansion is concentrated in network-prominent host plants

Introduced host distributions expanded reconstructed resource opportunity for 206 of 239 species (86.2%). Aggregate species × WGSRPD3 coverage increased from 26,530 under native-only host distributions to 41,083 under contemporary distributions, an addition of 14,553 units (+54.9%).

A simple same-family matched null initially suggested that actual host identities were unusually expansion-prone. In the 215-species conservative panel, observed host portfolios produced 13,529 introduced-added species × region units versus a null median of 8,046, and no value among 1,999 randomizations was as large as observed for total, mean, or median expansion.

That inference did not survive the stricter plant-level bias test. Among 207 species with sufficient alternative host pools, observed total added units were 12,853. Matching candidate hosts on native-range breadth alone still produced a much smaller null median of 7,347 (0/999 randomizations as large as observed). However, weighting candidates by use by other Lepidoptera produced a null median of 15,236, exceeding the observed total. When native-range matching and other-Lepidoptera-use weighting were combined, the null median was 12,971 (95% interval 12,084–13,843; p = 0.584). Mean log expansion was 0.439 versus a combined-null median of 0.409 (p = 0.073), and median log expansion was 0.333 versus 0.309 (p = 0.079).

Thus actual butterfly host identities are associated with greater anthropogenic expansion than arbitrary same-family plants, but this excess is largely captured by plants that are already geographically broad and prominent across the Lepidoptera–host network. The plant-level analysis confirmed that this prominence axis is itself associated with redistribution: across 8,909 host plants in 278 families, log network degree correlated with log expansion at rho = 0.313; the partial rank association controlling native breadth was 0.267, and the within-family/native-breadth-stratum correlation was 0.295 (permutation p = 0.0002). Because HOSTS usage also reflects recording intensity, these associations cannot separate ecological host commonness from database visibility.

### 3.2 Expansion is not preferentially greater in broad family-level generalists

Host-family breadth was nearly unrelated to proportional expansion in the full panel (Spearman rho = 0.008) and conservative subset (rho = 0.015). The stricter host-bias nulls did not reveal a hidden positive specialization gradient.

The near-zero relationship persisted after accounting for finite map support and taxonomic structure. Controlling native resource breadth gave a rank association of 0.012. A butterfly-Family random-intercept model estimated a standardized rank effect of 0.021 (95% CI -0.111 to 0.154; p = 0.752). The Kawahara et al. (2023) tree contained exact matches for 124/239 panel species: Brownian rank-PGLS gave beta = 0.140 (95% CI -0.036 to 0.317; p = 0.119), while estimated Pagel lambda was 0.033 and the corresponding rank-PGLS effect was -0.051 (-0.234 to 0.131; p = 0.581). Restricting analyses to progressively narrower native envelopes produced correlations of 0.028–0.067.

Regional stratification likewise showed that the global result was not solely driven by Northern America. Within dominant native-resource regions, rho was -0.005 in Northern America (n = 102), 0.013 in Africa (n = 38), 0.029 in Temperate Asia (n = 24), 0.051 in Tropical Asia (n = 9) and 0.089 in Southern America (n = 32). Europe showed a modest positive association (rho = 0.233; n = 34).

### 3.3 Introduced host geography aligns with separate-panel butterfly occurrence

Thirty-one of 32 independent-panel species yielded contemporary GBIF records. Across the panel, 115 butterfly species × WGSRPD3 observations occurred outside native host-resource envelopes. Adding introduced host distributions recovered 66 of these units (57.4%).

Twenty-three species had at least one outside-native occurrence, and 18 recovered at least one unit after introduced hosts were added. The median species-level recovery fraction was 0.60; eight species recovered all outside-native observed units. Seven of those eight complete-recovery cases involved only 1–4 outside-native units; the exception was *Pyrgus communis*, with 22/22 units recovered.

Recovery exceeded geometric expectation at both pooled and species levels. A uniform envelope-enlargement null had median recovery of 38 units (95% interval 30–46), while a null preserving each species' broad regional distribution of added host geography had median 49 (43–56); none of 199,999 draws reached the observed 66.

Treating species rather than species × region units as the replication level gave the same conclusion. Mean species-level recovery was 0.583 versus a region-matched null median of 0.385 (95% interval 0.308–0.470), with no exceedance in 99,999 draws. A species-cluster bootstrap placed the pooled recovery fraction between 0.349 and 0.774 (95% percentile interval).

*Pyrgus communis* contributed 22 of the 66 recovered units, but removing it still left 44/93 units recovered (47.3%) versus a region-matched null median of 30 (24–36), again with no exceedance in 99,999 draws. Removing *Pieris brassicae* left 62/111 (55.9%) versus a null median of 46 (40–52). Recovery was heterogeneous rather than universal: *Erynnis tristis* recovered 0/11 outside-native units and *Historis acheronta* recovered 5/13.

### 3.4 The original portfolio-concentration gradient was largely structural

Before structural adjustment, host-family breadth correlated with effective contributor number (rho = 0.492) and maximum single-host share (rho = -0.486) among 191 expanded conservative-panel species. After exact resolved host-species richness was held fixed, these associations fell to 0.059 and -0.071; conditional permutation p-values were 0.502 and 0.415. Joint adjustment for host-species richness and butterfly Family gave -0.039 and 0.015.

Thus the striking marginal specialist–generalist architecture contrast largely reflects finer-scale host-species portfolio richness and is not treated as an independent mechanism.

### 3.5 Secondary climate analysis

Climate-associated filtering remained evident in the 24 climate-informative species. Median filtering scores were 0.801 in the original analysis and 0.714–0.750 after 250–1,000-km distance matching; 21/24, 23/24 and 22/24 species remained above the neutral score of 0.5 at the three matching scales.

The pre-specified prediction that broader host-family diets weaken filtering was not supported (partial Spearman rho = -0.166, one-sided p = 0.2237) and was imprecisely estimated (bootstrap 95% interval -0.583 to +0.261). Full diagnostics are shown in Supplementary Figure S1.

## 4. Discussion

Human redistribution of larval host plants has altered butterfly resource geography at a scale that is visible across hundreds of species. Introduced host distributions increased aggregate reconstructed resource opportunity by 54.9% and expanded the envelope of more than 86% of resource-eligible butterflies. This expansion was not concentrated in family-level generalists, survived finite-geography and butterfly-Family sensitivities, and—most importantly—was reflected in independent butterfly occurrences: introduced hosts recovered 57.4% of contemporary species × region observations that fell outside native host-resource envelopes.

The central ecological result is therefore not a specialist–generalist contrast in portfolio concentration. It is a separation between **taxonomic diet breadth, resource geography and realized distribution** under global change.

### 4.1 Plant-level host prominence, not broad butterfly diets, structures resource globalization

A species can be taxonomically specialized yet geographically resource-rich when its few host lineages are widespread. Human transport amplifies this decoupling by moving host plants beyond their native biogeographic ranges.

The null-model sequence identifies an important boundary around the apparent host-identity signal. Actual butterfly hosts generated much more anthropogenic expansion than uniformly sampled same-family alternatives, and the difference remained after native-range matching alone. But the excess disappeared once alternatives were also weighted by use by other Lepidoptera. The complementary plant-level analysis shows why: network degree was positively associated with geographic expansion even after native breadth and plant-family structure were controlled. Previous work has shown that Lepidoptera interactions are strongly concentrated on a minority of high-value host lineages (Narango et al. 2020) and that repeated use of particular host clades shapes butterfly–plant network structure (Braga et al. 2018). Our result adds a geographic-global-change dimension: the same plant-level prominence axis is associated with how strongly host ranges have expanded beyond their native geography.

This does not show that network prominence itself causes plant redistribution. Widely used hosts may be common, weedy or cultivated, but they may also simply be better documented in HOSTS. The defensible result is therefore that **plant-level host prominence and starting biogeography account for the apparent butterfly-specific host-identity excess**. Network-prominent plants are disproportionately redistributed, but the data do not support the stronger claim that butterfly host choice independently targets unusually globalized plants.

At the same time, broad taxonomic diets do not translate into disproportionate benefit from host globalization. The near-zero family-breadth association survived ceiling, regional, random-effect and species-level phylogenetic sensitivities. Anthropogenic opportunity is therefore better understood through the identities and biogeography of resource species than through a one-dimensional specialist–generalist axis.

### 4.2 Secondary occurrence validation shows that added resource geography is ecologically relevant

A resource-envelope analysis can otherwise remain purely potential. This secondary occurrence comparison provides a bridge to realized biogeography.

More than half of species × region observations lying outside native host-resource envelopes were brought inside the envelope by adding introduced host ranges. Crucially, the result persisted when butterfly species were treated as the replication unit and when the largest-contributing species was removed. The validation therefore cannot be reduced to the fact that contemporary envelopes are larger or to a single high-leverage butterfly. This does not prove that the recorded butterfly used the introduced host at that locality, nor that host introduction caused colonization, but it demonstrates that a native-only view of larval resources systematically misses contemporary geographic opportunities that coincide with butterfly presence.

That distinction matters for macroecological analyses that combine consumer distributions with resource distributions. When resources themselves have been redistributed, native resource maps can create apparent consumer–resource mismatches that are partly artifacts of treating present-day interaction opportunity as if plant geography were still native.

### 4.3 Portfolio architecture should not be mistaken for an independent family-breadth mechanism

The original descriptive analysis suggested a striking contrast: one-family butterflies appeared to obtain added opportunity from more concentrated host portfolios than broad generalists. The structural diagnostics show why that pattern must be interpreted cautiously.

Effective contributor number cannot exceed the number of contributing hosts, and the minimum possible largest-host share falls as more hosts are available. Once exact resolved host-species richness was fixed, the apparent family-breadth gradient essentially vanished; jointly controlling butterfly Family did not restore it. Thus the marginal architecture result mostly expresses the finer-scale number of known hosts rather than an independent property of deep taxonomic generalism.

This is still biologically informative because it shows that "host breadth" contains nested levels. But the relevant conclusion is not that family specialists and generalists intrinsically assemble anthropogenic opportunity differently. Rather, species-level host richness determines much of the portfolio geometry that family-level categories only imperfectly summarize.

### 4.4 Resource opportunity is not realized geography

The occurrence validation shows that contemporary host geography matters, but it does not imply that all available resource geography is realized. In the smaller climate panel, climatically mismatched resource units were less often observed even after geographic matching. At the same time, broader host-family diets did not detectably weaken this filtering, and the effect estimate was too imprecise to exclude moderate associations.

We therefore treat climate as a secondary boundary on resource realization rather than a co-equal mechanism in this paper. Dispersal history, habitat, phenology, adult resources, biotic interactions and imperfect detection can also generate unrealized resource opportunity.

### 4.5 Implications for global change biogeography

Global change redistributes interacting species as well as climate. Consumer biogeography therefore cannot always be interpreted against static or native-only resource templates.

For butterflies, introduced plants can enlarge potential larval-resource geography across both specialists and generalists. Some of that added opportunity coincides with contemporary butterfly occurrence, while much remains unrealized and environmentally filtered. Similar logic should apply to other consumers whose resources are transported, cultivated, invaded or otherwise redistributed by humans.

A useful next step is to replace coarse opportunity envelopes with local realized host-use data and dated introductions. That would allow direct tests of whether introduced resources facilitate colonization, whether consumers switch among hosts after arrival, and whether resource redistribution changes persistence rather than merely geographic potential. 

## 5. Limitations

HOSTS is incomplete and geographically uneven, and the conservative filter cannot make host-interaction sampling globally uniform. Regional stratification reduces concern that the main breadth–expansion result is solely a North American artifact, but tropical strata remain small. Species-level PGLS was possible for only 124/239 butterflies directly represented in the Kawahara et al. (2023) tree; unmatched species were excluded rather than phylogenetically imputed.

WCVP introduced status identifies contemporary distribution status rather than the timing or pathway of introduction. We therefore infer changed **resource opportunity**, not a causal historical effect of host introduction on butterfly range expansion.

Occurrence validation is presence-side and opportunistic: the 32-species panel was constructed for a stratified climate test, not specifically for validating introduced-host opportunity. A butterfly observation in a region covered by an introduced known host does not prove larval use of that host population, and GBIF non-observation is not confirmed absence. The panel also includes species with documented non-native populations, including *Pieris brassicae* (Phillips et al. 2020); excluding that species did not change the occurrence-overlap conclusion.

WGSRPD3 is coarse and imposes a finite geographic ceiling. Sensitivity analyses indicate that this ceiling does not explain the near-zero specialization gradient, but finer spatial data would improve inference.

The climate test contains only 24 informative species. Its broad bootstrap interval means moderate effects remain plausible despite the unsupported directional prediction.

Finally, the matched-host nulls, occurrence-overlap nulls, ceiling, phylogenetic and spatial sensitivities were added after manuscript review and are therefore post-hoc robustness analyses. Other-Lepidoptera consumer count is an imperfect joint proxy for ecological host prominence and HOSTS recording intensity, so the stricter null can show that this axis explains the simple host-identity excess but cannot distinguish biological commonness from database bias.

## 6. Conclusions

Human redistribution of host plants has substantially expanded reconstructed butterfly resource geography: more than 86% of resource-eligible species gained opportunity and aggregate species × region coverage increased by 54.9%. The apparent excess expansion of actual host identities relative to a simple same-family null disappeared after accounting jointly for host native-range breadth and network-wide Lepidoptera use. Across the wider host-plant pool, network-prominent plants were themselves disproportionately redistributed after native breadth and plant-family structure were controlled. Anthropogenic resource gain is therefore better understood through plant-level prominence and biogeography than through broad butterfly diets or a butterfly-specific preference for unusually globalized hosts.

The secondary occurrence analysis shows that the added geography is not merely cartographic potential. Introduced host distributions recovered 57.4% of contemporary butterfly species × region observations that lay outside native host-resource envelopes, significantly more than expected after preserving both envelope size and broad regional placement. Yet resource availability alone did not reproduce realized geography: climate-associated filtering persisted after regional and distance controls.

The resulting picture is layered rather than one-dimensional: **host taxonomy describes interaction breadth; plant biogeography determines where resource opportunity exists; human redistribution changes that opportunity; and climate-associated plus other filters shape how much is realized.** 

## Figure legends

**Figure 1. Anthropogenic resource expansion is widespread and concentrated in network-prominent host plants.** **a**, Native versus contemporary host-resource breadth across 239 butterflies; 206 species expanded and aggregate coverage increased by 54.9%. **b**, Observed total introduced-added resource units in the 207-species strict-null subset compared with native-range-matched, other-Lepidoptera-use-weighted, and combined native-range-plus-use null expectations. The observed total exceeded the native-range null but not the combined null (12,853 observed vs median 12,971; p = 0.584). **c**, Across 8,909 host plants, Lepidoptera consumer degree was positively associated with plant geographic expansion; the association persisted after native-breadth adjustment and within plant-family × native-breadth strata (r = 0.295, permutation p = 0.0002). Network degree is treated as a joint proxy for ecological host prominence and recording intensity.

**Figure 2. Introduced host geography recovers butterfly occurrences beyond structural overlap expectations.** **a**, Fraction of outside-native occurrence units recovered for each of 23 informative species. **b**, Observed aggregate and species-level recovery versus structural nulls, with leave-one-species-out sensitivity. Removing *Pyrgus communis* left 44/93 units recovered versus a region-matched null median of 30; no region-matched null draw equalled or exceeded the observed value in 99,999 replicates.

**Figure 3. The absence of a broad-generalist advantage is robust to finite geographic support and regional composition.** **a**, Host-family breadth versus proportional expansion after progressively excluding broad native resource envelopes. **b**, Within-region associations for species grouped by dominant native host-resource region.

**Supplementary Figure S1. Climate-associated filtering persists after geographic controls, whereas the predicted host-breadth release is unsupported.** **a**, Filtering scores under the original analysis, within-region restriction and three distance-matching calipers. **b**, Primary partial host-family-breadth effect with bootstrap 95% interval and approximate 80%-power detectable-effect threshold.

## Data and Code Availability

For double-anonymous review, analysis code and the inputs required to reproduce the reported results and figures will be supplied through an anonymized reviewer-access link:

**Reviewer link:** [ANONYMIZED REVIEW LINK]

Public source datasets include LepTraits, HOSTS, WCVP, GBIF, CHELSA and WGSRPD. Exact source identities and retrieval rules are documented in the anonymized reproducibility materials. A permanent public archival snapshot and DOI will replace this statement in the final version.

## References (working)

- Braga, M. P., Guimarães, P. R., Wheat, C. W., Nylin, S. & Janz, N. 2018. Unifying host-associated diversification processes using butterfly–plant networks. *Nature Communications* 9: 5155. https://doi.org/10.1038/s41467-018-07677-x
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
- Narango, D. L., Tallamy, D. W. & Shropshire, K. J. 2020. Few keystone plant genera support the majority of Lepidoptera species. *Nature Communications* 11: 5751. https://doi.org/10.1038/s41467-020-19565-4
- Phillips, C. B., Brown, K., Green, C., Toft, R., Walker, G. & Broome, K. 2020. Eradicating the large white butterfly from New Zealand eliminates a threat to endemic Brassicaceae. *PLoS ONE* 15: e0236791. https://doi.org/10.1371/journal.pone.0236791
- Rashid, S., Wessely, J., Hausharter, J., Moser, D., Gattringer, A., Fiedler, K., Hülber, K. & Dullinger, S. 2026. Food Plant Availability Constrains Climatic Niches of Host-Specialized Europe-Centred Butterflies. *Diversity and Distributions* 32: e70245. https://doi.org/10.1111/ddi.70245
- Robinson, G. S., Ackery, P. R., Kitching, I., Beccaloni, G. W. & Hernández, L. M. 2023. HOSTS - a Database of the World's Lepidopteran Hostplants [Data set]. Natural History Museum. https://doi.org/10.5519/havt50xw
- Shirey, V., Larsen, E., Doherty, A. et al. 2022. LepTraits 1.0: A globally comprehensive dataset of butterfly traits. *Scientific Data* 9: 382. https://doi.org/10.1038/s41597-022-01473-5
- Slove, J. & Janz, N. 2011. The relationship between diet breadth and geographic range size in the butterfly subfamily Nymphalinae: a study of global scale. *PLoS ONE* 6: e16057. https://doi.org/10.1371/journal.pone.0016057
