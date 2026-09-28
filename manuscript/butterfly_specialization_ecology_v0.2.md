# Anthropogenic host redistribution expands butterfly resource geography, but climate-associated filtering persists

**Working manuscript v0.2 — 2026-09-28**

## Abstract

**Aim:** To quantify how anthropogenic redistribution of larval host plants changes butterfly resource geography, test whether the expanded resource envelope is reflected in contemporary butterfly occurrences, and evaluate whether climate-associated filtering persists within that opportunity.

**Location:** Global.

**Time period:** Contemporary host and butterfly distributions; butterfly occurrences 2010–2026; CHELSA climatology 1981–2010.

**Major taxa studied:** Butterflies and their larval host plants.

**Methods:** For 239 butterflies, we reconstructed native and contemporary larval-resource geography from LepTraits, HOSTS and WCVP and quantified the resource opportunity added by introduced host distributions. We diagnosed finite-geography, host-richness and butterfly-family structure in these comparisons. Using a separately assembled 32-species occurrence panel, we then asked whether introduced hosts recover butterfly occurrences lying outside native host-resource envelopes. Finally, 24 quality-qualified species entered a prospectively frozen climate cross-fit, with post-hoc within-region and distance-matched sensitivities addressing geographic accessibility.

**Results:** Introduced host distributions expanded reconstructed opportunity for 206/239 species (86.2%) and increased aggregate resource geography from 26,530 to 41,083 species × region units (+54.9%). Proportional expansion was nearly unrelated to host-family breadth (Spearman rho = 0.008; butterfly-family-adjusted rank correlation = 0.021) and remained near zero after excluding the broadest native resource envelopes. In the independent occurrence panel, 115 species × region observations lay outside native host-resource envelopes; introduced host ranges recovered 66 (57.4%) of them, affecting 18 of the 23 species with such observations. Climate-associated filtering remained strong within contemporary resource opportunity (median score = 0.801; 23/24 species >0.5) and persisted under 250–1000 km distance matching (median 0.714–0.750). The prospectively predicted weakening of climate filtering with broader host-family diets was not supported (partial rho = -0.166, one-sided permutation p = 0.2237).

**Main conclusions:** Human redistribution of host plants has substantially enlarged butterfly resource opportunity across the taxonomic specialization spectrum, and independent occurrences show that this added geography is ecologically relevant for contemporary distributions. Yet expanded resource opportunity is not equivalent to realized butterfly geography: climate-associated filtering persists even after coarse accessibility controls. Family-level diet breadth is therefore a poor stand-alone predictor of either proportional anthropogenic opportunity gain or release from abiotic filtering.

**Keywords:** butterflies; climate filtering; ecological specialization; global change; host breadth; introduced plants; resource geography 

---

---

## 1. Introduction

Specialization is one of the most familiar axes used to compare herbivorous insects. Butterflies are routinely described as specialists or generalists according to the number or phylogenetic breadth of larval host plants they use, and diet breadth has been linked to geographic range size, range dynamics, diversification and environmental gradients (Slove & Janz 2011; Lancaster 2020; Gross et al. 2026). Yet taxonomic diet breadth and the geography of those resources are not the same ecological quantity.

A butterfly feeding within one plant family may use one host species or dozens, and those hosts may themselves range from local to nearly cosmopolitan. Human transport adds a second source of decoupling because plant geography is no longer determined only by native biogeographic history. Regional studies have shown that introduced plants can become larval hosts and can accompany butterfly persistence or range change (Graves & Shapiro 2003), while recent global work documents extensive human-mediated redistribution of Lepidoptera themselves (Couto et al. 2026). What remains less clear at global scale is how strongly host redistribution has changed the geographic template of potential butterfly resources.

A resource-side reconstruction is biologically useful only if the added geography bears some relationship to contemporary butterfly distributions. A native-only resource envelope may classify an observed butterfly region as host-unavailable even when a known host has been introduced there. Conversely, the presence of a known host does not guarantee butterfly occurrence because climate, dispersal, habitat, phenology and biotic interactions can filter otherwise available resource geography. These processes motivate a layered view in which plant distributions define resource opportunity and additional filters shape its realization.

Here we ask four linked questions. First, how strongly do introduced host distributions expand butterfly resource geography, and is proportional expansion greater in taxonomic generalists? Second, does adding introduced host geography recover contemporary butterfly occurrences that fall outside native host-resource envelopes? Third, how closely does family-level diet breadth correspond to geographic resource breadth, and do apparent specialist–generalist differences in host-contribution architecture persist after the mechanically relevant number of resolved host species is held fixed? Fourth, within contemporary host-resource opportunity, is butterfly geography climatically filtered, and does broader host-family breadth weaken that filtering?

The resource reconstruction is descriptive and uses fixed LepTraits, HOSTS and WCVP inputs. The occurrence validation uses a separately assembled 32-species panel. The host-breadth/climate prediction was generated from an earlier exploratory pilot and then evaluated on a separately frozen panel with prospectively fixed quality gates and a response-blind statistical correction before independent ecological responses were opened. Post-hoc ceiling, taxonomic-family and geographic-distance sensitivities are reported explicitly as robustness analyses rather than as preregistered tests.

We expected anthropogenic host redistribution to enlarge resource opportunity across much of the specialization spectrum, but did not assume that a larger number of host families must produce a larger proportional gain. We further expected the contemporary envelope to recover at least some butterfly occurrences that a native-only envelope misses. For the independent climate test, we prospectively predicted that broader taxonomic diets would weaken climatic filtering within contemporary host-resource opportunity. 

## 2. Methods

## 2. Methods

### 2.1 Study overview

The analysis separates four ecological layers:

1. **Taxonomic host breadth**: LepTraits number of larval host-plant families.
2. **Geographic resource breadth**: the union of WGSRPD level-3 units occupied by resolved host plants.
3. **Anthropogenic resource expansion**: additional resource units obtained when introduced host distributions are retained.
4. **Realized-geography filters**: contemporary butterfly occurrence relative to native and contemporary resource envelopes, followed by climate-associated filtering within contemporary resource opportunity.

The S1 descriptor set contained 339 species, of which 239 had positive host-family breadth and a reconstructed native resource envelope. A conservative host-taxonomy lower-bound subset contained 215 species. The independent occurrence panel contained 32 species; 31 yielded contemporary occurrence records and 24 passed the frozen pre-climate quality gate and were climate-informative.

Host-contribution concentration was retained as a secondary diagnostic rather than a headline mechanism. Because effective contributor number and maximum single-host share are mathematically constrained by the number of contributing hosts, we explicitly tested whether their marginal association with host-family breadth remained after exact resolved host-species richness, and separately butterfly Family, were held fixed. 

### 2.2 Butterfly trait panel

### 2.2 Butterfly trait panel

We used the exact S1 butterfly reconstruction derived from the frozen LepTraits 1.0 input (Shirey et al. 2022). The resource descriptor table contained 339 species with major trait axes, including host-family breadth, wing-size proxies, voltinism, and habitat-affinity traits. Resource analyses were restricted to species with positive host-family breadth and at least one reconstructed native host-resource WGSRPD3 unit, yielding 239 species.

Because host-interaction databases are incomplete, we also defined a conservative host-taxonomy lower-bound subset in which the number of resolved HOSTS-WCVP host species was at least as large as LepTraits host-family count. This condition is necessary but not sufficient for interaction completeness. It yielded 215 species.

### 2.3 Host interaction and plant distribution reconstruction

Larval host records came from the HOSTS database of lepidopteran host plants (Robinson et al. 2023), using the fixed mirror commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`. Plant taxonomy and distribution came from the World Checklist of Vascular Plants (WCVP; Govaerts et al. 2021) through the fixed rWCVPdata v13 snapshot (commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`).

Species-level host names were resolved to accepted WCVP plant identifiers. The **native resource envelope** for each butterfly was the union of WGSRPD3 units in which any resolved host was recorded as native, extant, and non-doubtful. The **contemporary resource envelope** retained both native and introduced records while continuing to exclude extinct and location-doubtful records.

Geographic units followed level 3 of the World Geographical Scheme for Recording Plant Distributions (WGSRPD; Brummitt et al. 2001). WGSRPD3 is intentionally coarse. We therefore interpret these envelopes as regional resource opportunity, not as local host occupancy or confirmed butterfly habitat.

### 2.4 Taxonomic versus geographic specialization

For each of the 239 resource-eligible butterflies, we compared LepTraits host-family count with native host-resource WGSRPD3 breadth. We calculated Spearman rank correlations for the full resource-eligible panel and the 215-species host-taxonomy lower-bound subset.

To quantify discordance between taxonomic and geographic specialization, both variables were converted to percentile ranks. Large discordance was defined a priori in the reconstruction as an absolute rank difference >= 0.5. We also identified conservative "specialist-wide" species: one-family butterflies whose geographic resource breadth was in the upper quartile of the host-taxonomy-adequate panel.

### 2.5 Anthropogenic expansion of resource opportunity

For each butterfly we calculated native resource units, contemporary resource units, introduced-added units, contemporary/native ratio, log resource expansion = log(1 + contemporary) - log(1 + native), and the introduced share of contemporary geography. We summarized these metrics across four host-family strata and calculated descriptive rank associations with host-family breadth.

Three post-hoc structural sensitivities were added after manuscript review. First, because WGSRPD3 contains a finite 369-unit global support, we repeated the family-breadth association after progressively excluding species with the broadest native resource envelopes and calculated rank-adjusted associations controlling native resource breadth. Second, we residualized rank variables by butterfly Family using the Family field in the exact frozen LepTraits table as a coarse correction for taxonomic non-independence. Third, we separated absolute added units from proportional expansion so that large starting resource envelopes could not masquerade as a specialization effect.

These analyses quantify potential geographic resource opportunity generated by host redistribution. They do not establish that butterflies use every introduced host population or that host introduction caused butterfly range expansion.

### 2.6 Independent occurrence validation of introduced-host opportunity

We used the separately assembled 32-species independent panel to ask whether introduced host ranges improve agreement between resource envelopes and contemporary butterfly occurrences. For each species, GBIF records were mapped to WGSRPD3 units. We counted observed units falling outside the native host-resource envelope and then asked how many of those units entered the resource envelope when introduced, extant, non-doubtful WCVP host distributions were added.

The validation endpoint was the fraction of native-envelope-external butterfly species × WGSRPD3 observations recovered by introduced host ranges. We also summarized the number of species for which at least one such unit was recovered and the species-level recovery fraction. This is presence-side validation: it tests whether contemporary host geography can account for observations missed by native host geography, but does not verify larval use at each observed locality.

### 2.7 Host-contribution architecture as a structural diagnostic

For butterflies whose resource envelopes expanded, introduced-added units were fractionally assigned among contributing host species. We calculated maximum single-host fractional share, top-two-host share, effective contributor number and the fraction of added units supported by multiple hosts.

The original descriptive analysis related these quantities to host-family breadth. Because portfolio-concentration statistics are bounded by host number, we then performed two post-hoc diagnostics. We first removed exact resolved host-species-count group means from rank-transformed host-family breadth and architecture metrics. We then repeated the adjustment within joint strata of exact resolved host-species count and butterfly Family. These analyses ask whether family-level diet breadth contains architecture information beyond the finer-scale number of resolved hosts; they do not turn the exploratory decomposition into a causal test.

### 2.8 Taxonomic versus geographic specialization

For each of the 239 resource-eligible butterflies, we compared LepTraits host-family count with native host-resource WGSRPD3 breadth. We calculated Spearman rank correlations for the full panel and the 215-species host-taxonomy lower-bound subset. Both variables were also converted to percentile ranks to identify strong taxonomic–geographic discordance. 

### 2.9 Independent climate-filtering panel

### 2.8 Independent climate-filtering panel

An exploratory 10-species pilot was used only to generate the later host-breadth/climate prediction. Those ten species were excluded from the independent panel.

Before independent occurrence or climate responses were opened, we froze a 32-species panel balanced across four host-family strata (8 species each: 1, 2, 3-5, and 6+ families). Within strata, species were selected deterministically to span native resource breadth, resolved host-species richness, and wing size.

The frozen hypothesis was:

> Within larval host-resource opportunity, climatic filtering of realized butterfly distribution weakens as taxonomic host breadth increases.

### 2.10 Butterfly occurrence acquisition and quality gates

Butterfly occurrences were acquired through the GBIF Occurrence API (GBIF.org 2026) for 2010-2026 using fixed filters requiring coordinates, no flagged geospatial issue, and occurrence status PRESENT. Each species had six deterministic ordinal windows of up to 300 records. Transport failures were handled through checkpointed technical recovery without species replacement or changes to scientific filters.

The initial independent execution was not evaluable because only 10 species passed the frozen pre-climate gate, below the required minimum of 12. Fourteen species had incomplete transport. A post-gate technical recovery was therefore applied uniformly to all and only those 14 partial-transport species. The original quality thresholds, species panel, resource definition, and climate analysis were unchanged. After this completion, 31/32 species had complete occurrence transport; one species remained rejected at the GBIF taxon-resolution stage.

A species passed the pre-climate quality gate only when all of the following held:

- occurrence transport was complete;
- resolved host-species count was at least host-family count;
- no more than 10% of observed WGSRPD3 units lay outside the contemporary host envelope;
- at least five contemporary host units had >=10 GBIF records from other independent-panel butterflies;
- among those effort-supported units, at least two contained the target butterfly and at least two did not.

Twenty-four species passed after uniform transport completion, exceeding the frozen minimum of 12.

### 2.11 Climate cross-fitting

Climate was represented by CHELSA v2.1 1981-2010 variables BIO1, BIO7, BIO12, and BIO15 (Karger et al. 2017, 2021).

Within each quality-qualified species, observed WGSRPD3 units inside the contemporary host-resource envelope were deterministically split 50:50 into training and evaluation sets using a SHA256 species-by-unit rule.

The species climate center and scale were estimated only from occurrence-level climate values in training-observed units. Climate mismatch for each WGSRPD3 resource unit was calculated as the root mean squared standardized deviation from that training climate.

Evaluation compared effort-supported held-out observed units with effort-supported contemporary host units in which the butterfly was never observed. The species-level **climate-filtering score** was

`P(mismatch_never-observed > mismatch_held-out-observed) + 0.5 * P(tie)`.

A value of 0.5 is neutral; larger values indicate that butterfly-unobserved parts of the reconstructed resource envelope are climatically farther from the species' training niche than held-out observed portions.

Species required at least 30 training occurrence records with climate, at least two effort-supported held-out observed units, and at least two effort-supported never-observed units. All 24 pre-climate-qualified species were climate-informative.

### 2.12 Independent test of host-breadth climate release

The primary response was the species climate-filtering score. The predictor was LepTraits host-family count, and contemporary host-resource WGSRPD3 breadth was the frozen control.

The observed statistic was partial Spearman correlation on ranks. The null distribution used 9,999 deterministic, species-identity-invariant Freedman-Lane-style reduced-response residual permutations. The one-sided alternative predicted a negative association: broader host-family breadth should reduce climate filtering. Alpha was 0.05.

The residual-permutation implementation was corrected before independent ecological responses were opened to preserve predictor-control structure and to ensure row-order invariance. The effect definition, panel, direction, alpha, and permutation count were unchanged.


### 2.13 Geographic-accessibility and precision sensitivities

Because climatically mismatched resource units may also be geographically remote, we performed two post-hoc accessibility sensitivities. First, comparisons were restricted to WGSRPD level-1 regions containing a training or held-out observed unit for the focal species. Second, held-out observed and never-observed resource units were greedily matched by distance to the species' nearest training-observed WGSRPD3 region using calipers of 250, 500 and 1000 km. Filtering scores were recomputed on these matched pairs.

We also quantified the precision of the n = 24 host-breadth test. A Fisher-z approximation for a partial correlation with one control was used only as an interpretable precision diagnostic alongside the preregistered permutation p-value.

---

## 3. Results

### 3.1 Introduced host ranges broadly expand resource opportunity without a generalist advantage

When introduced host distributions were retained, 206 of 239 species (86.2%) gained reconstructed geographic resource opportunity. Total species × region units increased from 26,530 under native-only host distributions to 41,083 under contemporary distributions, an addition of 14,553 units (+54.9%).

Proportional resource expansion was almost unrelated to host-family breadth (Spearman rho = 0.008). The conservative 215-species subset gave the same result (rho = 0.015), with 191 species expanding. Median contemporary/native ratios were 1.390, 1.442, 1.411 and 1.361 across the one-, two-, three-to-five- and six-plus-family strata, respectively.

This near-zero relationship was not explained by the finite WGSRPD3 ceiling. The maximum native resource breadth was 305 of 369 units; after restricting the analysis to species below native-breadth thresholds of 250, 200, 168, 150 and 100 units, rho remained between 0.028 and 0.067. Rank-adjustment for native breadth gave a host-family association of 0.012 with log proportional expansion. Although host-family breadth was weakly associated with the absolute number of added units (rho = 0.127), that relationship fell to 0.017 after controlling native breadth. A coarse butterfly-Family correction likewise left the proportional result essentially unchanged (family-adjusted rank correlation = 0.021).

### 3.2 Introduced host ranges recover contemporary butterfly occurrences missed by native host geography

The independent occurrence panel provided direct presence-side validation of the resource reconstruction. Thirty-one of 32 species yielded contemporary GBIF records. Across the panel, 115 butterfly species × WGSRPD3 observations occurred outside native host-resource envelopes. Adding introduced host distributions recovered 66 of those units (57.4%).

Twenty-three species had at least one observed unit outside their native host envelope, and introduced hosts recovered at least one such unit in 18 species. The median species-level recovery fraction among those 23 species was 0.60; eight species had all native-envelope-external observed units recovered. For example, all 22 outside-native observed units of *Pyrgus communis* and all four of *Pieris brassicae* entered the contemporary resource envelope after introduced hosts were included.

Thus the anthropogenic expansion is not only a change in reconstructed plant geography: for a substantial fraction of independent contemporary butterfly observations, introduced host distributions convert a native-host mismatch into host-available geography.

### 3.3 Family-level host breadth only partly tracks native resource geography

Across 239 resource-eligible species, host-family breadth and native geographic resource breadth were positively but weakly associated (Spearman rho = 0.276). In the 215-species host-taxonomy lower-bound subset, the association was rho = 0.398.

Within the conservative subset, 33 species differed by at least 0.5 between their percentile ranks for taxonomic and geographic specialization. Eighteen one-family butterflies nevertheless occupied the upper quartile of geographic resource breadth, including *Hylephila phyleus* (305 native resource units), *Atalopedes campestris* (303) and *Pararge aegeria* (241). Family-level diet breadth therefore does not uniquely determine the geography of larval-resource opportunity.

### 3.4 The apparent specialist–generalist architecture gradient is largely structural

In the unadjusted 191-species expanded subset, host-family breadth correlated with effective contributor number (rho = 0.492) and maximum single-host share (rho = -0.486). However, these concentration metrics are mechanically constrained by the number of available hosts.

After exact resolved host-species count was held fixed on the rank scale, the residual associations fell to 0.059 for effective contributor number and -0.071 for maximum single-host share; conditional permutation tests gave two-sided p = 0.502 and 0.415, respectively. Joint adjustment for exact host-species richness and butterfly Family reduced the corresponding correlations further to -0.039 and 0.015.

The large marginal architecture contrast therefore should not be interpreted as an independent effect of family-level specialization. Instead, it mainly reflects finer-scale host-species portfolio richness. We retain the decomposition as a descriptive account of how added opportunity is assembled, but not as the central biological mechanism.

### 3.5 Climate-associated filtering persists after geographic-distance controls

Twenty-four independent-panel species passed the frozen pre-climate gate and were climate-informative. The original cross-fit showed widespread climate-associated filtering: the median filtering score was 0.801, 23/24 species exceeded the neutral value of 0.5, and 22/24 had positive median mismatch differences between never-observed and held-out observed resource units.

This pattern persisted when coarse accessibility was restricted. Within only WGSRPD level-1 regions containing observed units, the median score was 0.846 and 23/24 species remained above 0.5. Matching held-out observed and never-observed units by distance to the nearest training-observed region produced median scores of 0.714, 0.750 and 0.750 at 250, 500 and 1000 km calipers, with 21/24, 23/24 and 22/24 species above 0.5, respectively.

Geographic distance therefore explains part of the original separation but does not account for the general tendency of climatically mismatched host-resource regions to remain unobserved.

### 3.6 The predicted release from climate filtering in broad generalists was not supported

The prospectively frozen prediction that broader host-family diets weaken climate filtering was not supported. After controlling contemporary resource breadth, the partial Spearman association was -0.166; 2,236 of 9,999 residual permutations were at least as negative as observed (one-sided p = 0.2237).

The post-hoc spatial sensitivities did not rescue the predicted negative direction. Partial associations controlling resource breadth were +0.281 within occupied WGSRPD level-1 regions and +0.319, +0.347 and +0.123 under 250, 500 and 1000 km distance matching, respectively.

The independent panel was nevertheless modest in size. A Fisher-z approximation gives an approximate 95% interval of -0.54 to +0.26 around the observed partial correlation, and an effect of roughly |partial rho| = 0.50 would be required for about 80% power at one-sided alpha = 0.05. The result therefore rejects neither all negative effects nor all biologically moderate effects; it specifically fails to support the preregistered directional prediction at the observed precision. 

## 4. Discussion

Human redistribution of larval host plants has altered butterfly resource geography at a scale that is visible across hundreds of species. Introduced host distributions increased aggregate reconstructed resource opportunity by 54.9% and expanded the envelope of more than 86% of resource-eligible butterflies. This expansion was not concentrated in family-level generalists, survived finite-geography and butterfly-Family sensitivities, and—most importantly—was reflected in independent butterfly occurrences: introduced hosts recovered 57.4% of contemporary species × region observations that fell outside native host-resource envelopes.

The central ecological result is therefore not a specialist–generalist contrast in portfolio concentration. It is a separation between **taxonomic diet breadth, resource geography and realized distribution** under global change.

### 4.1 Plant redistribution changes the geographic meaning of specialization

A species can be taxonomically specialized yet geographically resource-rich when its few host lineages are widespread. Human transport amplifies this possibility by moving known host plants far beyond their native biogeographic ranges. Consequently, the number of host families does not specify how much geographic opportunity a butterfly gains when plants are redistributed.

The near-zero relationship between family breadth and proportional expansion is robust to the obvious finite-area objection. Species with large native envelopes have less room to expand on a 369-unit map, but excluding progressively broader native envelopes did not reveal a hidden generalist advantage. Likewise, the small positive association with absolute added units disappeared after starting resource breadth was controlled. The result is best read as a property of **relative opportunity gain**: broad taxonomic diets do not automatically translate into disproportionate benefit from host globalization.

### 4.2 Independent occurrences show that added resource geography is ecologically relevant

A resource-envelope analysis can otherwise remain purely potential. The independent occurrence comparison provides a stronger bridge to realized biogeography.

More than half of species × region observations lying outside native host-resource envelopes were brought inside the envelope by adding introduced host ranges. This does not prove that the recorded butterfly used the introduced host at that locality, nor that host introduction caused colonization. But it demonstrates that a native-only view of larval resources systematically misses contemporary geographic opportunities that coincide with butterfly presence.

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

First, HOSTS is incomplete and geographically uneven, with strong representation from Europe and North America. The conservative lower-bound filter removes obvious host-count inconsistencies but does not make interaction sampling globally uniform. The butterfly-Family sensitivity addresses taxonomic non-independence only coarsely; it is not a substitute for a dated species-level phylogeny or PGLS.

Second, WCVP introduced status describes distribution status rather than the date, pathway or local abundance of introduction. We therefore infer changed **resource opportunity**, not a historical causal effect of plant introduction on butterfly range expansion.

Third, the occurrence validation is presence-side. A butterfly observation in a region newly covered by an introduced known host does not demonstrate larval use of that introduced population. It shows that the contemporary resource envelope resolves a geographic mismatch present under the native-only envelope.

Fourth, WGSRPD3 is coarse. Availability can be heterogeneous within regions, and the finite 369-unit support creates a formal upper bound on proportional expansion. Ceiling sensitivities indicate that this bound does not explain the near-zero family-breadth association, but finer spatial data would improve inference.

Fifth, GBIF non-observation is not confirmed absence. The climate analysis used an effort proxy based on other panel butterflies, cross-fitting, within-region restriction and distance matching, but remains a macroecological filtering analysis rather than a standardized occupancy model.

Sixth, the independent climate test contains only 24 informative species. The approximate 95% interval around the observed partial association is broad (-0.54 to +0.26), so moderate effects remain plausible even though the prospectively predicted negative association was not supported.

Seventh, several robustness analyses—including ceiling diagnostics, exact host-count adjustment, butterfly-Family adjustment and spatial matching—were motivated after inspection of the manuscript and should be treated as post-hoc sensitivities. They strengthen or delimit existing claims but are not independent confirmatory tests.

Finally, host-contribution concentration metrics are structurally bounded by portfolio size. Our reanalysis shows that the original family-breadth architecture gradient largely disappears once exact resolved host-species richness is held fixed. We therefore do not interpret that marginal gradient as an independent specialization mechanism. 

## 6. Conclusions

Human redistribution of host plants has substantially expanded reconstructed butterfly resource geography: more than 86% of resource-eligible species gained opportunity and aggregate species × region coverage increased by 54.9%. This proportional gain was essentially unrelated to family-level diet breadth and remained so after ceiling and butterfly-Family sensitivities.

Independent occurrences show that the added geography is not merely cartographic potential. Introduced host distributions recovered 57.4% of contemporary butterfly species × region observations that lay outside native host-resource envelopes. Yet resource availability alone did not reproduce realized geography: climate-associated filtering persisted after regional and distance controls.

The resulting picture is layered rather than one-dimensional: **host taxonomy describes interaction breadth; plant biogeography determines where resource opportunity exists; human redistribution changes that opportunity; and climate-associated plus other filters shape how much is realized.** 

## Figure legends

**Figure 1. Introduced host distributions broadly expand reconstructed butterfly resource opportunity.** **a**, Native versus contemporary host-resource breadth for 239 resource-eligible butterflies. Overall, 206/239 species expanded and aggregate species × WGSRPD3 units increased from 26,530 to 41,083 (+54.9%). **b**, Log proportional resource expansion by host-family breadth class. Host-family breadth was nearly unrelated to proportional expansion (Spearman rho = 0.008), and this result was stable to native-breadth and butterfly-Family sensitivities.

**Figure 2. Introduced host geography recovers independent butterfly occurrences missed by native host envelopes.** For the 23 independent-panel species with at least one observed WGSRPD3 unit outside the native host envelope, bars show the fraction of those units recovered when introduced host ranges are retained. Across species, 66/115 native-envelope-external species × region observations (57.4%) were recovered; 18/23 species recovered at least one unit and the median species-level recovery fraction was 0.60.

**Figure 3. Taxonomic host breadth only partly predicts geographic resource breadth.** Host-family count is plotted against native host-resource WGSRPD3 breadth for the 215 species meeting the conservative host-taxonomy lower-bound criterion. Eighteen one-family butterflies fall in the upper quartile of geographic resource breadth. Spearman rho = 0.398.

**Figure 4. Climate-associated filtering persists after geographic accessibility controls.** Species filtering scores for the 24 climate-informative butterflies are shown under the original cross-fit, same-WGSRPD-level-1 restriction and distance matching at 250, 500 and 1000 km. Median scores were 0.801, 0.846, 0.714, 0.750 and 0.750, respectively; the neutral expectation is 0.5.

**Figure 5. The preregistered prediction that broader diets weaken climate filtering was not supported.** Climate-filtering score is plotted against host-family count, with contemporary resource breadth represented separately. The frozen primary test gave partial Spearman rho = -0.166 (one-sided permutation p = 0.2237; n = 24). Post-hoc spatial sensitivities did not recover the predicted negative direction, and the approximate 95% interval for the primary effect was -0.54 to +0.26. 

## Data and Code Availability

## Data and Code Availability

All analysis code, frozen scientific protocols, input identities, archived workflow definitions, and result receipts supporting this manuscript are versioned in the `chocho` ecology repository. The manuscript figures are rendered with `scripts/render_butterfly_specialization_manuscript_figures.py`; the historical workflow definition used for the audited figure run is preserved at `provenance/workflows/butterfly-specialization-manuscript-figures-v01.yml`, with exact execution provenance retained in the repository.

The frozen LepTraits input used by the reconstruction is vendored in the repository and hash-pinned. Larger external sources (HOSTS, WCVP, GBIF, CHELSA and WGSRPD) are reconstructed or retrieved from their original providers using pinned repository commits, file hashes, query rules and workflow-run provenance recorded in the repository. A permanent archival snapshot of the code, frozen inputs required for reproducibility, source-artifact receipts and result receipts will be deposited before submission. Its public DOI should be provided on the separate title page and in the final public manuscript; the double-anonymous review manuscript should instead use an anonymized reviewer-access link whose landing page and metadata do not identify the authors.

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

Integrated synthesis:

- `benchmarks/exploratory/butterfly_specialization_ecology_synthesis_v0.1.json`
- `docs/exploratory/BUTTERFLY_SPECIALIZATION_ECOLOGY_SYNTHESIS_V0_1.md`
