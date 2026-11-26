# Supporting Information

**Port-associated introduced floras reveal a weaker oceanic isolation gradient in mapped butterfly host resources**

Supporting methods and results are organized by analysis. Tables S1–S9 and Figures S3.1, S3.2 and S5.1 provide additional results, sensitivity analyses and geographical-transfer checks. Appendix S8 records data and code provenance.

## Appendix S1. Data sources and study populations

**Data sources and study design**

Butterfly–host associations were obtained from HOSTS (Robinson et al., 2023) and matched to native and introduced plant distributions in WCVP v13 (21 May 2024; Govaerts, 2024). Geographic units follow WGSRPD3 (Brummitt, 2001). The 155-species prediction and 224-species community analyses use the same GBIF occurrence download (GBIF.org, 2026; DOI 10.15468/dl.umqmsf) but different samples and outcomes. HOSTS does not exhaust all actual interactions, and WCVP regional origin assignments do not establish introduction dates. The separate historical and port-buffer GBIF acquisitions are not covered by this DOI. Versions and access limitations are listed in Appendix S8.

**Table S1.** Source populations, taxonomic identities and geographic grain.

| Analysis | Cohort | Measured response |
| --- | --- | --- |
| Fixed known-host opportunity | 239 HOSTS butterflies × 369 WGSRPD3 botanical regions | Presence of ≥1 reported native or native+introduced larval host plant, not butterfly occurrence |
| Island–mainland host opportunity | 126 insular and 237 mainland botanical regions (6 mixed separate) | Fixed-source host availability and new gap filling |
| Port × whole plant flora | The same 126 botanical islands, WPI ports × WCVP accepted species | Adjusted contemporary native/introduced plant richness and 7,875 paired island floristic similarities |
| Marine resistance | Source-matched 104 island botanical regions × 5,355 reused flora pairs | Native vs introduced plant Jaccard and undirected sea electrical resistance (not measured currents) |
| Independent later butterfly indexing | 155 butterflies, 208,482 species×cell test cases | Retrospective later indexed GBIF records, not population change |
| Independent GBIF community | 224 butterflies, 169 original primary regions, 14,196 reused region pairs | 2018–2024 recorded butterfly Jaccard conditional on earlier records and observation proxies |
| External island endemic | Independent 41 extant Macaronesian butterflies incl. 18 endemic | Named host records, independent wild/garden larval use, catalogue coverage |

## Appendix S2. Host-resource reconstruction and resource gains

**Methods**

For each geography class, native coverage B is the proportion of all butterfly–region combinations with a mapped native host, and gap filling F is the proportion of initially empty combinations newly covered by introduced hosts. Relative gain over native coverage is F(1−B)/B. The area × Level1 matched comparison uses an introduced-attributable contemporary share, a distinct denominator from gap filling.

Butterflies were ordered by global native-resource extent (number of WGSRPD3 regions with at least one mapped native host) and divided into groups of 60, 60, 60 and 59. Alphabetical species-name order resolved the tie at the upper boundary. Each butterfly was evaluated across all 126 island regions. Within groups, pooled gap filling divides added combinations by initially empty combinations; proportional gain divides additions by native-covered combinations. Each newly covered combination was allocated equally among its contributing host plants.

**Results**

**Table S2.** Native and contemporary host coverage in island and mainland regions.

| Source fixed host coverage | Island (126 botanical regions) | Mainland (237) |
| --- | --- | --- |
| Butterfly×region potential cells | 30,114 | 56,643 |
| Native known-host covered | 4,829 | 21,215 |
| Added by introduced known hosts | 5,741 | 8,686 |
| Current known-host covered | 10,570 | 29,901 |
| Native baseline coverage | 16.04% | 37.45% |
| Fraction of initially missing resource combinations newly filled | 22.71% | 24.52% |
| Relative gain over native host coverage | +118.89% | +40.94% |
| Introduced-only source coverage as fraction of current covered combinations | 54.31% | 29.05% |

For the full 369-region panel, native-only coverage was 26,530 and contemporary coverage 41,083 butterfly–region combinations, an increase of 14,553 (+54.9%); 206/239 species gained at least one region. The absolute island–mainland coverage gap narrowed from 21.4 to 17.7 percentage points (contemporary coverage 35.1% and 52.8%, respectively). Host-family count versus proportional gain gave Spearman ρ=0.008 (bootstrap 95% −0.111 to +0.128; strict equivalence not confirmed). Species-equal, spatial and crop-exclusion sensitivity analyses are documented in the reproducibility records (Appendix S8). The matched area × Level1 introduced-attributable contemporary share difference was +10.7 percentage points; raw gap filling was 22.7% on islands and 24.5% on mainlands.

Among **670** source-accepted introduced host contributors, **38 (5.7%)** provided half of all 14,553 added butterfly×region opportunities using equal split of credits within each newly covered combination (plant Gini 0.7774, top-ten cumulative contribution 25.1%). The original 239 butterfly source native-host geographical breadth predicted absolute island additions (tied-rank Spearman ρ=+0.578896); broadest 59 species gained 2,271/5,741 combinations (39.56%), narrowest 60 gained 644/5,741 (11.22%). **82/97** one-host-family butterflies gained new global potential regions; excluding FAO crop hosts still left +46–47% global resource increase. These are mapped opportunities, not larval abundance, food quality or individual fitness.

**Table S3.** Island resource-gap filling by global native-resource extent.

| Global native-resource extent group | Butterflies | Native island units | Initially empty island units | Added units | Gaps filled | Relative gain |
| --- | --- | --- | --- | --- | --- | --- |
| Narrowest | 60 | 339 | 7,221 | 644 | 8.9% | +190.0% |
| Second | 60 | 710 | 6,850 | 1,077 | 15.7% | +151.7% |
| Third | 60 | 1,238 | 6,322 | 1,749 | 27.7% | +141.3% |
| Broadest | 59 | 2,542 | 4,892 | 2,271 | 46.4% | +89.3% |

The broadest-group gap-filling estimate remains 46.0–46.9% when both species tied at 169 global native regions are assigned to the upper or lower group, respectively. The difference in mean island opportunities between the broadest and narrowest groups grows in absolute units (37.4 native to 65.2 contemporary regions per butterfly); the Gini coefficient of island opportunities across all 239 butterflies nevertheless declines from 0.462 to 0.427. We therefore distinguish unequal absolute gains from an increase in overall relative inequality, which is not observed here. The low native baseline explains proportionally greater island than mainland gains despite similar gap filling, whereas broadly supplied butterflies fill a larger fraction of their island gaps but the narrowest group has a greater proportional gain. These are post-result descriptive accounting comparisons, not independent tests of causal selection by ports, dispersal ability, dietary specialization, population growth or colonization.

Source-constructed butterfly-resource island-pair Jaccard increased +0.23667, compared with mainland +0.15044, but the Level1×island-class × fixed species/region gain-margin null left island-specific excess only **+0.00119**. This is not observed butterfly community homogenization.

## Appendix S3. Island isolation, port-associated plant supply and marine geography

### S3.1. Plant richness, composition and transport-related covariates

**Methods and results: port associations**

We related log(1 + native or introduced whole-flora richness) to log₂(1 + WPI port count) on the same 126 botanical island regions, controlling for log area and TDWG Level1. The introduced coefficient was +0.2854 (HC3 95% interval +0.0622 to +0.5087); the native coefficient was +0.0954 (−0.0774 to +0.2683). The paired introduced-minus-native contrast was +0.1900 (−0.0187 to +0.3987). Nine overlapping leave-Level1-out models retained positive introduced slopes (+0.202 to +0.380), as did the 91 port-positive regions. These associations do not establish historical ship transport or an origin-selective port effect.

For the 98 regions with botanical recording-volume data, the introduced-flora model used unstandardized log₂(1 + WPI port count), log area, log(1 + native richness), log(1 + total indexed GBIF Plantae records, 1700–2026) and Level1 effects. The separate documented introduced-host model used the native-host baseline; an additional specification conditioned on total introduced richness and used native whole-flora richness instead. Permutation inference used 4,999 port-label permutations within Level1 × area halves. We also fitted the direct within-region origin contrast, log(1 + introduced richness) − log(1 + native richness), using the same unstandardized port predictor and log area, Level1 and additional geographic/effort covariates as indicated in Figure S3.1.

**Results: recording-volume and origin-contrast sensitivity**

In the 98-region analysis, the port coefficient for introduced flora was +0.24448 (permutation p = 0.0012); using pre-1900 botanical records yielded +0.22250 (p = 0.0014). All nine leave-Level1-out estimates were positive (+0.187 to +0.283). The port coefficient for documented introduced-host richness was +0.16391 (p = 0.011); conditioning on total introduced richness yielded −0.04153 (p = 0.093). These are distinct specifications rather than a formal mediation test. For the direct origin contrast, the 98-region effort-adjusted estimate was +0.24291 (permutation p = 0.0036; HC3 p = 0.106; nine-division clustered p = 0.104). The 126-region estimates were +0.19001 (clustered p = 0.137) and +0.19345 after adjustment for mainland distance and latitude (p = 0.140). The corresponding 98-region geographic models without and with recording volume yielded +0.25997 (p = 0.0904) and +0.25331 (p = 0.0898; 95% interval −0.0494 to +0.5560). All clustered intervals include zero (Figure S3.1). The 28 excluded regions were disproportionately large and port-rich. The permutation and clustered analyses use different uncertainty assumptions; neither establishes that ports mediate the introduced-only host-share gradient or carried particular plants.

![Figure S3.1](figS33.png)

Figure S3.1. Sensitivity of the direct port–origin contrast. The response is log(1 + introduced plant richness) − log(1 + native plant richness), with unstandardized log₂(1 + port count) as predictor. All models include log area and TDWG Level1. Geographic adjustment adds mainland distance and absolute latitude; effort adjustment adds indexed plant-record volume. Points and whiskers are estimates and nine-Level1 clustered 95% t intervals. Blue indicates 126 regions; orange indicates the common 98-region effort cohort. The five estimates are not independent replications; all intervals include zero.

**Methods: marine predictors and reciprocal origin comparisons**

The coastline mask uses GSHHG version 2.3.7 (Wessel & Smith, 2017). Separate OSCAR checks accessed legacy 1/3-degree, approximately five-day gridded velocity estimates from a NOAA ERDDAP mirror (ESR, 2009). That mirror labels its version 2017.0; byte equivalence to the NASA product remains unverified. OSCAR Version 2 access retrieved catalogue metadata only. These velocity estimates are separate from the undirected resistance geometry.

Reciprocal origin comparisons and competing transport predictions. A stronger seawater-related signal in native flora and a stronger port association in introduced flora are contrasting predictions of the two-layer hypothesis, not an assumption of disjoint pathways. We applied each predictor to both origins to evaluate origin specificity. Native and introduced floras are reciprocal observational comparators, not strict causal negative controls: human transport may redistribute native species within their native ranges, and introduced species may subsequently disperse naturally.

**Results: plant composition and marine resistance**

Among 7,875 native-flora comparisons and 7,874 defined introduced-flora comparisons from the same 126 regions (one empty-union introduced pair excluded), median introduced-flora Jaccard was **0.046030**, native **0.004487** (~10.26×); this is contemporary composition of different species pools, not tenfold change through time. Separate 104-island original single-island holdout electrical resistance yielded introduced **+4.239%**, native **−0.110%** plant-Jaccard predictive MSE reduction. Strict both-island-endpoints-unseen 1,089 dyads in five source folds yielded introduced **+2.735%** (descriptive −1.649,+5.508%) versus native **−0.044%** (−0.984,+0.431%); paired difference interval crosses zero. A *different* shortest water path helps both origin classes. Alternative Simpson similarity can reverse origin ordering. The electric resistance is undirected model geometry, not observed currents.

![Figure S3.2](figS32.png)

Figure S3.2. Undirected marine electrical resistance and floristic prediction. The 104-region source comprises 5,355 island pairs. Relative MSE improvements beyond the original geographic/richness baseline are shown for one island held out and for both endpoints unseen (1,089 test pairs in five folds). Whiskers are descriptive resampling intervals; for the original one-island test the reported interval endpoints are approximate. In the strict test, both origins improved in three of five folds; the introduced-minus-native interval crossed zero. Resistance represents uniform-ocean geometry, not measured currents. A distinct shortest-water-path model provided information for both origins (Appendix S3).

**Methods: crossed geographic validation**

In the original 104-island crossed validation, both origin models used the same 5,355 eligible dyads, a baseline (B) containing geodesic distance, same botanical Level1, latitude gap and accepted-richness margins, and additional shortest-water graph geometry (G) or the product of contemporary log1p WPI port counts (P). GP includes both. Positive values below indicate percentage reductions in prediction MSE relative to the stated reference; the first two rows of contrasts use the original one-island-held-out source, the next two use the stricter five-fold test with both island endpoints excluded from training (1,089 test dyads). Intervals are the original descriptive resampling intervals, not independent natural experiments.

**Results: crossed validation**

**Table S4.** Crossed marine-geometry and port prediction: origin-specific relative MSE improvements.

| Validation and incremental predictor | Native MSE improvement (%) | Introduced MSE improvement (%) |
| --- | --- | --- |
| One island unseen: G after P | +2.2298 [1.0881, 3.3812] | +1.8941 [0.6788, 3.1342] |
| One island unseen: P after G | −0.0817 [−0.1171, −0.0528] | +0.1640 [−0.3388, 0.6744] |
| Both islands unseen: G after P | +1.8556 [−1.9402, 4.1334] | +1.9511 [−1.5178, 3.8374] |
| Both islands unseen: P after G | −0.2683 [−0.9960, −0.0523] | −0.0521 [−1.2901, 0.4880] |

Paired descriptive contrasts used identical held-island resamples for both origins (9,999 draws) or all 3,125 ordered five-fold resamples in the strict test. The native-minus-introduced G-after-P contrast was +0.3357 percentage points (resampling interval −1.0314 to +1.6456); the introduced-minus-native P-after-G contrast was +0.2457 points (−0.2800 to +0.7706). Both remained positive in every leave-one-held-island-out loss aggregation, but the resampling intervals spanned zero. Strict double-unseen contrasts were −0.0955 (−4.3838 to +2.8962) and +0.2162 (−0.7947 to +1.0607), respectively. These post-result paired diagnostics compare relative MSE improvements within each origin; shared dyad endpoints and overlapping training sets preclude interpreting them as independent causal tests. Thus shortest-water geometry provides nominal predictive information for both origin classes, whereas port counts do not establish additional introduced-flora pair-similarity prediction. The port result for island-level introduced species richness is a different response and remains positive; It does not identify pairwise shipping connectivity. Direct origin contrasts are not established by comparing whether separate intervals exclude zero. The actual 2020 AIS overlap comprised just seven verified island pairs, insufficient for global origin–destination shipping inference. Shortest-water geometry and uniform-ocean circuit resistance are separate predictors, neither measuring currents.

**Results: directional-current diagnostics and dispersal traits**

The independent legacy NOAA 2005 directional reversal diagnostic also tests whether apparent current–flora correspondence distinguishes direction. Of 12 model-positive directions, only ten had evaluable reverse physics; none of those ten had reverse model-positive physics, yet source-native/recipient-introduced taxon co-listing occurred in nine directions each way. The mean forward-minus-reverse plant Jaccard difference was +0.00181 (post-hoc exact sign-flip p = 0.186). In the matched-control diagnostic only two of 12 positive directions had at least three source-specific distance/richness controls; observed co-listing in both did not exceed the control expectation. The original launch selection and model censoring limit these comparisons. The subsequent higher-resolution HYCOM test (S3.5) likewise did not provide robust directional specificity. These diagnostics did not establish directional specificity; they do not rule out marine dispersal.

Among 22 regions with EuDiS trait information (Vargas et al., 2023), the native-minus-introduced share of plants with recorded seawater-dispersal morphology was +2.9515 percentage points (descriptive region-resampling interval +2.0395 to +4.5323; 17 positive, four negative, one tie). Coverage was predominantly European. The broader water-dispersal category in Tonga (Fall & Drezner, 2011, Table 1) showed the opposite direction (native 137/364; post-European introduced 137/309). These are different trait definitions, not matched replications.

### S3.2. Isolation gradients of fixed-239 butterfly known-host opportunities

**Methods**

For the original 126 WCVP/WPI botanical islands and 239-species HOSTS×WCVP region counts (**4,829 native-only**, **10,570 contemporary**, **5,741 introduced-only additional** butterfly×region opportunities), approximate mainland distance was the nearest WGSRPD3 boundary separation to 237 botanical mainland regions (chosen in longitude–latitude and converted to WGS84 geodesic endpoint distance with ±360 longitude wrapping). This is **neither** the 104-island undirected ocean circuit **nor** an exact ellipsoidal polygon-distance optimization.

Each log1p plant/host count was standardized and regressed against standardized log2(1+modern ports), log1p(distance), log region area, absolute centroid latitude and TDWG Level1 intercepts (126 islands).

**Results**

Native whole-flora distance β=−0.083 (HC3 p=0.106) was uncertain, whereas introduced whole-flora port β=+0.350 (HC3 p=0.00865; nine-cluster p=0.049). Native-host distance β=−0.1676 (HC3 p=0.000095; nine-Level1 cluster p=0.014) and contemporary-host distance β=+0.0192 (HC3 p=0.751) produced paired standardized slope shift +0.1868 (HC3 95% +0.0749,+0.2986). Critically, the nine-Level1 CR1/t8 interval for that shift crossed zero (95% −0.0172,+0.3907; p=0.0677). A non-significant contemporary slope does not establish exact absence of an isolation relationship.

For **123** host-covered islands, introduced-only potential share was fitted as log((introduced-only+0.5)/(native-covered+0.5)), one row per island: isolation β=+0.4251 (HC3 95% +0.2446,+0.6055; nine-cluster t8 95% +0.1121,+0.7381; p=0.0140). Equal-island fitted shares were **55.2%, 65.6%, 71.5%** at 100, 1,000 and 4,000 km; these are **not** measured butterfly dietary dependence or ratios of the separate fitted host-count curves. Within-Level1 999-draw bootstrap interval for paired distance shift was +0.095,+0.295, all 126 leave-one-island-out slopes were positive (+0.170,+0.203), and all nine Level1 leaveouts remained positive (+0.133,+0.220). Log port counts and log area correlated r=+0.807 (raw counts/area r=+0.271); geographical and historical dependence remain.

This analysis concerns mapped known-host opportunities for 239 butterflies across 126 botanical island regions. It does not measure butterfly presence, ship-borne plant arrival or demographic rescue. An earlier exploratory analysis with a different sample is not used here (Appendix S8).

### S3.3. Additional spatial autocorrelation and MEM10 sensitivity

**Methods and results**

Source-island WGSRPD3 region centroids (not necessarily individual island centroids) were linked using a prespecified 2,000 km great-circle binary neighbourhood, with 350 undirected edges and five unconnected regions. Plain model residual Moran's I used 1,999 global and within-TDWG-Level1 permutations.

**Results: residual spatial dependence**

The native-host model had residual I=+0.132 (within-Level1 p=0.001), while the introduced-only host opportunity log-odds model had residual I=+0.123 (within-Level1 p=0.003), indicating that the plain HC3 regressions do not fully absorb geographic dependence.

As an additional exploratory spatial-filtering sensitivity, the top **10 positive eigenvectors** of the centered spatial weight matrix were added to the original fixed covariates without choosing vectors by outcome association. Under MEM10, the native-host mainland-distance coefficient remained negative (**β=−0.1330, HC3 p=0.0204**) and the contemporary-host distance coefficient was uncertain (**β=+0.0709, p=0.3170**). The paired standardized current-minus-native slope shift was positive (**β=+0.2039, HC3 p=0.0043**) and the introduced-only share log-odds distance association remained positive (**β=+0.3949, HC3 p=0.000786**). The port coefficient for introduced whole flora was also positive (**β=+0.3877, p=0.00560**). The *nine-Level1 cluster t8* paired-slope test remains uncertain (p=0.0677): MEM10+HC3 and nine-cluster t8 are distinct sensitivity designs and **do not constitute independent biological confirmations**. Five disconnected regions, broad multi-island centroids, spatial residual structure and correlated ports/area remain limits.

The spatial weights, permutation results and model outputs are documented in Appendix S8. No historical introduction or butterfly population data enter these models.

### S3.4. Mainland isolation and botanical observation-effort index sensitivity

**Methods**

We separately tested whether the original 126-island Figure 3a source-defined HOSTS×WCVP potential-resource association was solely a correlate of uneven botanical recording. The fixed 239-butterfly responses and approximate WGSRPD3 nearest-mainland distances were independently reconstructed from the same WCVP/port table, 369-region HOSTS panel and pinned TDWG polygons. We then linked the pre-existing GBIF all-Plantae index counts (1700–2026; historical pre-1900 subset separately) to the **98 of 126 botanical island regions with a successfully indexed plant-effort record**. For transparent cohort control, **all models below are re-estimated on the identical 98 islands** (all 98 resource-active), with fixed standardized log2(1+ports), log island region area, absolute centroid latitude and nine TDWG Level1 intercepts; the added covariate is standardized log1p(all indexed plant records). Response counts were standardized by log1p *within the common 98-island panel*, while the introduced-only host-share response remained region-level log odds with a 0.5 continuity adjustment.

**Results**

**Table S5.** Island resource isolation gradients with and without botanical recording effort.

| Same 98-island source | No plant-record effort: distance β, 9-Level1 t8 p | All-plant-record effort adjusted: distance β, 9-Level1 t8 p |
| --- | --- | --- |
| Native documented host opportunities | −0.1957, p=0.0310 | −0.2162, p=0.0147 |
| Contemporary documented host opportunities | +0.0062, p=0.9522 | −0.0201, p=0.8406 |
| Paired contemporary-minus-native standardized slope | +0.2019, p=0.1207 | +0.1961, p=0.1390 |
| Introduced-only share of currently covered known-host opportunities (log odds) | +0.5091, p=0.0111 | +0.5247, p=0.0096 |

Using only pre-1900 indexed botanical records instead of the all-period count yielded native distance **β=−0.2104, nine-cluster p=0.0215** and introduced-only share **β=+0.5105, p=0.0088**, with the paired slope difference still uncertain (β=+0.1809, p=0.1206). The corresponding HC3 results and nine-cluster t8 intervals are included in the aggregate source receipt.

This is an observation-index proxy sensitivity, not an identified or full-population recording-effort adjustment. The 28 omitted botanical regions have no usable original indexed plant-effort measure rather than observed zero effort. They are not missing at random: the separate data showed disproportionately large and port-rich omitted regions. Results on the common 98 cannot be extrapolated to all 126 or described as independent biological replication. The original full-126 evidence and its conservative paired-slope uncertainty (nine-cluster p=0.0677) remain primary; no new historical ship/cargo, butterfly occurrence or recruitment evidence is inferred.

### S3.5. Exploratory directional-current evidence and limits

**Methods**

Using the same WCVP origin assignments and EuDiS trait coding, we identified 17 source regions with adequate within-family trait-positive and trait-negative native plant pools. The historical OSCAR pilot lacked sufficient coastal and day-30 support; we therefore used the independent HYCOM+NCODA GLBu0.08 experiment 19.1 reanalysis (source product; 2005 daily archive), approximately 8-km horizontal resolution, 15-m depth, six release dates in May–October 2005. The fixed roster comprised three separated wet-grid starts per region/date within 50 km of original GSHHG wet anchors (306 trajectories). Midpoint integration used 6-hour steps, a 30-day horizon and a land-mask check along segments of at most 2 km. Missing velocities and land contact terminated integration; no missing trajectory was imputed as zero transport.

**Results**

All 306 starts were reconciled; 52 reached day 30, from 12 source regions. The other 254 trajectories ended at masked/missing velocities (247) or land contact (7); source retrieval/boundary errors were zero. Across 1,236 directed comparisons with 104 recipient regions, we related family-standardized trait-positive-minus-negative listing excess to modeled progress toward the recipient representative point. Source and recipient fixed effects, cubic log initial distance, same Level1 and absolute latitude difference were included, with equal source weights and equal dates within source. For each 100 km of modeled progress, the native excess increased by 0.297 percentage points (full-design two-way clustered descriptive 95% interval 0.112–0.481), introduced excess was −0.101 (−0.223–0.022), and their paired difference was +0.397 (0.222–0.572). Requiring two valid dates retained the paired estimate +0.406 (0.236–0.576; 11 sources, 1,133 comparisons). These are exploratory spatial associations conditional on model-supported trajectories.

Sensitivity analyses constrain the interpretation. Omitting Japan changed the paired estimate to +0.006 and native estimate to −0.215; Japan accounted for 77.9% of residual predictor sum of squares. The paired estimate was −0.059 within 5,000 km (220 comparisons), −0.374 within 3,000 km (105 comparisons), and −0.752 within Level1 (121 comparisons); all intervals included zero. The 1,500-km cohort failed the minimum support gate. In 4,096 source-wise geometric orientation assignments, 33.8% had a partial correlation at least as large as the observed positive paired contrast. This placebo is geometric, not an exchangeable causal null or reversed-current simulation. The positive global association therefore provides no direction-specific support; the result is sensitive to the geographic composition of the comparisons, and the geometric diagnostic does not distinguish the observed direction from alternative orientations.

A separate fixed-anchor diagnostic identified 15 directed comparisons compatible with the computed path lengths among completed trajectories, involving eight sources. No completed trajectory newly ended within 50 km of another fixed recipient anchor. One incomplete trajectory, released from Sardinia on 26 September 2005, ended about 22.9 km from a Corsica anchor (initial distance 149.4 km, computed path length 208.2 km). This is model-supported proximity, not coast-wide arrival, viable seed dispersal or historical establishment. The narrow support prevented a new ecological fit. These exploratory directional and proximity diagnostics provide no additional direction-specific evidence for the two-layer interpretation, which rests on the separate native–introduced morphology, port association and undirected marine-geometry results.

The 2005 model, present-day floras, nonrandom trait coverage and conditioning on complete water tracks limit historical inference. EuDiS is centered on European flora and only 17 source regions meet the trait-coverage gate; additional releases would not resolve this source-coverage constraint. No further threshold or model searches were used for these results. Public-input redistribution conditions and reviewer access remain unresolved.

### S3.6. Taxonomic resolution of the native-host isolation gradient

**Methods**

We reconstructed the original fixed 239 butterfly × 126 island-region panel from the pinned HOSTS and WCVP v13 inputs used in S3.2–S3.4. All 239 per-butterfly global native and contemporary opportunity totals and all 126 island-level native and contemporary counts matched the original Q1 sources exactly (4,829 native and 10,570 contemporary island butterfly–region units). Exact documented host coverage used the dedicated original resolved host distributions. The host edges contained 1,706 resolved plant IDs and 2,596 distinct butterfly–plant pairs; 37 resolved IDs were outside the strict accepted-species whole-flora background and retained with their original name and family metadata. Genus and family presence used the accepted-species background plus those original resolved hosts occurring locally. A butterfly counted as covered at a given rank if at least one of its documented hosts or corresponding higher-rank taxa was mapped in that region. These are nested potential taxonomic matches, not evidence that every congener or family member is edible. All 126 regions were included for counts, including Bouvet Island with no mapped flora; composition and native range summaries with undefined denominators used 125 regions.

We applied the S3.2 standardized log1p count model, HC3 intervals and nine-Level1 t8 inference, the same ten geography-selected Moran eigenvectors, and a combined MEM10 plus clustered sensitivity. Paired standardized slope contrasts were fitted on the same islands, not inferred from differing significance. We also modeled continuity-corrected log ratios of covered butterflies at species versus higher ranks, mean fraction of each butterfly's documented family set represented, and the log odds of native plant species belonging to families containing any documented focal host. The latter is a whole-flora composition measure distinct from whether a butterfly's host family occurs. These were exploratory diagnostic endpoints, not independent tests of actual feeding. All coefficients below use standardized isolation distance; responses are standardized log1p counts for count/paired rows, log ratios or log odds where stated.

**Results**

**Table S6.** Taxonomic-resolution and spatial-filtering sensitivity of mainland-distance coefficients.

| Endpoint | Base β | HC3 p | Level1 t8 p | MEM10 β | MEM10 HC3 p | MEM10 + t8 p |
| --- | --- | --- | --- | --- | --- | --- |
| Native species coverage | −0.1676 | 0.000095 | 0.0140 | −0.1330 | 0.0204 | 0.1095 |
| Native genus coverage | −0.0928 | 0.1087 | 0.2132 | −0.0628 | 0.3169 | 0.4038 |
| Native family coverage | −0.0126 | 0.7793 | 0.7869 | +0.0051 | 0.9236 | 0.9245 |
| Native species minus family slope | −0.1550 | 0.000237 | 0.00809 | −0.1382 | 0.00879 | 0.0587 |
| Native species minus genus slope | −0.0748 | 0.0666 | 0.0881 | −0.0702 | 0.1038 | 0.1241 |
| Native log((species+0.5)/(family+0.5)) | −0.2565 | 0.000023 | 0.00892 | −0.2215 | 0.00721 | 0.0793 |
| Native log((species+0.5)/(genus+0.5)) | −0.1881 | 0.000145 | 0.00713 | −0.1678 | 0.0114 | 0.0748 |
| Native plant species in host families (log odds; 125) | −0.2639 | 0.000195 | 0.00967 | −0.1726 | 0.0157 | 0.0281 |
| Contemporary species coverage | +0.0192 | 0.7508 | 0.7847 | +0.0709 | 0.3170 | 0.3552 |
| Contemporary genus coverage | +0.0573 | 0.3487 | 0.3022 | +0.0895 | 0.1771 | 0.2418 |
| Contemporary family coverage | +0.0340 | 0.4272 | 0.3827 | +0.0454 | 0.3689 | 0.3788 |
| Species contemporary minus native slope | +0.1868 | 0.00107 | 0.0677 | +0.2039 | 0.00431 | 0.0657 |
| Genus contemporary minus native slope | +0.1501 | 0.00008 | 0.0919 | +0.1523 | 0.00205 | 0.1068 |
| Family contemporary minus native slope | +0.0466 | 0.0458 | 0.2452 | +0.0403 | 0.2420 | 0.4015 |

The nearest and farthest distance quartiles each contained 32 regions. Median coverage counts were 74.5 versus 6.5 butterflies at species level, 158 versus 54.5 at genus level and 227 versus 189.5 at family level. These are unadjusted summaries. The paired native species–family slope contrast had a nine-Level1 interval of −0.257 to −0.053; the combined MEM10-plus-cluster interval was −0.283 to +0.006. The direct species–genus contrast was uncertain. Family coverage may saturate: the standardized mean documented-family fraction had a weaker negative association (β −0.110, cluster p = 0.125). Excluding all butterflies with any Fabaceae or Poaceae host (120 retained) left native species coverage negatively associated with isolation (β −0.218, cluster p = 0.00480) and family coverage uncertain (β −0.013, p = 0.849); MEM10 plus clustering yielded p = 0.0586 for species coverage.

On the same 98 islands as S3.4, adding the plant-record-volume proxy retained negative species coverage (β −0.216, cluster p = 0.0147), species/family log ratio (β −0.298, p = 0.0140) and representation of native species belonging to host families (β −0.246, p = 0.00640). Family and genus coverage also had negative point estimates (β −0.156 and −0.193; HC3 p = 0.032 and 0.003; cluster p = 0.121 and 0.074). Thus absence of a full-cohort family effect does not establish that family-level disharmony is absent. Botanical indexing volume is not direct effort in documenting butterfly–host interactions.

Across the 125 nonempty native floras, the standardized logit share of native plants recorded in only one WGSRPD3 region increased with isolation (β +0.328, cluster p = 0.00988; MEM10 + t8 p = 0.0493). This is a geographic-range proxy rather than verified island endemism. On identical rows, adjustment for that share attenuated the species/genus log ratio slope from −0.179 to −0.113 and the species/family slope from −0.244 to −0.216 (adjusted cluster p = 0.156 and 0.067). These comparisons do not establish mediation. The stronger association at species than family rank, the composition signal, and partial range-size attenuation suggest that both identity and floristic composition structure the known-host gradient. They cannot distinguish missing documentation from ecological incompatibility. Contemporary slopes and direct origin contrasts do not establish exact removal of isolation or a uniquely species-specific effect of introductions. This analysis is source-coupled to the main resource map, not independent larval evidence.

The reconstructed counts agree with the fixed 239-species resource panel; this analysis does not establish novel feeding associations, establishment or demographic effects. Input and code limitations are recorded in Appendix S8.

## Appendix S4. Host-specific prediction controls and catalogue coverage

**Methods and results: shared-input controls**

On the original 104 botanical islands, H=668 introduced WCVP plant taxa registered as hosts of the 239 butterflies; O=7,219 introduced taxa **not registered** as hosts (not confirmed true non-hosts). Under original geography E controls, H-after-O improved prediction of *WCVP/HOSTS-derived potential butterfly resource composition* by +17.609% MSE. With 199 plant-family×island-occupancy-matched label draws, observed source H-after-O +16.365% versus median matched null +7.699% (9/199 null scores at least as large, descriptive rank .050). Contemporary whole alien flora similarity improved source-potential butterfly-resource prediction +26.363% against a 225 serial-state fixed source-gain-margin null median +18.974% (excess +7.389 points). These models share WCVP/HOSTS taxonomic input between predictors and potential response and **are not independently observed butterfly evidence**.

### S4.1. Figure 4b: known-host documentation-richness proxy sensitivity

**Methods**

The original Figure 4b association on the fixed **239 butterflies** was between global native known-host resource breadth and the number of newly source-supported island butterfly–region opportunities (Spearman **ρ=+0.578896**). An important ascertainment alternative is that better-studied butterflies have more documented host plants, creating larger reconstructed native and introduced resource maps. The source-frozen 239-butterfly HOSTS/WCVP metrics and corresponding 126-island species-level source counts were therefore paired by identical accepted butterfly species names. We used the **number of distinct resolved documented host-plant species** and **number of documented host-plant families** as imperfect proxies for catalogue coverage; neither equals the number of independent field observations or publications.

**Results**

Across the 239 species, the documented host-plant species count correlated with native global resource breadth (**ρ=+0.626**) and added island resources (**ρ=+0.573**). Residualizing tied midranks of the predictor and outcome against tied midranks of documented host-plant count alone retained a partial rank correlation **+0.344**; additionally controlling documented host-family count retained **+0.348**. A fixed-seed (20261121) 2,000-resample *species* bootstrap returned a descriptive percentile interval **+0.219 to +0.465**. A 2,000-draw permutation of island additions **within five predefined documented-host-count tiers** (1, 2, 3–4, 5–9, 10+; sample counts 29, 35, 43, 51, 81), recomputing the rank adjustment each draw, gave a descriptive two-sided tail probability **p=0.0005**. Omitting each butterfly separately gave adjusted correlations **+0.334 to +0.368**. The source-taxonomically more adequate subset (**n=215**) gave **+0.377** after both controls; butterflies with exactly one documented host family (**n=97**) gave **+0.463** after host-count control; butterflies with positive island additions (**n=198**) retained **+0.247** after both controls.

The follow-up used the same pinned 140,485-row HOSTS catalogue and WCVP name-resolution code. All 239 butterfly names and resolved host counts matched the frozen cohort. These butterflies had 6,672 raw rows (6,649 unique after excluding identifier fields and deduplicating identical content), including 4,695 complete plant-binomial rows and 3,898 resolved rows (3,890 unique). Neither catalogue identifier denotes a verified independent publication; the exported schema contains no explicit publication/reference field. Adding raw row count to both tied-rank controls retained partial ρ=+0.357 (2,000 species-bootstrap descriptive 95% interval +0.238 to +0.474; all 239 leave-one-out estimates +0.343 to +0.375). The four fixed alternative row-depth definitions gave partial ρ=+0.357 (deduplicated all rows), +0.348 (resolved rows), +0.348 (deduplicated resolved rows), and +0.348 (resolved rows per resolved host), with descriptive bootstrap intervals excluding zero in each case. The adequate-taxonomy and one-family subsets retained +0.380 (215 butterflies) and +0.464 (97 butterflies, controlling host species and raw rows), respectively. Raw row count was strongly correlated with resolved host count (ρ=+0.912), so these are correlated catalogue-depth controls, not an independent measurement of research effort. All predefined variants are retained in the accompanying source receipt.

These results show that controlling for catalogued host species, families and recorded row depth does not eliminate the association in Figure 4b. They do not demonstrate that the relationship is independent of true research effort, host breadth, phylogenetic non-independence, species selection, or the common HOSTS/WCVP data construction of both axes. The bootstrap is across the fixed selected species, and the tier-restricted permutation depends on exchangeability assumptions; neither is an independent ecological replication or a field-derived population confidence interval. Figure 4b therefore describes a mapped resource relationship, not as evidence of butterfly population growth or superior host quality.

## Appendix S5. Independently recorded butterflies and geographic transfer

### S5.1. Retrospective prediction of later records (155 species)

**Results**

Across 155 source-selected butterflies / 208,482 species×cell evaluations, adding known host-resource and host-family predictors lowered species-equal retrospective later GBIF record MSE by **+4.120%**, spatial block **+2.995%**, nested ≥100 km subset **+2.736%**; **119/155** and **115/155** benefited in temporal and spatial source tests. Contemporary plant maps may postdate the observation windows. These are **indexed records, not new colonies or biological abundances**.

Reproducibility and limitations. The 842,964-row prediction panel was recovered from a private archive and verified by checksum. Using the archived random-forest script and locked model settings, we independently refitted the temporal, spatial-block and unseen-species models (193,703 evaluable training rows and 208,482 test rows across 155 species before fold exclusion). For the temporal and spatial-block fits, all seven score groups, per-species errors and row-level predictions agreed with the archived outputs within prespecified tolerances (relative 1e−6, absolute 1e−9). The ≥100-km result used 102,728 predictions for 1,788 cells from the spatial fit. For unseen species the original and fresh relative MSE reductions were −0.20852% and −0.20893%, respectively. Both indicate higher error with host predictors. Numerical validation nevertheless failed for this holdout: 1,062 of 208,482 host-composition model predictions differed from their archived values beyond tolerance. These checks start from the original panel. Of four sensitivity fits, effort ≥20 and the targeted 250-km radius rerun agreed numerically; excluding the largest data provider and the 100-m precision filter retained positive gains but failed numerical equivalence (78 of 193 comparison groups in each). All 13 input files were recovered, but rebuilding the 842,964-row panel did not yield byte-identical Parquet: the only exact numeric discrepancy was distance_km. A targeted reconstruction of this distance column matched all 280,988 species–cell keys and 175 temporal comparisons when the other 95 archived columns were retained. The raw occurrence-to-day aggregation has not been rerun; full raw-to-model reproduction and public access to the original panel remain unresolved. The separate 224-species community refit covered 169 regions and 676 region–model score comparisons.

### S5.2. Recorded community composition and recording-intensity sensitivity

**Methods and results**

Independent GBIF DOI **10.15468/dl.umqmsf**: 14,637,152 original occurrence records; derived source 91,882 positive species×region×year records and 5,390 region×year effort cells. Mean indexed butterfly Jaccard was originally **0.114213** (2001–2010) and **0.167880** (2018–2024): raw delta +0.053667. Under 499 within-region recording-intensity-matched draw simulations, median raw delta was reduced to **+0.015064** (datasetKey×calendar-day proxy) or **+0.002914** (indexed record-count proxy). This does not quantify the proportion of biological homogenization due to observation effort.

The regional held-out MSE improvement from adding documented host flora after other introduced flora was +0.0913%, 224 butterflies/169 botanical regions/14,196 dependent pairs. After 29 repeated proxy-thinning draws of both prior and late butterfly similarity, median source MSE increments were +0.215% (conditional Monte Carlo 0.0659–0.4706) and +0.438% (0.1543–0.7362), positive 29/29 draws each. These 58 draws are not independent populations or true detection-corrected occupancy.

### S5.3. Recording-eligibility sensitivity

**Results**

**Table S7.** Recorded-community prediction across alternative recording-eligibility cohorts.

| Original effort cohort | Regions | Dependent region pairs | H-after-O recorded-community MSE improvement | Descriptive regional bootstrap central 95% |
| --- | --- | --- | --- | --- |
| Relaxed: at least 5 butterflies, 5 early dataset-days, 10 later | 194 | 18,721 | +0.3325% | +0.0712 to +0.5770% |
| Main: 5 butterflies, 10 early, 20 later | 169 | 14,196 | +0.0913% | −0.0915 to +0.2621% |
| Stringent: 10 butterflies, 20 early, 50 later | 133 | 8,778 | −0.0764% | −0.1554 to −0.0007% |

The stringent cohort reverses the sign. Region composition also changes with original selection threshold, so this comparison is an observation/region-source sensitivity, not a controlled observation-effort experiment. Later 169-region proxy matching does not erase sensitivity to the 194- and 133-region samples.

### S5.4. Geographic transfer and South American holdout

**Results**

Eight entire Level1 divisions removed from source fitting: original world-equal mean H-after-O MSE gain **+0.0971%**, proxy thinning medians **+0.6945% / +0.7793%**, with only 5–7 of eight world groups improving per draw. Level1=8 South America (20 regions, 3,170 heldout pairs) was adverse primarily for **2,980 South America–other-world pairs** rather than 190 within-SA pairs. Cross-world gain was **−4.036% / −2.697%** under the two original matching proxies. The unthinned source observed butterfly Jaccard mean was 0.06447 versus M2 prediction 0.08281 and after H added M3 prediction 0.08380. Registered H plants represented 24.1% of alien flora in South America and 23.6% elsewhere (regional medians). This indicates source prediction/calibration shift, not continent-specific biological rejection of introduced larval hosts.

The regional holdout sample comprised 169 WGSRPD3 regions and 14,196 dependent region pairs. Model M3 adds documented larval-host introduced-flora similarity (H) to M2, which already includes other introduced plants (O), geography, climate, earlier butterfly similarity, plant geography and observation-effort covariates. The unthinned improvement was +0.0913%; median improvements after the two proxy-thinning procedures were +0.215% and +0.438%.

For the eight held-out Level1 divisions, the world-equal improvements were +0.097% unthinned and +0.695% or +0.779% after thinning. Only five to seven divisions improved in each simulated draw, so the average does not imply universal transfer.

(c) South America: geographic transfer counterexample. The entirely unseen South American division comprised 20 regions, 190 within-division and 2,980 cross-world dyads. Median H-after-O MSE gains for cross-world dyads were −4.036% (datasetKey×calendar-day, conditional Monte Carlo 95% interval −7.268 to −1.473%) and −2.697% (record-count proxy, −8.043 to −0.137%). Within-South-America medians were −0.496% and −0.126%, respectively. All four estimates and their conditional Monte Carlo quantiles are shown in Figure S5.1c. This is a documented prediction-calibration failure, not evidence that South American butterflies cannot feed on introduced plants.

![Figure S5.1](figS51.png)

Figure S5.1. Recorded community prediction and geographic transfer (224 butterflies). Relative held-out MSE improvement from adding documented host-flora similarity to a model already including other introduced flora. (a) Regional holdouts: the unthinned point has the descriptive regional bootstrap interval (Table S7); the thinned estimates have conditional Monte Carlo 2.5–97.5% quantiles across 29 draws. (b) Entirely unseen TDWG Level1 divisions; no interval is available for the unthinned average. (c) Within-South-America and cross-world dyads for the South American holdout, with conditional Monte Carlo quantiles. Positive values indicate improved prediction. Conditional simulation quantiles are not population confidence intervals.

### S5.5. Recorded island butterflies and mapped introduced-only host opportunities

**Methods**

The original independent GBIF annual source (DOI 10.15468/dl.umqmsf; 91,882 positive species–region–year cells and 5,390 region–year effort cells, 13,587,452 mapped records) was joined by exact butterfly species and WGSRPD3 identity to the frozen 239 × 369 host-opportunity matrix. All 224 indexed species belonged to the fixed panel. All 126 island-region native/current counts matched the original 239-species matrix; the common 224-species island subset contained 4,695 native-host and 5,617 introduced-only opportunity units. For each species–island combination, recording was defined as at least one indexed year within the stated interval. The mutually exclusive mapped resource categories were native-host present (including possible additional introduced hosts), introduced-only host present, and no documented mapped host. No recording was not treated as absence. No mapped host does not exclude an unknown native host or a transient adult.

We fitted equal-region OLS models to the log((introduced-only recorded + 0.5)/(native-host recorded + 0.5)) among regions with recorded host-mapped butterflies. Covariates were standardized log mainland separation, log area, log2(port count + 1), absolute centroid latitude, TDWG Level1 effects and standardized log(datasetKey-calendar-day activity + 1) for the recording period. Minimum support was 40 contributing regions and five Level1 divisions. We report HC3 and Level1-clustered t inference (eight divisions in the supported models). Results are exploratory source-conditional associations, not independent replicated populations.

**Results**

Across the 126 regions, the 2001–2024 recorded species–region combinations comprised 916 native-host, 193 introduced-only and 62 no-mapped-host units. The introduced-only units occurred in 69 regions; 143 were recorded in at least two distinct years, across 57 regions. The introduced-only fraction among recorded combinations with mapped hosts was 17.4% (193/1,109), versus 54.5% among all 10,312 mapped opportunities for these butterflies. The corresponding fractions indexed were 3.44% of introduced-only opportunity units and 19.5% of native-host units, which must not be interpreted as detection-corrected occupancy or establishment probabilities. In the nearest and farthest distance quartiles, recorded introduced-only shares among host-mapped units were 11.4% (70/612) and 61.4% (27/44), respectively; all 32 near regions but only 12 of 32 far regions contributed to these fractions. In the original 169-region recording-eligibility cohort, only 24 island regions remained, with only one far-quartile island; thus we did not fit the multivariable isolation model to that cohort.

**Table S8.** Recorded introduced-only host share: isolation associations across periods and recording thresholds.

| Period and recording support | Regions | Isolation β | Eight-division clustered p |
| --- | --- | --- | --- |
| 2001–2024, all active | 84 | +0.439 | 0.0396 |
| 2001–2024, at least two recorded years per species–region | 74 | +0.499 | 0.0446 |
| 2018–2024, all active | 82 | +0.556 | 0.0270 |
| 2018–2024, at least two recorded years | 73 | +0.487 | 0.0589 |
| 2001–2024, regions with ≥20 activity days in 2001–2024 | 54 | +0.270 | 0.1093 |
| 2018–2024, same 54-region effort threshold | 54 | +0.337 | 0.0347 |
| 2018–2024, regions with ≥20 activity days specifically in 2018–2024 | 52 | +0.231 | 0.2239 |
| 2001–2024, excluding monarch | 83 | +0.412 | 0.0576 |
| 2018–2024, excluding monarch | 81 | +0.531 | 0.0363 |

For the 84-region model, the clustered 95% interval was +0.028 to +0.851; conditioning on potential resource share reduced the estimate to +0.285 (clustered p = 0.152). The observed-minus-potential log-odds gradient was uncertain (β +0.182, p = 0.324). These comparisons are consistent with resource geography but cannot establish a further independent isolation effect or a causal resource-mediated mechanism. Restricting to the most intensively recorded islands weakened inference. Neither repeated records nor historical resource maps establish breeding.

The monarch had 22 island regions with native-host opportunity and 42 additional regions with only mapped introduced hosts. Of the latter, 20 had indexed records and 18 had records in two or more years. These counts do not identify new colonizations; the contemporary host map does not date historical milkweed introductions. Data attribution and access limitations are detailed in Appendix S8.

## Appendix S6. Island-endemic butterflies and documented larval feeding

**Checklist coverage and reconstructed opportunities**

An independent Macaronesian checklist identified **41 living butterflies including 18 endemics** (13 Canary, three Madeira, two Azores). The fixed global 239 source selected **0/18** endemics (13/23 nonendemics). Among eight *nonrepresentative* endemic butterflies with accepted named foodplant species, four had added **global** potential host regions and four did not; none gained a newly host-present original *home archipelago* botanical region. Source-incomplete native hosts can change either side of the reconstructed envelope and make 4/8 nonrepresentative of all 18.

**Documented feeding observations**

Independently verified source photographs document *Pieris cheiranthi* larvae on **introduced Tropaeolum majus**, Tenerife (wild, 2015) and La Palma (garden, 2011). Independent native feeding in *Euchloe grancanariensis* on native *Hirschfeldia incana* further clarifies that documented exotic host use is not evidence that native food options have disappeared. The first case identifies **individual local larval feeding**, not new establishment or population benefit. Third-party photographs are not reproduced in journal files without rights clearance.

## Appendix S7. Local host records, port inventories and exploratory mechanisms

**Results: geographically localized host reports**

**Geographically localized HOSTS reports.** Island locality median documented introduced-host use **42.4%** (13 localities) versus mainland **9.1%** (15); exact FAO crop exclusion **20.9% vs 8.1%**. In 67 butterflies reported in both types, plant-family-matched and floristic supply-adjusted *extra island introduced-host use* was **−3.87 percentage points** (descriptive −15.85,+8.13): raw documented differences do not identify preference beyond plant availability.

**Results: port-associated records**

**Port-indexed butterfly/plant co-records.** In the original **17** sampled 5-km port buffers, **31** butterfly species were indexed, **16** with co-indexed compatible candidate-host plants from the same port neighborhood, yielding **8,525** original butterfly records and **59** species×port co-recording combinations (2001–2024 pooled; no synchronized larval feeding or matched nonport controls). Separately, the original 3,630-port roster had **1,596** complete botanical queries across 142 country codes, of which **1,150** contained a source-indexed regionally introduced host candidate. These comprised 687 candidate plant species and 25,742 plant–port pairs; 11,213 other plant–port candidate pairs had unresolved origin. In the original 17-island, 38-port early botanical sample (2001–2010), accepted-species-key matching identified 8/10 leading global host contributors. Their global fractional contribution was 20.7014% of the 14,553 new potential opportunities; the top ten together contributed 25.1011%. This is a contribution of the shared plant identities to the global resource map, not a fraction of global gains realized by the sampled ports. Buffers can overlap; the inventory lacks matched nonport controls and standardized recording effort. Regional contemporary plant origin does not establish site-specific origin in the earlier period. These are not vessel cargo observations or field interactions.

Table S9. Ten leading global host-plant contributors and port-buffer records in the 2001–2010 sample (38 buffers on 17 landmasses).

| Host plant | Share of global additions (%) | Port buffers with records |
| --- | --- | --- |
| *Medicago sativa* | 4.55 | 7 |
| *Poa pratensis* | 3.29 | 8 |
| *Oxalis corniculata* | 3.01 | 2 |
| *Cynodon dactylon* | 2.60 | 0 |
| *Dactylis glomerata* | 2.25 | 9 |
| *Zea mays* | 2.03 | 2 |
| *Avena sativa* | 2.03 | 5 |
| *Urtica urens* | 1.82 | 8 |
| *Oryza sativa* | 1.80 | 0 |
| *Digitaria sanguinalis* | 1.72 | 1 |

The port-buffer counts are retained from the earlier inventory display; they represent recorded presences, not standardized abundance or historical ship arrivals. Occurrence-level attribution and independent recreation of these counts remain outstanding.

Association between host-resource gains and port records (exploratory, same-source host catalogue). Among 225 butterfly species shared between the frozen 239-resource panel and original island-port inventory, 197 gained global mapped resource units and 159 had at least one island port with a recorded regionally introduced host candidate during 2001–2024. Spearman correlation between global added butterfly–region units and the number of sampled island ports with recorded introduced-host candidates was +0.65243; the partial rank association after host-count control was +0.56121 and after the source log1p host-count, native resource-unit and prior recorded-island controls was +0.54958. This retrospective comparison is not independent of the HOSTS/WCVP inputs used for the resource map, does not establish spatial or temporal precedence, and is subject to selection and heterogeneous botanical coverage. It is a taxon-resolved correspondence with the port-associated plant supply, not a port-specific enrichment or a ship-borne pathway.

**Exploratory nutrient and historical-arrival proxies**

**Port nutrient proxies (not measured terrestrial nitrogen).** Original near-port crop-N plant-composition proxy β **+0.08617** (descriptive +0.02375,+0.14859) and introduced-host fractional flora proxy β **+0.10144** (+0.03770,+0.16519) had positive source fits; stricter raster, botanical coverage and input-intensity sensitivities weakened some. These are botanical **record composition** correlations, not an identified sewage/fertilization-to-larval-host mechanism.

**First-index and physical route limits.** Among 317 apparently new 2001–2010-empty butterfly×island source cells, **145 already had pre-2001 GBIF records**. All 24/24 selected established-positive alien Lepidoptera cases had source-listed introduced hosts in recipient regions, but no controlled failed recipient comparison. Original NOAA pilot simulated 280 paths, 12 other-island coastal-nearshore directions and only two matched destinations, insufficient to identify a historical transport route. Colonisation, ships and individual propagule arrivals are not measured by the present paper.

## Appendix S8. Reproducibility and data access

Dataset versions and provenance. HOSTS DOI 10.5519/havt50xw (mirror commit 808e0b869f9ec1adf8efff87cf6a395adda103e0); WCVP v13 DOI 10.34885/0yex-xv26 (21 May 2024); WGSRPD3 geometry commit 52da7828aba9d461dd133c27b3bd7a4407161f54; GBIF DOI 10.15468/dl.umqmsf; versioned and checksum-verified 224-species annual regional positive and effort aggregates. The 155- and 224-species analyses share the GBIF download but not identical endpoints. The source-evidence register and archived analysis versions preserve exploratory, positive and adverse results.

Excluded exploratory results. The Lepidoptera-wide HOSTS network-degree analysis (8,909 plant species; degree-one consumers 26.2% introduced-range increase, ≥6 consumers 65.9%; stratified ρ ≈ +0.295) cannot distinguish catalogue coverage from interaction degree. A prior 51.4% first-index result lacks a fixed cohort, and a 32-species comparison (66/115 versus 60/109) has unresolved denominator scope. An earlier 125-island/2,517-butterfly isolation analysis and preliminary taxonomic-rank input files were unavailable for direct reconstruction. These results are excluded from the paper’s confirmatory claims.

Analysis reconstruction. The fixed 239-species resource accounting and the isolation regressions were replayed on the archived derived regional and species-level tables. The 126-region/239-species counts and 239 global species totals matched the fixed tables. The spatial-filtering analysis was previously reproduced with the pinned botanical geometry; its published estimates were checked against the archived results in the present figure audit, rather than independently refitted from geometry. The earlier exploratory sample is not represented as replicated.

Port and marine analyses. The crossed port–water models were refitted from fixed dyad and area tables; all 1,726 report fields agreed within 7.3e−14. The exact model outputs and code versions are retained in the internal reproducibility record.

The plant-recording-effort and host-catalogue-depth sensitivity outputs and analysis scripts are identified in the reproducibility record.

The port inventory and taxon-association outputs are identified in the reproducibility record. Occurrence-level attribution for the historical and port-buffer acquisitions remains incomplete.

Figure inputs. The figure values and row-level plotting arrays were checked against the archived fixed data and model results. The 239 species pairs reproduce Spearman ρ = 0.5788964068803357 and 5,741 island additions; the 84 recorded regions sum to 916 native-host and 193 introduced-only combinations. The three-distance predictions are verified model outputs shown as point estimates, not newly fitted curves. The recorded-share curves use the archived model coefficients, also verified by refitting the archived derived regional inputs. All display inputs and plotting code are retained for reproducibility; release and reviewer access remain to be resolved. This check is not a reconstruction of all raw occurrence data.

**Data-source limits.** HOSTS, WCVP v13, GSHHG and the original 224-butterfly GBIF download have fixed source identifiers. The separate historical first-index and port-buffer GBIF queries did not preserve original occurrence-level contributing publisher identifiers. Reproduction therefore depends on archived query receipts, permitted derived tables and source-specific access terms; these source limits are distinct from ecological evidence. The official World Port Index citation describes the 2019 27th edition. The source GeoJSON matches a comparison mirror to six decimal places, but its original download, conversion and edition remain unverified.

### Supporting Information references and source citations

Brummitt, R. K. (2001). *World geographic scheme for recording plant distributions* (2nd ed.). Pittsburgh, PA: Hunt Institute for Botanical Documentation, Carnegie Mellon University. Retrieved from https://rs.tdwg.org/wgsrpd/doc/data/

ESR. (2009). *OSCAR third degree resolution ocean surface currents* (Version 1) [Data set]. NASA Physical Oceanography Distributed Active Archive Center. https://doi.org/10.5067/OSCAR-03D01

Fall, P. L., & Drezner, T. D. (2011). Plant dispersal, introduced species, and vegetation change in the South Pacific Kingdom of Tonga. *Pacific Science, 65*(2), 143–156. https://doi.org/10.2984/65.2.143

GBIF.org. (2026, October 29). *GBIF occurrence download* [Data set]. https://doi.org/10.15468/dl.umqmsf

Govaerts, R. (2024). *World Checklist of Vascular Plants (WCVP)* (Version 13, 21 May 2024) [Data set]. Board of Trustees of the Royal Botanic Gardens, Kew. https://doi.org/10.34885/0yex-xv26

National Geospatial-Intelligence Agency. (2019). *World port index* (27th ed., Publication 150). Springfield, VA: Author. Retrieved from https://msi.nga.mil/api/publications/download?key=16694622%2FSFH00000%2FPub150bk.pdf&lang=en

Robinson, G. S., Ackery, P. R., Kitching, I., Beccaloni, G. W., & Hernández, L. M. (2023). *HOSTS - a database of the world's lepidopteran hostplants* [Data set]. Natural History Museum. https://doi.org/10.5519/havt50xw

Vargas, P., Heleno, R., & Costa, J. (2023). *EuDiS - A comprehensive database of the seed dispersal syndromes of the European flora* (Version 5) [Data set]. Figshare. https://doi.org/10.6084/m9.figshare.22251028.v5

Wessel, P., & Smith, W. H. F. (2017). *GSHHG: Global self-consistent hierarchical high-resolution geography* (Version 2.3.7) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7007502

The main text contains Figures 1–5. Supporting figures are numbered by appendix; Tables S1–S9 retain the numerical results and sensitivity analyses.
