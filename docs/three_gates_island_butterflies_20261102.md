# Three gates to island butterfly colonization: hosts × climate × access

Status: **exploratory framework and empirically audited source feasibility**, 2026-11-02. Analysis branch only; NOT integrated into frozen GEB submission.

## Ecological question

**When butterfly colonization of an island is first documented, which independently varying constraints—larval hosts, climate, or arrival opportunity—were compatible with the observed change, and does human maritime access compensate for low natural dispersal potential?**

Three gates:
- **Resource (H):** established suitable larval host plants, including introduced plants.
- **Climate (C):** season-specific development and overwintering suitability, not just mean temperature.
- **Access (A):** actual chance of propagules reaching the island through natural flight/winds or human transport (plant stock, vessels, airports).

Colonization requires available hosts, adequate climate, arrival, and subsequent reproduction. These constraints are **not necessarily binary**, and host/climate/access can co-vary; source records do not identify which was **last** to appear without independently dated evidence. The biologically strongest inference is therefore a **comparative hazard/limitation model**, not deterministic last-gate classification.

### Source-validated initial findings

| Component | Source-audited state | Limits |
|---|---|---|
| Fixed GEB resource panel | 239 butterfly species, 26,530 native and 41,083 contemporary species×WGSRPD3 resources (+54.9%) | Potential host opportunity, not butterfly establishment |
| Q1 islands | 126 strict island WGSRPD3 regions; 4,829→10,570 (+118.9%) versus mainland 21,215→29,901 (+40.9%) | Area and biogeography confounding |
| New island host opportunities | 5,741 butterfly×island cells; 11,573 butterfly–introduced-host–island edges; 564 host taxa | Same plant/host source used for all resource construction |
| Q4 deduplicated plant dates population | 5,611 introduced-host×WGSRPD3 island pairs across 123 islands (includes 66 Canary CNY host–island pairs) | Modern WCVP alien layer, not historical introduction dates |
| GBIF pilot | Fixed original 500 response-blind pairs with **219 earliest indexed plant years** (median 1975, range 1805–2026), plus 3 geographically selected Canary CNY pairs each with a dated GBIF plant index (1964, 2020, 2019), giving **222/503 dated**. Original 500 retain **163** API-success nonrecords; **94** failed queries (69 HTTP429; 25 HTTP400); **24** other taxon/geometry exclusions | First indexed occurrence is only an upper limit on true date of plant presence. 69 failed calls require throttling; do not count as plant absences |
| Literature event seed | 47 candidate historical events, **13 butterfly species**, 29 physical islands/localities, 25 mapped species×WGSRPD3 region cells, 26 events from 1981 onwards, 20 with a subsequent established/breeding upper year | Not 47 independent establishment years; 3 events without WGSRPD mapping |
| Literature source bridge, fixed pilot | 25/47 lead events are at corrected modern host resource cells; 10 had dated suitable hosts in the independent response-blind 500-pair pilot, 8 showed plant index before butterfly documentation (3 species×WGSRPD3 cells), and 1 had all proxy fields | Historical baseline only; record order is not arrival or local larval feeding |
| LepTraits | 239/239 species wing measurement, 129 diapause/overwintering stage, 148 oviposition style; 239 disturbance affinity | Wingspan is only a noisy flight proxy; no direct migration variable in this frozen CSV |
| Butterfly presence | Existing original 224-species / 14,425,536-record island/year file was previously source-verified, but source Parquet is private and not mounted in this cloud workflow; public 31-species GBIF DOI provides 13 recorded new-resource cells among 920 potential new-resource cells in the corrected island definition | Cannot pretend to possess the full original 224-species time series |
| NASA POWER thermal pilot | 15/15 literature WGSRPD3 regions with complete **1981–2025** annual T2M series from monthly NASA POWER; 21 literature cases have sufficient ten-year pre-first-record temperature history. Earlier 5-region 1981–2020 pilot is retained only as a source test | Not a species-specific thermal threshold. Historical pre-1981 establishment dates are outside POWER chronology |
| Modern World Port Index | 3,630 published WPI points, 3,274 uniquely mapped to WGSRPD3 (direct/within 5 km). 982 mapped to 91/126 corrected strict island regions; Canary Islands were initially omitted and contributed four modern indexed ports | Public mirror last edited 2023; count is neither historical ship traffic nor date of propagule introduction |

| Targeted host chronology | After a literature-case-selected search, **67/84 introduced-host×region pairs** have an indexed first plant year, and **21/47 literature lead rows** have at least one region-level compatible host record older than the butterfly first report; only **11 distinct butterfly×WGSRPD3 cells**. Full observed proxies (host lead, climate context, ports, wings) now overlap in **11 lead rows / 2 species / 7 region cells** | Selected on butterfly detections, regional botanical index dates are not actual naturalization or physical-island host presence |
| GBIF-wing-port exploratory screen | In the corrected nonrandom 31-butterfly DOI, **13/920** new-resource island cells had butterfly observations; all 13 indexed cells occur in regions with at least one modern WPI port; 736/920 source opportunities have a recorded WPI port. Permutation within species×Level1: port score p(two-sided)=0.246, wing×port score p=0.218 | Strong GBIF reporting selection, only 13 positive observation cells; no vessel-transport effect or colonization regression identified |

Primary data receipts: `results/q4-host-first-year-v02/`, `results/q4-gbif-first-year-pilot/`, `results/q4-canary-plant-firstyear-v02/`, `results/colonial-host-native-source-q2-v02/`, `results/q4-literature-hostyear/`, `results/three-gate-island-integration-v02/`, `results/q3-observed-butterfly-island-opportunity-v02/`, `results/three-gate-nasa-power/`, `results/three-gate-wpi-port-v02/`.

### Geographic source repair and study readiness

- **2026-11-02 correction:** WGSRPD3 CNY Canary Is. was mistakenly classified as mainland in an early manual 125-island list. The official area-name audit confirms it is a strict island; RHO Rhode I. refers to the US state, not an oceanic island. Q1 is now **126 island / 237 mainland** regions. The original 239-species GLOBAL count (26,530→41,083) remains exactly the same; independent rebuild of the Q2 introduced-host edges yields **5,741** island butterfly-resource opportunities and **11,573** distinct butterfly×plant×island edges. Older 125-island Q2 colonial-network tests remain **legacy exploratory sensitivity**. The 126-island Q1 area-covariate result has been rebuilt independently (island coefficient +0.0672; within Level1×area-quarter label permutation p=0.0002, descriptive).
- LepTraits wingspan columns WS_* are recorded in **centimeters** (Shirey et al. 2022); the case matrix provides separate cm and converted mm columns and never merges raw cm with the GEB pre-scaled wing proxy.
- The fixed 500-plant-pair pilot supported complete proxy coverage for **1** of 47 literature leads. A later case-selected GBIF enrichment (84 candidate introduced host×region pairs, 67 dated; 6 query failures) increased documented pre-butterfly host-index overlaps to **21 of 47 literature events**, representing **11 distinct butterfly×WGSRPD3 cells**. It increased complete available source-proxy coverage to **11 events**, but these represent only **2 butterfly species and 7 distinct butterfly×WGSRPD3 cells**. The last biological gate remains unidentified in **all 47**. **Do not fit a wing×port colonization regression from these positive-only selected cases.**

### Core tests, when observation source is actually available

1. **Denominator:** island×species×calendar year eligible risk sets, including never indexed species, with independently classified butterfly nonnative/native status and no historical records before left-truncation. Treat missing GBIF record as *not indexed*, not absent. Validate 224-species source and source taxon keys before joining to 239 HOSTS panel.
2. **Observation model:** estimate butterfly detectability from native-region effort, focal island total butterfly reporting, sampling source/provider and year; account for iNaturalist-era recording jump. Repeated records, duplicate datasets, specimens and observation days should not create independent colonization events.
3. **Host calendar:** plant year-facet history with all-island GBIF botanical effort and exact species keys, preferably dated naturalization or nursery-import first reports. GBIF first indexed plant date is a right-end observation bound (`true_arrival <= indexed_year`), not arrival time. Host availability prior to butterfly source event is uncertain unless independently documented.
4. **Climate calendar:** use NASA POWER 1981+ island time series to compute species-dependent seasonal cold limitation and degree-day thresholds calibrated without island outcome leakage; use independent historical reconstructions for earlier events or exclude those events from climate last-gate classification.
5. **Arrival/access:** compare port exposure to self-dispersal proxies through **wing-size × port**, **wing-size × mainland sea distance**, **development-stage transportability × port**, with island area, host richness, biogeography, geography and sampling controlled. Overwintering stage is not automatically the transported stage; port count is not a vessel route. Air transportation, horticultural imports and winds are alternatives.

Suggested comparisons: pre-registered baseline with island/species effects and reporting effort; host-only, climate-only, access-only; all additive; interaction model; leave-island and leave-species-out checks. Predict a **weakening of isolation effect with large wing size**, and **stronger positive port association for small-wing species** only conditionally on reporting effort and resource/climate. Inferences about ships need independent route/plant-shipment data; a port interaction alone is not direct proof of vessel-mediated introduction.

### Boundary of inference

Avoid statements that "98.6% of potential resources remain unused" (corrected Q3: 13/920 GBIF-recorded overlap from a nonrandom 31-species partial sample). A host recorded before a butterfly does **not** imply host arrival before actual butterfly arrival. Both observations are right-censored proxies. Port number and adult wing size are **not sufficient** to deduce arrival mode. A novel host opportunity for a native butterfly is not inherently an alien-butterfly colonization.

### Literature leads (frozen CSV)

`data/frozen/exploratory_island_butterfly_invasion_literature_seed_20261102.csv`: sources linked per event. Examples include *Pieris rapae* 1930 New Zealand; *Cacyreus marshalli* Mallorca/Canaries and Madeira (Madeira establishment uncertain); *Papilio demoleus* Jamaica 2006→established by 2008, Mahé Seychelles 2016→established by 2020; *Brephidium exilis* Oahu 1978→established by 1980; Hawaiian skipper and monarch first records. Dominican Republic/Puerto Rico *P. demoleus* historical first-year sources conflict and remain flagged, not resolved by averaging.

Study documents: https://www.nature.com/articles/s41597-022-01473-5 (LepTraits), https://pubmed.ncbi.nlm.nih.gov/21988561/ (wingspan meta-analysis), https://power.larc.nasa.gov/docs/services/api/temporal/monthly/ (NASA POWER), https://www.redalyc.org/pdf/455/45528755006.pdf (Canary colonizations), https://gd.eppo.int/taxon/PAPIDE/distribution/JM (Jamaica), https://www.biotaxa.org/em/article/view/em.2020.28.7 (Seychelles), https://teara.govt.nz/en/1966/butterfly-white (NZ).

**Readiness:** source feasibility established; biological last-gate attribution and wing×port establishment effect **not yet identified**. The immediate bottleneck is source-backed island/year observation risk sets and independently known establishment/native status, not lack of potential-resource or trait data.


### New exploratory evidence: indexed host chronology and the wing×port screen (2026-11-02)

The case-targeted botanical search extends a **different** sample from the original response-blind 500 plant×island pilot. It finds compatible *GBIF indexed* plants before the butterfly literature report for 21 of 47 leads (11 distinct butterfly×WGSRPD3 cells). This is a **documentary chronology**, not true first introduction or butterfly arrival. Source files are `results/q4-targeted-literature-hostyear/` and `results/three-gate-integrated-targeted-v03/`. The latter source-backed matrix contains 11 rows with all *available proxy fields*, but only 2 species (Leptotes pirithous and Papilio demoleus), making causal decomposition across three gates nonidentifiable.

The nonrandom 31-species GBIF download screen includes 920 source-defined newly enabled butterfly×island opportunity cells and only 13 with an indexed butterfly record. All 13 happen to be in 2023-WPI-port-bearing WGSRPD3 regions, but the species- and region-matched port score is not distinguishable from the conditional randomized distribution (two-sided p=0.2458); neither is wing×port (p=0.2181). Requiring an additional other-butterfly-indexed observation stratum leaves port p=0.574 and wing×port p=0.267. The observation activity surrogate itself is related to where the 13 records occur; it is not an independent ecological supply mechanism. Source: `results/q3-wing-port-index-screen-v02/wing_port_observation_screen_audit.json`.

**Primary empirical claim maintained:** human-mediated host plant redistribution expands **potential** butterfly resource geography disproportionately on islands. **Hypothesis still open:** if matched source-year risk sets become available, assess which among resource, species-relevant climate and access histories best predicts *subsequent indexed butterflies*, and separately whether the result generalizes to independently documented naturalization/establishment. No original positive results are invalidated by a failure to identify a causal arrival mechanism.


### Botanical effort-controlled port–resource evidence

A new complete-case sensitivity (98/126 islands with GBIF indexed all-plant reporting effort) independently queries WCVP introduced species, HOSTS-selected introduced host taxa, WPI ports and GBIF botanical collection intensity. All source files are audit-pinned. **Modern ports are positively associated with the quantity of introduced flora, even after conditioning on native flora richness, region area, broad Level1 biogeography and indexed botanical records**:

- WCVP all introduced accepted plant species: β=+0.2445 in log(1+species) per log2(1+WPI ports), corresponding to +27.7%, 4,999 Level1×area-half permutations, two-sided centered p=0.0012.
- HOSTS-selected introduced larval-host taxa: β=+0.1639, +17.8%, two-sided centered p=0.0110; the native butterfly-resource baseline is controlled.
- WCVP introduced/native species log-ratio still has a positive port association (β=+0.2429; p=0.0036).
- WPI itself also positively correlates with indexed botanical reporting effort: β=+0.5239, +68.9%, p=0.0002; this is a source of confounding, not an independent access mechanism.
- Sensitivities: restricting to indexed-port islands only (n=69) β=+0.2354, p=0.0222; botanical effort ≥1,000 records (n=71) β=+0.2193, p=0.0466; adjusting pre-1900 botanical effort instead of all-years β=+0.2225, p=0.0014.

The 28 unqueried/non-comparable islands have median area ≈25,120 km² versus ≈956 km² among the 98 included. This excluded cohort includes Japan, Great Britain and Cuba; geographic missingness and human-population/port/recording confounding prevent a causal or all-islands generalization. GBIF reporting itself may be affected by historic transport and settlement, so recording-effort adjustment is **sensitivity**, not causal control.

**The three-gate inference has therefore advanced most strongly in the *upstream plant-supply* link, not the hypothesized port-mediated butterfly arrival, wing-size interaction or colonization-delay mechanism.** These remain distinct open tests and require genuine island×species×year observation/risk sets or independent establishment dates.

Full source and audit: [effort-adjusted WPI/flora model](../results/three-gate-port-botanical-effort/botanical_effort_adjusted_port_flora_audit.json), [case-island model input](../results/three-gate-port-botanical-effort/botanical_effort_flora_island_inputs.csv), and [ecological interpretation and selection caveat](./three_gate_port_resource_source_results_20261102.md).
