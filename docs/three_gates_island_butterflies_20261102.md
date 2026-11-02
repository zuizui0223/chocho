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
| Literature source bridge | 25/47 events fall in corrected Q2 resource-added butterfly×WGSRPD3 cells; 10 events had a dated host in the fixed original 500-pair pilot, 8 showed host index before the butterfly literature record (3 unique species×WGSRPD3 cells). Only **1/47 cases** presently has all recorded host lead, NASA preceding-decade T2M, modern WPI port count, and LepTraits wing measure | Same archipelago's plant index can be on a different physical island; no causal lead or local feeding shown |
| LepTraits | 239/239 species wing measurement, 129 diapause/overwintering stage, 148 oviposition style; 239 disturbance affinity | Wingspan is only a noisy flight proxy; no direct migration variable in this frozen CSV |
| Butterfly presence | Existing original 224-species / 14,425,536-record island/year file was previously source-verified, but source Parquet is private and not mounted in this cloud workflow; public 31-species GBIF DOI provides 13 recorded new-resource cells among 920 potential new-resource cells in the corrected island definition | Cannot pretend to possess the full original 224-species time series |
| NASA POWER thermal pilot | 15/15 literature WGSRPD3 regions with complete **1981–2025** annual T2M series from monthly NASA POWER; 21 literature cases have sufficient ten-year pre-first-record temperature history. Earlier 5-region 1981–2020 pilot is retained only as a source test | Not a species-specific thermal threshold. Historical pre-1981 establishment dates are outside POWER chronology |
| Modern World Port Index | 3,630 published WPI points, 3,274 uniquely mapped to WGSRPD3 (direct/within 5 km). 982 mapped to 91/126 corrected strict island regions; Canary Islands were initially omitted and contributed four modern indexed ports | Public mirror last edited 2023; count is neither historical ship traffic nor date of propagule introduction |

Primary data receipts: `results/q4-host-first-year-v02/`, `results/q4-gbif-first-year-pilot/`, `results/q4-canary-plant-firstyear-v02/`, `results/colonial-host-native-source-q2-v02/`, `results/q4-literature-hostyear/`, `results/three-gate-island-integration-v02/`, `results/q3-observed-butterfly-island-opportunity-v02/`, `results/three-gate-nasa-power/`, `results/three-gate-wpi-port-v02/`.

### Geographic source repair and study readiness

- **2026-11-02 correction:** WGSRPD3 CNY Canary Is. was mistakenly classified as mainland in an early manual 125-island list. The official area-name audit confirms it is a strict island; RHO Rhode I. refers to the US state, not an oceanic island. Q1 is now **126 island / 237 mainland** regions. The original 239-species GLOBAL count (26,530→41,083) remains exactly the same; independent rebuild of the Q2 introduced-host edges yields **5,741** island butterfly-resource opportunities and **11,573** distinct butterfly×plant×island edges. Older 125-island Q2 colonial-network tests and v0.1 area-covariate results are **legacy exploratory sensitivity**, not contemporary primary 126-island estimates until their source geometry is rebuilt.
- LepTraits wingspan columns WS_* are recorded in **centimeters** (Shirey et al. 2022); the case matrix provides separate cm and converted mm columns and never merges raw cm with the GEB pre-scaled wing proxy.
- Three-gate evidence matrix has **47 literature lead rows**, **25** with corrected Q2 introduced host opportunity, **21** with pre-first-record decade NASA temperature, **44** with mapped modern WPI port counts, **47** with a LepTraits wing measure, but only **1** with all available proxy fields including indexed host lead. That sole candidate is *Papilio demoleus* Puerto Rico; its butterfly literature first date is disputed and none of the three bottlenecks has an independently determined onset. **Do not fit a wing×port colonization regression from these 47 positive-only selected cases.**

### Core tests, when observation source is actually available

1. **Denominator:** island×species×calendar year eligible risk sets, including never indexed species, with independently classified butterfly nonnative/native status and no historical records before left-truncation. Treat missing GBIF record as *not indexed*, not absent. Validate 224-species source and source taxon keys before joining to 239 HOSTS panel.
2. **Observation model:** estimate butterfly detectability from native-region effort, focal island total butterfly reporting, sampling source/provider and year; account for iNaturalist-era recording jump. Repeated records, duplicate datasets, specimens and observation days should not create independent colonization events.
3. **Host calendar:** plant year-facet history with all-island GBIF botanical effort and exact species keys, preferably dated naturalization or nursery-import first reports. GBIF first indexed plant date is a right-end observation bound (`true_arrival <= indexed_year`), not arrival time. Host availability prior to butterfly source event is uncertain unless independently documented.
4. **Climate calendar:** use NASA POWER 1981+ island time series to compute species-dependent seasonal cold limitation and degree-day thresholds calibrated without island outcome leakage; use independent historical reconstructions for earlier events or exclude those events from climate last-gate classification.
5. **Arrival/access:** compare port exposure to self-dispersal proxies through **wing-size × port**, **wing-size × mainland sea distance**, **development-stage transportability × port**, with island area, host richness, biogeography, geography and sampling controlled. Overwintering stage is not automatically the transported stage; port count is not a vessel route. Air transportation, horticultural imports and winds are alternatives.

Suggested comparisons: pre-registered baseline with island/species effects and reporting effort; host-only, climate-only, access-only; all additive; interaction model; leave-island and leave-species-out checks. Predict a **weakening of isolation effect with large wing size**, and **stronger positive port association for small-wing species** only conditionally on reporting effort and resource/climate. Inferences about ships need independent route/plant-shipment data; a port interaction alone is not direct proof of vessel-mediated introduction.

### Boundary of inference

Avoid statements that "98.6% of potential resources remain unused" (13/912 was GBIF-recorded overlap from 31-species partial sample). A host recorded before a butterfly does **not** imply host arrival before actual butterfly arrival. Both observations are right-censored proxies. Port number and adult wing size are **not sufficient** to deduce arrival mode. A novel host opportunity for a native butterfly is not inherently an alien-butterfly colonization.

### Literature leads (frozen CSV)

`data/frozen/exploratory_island_butterfly_invasion_literature_seed_20261102.csv`: sources linked per event. Examples include *Pieris rapae* 1930 New Zealand; *Cacyreus marshalli* Mallorca/Canaries and Madeira (Madeira establishment uncertain); *Papilio demoleus* Jamaica 2006→established by 2008, Mahé Seychelles 2016→established by 2020; *Brephidium exilis* Oahu 1978→established by 1980; Hawaiian skipper and monarch first records. Dominican Republic/Puerto Rico *P. demoleus* historical first-year sources conflict and remain flagged, not resolved by averaging.

Study documents: https://www.nature.com/articles/s41597-022-01473-5 (LepTraits), https://pubmed.ncbi.nlm.nih.gov/21988561/ (wingspan meta-analysis), https://power.larc.nasa.gov/docs/services/api/temporal/monthly/ (NASA POWER), https://www.redalyc.org/pdf/455/45528755006.pdf (Canary colonizations), https://gd.eppo.int/taxon/PAPIDE/distribution/JM (Jamaica), https://www.biotaxa.org/em/article/view/em.2020.28.7 (Seychelles), https://teara.govt.nz/en/1966/butterfly-white (NZ).

**Readiness:** source feasibility established; biological last-gate attribution and wing×port establishment effect **not yet identified**. The immediate bottleneck is source-backed island/year observation risk sets and independently known establishment/native status, not lack of potential-resource or trait data.
