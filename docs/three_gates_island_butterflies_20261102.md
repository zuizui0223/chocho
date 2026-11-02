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
| Q1 islands | 125 strict island WGSRPD3 regions; 4,743→10,417 (+119.6%) versus mainland 21,301→30,054 (+41.1%) | Area and biogeography confounding |
| New island host opportunities | 5,674 butterfly×island cells; 11,455 butterfly–introduced-host–island edges; 560 host taxa | Same plant/host source used for all resource construction |
| Q4 deduplicated plant dates population | 5,545 introduced-host×WGSRPD3 island pairs across 122 islands | Modern WCVP alien layer, not historical introduction dates |
| GBIF pilot | Fixed, response-blind 500 pairs; **219 earliest indexed plant years** (median 1975, range 1805–2026); **163** API success but no indexed observation; **94** failed queries (69 HTTP429;25 HTTP400); **24** other taxon/geometry exclusions | First indexed occurrence is only an upper limit on true date of plant presence. 69 failed calls require throttling; do not count as plant absences |
| Literature event seed | 47 candidate historical events, **13 butterfly species**, 29 physical islands/localities, 25 mapped species×WGSRPD3 region cells, 26 events from 1981 onwards, 20 with a subsequent established/breeding upper year | Not 47 independent establishment years; 3 events without WGSRPD mapping |
| Literature source bridge | 21/47 events fall in Q2 resource-added butterfly×WGSRPD3 cells; 10 events had dated introduced host in the fixed pilot, 8 showed host index before butterfly literature record, but these reduce to only **3 independent species×WGSRPD3 cells** | Same archipelago's plant index can be on a different physical island; no causal lead or local feeding shown |
| LepTraits | 239/239 species wing measurement, 129 diapause/overwintering stage, 148 oviposition style; 239 disturbance affinity | Wingspan is only a noisy flight proxy; no direct migration variable in this frozen CSV |
| Butterfly presence | Existing original 224-species / 14,425,536-record island/year file was previously source-verified, but source Parquet is private and not mounted in this cloud workflow; public 31-species GBIF DOI provides only 13 recorded new-resource cells | Cannot pretend to possess the full original 224-species time series |
| NASA POWER thermal pilot | 5 complete region-point **1981–2020** T2M histories (JAM, PUE, NZN, CUB, SEY); 2001–2020 means rise relative to 1981–2000 by +0.18 to +0.55 °C | Not a species-specific thermal threshold. Historical pre-1981 establishment dates are outside POWER chronology |
| Modern World Port Index | 3,630 published WPI points, 3,274 uniquely mapped to WGSRPD3 (direct/within 5 km). 978 mapped to 90/125 strict island regions | Public mirror last edited 2023; count is neither historical ship traffic nor date of propagule introduction |

Primary data receipts: `results/q4-host-first-year/`, `results/q4-gbif-first-year-pilot/`, `results/q4-literature-hostyear/`, `results/three-gate-island-butterfly/`, `results/three-gate-nasa-power/`, `results/three-gate-wpi-port/`.

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
