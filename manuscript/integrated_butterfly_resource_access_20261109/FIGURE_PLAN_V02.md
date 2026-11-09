# Figure system v0.2 — one butterfly paper, islands and ports in the foreground
**Canonical manuscript:** ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V02.md | 2026-11-09
**Status:** source-designed figure programme; existing separate GEB and JBI figures are reproducible references, not yet a completed new set.

## What each figure must prove — without borrowing unmeasured causal arrows

| Figure | Scientific question | Primary source grain | Visual evidence | Adverse / scope control |
|---|---|---|---|---|
| **1. Island-port resource gateway?** | Where do introduced plants enter the story of butterflies, and which relationships are measured? | 239 butterfly×WGSRPD3; 98 island botanical regions; 225×51 record sample | Butterfly → fixed compatible plant list → native and introduced host geography; island amplification; broad port flora; actual butterfly index records | Dashed host-to-record associations, no ship-transport or butterfly-colonisation arrow |
| **2. Island amplification of resources** | Why are islands disproportionately affected by introduced host geography? | 126 insular, 237 mainland WGSRPD3; 15 matched Level1×area strata; 239 butterfly host envelopes | Native→contemporary +118.9% island vs +40.9% mainland, matched +10.7 pp; global +54.9% context; known hosts 38/670 deliver half, 82/97 one-family butterflies also gain | Source unit is potential butterfly×region opportunities, not butterfly population or dispersal |
| **3. Ports, alien plants and named larval hosts** | Do port-rich island regions concentrate butterfly-compatible alien plants separately from general introduced flora? | 1596 source-complete ports (142 country codes) and **distinct** 98 WGSRPD3 botanical-effort island regions; 17 landmasses/38 port buffers with 2001–2010 indexed accepted plant names | 1596/3630 complete ports, 1150 with introduced-host candidates, 687 plant species and 25742 candidate host×port records; three distinct beta estimates +0.2445, +0.1639, −0.0415 after whole-flora control; 8/10 major global host-contributor species seen pre2011 in at least 1 studied port buffer | 1596 complete port locations NOT 1596 matched independent island controls; 8 species NOT shipped/verified at that time as non-native |
| **4. Butterflies recorded where host opportunities are indexed** | Which butterflies and host matches are observable? | 225 candidate butterflies species bridge, 51 landmass 1275 dependent pairs, small 17-port 5km source; separate 155 butterfly forecasting population | Species added potential resource vs number of indexed candidate-host ports ρ +.652 (partial+.550); candidate host- and later butterfly-composition ρ +.300 (adjusted +.215), earlier butterfly +.418; 208/1163 = 17.9% compatible host coindex vs 11.9% host-family/degree-only null; name 16 butterfly species with 59 same-port butterfly×compatible-host candidate coindexed combinations; separate 155 species forecasting host increment +4.12% temporal and +3.00% spatial | 225 port coindex and 155 forecast population are not one matched trial; host records pooled 2001–2024 may postdate butterfly; 1275 pairs not independent |
| **5. What can the island–port evidence actually identify?** | Do earlier compatible plant records precede butterfly indexing; do access/traits support beyond-geography mechanism? | 113 historical butterfly first indexes/17 islands, 225×51 new public GBIF source, 199 traits and 111 indirect migrant candidates | Earlier compatible plant record 34/113; basic family/degree-null 25.2% vs observed30.1%, but geographically preserved null 30.26%; independent 225×51 butterfly first-index records 221 early / 96 late; within-region source proximity p=.00033 and residual sea access p=.00233; geographic holdout adds sea → Brier −6.48%; indirect migrant 49-island heavy-port exclusion ocean residual p=.229 | Earlier compatible host does not boost five-year reindex (26/34 vs 70/79); fixed prior island not true origin; ocean graph no relay, wing/migrant status no robust independent effect |

## Butterfly names and original plant names to show

- Plant species with globally disproportionate contribution and verified 2001–2010 acceptance at some port buffers: *Medicago sativa* (7 port buffers/3 islands), *Poa pratensis* (8/4), *Dactylis glomerata* (9/5), and other taxonomically verified contributors. Avoid depicting *Cynodon dactylon* as historically indexed in this specific 38-port source; its query returned zero indices here.
- Observational butterfly × port-buffer catalogue matches in the small original 17-port source: *Pieris rapae*, *Pieris brassicae*, *Vanessa cardui*, *Anthocharis cardamines*; these are not proven caterpillar feeding or the actual body of a migrant.
- Species-specific 225×51 first-index event profiles: *Zizula hylax* (5 new landmass indexes in 2018–2022), *Nymphalis antiopa* (4), *Danaus plexippus* (3) and *Vanessa cardui* (3). These are first GBIF indexing events, NOT their island arrival or ship-travel distances.
- Contrasting named plants and butterflies require explicit source/sample labels, because the 2001–2010 38-port plant subset, 2001–2024 17-port butterfly co-index, 51-landmass original pair and 225×51 public temporal panel are not one identical joint field survey.

## Source files for plots (all locked and already committed)

**Fig 2** 
- results/island-mainland-resource-q1/island_mainland_summary.json
- results/island-mainland-resource-area-v02/area_adjusted_summary.json
- benchmarks/exploratory/butterfly_anthropogenic_resource_expansion_result_v0.1.json
- provenance/reviewer_defenses/results/butterfly_host_contribution_concentration_v0.1.json

**Fig 3**
- provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_validation.json
- results/three-gate-port-botanical-effort/botanical_effort_adjusted_port_flora_audit.json
- results/three-gate-port-flora-host-conversion/island_port_flora_host_conversion_audit.json
- provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_top10_pre2011_port_plant_overlap.json

**Fig 4**
- provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_to_port_species_bridge.json
- provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_audit.json
- provenance/reviewer_defenses/results/island_shared_host_match_20261030/island_shared_host_match_audit.json
- provenance/reviewer_defenses/results/island_host_null_20261030/island_host_family_null_audit.json
- provenance/reviewer_defenses/results/port_host_corecording_partial_20261030/port_alien_host_butterfly_corecording_audit.json
- provenance/reviewer_defenses/results/host_island_synthesis_20261030/forecast_temporal_results.json
- provenance/reviewer_defenses/results/host_island_synthesis_20261030/forecast_spatial_results.json

**Fig 5**
- provenance/reviewer_defenses/results/original230_dated_host_join_20261031/original230_prior_host_recurrence_audit.json
- provenance/reviewer_defenses/results/host_geographic_prevalence_null_20261031/host_geo_frequency_null_results.json
- data/exploratory/public51_butterfly_index_20261109/full51_source_validation.json
- results/public51_butterfly_trait_access_20261109/public51_trait_summary.json
- results/public51_butterfly_trait_access_20261109/public51_butterfly_conditional_trait_marine.json
- results/public51_author_selected_migrant_access_20261109/author_selected_migrant_no9_15_source_null.json

**No figures have been newly rendered from this updated v0.2 plan yet.** Existing earlier JBI results are source-audited but different manuscript figures. Integrate them only after numerical and visual source-denominator checks.

## Figure caption discipline

- Report native→contemporary percentages with denominator native resource combinations.
- Report matched island resource as **percentage-point difference in contemporary added-resource share**, not percentage increase.
- Report port betas as three **separate log-outcome models on 98 regions**, not a common causal decomposition or percent butterfly benefit.
- Report geographic source-distance tests as observational *GBIF index* and *not* measured butterfly flight.
- Report 208/1163 source pairs as **observed compatibility** only; its family/degree null does not preserve botanical geography; do not reuse the 34/113 geography-preserved null as though applied to the 208/1163 cohort.
- Draw the null/negative test adjacent to the strongest positive interpretation in Fig 5, without dismissing the genuine earlier compatible host observations.
- Use 5 final main figures and keep messy fully source-audited but mechanistically unidentifiable AIS/wind/shipping current maps in Supplementary.

**Journal direction:** Global Ecology and Biogeography for a global cross-scale butterfly resource-geography argument with islands and ports as the headline natural experiment analogue; Journal of Biogeography is fallback for the same manuscript only. Both require a clear boundary between association and a causal arrival pathway.
