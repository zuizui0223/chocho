# Four original WPI ports: PortWatch static network edges versus 2020 real source triplets

**Complete bounded source-only empirical comparison, 2026-11-07.** The four original WPI/PortWatch port identities were verified in parent [PR #144](https://github.com/zuizui0223/chocho/pull/144) using exact normalized name, ISO country and ≤25 km distance. The present twelve link conditions are **fixed before the new link API queries** in [PROTOCOL.md](PROTOCOL.md), without butterfly or plant response inspection.

## New result

Across the **12 ordered pairs** formed from four preverified original WPI ports, PortWatch's officially hosted **static Spillover Simulator Maritime Connections model** includes **3 directed link records**. The separate [Duvallet et al. 2020 Zenodo origin/destination source](https://doi.org/10.5281/zenodo.10783698), previously matched conservatively to these same WPI ports in PR #141, contains triplet rows on **5 directions**. All **3** provider modeled directions have an independently matched Zenodo 2020 origin/destination record, while **2 Zenodo-positive directions have no PortWatch static model record** in the current bounded source query.

| Original WPI origin → destination | PortWatch static directed modeled edge records | Previous independently matched Zenodo 2020 triplet rows |
|---|---:|---:|
| Belfast → Santa Cruz de Tenerife | 0 | 1 |
| Belfast → Limassol | 0 | 0 |
| Belfast → Colombo | 0 | 0 |
| Santa Cruz de Tenerife → Belfast | 1 | 1 |
| Santa Cruz de Tenerife → Limassol | 1 | 1 |
| Santa Cruz de Tenerife → Colombo | 0 | 0 |
| Limassol → Belfast | 1 | 1 |
| Limassol → Santa Cruz de Tenerife | 0 | 2 |
| Limassol → Colombo | 0 | 0 |
| Colombo → Belfast | 0 | 0 |
| Colombo → Santa Cruz de Tenerife | 0 | 0 |
| Colombo → Limassol | 0 | 0 |

A **PortWatch record count of 1 is a model-network edge**, not one AIS vessel voyage or one ship carrying butterfly larvae. A model-network count of 0 means no row in the presently hosted simulator for that specific pair, **not** that shipping did not occur. Conversely, a Zenodo 2020 zero means no previously matched source triplet among the 2020 Jan–Oct source records, not zero actual voyages. Source systems have different observation windows, coverage, estimands and missingness. This four-port enriched nonrandom pilot is **NOT a statistical validation of a global source network**.

## The provider's four port IDs have many model-network records in the complete source layer

The following are the provider's **count-only** query results for all source model edges where each preverified PortWatch port is an origin or destination. They are counts of static source **edge records**, NOT unique ships, departures or distinct biological colonization corridors.

| Original matched WPI port | Provider modeled origin records | Provider modeled destination records |
|---|---:|---:|
| Belfast (WPI 28) | 222 | 252 |
| Santa Cruz de Tenerife (WPI 428) | 161 | 157 |
| Limassol (WPI 817) | 109 | 116 |
| Colombo (WPI 1416) | 194 | 165 |

The official source [PortWatch ports registry](https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/PortWatch_ports_database/FeatureServer/0) reports 2,065 source ports, and the [Maritime Connections source layer](https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/Spillover_Simulator_Maritime_Connections/FeatureServer/0) reports **155,521 static model connection records** (at source audit time). Original WPI provider pilot matches: WPI28 Belfast ↔ port141, WPI428 Santa Cruz de Tenerife ↔ port1153, WPI817 Limassol ↔ port647, WPI1416 Colombo ↔ port254. Original WPI1525 Kahului did **not** uniquely match the pilot and was excluded without imputation.

## Fixed numerical method and independent audit

1. Verify parent WPI194 source Git blob, original 51 landmass IDs and 182 marine-qualified terminals; require all four preverified PortWatch provider source IDs and original WPI names, countries and original-port coordinates.
2. Issue exactly **12 provider count-only queries**, one for each ordered pair of different original four WPI ports. On each positive model count, request at most **one provider network example** (three rows total observed here), verifying actual from/to source IDs and coordinates within 25 km of original WPI geography.
3. Issue exactly **eight additional count-only requests** for model out/in edge-record counts of the four PortWatch source IDs. No bulk pagination, all-network export or automated collection of 155,521 rows.
4. Read the existing original 94-edge verified Zenodo 2020 WPI port-to-port table and compare previous individual source triplet-count support for each of the same 12 directions.
5. Independently re-evaluate all 12 original port pairs and all eight degree queries with a separate Python-standard-library auditor. No ecological regression or 51-island original butterfly response has been fitted or modified.

The [source receipt](source_four_port_directed_model_edge_pilot.json) permanently records query result values, original WPI IDs, original island groups, source server response attempts, and the three model-link example fields. Additional machine-readable compact source tables: [12 directed-pair comparison rows](source_four_port_directed_model_edges_vs_2020_vendor_triplets.csv), [4 port in/out degree model counts](source_four_WPI_model_in_out_edge_record_counts.csv), and [independent source auditor](verify_four_directed.py).

## Inference, legal and study positioning

The [IMF Copyright and Usage terms](https://www.imf.org/en/About/copyright-and-terms) restrict automated bulk downloads without permission and distinguish statistical Data from other site content. This study made **20 count-only and 3 single-record queries only**, and attributes the hosted data to **Source: International Monetary Fund, IMF PortWatch**. Any systematic extraction of all 155,521 records or all original 194 WPI candidate matches must first resolve provider-specific permission and proper attribution/redistribution rights; this four-port sample is not a right to bulk harvest.

**Biological limitations:** contemporary/global static PortWatch model edges are neither the original 2020 Zenodo source ship triplets, historical 2001–2010 actual voyages, cargo of host plants/seed, larval host plant establishment, butterfly movement/feeding, nor confirmed island colonization. Sample was enriched from successfully name-matched ports and is subject to selection bias. Earlier global GEB +54.9% potential host portfolio opportunity, the original PR #137 750-km AIS node-placebo specificity failure, PR #138 194-port historical general indexed plants and PR #140 990-pair cross-community predictive results remain separate and unchanged.
