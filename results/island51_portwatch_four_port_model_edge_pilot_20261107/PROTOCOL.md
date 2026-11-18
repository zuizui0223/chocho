# Four source-verified WPI ports: bounded PortWatch directed model-edge pilot

Date 2026-11-07. Source-only exploratory diagnostic. Parent PR144 confirmed IMF PortWatch registry of 2065 ports and Spillover Simulator Maritime Connections table with 155521 modeled edge records, containing from_id, to_id, from/to geographic coordinates, distance, transit_days and input_days_* simulation parameters. The earlier pilot had four source-verified WPI port identities out of five enriched locations, and earlier actual 2020 vendor triplets. These were known before this new directed-link query; results for the current four-port directed-links question were not.

## Four source-verified original WPI and PortWatch IDs (freeze)
- WPI28 BELFAST island26 GB -> PortWatch port141, matching original WPI coordinate within 0.856 km.
- WPI428 SANTA CRUZ DE TENERIFE island206 ES -> PortWatch port1153, within 1.258 km.
- WPI817 LIMASSOL island85 CY -> PortWatch port647, within 3.646 km.
- WPI1416 COLOMBO island31 LK -> PortWatch port254, within 1.128 km.
Original fifth WPI1525 KAHULUI was unverified by provider registry pilot; do not impute.

The source pinned prior receipt is results/island51_portwatch_global_connections_20261107/WPI5_bounded_registry_source_pilot_receipt.json. The original 51 landmass/194 biological port/182 marine terminal cohort and frozen source git blob a4681c3e3c1406d3a37f51464a298286e9d0ee4f are unchanged. No butterfly/plant response enters the source selection.

## Exact bounded public query
Official provider network URL:
https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/Spillover_Simulator_Maritime_Connections/FeatureServer/0

For exactly 4*(4-1)=12 ordered pairs of verified source PortWatch port IDs: one public count-only (returnCountOnly=true) query per ordered condition from_id='portA' AND to_id='portB'. For each condition with positive count, request AT MOST one source row, returning provider fields from_id, to_id, from_name, to_name, from_iso3, to_iso3, from_lat, from_lon, to_lat, to_lon, distance, transit_days. Absolutely no bulk downloads, no all-155521 rows, no source vessel identifiers or time-series harvest. Queries that fail must be marked ERROR, not interpreted as zero edges. In addition perform one count-only outdegree (from_id) and one count-only indegree (to_id) per each of four ports (8 count requests). Budget at most 20 count-only plus 12 one-row samples. Keep all 12 ordered pairs even if zero.

Verify each returned model-edge sample includes identical from_id and to_id and provider source coordinates are consistent with the original WPI and registry pilot within the original <=25km identity threshold. No network model value is called the number of actual voyages or cargo, and model simulator 'input_days_*' values represent disaster disruption scenarios.

## Separate previous 2020 real-OD comparison
Use source-archived old Zenodo 2020 WPI origin/destination list results/island51_zenodo2020_portvoyage_audit_20261107/observed2020_verified_WPI_port_links.csv, which has 94 source-matched directed port edges. For same exact 12 ordered WPI pairs report the independently verified matched 2020 vendor triplet count; zero means no matched 2020 source rows, NOT no actual shipping. Compare provider network STATIC MODELED EDGE PRESENCE with earlier real-2020 matching without assuming equivalence, agreement proof or historical time consistency. Report source-only pair comparison without ecology.

## Provider restrictions and inference
IMF Copyright and Usage https://www.imf.org/en/about/copyright-and-terms distinguishes licensed statistical data from general content and prohibits automated BULK harvesting without permission. The above selected 12-pair and four-port small sample is deliberately not a systematic extraction. Permission required before all 194 original ports or all 155521 model rows are harvested or redistributed where covered by the IMF general terms. Cite 'Source: International Monetary Fund, IMF PortWatch' and official endpoints. No 51-island butterfly regression or ship-borne host plant/larval colonization mechanism. Model edges are STATIC spillover network, not confirmed AIS vessel trips, and from/to ports do not independently verify cargo movement.
