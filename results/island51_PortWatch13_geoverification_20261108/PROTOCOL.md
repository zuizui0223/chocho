# Frozen independent source-based PortWatch verification of 13 new GLSN+UNCTAD original WPI ports

Date 2026-11-08. This follows Draft PR147 in which exactly 17 original WPI ports were selected because they appeared in the source country/name matched GLSN2018 liner model, historical UNCTAD PLSCIQ and independently verified Zenodo2020 AIS inferred OD; 15 of those 17 passed independent IMF PortWatch registry name+ISO3+within 25km checks. PR147 data and all prior positive and adverse ecological results remain unchanged.

## Locked 13 new source-only ports
Original WPI194 geography Git blob a4681c3e3c1406d3a37f51464a298286e9d0ee4f (51 original landmasses / 182 marine-qualified terminals) and original PR146 194-row cross-source candidate CSV are immutable. Select original ports with a unique source GLSN2018 port name/ISO country AND an old UNCTAD 2006–2021 name/country series AND no previously matched original Zenodo2020 WPI port endpoint. The last condition concerns source observability, not absence of ships or butterflies.

Exactly **13** original WPI IDs in numeric order are locked: 96 Korsakov (RU), 98 Kholmsk (RU), 146 Miri (MY), 284 Benoa (ID), 402 Samarinda (ID), 408 Banjarmasin (ID), 410 Pontianak (ID), 415 Kuching (MY), 416 Sibu (MY), 576 Semarang (ID), 656 Palembang (ID), 720 Panjang (ID), 1516 Antsiranana (MG). Five landmass groups: 9, 10, 12, 19, 28. The complete GLSN2018∩UNCTAD country/name candidate set is 30 ports: 17 old triple-source candidates + these 13. No butterfly/plant observation or outcome informed this selection.

## Strict and bounded official registry requests
Official IMF PortWatch 2065-port registry:
https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/PortWatch_ports_database/FeatureServer/0

Only the exact 13 original WPI coordinate boxes are queried, with 0.35 degrees on each latitude/longitude axis. Two GET requests per box: returnCountOnly and at most 15 source rows restricted to the provider fields portid, portname, country, ISO3, lat, lon, LOCODE. Maximum **26 official public API requests, 195 provider source example rows**, no pagination, no total registry copying or any of the 155521 maritime model edge records. Provider reuse and attribution are subject to IMF terms. Do not silently treat incomplete/truncated boxes as absence.

A source WPI port is independently matched only when a complete source sample has exactly one original country ISO2 mapped via pycountry to the same provider ISO3, the exact punctuation/whitespace-normalized NFKD name agrees, and WPI centroid to provider registry point has spherical haversine distance no greater than **25 km**. Retain original WPI ID/island, original WPI name/country/coordinates, provider ID/LOCODE if verified, geodesic and reasons for every rejected/ambiguous/incomplete row. The PortWatch provider port ID is not the original WPI port ID; the reported LOCODE is not independently certified against its issuing authority.

## Independent reproducibility and model stop
Publish complete 13-row original WPI source outcome and compact source receipt only. Full bounded provider example response/query audit is a 30-day CI artifact rather than a bulk PortWatch dataset in Git. A separately authored Python-standard-library auditor independently recomputes every original source-selection identity, provider query bounds, all source name/ISO3/haversine checks, and failure statuses.

A separate follow-up may merge these verified 13 with previous 15 geoverified old source WPI IDs, preserving the 2 prior nonmatches, and calculate source-only 2018 GLSN candidate model edges and UNCTAD historical index coverage. NEVER treat a missing Zenodo2020 source endpoint as no ship traffic, a GLSN model co-service edge as observed voyage, a historical connectivity index as cargo volume, or any of them as actual butterfly colonization.

No butterfly model and no journal submission. Preserve global GEB +54.9% potential host-resource expansion, all original AIS predictive results and PR137 local geography placebo limitation.
