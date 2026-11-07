# Frozen IMF PortWatch static maritime-connections discovery, bounded metadata only

Date 2026-11-07. The previous public Mendeley Data API for dataset DOI 10.17632/vzzy3b9gg4.1 returned HTTP 401 on requested metadata/files/folder snapshot, so no network file bytes or original-WPI port matches were retrieved. Existing PR141 source 2020 observed voyage original 51-island links remain 12 verified and PR142 all-27-unverified-identity what-if remains 28 hypothetical links (NEITHER a complete transport map). Keep all earlier positive AIS potential prediction and negative early group response intact.

## Better direct-source candidate

International Monetary Fund PortWatch (with University of Oxford/UN Global Platform). Official PortWatch https://portwatch.imf.org/ and public ArcGIS FeatureServer organization https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services. IMF PortWatch port registry layer publicly documented at PortWatch_ports_database/FeatureServer/0 with portid, portname, country/ISO3, lat/lon, UN/LOCODE and vessel_count columns. PortWatch Port Monitor shows top five up/downstream port connections derived from AIS 2019–2025. A separate peer-reviewed source describes a static PortWatch Maritime Connections table with approximately 155,521 complete rows and 11,924 directed subnetwork links among 202 selected study ports (Al-Azazi 2026 online, Reliability Engineering & System Safety 2027 issue, DOI 10.1016/j.ress.2026.113682); this establishes potential source existence, NOT our own dataset retrieval, coverage or verified direct passages. The model network may be static aggregation, not vessel-level trips. User's original WPI194 source is unchanged, Git blob a4681c3e3c1406d3a37f51464a298286e9d0ee4f.

## Access, permissions and bounded metadata workflow

The official IMF Copyright and Usage terms effective October 11, 2024 allow certain reuse of published statistical Data with citation, but prohibit automated bulk copying without explicit permission where applicable. Because the specific hosted global maritime-connections feature layer's classification, copyright and redistribution terms have NOT been separately inspected, DO NOT run a systematic bulk download of 155k rows or republish provider records as our own. This stage is a bounded public ArcGIS layer discovery/source-schema probe, not full extraction.

Perform at most:
- one service-directory metadata query (f=json), and one official ArcGIS public search API metadata query for maritime connections to find exact source service/item names;
- one official port registry layer schema (known service URL), and at most four source-discovered Maritime Connections/Network or similarly named layer schema documents;
- if layer genuinely exists and is publicly readable, request returnCountOnly=true, and at most FIVE example rows, only source port IDs, directional endpoint keys, link aggregate/period names and no personal vessel identifiers. Treat fields as names/metadata, not proof of voyage observation. Capture layer page provider copyrightText, editing time, source ID namespace, published attribution and any errors.
- only public nonauthenticated GET with bounded response bytes/timeouts. No login bypass, no brute-force enumeration. No source full data extraction, no 194 individual port-ID queries yet. Archive source metadata (not large provider rows), source URL, retrieval date, id/field/schema detail and provider's explicit Terms of Use.
- verify original WPI194 source count, 51 landmasses and 182 original marine qualified terminals by original source Git blob.
- no matched WPI pairs, AIS biological model, inferred cargo transport, colony arrival or causal conclusion until a separately frozen original port matching and provider rights plan is accepted.

## Next conditional step if a source exists and access conditions permit
A scientifically usable global connection table must provide resolvable port origin/destination keys, distinct port registry with lat/lon and verifiable original WPI matching, and at least source-only 80 unique original landmass pairs with verified provider links. All missing original islands/pairs must be flagged unobserved, not absent. Provider edges may represent aggregate AIS 2019-2025 shipping CAPACITY, not the 2020 single-voyage trips or a confirmed 2001–2010 dispersal opportunity, so never merge those estimands silently.

Source-only stage may end in no accessible network or rights restriction. Report that concrete result rather than claiming model success. Prior GEB +54.9% potential larval host opportunity, 51-island host/AIS predictive studies, and source temporal leakage limitations remain unchanged.
