# Actual 2020 port-to-port AIS source, conservatively matched to the original WPI194

**Status: source-matching and independent audit completed, 2026-11-07. No butterfly/plant-route prediction has been fitted.** This adds direct AIS-derived origin/destination events to the existing chocho maritime analyses without altering the established global GEB resource-expansion and island-community findings.

## Source and original cohort

Duvallet et al. (2024), *Date-Ship-Port triplet data based on AIS data*, [Zenodo DOI:10.5281/zenodo.10783698](https://doi.org/10.5281/zenodo.10783698), CC BY 4.0. Original archive triplet_2020.zip: 32,453,859 bytes, publisher MD5 af27e1bef1f706e2e29ed53aee47ff10, source SHA256 6d5f1273ba98024ea485dc639f8a1c1ea425b83bb95e01d725cb54a8c4aa5ec0. It contains a 113,764,323-byte CSV with from/to dates, port identifiers, names and positions. No vessel-identifying field is republished here.

**Available source departures span January 1 to October 7, 2020**; do not report this as a complete calendar year. Source records are AIS-derived inferred port-origin/destination episodes, not independent proof of nonstop cargo voyages.

The original pinned WPI2019 cohort consists of **194 biological ports, 182 marine-qualified terminals and 51 island/landmass groups**, source blob a4681c3e3c1406d3a37f51464a298286e9d0ee4f. The external Zenodo port identifier and original WPI OBJECTID are different namespaces and were never equated numerically. Primary linkage required normalized **exact port name plus no more than 25 km original source-coordinate separation**, uniquely identifying one WPI port. No source-only location match without name identity was treated as verified.

## All published 2020 voyage source and WPI crosswalk

| Observation support | Directly audited records |
|---|---:|
| Original 2020 source origin/destination rows scanned | **645,581** |
| Source rows with usable 2020 departure/arrival chronology and positions | **645,581** |
| Unique named/unnamed external port-position source nodes | **1,923** |
| External nodes with one exact-name and ≤25 km original WPI match | **46** |
| Distinct original WPI ports matched | **46 / 194** |
| Original island groups with at least one matched WPI endpoint | **26 / 51** |
| Source triplets with both original WPI endpoints matched | **442** |
| Source triplets with different original WPI ports in the SAME island group | **309** |
| Source triplets between different original island groups | **133** |
| Different directed original WPI port-to-port links (includes within-island) | **94** |
| Directed original island-to-island links | **20** |
| Distinct unordered original island pairs with at least one source-observed movement | **12 / 1,275** |

The original source 1,923 nodes include 391 lacking reliable names and 1,486 failing the exact-name-plus-distance WPI mapping. These are **unmatched source coverage**, not biological evidence that an island has no voyages.

### Every source-observed original island-pair link

| Original landmass IDs | Matched 2020 origin/destination source records |
|---|---:|
| 15–26 (e.g., British mainland WPI ports ↔ Belfast-area ports) | **77** |
| 206–397 (Santa Cruz de Tenerife ↔ Santa Cruz de La Palma) | **20** |
| 15–85 (UK mainland-side ports ↔ Cyprus) | 10 |
| 206–234 (Tenerife ↔ Puerto del Rosario/Fuerteventura) | 10 |
| 54–85 | 4 |
| 85–206 | 3 |
| 26–206 | 2 |
| 54–206 | 2 |
| 214–242 (Kahului/Maui ↔ Barber's Point/Oahu) | 2 |
| 10–23 | 1 |
| 15–21 | 1 |
| 26–85 | 1 |
| **Total source-observed inter-island rows** | **133** |

The dominant original group link 15–26 has **77/133** of all matched cross-island observations. The observed episodes do not independently establish cargo, ship-borne host plants or butterfly movement. Their uneven distribution and selective original WPI matching cannot be extrapolated to the remaining original island pairs.

## Precommitted source-only feasibility gate failed

The [frozen port identity and coverage protocol](CROSSWALK_PROTOCOL.md) required ≥15 original island groups with a verified 2020 endpoint **AND** ≥80 original unordered island pairs with matched OD episodes before fitting a 51-island route-based biological predictor. The source-island endpoint gate is met (**26**); the link coverage gate fails (**12**). **No direct-route butterfly/plant prediction or formal shipping-causality inference is attempted.** The other 1,263 original island pairs have *zero matched episodes in this limited source*, not zero actual voyages.

## Pending source identity review, no inferred routes

Of 1,486 external nodes with a reliable but unmatched source name, **27** are no more than 10 km from a WPI port and **18** are no more than 5 km away. The full [27-row geographical-neighbor candidate ledger](source_only_10km_unverified_port_candidates.csv) explicitly says identity_verified=False for every row, and none of these locations contributes to the matched 133 source episodes. They include possible name aliases (Hull versus KINGSTON UPON HULL, Rodbyhavn versus RODBY HAVN), but also clearly unsafe automatic assignments (Felixstowe versus HARWICH, London versus GRAVESEND, Koh Si Chang Transhipment Area versus MANILA). **Proximity is NOT port identity**. Each needs independent port authority/UNLOCODE and positional validation before promotion; do not change the primary counts by silently assuming a synonym or using only a nearest neighbor.

## Auditable source records and limits

The [Zenodo ZIP/MD5 inspector](../../scripts/audit_island51_zenodo2020_voyage_source.py) and [full source-to-WPI OD reconciliation](../../scripts/reconcile_island51_zenodo2020_WPI_port_voyages.py) can rerun against the original DOI archive. The [independent standard-library port identity auditor](verify_WPI_crosswalk.py) checks all original WPI rows, normalized exact-name/geodesic identity, aggregate OD links and all 1,275 pair-level matched-observation flags. Durable files: [original 2020 ZIP source receipt](original2020_AIS_triplet_schema_receipt.json), [1,923 endpoint identity decisions](external2020_endpoints_to_WPI194_crosswalk.csv), [194 original port coverage rows](original194_source_matched_port_coverage.csv), [94 WPI port-to-port links](observed2020_verified_WPI_port_links.csv), [all original 1,275 island pairs](original51_all1275_direct2020_source_pair_observation.csv), [final source receipt](source_matched_voyage_receipt.json) and [source independent audit](independent_crosswalk_audit.json). Full third-party raw ZIP is referenced by DOI instead of being committed.

**Ecological limits:** 2020 vessel events postdate 2001–2010 plant/butterfly indexing. Earlier indexed general plants are not confirmed butterfly hosts or dated introductions; original candidate host predictors pool 2001–2024 and may postdate later butterflies. Butterfly first indexing is not arrival, reproduction or colonization. AIS source triplets do not document seed/larval cargo. Previous GEB +54.9% potential host-resource opportunity, original AIS+wind predictive correlations, and strict local-750 km placebo limits retain their original distinct inferential scopes.
