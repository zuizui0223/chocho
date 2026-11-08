# Geographic verification of 30 original WPI GLSN–UNCTAD ports: 24 real registry identities, 276 source-only port-pair comparisons

**Completed 2026-11-08. Original 51-island/194-WPI-port macroecology datasets and their previous positive/adverse results have not been changed. No butterfly/community regression or journal submission was performed.**

## Main contribution

Original PR146 recovered exact publisher 2018 GLSN weighted liner co-service source (887 ports, 14,052 undirected modelled port-pair edges), and 2006-Q1 to 2021-Q3 UNCTAD PLSCIQ port connectivity series (1,300 historical port index series; 81,900 port×quarters). In the ORIGINAL WPI2019 194-source port cohort these archives have **30 ports with globally unique source country/name matches to BOTH GLSN2018 and historical UNCTAD PLSCIQ**. Neither GLSN nor UNCTAD original table gives port centroid coordinates, so country/name agreement alone was not an independent physical harbor identity.

PR147 independently checked 17 of the 30 ports additionally WPI identity matched to Zenodo's selected 2020 AIS OD source: 15 passed strict PortWatch geographic check and 2 (Sheerness and Harwich) were not matched. In this PR, the **remaining 13 original WPI double-source candidates with NO existing Zenodo2020 WPI endpoint** were selected BEFORE new public registry queries, and **9 of 13 independently matched** to an original IMF PortWatch registry identity; 4 were unmatched. No butterfly observations, host plant records or retrospective ecological scores informed either source sample.

**Overall: 24 of 30 original GLSN∩UNCTAD country/name candidate ports are independently geoverified, spanning 17 of 51 original island/landmass groups, each with an IMF PortWatch registry portid, independently recorded lat/lon, and a provider-reported UN/LOCODE.** The remaining **6 original ports remain unverified**: WPI758 Sheerness, WPI764 Harwich, WPI146 Miri, WPI284 Benoa, WPI416 Sibu, WPI1516 Antsiranana. An unverified provider match is NOT evidence that the port is closed, lacks actual shipping or has no butterfly resource opportunities.

## New cohort: exact 13 original port identities and independent result

Official IMF PortWatch registry: https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/PortWatch_ports_database/FeatureServer/0 . Full registry 2,065 ports was NOT harvested. The [source-frozen acquisition protocol](PROTOCOL.md) limited queries to exactly **13 prespecified WPI centroid boxes** (0.35° lat/lon half-width), one count-only and one at-most-15-feature sample per box: **26 official queries / 11 provider example rows**. A match required EXACT normalized port name, ISO3 via pycountry and ≤25 km original WPI-to-provider great-circle distance, with exactly one qualifying provider record and a complete nontruncated box. A separately written Python-standard-library auditor reconstructed all requests/source features and status decisions, with no guessed source locations.

| New original WPI port ID and name | Original island group | Independently geoverified PortWatch ID | Provider LOCODE | WPI-to-PortWatch distance |
|---|---:|---|---|---:|
| 96 Korsakov | 28 | port594 | RU KOR | 6.47 km |
| 98 Kholmsk | 28 | port2444 | RU KHO | 1.00 km |
| 402 Samarinda | 9 | port1137 | ID SRI | 9.24 km |
| 408 Banjarmasin | 9 | port112 | ID BDJ | 5.52 km |
| 410 Pontianak | 9 | port925 | ID PNK | 5.04 km |
| 415 Kuching | 9 | port604 | MY KCH | 7.97 km |
| 576 Semarang | 19 | port1179 | ID SRG | 0.82 km |
| 656 Palembang | 12 | port875 | ID PLM | 6.51 km |
| 720 Panjang | 12 | port881 | ID PNJ | 3.29 km |

**The 4 new unverified ports:** Miri, Benoa, Sibu, Antsiranana have no unique port-name/ISO3/position match inside their completely sampled source boxes. A provider-reported LOCODE is a useful port key, but its presence is not independently authenticated against the UN/LOCODE issuing body. Full rejected/unmatched source IDs and reasons are retained.

Public evidence: [all 13 original WPI port rows](WPI13_independent_PortWatch_source_geography_LOCODE.csv), [provider-count/selection and checksum receipt](WPI13_PortWatch_source_identity_receipt.json), [independent per-source-feature audit](independent_WPI13_source_audit.json). Original bounded 26 request details and 11 provider feature records are stored in a time-limited GitHub Actions source artifact rather than published as an IMF global registry copy. Source-only scientific [13-port registry CI](https://github.com/zuizui0223/chocho/actions/runs/42436460525) succeeded; its audit request-count metadata was clarified from an inherited old-cohort display label, not by changing source data or any port-match decisions.

## Consolidated 24 source-matched ports, original 276 port pairs and all 1,275 island pairs

The 15 independently verified ports from [PR147](https://github.com/zuizui0223/chocho/pull/147) and 9 newly verified ports were combined with exact original WPI IDs, separately frozen [combined 30/24-port protocol](COMBINED30_PROTOCOL.md), and original publisher-sourced GLSN2018 modelled edges/GC/weights, original Zenodo2020 AIS-inferred WPI port-to-port OD source rows and UNCTAD 2006–2010 historical port-index coverage. A port contributes only if its actual independent PortWatch registry identifier, name, country and original WPI centroid have been checked; the six failed port identities remain outside the 24-port pair denominator.

| Outcome, strictly source-only | Original 15 verified WPI ports | Expanded 24 independently verified WPI ports |
|---|---:|---:|
| Original island groups represented | 14 | **17** |
| Possible unordered matched WPI port pairs | 105 | **276** |
| **GLSN2018 modelled liner co-service** port-pair edges within this source subset | 6 | **20** |
| Original unordered island-pair types touched by these GLSN **modelled** edges | 6 | **9** |
| Actual 2020 Zenodo AIS-inferred OD source matching: distinct unordered WPI port pairs | 10 | **10** |
| AIS2020 directed WPI port OD source links | 17 | **17** |
| AIS2020 matched source triplet rows across the matched ports | 64 | **64** |
| Original island-pair types with matched AIS2020 source OD rows | 7 | **7** |
| Distinct original WPI port pairs with BOTH 2018 modelled co-service and 2020 recorded AIS-inferred OD | 2 | **2** |
| Source WPI port indices with ≥16/20 nonmissing UNCTAD 2006–2010 quarters | 14 | **22** |

The two original WPI port-pair identities shared by GLSN2018 model and Zenodo2020 source remain **WPI28 Belfast ↔ WPI1143 Southampton** and **WPI428 Santa Cruz de Tenerife ↔ WPI430 Puerto del Rosario**. Source original GLSN co-service weights 800 and 4,606 **are NOT number of vessel passages**; the separate Zenodo2020 AIS source has 15 and 10 inferred OD triplet rows on these port pairs. The fact that there are NO **additional** AIS recorded port-pair observations after admitting the 9 new geographic source port matches is expected: the new 13 original WPI candidate ports were *selected precisely because they were absent from the existing original-WPI-matched 2020 AIS excerpt*. That absence is an observation coverage gap, not evidence that marine connectivity was absent.

The added 14 GLSN modelled source port-pair edges reveal additional *model co-service* structure, especially among original source island groups 9, 12 and 19, as well as Korsakov–Kholmsk within original island group 28. These are **source-modelled route opportunities**, not observed direct passages or inferred plant/butterfly ship cargo. The entire original 51 island/landmass universe has **1,275 unordered original island pairs**, many with NO verified-port source observation; those source zeros are NOT observed biological absences and must not be treated as no shipping.

Complete verified-data deliverables:

- [All 276 matched original WPI24 port pairs](WPI24_geoverified_all276_GLSN_UNCTAD_AIS_port_pair_comparison.csv): both WPI/PortWatch IDs and provider LOCODE, source GLSN modelled edges and original GC/weight, original AIS2020 two-direction inferred OD counts, original PLSCIQ historical observation quarters/index means, explicit no-traffic-inference flags.
- [All 1,275 original WPI51 island pairs](original51_all1275_geo_verified24_source_island_pair_observation.csv): source port coverage, number of original GLSN *modelled* edges, AIS2020 *recorded* source triplet counts, observed-versus-modelled match status and no biological claim.
- [Complete source pair/audit receipt](WPI24_3_source_verified_port_graph_audit.json) and [independent numerical audit](independent_WPI24_3source_graph_audit.json).
- [Source-only reproduction](../../scripts/summarize_original_WPI24_geoverified_maritime_sources.py) and [independently authored Python-stdlib auditor](verify_WPI24_cross_source.py).

The independent numerical audit was run after repairing an **outdated metadata key** inherited from the earlier 15-port source-audit script. It then reconstructed all **276 distinct port-pair** and **1,275 original island-pair** values, every GLSN model edge and exact weight, original UNCTAD port index source availability, and all original 2020 AIS directed OD counts. Crucially it compared all original **105 prior verified 15-port pair rows** against their corresponding subset of the new 276-pair table: **all historical source metrics reproduced exactly**. [Recomputed 24-port scientific CI](https://github.com/zuizui0223/chocho/actions/runs/42436864992) succeeded, including independent source audit.

## Limits and next ecological decision

This is a deliberate **data-quality and cross-source coverage** advance, not newly demonstrated global butterfly movement or a causal result. The 30-port source cohort was chosen because GLSN and UNCTAD already contained country/name candidates; it is NOT a representative random 194-port sample. GLSN original 2018 edges represent modeled liner-service co-membership, NOT individual direct observed vessel trips. UNCTAD historical PLSCIQ values are a port connectedness index (available from 2006), NOT bilateral 2006–2010 voyages. Zenodo AIS2020 OD source episodes are selective/inferred and cover Jan 1–Oct 7, 2020; no source 2020 event cannot be treated as physical shipping absence or as evidence of no larvae/seed transport.

The full 51-island real-route study's previous minimum 80 *verified observed island-pair links* is NOT met by counting 20 modelled GLSN port-pair edges or 9 modeled island-pair types. Model edges cannot silently be converted into confirmed 2020 voyages. No butterfly distribution prediction has been fitted using any newly joined source panel; no 2001–2010 plant introductions, true butterfly feeding/breeding, successful establishment or causal ship/cargo transport has been inferred.

The previous **global GEB +54.9% potential larval host portfolio opportunity**, prior positive AIS later-community predictive increments, negative early-group response metrics and the PR137 750km geography-node AIS placebo limitation remain independently preserved. No journal submission or original GEB paper alteration has occurred.
