# Independently source-verified WPI ports shared across maritime datasets

**Source-only ecological infrastructure result, 2026-11-08.** Original study 51 island/landmass groups, 194 WPI botanical ports, 182 marine-qualified terminals. **NO butterfly model, colonization analysis or journal submission.** This is a separate independent port-identity validation for maritime data already collected in [PR #146](https://github.com/zuizui0223/chocho/pull/146).

## What became independently verifiable?

Previously GLSN2018 liner co-service network ports and the 2006–2021 UNCTAD PLSCIQ port indices were linked to original WPI194 **only by country and normalized name**. Their source sheets had no original port latitude/longitude, so 34 GLSN and 40 UNCTAD WPI entries were still mere possible harbor identities. From these, exactly **17 original WPI ports** were already separately matched to Zenodo 2020 AIS port-origin/destination positions and shared both the GLSN and UNCTAD source-name/country keys. This 17-port cohort was fixed entirely from source information, prior to new PortWatch queries, without examining butterfly outcomes.

The official [IMF PortWatch port registry](https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/PortWatch_ports_database/FeatureServer/0) contains independently recorded portid, port name, ISO3 country, latitude, longitude and a registry-reported UN/LOCODE. Following a [source-frozen 17-port protocol](PROTOCOL.md), we queried **only 17 limited original-port geographic boxes**, each with one count-only and one maximum-15-row registry sample: **34 public queries, 31 returned example records**, NOT a bulk download of PortWatch's 2,065 global ports or its 155,521 model-network records. Provider rights and attribution remain applicable.

A source candidate is accepted only if exactly ONE provider port within the WPI search box has **canonical exact name, identical ISO3 country and spherical haversine distance ≤25 km** from the original immutable WPI port. A source box with incomplete records must fail, not claim a missing biological harbor.

**Outcome: 15 of 17 original WPI candidates were independently position-confirmed, representing 14 different original landmass groups.** All 15 have a PortWatch-reported LOCODE; provider registry LOCODE presence is NOT independent UN/LOCODE issuance certification. Exact four original provider port-ID matches from the previous WPI five-port pilot were reproduced. Sheerness (WPI 758) and Harwich (WPI 764) returned no qualifying exact name/country/position PortWatch match; they remain **UNVERIFIED**, not absent or closed ports.

### All fifteen verified source ports

| Original WPI ID | Original port | IMF PortWatch portid | Registry LOCODE | WPI-to-provider distance |
|---|---|---|---|---:|
| 28 | Belfast | port141 | GB BEL | 0.86 km |
| 274 | Cagliari | port203 | IT CAG | 3.21 km |
| 404 | Balikpapan | port102 | ID BPN | 5.76 km |
| 428 | Santa Cruz de Tenerife | port1153 | ES SCT | 1.26 km |
| 430 | Puerto del Rosario | port1046 | ES FUE | 0.82 km |
| 668 | Manila | port694 | PH MNL | 10.19 km |
| 682 | Haikou | port437 | CN HKO | 1.98 km |
| 767 | Port Louis | port965 | MU PLU | 0.44 km |
| 817 | Limassol | port647 | CY LMS | 3.65 km |
| 876 | Immingham | port493 | GB IMM | 1.98 km |
| 916 | Corner Brook | port2437 | CA CBK | 2.09 km |
| 1143 | Southampton | port1216 | GB SOU | 2.08 km |
| 1416 | Colombo | port254 | LK CMB | 1.13 km |
| 1515 | Toamasina | port1301 | MG TOA | 1.41 km |
| 1546 | Napier | port794 | NZ NPE | 0.89 km |

The complete [17 original port source identity and outcome table](WPI17_independent_PortWatch_source_geography_LOCODE.csv), including the TWO nonmatches, [PortWatch provenance and source-count receipt](WPI17_PortWatch_source_identity_receipt.json) and [independent source-row audit](independent_WPI17_source_audit.json) are retained. Every one of the bounded 31 source registry sample records, full source response status and geodesic calculation was independently rechecked, with the detailed provider rows archived ONLY in the temporary CI artifact (not republished as an IMF global port registry). All 17 WPI source identities have been kept even when provider matching fails.

## Original data-source contrasts restricted to these fifteen actual geographically verified original ports

We then joined the same **15 independently geoverified original WPI IDs** to their already archived **GLSN2018 14,052 modelled undirected co-service edges**, **UNCTAD historical port indices** and **Zenodo2020 separately verified AIS origin/destination records**. This is a source coverage crosscheck, not a causal butterfly outcome comparison. No geographically unverified source name was promoted.

| Support within the 15 geoverified original ports | Source-count result |
|---|---:|
| Original possible WPI port-pairs | **105** |
| 2018 GLSN modelled liner co-service candidate pairs | **6** |
| 2020 Zenodo AIS-inferred OD port pairs with source records | **10** |
| Distinct directed 2020 AIS original WPI port links | **17** |
| Indexed 2020 original OD triplet rows across all 15 source ports | **64** |
| Of those source OD rows connecting **different original island groups** | **47** |
| GLSN-modelled port pairs with at least one independent 2020 source OD record | **2** |
| Distinct GLSN original island-pair types | **6** |
| Distinct 2020 source AIS original island-pair types | **7** |
| Original island-pair types represented in both sources | **2** |
| Source WPI ports with ≥16 nonmissing UNCTAD 2006–2010 quarterly indices | **14 of 15** |

**The two port-pair identities appearing in both GLSN and 2020 source data** are:

- **Belfast ↔ Southampton, WPI 28–1143:** GLSN2018 modelled co-service source weight 800 (NOT ship count), with **15** later Zenodo2020 AIS-derived port OD source rows.
- **Santa Cruz de Tenerife ↔ Puerto del Rosario, WPI 428–430:** GLSN2018 source weight 4,606 (NOT ship count), with **10** later 2020 source OD rows.

Among the six 2018 modelled links there are also Manila–Haikou, Port Louis–Colombo, Port Louis–Toamasina, and Southampton–Colombo with no matching 2020 source rows in this strictly selected original port cohort. That does **NOT** show these links disappeared, since modelled network co-service and a selective Jan–Oct 2020 AIS inferred voyage source are different statistical observables, subject to different windows, vessel coverage and port matching.

All **105 original WPI verified-port pairs** are archived in [WPI15_geoverified_all105_GLSN_UNCTAD_AIS_port_pair_comparison.csv](WPI15_geoverified_all105_GLSN_UNCTAD_AIS_port_pair_comparison.csv), with original WPI/PortWatch IDs, LOCODE, original island IDs, historical UNCTAD quarters/index means, exact GLSN source weight and GC only when present, and observed 2020 source rows *by direction*. Missing/zero link entries explicitly say **unobserved source is NOT no shipping**. All original **1,275 island-pair** comparisons are in [original51_all1275_geo_verified15_source_island_pair_observation.csv](original51_all1275_geo_verified15_source_island_pair_observation.csv), keeping the original 51-island denominator rather than discarding unmatched islands or inventing absent routes.

The [source graph receipt](WPI15_3_source_verified_port_graph_audit.json) and [independent source-only numerical audit](independent_WPI15_3source_graph_audit.json) independently reconstructed all **105** port pairs, original GL/UNCTAD/AIS rows, 17 observed source directions and all **1,275** original island pairs. The standalone [reproduction code](../../scripts/summarize_original_WPI15_geoverified_maritime_sources.py) requires only archived provider data; it does not require another API call.

## What these data establish — and what remains untested

**Established:** the original selected 17 WPI candidate harbor IDs can be checked against a fourth, geographically independent port registry. A source-verified subset of **15 ports / 14 original islands** has independent port IDs, recorded official registry coordinates and LOCODE strings. Two original port-pair identities have evidence from distinct GLSN2018 modeled co-service and 2020 AIS-inferred OD sources, with the historical UNCTAD port connectedness series available for most of these harbors.

**Not established:** AIS-derived direct plant seed or butterfly transport, 2001–2010 shipping routes, true butterfly arrivals/feeding/reproduction/colonization, equivalence between the 2018 GLSN model weights, 2020 recorded AIS trips and UNCTAD's index values, or a completed global transport network. This source set is **enriched for ports already present in all three datasets**, far from a representative sample of 194 original ports/51 islands. For comparison, the prior 2020 real OD source identified only 12 out of 1,275 original island pairs, and the 27-nearby-port overinclusive sensitivity only 28 hypothetical pairs; no full 51-island causal real-route model is justified.

No original GEB manuscript or prior macroecological result was changed. Earlier positive global +54.9% *potential* butterfly host resource opportunity, later butterfly AIS prediction increments, adverse earlier plant/butterfly source metrics and PR137 local AIS graph-reassignment falsification remain separate and intact. **No paper submission.**
