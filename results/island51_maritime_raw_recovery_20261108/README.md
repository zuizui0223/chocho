# Original GLSN2018 and historical UNCTAD PLSCIQ sources recovered

**Data acquisition and source-identification audit completed on 2026-11-08. No butterfly regression fitted. No manuscript submitted.** The positive and negative biological results already in chocho are unchanged.

## 1. GLSN2018 actual original publisher files recovered

Source: [Network-Maritime-Complexity/GLSN-efficiency-vulnerability](https://github.com/Network-Maritime-Complexity/GLSN-efficiency-vulnerability), MIT license. The source publisher Git blob IDs were verified against bytes downloaded directly on GitHub Actions:

- Original port 2018 XLSX: Git blob SHA-1 5486d8c204a7c49b375a14762fe5cb36c7877d96, 34,675 bytes.
- Original Edge 2018 XLSX: Git blob SHA-1 18098cab432ddbcd7e0a7091fc1c4cedbda40fb7, 1,018,504 bytes.

Both exact original XLSX sheets were extracted via stdlib ZIP/XML to permanent, MIT-attributed source CSV tables in this PR.

| Publisher source sheet | Header columns | Data rows |
|---|---|---:|
| [2018 global ports](GLSN2018_ports_original_sheet_1.csv) | port, ISO3, Country Name | **887** |
| [2018 modelled weighted network edges](GLSN2018_edges_original_sheet_1.csv) | port1, port2, GC, weight | **28,104** |

**Every edge is present in both directions.** The complete publisher edge list contains **14,052 unique unordered port-name pairs**, each represented exactly twice with **no inverse-row mismatches in GC or weight** and no unknown source port-name endpoints. The edge is a *liner-service modelled co-membership connection*, NOT a proven direct physical vessel voyage or a number of ship trips. The original weight cannot be treated as butterfly transport volume.

The original Excel bytes, extracted CSV source rows, cryptographic hashes and environment are retained in the [source acquisition GitHub Actions artifact](https://github.com/zuizui0223/chocho/actions/runs/42430381585), subject to artifact expiry; the extracted source tables remain in the repository with [original-source receipt](GLSN2018_recovered_source_receipt.json).

## 2. Matched to the original WPI194 island-port universe, without false geographic confirmation

Original source: 194 WPI2019 biological ports, 182 marine-qualified terminals, 51 original landmass groups (original WPI Git blob a4681c3e3c1406d3a37f51464a298286e9d0ee4f). The GLSN port sheet contains a country ISO3 code and a port NAME but **no original port coordinate**. Only globally unambiguous GLSN port names and one-to-one country/name WPI matches were counted as source **candidate identities**, explicitly NOT independent location-verified port identities.

| Cross-source finding | Audited count |
|---|---:|
| Original WPI194 ports with one-to-one GLSN ISO3→ISO2 country+name candidates | **34** |
| Original island/landmass groups among those candidate ports | **18** |
| GLSN modelled unordered port pairs among those 34 WPI candidates | **37** |
| Original unordered island pairs touched by those modelled candidate edges | **10 of 1,275** |
| Original WPI194 ports with a one-to-one historical UNCTAD country+name candidate | **40** |
| Original WPI ports appearing in GLSN candidate, UNCTAD candidate AND **Zenodo2020 independently verified WPI port endpoint** | **17** |

Source-backed [all-194 port/country/name/UNCTAD/AIS joint ledger](GLSN2018_UNCTAD2006_2010_Zenodo2020_originalWPI194_candidate_join.csv), [37 GLSN candidate WPI modelled weighted edges](GLSN2018_originalWPI_source_candidate_model_edges.csv) and [all 1,275 original island-pair modelled candidate edges](GLSN2018_original51_all1275_modelled_island_pair_candidates.csv) are preserved. All rows explicitly mark GLSN's missing geographic identity confirmation. An absent modelled edge does NOT mean no actual shipping, and no additional voyages are added to the previously independently verified 2020 Zenodo AIS source.

### Independent 2018-model versus 2020-AIS source overlap

Within the **17 WPI ports** shared among GLSN name+country candidates, UNCTAD PLSCIQ name+country candidates and independently AIS2020-verified WPI ports, the original GLSN model has **6** candidate unordered port-pair edges. Zenodo2020 has **15 different unordered WPI port-pair original OD links** within these 17 ports (25 directed links). The overlap is exactly TWO WPI port-pair identities: **28–1143 Belfast–Southampton** and **428–430 Santa Cruz de Tenerife–Puerto del Rosario**. These are cross-source compatibility examples, not independent proof that a modelled liner connection is a non-stop observed cargo voyage; the original source windows and coverage differ.

All 14,052 original GLSN edges, 194 WPI joins, 37 candidate links and 1,275 island pair entries were reconstructed in [independent source audit](verify_GLSN2018_original_WPI_source.py); [source run](https://github.com/zuizui0223/chocho/actions/runs/42430798011).

## 3. Exact historical UNCTAD/DBnomics PLSCIQ source vintage recovered

The 2026-10-30 archived port-level liner shipping connectivity index source was recovered from the public DBnomics / UNCTAD API. All THREE independently re-fetched original response SHA256s match the previous acquisition audit:

| Original source series offset | Exactly matching original and newly recovered SHA-256 |
|---|---|
| offset 0: 500 series | 29a343bb6b1113b1699a9ca1821952cce1cb8f150831921da68e248e78fb3e1b |
| offset 500: 500 series | c82e2d2119b0f4b03e4afb18c04f9c4442f19a2bfce96efbc8d30a87aed42933 |
| offset 1000: 300 series | 7eb1eaf4349008b784d9e6d2bd834cadf29c51fa493c2cdd2f83e2d767978255 |

Publisher vintage directory hash: 5cfa0ae00d0c42645dae30be6a5fad068dc35c5f. Confirmed **1,300 original port series, 81,900 port-quarter values, 57,312 nonmissing**, 2006-Q1 through 2021-Q3, historical index base Maximum Q1 2006=100. **This is exactly the OLD version**, not a changed base or partial reconstruction.

All three original source JSON page responses are in the [30-day Actions source artifact](https://github.com/zuizui0223/chocho/actions/runs/42430381585), not redistributed as a complete public provider dataset while upstream rights remain to be independently established. Source [historic version/hash audit](UNCTAD_PLSCIQ_source_recovery_receipt.json) and [WPI194 original country/name index history candidates](UNCTAD_PLSCIQ_originalWPI194_name_country_candidates.csv) are permanent in the PR. Exactly **40** original WPI ports have unambiguous country/name candidate series; **33 ports from 21 original island groups** have 16 or more indexed quarters across 2006–2010. The provider has no independent original port coordinates, so matches remain *source country/name candidates*. Missing values remain missing rather than zero.

## Scientific boundaries and the next actual data need

The GLSN2018 network is a modelled liner-service co-membership graph, UNCTAD PLSCIQ is a port-level index and Zenodo2020 contains inferred AIS departure/arrival records. These sources represent three different measurement types, not three identical independent ship-voyage inventories. The retrieved data do NOT determine historical shipborne introduction of plants, host ingestion, butterfly actual arrival, successful colonization or causal stepping-stone movement.

The new historical index is dated **2006–2010** and is therefore more useful for temporal source-hypothesis work, but its 33 source port names/21 original island groups lack independent port-area confirmation. The original candidate-host Jaccard pooled 2001–2024 may postdate the later butterfly response. No original GEB manuscript was changed and no biology model has been fit in this recovery branch. Earlier positives (global +54.9% *potential* host portfolio opportunity; later butterfly AIS predictive increments), unfavorable early plant/butterfly metrics and PR137 local-AIS relabeling sensitivity remain unchanged.

**Current remaining task:** independently validate GLSN and UNCTAD source port identities using a port-ID/coordinate authority, before promoting any modelled network or connection index to a source covariate in the original 51-island response models.
