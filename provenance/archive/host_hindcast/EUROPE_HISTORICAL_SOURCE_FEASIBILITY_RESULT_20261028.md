# European 39-butterfly hindcast: historical GBIF availability PASS, observation process confounded

**2026-10-28 — complete original-source source-feasibility audit. No climate/host hindcast fitted yet; no ecological expansion or warming effect shown.**

## Precommitted source gate

Original chocho 239 accepted GBIF source identity manifest `sha256:18f84371cf2a9a71ca06202fac2246d9aa3df0465bd973d06b4e9e8bfa58e8d2`; 39 preselected European butterfly species. *Lasiommata megera* was `TAXON_UNRESOLVED` in the ORIGINAL manifest and remains in the 39 denominator but **not** reidentified using a newer live key; 38 original keys were queryable. A 1970–1990, 1991–2005 and 2010–2024 3-period × 38-key source-only GBIF pilot **completed 114/114 queries**, each requiring both `COUNTRY` and `BASIS_OF_RECORD` facets from a `limit=0` response (zero occurrence bodies retrieved).

**Original fixed information gate PASS**: a species qualified in a country only if ≥20 original coordinate-valid records dated 1970–1990 and ≥50 dated 2010–2024 in that SAME country. **37/39** fixed original names qualified in ≥1 country; all **8/8** countries had some qualifying species; fixed minimums were 25 species, 5 countries. The two nonqualifiers: unresolved original GBIF taxon *Lasiommata megera*, and *Boloria napaea* (no qualifying paired-country counts). These are source qualification outcomes, **not** biological absences.

## Sources alone reveal severe changing observation geography

The counts are **GBIF original indexed occurrence records**, often multiple visits of one individual/site, with different era lengths (21 vs 15 calendar years). They must **never** be read as butterfly abundance. Country is derived from the occurrence source; the records are not constrained to the identical monitored site in both time periods.

| Country | Eligible species | 1970–1990 GBIF index records | 2010–2024 index records | Share of early 8-country records | Share of late 8-country records | Raw records per-calendar-year multiplier (not biological) |
|---|---:|---:|---:|---:|---:|---:|
| GB | 28 | 1,130,062 | 11,635,336 | 91.81% | 56.33% | 14.4× |
| FR | 37 | 52,059 | 2,874,905 | 4.23% | 13.92% | 77.3× |
| ES | 33 | 25,955 | 341,863 | 2.11% | 1.66% | 18.4× |
| NL | 24 | 9,417 | 3,904,965 | 0.77% | 18.91% | 580.5× |
| SE | 29 | 6,746 | 719,748 | 0.55% | 3.48% | 149.4× |
| DE | 32 | 4,305 | 1,045,594 | 0.35% | 5.06% | 340.0× |
| IT | 18 | 1,379 | 111,795 | 0.11% | 0.54% | 113.5× |
| PL | 18 | 952 | 21,818 | 0.08% | 0.11% | 32.1× |

Fixed 8-country source records: **1,230,875** (1970–1990) versus **20,656,024** (2010–2024), much of this change reflecting sampling/digitization/reporting differences. GB accounts for **91.81%** of all early and **56.33%** of late fixed-cohort source records. Aggregating these eight countries would conflate a major shift in where people recorded butterflies with any putative distribution front change. Netherlands and Germany also show enormous reporting changes; do not claim they acquired 580× or 340× more butterflies.

### Original worldwide recording basis, same 38 keys

| GBIF source-recording era | Worldwide indexed records for 38 species | Human observation % | Preserved specimen % |
|---|---:|---:|---:|
| pre_1970_1990 | 1,570,100 | 91.75% | 7.19% |
| mid_1991_2005 | 7,439,342 | 98.91% | 0.92% |
| post_2010_2024 | 26,204,973 | 99.13% | 0.17% |

These recording-basis facets describe **worldwide** original 38-species records, not the UK/EU subset and do not give country×basis joins. Museum-versus-human classification is a necessary diagnostic but insufficient to infer equal survey effort.

## Meaning and next data gate

**Positive:** the original 1970–1990 butterfly record availability gate is *not* the limiting factor. We can select a historically paired survey panel without searching for more European butterfly names. **Not yet positive:** relative predictive skill of host plants, botanical range tracking, climate-front shift, actual colonization, and human-managed host contribution are all untested.

The next step must establish **historically repeated actual monitoring effort within the same spatial units**, especially in UKBMS transects. A GBIF grid with records for other butterflies is only an observation proxy; true absence or local colonization needs survey visits with zero detections and fixed effort. eBMS detailed transect histories have access terms and may require an explicit licensed request. The UKBMS GBIF occurrence source has earlier historical coverage, but its field-level absence and land-management history must still be verified before using it as a full detection history.

If enough repeated UK monitoring units exist, score a time-forward new-grid occurrence/occupancy prediction using pre-1991 or other predeclared early-window climate from E-OBS (1950–2025) and source-dated host plant observations. Compare climate/effort/land-use-only versus source-verifiable native wild hosts, wild introduced hosts, and independently verified crop/ornamental/grass-maintained hosts. **WCVP native/introduced is not cultivated status**, and pre-2010 host presence must not be reconstructed from post-2010 observations. Include matched nonhost plant, land cover and survey-effort negative controls. Only if out-of-time host models improve against these controls can future CMIP6 be interpreted mechanistically.

## Exact source/provenance

- [Source-run GitHub Actions #40631110936](https://github.com/zuizui0223/chocho/actions/runs/40631110936) **SUCCESS**, 114/114 real source queries, 5 source-independent tests; [source artifact](https://github.com/zuizui0223/chocho/actions/runs/40631110936/artifacts/12790973591) original ZIP SHA256 `a1b47cf7990dbec41e5d9eb3b98222f8c0fc7d97c01f5941ce0523e97a69e1d6`.
- Original complete derived JSON SHA256 `3f0a819752612211158dc9f0031181254ae35f14fa67fdbfd993333b7e0ba467` (individual per-query source checksums, source URLs and exact fixed GBIF taxon IDs, no original butterfly coordinate records).
- `provenance/archive/host_hindcast/EUROPE_HOST_MEDIATED_HINDCAST_PROTOCOL_20261028.md` was frozen before all historical source outcomes; one original unresolved GBIF name retained as explicitly unscorable after a source-identity failure with no historical outcomes inspected.
- Prior: [Rashid et al. 2026](https://doi.org/10.1111/ddi.70245), 118 Europe-centred butterflies and static-host future climate range effects; that broad idea is prior art. The possible new question is demonstrated historical **human-maintained host** assistance vs wild-host lag under independent effort-adjusted validation.

Keep exploratory Draft; original GEB `main` and manuscript remain unchanged.
