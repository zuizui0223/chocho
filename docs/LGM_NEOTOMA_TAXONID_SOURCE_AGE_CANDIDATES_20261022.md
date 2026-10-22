# Neotoma: taxon-ID-specific fossil source and exact 18–24ka sample metadata (2026-10-22)

**Scientific state: candidate fossil occurrence metadata only.** An exact botanical taxon-ID match, dataset type and numeric age field are NOT independent confirmation of calibrated sample age, fossil pollen species identification or 21ka plant occupancy. Neither source insect genotypes nor older negative true-host-vs-other-plant climatic comparisons were reopened.

## Source-provenance failure and repair

The previous 15 plant names ×2 fossil proxy API queries had 12 HTTP 200 returned payloads but only **two distinct payload SHA256s**, irrespective of requested species within pollen or plant-macrofossil categories. [v11 source negative control #39722583906](https://github.com/zuizui0223/chocho/actions/runs/39722583906) verified that generic `/v2.0/data/occurrences?taxonname=...` ignored even the *impossible* botanical name; **never count those as plant fossil occurrences**.

Official taxon registry `/v2.0/data/taxa` accurately distinguishes *Alnus glutinosa* (Neotoma ID 1299), *Betula pendula* (3461) and the impossible sentinel. ID-specific `/v2.0/data/taxa/{taxonid}/occurrences` returns distinct source occurrence records; verified by [v12 #39723327745](https://github.com/zuizui0223/chocho/actions/runs/39723327745). Its two 20-item/pages per species were **not age-targeted**, and their 80 returned rows' central ages did not lie in 18–24ka; this does NOT imply old fossils absent.

## Validated *source metadata* window using numeric taxonid query

[v13 #39723924751](https://github.com/zuizui0223/chocho/actions/runs/39723924751) ran only two fixed plant taxa ×2 source proxy types, `taxonid={1299,3461}`, `datasettype={plant macrofossil,pollen}`, `ageyoung=18000`, `ageold=24000`, `limit=30`, plus sentinel taxon ID 9999999 and 100–200 BP age controls. All source classifier tests and original v11/v12 SHA checks passed.

| Source request | Source metadata records returned | Validation scope |
|---|---:|---|
| *Alnus glutinosa* pollen | **27** | 27 rows match Neotoma ID1299, exact database taxonname, pollen datasettype, 18–24ka numeric `age.age` |
| *Alnus glutinosa* plant macrofossils | 0 | No records **in this bounded query**, NOT fossil absence |
| *Betula pendula* pollen | 0 | Same limitation |
| *Betula pendula* plant macrofossils | 0 | Same limitation |
| Impossible taxonid 9999999, pollen, 18–24ka | 0 | Negative filter integrity control |
| Alnus ID1299, pollen, 100–200 BP | 2, distinct source IDs | Temporal source query response changes under different ages |

Among **27 Alnus pollen rows**, there are **26 unique occurrence IDs** at **6 distinct site IDs**, each corresponding to one of **6 dataset IDs**. Dates of the reported central `age.age` values range **18,241–23,611 BP**. One duplicate occurrence ID must not be counted as a separate specimen. Only **1/27** occurrence rows supplies non-null `ageolder` and `ageyounger`; other central age values have no source uncertainty interval available at this API level.

Site contributions, not six independent inferred colonization events:

| Neotoma site ID | Rows | Unique occurrence IDs |
|---|---:|---:|
| 9834 | 9 | 9 |
| 23855 | 7 | 7 |
| 15750 | 6 | 5 |
| 3392 | 3 | 3 |
| 14474 | 1 | 1 |
| 14485 | 1 | 1 |

These records' `site.location` fields are polygon or point geometries, but their coordinates are **sediment/sample-site metadata**, not source-confirmed plant organism locations. Transported or morphologically genus-level pollen might lead to inappropriate binomial assignment. The 18–24ka central age is a **database-model age field**, and `ageolder`/`ageyounger` are mostly missing. No independent age–depth model, sample chronology, proxy taxonomic specificity, publication or depositional origin has been reviewed.

Thus: **27 source-level candidate age/taxon/proxy records, 0 confirmed independent species-level LGM plant occupied locations at present.** The two taxa are selected for API technical checks, not a representative 416-host panel. Queries are bounded and do not establish fossil absence for taxa returning zero.

## Next decisive ecological dependency

Resolve the six source datasets' site chronologies, age-depth controls and original pollen taxonomic identification through the official Neotoma `/v2.0/data/datasets/{datasetid}/chronologies`, individual dataset metadata/downloads and publication records. Classify source confidence at **dataset-site**, not pseudoreplicated occurrence-row grain, and never include records whose pollen identification cannot support actual species identity in fossil-constrained LGM native-host distribution models.

Only independently validated dated **plant macrofossils/species-diagnostic pollen**, with depositional/geographic provenance, could begin constraining glacial host occupation. That remains different from proving ancestral host-mediated insect genetic divergence. The original TTF v0.5 and chocho v0.8 climate-matched alternative-plant **negative** results stay frozen.
