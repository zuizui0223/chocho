# Subnational GBIF geography and developmental-stage labels are available, but ecological realization remains unresolved

**20 October 2026; all 48 positive species × country cells, no selection by stage outcomes.**
This source-locked macroecological audit retrieved compact facets, not
individual records or raw coordinates; all 48 results exactly match the
previously frozen GBIF country counts. Successful
[acquisition](https://github.com/zuizui0223/chocho/actions/runs/39414408960)
and [independent archive replay](https://github.com/zuizui0223/chocho/actions/runs/39415280084).

## Geographic resolution and source stage information

| Field | Validated quantity |
| --- | ---: |
| Source-positive butterfly species × countries | **48** (20 species; 23 countries) |
| Indexed GBIF occurrence records | **756** |
| Assigned to GADM level-1 administrative divisions | **736/756 (97.35%)** |
| Unassigned to GADM1 | **20/756 (2.65%)** |
| One observed administrative division | **18/48 positive cells** |
| Median GADM1 units, all cells | **2** |
| Cells with ≥5 indexed records | **20** |
| Median GADM1 units among ≥5-record cells | **5.5** |
| Provider-submitted lifeStage strings | **500/756** |
| establishmentMeans values | **0/756** |

The two previous candidate cases each span **six GADM1 units**:
*Aguna asander* in the Dominican Republic (9 source records) and
*Papilio demoleus* in Israel (86). Counts from different administrative
units are not independent populations. The 18 apparently single-division
cases disproportionately include sparsely recorded units; a single
record provides no evidence of narrow range occupancy.

## Exactly what the 500 life-stage labels mean

| Raw label | GBIF records |
| --- | ---: |
| Adult | 423 |
| Imago | 17 |
| Larva | **33** |
| Caterpillar | **1** |
| Pupa | **4** |
| Egg | **1** |
| Unknown | 21 |
| **Total** | **500** |

Exactly **39 raw early-stage labels** are present in **11 species×country
cells and only three butterfly species**: *Papilio demoleus* **27**,
*Colobura dirce* **11**, *Anteos clorinde* **1**. They have **not**
been independently reviewed against individual images, original provider
observations or identified larval host plants; they are **not** 39
independent successful feeding, survival or recruitment events.

| Butterfly | Country | Introduced-added host resource at country scale | Indexed records | Early-stage labels | Mapped administrative divisions |
| --- | --- | --- | ---: | ---: | ---: |
| *Anteos clorinde* | CU | あり | 34 | 1 | 7 |
| *Colobura dirce* | CU | なし | 8 | 1 | 4 |
| *Colobura dirce* | DO | なし | 50 | 6 | 16 |
| *Colobura dirce* | HT | なし | 2 | 1 | 2 |
| *Colobura dirce* | SR | なし | 16 | 2 | 6 |
| *Colobura dirce* | VE | あり | 56 | 1 | 9 |
| *Papilio demoleus* | CU | あり | 35 | 1 | 11 |
| *Papilio demoleus* | DO | あり | 199 | 19 | 27 |
| *Papilio demoleus* | HT | あり | 8 | 2 | 3 |
| *Papilio demoleus* | IL | なし | 86 | 3 | 6 |
| *Papilio demoleus* | SY | なし | 26 | 2 | 5 |

The added-host-country group has **24 such source labels across 467
records** (327 records with any life-stage label); the non-added group
has **15 across 289** (173 with any life-stage label). Raw rates (~5.1%
in either group) are **not a null-effect finding**: different species,
observers, platforms and habitat/host catalogue coverage are mixed
and the records are not independent sampled juveniles. The underlying
country host flags do not establish local host species, host plant
arrival time or larval success.

## Implication for the macroecological study

The new result is a **real, source-verified observation-data constraint**:
subnational administrative geography is generally populated, but
biological establishment metadata is wholly missing and potentially
juvenile-stage evidence is narrowly concentrated. The prior 229-taxon
country association (+3.91 pp; 95% interval −6.42 to +14.94) still
has only **21 informative strata** against a fixed minimum of 50.
The source cannot demonstrate that newly relocated host plants
produced butterflies.

A meaningful next test must match plant and butterfly locations **within
subnational areas**, identify which plant was actually used, then
verify host-dependent development or population recruitment. Further
retesting these same country detections without new local host data
will not resolve the biological mechanism.

## Reproducibility and immutable separation

- [Complete safe GADM1/lifeStage source facets](gbif_positive48_gadm_lifestage_facets_20261020.json)
- [Source-verified derived tables](gbif_positive48_subcountry_stage_summary_20261020.json)
- [Frozen source-blind protocol](../../protocols/GBIF_SUBCOUNTRY_POSITIVE48_QC_20261020.md)
- [Executable source replay](../../../scripts/summarize_gbif_positive48_subcountry_20261020.py)
- GEB frozen manuscript remains unchanged; all files belong to the
  separate exploratory PR #51 and not its anonymous submission snapshot.
