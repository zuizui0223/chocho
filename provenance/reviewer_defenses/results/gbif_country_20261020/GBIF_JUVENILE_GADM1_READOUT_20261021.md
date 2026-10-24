# Indexed juvenile-stage butterfly observations can be localized within countries, but not linked to larval hosts

**21 October 2026; complete exploratory, source-conditional GADM1 result.**
We used **all 15 exact nonzero raw life-stage label × butterfly × country combinations** from the previous, fixed 48-positive-cell GBIF index snapshot. No outcomes were selected after viewing the mapped administrative subdivisions.

Every stage-filtered GBIF count equals the original marginal count, and each returned GADM1 category is nested within the same butterfly's all-stage GADM1 count. The first complete retrieval saved source aggregate queries but exited due to a misplaced `flush` keyword during console logging. **No ecological computation failed**. [The second Actions run](https://github.com/zuizui0223/chocho/actions/runs/39449326117) downloaded and independently rebuilt **the original artifact** without another GBIF request; both source and replay are preserved.

## Main observational findings

| Verified indexed facet metric | Quantity |
| --- | ---: |
| Exact lifeStage-conditioned GBIF source queries | **15/15** |
| Source-positive species×country cases carrying juvenile-stage strings | **11** |
| Source-confirmed recorded life-stage strings (Egg, Larva, Caterpillar, Pupa) | **39** |
| Of those with first-level administrative code | **39/39 (100%)** |
| Distinct GADM1 subdivisions across all 11 species-country cases | **22** |
| Sum of distinct GADM1 subdivisions per species-country (allows repeated administrative codes across taxa) | **24** |
| Added-host-country source-flag juvenile strings | **24** across 5 source-positive cases |
| Non-added-country source-flag juvenile strings | **15** across 6 source-positive cases |

The 39 are raw provider lifeStage metadata and may include duplicate sightings, non-host observations and misidentifications. **They do not demonstrate successful larval feeding, recruitment, establishment or demographic changes**.

| Butterfly species | Country | Original country added-host flag | Juvenile-stage source labels | GADM1 units with juvenile labels | GADM1 units with any focal record |
| --- | --- | --- | ---: | ---: | ---: |
| *Anteos clorinde* | CU | Added | 1 | 1 | 7 |
| *Colobura dirce* | CU | Not added | 1 | 1 | 4 |
| *Colobura dirce* | DO | Not added | 6 | 3 | 16 |
| *Colobura dirce* | HT | Not added | 1 | 1 | 2 |
| *Colobura dirce* | SR | Not added | 2 | 1 | 6 |
| *Colobura dirce* | VE | Added | 1 | 1 | 9 |
| *Papilio demoleus* | CU | Added | 1 | 1 | 11 |
| *Papilio demoleus* | DO | Added | 19 | 10 | 27 |
| *Papilio demoleus* | HT | Added | 2 | 1 | 3 |
| *Papilio demoleus* | IL | Not added | 3 | 3 | 6 |
| *Papilio demoleus* | SY | Not added | 2 | 1 | 5 |

The strongest source-localized evidence consists of *Papilio demoleus* in Dominican Republic (**19 labels distributed across 10 GADM1 units**, from a country-level added-host flag) and Israel (**3 labels across 3 GADM1 units**, with no added-host flag in the original source catalogue). These are **not** matched treatment and control populations. An entire GADM province may still contain unrelated plants and butterfly locations, and the source life-stage annotations have not been independently adjudicated.

## Why the macroecological question remains unresolved

The original country contrast uses **229 butterflies** and 33,817,891 accepted indexed records but only **21 informative species × large-region strata**. Querying juvenile labels at the subcountry level improves the *location of existing observed metadata* but cannot add independent non-detection effort controls, plant locations or true host-use data. It therefore cannot lift the primary information gate or distinguish allopatry from finer-scale mismatch between butterflies and relocated hosts.

The right next ecological test is **an original HOSTS/WCVP accepted-plant-ID matched subcountry observation analysis**, first reconstructing which introduced hosts are proposed for the five added-host juvenile cases, then independently querying GBIF **plant** observations in those candidate GADM1 units. Apparent GADM-level co-detection still does not prove feeding or host-mediated recruitment. A source-matched absence/control background and plant/animal recording effort remain essential before fitting a macroecological effect model.

## Reproduction and provenance

- [Complete archived 15 conditional raw facet summaries, URL and SHA-256](gbif_juvenile_gadm1_conditional_20261021.json)
- [Compact cross-tab and species-country table](gbif_juvenile_gadm1_readout_20261021.json)
- [Frozen before querying 15 facets](../../protocols/GBIF_JUVENILE_GADM1_CONDITIONAL_20261021.md)
- [Source-exact Python implementation](../../../scripts/audit_gbif_juvenile_gadm1_20261021.py) and the [no-network replay](https://github.com/zuizui0223/chocho/actions/runs/39449326117)
- The original frozen GEB manuscript and blind review snapshot are unchanged.
