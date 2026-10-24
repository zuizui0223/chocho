# GloBI source-provenance correction for the apparent georeferenced larval-feeding records

**24 October 2026 — source definition audit, separate from fixed taxon sampling.**

## Evidence verified from primary data-source repositories

1. **`globalbioticinteractions/menares2025-herbivory` is literature-derived and geographically coarsened.** Its original [GloBI configuration](https://github.com/globalbioticinteractions/menares2025-herbivory/blob/main/globi.json) describes historical German butterfly–plant relationships compiled across whole federal states and long time windows. Its two source tables explicitly assign **one locality coordinate per broad region** (Baden-Württemberg `9,48.5`, Brandenburg `13.01582,52.45905`) and coarse `eventDate` periods **1975/2005** or **1960/2018**, while labelling insect larval stage. These are source mapping rules, **not independently measured exact caterpillar encounter coordinates/dates**. The source [README](https://github.com/globalbioticinteractions/menares2025-herbivory/blob/main/README.md) describes the compilation as literature based. Consequently any GloBI query row from this dataset that passes a numerical `latitude/longitude` and `source_specimen_life_stage=larva` filter is **NOT proven a contemporary site-specific observed feeding event**.
2. **`globalbioticinteractions/inaturalist` does contain citizen-contributed interaction claims.** Its [repository](https://github.com/globalbioticinteractions/inaturalist/blob/main/README.md) states that GloBI indexes observations from iNaturalist contributors. Its [interaction type map](https://github.com/globalbioticinteractions/inaturalist/blob/main/interaction_types_mapping.csv) maps the free-submitted observation fields *Eating*, *Feeding on*, *Food Source* and other related labels into GloBI `eats`. Thus an `eats` record is a **source-contributed association** and may be useful, but requires the source observation, photographed life stage, exact plant match, and temporal provenance verification before it can be upgraded to directly verified field feeding.
3. GloBI source labels distinct from `HOSTS` do **not** guarantee independent raw primary observations. Citations can overlap historical HOSTS literature, and studies can use GloBI-supplied dates or centroids. Independent source is a **necessary but not sufficient** condition for new biological confirmation.

## Impact on interpretation of first 20-species receipt

The initial [v0.1 machine result](globi_larval_host_feasibility_20261024.json) captured **14** nominally strict plant + `eats` + explicit larva + numeric-coordinate + non-HOSTS-source rows in **5/20** response-blind focal butterfly taxa, but five taxa had truncated pages. The distinct source IDs behind those 14 rows were only **iNaturalist** and **Menares 2025 herbivory**. Since one of these sources maps regional centroids to literature links, the 14 rows **cannot be called 14 independently verified local larval feeding events**.

**Source-level gate verdict:** no new individual-level larval or adult recruitment effect is established, regardless of whether the fixed taxon pagination finds further GloBI associations. Page-corrected v0.2 counts, if obtained, will decide only the *pre-specified metadata-availability* gate; their interpretation ceiling remains potential candidate feeding claims.

## Exact next biological threshold

For a publishable macroscopic realized-use study, require occurrence-anchored interactions with source record IDs, actual dated and georeferenced juvenile encounters, explicit plant identification, and adequate independent species-region controls. Site-level host presence without verified caterpillar feeding is not an observed interaction; caterpillar feeding without survival-to-adult is not population establishment. Do not substitute artificial `date` midpoint or region centroids for actual field encounters.

The frozen `chocho/main` GEB paper and earlier exploratory disturbance analysis remain unchanged.
