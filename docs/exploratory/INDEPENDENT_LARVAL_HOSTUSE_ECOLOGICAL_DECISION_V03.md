# Independent larval host-use evidence: ecological decision v0.3

**Date:** 2026-10-08. **Status:** Exploratory. Main GEB manuscript unchanged.

## Decision

The 11-species GloBI pilot independently documents potential larval feeding on plants whose distributions are introduced in the focal WGSRPD3 regions. It **does not** establish that butterflies colonize areas lacking native host plants. The original prespecified coverage gate FAILED (45 records but only two butterfly species in introduced-only resource cells; threshold three).

## Evidence and source-level correction

Sources: [GloBI v0.2 run 37712115731](https://github.com/zuizui0223/chocho/actions/runs/37712115731), [45 iNaturalist originals audit run 37712707184](https://github.com/zuizui0223/chocho/actions/runs/37712707184).

- 623 deduplicated geographically mapped, dated potential larval feeding records, of which 153 involve known HOSTS food plants that are introduced in their WGSRPD3 region.
- 45 of these fall in a butterfly-region cell where the original frozen HOSTS-WCVP union has no known native host opportunity. The 45 have unique iNaturalist observation IDs, all retrievable and with photographs; all marked captive=false.
- **44/45** match the original butterfly taxon, year, explicit larval annotation and an iNaturalist field asserting feeding (Eating, Feeding on, or Interaction->Herbivore of). **Observation 338738** is presently identified as a *Selasphorus* bird instead of the GloBI *Danaus plexippus* label; the event year also disagrees. Exclude it: https://www.inaturalist.org/observations/338738
- Of 44 source-consistent records, 3 report geographic uncertainty >1km and 30 have missing accuracy estimates. Photo presence and research-grade taxonomy do not independently prove consumption; field entries are participant assertions. No host quality, density or fitness is measured.

## Origin matters

Of the 44 source-consistent records, **43 concern butterflies beyond their original biogeographic range**: *Pieris rapae* in New Zealand and *Danaus plexippus* in New Zealand, Australia or Hawaii. Some monarch colonization may involve natural long-distance dispersal; do not label every such butterfly as a human introduction.

The one original-native-range case is the Puerto Rican monarch subspecies *Danaus plexippus portoricensis*, [observation 109870256](https://www.inaturalist.org/observations/109870256). Its 2022-03-29 iNaturalist record is research grade, wild, 43m positional accuracy, and two photographs, with a larval annotation and a Feeding on field referencing *Calotropis procera* taxon ID 120917. The associated separate plant record [109865691](https://www.inaturalist.org/observations/109865691) is research-grade *C. procera*, photographed on the same date at a nearby location (8m stated accuracy). Both records reciprocally link each other. The provenance was verified through [GitHub Actions 37712900417](https://github.com/zuizui0223/chocho/actions/runs/37712900417). This verifies linked taxon claims, **not photographically adjudicated ingestion**.

This is NOT a new local host-use discovery. A pre-existing Puerto Rican study by López, Santiago and Puente-Rolón compared monarch larval response to *C. procera* and *Asclepias curassavica*, reporting a nonsignificant tendency toward *Calotropis*: https://www.arecibo.inter.edu/wp-content/uploads/portal/pdf/inter_scientific_05.pdf

Kew POWO classifies *Calotropis procera* **introduced** to Puerto Rico and *Asclepias curassavica* **native** to Puerto Rico:
- https://powo.science.kew.org/taxon/1004515-2
- https://powo.science.kew.org/taxon/94213-1

However, the **frozen 42 accepted HOSTS records for Danaus plexippus omit the documented Asclepias curassavica host link** while including Calotropis procera. Thus the Puerto Rico cell's introduced-only classification is an inventory artifact: a known native host exists but the corresponding butterfly-host link is missing from the frozen reconstruction. Do not treat an empty reconstructed native host union as a biological absence.

## Broader ecological implications

The strict within-butterfly-family and within-host-plant-family identity null remains a modest structural association (+0.0105 Jaccard over its null, one-sided p=0.030, independent-seed p=0.026). Most host-linked climate selectivity was already present in the native reconstruction (+0.0405 residual); the additional globalization signal is much smaller.

These observations support that some potential host links are used in the field, but **do not independently validate the broad claim that plant introductions cause native butterfly range filling or improve performance**. The 11-species independent observational coverage gate failed, and one apparent exception is misclassified due to native host-link incompleteness.

See Allf et al. (2026) PLOS Biology https://doi.org/10.1371/journal.pbio.3003988 regarding iNaturalist feeding annotations versus directly confirmed feeding.

**Decision:** Do not promote an independent realized butterfly ecological paper from these 45 cases. Do not loosen thresholds post hoc. A future independent panel should sample novel butterfly species before seeing observational outcomes and verify all native host alternatives. To test ecological preference or fitness, obtain experimental larval performance or observation denominators, not interaction-positive presences alone.

Main GEB manuscript unchanged.
