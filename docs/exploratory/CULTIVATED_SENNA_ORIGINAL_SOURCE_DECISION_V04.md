# Cultivated Senna original-source audit: completed decision v0.4

**2026-10-08. Status: EMPIRICALLY OBSERVED CULTIVATED HOST REPORTING; DIRECT ECOLOGICAL GATE NOT MET.** Branch `analysis/resource-homogenization-v01` only, main GEB manuscript unchanged.

## Exact reproducible original-source observations

[GitHub Actions run 37728156494](https://github.com/zuizui0223/chocho/actions/runs/37728156494) completed successfully. Artifact `butterfly-cultivated-senna-observer-larva-v02` contains `original_observer_and_larval_audit_v02.json` and source-assertion candidates.

- All **573/573** original taxon/date/captive/photo/coordinate source records from the 4-Senna Florida pilot were retrieved; zero metadata mismatches in the defined checks.
- **59/59** Florida *Senna polyphylla* photo-supported records marked cultivated remained source-consistent.
- The 59 records involved **56 unique iNaturalist observers**, maximum **two** records by the same observer (**3.39%**), **46 approximate 5-km grids**, **58 unique observer×grid pairs** and **10 distinct observation years** (2010–2025). The predetermined observer-diversity gate passed.
- A photographed/cultivated original observation is **not** a biologically validated individual planted shrub, garden, naturalized stand, larval feeding event, or fitness observation.
- The frozen WCVP geographic layer has no *S. polyphylla* Florida row; this remains a valid **horticultural exposure not captured in a naturalized-plant range database**.

## Independent butterfly use endpoint

The frozen species, date and life-stage query investigated original iNaturalist Florida photo observations annotated Larva:

| Butterfly | All queried exact larval records | Source fields on records | Explicit feeding or partner field linking to S. polyphylla |
| --- | ---: | ---: | ---: |
| *Phoebis sennae* | 338 | 28 | **0** |
| *Phoebis philea* | 577 | 70 | **0** |
| **Total** | **915** | **98** | **0** |

Both species' source queries were exhausted (not capped). **17** exact annotated larval observations were reported within 2 km and ±365 days of one of the cultivated *S. polyphylla* records; these are strictly incidental spatiotemporal co-records and **not trophic links**. Absence of a feeding-field link cannot establish absence of actual use. The original photos were not adjudicated for feeding.

The separate [Koptur et al. (2024), Insects](https://doi.org/10.3390/insects15020123) rearing study demonstrates larvae can emerge as adult butterflies on cultivated *S. polyphylla* at its study sites. Both its native and nonnative *Senna* were cultivated. The published source's parasitism table has internally conflicting cell and grand totals; survival/whole-cohort effect signs are not identified with 176 unresolved fates. These already-published empirical patterns must not be relabeled as a new ecological discovery.

## Decision

**Stop the garden-observation shortcut to a new butterfly ecological mechanism paper.** Source independence is demonstrated for the **botanical occurrence reporting**, not for consumer response. The direct-field larval-usage gate failed and there is no replicated fecundity, confirmed larval-to-adult success, or enemy-exposure contrast for cultivated versus naturalized host plants.

Retain the horticultural inventory blind spot as a valid limitation and potentially a separate public-data/design contribution. For genuinely new biological inference, prioritize independently observed *host performance* and its association with different naturalized host redistribution profiles across multiple plant species, with taxonomy, native range and evidence class controlled. The independent monarch 127-plant quality × pinned WCVP globalization test is frozen under `MONARCH_HOST_QUALITY_GLOBALIZATION_PROTOCOL_V01.json`; its workflow is `butterfly-monarch-host-quality-globalization.yml`. Neither test demonstrates local butterfly colonization, realized competition or fitness.

No change to `main` or GEB submission.
