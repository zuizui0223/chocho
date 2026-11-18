# Original-island GBIF Larva/Egg query for port-indexed candidate hosts

**14 November 2026 | Completed evidence-bounded exploration, not confirmation of establishment.**

## New source-verified finding

From **15** earlier-introduced-host-indexed butterfly × original-island cells, covering **11 GSHHG 2.3.7 landmass IDs**, a fresh 2011–2017/2018–2022 GBIF stage-facet acquisition returned **one positive species × island × period × stage cell**. There were **152 of 152** polygon-fragment × period × stage facet searches (all successful), producing **60 fully covered case-stage combinations**.

The positive cell was:

| Island source key | Butterfly | Time window and stage | Indexed stage occurrences | Raw GBIF evidence |
|---|---|---|---:|---|
| GSHHG landmass 951 (original source ID, not independently geocoded island name) | **Pieris rapae** (GBIF speciesKey 1920496) | **2018–2022, Larva** | **2** | GBIF 2814126623 (**2020-06-06**), GBIF 3772636962 (**2022-04-29**) |

The source API returned both as **HUMAN_OBSERVATION, lifeStage=Larva**, with exact speciesKey identity. The 2020 record reports coordinate uncertainty **65 m**, while the 2022 record has **no numerical uncertainty reported**. Both were indexed within the frozen GSHHG original-island polygon. They share GBIF contributing datasetKey **50c9509d-22c7-4a22-a47d-8c48425ef4a7**.

This establishes **indexed immature-stage observations in two separate calendar years** within that source island polygon. It is more specific than simply two later adult occurrences, but is **NOT proof of multi-generation population persistence, successful recruitment, or independent founding**. The observation data do not say which plant the larvae were feeding on. The original 2001–2010 port inventory for this *P. rapae*–island cell indexed **both regionally native and candidate introduced known larval host taxa**; attribution of the larvae to an introduced host is consequently unsupported.

### Entire 15-case source and observation boundaries

- All 15 had a pre-2011 indexed candidate introduced host near an original WPI port; 14 of 15 were in the **indirect author-selected review-related candidate migrant positive list**, but the one remaining candidate has **UNKNOWN migration status**, not evidence of nonmigration. No actual flight path, independent dispersal distance, or source population is verified.
- The other **14** case cells did not yield a stage facet under the selected query years and standardized GBIF lifeStage=Larva/Egg. This means **no qualifying index in the queried source**, not true ecological absence of immature stages, oviposition or breeding. GBIF lifeStage metadata are incomplete and vary strongly among data providers.
- Some islands required multiple nonoverlapping GSHHG fragments to preserve their source geometry. Each query preserved source checksum, fragment WKT hash, query-period, stage and speciesKey facet response hash.
- The first butterfly indexing in 2011–2017 and subsequent 2018–2022 re-indexing are GBIF observation windows, **not known introduction and colonization dates**.
- No source-data evidence traces the 2020 or 2022 larvae to named introduced host plants; no evidence establishes maritime transport or island-to-island stepping stones.

## Reproduce and audit exact provenance

Executed successfully in [source-verified workflow #43579732967](https://github.com/zuizui0223/chocho/actions/runs/43579732967). [Artifact of all 60 stage combinations, 152 individual query receipts, full source SHA256/geometry audit and GBIF record metadata](https://github.com/zuizui0223/chocho/actions/runs/43579732967/artifacts/14107767109).

The code is in scripts/acquire_port_host_immature_gbif_20261114.py and the workflow in .github/workflows/chocho-port-host-immature-facet.yml. GSHHG source ZIP digest was independently verified as 8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41 before querying.

**Scientific verdict:** One named butterfly on one selected island now has temporally separate original larval observations; the evidence does *not* identify the larval food plant or whether human-redistributed hosts enabled this species' island establishment. This positive case should be shown transparently as an illustrative stage witness and not counted as a proven mechanism or global replicated result.
