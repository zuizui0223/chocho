# Two missing *Castilleja* hosts: no change in broad-scale resource envelope

**Date:** 2026-10-08. **Status:** COMPUTED SOURCE-FROZEN, POST-HOC ONE-BUTTERFLY COMPLETENESS SENSITIVITY; NO NEW DEMOGRAPHIC FINDING. Primary GEB manuscript and PR #38 unchanged.

## Study and frozen execution

A previous exact-pinned source audit of *Euphydryas editha* found that two *Castilleja* spp. used by Taylor's checkerspot butterflies (Haan, Bowers & Bakker 2021, doi:10.1038/s41598-020-80413-y) were **not exact butterfly–plant interaction rows** in the original frozen HOSTS inventory, despite being legitimate locally documented hosts. We separately froze `EUPHYDRYAS_MISSING_CASTILLEJA_AUGMENTATION_PROTOCOL_V01.json` after observing their absence and the original 52/114 *Plantago*-specific geographic dependence. The added host set was fixed to *Castilleja hispida* and *C. levisecta*, with no genus expansion.

**Source-frozen successful workflow:** [GitHub Actions 37770418250](https://github.com/zuizui0223/chocho/actions/runs/37770418250), triggering commit `a93ff6f0769be6cd1254cf8e211b8a54a34cb82a`, artifact `chocho-euphydryas-two-host-completeness-sensitivity-v01` (ID 11548281153). Original HOSTS commit `808e0b869f9ec1adf8efff87cf6a395adda103e0` and original WCVP commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`.

## Exact botanical identities

| Independently documented missing host | Unique exact WCVP accepted ID | Native WGSRPD3 regions | Contemporary WGSRPD3 regions | Introduced-only |
| --- | --- | ---: | ---: | ---: |
| *Castilleja hispida* | `2705047` | 6 | 6 | 0 |
| *Castilleja levisecta* | `2705113` | 3 | 3 | 0 |

Both accepted taxa exist botanically, but their species-level links with *E. editha* were absent from the frozen butterfly host list. They were added solely in this diagnostic, not to global paper inputs. All of their WCVP level-3 resource regions were **already inside the union of the other 30 recorded host distributions**.

## Main numerical result

| Structural quantity for *E. editha* | Frozen 30-host map | 32-host augmentation |
| --- | ---: | ---: |
| Accepted known host species | 30 | 32 |
| Native resource WGSRPD3 regions | 131 | 131 |
| Contemporary resource WGSRPD3 regions | 245 | 245 |
| Contemporary minus native ('introduced-added') resource regions | 114 | 114 |
| Added resource regions dependent solely on *Plantago lanceolata*'s introduced distribution | 52 | 52 |
| Fraction of introduced-added regions uniquely mapped to introduced *Plantago* | 45.6% | 45.6% |
| Previously introduced-added regions reclassified native after Castilleja addition | — | 0 |
| Previously Plantago-unique regions additionally covered by Castilleja | — | 0 |

**Outcome:** The two documented omissions affect **interaction-list completeness** but, at the WGSRPD3 regional scale, do **not** affect **geographic envelope completeness** in this single example. The result is compatible with the original 131/245/114 geographic calculation. This is a nontrivial audit outcome (zero difference was not imposed), but **not** evidence that host incompleteness is innocuous across 239 butterflies.

## Region-level host composition despite unchanged envelopes — completed diagnostic

The additional post-hoc source-frozen composition protocol `EUPHYDRYAS_HOST_COMPOSITION_WITHIN_FIXED_ENVELOPE_PROTOCOL_V01.json` and original-data workflow [37772294904](https://github.com/zuizui0223/chocho/actions/runs/37772294904) completed successfully. Its triggering source [commit `343f8278`](https://github.com/zuizui0223/chocho/commit/343f8278f5a957ffcc494ea5ca752b0fcdcecce5) measures exact accepted known-host species counts per WGSRPD3 botanical region, not field larval food suitability.

**Nine newly documented host-species × region incidences affected six already-covered regions; 3 regions gained one and 3 gained two recorded known host species.** The exact native versus contemporary reconstructed host-count changes are:

| Botanical region code | Native 30-host → 32-host count | Contemporary count |
| --- | --- | --- |
| ABT | 2 → 3 | 3 → 4 |
| BRC | 2 → 4 | 4 → 6 |
| IDA | 7 → 8 | 9 → 10 |
| MNT | 5 → 6 | 7 → 8 |
| ORE | 14 → 16 | 16 → 18 |
| WAS | 5 → 7 | 7 → 9 |

- **Area remains unchanged:** native 131 → 131, contemporary 245 → 245, introduced-added 114 → 114.
- **No newly converted single-host regions:** native regions with exactly one *recorded* host stay **39**, contemporary exactly-one host stay **82**; **zero** regions moved from exactly one to two or more recorded hosts.
- **No overlap with baseline introduced-added regions:** all six regions lie in the already-present native resource envelope, and zero of the 114 added regions gain a *Castilleja* incidence. Original *Plantago* unique introduced-range dependence remains 52/114.
- These counts are species-level **network projection** diagnostics within the published, frozen botanical region scheme. WGSRPD3 covers extensive heterogeneous habitats; plants sharing a region do not establish shared individual butterfly habitats, stage-appropriate phenology, nor functional fallback.

The zero regional-area and unchanged single-host results are useful **negative completeness tests**. The new source confirms locally documented host connections omitted from the original HOSTS list, but does not provide evidence that the original region-level analysis missed area expansion, nor that global host-list incompleteness is unimportant. Do not elevate this deterministic post-hoc change in host richness to an independent ecological discovery or a population-resilience test.

## Why this is interesting biologically but not a second result about fitness

Two different types of incomplete knowledge must be separated:

1. **Taxonomic edge missing:** a butterfly's observed ability to feed/oviposit on a plant is omitted from the interaction inventory. The new source study repairs that known relationship locally.
2. **Map area unchanged:** when adding that relationship at its WCVP geographic range, all regions are already covered by other mapped hosts. A union-based regional envelope cannot distinguish the new plant from already-present mapped supply.

The second result **does not validate local host substitution**. In the source experiment, larval growth on *Castilleja* depended on **bracts versus leaves**, and field plant phenology and tissue availability differed from *Plantago* (Haan et al. 2021). Adding a botanical range cannot recover that tissue/season effect. Likewise, Singer & Parmesan (2018), doi:10.1038/s41586-018-0074-6, documented an eco-evolutionary trap mediated by host choice and grazing-induced microclimate change, not a disappearance of botanical WGSRPD3 plant records.

This is thus a particularly explicit case of **regional spatial saturation hiding functional host diversity**. Do not infer a global 'no effect' or correct 239 species by the number of missing pairs in this handpicked example. Neither local experiment measured global host use or the geographic response to reintroducing *Castilleja*.

## Decision

- **Accepted:** Original four-host crosswalk and two-host addition both reproduce, pass their fixed source checks and have audit trails.
- **Rejected:** Treating 32/30 host count change as a change in the geographic opportunity claim; using an unchanged envelope as proof of adequate larval food throughout the butterfly's range; converting a structural *Plantago* 52/114 into demographic dependence; claiming host-list gaps are globally negligible.
- **Stop:** Do not perform further a posteriori host-choice additions or region-cutoff sweeps to force a nonzero area result. Next biological evidence needs stage-specific local host availability, feeding, cohort survival and population responses at independent sites.
- **Paper:** Primary GEB v0.2 resource geography, PR #38 and manuscript figures are unchanged.
