# Evidence-filtered European butterfly foodplants: reproducible source pilot

**Date:** 2026-10-08. **Status:** PILOT POSITIVE SOURCE ACCESS / TAXONOMY UNRECONCILED / NO CAUSAL ECOLOGICAL EFFECT. Analysis branch only; GEB manuscript and PR #38 unchanged.

## Frozen protocol and source identity

After the original Clarke (2024) Dryad / publisher / Europe PMC source downloads failed in runs [37774669703](https://github.com/zuizui0223/chocho/actions/runs/37774669703) and [37775172709](https://github.com/zuizui0223/chocho/actions/runs/37775172709), we froze the distinct `BCE_CLARKE_HIGH_EVIDENCE_HOST_PILOT_PROTOCOL_V01.json` before inspecting the 12 chosen butterfly–host comparisons.

The Butterfly Conservation Europe (BCE) [taxonomy](https://www.bc-europe.eu/taxonomy.php) and [butterfly species pages](https://www.bc-europe.eu/butterfly.php?genus=Iphiclides&species=podalirius) explicitly state their foodplant names derive from **Clarke (2024)**, DOI [10.1002/ece3.10834](https://doi.org/10.1002/ece3.10834), and **exclude evidence statuses 4–6**. It is an **independently compiled secondary/derived evidence list**, not original raw Dryad, field observations, country-specific host preference, or independently generated studies. Original comparison is pinned `globalbioticinteractions/HOSTS` commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`.

The fixed selection used the full 239-species chocho panel and all 501 BCE butterfly taxa, intersected **exact scientific binomials**, then selected the first 12 alphabetically **before comparing any host names**. The 92-species exact butterfly taxonomic intersection was an outcome of the source-access run, not a retroactively edited panel.

## Successful exact-name source audit

[GitHub Actions run 37776666403](https://github.com/zuizui0223/chocho/actions/runs/37776666403) completed SUCCESS; source hashes, page URLs, excluded generic/infraspecific plant names, and full per-butterfly name lists are in artifact `chocho-bce-clarke-evidence-host-pilot-v01`.

- 239 frozen chocho butterflies; 501 BCE European species rows; **92 exact binomial intersections**.
- Twelve alphabetically chosen European butterfly species were attempted, **11** with at least one exact species-rank BCE foodplant.
- In the 11 eligible species: **213** distinct BCE exact binomial butterfly–plant associations; **81** HOSTS exact binomial associations; **42** overlap in exact plant spelling; **171 BCE-only exact spellings** (213 minus 42).
- The twelfth, *Boloria napaea*, had no species-rank foodplants recoverable from the site but might have genus/infraspecific data; its zero is **not hostless-in-nature evidence**, and it is excluded from these 11-species totals.
- The above is an intentional small **alphabetic/European** sample, **not** representative of all 239 butterflies or all of Clarke's 464 species with foodplant information.

Per-butterfly source record counts:

| Butterfly | BCE exact plant species | HOSTS exact raw plant species | Matching exact binomials |
| --- | ---: | ---: | ---: |
| *Aglais urticae* | 2 | 6 | 2 |
| *Anthocharis cardamines* | 61 | 11 | 7 |
| *Apatura ilia* | 8 | 10 | 6 |
| *Aphantopus hyperantus* | 22 | 7 | 6 |
| *Aporia crataegi* | 27 | 14 | 6 |
| *Araschnia levana* | 2 | 4 | 2 |
| *Argynnis paphia* | 7 | 11 | 3 |
| *Boloria pales* | 5 | 3 | 1 |
| *Brenthis daphne* | 6 | 4 | 2 |
| *Callophrys rubi* | 59 | 10 | 6 |
| *Carcharodus alceae* | 14 | 1 | 1 |

The overlap 42/213 is **not sensitivity for biologically confirmed feeding**: HOSTS matches have not yet been synonym-reconciled to BCE/Clarke taxonomic concepts, BCE includes only a partial geographical scope, and their primary-literature sources can overlap. Large BCE-only counts are a cue for quality-controlled taxon joins and original reference adjudication, not a robust global false-negative estimate.

## Next frozen scientific decision

Expansion from 12 to all **92** already defined exact-butterfly BCE–chocho overlaps would remove the small alphabetic-subset bias, but still remain a geographically bounded European source. Its design must be explicitly frozen after inspecting the pilot result, and includes every matched butterfly whether or not BCE has species-level foodplants. Compute exact names only first; **stop** before calling BCE-only names biological missing edges. Follow with exact taxon acceptance/synonym audit and source references before any map changes. No new results about butterfly colonization, interaction strength, local demography, or fitness.

**Primary GEB paper remains unchanged.**
