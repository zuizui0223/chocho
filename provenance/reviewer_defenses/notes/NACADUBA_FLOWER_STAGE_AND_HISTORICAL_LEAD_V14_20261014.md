# v14: Historical plant lead and floral-stage larval opportunity in Nacaduba

Date: 2026-10-14. **Exploratory retrospective natural-history audit, not a study of causal colonization or estimated butterfly invasion lag.** The GEB submission PR #38 and its frozen numerical estimates are unchanged.

## Source chronology
Nghiem, Tan & Corlett (2015, [Table 2](https://doi.org/10.1177/194008291500800116)) list **1890** as the first official record of the invasive Australian *Acacia auriculiformis* in Singapore. Local *Nacaduba biocellata* was retrospectively photographed in **2004**, added to the local species checklist in **2008**, and documented as feeding on alien *A. auriculiformis* in immature stages and then eclosing under rearing conditions in **2010**.

The arithmetic record separations are **114, 118 and 120 years**. These are *documentation intervals*, not plant introduction–to–butterfly arrival, evolution, colonization, population establishment, or insect response delays. The 2004 photograph was only publicly identified retrospectively, and neither the exact introduction/establishment date nor plant flowering status at each butterfly observation is known.

Sources: [2008 checklist](https://butterflycircle.blogspot.com/2008/02/new-lycaenidae-species-for-singapore_10.html), [2014 review of first 2004 photo and seasonal pulses](https://butterflycircle.blogspot.com/2014/01/butterflies-galore-two-spotted-line-blue.html), [2010 detailed life history](https://butterflycircle.blogspot.com/2010/01/life-history-of-two-spotted-line-blue.html). A 2015 Singapore flora survey establishes that *A. auriculiformis* is a common invasive tree, not necessarily a direct conservation benefit.

## Life-stage constraint based on a genuine field-host sample
2010 Singapore observers collected **20+ immature stages** from flowering *A. auriculiformis* in two visits within one week and described flower-bud and flower feeding, eggs near buds, and captive eclosion. Reported egg, four instar, prepupal and pupal durations yield:

- **9.5–11 days** of larval bud/flower feeding;
- **11.5–13.5 days** from egg laying to the end of larval feeding;
- **17.5–19.5 days** from egg to adult **under the described rearing conditions**.

These are descriptive sums of approximate timings, **NOT confidence intervals, field survival rates, the life of a single individual blossom, or a measured climate or flowering mismatch**. Resources can be supplied by sequential buds or plants rather than one lasting bloom.

## Triangulation with natural-range larvae and other introduced hosts
Braby (2015, [Western Australia field-voucher study](https://doi.org/10.18195/issn.0312-3162.30(2).2015.073-097)) collected **two *N. biocellata biocellata* larvae from flowering *Acacia plectocarpa*** on May 16, 2011. Reared adults emerged on May 25 and May 29. These 9/13 days are **from larva collection to adult**, not egg-to-adult. The Singapore wild larvae were confirmed on alien *A. auriculiformis*; another alien, *A. mangium*, was accepted **in captivity only**. Thailand *N. b. baliensis* was reported on alien *Pithecellobium dulce* in Jennarongsak (2026, *Butterflies* 98:52–56), with an identifiable bred specimen dated Jan 22, 2025 ([original specimen compilation](https://yutaka.it-n.jp/lyc4/81285001.html)). The Thai larval feeding organ is **not identified** from the reviewed open metadata.

All these named *Acacia*, *Pithecellobium* and frozen HOSTS *Vachellia* belong to Fabaceae, specifically the [mimosoid clade nested within Caesalpinioideae](https://doi.org/10.12705/661.3). Thus the observed association is **geographically novel local access to related plants** rather than evidence of evolutionary leaps to unrelated host lineages. Thai and Singapore butterfly subspecies identity was not independently proven to be identical, and no ancestral state was reconstructed.

## Reconstructed geographic sensitivity (from pinned v13 input)
The one-butterfly HOSTS+BCE host mapping yields **50 native / 131 current / 81 added** butterfly×WGSRPD3 cells. Hypothetically extending the two observed local wild host plants globally for this butterfly yields **54 / 157 / 103** (+4 native, +26 contemporary, +22 introduced-only). THA and MLY are classed introduced-only both before and after. **This is not 22 real colonizations.** It is a test of host-list completeness and demonstrates that knowing a host-bearing region exists is distinct from correctly identifying which host plants can support breeding there.

## Testable ecological prediction
For a focal butterfly lineage with a roughly ten-day bud/flower-feeding phase, local reproductive opportunity will depend on the **density and continuity of suitable flowering tissue**, *in addition to botanical occupancy* and host lineage. The Singapore butterfly has repeated seasonal population bursts in natural-history reports, while NParks labels the host evergreen and [free-flowering](https://www.nparks.gov.sg/florafaunaweb/flora/2/6/2693). **This apparent contrast is not a tested relationship:** budding intensity through time was not tracked; larval natural enemies, ant attendance, weather, dispersal, mowing/planting and detection may all explain the pattern.

More spatial permutations or more post-hoc GBIF first-year checks cannot identify this mechanism. A distinguishing design would record sequential flowering-tissue density, egg dates, wild larval cohorts, pupation and observed eclosion across multiple host trees and regions, plus a comparison against periods when host trees are present but lack buds.

## Verification and scope
A reproducible source-staged local v14 full verification package exists with pinned v13 numeric sources, verified citations, six regression tests, an independent Decimal/CSV audit and twelve passed parent v13 tests. The self-contained GitHub reproduction under `research/nacaduba_stage_v14/` stores precisely labelled source facts, a standard-library timeline/stage calculator and tests. GitHub's smaller reproduction checks the **source-ledger arithmetic**; the more complete v13 source archive in the local bundle is needed to audit the original spatial data reconstruction. The publication CI and anonymous review bundle were not rerun or modified.

Local archive SHA256: `f5c3bddc9fdea5bce81cf9e019a674ed5ee8af5119f127dd80f17b1e24a2f70a`.
