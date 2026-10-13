# v12 — Older occurrence and historical vagrancy falsify simple GBIF first-year colonization claims

2026-10-13 JST. **Post-hoc source-history audit; no claim of new field observations, reproductive success or independent butterfly colonization estimation.** Frozen GEB manuscript PR #38 remains unchanged at `88209b8a215ac1bc9c4c5e67622fa104179c8895`.

## Four additional v11-unreviewed candidates

- **Papilio memnon × AND (Andaman Islands)** — First retained GBIF 2025. Evans' 1931 Census of India Andaman & Nicobar butterfly list reports *P. memnon agenor* as a **straggler** on the Andamans and Car Nicobar. The fact of historical occurrence **does not demonstrate a resident breeding population**, much less a new colonization in 2025. Source: https://censusindia.gov.in/nada/index.php/catalog/32394/download/35575/40153_1931_REP.pdf . The distinct *Papilio mayo* is documented separately.
- **Pyrgus communis × MSO (Missouri)** — First retained GBIF 2019. BAMONA's January 6, 2010 [Boone County retrospective records](https://www.butterfliesandmoths.org/blogs/829/new-boone-county-missouri-records) includes *P. communis* at Big Muddy NWR. This proves an earlier regional report, not larval development or establishment.
- **Pyrgus communis × MIN (Minnesota)** — First retained GBIF 2024. Minnesota DNR's 2014 conservation species profile describes *P. communis* as a species occurring in Minnesota. This is a dated **agency state-level species account**, not a voucher- or date-specific record (https://www.lrl.mn.gov/docs/2015/other/150681/PFEISref_2/MDNR%202014b.pdf).
- **Pyrgus communis (now Burnsius communis) × WDC (District of Columbia)** — First retained GBIF 2025. [BAMONA record 1144419](https://www.butterfliesandmoths.org/explore-data?region=45223) has **observation date October 1, 2017**, but **verification date November 30, 2021**. An event date before the train/test split does not mean a source was published/verified before 2018; any true prospective hindcast must respect when this information actually became available.

## Cumulative v12 source-audited candidate panel

| Outcome from same 26 introduced-only, late-first source-stable species×regions | n cells | n retained GBIF occurrence rows |
|--|--:|--:|
| Earlier regional presence, vagrant account or early-dated verified observation | **12** | **188** |
| Taxonomy identification hold (*Pyrgus communis* × TEN) | 1 | 5 |
| Still not host/history audited | **13** | **22** |
| Total | **26** | **215** |

Eleven of the 12 were published by 2017; one (*Pyrgus* × WDC) was observed in 2017 but verified in 2021. All sources were retrieved post-hoc and may overlap the original GBIF records. None of these counts forms a fair estimate of false-positive colonization claims in the global population.

Historical vagrants and established populations are distinct ecological categories. No new confirmed breeding episode was found during v12; larval and organism-stage records from earlier v11 remain separately reported.

## Biological next step, not an unearned claim

Historic *Troides aeacus* × WHM work (Smetacek 2011; DOI https://doi.org/10.18473/lepi.v65i2.a5) describes the larval host *Aristolochia dilatata* and mature seed-bearing plants, providing a plausible test of **host developmental stage and quality** rather than mere regional host presence. Observational rearing or host-stage preference in one species cannot be extrapolated to 239 butterflies without host availability/phenology and matched field survival.

Avoid repeating an undifferentiated post-2017 presence/absence hindcast using the 28/31 species 1,800-record capped sample: it is already known to conflate recording, source update date, transient sightings, identification, and genuine breeding.

## Checks and files

The standalone v12 package contains the full unchanged v11 archive, its original 26-case table, four case additions with source URLs/date type, a deterministic merge script, independent no-shared-code CSV audit, four regression tests, all outputs and SHA256 manifest. Fresh v11: 5 tests passed; v12: 4 tests passed, independent counts rechecked. Final v12 ZIP SHA256: `0152c5a5d7f9d51ec3d64b804b3bf815ea38d288c45170e92a4382246e62870f`.

Receipt: `provenance/reviewer_defenses/results/butterfly_temporal_historical_addendum_v12_20261013.json`. No frozen paper, author-facing declaration or blinded bundle was changed, and cumulative code was not merged into the remote analysis branch.
