# Multi-source GBIF butterfly time trends cannot be reduced to a newly created provider dataset

**20 October 2026.** This posthoc provenance audit is complete and returns all
**11 dataset-conditioned year facets across the seven locked non-singleton
butterfly species × country cells**. All 11 dataset totals, seven country totals,
and seven year marginals **exactly match** the frozen earlier GBIF source. No
individual occurrence records, observer details or coordinates were obtained.

## Main nontrivial source finding

The striking **86 `Papilio demoleus` records in Israel (2023: 2; 2024:
23; 2025: 61)** are distributed across **three distinct indexed datasets**.
The 2025 total splits as **33 iNaturalist (54.1%), 27 INPA - Observations
(44.3%), and 1 eButterfly (1.6%)**.

| Dataset | Total 2018–2025 | 2023 | 2024 | 2025 |
| --- | ---: | ---: | ---: | ---: |
| iNaturalist Research-grade Observations | 43 | 1 | 9 | 33 |
| INPA - Observations | 42 | 1 | 14 | 27 |
| eButterfly Surveys | 1 | 0 | 0 | 1 |
| **Total** | **86** | **2** | **23** | **61** |

The INPA GBIF dataset was **created in January 2025** but includes
**one indexed observation dated 2023 and 14 dated 2024**, demonstrating why
provider publication timing is not equivalent to observation timing.
Removing all INPA records still leaves an increase from **9 (2024)
to 34 (2025)** among the remaining platforms. Within iNaturalist,
it is **9→33**. Thus **a newly registered INPA dataset alone does not
explain the 2025 marginal**, although source-specific reporting effort,
backfilling, duplicate observers, or increasing population abundance
could all contribute. We did not measure any of those separately.

The original country-scale host catalogue classifies Israel for
`Papilio demoleus` as **without introduced-added host opportunity** under
its specific plant-range and climate definitions. It is not a verified
host-free country, nor evidence of host-independent butterfly expansion.
The sign difference between this case and `Aguna asander` in the Dominican
Republic warns against attributing raw record growth to new hosts.

For `Aguna asander` in the Dominican Republic, all **9** records come
from iNaturalist across **six** observation years: 2019: 2, 2020: 2,
2022: 1, 2023: 1, 2024: 2 and 2025: 1. Repeated years demonstrate
a temporally distributed *recording trail*, not a self-sustaining
population, independent generations, or actual local use of introduced hosts.

## Full source-audited cohort

| Butterfly | Country | Indexed records | Dataset-specific totals | Complete year marginals |
| --- | --- | ---: | --- | --- |
| *Aguna asander* | DO | 9 | iNaturalist Research-grade Observations: 9 | 2019: 2, 2020: 2, 2022: 1, 2023: 1, 2024: 2, 2025: 1 |
| *Aguna asander* | EC | 3 | iNaturalist Research-grade Observations: 3 | 2019: 1, 2020: 1, 2022: 1 |
| *Aguna asander* | HT | 3 | iNaturalist Research-grade Observations: 3 | 2020: 1, 2023: 2 |
| *Danaus plexippus* | DE | 2 | iNaturalist Research-grade Observations: 1; NABU|naturgucker: 1 | 2023: 1, 2025: 1 |
| *Melanargia galathea* | RW | 4 | JRS Wetands data Rwanda, 2019: 4 | 2019: 4 |
| *Papilio demoleus* | IL | 86 | iNaturalist Research-grade Observations: 43; INPA - Observations: 42; eButterfly Surveys: 1 | 2023: 2, 2024: 23, 2025: 61 |
| *Papilio demoleus* | LB | 3 | iNaturalist Research-grade Observations: 2; Observation.org, Nature data from around the World: 1 | 2021: 1, 2024: 2 |

All seven cells are still exactly the seven initially identified in the
frozen temporal panel; **none is an independent replicate of an ecological
host-plant treatment**. The four other countries, their provider combinations
and full year counts are retained above. The parent response categories were
all HUMAN_OBSERVATION. These facet breakdowns cannot validate individual
organisms, taxonomic misidentifications, larvae, egg survival, recruitment
or settlement.

## Audit protocol, validation and limits

The [fixed source-conditioned plan](../../protocols/GBIF_COUNTRY_DATASET_YEAR_CONFOUND_20261020.md)
used exactly the eleven nonzero dataset × species × country strata already
defined by the [seven-cell source receipt](gbif_country_provenance_facet_qc_20261020.json)
and the [six-key dataset registry](gbif_country_dataset_registry_qc_20261020.json).
All occurrence facets retained the same accepted basis filters, geospatial
quality rules, taxonKey, country and observation-year window. Every complete
dataset-conditioned year series was summed back to both the original
dataset totals and original country-level year marginals with **zero drift**.

- **Source evidence:** [full eleven cross-tabs, original GBIF query URLs,
  response hashes and timestamps](gbif_country_dataset_year_matrix_20261020.json).
- **Reproduction:** [audit code](../../../scripts/audit_gbif_country_dataset_year_20261020.py),
  [no-network tests](../../../tests/test_gbif_country_dataset_year_20261020.py),
  [successful GitHub workflow](https://github.com/zuizui0223/chocho/actions/runs/39409318315).
- **Snapshot:** result source SHA-256 for the seven-cell QC
  `69f5f1092fe853ac44080e81247f82715087a6ee13c250a365734bf37e91294e`; registry
  `3d714f4cde971fb6d88e957789e5ec44757d5117ec92b24a7f18e39084c0744c`.

This does **not** confirm biological range expansion, butterfly-native-range
limits, wild breeding, host-mediated recruitment or increased abundance.
The baseline country host-association result (+3.91 pp) and temporal result
(+5.65 pp) **retain their zero-spanning uncertainty** and still fail the
minimum information gate. Source provenance helps distinguish possible
reporting artifacts, but does not remove the need for precise matched
locations, host arrival chronology or population evidence.
