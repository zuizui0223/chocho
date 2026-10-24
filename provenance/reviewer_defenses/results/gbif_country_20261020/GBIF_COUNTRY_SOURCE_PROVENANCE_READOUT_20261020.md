# Source-provenance audit of the seven non-singleton butterfly country detections

**Completed 20 October 2026** from full-index compact GBIF occurrence facets
[reproducible audit run](https://github.com/zuizui0223/chocho/actions/runs/39401766071). All 7 GBIF live-index
counts **exactly matched** the frozen country facets, and all 7 response
triplets (datasetKey, year, basisOfRecord) were complete, unique and internally
consistent. Data retrieval was prospectively bounded to the cohort selected
by the [source-blind audit protocol](../../protocols/GBIF_COUNTRY_SOURCE_INDEPENDENCE_GATE_20261020.md).

| Focal butterfly | Country | Later records (2018–2025) | Years indexed | GBIF datasets | Largest dataset share | Added-host flag in original country model |
| --- | --- | ---: | --- | ---: | ---: | --- |
| *Aguna asander* | DO | 9 | 6 (2019, 2020, 2022, 2023, 2024, 2025) | 1 | 100% | yes |
| *Aguna asander* | EC | 3 | 3 (2019, 2020, 2022) | 1 | 100% | no |
| *Aguna asander* | HT | 3 | 2 (2020, 2023) | 1 | 100% | yes |
| *Danaus plexippus* | DE | 2 | 2 (2023, 2025) | 2 | 50% | yes |
| *Melanargia galathea* | RW | 4 | 1 (2019) | 1 | 100% | no |
| *Papilio demoleus* | IL | 86 | 3 (2023, 2024, 2025) | 3 | 50% | no |
| *Papilio demoleus* | LB | 3 | 2 (2021, 2024) | 2 | 67% | no |

**All seven positive country cells have `HUMAN_OBSERVATION` only**, with no
preserved specimens in these very specific filtered source queries. This does
not mean no specimens exist anywhere; only that none occur within the
seven source-locked filtered country/time combinations.

## Two >=5-record comparisons

1. **Aguna asander / Dominican Republic:** exactly 9 indexed records in
   2019, 2020, 2022, 2023, 2024 and 2025 (2, 2, 1, 1, 2, 1 by year).
   All nine are in one dataset, GBIF key
   `50c9509d-22c7-4a22-a47d-8c48425ef4a7`, the
   [iNaturalist Research-grade Observations](https://www.gbif.org/dataset/50c9509d-22c7-4a22-a47d-8c48425ef4a7)
   source. The original country host catalogue flags introduced-added larval
   resource opportunity. Multiple calendar years *suggest repeatedly indexed
   observations*, but do not establish independent individuals, larval
   recruitment, breeding habitat, or a host-mediated range expansion.

2. **Papilio demoleus / Israel:** exactly 86 indexed records spanning
   2023 (2), 2024 (23) and 2025 (**61**). Three datasets contribute
   **43, 42 and 1** records, with the largest fraction 50.0%.
   The 43-record component is iNaturalist; the other two are catalogued
   in the [source facet receipt](gbif_country_provenance_facet_qc_20261020.json)
   but not yet independently validated for duplicate observations, publishing
   provenance or wild establishment. The original country catalogue **does
   not** flag an added-host-resource country. This is a useful recording
   counterexample to interpreting higher record counts as host-enabled
   recruitment.

## Remaining five cells and stronger independence requirements

*Aguna asander* / Ecuador and Haiti have 3 records each, both entirely
from iNaturalist, spanning 3 and 2 years respectively. *Danaus plexippus*
/ Germany has two records in 2023 and 2025 across two dataset keys.
*Melanargia galathea* / Rwanda has four records in **2019 only**, from
one dataset. *Papilio demoleus* / Lebanon has three records (2021 and 2024)
from two datasets.

Four of seven cells each draw **100% from a single GBIF dataset**, and six
of seven rely in part on iNaturalist. Dataset count is a provenance
heterogeneity metric, **not a count of independent field surveys**.
Record count and observation year are subject to sampling effort, repeated
visits, occasional vagrants and potentially duplicate indexing.

None of the seven has explicit larval host-use, local recruitment,
establishment, or dated introduction evidence in these three compact facets.
The macro effect remains unresolved. Do not promote these seven observations
or the 86-record maximum to a global conservation conclusion.

## Dataset registry identities validated

The six unique dataset keys in these seven late-only candidate cells were
subsequently resolved directly against the GBIF dataset registry (all 6/6
successful). This is **posthoc provenance metadata** only, not an
independence or successful breeding verification.

| GBIF dataset registry title | DOI | Key (prefix) |
| --- | --- | --- |
| iNaturalist Research-grade Observations | 10.15468/ab3s5x | `50c9509d` |
| NABU|naturgucker | 10.15468/uc1apo | `6ac3f774` |
| INPA - Observations | 10.15468/nqx5dj | `6c8d87be` |
| Observation.org, Nature data from around the World | 10.15468/5nilie | `8a863029` |
| JRS Wetands data Rwanda, 2019 | 10.15468/k5weqa | `8fbdc3e6` |
| eButterfly Surveys | 10.15468/ykxm8x | `cf3bdc30` |

In the most repeatedly indexed *Papilio demoleus* / Israel comparison, the
86 observations split as **43 iNaturalist Research-grade Observations,
42 INPA - Observations and 1 eButterfly Surveys**, representing three
different GBIF publishing-organization keys. The three platforms provide
source heterogeneity, but GBIF facet totals cannot establish that their
underlying observations are deduplicated, independent visits or breeding
records. The 61 records indexed in 2025 are only a marginal year total;
this audit did not measure which of those years came from which platform.

The two German *Danaus plexippus* observations split across iNaturalist and
**NABU|naturgucker**, while all four 2019 Rwandan *Melanargia galathea*
records belong to the single sampling-event resource
**JRS Wetands data Rwanda, 2019**. These provenance differences show why
record counts should not be substituted for validated demographic evidence.

[Complete registry titles, DOI, publishing keys and source hashes](gbif_country_dataset_registry_qc_20261020.json)
are saved independently of the ecological effect calculation; the
[successful metadata workflow](https://github.com/zuizui0223/chocho/actions/runs/39402802072)
and exact six-ID source receipt are retained. The frozen country and GEB
manuscript results remain unchanged.

## Dataset × year cross-tab validation (posthoc)

All 11 dataset-conditioned year facets for the seven fixed cells were
fetched and reconciled with the original marginal totals, without any drift.
In Israel the `Papilio demoleus` **61 records dated 2025** split as
**33 iNaturalist, 27 INPA, 1 eButterfly**; excluding all INPA observations
still leaves **9 indexed records in 2024 and 34 in 2025** from other
providers. The INPA dataset was registered in 2025 but contributes
observations dated 2023 and 2024 as well. This establishes that a single
newly published provider dataset does not by itself explain the later-year
distribution; it does **not** establish independent abundance growth.
The Dominican Republic `Aguna asander` 9 records remain distributed across
six years but belong to a single dataset. See the
[complete dataset-year readout](GBIF_COUNTRY_DATASET_YEAR_READOUT_20261020.md)
and [source-validated JSON](gbif_country_dataset_year_matrix_20261020.json).

## Reproducibility and verification

- Original country aggregate (unchanged): 23,433 rows, 229 noncolliding taxa,
  33,817,891 indexed records within accepted basis types.
- Prospective selection: **zero** 2010–2017 records, **>=2** 2018–2025
  records, other recorded butterfly richness >=3 in **each** 8-year period.
  Seven positive cells (three added-host-resource, four control).
- Exactly seven live facet queries; all returned HTTP 200, complete facet
  counts, three requested facet categories and unchanged frozen record
  totals. UTC retrieval ended at `2026-10-20T09:35:38.177233+00:00`.
- [Saved complete facet summary, query URLs and response SHA-256](gbif_country_provenance_facet_qc_20261020.json).
- [Reproduce via checked source-aware script](../../../scripts/audit_gbif_country_record_provenance_20261020.py)
  and [workflow 39401766071](https://github.com/zuizui0223/chocho/actions/runs/39401766071).
- This audit does **not** change the original statistical effect, source
  benchmark, or GEB manuscript. No individual observation records were
  downloaded or redistributed.
