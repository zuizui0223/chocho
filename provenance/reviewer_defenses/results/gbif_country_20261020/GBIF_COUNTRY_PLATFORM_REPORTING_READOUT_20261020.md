# The 2025 Israel Papilio record increase persists after matching each data publisher's reporting volume

**20 October 2026 — complete posthoc reporting-effort diagnostic**, using
the **same previously selected seven butterfly country cells** and eleven
dataset×country cohorts, not a new ecological treatment or independent
global replication.

All **22** compact GBIF year-facet queries succeeded; 11/11 all-taxa and
11/11 legacy-key Lepidoptera denominators have exact valid year totals and
contain all focal counts for the same dataset/country/year. The
legacy-order-key nesting check is necessary but **not sufficient** to prove
taxonomy remains perfectly comparable across 2026 indexing. No individual
records, observers, coordinates, or breeding evidence were requested.

## Matched 2024–2025 denominator results for Papilio demoleus / Israel

| Publishing source | Denominator | Focal records 2024→2025 | Denominator records 2024→2025 | Focal per 10,000 records 2024→2025 | Proportion change |
| --- | --- | --- | --- | --- | --- |
| iNaturalist | 同一提供元・全分類群 | 9→33 | 34,270→46,687 | 2.63→7.07 | **2.69×** |
| iNaturalist | 同一提供元・旧GBIFチョウ目キー | 9→33 | 3,056→4,436 | 29.45→74.39 | **2.53×** |
| INPA | 同一提供元・全分類群 | 14→27 | 143,482→161,990 | 0.98→1.67 | **1.71×** |
| INPA | 同一提供元・旧GBIFチョウ目キー | 14→27 | 6,777→1,795 | 20.66→150.42 | **7.28×** |

The original **86** Israeli `Papilio demoleus` observations comprise 43
iNaturalist, 42 INPA and 1 eButterfly records over 2023–2025.
The **61 records dated 2025** split 33 iNaturalist, 27 INPA,
1 eButterfly; the latter has no 2024 focal case and cannot support a
year-over-year within-source ratio.

**Key interpretation:**
- For iNaturalist, *all-taxa reporting volume rose 36.2%* (34,270→46,687)
  while focal records rose from 9 to 33. Thus the focal share rose
  **2.69-fold**. The more target-aligned Lepidoptera denominator rose
  45.2% (3,056→4,436), giving a **2.53-fold** share increase.
- For INPA, all-taxa volume rose **12.9%** (143,482→161,990) while
  focal records rose 14→27, giving **1.71-fold** share increase.
  But the legacy-key Lepidoptera denominator *fell sharply* from
  6,777 to 1,795 between 2024 and 2025, producing a **7.28-fold**
  nominal focal share increase. This unstable taxon-specific denominator
  could reflect reporting mix, lag, taxonomic processing or survey
  priorities and should **not** be interpreted as stronger population growth.
- The source-conditional 2025 increase is not merely the one-time
  introduction of INPA to the GBIF dataset registry in January 2025:
  both iNaturalist and INPA contain records already dated 2023 and 2024.
  Nevertheless, neither broad nor Lepidoptera-normalized records
  equal standardized butterfly abundance, larval host use or viable
  wild populations.

## Dominican Republic Aguna asander

The previously source-audited `Aguna asander` / Dominican Republic case
has **9 iNaturalist observations spread over six years**, 2019, 2020,
2022–2025. These were 2,2,1,1,2,1 focal records while the iNaturalist
all-taxa country reporting activity grew markedly over the same period
(4,603 in 2019; 13,935 in 2025). Its Lepidoptera records increased
from 879 to 1,395 over these two endpoints. Such very low focal counts
cannot support a demographic trend.

## All seven cases and source-bounded population

The [full saved 22-query result](gbif_country_platform_reporting_20261020.json)
contains every pair's year-by-year broad all-taxon counts, legacy
Lepidoptera counts, focal counts, nesting statuses, original dataset keys,
GBIF query URLs, response SHA-256, and UTC fetch receipt. Two examples
were reported prominently **because the prior fixed source selected
the >=5-record subset**, not because the outcome favored a host effect.

The initial country-scale host-resource model labels **Israel/Papilio as
not having introduced-added-host opportunity** in that narrowly defined
catalogue, whereas **Dominican Republic/Aguna has added-host opportunity**.
Neither country-scale flag proves absence of host plants at local encounter
sites. A within-platform observation-share increase in one *non-added*
country undermines simplistic attribution of all butterfly records to
documented introduced hosts, but does not identify an alternative
biological cause.

## Reproducibility and boundary

- [Fixed posthoc protocol](../../protocols/GBIF_COUNTRY_PLATFORM_EFFORT_DENOM_20261020.md)
- [22-query source and verification receipt](gbif_country_platform_reporting_20261020.json)
- [Re-executable Python code](../../../scripts/audit_gbif_country_platform_reporting_20261020.py)
- [Completed GitHub workflow](https://github.com/zuizui0223/chocho/actions/runs/39410608472)
- Data source SHA-256: `ceb7571e8e2d0fb346e60d0c8bec7b46a890d870d2fa8d5a3f9af460933cb2f5`
- 2018–2025 only; each provider's records are indexed occurrences, not
  independent sampled insects. Different platforms may share individual
  sightings. The 2026 GBIF taxonomy transition limits historical comparability
  even though every legacy Lepidoptera denominator passed the nested-count
  check.
- The full 229-taxon country comparison remains **+3.91 pp with 95%
  genus-bootstrap interval −6.42 to +14.94** and only **21 informative
  strata**. The paired-year comparison remains **+5.65 pp, −15.02 to +22.04**
  with **10 informative strata**. Neither passes the original gate
  of 50 informative strata, and neither establishes host-driven range
  expansion, population persistence, or recruitment.

**Further evidence required:** occurrence-level locality clusters, deduplication,
standardized recording effort and direct host use/recruitment or dated
arrival of host plants; the aggregate facet source cannot provide these.
