# Full 26 positive added-host butterfly country source census

**23 October 2026.** Completed from the unchanged checksum-verified original
HOSTS/WCVP/BCE/CHELSA source. **No new GBIF plant outcomes were used**.
The five earlier juvenile-labelled cases were not used to filter these 26.

| Input/output | Verified source quantity |
| --- | ---: |
| Butterfly species × country positive record cells | 48 |
| Original resource-added butterfly×country cells | **26** |
| Non-added positive cells retained as source context | 22 |
| Source butterfly×country×accepted host candidate links | **46** |
| Distinct source WCVP plant IDs | **24** |
| Unique plant×country candidates for observational GBIF facets | **44** |
| Unique butterfly taxa in the 26 cells | **16** |
| Distinct countries | **11** |
| Candidate-plant names already passing strict WCVP binomial syntax | **44/44 source query pairs** |

## Fragility of source-recognized plant opportunity

The number of original `WCVP-accepted introduced-source hosts`
per *positive* butterfly × country cell is sharply uneven:

| Number of accepted introduced host candidates | Butterfly×country cells |
| ---: | ---: |
| 1 | 19 |
| 2 | 3 |
| 3 | 1 |
| 6 | 3 |
| **Total** | **26** |

**19/26 cells have just one introduced host plant in the source
catalogue**. This points to low *documented* functional substitution
in the mapped host catalogue, NOT verified exclusive host dependence
or ecological fragility in the field. Source host links may be
taxonomically or biologically incomplete.

All 48 original source resource flags (26 added, 22 nonadded)
reconstruct correctly from the author-only accepted plant IDs,
country/material WGSRPD3 region intersections and native-host
guard. All 20 of the earlier five-case botanical links and the
eight accepted plant IDs have been compared against the original
saved five-case source; no source-identity drift occurred.

## All 26 pre-outcome source cases

| Butterfly | ISO2 country | Introduced host candidate count | Original GBIF focal record count |
| --- | --- | ---: | ---: |
| *Aguna asander* | CU | 1 | 2 |
| *Aguna asander* | DO | 1 | 9 |
| *Aguna asander* | HT | 1 | 3 |
| *Anteos clorinde* | CU | 1 | 34 |
| *Anteos clorinde* | DO | 1 | 15 |
| *Anteos clorinde* | HT | 1 | 1 |
| *Battus philenor* | GT | 1 | 1 |
| *Calycopis cecrops* | VE | 1 | 2 |
| *Colias philodice* | GT | 3 | 4 |
| *Colias philodice* | HN | 2 | 1 |
| *Colobura dirce* | VE | 1 | 56 |
| *Danaus plexippus* | DE | 1 | 2 |
| *Danaus plexippus* | SE | 1 | 1 |
| *Epargyreus clarus* | GB | 2 | 1 |
| *Euploea core* | DE | 1 | 1 |
| *Historis acheronta* | CU | 1 | 23 |
| *Historis acheronta* | VE | 1 | 14 |
| *Kricogonia lyside* | HN | 1 | 3 |
| *Papilio demoleus* | CU | 6 | 35 |
| *Papilio demoleus* | DO | 6 | 199 |
| *Papilio demoleus* | HT | 6 | 8 |
| *Papilio machaon* | NA | 1 | 1 |
| *Phoebis agarithe* | DO | 1 | 24 |
| *Phoebis agarithe* | HT | 1 | 1 |
| *Pieris brassicae* | NA | 2 | 25 |
| *Vanessa atalanta* | GA | 1 | 1 |

The full [source-derived plant name/ID and deduplicated 44-pair
inventory](gbif_all26_source_host_country_inventory_20261023.json)
is intended as the fixed query frame for the next *observation*
stage. A source plant's range in a WGSRPD3 unit is not evidence that
the plant was actually growing near an observed butterfly.
A country with no indexed candidate-host plant records must be
labelled **not recorded**, never biologically absent.

## Verification and execution

- [Source-locked pre-outcome design](../../protocols/GBIF_ALL26_INTRODUCED_HOST_INVENTORY_20261023.md)
- [Completed reconstructive GitHub Actions run](https://github.com/zuizui0223/chocho/actions/runs/39878952882)
- [Source reconstruction script](../../../scripts/derive_gbif_all26_source_host_country_inventory_20261023.py)
- Author-only source zip SHA-256 `408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323`;
  source manifest SHA-256 `06ac6ec3ce506421a6006f994dba7549a0daa74caf3da5318af7a889cb5ff819`.
- The original GEB manuscript and anonymous reviewer archive are
  unchanged and must exclude the independent exploratory source.
