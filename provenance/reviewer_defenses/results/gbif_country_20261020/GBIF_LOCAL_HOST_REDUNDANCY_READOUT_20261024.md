# Plant resource portfolios lose observed redundancy at subcountry grain

**24 October 2026 — successful explicitly posthoc source-index macroecological audit.**
This is an extension of the **26 pre-existing positive butterfly×country cells**,
not an independent sample of 26 butterfly populations or a measured
butterfly recruitment rate. It uses only checksum-audited original
GADM1 GBIF butterfly and accepted candidate plant species recording
facets, no new occurrences, no raw localities, and no re-estimated
country-level source-host opportunity.

## Exact local recording result

Of the 26 positive source-added-host butterfly×country cells,
**19 have only one recognized introduced source host and seven have
multiple source candidate plants** (three of the latter have six).
The source has 24 distinct accepted plant taxa, 46 focal
butterfly×country×plant links, 44 unique GBIF plant×country
queries, and 3,595 indexed plant records (3,203 of those from
*Robinia pseudoacacia* / Great Britain).

Counting **distinct source candidate plant taxa recorded by GBIF in
the SAME first-level administrative division as a recorded butterfly**,
the full 453 mapped butterfly GBIF record indices are distributed as:

| Number of candidate plant species with an indexed record in that administrative division | Butterfly GBIF records | Butterfly-positive administrative divisions |
| --- | ---: | ---: |
| **0** | **222** | **66** |
| **1** | **226** | **38** |
| **2 or more** | **5** | **2** |
| **Total** | **453** | **106** |

Thus **226/231 (97.84%)** of butterfly GBIF record indices sharing a
GADM1 unit with **any** documented candidate host plant share their
GADM1 with **only one recorded candidate plant species**. In the
seven cases where the plant catalogue actually provides two or more
candidate introduced host species, there are **129** GADM1
co-detected butterfly records, and only **5/129 (3.88%)** have
two or more candidate host species co-recorded.

**This is an observation-level pattern, not 97.84% biological
dependence on a single actual larval host.** Four of the 26 source
cases had no candidate plant GBIF country record at all, and low
plant-recording depth reduces the probability of demonstrating
within-unit redundancy.

## Preserved source-count sensitivity

The focal GBIF butterfly record cohorts (453 indices across 106
butterfly-positive GADM1 units) are *held exactly fixed*. Only
the minimum number of indexed GBIF records **of each distinct
candidate plant species within the same GADM1 unit** changes.

| Botanical minimum per plant species and division | Butterfly indices without qualifying plant | With exactly one qualifying candidate plant species | With >=2 qualifying candidate species | With >=1 plant / 453 | Multi-candidate-case butterfly records with >=2 plant species |
| --- | ---: | ---: | ---: | ---: | ---: |
| ≥1 | 222 | 226 | 5 | 231/453 | 5 |
| ≥2 | 351 | 101 | 1 | 102/453 | 1 |
| ≥3 | 408 | 44 | 1 | 45/453 | 1 |
| ≥5 | 425 | 28 | 0 | 28/453 | 0 |

When requiring ≥5 plant records per plant and division, **no**
butterfly-indexed GADM1 has two or more separately documented
candidate host taxa under these source filters; only **28/453**
butterfly indices meet even one qualifying plant record. Tightening
this observation threshold must reduce the documented overlap
by construction. The sensitivity measures **data support**, not true
availability, botanical abundance or an ecological interaction
threshold.

## Six documented plant candidate names do not guarantee local overlap

In the three *Papilio demoleus* cases (Cuba, Dominican Republic,
Haiti), the original HOSTS/WCVP source lists the same six
introduced-range botanical candidate species. Yet:

- **Dominican Republic:** 107 butterfly records fall in
  candidate-plant-positive GADM1 units, **zero** have two
  candidate plant species indexed in the *same* division.
  Removing only *Triphasia trifolia* from the index would
  reduce the 107 co-detected butterfly indices by **54**
  (~50.5%); this is a hypothetical *recording coverage
  change*, not a plant extinction scenario.
- **Cuba:** 9 butterfly indices co-detected with one candidate
  plant species in their administrative unit; none with two.
- **Haiti:** 8 butterfly indices co-detected with a candidate
  plant; only 4 share the division with two indexed candidate
  host species. No other 26-case source country except the
  *Epargyreus clarus* / Great Britain single-record
  source case contains ≥2 candidate plant species at the
  recorded butterfly administrative division level.

These are **source-host network and GBIF observational geography**
patterns. Plant GBIF observations might be cultivated, from
different habitats or years, incomplete or unrelated to butterfly
feeding. Original accepted host links may contain outdated
biological assignments. The same GADM1 is neither same-site
co-occurrence nor verified oviposition/feeding.

## Immature-stage source cross-check

For the separately fixed **five** juvenile-labelled,
added-host-country cases, all **24 immature-stage GBIF source
labels** occur in GADM1 units with **0**, **1**, or **≥2**
indexed candidate plant species in counts **8 / 14 / 2**.
The original 39 raw juvenile-stage records (three species, 11
country pairs) are not 39 independently validated feeding or
recruitment events.

A separate audit reviewed only the publicly available selected
GBIF record annotation fields for all **39** juvenile-labelled
records. **12 records had any nonempty selected text field;
none had an exact source-accepted candidate host plant binomial
mention or an explicit host-use/feeding keyword under the fixed
parser.** Missing structured host-use evidence is **not** a
test of whether larvae actually used host plants. See the
[privacy-safe missingness summary](gbif_juvenile_host_annotation_coverage_20261023.json).

## Ecological inference gate

One may hypothesize that regional host-resource portfolios
are spatially complementary instead of locally redundant,
or that recorded butterflies disproportionately encounter
few source-listed candidates. **Neither can be established as
a field ecology conclusion from these data**, because the
botanical index is sparse and the butterflies were selected
on recorded presence and resource opportunity. There is no
species×GADM matched biological *absence* background,
direct plant association, wild establishment or larval
development outcome.

More country-only reanalysis cannot fill this gap. An empirical
host-realization study would need the original observation
identities with photo/taxon adjudication, flowering/leafing
and cultivation context, validated host plant contact,
standardized survey effort, and linked juvenile recruitment.

## Source verification and separation

- [Exact data-only numerical summary](gbif_local_host_redundancy_summary_20261024.json);
  [fully reproducible code](../../../scripts/analyze_gbif_all26_local_host_redundancy_20261024.py).
- [Source-only workflow run 40073886529](https://github.com/zuizui0223/chocho/actions/runs/40073886529)
  (offline source replay with five regression tests).
- [Original 44 GBIF botanical-country facet audit](gbif_all26_host_gadm_source_20261023.json);
  [48 butterfly country/GADM original source](gbif_positive48_gadm_lifestage_facets_20261020.json);
  [five stage-label original cohort](gbif_source_introduced_hosts_five_cases_20261021.json).
- [Host-use annotation availability source run 39919887018](https://github.com/zuizui0223/chocho/actions/runs/39919887018),
  persist **aggregated numeric privacy-safe output only**.
- Fixed original GEB v0.2 manuscript and its anonymous submission
  package remain unchanged.
