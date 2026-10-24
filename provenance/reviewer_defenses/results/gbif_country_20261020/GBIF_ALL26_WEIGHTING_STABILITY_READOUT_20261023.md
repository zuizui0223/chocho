# When host co-recorded geography changes meaning under alternative macroecological weights

**23 October 2026 — exact source-replay and explicitly POSTHOC sensitivity.** An
earlier [source-verified exact all26 host-identity reference](https://github.com/zuizui0223/chocho/actions/runs/39880996039)
found **231 of 453 butterfly records** in GADM1 units containing a
recorded candidate introduced host plant, versus a conditional mean
**230.0833/453** under country-wise rearrangement of accepted host
identity with both host-link degree sequences fixed. The original
all-record difference **+0.2024 pp** is not compelling on this
artificial, nonbiological reference.

**New finding:** the same source geography and exact 13,826,598,912
possible country-degree graphs produce markedly different
*descriptive departures* when choosing a different ecological sampling
unit. All following reference tails are **not ecological p-values**.

| Pre-existing data summarized as | Unit denominator | Observed GADM co-recording | Artificial reference expectation | Difference | Exact source-label upper tail, NOT biological p |
| --- | ---: | ---: | ---: | ---: | ---: |
| Record-weighted 453 GBIF indices | 453 | 50.99% | 50.79% | +0.20 pp | 0.6196 |
| 26 butterfly×country equal-weight | 26 | 49.41% | 42.55% | +6.86 pp | 0.02771 |
| 16 butterfly taxa equal-weight | 16 | 49.47% | 44.29% | +5.18 pp | 0.09932 |
| 11 country equal-weight | 11 | 52.96% | 48.47% | +4.49 pp | 0.04555 |
| 17 cases with ≥2 records | 17 | 46.16% | 38.78% | +7.38 pp | 0.04742 |
| 11 cases with ≥5 records | 11 | 48.61% | 38.33% | +10.28 pp | 0.0002137 |
| 9 cases with ≥10 records | 9 | 47.07% | 36.00% | +11.06 pp | 0.001058 |

## The nontrivial scale problem

Treating **a randomly chosen indexed butterfly record** as the target
produces almost no additional alignment to the documented candidate host
identity (+0.20 pp). Treating **a randomly chosen butterfly×country
source case** as the target produces +6.86 pp. Both are technically
correct for their **different** estimands on the **same records**.
Choosing the estimand after looking at those results can create
false assurance that a biological mechanism has been discovered.

Neither species nor country weighting is stable justification for a
host-driven expansion claim. Species-balanced conditional label tail
is 0.0993; country-balanced is 0.0456; the latter gives very sparse
and sometimes source-determined countries the same weight as countries
with many records. Filtering to >=5 records (eleven source cases)
still leaves a large positive *artificial label reference* excess
(+10.28 pp), but those eleven were chosen posthoc by their own
observation depth and have no wild population independence.
The source-only tails cannot establish actual butterfly acceptance
or larvae surviving on recorded plants.

The fixed country host-degree graph is **unchanged across all seven
metrics**, including after recording-depth restrictions. Graphs may
connect a butterfly to plant species that were never its host, so
this reference concerns structure in the **recorded host identity
network** and spatial observation coverage, not ecologically plausible
transfers of larval feeding.

## Dependence on geography and source structure

The equal-case +6.86 pp conditional mean departure is concentrated
in some national settings. Without Cuba (four cases), the difference
is **+5.17 pp** and the artificial reference upper tail **0.0995**;
without the Dominican Republic (four cases), the corresponding
difference is **+5.92 pp** and upper tail **0.0762**.
The source-country examination was posthoc. The previously
completed record-weighted country omissions show opposing
geographic contributions: Cuba +13.25 indexed-record units
against its degree-fixed reference, Dominican Republic −9.28.
This can cause case-balanced and record-weighted conclusions
to diverge while both reproduce exactly.

Four of eleven national graphs have only one allowable degree-fixed
host network, and only six countries vary the record-weighted
statistic under rearrangement. These are not eleven independent
natural experiments.

## Biological boundary and next defensible target

The original 229-taxon country comparison still has **21
informative large-region strata versus a prespecified minimum 50**.
The 26-case subcountry analysis is selected on observed butterfly
record presence and an original source-added-host classification.
A host GBIF record in the same GADM1 region need not be a plant in
the same plot, established outside cultivation or accessible during
larval development. The host identity reference can assign genuine
non-hosts and the source records may share observers or multiple
records of one individual.

The **testable ecological question** remains when a newly distributed
host turns *potential access* into **developmental recruitment**.
Resolving it requires GBIF record-level host identities, dated
local plant establishment, independent larval use/development records
and suitable matched non-detection sampling effort. This archive
shows why reporting-grain and case weighting must be determined
before fitting any biological global analysis.

## Verifiable provenance

- [Full seven-metric, leave-country-out numeric archive](gbif_all26_weighting_stability_20261023.json)
- [Posthoc methodological transparency](../../protocols/GBIF_ALL26_WEIGHTING_SENSITIVITY_20261023.md)
- [Exact rational enumeration and graph convolution](../../../scripts/analyze_gbif_all26_weighting_stability_20261023.py)
- [Successful independent same-source CI run](https://github.com/zuizui0223/chocho/actions/runs/39918369972), with offline regression tests
- Original 26-cell [source botanical GADM facets](gbif_all26_host_gadm_source_20261023.json)
  and [48 original butterfly GADM facets](gbif_positive48_gadm_lifestage_facets_20261020.json)
- Frozen original GEB manuscript and anonymous submission supplement
  unchanged. The new audit is explicitly excluded from that snapshot.
