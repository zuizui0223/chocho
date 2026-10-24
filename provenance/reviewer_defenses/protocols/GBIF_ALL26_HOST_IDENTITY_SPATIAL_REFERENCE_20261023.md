# Conditional country botanical-host-identity spatial reference — frozen before enumeration

**23 October 2026; posthoc and intentionally nonbiological label reference.**
The observed 26 source-added butterfly×country candidate-host GADM1
panel (44 botanical country queries, 24 accepted WCVP plant taxa) is
already known. This plan fixes an **additional spatial-label reference**
before calculating its expected statistic or conditional tail.

## Exact fixed source arrays

Use ONLY source-verified
`provenance/reviewer_defenses/results/gbif_country_20261020/gbif_all26_host_gadm_source_20261023.json`, the earlier source-only
`gbif_all26_source_host_country_inventory_20261023.json`, and
the original all-stage 48 positive butterfly GADM1 facets. There
are 26 source-added butterfly×country positive cells from 16
butterfly species across 11 countries, 46 unique source
butterfly×country×plant links and 44 unique plant×country
queries. Original butterfly records in 26 cells have 453 mapped
GADM1 observations and **231** in GADM1 units with a candidate
source plant observation. All GBIF facet records are source-only,
not independent population samples.

## Matched country source-label null

Within each country, keep fixed:
1. The list of butterfly taxa with records, all per-species GADM1
   butterfly record counts, and its number of **source-defined
   candidate introduced hosts** (butterfly row degrees).
2. The list of source-listed candidate host taxa in the same
   country, each plant's GBIF GADM1 observation footprint and
   number of butterflies that list it (plant column degrees).
3. All no-GBIF-plant-index records as zero *source detections*
   without assuming botanical field absence.

**Enumerate all simple bipartite butterfly×plant graphs with the
same exact row and column degrees** in each country. Assign each
such graph equal reference mass. All countries are conditionally
independent under this expressly artificial label reference.
Do not shuffle a plant to another country, change the focal
butterfly or plant GADM1 counts, make a duplicate host edge,
or change native/introduced source classification.

For a graph, mark each focal butterfly GADM1 unit as
`source-host-indexed` if ANY plant currently connected to
that butterfly has >=1 original GBIF plant record in that
same country GADM1. Let `T` be the sum of focal butterfly
GBIF indexed records in all such marked GADM1 units
over all 26 butterfly×country source-positive cases.
The fixed denominator is **453** mapped butterfly records.
The observed `T` is **231**; this was known before
defining the label reference. Compute exact conditional
expected `T`, full null distribution by country-wise
integer-PMF convolution, and upper-tail P(T>=231).
These are **NOT biological p-values**.

Also report the exact number of possible degree-matched
graphs in each country, the count with >1 possible graph,
the fraction of 26 butterfly cases in those countries,
source-record weighting sensitivity, and leave-one-country
descriptives. Do not select favorable countries after
viewing the output. A country with one possible graph
adds constant T but no reference information.

## Interpretation and stop rules

The swaps may connect a butterfly to a botanical plant that
is NOT its documented host. Thus the reference asks only:
**does the source-listed host identity align with the
GBIF observation geography more than arbitrary link
identity under these margins?** It does not ask whether
true butterflies choose or survive on host plants.
Host associations can share plant family and source
quality, with GBIF cultivated-plant and observer bias,
strong human attention and repeated records.
There is no host-use, larval survival, demographic
establishment or native-versus-introduced paired
comparison in these data.

If enumeration of any country's graphs is computationally
intractable (>200,000 states), report NOT_IDENTIFIABLE
and do not silently use an unverified MCMC replacement.
If all country distributions have one state, do
not report an inferential tail.

Results are an exploratory data-geometry negative control
only. The 229-taxon country macro effect and frozen GEB
submission/review data are not changed.
