# Posthoc weighting of exact source host-identity spatial reference — 23 October 2026

**This audit is POSTHOC.** The 26-case original record-weighted
result (+0.2023547 pp), original case-balanced expectation (+6.8580651 pp)
and the shared 13,826,598,912-degree-constrained possible source-label graphs
were **known before** this analysis. No result should be described as a
pre-registered ecological null test. The exact tails are probabilities over
artificial, biologically unconstrained host-label reallocations.

## Fixed evidence and computational contract

Reuse the completed **26 butterfly×country, 16-species, 11-country**
source-eligible observations, the 44 GBIF plant-country GADM facets and
original 48 butterfly-country GADM facets. Keep original butterfly
counts and plant footprints, butterfly row candidate-host degrees, plant
column host degrees and country identities unchanged. Exhaustively enumerate
each country's simple bipartite graph; verify 13,826,598,912 possible joint
graph combinations, then compute its exact conditional score distribution
through integer-lattice country convolution, not MCMC.

Compute all seven **descriptive, different-estimand** summaries on the
same source rows/graphs: record-weighted, equally weighted butterfly×country,
equally weighted butterfly species (average countries within species),
equally weighted countries (average butterfly cases within country),
and equally weighted source cells with >=2, >=5, >=10 indexed butterfly
records, while keeping all source graph degrees fixed even when a case
is not included in a threshold summary.

For each display observed share, exact mean under artificial label reference,
percentage-point difference, support size, graph count and upper tail
P(reference score >= observed). Exact tails **must never** be called real
ecological/causal p-values or a measure of evidence for actual host use.
Posthoc leave-one-country descriptive effect for all 11 countries; exact
country-deletion tails for CU and DO, with retained-case denominator
recomputed (not the original 26) and no preferential dropping of cases.

### Explicit interpretation limit

The null can assign non-host plants as larval hosts; a positive departure
can reflect botanical trait groupings, recording coverage and spatial
data-geometry. Lack of difference does not mean genuine host choice is
unimportant. GBIF occurrence records are not independent insect populations,
and the original 21 informative country-level strata remain short of 50.
Nothing here establishes locally eaten plants, survival, reproduction,
host arrival chronology, introduced-plant benefits or population expansion.

Preserve this outside the frozen GEB v0.2 anonymized reviewer package,
and do not edit the original GEB submission manuscript.
