# Island and port measurement continuation — 2026-10-29

**Goal:** turn the source-verified island/port inventory into measured geographic exposures and determine whether the existing butterfly panel can identify the proposed host × island interaction.

**Spec:** `ISLAND_ISOLATION_PORT_SOURCE_GATES_20261029.md`. This is a continuation of the user's requested analysis; the original GEB manuscript remains a separate frozen deliverable.

## Decisions before calculation

1. Reuse the SHA-locked GSHHG 2.3.7 high-resolution L1 source and all 3,689 source-locked WPI 2025 entries. Use source latitude/longitude attributes, not the problematic transformed ArcGIS geometries.
2. Retain the existing ≥5 km² island universe and ≥20 km² subset. Check five continental seeds and named island seeds against dissolved source geometry. Match the five seeded continental landmasses to the prior >7 million km² census; fail on disagreement. Antarctica is absent from this L1 source, so distance means distance to the five mapped continental landmasses, not every continental coast. Flag islands for which the south-of-60°S cap could be nearer; exclude these from the primary support summary.
3. Compute minimum **spherical great-circle coastline separation**, sampling each original coast segment at ≤1 km. Report the sampled upper estimate and a conservative discretization lower bound (1 km for two sampled coastlines; 0.5 km for a point and one coastline). These bounds cover sampling only, not source shoreline or spherical-Earth error. Remove artificial ±180° cut edges. Never substitute island centroids.
4. Count WPI registered entries within 50 and 100 km of island coasts, carrying threshold ambiguity through certain/possible counts. A nearby entry may be on another island or the mainland. The exposure is proximity to a mapped port location, not traffic, an on-island port, or a dated arrival pathway.
5. Audit the frozen biological panel at its actual country grain. Preserve species × Level1 strata and target-group recording effort. Record zero informative overlap rather than fitting an unsupported interaction. Country codes do not define independent physical islands.

## Work and verification

- [x] Add and run independent geometric fixtures: dateline distance, omitted cut edges, long-segment coverage, known coast-gap interval, and port threshold ambiguity.
- [x] Implement `scripts/measure_island_port_geography_20261029.py`; write derived island-level exposure table, aggregate support table and source/QA manifest under `data/exploratory/island_port_20261029/`.
- [x] Run the independent country-panel support audit and record source checksums, exclusions and exact denominators.
- [x] Compare seed identities and previous 8,260 / 3,696 island counts; summarize geographic combinations useful for subsequent data acquisition. Do not label this geographic result a butterfly effect.
- [x] Exclude new exploratory scripts/tests from the existing GEB reviewer snapshot, run the relevant tests and the repository suite, independently review code/results, then prepare the existing exploratory branch update.

## Rejected shortcuts

Country-level port averages cannot stand in for distances of butterfly observations. Positive alien-establishment cases do not provide establishment failures. Present port locations cannot be silently backdated. Regional host status cannot be copied to individual islands as verified local plants. A marginal island-country detection enrichment cannot estimate a within-species host interaction.

## Remaining biological test

Use island-resolved, dated butterfly recording and target-group backgrounds, independently verified host occurrence and an explicit mainland reference. For port analyses also require within-island locations and a coastal pseudoport comparison. Hold out physical islands/archipelagos and do not count countries sharing a landmass as independent islands.

## Bounded island recording acquisition, fixed before outcome queries

Open a new broad Papilionoidea recording panel using the verified GBIF COL XR namespace (`checklistKey=7ddf754f-d193-4cc9-b351-99906754a03b`, `taxonKey=5G9`). This is distinct from the older numeric-backbone 239-species panel; do not merge species keys without a verified taxonomy bridge.

- Candidate physical islands: source area 20–20,000 km², source representative latitude within ±55°, and no unresolved Antarctic comparison.
- Twelve prespecified geographic combinations: area 20–100 / 100–1,000 / 1,000–20,000 km² × direct continental gap certainly <100 / ≥500 km × WPI entry certainly within 50 km / all listed entries certainly beyond 100 km.
- Select up to two islands per combination: largest island first, then the spatially farthest candidate from that island. Selection uses geography only; freeze the complete roster and source checksums before fetching butterfly outcomes. Do not replace low-record islands after seeing counts.
- Query exact source island polygons for 2010–2017 and 2018–2025, including descendant Papilionoidea, coordinates present, no indexed geospatial issue and occurrence status PRESENT. Request aggregate species, basis-of-record, dataset and year facets, with complete independent pagination and index-stability checks.
- Retain query failures as missing, true successful empty responses as zero indexed records, and higher-taxon identification residuals explicitly. No simplification or rectangular facet substitution is allowed.
- Current retrospective GBIF indexing is not a reconstruction of data available in 2017. The 2025 WPI exposure is not established as pre-2018 history. This acquisition opens an independent observation panel; it does not estimate colonization, host mediation or port causality.

## Targeted barrier-definition check

For the eight named sanity-check islands, additionally compare distance to a continental coast with distance to the nearest other source landmass ≥5 km². Keep this local-neighbor measure distinct: the neighbor is not necessarily a butterfly source or a functioning stepping stone. Use spherical cap lower bounds only for safe candidate pruning; calculate final distances from coast samples, not landmass centers.

## Execution notes

The initial Okinawa sanity-check point (128.0°E, 26.5°N) was outside the source shoreline. The check failed before measurement. Direct source inspection showed both (127.9°E, 26.5°N) and (127.8°E, 26.35°N) inside the same 1,202 km² source island; the first is now used. This changes the check location, not the island geometry or analysis universe.

Independent review identified two requirements now covered by regression tests: verify the complete WPI page digest against prior successful Actions artifact 12890257067 (not just metadata/IDs), and bound Antarctic proximity using the actual spherical coast samples minus half the sampling step. Inspection of that prior artifact also exposed one numeric (0,0) placeholder. All 3,689 original entries are audited; only 3,688 enter distance calculations. The WPI source digest independently matched the prior Actions output exactly.

The bounded 24-island pilot subsequently completed all 48 windows without replacing cases. Only 10 islands have ≥3 recorded species bins in each window, and only two area×gap blocks retain such coverage in both port classes. The independent named-neighbor measurement also completed all eight islands. See the result document for exact values and the remaining biological inference limits.
