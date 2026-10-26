# Source-verified butterfly / introduced-host point proximity: 26 original cases

**Date:** 2026-10-26  
**Scope:** Exploratory **indexed butterfly–accepted candidate larval-host plant point co-recording**, not observed larval resource use, demographic range expansion, ecological competition, realized climatic niche expansion or causal effects of introductions. The GEB submission PR #38/main and the distinct climate-novel country analysis PR #72 remain unchanged.

## Scientific finding and exact acquisition

The earlier WCVP/HOSTS/BCE plant-range reconstruction shows potential resources in country/botanical-region units. The decisive new question here is whether actual public indexed observation POINTS of the butterfly and an accepted source candidate plant are near each other in both space and the index year.

**Fully passed independent real-source workflow:** [GitHub Actions #40344793069](https://github.com/zuizui0223/chocho/actions/runs/40344793069), seven tests passed. [Full original numerical receipt and 70 query hashes, aggregated only](https://github.com/zuizui0223/chocho/actions/runs/40344793069/artifacts/12661970857).

- **70/70 original source-accepted taxon×country GBIF record retrievals COMPLETE**, 26 butterfly queries + 44 unique accepted source plant queries. Page size 300; 80 requests; no transport retries; no original vs new facet total-count drift across all 70. The new point index 2010–2025 source total is **4,062**: exactly **467 butterfly records and 3,595 unique plant×country query records**.
- Exact original source cohort: 26 butterfly-positive introduced-host opportunity country cases, 16 butterfly species / 11 countries, 24 source candidate plant species, 46 butterfly×plant×country links. The earlier 48 country-positive contrast's 22 non-added cases are *not* known nonhost or butterfly absence controls.
- All 26 cases have a complete source snapshot and all 467 butterfly indices have validated accepted source species keys, country, year and valid point coordinates under the frozen filter.
- **75/467 (16.06%)** butterfly indices are <=5km from an indexed source candidate accepted host plant within **±3 source observation years**. The other 392 are NOT verified negative host-use outcomes.
- **14/467 (3.00%)** indices <=1km and ±3 years; **160/467 (34.26%)** <=25km and ±3 years. Same-year <=5km only **16/467 (3.43%)**; ignoring observation year <=5km is **98/467 (20.99%)**. These are alternative *detection definitions*, NOT seasonal resource or actual plant occupancy contrasts.
- **288/467** butterfly indices have reported coordinate uncertainty <=1km. Restricting *both* the butterfly and botanical source point to documented uncertainty <=1km produces **34/288 (11.81%)** <=5km and ±3 years. The strict denominator is only qualifying butterfly indices, not comparable survey effort across all countries.
- **22/26** source cases have >=1 candidate botanical indexed point; **4/26** have no accepted candidate plant GBIF point record under frozen filters. Missing botanical GBIF indices do NOT prove absence of the host plant.
- The source-dependent text screen found **0 literal juvenile-stage + exact accepted candidate-host binomial + action phrase annotations**. This means *no source-confirmed larval host use from that restrictive parser*, NOT zero actual feeding.

## Taxonomic and geographic concentration of the observed point matches

**Two butterfly×country source cases contribute 60 of the 75 primary near-pairs (80%):**

- *Papilio demoleus*, Dominican Republic: **37/199** butterfly indexed records have a candidate accepted plant point within 5km/3years.
- *Colobura dirce*, Venezuela: **23/56** indexed records meet the same threshold.

Those two cases already account for 255 of the 467 butterfly indices (54.6%). After removing both, near-host index co-recording falls to **15/212 = 7.08%**. These are not independent within-case butterfly populations, and this descriptive sensitivity is not a source-preserving geographic null test.

Across the whole fixed 26-case cohort, only **8 cases** contain >=1 positive <=5km/3year co-recording. This heterogeneity cautions against treating the 16.1% pooled point-index rate as a global ecological species response.

### A fully prespecified sparse-source pilot that cannot substitute for the full cohort

The separate [12-case low-volume source-gate](BUTTERFLY_POINT_LOWVOLUME_DIAGNOSTIC_PROTOCOL_20261026.md) preserved all frozen butterfly-positive country cases with <=20 original butterfly record indices and <=20 total source candidate botanical record indices *before inspecting point distance*. Its [source-complete result](BUTTERFLY_POINT_LOWVOLUME_DIAGNOSTIC_RESULT_20261026.md) was **0/46** within 5km/3years, but **10/46** within 25km/3years, across 30 complete original query identities. Because it was deliberately selected on low record volume, **0/46 is not an independent counterexample to the full 75/467 nor a population ecological absence result**. Instead, it shows that the source/index density strongly determines the resolvability of potential host access at local scales.

## Claims retained and rejected

**Retained:** the original macro host-resource region expansion is real at its WCVP-defined potential-resource scale; it does not automatically establish a close source-verified plant point at a recorded butterfly point. A complete, checksum/identity-controlled point comparison demonstrates a measurable loss of ecological identifiability when moving from broad regions to local plant co-recording. The source point-proximity distribution is strongly concentrated by species×country and indexing intensity.

**Rejected:** actual butterfly climate niche entry, proven larval use of a local introduced host, an increase or decrease in butterfly populations, a species-level host-plant benefit, a competition/interaction strength shift, or a claim that botanical record scarcity means true local host absence. GBIF co-recordings can reflect shared observer preferences, garden/cultivation sources, platform coverage and multiple indexed records of nonindependent sightings.

**Next decisive biological comparison:** apply validated butterfly point geography to the *distinct* 73 positive source country observations that were included in the prespecified climate-analogue vs climate-novel PR #72 test; classify point-level climate (not only host-range country climate), record known host geography and temporal overlap, and verify host-linked juvenile observations or repeated local recruitment. Correct for sampling effort and compare against non-host/other-plant geography using a source-locked structural null. If such endpoints remain unavailable, publish potential host-resource *niche opportunity*, not butterfly realized niche evolution.

## Source code, audit and publication firewall

- [Prespecified protocol](../../provenance/reviewer_defenses/protocols/BUTTERFLY_HOST_POINT_PROXIMITY_20261026.md).
- [Reproducible public GBIF point retrieval and source-only analysis](../../scripts/analyze_butterfly_host_point_proximity_20261026.py).
- [Regression tests](../../tests/test_butterfly_host_point_proximity_20261026.py).
- [Green source-run](https://github.com/zuizui0223/chocho/actions/runs/40344793069), [derived result artifact](https://github.com/zuizui0223/chocho/actions/runs/40344793069/artifacts/12661970857).

No raw point coordinates, GBIF occurrence record bodies, observers or locality identifiers are included in the retained derived artifact. The source remains a *live-index* retrieval rather than a DOI-stable full GBIF occurrence download; the archived SHA receipts and method make the audit inspectable but do not guarantee future point-exact replication.
