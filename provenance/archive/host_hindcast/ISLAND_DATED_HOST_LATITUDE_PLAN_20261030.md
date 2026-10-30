# Island dated-host and latitude-comparability implementation plan

> **For agentic workers:** execute the three independent, file-scoped tasks below using test-first implementation and independent review. The user explicitly authorized continued analysis without another design approval. Preserve the existing draft PR and all prior inputs/results.

**Goal:** determine whether the recorded-source distance association survives a basic latitude-zone comparison, and acquire actual pre-2018 island plant records for the fixed original host catalogue.

**Architecture:** preserve the frozen 96-island panel, 93 early-recorded butterfly concepts, 28 eligible targets and three risk masks. Repartition the existing exact comparison by longitude and latitude zone without changing outcomes or distances. Independently resolve the early host names to the same COL XR namespace and acquire complete count-only Plantae facets on the 28 targets for 2010–2017. Join those inputs with explicit unresolved/missing states.

**Tech stack:** existing Python 3.12, standard library, NumPy/SciPy, Shapely/geopandas, existing exact finite-table implementation and GBIF facet collector.

**Spec:** continuation requirements in `ISLAND_EXACT_NULL_HOST_ROSTER_RESULTS_20261030.md`, sections 3–5, plus this fixed specification. Starting commit `77d9e7295d10396254bcbf0994ae7a55a12b19d7`.

## Status before new measurements

The existing outcomes, ratios 0.488/0.550/0.588, 72.8% western-block concentration and an audit showing one target contributes about 46.2% of the primary distance contrast were already known. Its longitude block mixes tropical targets and targets at 46–54°N. This is explicitly exploratory ecological reanalysis, not a preregistration or independent replication. No focal plant taxonomic or occurrence results have been inspected. Only the COL XR kingdom match has been checked: Plantae key `P`, exact accepted KINGDOM.

## Global constraints and immutable sources

- Do not change prior source bytes, old exact/MCMC results, risk masks, source coast distances, original port/richness gates or the stopped 120-island baseline roster. Do not query later outcomes for that stopped roster.
- Original source keys are strings. Current COL XR checklist is `7ddf754f-d193-4cc9-b351-99906754a03b`; never join these to old numeric Backbone keys.
- All analyses concern indexed recording. Earlier focal nondetection is not absence; later recording is not first arrival, breeding or establishment. Static WPI2025 proximity is not historical shipping exposure.
- Plant queries use only 2010–2017, PRESENT, coordinates, no indexed geospatial issue and exact original island polygons; `limit=0` and empty results are mandatory. No biological coordinates or individual occurrence identifiers are acquired or persisted.
- Target roster: exactly the 28 rows with `fixed_target_recording_support=1` in the old target-support CSV. This roster depends on already known later all-butterfly richness; all conclusions remain conditional on that selection.
- Plant candidates: exactly the 1,040 `query_eligible=True` source plant IDs linked to at least one of the 93 early taxa. Keep all 1,062 early-linked plant IDs in coverage tables; the 22 infraspecific source names remain unresolved. Retain 93 butterflies and all 1,440 catalogue edges.
- Preserve partial taxonomic coverage. A matched host record confirms a catalogue-compatible recorded plant; a zero confirms only that no matched indexed host record was obtained under the query. Taxonomic or query failure never becomes a biological zero.
- Existing scientific reports remain in `provenance/archive/host_hindcast/`; all new source/derived products go to `data/exploratory/island_dated_host_20261030/` or `data/exploratory/island_latitude_null_20261030/`.
- No merge or reviewer comments. Publish verified continuation to existing draft PR #102 if the remote head remains the expected parent; stop on concurrent remote changes rather than overwrite.

| Source | SHA256 |
| --- | --- |
| `plant_taxon_query_candidates.csv` | `a1cc87b20fbd19415fa81ca4269b8414b0c32167bb6033ca0425354c89c1906f` |
| `original230_host_edges.csv` | `e55012d6b0eefd2b3e9e4b8f13fc4786ddc1e66c21f73e705bb20df92d0b5e22` |
| `baseline93_host_coverage.csv` | `f702ce482ca6b6411bb4530aa8b829400300cb5ef2464e367eef4e5a365ce9aa` |
| `recorded_source_target_support.csv` | `9bbf0d50f1a42b82879f23ef966a6de4a6dc36f5dc03cc92e5c511fc6bad3e74` |
| `expanded_island_panel_selection.csv` | `05c370d6fb63c623c48714ce14b377639c446d523bc742e4f93dfbb96d562538` |
| GSHHG 2.3.7 official full shapefile ZIP | `8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41` |

The original risk/margin/distance files must pass their existing pinned verification too.

## Review focus

1. Current accepted plant concepts, exact synonym handling, alternative homonyms and collisions must remain inspectable; unresolved names must not acquire species IDs by truncation.
2. Early-only plant filters and precise polygon hashes must be checked against actual outgoing requests, with complete facets, stable totals and first-page rechecks.
3. Tropical and hemispheric grouping must use geography only; old observed tables and source-distance scores must be unchanged.
4. Forced cells and finite conditional support must be reported; fewer possible tables cannot be described as ecological robustness.
5. Reused catalogue links and plant co-recording cannot establish introduced status or local feeding, and a port effect remains unestimated.

## Task 1 — conservative current plant taxonomy

Create `scripts/bridge_island_plant_taxonomy_20261030.py` and `tests/test_island_plant_taxonomy_20261030.py`.

- [ ] Test exact accepted species, exact synonym with accepted species, nonexact/higher-rank match, competing exact concept, wrong kingdom/family, failed request, source hash drift, and duplicate accepted concepts.
- [ ] Validate source pins; select 1,040 names by the fixed early93 link rule, sorted by source plant ID. Use official GET `/v2/species/match`, explicit checklist, scientificName, kingdom Plantae, source family, taxonRank=SPECIES, strict=true, verbose=true. At most three attempts per name, at least 1.05 s between global request starts. Preserve canonical request, raw taxonomy JSON and hashes in a deterministic archive; taxonomy responses contain no occurrence points.
- [ ] Accept only exact canonical matches, usage status ACCEPTED/SYNONYM, species rank for usage and accepted usage, accepted status ACCEPTED, a nonempty string accepted key, explicit Plantae classification and matching source family. Reject ambiguity diagnostics or a distinct exact competing accepted concept. Record all reasons. Multiple source plant IDs resolving unambiguously to one accepted concept are retained as source links and deduplicated per butterfly for occurrence joining, with a collision inventory.
- [ ] Export `plant_taxonomy_bridge.csv` (all 1,062 early-linked IDs, explicit unqueried source-name statuses), `plant_taxonomy_manifest.json`, and `plant_taxonomy_raw_responses.tar.gz`. Manifest counts all source/queried/failed/eligible rows and unique accepted concepts. Preserve infraspecific unresolved names without live matching in this continuation.
- [ ] Run meaningful tests, actual bounded acquisition, and source/count reconciliation. Partial queries remain explicit missing states.

## Task 2 — early plant recording on fixed target islands

Create `scripts/acquire_island_early_plants_20261030.py`, `tests/test_island_early_plants_20261030.py`, and `scripts/analyze_island_dated_hosts_20261030.py` with tests.

- [ ] Test the early-only filter, kingdom/checklist, count-only contract, target roster/hash checks, incomplete or unstable facet pages, taxonomy failure vs zero, and source-ID collision deduplication.
- [ ] Restore only the SHA-pinned official GSHHG archive, reconstruct original dissolved polygons and verify all target CCW WKT SHA256 values. WPI need not be reacquired.
- [ ] Acquire all Plantae species/basis/dataset/year facets for the fixed 28 targets in 2010–2017. Reuse the existing complete facet validation, 100 bins/page, maximum 100 pages/island, stable counts, zero-result contract and first-page recheck. At least 1.5 s between request starts, three bounded retries; preserve raw count-only request/response/receipt triples in a deterministic archive. No later plant queries or geometry simplification.
- [ ] Export complete/missing island query statuses, species counts, all-plant recording totals, unassigned fractions, and number of recorded years/datasets. Record the inherited botanical support screen (>=20 records, >=90% species assignment, >=3 positive years) only as a support descriptor; it does not replace the original target roster or make absence claims.
- [ ] Join actual current accepted plant keys to each original93 catalogue. On the old 2,385 primary cells export confirmed early host concepts and records, unresolved catalogue counts, complete-query status and later butterfly record status. Do not sum duplicate current concepts more than once. Retain complete no-match, taxonomically incomplete no-match and query-missing states distinctly.
- [ ] Summarize dated co-recording descriptively for all candidates and later-positive cells; no causal host-effect or introduced-host claim. Do not add a favorable subset, tune recording thresholds or create an unsupported absence endpoint. A new host-effect model requires an explicit later specification after inspecting available conditional support.

## Task 3 — latitude-zone-conditioned distance comparison

Create `scripts/analyze_island_latitude_null_20261030.py` and `tests/test_island_latitude_null_20261030.py`.

- [ ] Test zone boundaries, hemisphere distinction, identical risk masks/scores/outcomes, row-margin refinements, exact tiny-table enumeration, and independent table mask/margin validation.
- [ ] Fixed zones from each island's already archived representative latitude: south extratropical if latitude < -23.5 degrees; tropical if -23.5 <= latitude <= 23.5; north extratropical if latitude > 23.5. The conventional rounded tropical boundary is a geography control, not a measured climate classification.
- [ ] Refine each original species-by-longitude row into species-by-longitude-by-zone rows and retain target-island margins. Use the same three original scenarios and exact DP completion counts. Cap at 100,000 memoized states and 2,000,000 attempted transitions per block; incomplete runs yield no qualified tail.
- [ ] Preserve the original mean log(1 + upper coastline-gap) statistic including forced positives. Compute exact conditional means/variances with floating scores, every cell's exact occupancy and block contributions. Compare with the old ratios without treating attenuation as a causal proportion explained.
- [ ] Use exactly 100,000 independent uniform complete-table rank draws per scenario, with replacement, seeds 2026103011/2026103012/2026103013. No reruns or budget extensions after viewing results. Report inclusive closer-tail count and existing individual/simultaneous three-scenario binomial limits. If the full fiber is forced, label no randomization contrast and do not call it support.
- [ ] Preserve first-run lock, all sampled ranks/statistics, source/script/plan hashes and actual execution times. Report variable cells/positives/species/islands and forced cells. Do not pick a favorable zone or post-hoc climate cut point.

## Integration

- [ ] Review code and independent scientific counts; verify task-specific tests and one final full project suite.
- [ ] Exclude all new exploratory code/tests from the existing anonymous manuscript bundle and verify the unchanged 222-file SHA256.
- [ ] Save an integrated result report, a clear static scientific figure, source receipts and an execution ledger; publish verified files to draft PR #102 and verify checks on the exact published head.

**Stop conditions:** original source drift, missing/unstable source response, count-only contract failure, exact DP cap, or insufficient randomization support must remain an explicit limitation. Finish independent authorized work and report the precise unresolved input. Do not relax a frozen source or scientific gate to produce a positive result.
