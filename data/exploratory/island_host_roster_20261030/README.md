# Recovered original butterfly host roster

The verified original accepted HOSTS table supplies **2,524 butterfly–plant-ID links for all 230 conservative butterfly concepts**. Every original 239-species host count reconciles with the pinned source metrics. The 93 butterflies with early records in the fixed 96-island panel retain all **1,440** expected catalogue links. These are source-catalogue candidates, not island plant occurrences or observed host use.

| Catalogue scope | Butterfly taxa | Links |
| --- | ---: | ---: |
| Original count reconciliation |239|2,596|
| Conservative taxonomy bridge |230|2,524|
| Early 2010–2017 recorded subset |93|1,440|

There are 1,677 distinct source plant IDs; 1,642 pass the source-binomial name gate. The 35 withheld IDs all contain infraspecific rank markers (17 varieties, 16 subspecies, 2 forms). They remain in the catalogue, with explicit unresolved-query status. No name was promoted to a species by truncation. All 93 early butterflies have at least one eligible source query name, and 71 have no unresolved source names. The early subset contains 1,062 plant IDs, including 1,040 with eligible query names; 1,411 of its 1,440 links pass that source-name gate.

Current plant taxonomy has not been queried. Query eligibility does not establish current accepted GBIF identity, plant presence on an island, or larval host use. No new occurrence calls were made. The eligible butterfly roster is all 230, not the 105 species selected by early-or-later occurrence; early coverage uses only 2010–2017. No BCE-augmented country-case edges were used.

## Source and execution

The [successful source recovery](https://github.com/zuizui0223/chocho/actions/runs/40979986160) ran commit `ef86ac353cce571c2765fb5f8ebe49e18e06ca6e` and exported artifact `12943834119`. This derived ZIP was 124,608 bytes, SHA256 `414f5806034d1934ded8245efbe3a30edb81adaa1ac05c5f03956831a39c3200`. Its bytes, exactly five permitted filenames, and every manifest output hash were checked after download.

The original author artifact remained inside the Actions runner. The derivation verified the pinned inner ZIP, all 28 source-manifest members, host/metric/bridge/facet hashes, source metadata and each original species count. It did not extract the raw source members or execute archived code. The outer original artifact digest was verified by the download action; the derivation manifest accurately labels its own outer check as metadata-based. The raw author archive was not retrieved locally or included in derived outputs.

The initial import failed because a digit-only ID parser rejected an actual scientific-notation ID. A diagnostic-only retry found four such rows across the full 53,449-row source, including one among the original cohort's 2,633 source rows. The corrected parser checks exact positive integrality with Decimal and preserves the original ID string, without floating conversion, numerical normalization or merging. The all-239 count gate still applies. Both incomplete attempts are retained below.

## Files

| File | Role |
| --- | --- |
| `original230_host_edges.csv` | Every eligible butterfly–accepted-plant-ID source link, with butterfly bridge evidence and botanical name gates |
| `original239_host_count_alignment.csv` | All original per-species expected/recovered counts and the nine excluded butterfly concepts |
| `plant_taxon_query_candidates.csv` | Unique source plant IDs, retained accepted names/families, linked butterflies, and exact-name query eligibility |
| `baseline93_host_coverage.csv` | Early-only butterfly recording coverage and recovered/query-eligible/unresolved host counts |
| `original230_host_roster_manifest.json` | Original source identities, successful execution commit, counts, taxonomy audit and output hashes |
| `derived_artifact_retrieval.json` | Verification of the successful five-file artifact after download |
| `derived_host_coverage_summary.json` | Additional counts calculated solely from the recovered derived CSVs |
| `source_attempt_1_incomplete_receipt.json` | Original strict-parser failure, preserved from its one-file artifact |
| `source_attempt_2_identifier_diagnostic.json` | Aggregate-only identifier-format diagnostic from the retry's logged receipt |

Tests are in `tests/test_original230_host_roster_20261030.py` (59 passing tests). Replaying the source derivation uses `.github/workflows/island-original230-host-roster.yml` and the still-available pinned original artifact. Source expiry or any hash/cohort mismatch must produce an incomplete result, never a fabricated zero or silent source replacement. The whole continuation is described in the [integrated result report](../../../provenance/archive/host_hindcast/ISLAND_EXACT_NULL_HOST_ROSTER_RESULTS_20261030.md).
