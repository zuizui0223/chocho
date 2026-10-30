# Task 1 correction — independent review

## Current disposition

**Scientific correction and final full-output verification PASS; no blocking defect found.** Acquisition, corrected offline processing, full source/response/decision reconciliation and archive verification are complete. A minor request-timestamp verification limit is retained below. No shared repository file, acquisition cache or remote system was changed by this reviewer; no HTTP requests or subagents were used.

Reviewed repository: `/workspace/scratch/10d2ad986b7b/chocho`.

Reviewed sources: `scripts/bridge_island_plant_taxonomy_20261030.py`, its test file, the fixed Task 1 brief, the earlier independent review, actual cached taxonomy JSON/receipts, the retained original acquisition script, regression red output and final full-suite transcript.

## What changed and why it is correct

The only script difference between `/tmp/chocho-dated-workflow/task-1-acquisition-script-original.py` and the final repository script removes the requirement that an alternative's **accepted target** have rank `SPECIES` before that alternative is counted as a competing concept. The same-canonical alternative usage must still be an exact match of rank `SPECIES`, and its accepted key must differ from the selected primary accepted key.

The fixed plan says to reject a distinct exact competing accepted concept; it does not exempt alternatives whose species-level name usage resolves to a subspecies or variety. Such a response does not disambiguate the source binomial. Rejecting it as unresolved is therefore a correction to the intended conservative gate, rather than a new favorable specification. The final primary accepted match must still have species rank and all the previous status, namespace, kingdom, family and exact-name gates. The correction does not promote lower-rank target concepts to eligible species.

The request URL, primary query roster, cache contract, retry code, accepted-name verification request selection, source-edge preservation and output schema are unchanged. Exact synonyms that return to the same accepted key remain allowed. The scientific rejection state is preserved separately from HTTP success and must not become an occurrence zero.

The existing reason name `competing_exact_current_species_concepts` now includes a species-level **usage** whose accepted target has a lower rank. The stored alternative keys and raw diagnostics make this inspectable. In prose, describe it as competing exact name interpretations or distinct accepted concepts; not every case has to be called a botanical homonym. This wording consideration is not a blocker.

## Actual cached examples

The raw response bytes were checked against their receipt SHA256. Old and new `assess_match` results were compared on the same cached response without any request or cache mutation.

| Source ID / name | Primary accepted key | Distinct alternative accepted target | Old eligibility | Corrected eligibility |
| --- | --- | --- | --- | --- |
| `2700274` / Cardamine pratensis | `R4LD` | `5G7NW`, Cardamine pratensis atlantica, SUBSPECIES | True | False |
| `2853709` / Hosackia crassifolia | `3MQN7` | `9F343`, Acmispon glaber glaber, SUBSPECIES | True | False |
| `2744513` / Crotalaria saltiana | `ZP3K` | `5NCHV`, Crotalaria pallida pallida, VARIETY | True | False |
| `2337827` / Lepidium campestre | `6P99B` | `7K8RF`, Lepidium hirtum hirtum, SUBSPECIES | True | False |
| `114913` / Livistona chinensis | `3VP7Y` | Alternative usage returns to the same `3VP7Y` species concept | True | True |

All the competing alternative usages above are exact same-canonical `SPECIES` usages. Several are `AMBIGUOUS_SYNONYM`; the gate correctly does not treat lower ranking as disambiguation. The Livistona example confirms that usage multiplicity alone is not treated as multiple accepted concepts.

Selected source-body hashes:

- Cardamine pratensis: `84f623e8935897010b2f4876491ca19ef626d5d67be863b7b8021b7123889505`.
- Hosackia crassifolia: `95d2f81e022678dcd5b5094c203df4a2861b726150e4ece0886cb3cb938a3fe2`.
- Crotalaria saltiana: `e8096510fbbaf63e904623ac18b1f5b38d584e6ecd53dd73c2568a857562d73a`.
- Livistona chinensis: `161b833d8fc10a8e59f4666617ed7a47123e60bf0d6b477d51e1daf040a96cbe`.

At the first correction audit snapshot, 815 primary receipts and no accepted-verification receipts were present. Sixteen of those cached responses had expanded competing-key inventories under the final gate; these are an interim observation, not final eligible-count results. At the later example check, acquisition had reached 844 primary receipts and still no verification receipts.

## Tests and acquisition-versus-processing provenance

The regression transcript records the expected failures for both new lower-rank accepted-target cases: **2 failed, 18 passed**, specifically because the old gate returned eligible=True. The final full-suite transcript records **631 passed, 20 warnings, 5 subtests passed in 13.97s**. The tests were inspected and their retained evidence reused; the reviewer did not repeat the project suite.

Acquisition script SHA256: `5f26c60c1a2b743b4463e4e1744cee20f2ca8f4e923c1a7cd8f2bcb774862f75`.

Final processing script SHA256: `53410a0a99843a36967a96cac3a4c9a0ee4ff5fd6aeecb67dfe8f7f25745e9f6`.

The final output retains these as separate acquisition and processing identities in `plant_taxonomy_audit.json`, and the original acquisition script is archived as `plant_taxonomy_acquisition_script.py.txt`. The client's request-generation contract did not change. Completed corrected processing reused the successful primary/verification receipts, with cached URLs and body hashes verified. The following final checks replace the earlier pending state.

## Final full-output verification

The reviewer independently read the finished source tables, bridge, collisions, manifest, all taxonomy receipts and bodies, and every archive member. A separate decision implementation, without importing the production assessor, recomputed every primary response's final eligibility, rejection-reason set, competing accepted-key set, accepted status and accepted-status evidence. It also applied the corrected gates to each required accepted-name verification response. All fields agree with the final bridge.

| Check | Verified final result |
| --- | --- |
| Original cohort | 93 butterflies; all 1,440 catalogue edges retained |
| Bridge source IDs | All 1,062 distinct early-linked plant IDs, with exact set equality to the source roster |
| Primary queries | Exactly 1,040 SUCCESS responses; exact official URLs, checklist, source name/family, rank and strict/verbose parameters |
| Unqueried source names | Exactly 22; no corresponding primary receipts; explicit unresolved status |
| Accepted-name verification | Exactly 29 required requests, all SUCCESS; raw responses and URLs verified |
| Total HTTP acquisition attempts | 1,069; all corresponding receipts and raw-body hashes reconcile |
| Final conservative eligibility | 695 source rows and 695 unique accepted concepts |
| Eligible accepted-concept collisions | 0 groups; all per-row collision annotations independently reconcile |
| Accepted-status evidence across all primary rows | 1,011 primary responses; 24 independently verified accepted-name matches; 5 still unverified |
| Catalogue edges with eligible taxonomy | 865; remaining 575 unresolved; the sum remains 1,440 |
| Butterflies with at least one eligible concept | 84 of 93; the other 9 remain explicitly without resolved query concepts |
| Archive | 3,207 distinct sorted members, zero gzip/tar timestamps, all member bytes identical to the complete cache |
| Identity | All pinned source and manifest output hashes, processing script/test hashes, archived acquisition script and audit-manifest hash match |

The final rejection-reason counts independently reproduce: competing exact concepts 339, unqueried source names 22, accepted status unverified 5, non-species rank 2, source family unverified 2, nonexact match 1, submitted canonical mismatch 1. Reasons may overlap; these are not disjoint source-row denominators. HTTP success does not imply scientific eligibility.

The nine butterflies without an eligible accepted host concept are Araschnia levana, Carcharodus alceae, Colias croceus, Colias interior, Erynnis tages, Maniola jurtina, Polyommatus icarus, Pseudozizeeria maha and Wallengrenia otho. They must remain unresolved in the occurrence join, not be recoded as host absence. This set is unchanged by the correction.

## Verified acquisition and processing separation

Reviewed `/tmp/chocho-dated-workflow/reprocess-plant-taxonomy.py`, its strictness-delta receipt, saved original bridge/manifest and the final `plant_taxonomy_audit.json`. The replay explicitly replaced the client's network opener with a function that records and rejects network attempts. Its completed receipt records **zero attempted HTTP requests** and identical before/after cache file hashes.

Independently comparing the saved original and final bridge shows:

- All 1,062 source IDs remain in the same order.
- Every primary and accepted-verification request URL, response hash, query status and attempt count is unchanged.
- The competing-key inventory changes for 21 rows; eligibility changes for exactly 10 rows, always True to False.
- Eligible source rows decrease from 705 to 695. Eligible catalogue edges decrease from 878 to 865; all source edges remain present, so unresolved edges increase from 562 to 575.
- No previously rejected row becomes eligible, and the count of butterflies with at least one eligible concept remains 84.
- The acquisition and final raw archive hashes are identical: `24805b7161711538e30f443b5d784efad1b9038a6036a710b4eb482b0dca7039`.
- The archived acquisition script is byte-identical to the retained original script; final processing uses the distinct final hash recorded above.

The ten newly unresolved source names are Lepidium campestre, Passiflora biflora, Scabiosa atropurpurea, Boehmeria cylindrica, Cardamine pratensis, Crotalaria capensis, Crotalaria saltiana, Hosackia crassifolia, Petroselinum crispum and Prunus ilicifolia. Their original accepted keys and new competing keys remain in the bridge and correction receipt.

Actual acquisition timestamps remain in the original receipts, from `2026-10-30T02:06:03.698567+00:00` through `2026-10-30T02:24:45.298627+00:00`. Final manifest timing refers to cache processing and is not silently substituted for those original retrieval times. The audit explicitly preserves both identities and the correction history.

### Timestamp precision limit

The client spaces monotonic request-start markers by 1.05 seconds, but only the subsequent UTC timestamps were archived. The independently observed minimum UTC gap is 1.049838 seconds; the audit records six intervals up to about 162 microseconds shorter than 1.05 seconds. Consequently the strict monotonic threshold cannot be independently certified from UTC receipts alone. This limit is already stated in `plant_taxonomy_audit.json`; it does not alter taxonomy decisions or the cache-only correction and is not a scientific blocker. Do not rewrite those source timestamps or claim an exact all-UTC-intervals pass.

## Final disposition

**No required scientific or code fixes remain in this correction.** Its reduction in coverage is the appropriate conservative response to unresolved name interpretation. Carry the complete source roster, explicit unresolved rows and 695 eligible concepts into the subsequent host-record join; preserve the distinction between taxonomy eligibility, local indexed plant occurrence and demonstrated larval use.
