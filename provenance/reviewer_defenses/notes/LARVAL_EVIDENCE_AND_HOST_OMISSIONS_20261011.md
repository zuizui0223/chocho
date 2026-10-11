# Larval evidence changes the interpretation of introduced-only resource maps

Date: 2026-10-11. Status: completed bounded, post-hoc source audit. Not a representative validation sample, new field survey, causal test or corrected global resource estimate.

## Main decision

The investigation has moved from subtle residual spatial effects to observed larvae and explicitly reported host use. Public records provide positive immature-stage evidence in some predicted regions, but two documented native-host omissions also reclassify many candidate records. Geographic agreement alone must not be presented as validation of the exact larval host or introduced-host dependence.

## Fixed 32-species screen

The existing occurrence-validation panel was fixed before querying annotated larvae. The iNaturalist API returned 6,428 records dated 2010-01-01 through 2026-10-11. After excluding 1,062 obscured/private locations, 5,366 research-grade larval records had usable public coordinates. There were 122 records in regions classified as introduced-only by the frozen HOSTS/WCVP reconstruction, spanning six butterflies and 12 butterfly-region combinations. Ninety-two reported coordinate accuracy of 10 km or better; that threshold is not an uncertainty-buffer proof.

Twenty-eight species were retrieved completely under the query. Neptis sappho, Papilio polytes and Pieris brassicae reached the frozen oldest-first cap of 1,000 records each. Pyrgus communis did not pass strict name resolution and remains held, not absent. Counts are descriptive: repeated observations, uneven annotation effort, spatial clustering and capped sampling preclude an occupancy or successful-colonization rate.

## Two native-host omissions materially change the candidate interpretation

Museum Victoria documents Acacia aneura as a host associated with Nacaduba biocellata; Kew lists this plant as native in NSW, NTA, QLD, SOA and WAU. The frozen reconstruction retained only Vachellia farnesiana for that butterfly. Adding the documented native-host link separately adds five native-resource regions and changes its added-region count from 81 to 76.

The Florida Museum lists Passiflora multiflora as a larval host of Dryas iulia; Kew lists it as native in Florida and several Caribbean regions. This link is absent from the frozen reconstruction. Adding the published native range adds one previously unrepresented native-resource region, Florida, and changes Dryas added regions from 94 to 93.

These two targeted additions move 78 of the 122 larval candidates out of the introduced-only category: 68 Dryas records in Florida and ten Nacaduba records across three Australian regions. Forty-four candidate records remain. This is **not** a 64% global error estimate: the checks were targeted and many records share the same regional classification. It also does not prove that none of those larvae uses an introduced plant locally.

The same separate source supplement changes the earlier 32-species occurrence-overlap summary from 66 recovered / 115 outside-native cells to 60 / 109. The six reclassified occurrence cells are Dryas–FLA and Nacaduba–NSW/NTA/QLD/SOA/WAU. Only these two native host-range supplements were applied; other introduced ranges of these hosts were not incorporated. Do not silently substitute this illustrative sensitivity for globally representative correction.

Sources:
- https://collections.museumsvictoria.com.au/species/11366
- https://powo.science.kew.org/taxon/469715-1
- https://www.floridamuseum.ufl.edu/wildflowers/butterfly/julia/
- https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:675237-1

## Directly linked feeding report, with unresolved plant identification

In Kuala Lumpur, Malaysia (WGSRPD MLY), Dryas observations 361335797, 363585138 and 363586114 are linked to the same plant cluster on 14 and 21 May 2026. The observer explicitly reports feeding and identifies the plant as Passiflora suberosa, with a separate plant observation 363586246. That scientific name matches frozen host ID 2560211, which is introduced in MLY, and no native known host is mapped there in the fixed graph.

However, the linked plant observation is still **needs_id**. Record this as **observer-reported larval feeding with a provisional plant ID**, not an independently confirmed exact-host link. Reviewed images show a larva on a damaged leaf of a flowering vine; one photo does not verify completed development. Three linked observations on two dates are one same-plant cluster, not three independent sites or generations.

Two other named-host candidates do not resolve the host species: Papilio polytes observation 17789987 names Citrus at genus rank, and Colobura dirce observation 146491774 names Cecropia at genus rank. They cannot validate an exact WCVP host-species prediction.

## Published breeding confirms some resource opportunities, but reveals a wrong-host match

Phillips et al. (2020), doi:10.1371/journal.pone.0236791, report repeated breeding of Pieris brassicae in New Zealand, with eggs and larvae on hosts including Tropaeolum majus. The frozen graph contains that exact butterfly-host link and treats NZS as introduced-only. This is external published support for a real introduced larval resource, not a newly discovered population. Deliberate eradication is not evidence for natural demographic failure.

Burg et al. (2014), doi:10.1371/journal.pone.0104076, report Dryas iulia larvae feeding on introduced Passiflora foetida in Thailand. Yet the only fixed known host mapped in THA is Passiflora edulis; the Dryas–P. foetida link is absent. Thus even a correct geographic prediction can match the wrong documented host. The same study implicates human transport/release of the butterfly; it is not evidence that plants alone caused the butterfly introduction.

## Verification and preservation

Screen run 38124300207 at cb5f886092e549e2a8e5532946767c3506a8cae5; artifact 11696716253. Linked-source followup run 38124848902 at 635c943ae2545155b77a3eb1fd1b4babc7b8697e; artifact 11696319284. Archive SHA256s were verified.

An independent per-polygon vectorized scan and separately reconstructed plant-resource unions reproduced all 5,366 eligible records' geography and resource categories, with zero disagreements. The independent offline audit also reconstructed the original 66/115 cells and the two-link sensitivity. The followup runner passed all 101 repository tests with the full optional environment; the local suite, including four independent-audit tests, passed 105. The initial screen had 99 passed and one docx-dependent module skipped; the followup resolved that dependency gap.

The machine-readable result is ../results/butterfly_larval_evidence_audit_20261011.json. The downloadable evidence bundle contains the offline verification code, structured metadata, request receipts, fixed inputs and exact observation links. Source photographs and unnecessary user-profile fields are excluded from that reusable bundle.

## Publication boundary

Do not change or merge GEB PR #38 automatically. Its fixed +54.9% arithmetic remains conditional on its frozen graph; this audit is not a full corrected reconstruction. Nevertheless, the newly documented host omissions are scientific review items, not merely author-metadata tasks. The validation wording and source-completeness sensitivity should be reconsidered before asserting that all scientific submission issues are closed.
