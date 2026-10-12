# v08: Native-host triage and an occurrence-stage correction

Date: 2026-10-12. Targeted post-hoc source audit, not a new field study, global host correction, or independent biological validation. Frozen manuscript inputs were not modified.

## Findings

The v07b four-Acacia sensitivity was reproduced. After that prior sensitivity the source-stable introduced-only candidate panel contains 60 butterfly-by-WGSRPD3 cells and 1,475 mapped occurrence records. The current targeted audit partitions them as follows:

| Review state | Cells | Occurrence records |
|---|---:|---:|
| Pyrgus communis: no new native alternative from the specific hosts examined | 22 | 736 |
| Dryas iulia / Florida: HOLD for taxonomic-concept reconciliation | 1 | 167 |
| Remaining cells not host-audited in v08 | 37 | 572 |
| Total after the prior Acacia sensitivity | 60 | 1,475 |

Excluding the named HOLD leaves 59 pending cells with 1,308 records, **not 59 validated breeding cells**. No new field breeding confirmation is adopted in v08. The HOLD does not change any regional native/introduced/absent resource classification, and neither the four original global resource maps nor their estimates are changed.

## Pyrgus: a proposed missing-host explanation did not survive inspection

Pyrgus communis has 26 retained host taxon IDs. The [UC Davis field account](https://butterfly.ucdavis.edu/butterfly/pyrgus/communis) documents native Malvella leprosa and Sidalcea usage in western populations. The former and all five retained Sidalcea taxon IDs are already in the fixed data. Their pinned WCVP native ranges overlap none of the 22 eastern candidate cells. Re-adding these six IDs is exactly idempotent.

A separate raw HOSTS entry, Sida spinosa, is omitted by the exact-name resolver because its string has two canonical candidate IDs: 2601633 (Sida spinosa, Accepted) and 2596030 (Sida lancifolia, reached via an Invalid same-spelled record). The unqualified input lacks authorship. v08 does not automatically assign the ambiguous name. Even in a deliberately conditional sensitivity using ID 2601633, none of the 22 candidate cells gains native-host coverage; the plant is introduced in 20 of them. The conditional assignment would affect other global cells but is not adopted as a correction.

This is a bounded negative audit, not proof that native hosts are absent. The UC Davis source also describes lineage/voltinism complexity, so its western ecological statements should not be projected indiscriminately onto every record labelled communis.

## Dryas / Florida: exact spelling does not guarantee the same taxonomic concept

The [Florida Museum local account](https://www.floridamuseum.ufl.edu/wildflowers/flowers/corky-stemmed-passionflower/) lists Dryas iulia as a larval user of the plant called Passiflora suberosa. The fixed Dryas host list already includes P. suberosa ID 2560211, but that ID's pinned WCVP range does not include Florida, either native or introduced.

Current [POWO P. suberosa](https://powo.science.kew.org/taxon/321973-2) and [POWO P. pallida](https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:185104-2) are separate accepted taxa; Florida is listed as native for pallida. The primary revision, Porter-Utley (2014), DOI [10.3897/phytokeys.43.7804](https://doi.org/10.3897/phytokeys.43.7804), likewise distinguishes members of this species complex. Those facts support a taxonomic-concept review, not global synonymization, nor a claim that a particular butterfly specimen fed on pallida.

The 167 Florida occurrence records are therefore held out of any proposed introduced-only breeding validation pending a local host-name/voucher reconciliation. No global suberosa-to-pallida replacement or new host identity was introduced. Current web taxonomy was used to diagnose the conflict, not to overwrite frozen plant ranges.

## Correction to previous use of the word adult

Earlier v07b notes and column names called these records adult occurrences. This is not warranted by the retained source evidence. The archived acquisition script has two occurrence/search request dictionaries (lines 283 and 390), neither with a lifeStage filter. Its exported CSV retains species, gbif_key, latitude, longitude, year, basis_of_record and dataset_key, but no life stage. It contains 53,434 raw export rows; the original mapped panel contains 32 species, 1,081 species-region cells and 48,728 records.

All counts previously described as adult records should instead be read as **occurrence records whose life stage was not retained in this archive**. This does not mean every source record lacks life-stage information at GBIF, nor that no larvae occur in the source. Adult and larval counts cannot be identified from the retained representation. Counts and original ecological calculations are unchanged by this terminology correction. The new queue uses occurrence_records and an explicit not_retained stage label; the old v07b CSV is retained verbatim only as a historical input.

## Verification and deliverable

- The original v06 full suite was re-run: 169 passed in 4.32 seconds.
- Standalone v08: 8 tests passed, including duplicate-host idempotence, retaining ambiguous names, stage provenance, HOLD versus reclassification, and exact panel counts. The tests also pass after ZIP extraction.
- Independent NumPy/cross-tab implementation checked 4,324 source-variant-by-observed-cell classifications, all 64 original candidate ledger rows, 132 western-host-by-eastern-cell native flags and 22 conditional Sida flags. It does not import the main audit engine.
- The re-extracted ZIP passed all 30 file checksums. ZIP SHA256: 018bc80cd1ed487c0b6a90b38e7a014c5eac61437b5f740362a324c7348e77a5.

Complete standalone inputs, source provenance, code, tests, queue and independent audit are in chocho_native_host_triage_verified_v08_20261012.zip. This note does not claim that the standalone code or cumulative earlier local code has been integrated into the remote analysis branch. No changes were made to GEB PR #38.
