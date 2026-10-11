# Hybrid-marker audit: stable aggregate, nonuniform taxon-level omissions

Date: 2026-10-11. Retrospective source/method sensitivity, not a new feeding observation or independent ecological validation. The fixed 239-species panel and submission head 88209b8a215ac1bc9c4c5e67622fa104179c8895 remain unchanged.

## What was tested

Across the pinned HOSTS export (140,485 rows, 11,795 distinct nonempty host names), candidate WCVP species names were checked for a difference consisting only of the Unicode hybrid marker between genus and epithet. Exact matches and ambiguous names were retained. Only a unique canonical ID was assigned. No fuzzy matching, parent formulas, cultivars or intergeneric hybrids were inferred.

The main sensitivity preserves the original pipeline's accepted-ID policy, which can include Artificial Hybrid targets and species synonyms that resolve to Form. A separate restricted sensitivity only admits canonical targets with rank Species and status Accepted. Confounding the marker change with a new target-status/rank restriction would change two decisions at once.

## Global result

The original literal reconstruction was reproduced exactly, including 53,449 name/insect/accepted-ID rows, 148,254 native plant-region edges, 222,384 contemporary edges, and all 239 species' frozen metrics.

The original-policy correction rescues 63 raw names globally. In the focal panel, it adds 23 unique host-ID links involving 13 butterflies and ten plant taxonomic units. These are recovered recorded links, not newly discovered host use. Eight of the 23 links target forms, so taxon ID counts must not be called independent species counts.

| Fixed 239-species panel | Original | Marker reconciled |
|---|---:|---:|
| Native butterfly-region cells | 26,530 | 26,533 |
| Contemporary cells | 41,083 | 41,191 |
| Introduced-only added cells | 14,553 | 14,658 |
| Aggregate added/native | 54.8549% | 55.2444% |
| Butterflies with any added cell | 206 | 206 |

The aggregate difference is +0.3895 percentage points. The stricter Accepted-Species-only sensitivity gives 55.1954%; both round to 55.2%, but they are different reconstructions and must remain distinguishable.

## Changes are concentrated in Papilio

All 107 genuinely new introduced-only cells belong to five Papilio: P. demoleus +16 (83 to 99), P. machaon +6 (91 to 97), P. memnon +37 (56 to 93), P. polytes +23 (75 to 98), and P. zelicaon +25 (83 to 108). The recovered resources causing these additions are Citrus taxa. It is not valid to interpret these data-reconciliation increments as measured colonizations.

Two cells for Phalanta phalantha (NET, RUW) switch from introduced-only to native-host-covered through Populus × canescens. Pyrgus centaureae gains one native cell in ABT. Thus 107 newly added minus two reclassified gives a net gain of 105 added cells. Existing native and contemporary unions never decrease, but an introduced-only classification can disappear.

## Redundancy and rank boundary

The original 2,596 host links already included 30 subspecies, 20 varieties and two forms. Treating every accepted ID as a separate species can inflate resource redundancy. Under ID counting, single-host share changes from 58.9157% to 58.2208%. Grouping exact canonical species prefixes only within already-supported cells changes that comparison to 58.9638% versus 58.7051%. This grouping does not transfer feeding to unrecorded forms or import a parent species' other distribution records.

For P. demoleus/CUB, the single-lemon correction discussed previously would give seven taxonomic IDs. The full-marker correction gives ten IDs (including two forms) while CUB remains an already-added resource cell. This is not ten independent species and does not retrospectively identify the field collection tree as lemon.

## Verification and decision

The downloadable v03 code package passed 135 local tests. An independent pandas name join and SciPy integer sparse-matrix reconstruction reproduced 63 rescued names, 452 additional all-HOSTS insect-host ID links, all four fixed-panel before/after resource totals, and every one of 110 changed species-region classifications. The original and updated 239 by 367 region matrices contain 87,713 cells per version.

The global qualitative result is retained, while certain species-level maps require correction. The old +54.9% is still a valid description of the frozen literal-name pipeline, not the corrected output. Before replacing submission numbers with +55.2%, regenerate all affected tables/figures and ecological null tests together. No revised ecological null p-values, fitness inference or causal range-expansion claim is supplied here.

The full script, tests, complete name decisions, pair/cell deltas and independent verifier are in the downloadable v03 package and incremental patch. These receipt/note commits do not claim those local files were merged into the remote code branch.

Receipt: ../results/butterfly_hybrid_marker_global_v03_20261011.json. Pinned source acquisition: Actions run 38128030712, artifact 11697114593 (SHA256 ca92cbe98e8e76996c66b61249da0dc8402c31d3d1873e161909547bbe463b51).
