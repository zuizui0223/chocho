# After initial indexing: later butterfly re-indexing and candidate introduced hosts

**14 November 2026 | Executed source-locked follow-up in Draft PR #159**

## Decision

The re-index endpoint does **not establish that introduced larval hosts at ports increase butterfly establishment, population persistence or onward dispersal**. The small descriptive positive contrast has very wide island-bootstrap uncertainty; it changes sign under a sensitivity that requires more than one early butterfly record.

Previous positive **global resource expansion, island amplification, port-flora correlations, host-informed recording prediction and retrospective Circuitscape community gains remain intact**, and were not the endpoints re-tested here.

## Source and outcome

- The 225-butterfly × 51-landmass GBIF panel covers 2001–2010, 2011–2017 and 2018–2022.
- Candidate larval hosts were indexed in 2001–2010 WPI 5-km port buffers. Modern **regional** WCVP native/introduced status was assigned to those historical plant records. Local introduction dates, abundance and actual feeding were not observed.
- The new response is a **2018–2022 GBIF index** for the same butterfly and landmass after its first 2011–2017 index. This is NOT breeding or a self-sustaining population.
- Main botanical quality threshold: ≥10 indexed 2001–2010 plant records on the destination landmass. A missing plant record is UNKNOWN, never genuine plant absence.
- Optional source-condition sensitivity requires a conspecific indexed on another original island before 2011; it is not a proven source population.

## Reproduced results

| Sample | With indexed introduced candidate host | Without indexed introduced candidate host (true presence unknown) | Re-index difference | Island bootstrap 95% |
|---|---:|---:|---:|---:|
| **Main historical port plants ≥10 records** | **11/15 (73.3%)** | **100/144 (69.4%)** | **+3.89 pp** | **−23.34 to +29.84 pp** |
| Conspecific on another indexed island before 2011 | 10/13 (76.9%) | 94/127 (74.0%) | +2.91 pp | −28.53 to +28.43 pp |
| At least 2 early butterfly records | 5/8 (62.5%) | 88/115 (76.5%) | **−14.02 pp** | −56.63 to +21.57 pp |
| Native host also historically indexed | 4/7 (57.1%) | 13/24 (54.2%) | +2.98 pp | −30.98 to +35.71 pp |

These are **post-hoc sensitivity strata, not independent replications**. The original 159 first-index combinations consist of 15 introduced-host-indexed combinations (11 island IDs, 12 butterflies) and 144 without such a record. In the main panel, early butterfly observations average **4.60** for the former, **8.78** for the latter. Unequal early and later observation opportunity can affect re-indexing.

### Four historical host-record types (main n=159)

| Historical port index | First-index combinations | Later re-indexed |
|---|---:|---:|
| Candidate introduced host indexed, no native compatible host indexed | 8 | 7 |
| Both introduced and native compatible hosts indexed | 7 | 4 |
| Native compatible host indexed, no introduced compatible host indexed | 24 | 13 |
| Neither compatible host indexed; genuine presence **unknown** | 120 | 87 |

No entry in this table establishes that an introduced plant was the first usable larval host on the island.

## Source-check and executable reproduction

**GitHub Actions run #43578760844 succeeded**, including both first-index and later-reindex audit. It asserts the exact 159, 15, 11, 144, 100 main denominators and archives all source hashes, model outputs, subgroups, and the complete list of 15 introduced-host-indexed species–island cases.

- [Successful workflow](https://github.com/zuizui0223/chocho/actions/runs/43578760844)
- [Combined two-analysis JSON artifact](https://github.com/zuizui0223/chocho/actions/runs/43578760844/artifacts/14107612832)
- Reproduction: python scripts/audit_port_host_later_reindex_20261114.py --out-dir results/port_host_reindex_20261114

The audit uses 4,999 fixed-seed **whole-island bootstrap resamples** and 2-sided Fisher exact tests as unclustered diagnostics; it is not confirmatory. It is distinct from the previously executed species/island fixed-effect first-index probability model.

## Remaining evidence needed for ecological mechanism

1. Island-local named introduced host plant dated **before** credible butterfly first arrival; documented larval utilization and alternative native resource availability.
2. Repeated breeding or immature life stages across years, with independent observation-effort standardization.
3. Measured butterfly long-distance movement propensity; do not equate wing-size proxies with actual dispersal capacity.
4. Directional flight-relevant **wind vectors** or plant-propagule-relevant **ocean current vectors**, not undirected effective resistance described as a current/trajectory.
5. Source population and actual onward settlement evidence before describing stepping stones.

**Manuscript position:** Keep the GEB host-resource paper intact. The association to an establishment mechanism remains a testable hypothesis, not a supported causal result.