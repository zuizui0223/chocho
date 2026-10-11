# S9 after marker reconciliation: robust totals, changing cell identities

This completed v05 analysis crosses two choices: the conservative HOSTS hybrid-marker reconciliation and the original fixed BCE/Clarke S9 candidate-host augmentation. It is a retrospective source sensitivity, not new field observations or a globally representative correction.

## Same 239 butterflies, four networks

| Input scenario | Native butterfly-region cells | Contemporary cells | Introduced-only added cells | Added/native | Butterflies with expansion |
|---|---:|---:|---:|---:|---:|
| Original literal HOSTS | 26,530 | 41,083 | 14,553 | 54.8549% | 206 |
| Marker-reconciled HOSTS | 26,533 | 41,191 | 14,658 | 55.2444% | 206 |
| Original HOSTS plus BCE/Clarke | 29,852 | 44,819 | 14,967 | 50.1373% | 213 |
| Marker-reconciled HOSTS plus BCE/Clarke | 29,855 | 44,927 | 15,072 | 50.4840% | 213 |

The 1,027 extra candidate links involve 81 butterflies and 679 accepted plant IDs. No candidate link is already recovered by the marker correction. The list is held fixed, including the original unresolved-name exclusions; it is not an exhaustive reanalysis of the original Clarke relational data. IDs are not necessarily independent species. Full totals above include all regions; Jaccard and null comparisons retain the same original 355-region domain even though the augmented native envelope includes a 356th region.

## Why a stable aggregate is not a stable map

Relative to marker-reconciled HOSTS, BCE/Clarke changes 1,403 formerly introduced-only cells to native-host-covered cells, across 61 butterflies. This is 9.5716% of the original 14,658 introduced-only cells. At the same time, 1,817 new introduced-only cells appear, for a net increase of 414. A further 1,919 previously unsupported cells gain native hosts directly. Total native coverage increases by 3,322, and contemporary coverage by 3,736. The lower percentage is not a loss of contemporary resource opportunity; the native denominator also grows. These are source-driven reclassifications, not measured colonizations or ecological losses.

## Small conditional homogenization excess remains

Using the unchanged v04 fixed-margin sampler with butterfly totals, region totals, butterfly-by-continent additions and native structural exclusions, the marker+BCE regional Jaccard is 0.48608762. Its null mean is 0.48267321, leaving an excess of 0.00341441; the second seed gives 0.00341660. Marker-only excess is 0.00359422. All eight chains (four networks times two seeds, 499 states each) have upper-tail Monte Carlo p=0.002, the resolution floor. Chain states are autocorrelated: primary marker+BCE lag-1 is 0.515, approximate AR1 ESS 160/499. These are not independent field replicates or a causal effect. Legacy S9 observations were reproduced; the finite sampler schedule follows v04, not the legacy S9 scheduler.

## The previous selected cases remain candidates

Danaus plexippus/AZO, Pieris brassicae/NZS and Papilio demoleus/CUB remain introduced-only in all four reconstructions. This adds zero new field evidence. It checks only the robustness of previously selected cell classifications; actual absence of all native hosts and source-bibliography independence remain unproven.

## Verification and status

The cumulative local package passes 153 tests. Independent pandas/SciPy reconstruction checked 68 observational metrics, 3,824 species host/range sets, all 239 legacy S9 species-level count vectors and all changed-cell tables. Independent NumPy checks covered 3,992 trace rows, 160 null summaries and 24 saved randomized matrices. The four unaugmented traces reproduce v04 within 1e-12.

Pinned source recovery: run 38147019178, artifact 11705192851, SHA256 70e471beb9bfc74b2f61bfaea0efc92b8a13f874bb3fac918361e5e047bfbb31. The full code, inputs, traces, snapshots, audit and incremental patch are in the verified v05 download. Remote receipt/note commits do not merge accumulated local code.

The corrected occurrence/climate response, plant-level seasonal/geographic nulls, all publication figures and anonymous submission artifacts remain unupdated. PR #38 remains at 88209b8a215ac1bc9c4c5e67622fa104179c8895 and has not been modified or merged.

Machine-readable record: ../results/butterfly_marker_bce_s9_crossed_v05_20261011.json
