# Corrected host resources versus frozen butterfly occurrence and climate — v06

Date: 2026-10-11. Completed retrospective sensitivity, not a new independent observation study. The submission at 88209b8a215ac1bc9c4c5e67622fa104179c8895 is not changed.

## Result

The association between unrecorded butterfly resource regions and climate mismatch survives the hybrid-marker and BCE/Clarke corrections. The hypothesis that broader host-family diets weaken that separation remains unsupported. Host resources, climatic compatibility and realized breeding must remain distinct levels of evidence.

The occurrence comparison fixes the original 32-species panel and its 1,081 recorded butterfly-by-WGSRPD3 units (48,728 mapped records). Climate analysis fixes the original 24 quality-qualified species, not the 239-species global resource panel.

| Host source | Occurrences recovered beyond native resource geography | Median climate-separation score | Partial host-family association | One-sided residual-permutation p |
| --- | ---: | ---: | ---: | ---: |
| Literal HOSTS | 66/115 (57.39%) | 0.8014 | -0.1664 | 0.2237 |
| Marker reconciled | 66/115 (57.39%) | 0.7903 | -0.1680 | 0.2210 |
| Literal plus BCE/Clarke | 64/113 (56.64%) | 0.8076 | -0.1059 | 0.3148 |
| Marker reconciled plus BCE/Clarke | 64/113 (56.64%) | 0.8014 | -0.1124 | 0.3014 |

The climate score is the cross-pair probability that a never-recorded resource region has greater mismatch to the training-occurrence climate than a held-out observed resource region, with half credit for ties. It is not occupancy probability, breeding success or the fraction of regions made unsuitable by climate. In the combined correction 23/24 species have scores above 0.5. The nonsignificant diet-breadth result does not establish equivalence or a zero effect with only 24 species.

## Occurrence classification and new resource opportunity

Neptis sappho/GRC and Pieris brassicae/AZO move from introduced-only to native-host-covered. Both the recovered numerator and outside-native denominator decrease by two; the 49 occurrence cells outside all documented host geography do not change.

The combined correction adds 160 contemporary resource cells within the fixed 32-species panel. None overlaps the frozen butterfly occurrence sample; 77 have at least ten records of other panel butterflies and meet the original effort-proxy threshold. This is no additional occurrence confirmation, not an estimate that butterflies are absent from all 160 regions. The other 83 cells do not meet that proxy criterion. BCE/Clarke is a geographically biased source-completeness sensitivity, not a globally exhaustive correction.

The corrected continent-composition-matched uniform-placement null expects 48.0814 recovered cells versus 64 observed, with central 95% null range 42–54 and Monte Carlo upper-tail p=0.000005 (199,999 draws). Exact convolution gives approximately 5.43e-7 under the same independence and uniform-placement assumptions. Those nulls do not adequately condition on fine climate, recording effort, shared hosts or dispersal; they do not demonstrate local use or causation.

## Distance checks

The original distance-matching script was reused without scientific changes, with corrected resource breadth in the covariate metadata. In the combined correction the score medians are 0.7143 at 250 km, 0.75 at 500 km and 0.75 at 1,000 km; 21, 23 and 22 species respectively exceed 0.5. At 500 km the partial diet-breadth association is positive (+0.3407), not the hypothesized negative direction. This remains a post-hoc sensitivity, not evidence for a positive causal diet effect.

Distance is calculated between regional representative points and the nearest training-region representative point. Matching the distance-to-training-range does not establish local accessibility or exact matching of observed and unrecorded locations.

## Score provenance

The set of observed regions inside the contemporary host envelope is identical across all four variants. The original training/evaluation partition is therefore unchanged, and no new niche fit is required.

For previously unscored new resource cells, the original specified squared mismatch is a diagonal quadratic in four climate variables, with nine coefficients. Coefficients were algebraically recovered from even rows of saved regional climate/score pairs, required to have full rank, and checked against withheld odd rows; the maximum held-out error was 1.64e-14. Original saved scores were retained verbatim, and the recovered evaluator was only used for new eligible cells. Occurrence presence/absence labels were not used in that coefficient recovery.

An independent implementation recovered the same evaluator from the opposite rows and reproduced the new scores. This verifies an exact evaluator rematerialization from archived outputs, not a new raster extraction or independent verification of the original occurrence-level climate pixels. A separately labelled CPI exclusion check left the primary results unchanged.

## Verification and remaining scope

The parent v05 archive passed all 431 checksums. The old occurrence counts, old 24 climate scores and old primary p=0.2237 were reproduced. The old 24-row distance result table matched exactly; aggregate partial correlations matched within 5.2e-11.

The full updated local code passed 169 tests. Independent arrays, source joins, Mann–Whitney counts, combinatorial distributions, rank projections and great-circle matching checked 47,232 species-region cells, 896 species-level occurrence metrics, 96 climate scores, 79,992 conditional permutation statistics, eight exact overlap null distributions and 288 species-by-distance-caliper scores. Independent code is not independent biological evidence.

New code, fixed inputs, full traces, independent verification and an incremental patch are in the v06 downloadable package. Only this note and the result receipt are added remotely. Plant-level seasonal/geographic nulls on corrected host inputs, integrated manuscript/figure updates, anonymous bundles and cumulative code integration remain outside the completed scope.

Receipt: ../results/butterfly_occurrence_climate_corrected_v06_20261011.json.
