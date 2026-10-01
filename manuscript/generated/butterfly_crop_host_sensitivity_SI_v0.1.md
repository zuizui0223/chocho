# Supporting Information — Crop-host sensitivity

WCVP introduced distributions represent wild non-native distributions rather than cultivated acreage. As a sensitivity analysis, we additionally removed host taxa classified as crops by the FAO Indicative Crop Classification before reconstructing butterfly-level native and contemporary resource envelopes. The primary crop definition used exact botanical binomial matches; an intentionally over-conservative variant also removed all congeners when the FAO list specified a crop at genus level (spp.). Butterfly host-family breadth was kept fixed as the specialization predictor.

The FAO parser identified 179 exact crop binomials and 40 genus-level crop entries. In the butterfly host reconstruction, the exact rule excluded 142 resolved host taxa affecting 89 butterfly species; the conservative rule excluded 513 host taxa affecting 117 species.

## Table S1. Sensitivity of resource expansion and host-contribution architecture to crop-host exclusion

**A. Resource expansion in the original 239-species resource-eligible panel**

Analysis | Evaluable species | Expanded species | Native species-units | Contemporary species-units | Added species-units | Aggregate increase | Spearman rho (host-family breadth vs log expansion)
--- | ---: | ---: | ---: | ---: | ---: | ---: | ---:
Primary reconstruction | 239 | 206 (86.2%) | 26,530 | 41,083 | 14,553 | 54.9% | 0.008
FAO exact-binomial crop exclusion | 238 | 203 (85.3%) | 25,722 | 37,810 | 12,088 | 47.0% | 0.052
FAO exact + genus-level conservative exclusion | 233 | 197 (84.5%) | 24,386 | 35,614 | 11,228 | 46.0% | 0.057

**B. Host-contribution architecture in the host-taxonomy-adequate subset**

Analysis | Expanded species | Host taxa for 50% of added credit | Median maximum single-host share | Median effective contributors | Spearman rho (breadth vs effective contributors) | Spearman rho (breadth vs maximum share)
--- | ---: | ---: | ---: | ---: | ---: | ---: | ---:
Primary reconstruction | 191 | 37 | 0.546 | 2.421 | 0.492 | -0.486
FAO exact-binomial crop exclusion | 189 | 38 | 0.619 | 2.000 | 0.444 | -0.422
FAO exact + genus-level conservative exclusion | 184 | 37 | 0.659 | 1.890 | 0.422 | -0.408

Across the full resource-eligible panel, the number of host taxa required to account for half of all added fractional opportunity was 38 in the primary reconstruction, 38 after exact-binomial crop exclusion, and 37 under the conservative genus-expanded exclusion.

Excluding FAO-defined crop hosts did not change the qualitative conclusions: resource expansion remained widespread, host-family breadth remained only weakly associated with proportional expansion, and broader diets retained a more distributed host-contribution architecture.

The conservative genus-expanded analysis deliberately removes non-crop congeners in genera represented as spp. in the FAO source and is therefore treated as a stress test rather than an alternative primary definition.
