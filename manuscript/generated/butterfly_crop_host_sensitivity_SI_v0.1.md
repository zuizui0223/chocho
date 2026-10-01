# Crop-host sensitivity analysis receipt

Submission mapping: Supplementary Table S7 in the v0.2 manuscript package.

WCVP introduced distributions represent wild non-native distributions rather than cultivated acreage. The crop-host sensitivity removes FAO-classified crop taxa before reconstructing butterfly-level native and contemporary resource envelopes. The primary crop definition uses exact botanical binomials; an intentionally over-conservative variant also removes all congeners when the FAO list specifies a crop at genus level (spp.). Butterfly host-family breadth remains fixed as the specialization predictor.

The FAO parser identified 179 exact crop binomials and 40 genus-level crop entries. In the butterfly host reconstruction, the exact rule excluded 142 resolved host taxa affecting 89 butterfly species; the conservative rule excluded 513 host taxa affecting 117 species.

Analysis | Evaluable species | Expanded species | Added species-units | Aggregate increase | Spearman rho (host-family breadth vs log expansion) | Host species for 50% of added credit
--- | ---: | ---: | ---: | ---: | ---: | ---:
Primary reconstruction | 239 | 206 (86.2%) | 14,553 | 54.9% | 0.008 | 38
FAO exact-binomial crop exclusion | 238 | 203 (85.3%) | 12,088 | 47.0% | 0.052 | 38
FAO exact + genus-level conservative exclusion | 233 | 197 (84.5%) | 11,228 | 46.0% | 0.057 | 37

Excluding FAO-defined crop hosts did not change the qualitative conclusions: resource expansion remained widespread, host-family breadth remained only weakly associated with proportional expansion, and 37-38 host species still accounted for half of all added opportunity.

The conservative genus-expanded analysis deliberately removes non-crop congeners in genera represented as spp. in the FAO source and is therefore treated as a stress test rather than an alternative primary definition.
