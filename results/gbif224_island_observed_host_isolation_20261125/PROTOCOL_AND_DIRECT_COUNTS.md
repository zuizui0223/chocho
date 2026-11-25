# Independent observed butterfly records in the frozen island host-resource map

## Direct source join completed 25 November 2026 (prior to formal Python CI)
Read both exact gzip blobs directly from the two pinned GitHub source commits; verified gzip CRC and uncompressed sizes, 91,882 unique-source annual rows, 13,587,452 mapped records, 5,390 effort rows, 224 GBIF taxa all within the frozen 239 HOSTS panel, and global matrix totals 26,530 native, 14,553 introduced-only, 41,083 current units. The 126-region and distance-quartile rosters came from the source-verified S3.6 artifact #14997553306 (SHA256 4e665dd376606782a84fbb5e2517cc32cd54abcedf0a0a44f3591d09b001ceb5). The Python workflow separately verifies immutable SHA-256 sources, region-year aggregates, exact original region counts, and direct-count replay. This initial direct join is not an occupancy analysis.

| Cohort / period | Native-host recorded units | Introduced-only mapped host recorded units | No mapped host recorded units | Introduced-only fraction among mapped-host recorded units |
|---|---:|---:|---:|---:|
| 126 islands, 2001–2024 | 916 | 193 | 62 | 17.40% |
| Original 169 eligible global regions, 24 island subset, 2001–2024 | 551 | 69 | 22 | 11.13% |
| Same 24 islands, 2018–2024 | 521 | 65 | 22 | 11.09% |
| Same 24 islands, at least two distinct recording years 2018–2024 | 467 | 53 | 15 | 10.19% |
| Nearest 32 islands, 2001–2024 | 542 | 70 | 38 | 11.44% |
| Farthest 32 islands, 2001–2024 | 17 | 27 | 2 | 61.36% |

Across 126 islands, 69 had an introduced-only recorded unit and 84 had any recorded known-host unit. Potential introduced-only coverage was 5,617 out of 10,312 for the 224 common butterflies (54.47%); only 193/5,617 were indexed (3.44%), compared to 916/4,695 native-host units (19.51%). Thus potential gains are not realized distributions, and recording fractions are not occupancy or establishment rates. In the far quartile only 12/32 islands had recorded known-host butterflies versus 32/32 near. Under the original 169-region eligibility gate, 14 near islands but only one far island remain; the contrast is therefore descriptive and geographically under-supported. The planned models use explicit effort covariates and no significance-selected thresholds, keep insufficient-cohort status, repeat across recurrence and without monarch, and compare observed share with the potential share.

## Independent monarch bridge
Danaus plexippus: 22 native-host and 42 additional introduced-only botanical island regions (64 contemporary); 20/42 additional regions with 2001–2024 records, 18/42 with >=2 distinct recording years. These are not confirmed new breeding populations or dated before/after host changes. The 20 include Hawaii and Fiji. Original historical paper: Scudder (1875), Psyche 1, 81–84, DOI 10.1155/1875/85854, letter from L. H. Gulick reports Hawaiian monarch occurrence following milkweed introduction and describes subsequent Ponape occurrence. This is historical testimony, not an experiment or proof of independent flight. Hemstrom et al. (2022), Molecular Ecology 31, 4544–4557, DOI 10.1111/mec.16592, supports stepwise Pacific expansions including Hawaii to Marianas; genetic timing remains uncertain and transport mode unproven.

## Frozen analysis policy
No observed absence, true larval host selection, introduction causation or resident establishment inference. Contemporary host lists incomplete, record dates are not arrival dates. Do not equate island and botanical region. Report all adverse and sparse-cohort outcomes. Keep this new exploration outside PR179/192 manuscript until replay and model review. No raw GBIF occurrence redistribution or review-package authorization.
