# Butterfly host-list continuity versus introduced-area documentation: source-matched null and negative control

**2026-10-27 — completed source-locked, posthoc/exploratory fixed-margin audit.** Draft chocho PR #78 stacked on #77/#76, outside frozen GEB submission/main. Original 239 butterflies only, original exact HOSTS accepted pairs and their `Location` strings, original WCVP accepted taxonomic identifiers and geographical source. This report does **not** infer new larval host adoption or butterfly colonization.

## Exact empirical identity and successful independent replay

[Original-source GitHub Actions #40461539426](https://github.com/zuizui0223/chocho/actions/runs/40461539426): **success**, 12 source-independent regression tests passed; 11 original author-only WCVP/HOSTS constituent SHA256 checks, 140,485 HOSTS rows / 86,826 binomial rows, 239 focal butterfly species and 969 original disjoint country area atoms restored from their fixed source identities. [Full derived numerical result ZIP](https://github.com/zuizui0223/chocho/actions/runs/40461539426/artifacts/12714916855). No copyrighted raw accepted host interaction lists, individual `Location` triples or original geometries uploaded.

Under the original **95% whole-country botanical origin dominance**, from 3,586 unique `butterfly×accepted plant×original HOSTS Location` source triples:

- **156** different butterfly×accepted plant pairs had at least one report in the plant's WCVP **introduced** region, across 62 butterflies;
- **376** accepted butterfly×plant pairs had at least one geographically classifiable WCVP **native** region report, across 122 butterflies;
- **11** original butterfly×plant pairs had documentation of the same accepted association in **both** plant-native and plant-introduced source regions;
- **282** accepted butterfly×plant pairs had at least one report from WCVP-unlisted regions, across 107 butterflies; **21** were also documented in the plant's native source range.

The global accepted host-list universe for these source original 239 butterflies is **2,596 species×accepted plant pairs**. These numbers distinguish **unique pair overlap** from 179 introduced `butterfly×plant×Location` triples, 12 of which were also documented native. No record in a source does NOT mean the butterfly biologically rejects that host.

## First question: is 11 shared documented pairs unusual under simple source margins?

The artificial reference keeps each butterfly's complete fixed accepted host list size `U_s`, its number of host plants with WCVP-native-area reported use `N_s`, and its number of host plants with WCVP-introduced-area reported use `I_s`. It assigns the documentation labels uniformly across the butterfly's documented hosts and analytically convolves per-butterfly exact hypergeom distributions. This is **source-documentation exchangeability**, not a biological or source-geographic randomization test. A stronger conditional reference holds original accepted host **family** fixed too.

| Source pair overlap | Observed overlap | Source-margin null mean (central 95% reference interval) | Host-family-conditioned mean (central 95% reference interval) |
|---|---:|---:|---:|
| WCVP native ∩ introduced report | **11** | **28.25** (20–36) | **26.98** (20–35) |
| WCVP native ∩ WCVP-unlisted report, generic documentation control | **21** | **50.00** (40–60) | **47.20** (38–57) |

The native-introduced observed/source-margin expectation ratio is **0.389**, versus **0.420** for native-unlisted documentation; plant-family-adjusted ratios are **0.408** and **0.445** respectively. Both show **substantial SOURCE documentation segregation** relative to their artificial labels-exchangeable host pool; a similar deficit in the non-introduced negative control makes it unsafe to interpret the 11/156 as special evidence for the ecological origin of new host associations.

The exact artificial geography-conditioning left-tail probability for native-introduced is **4.10 × 10⁻⁶** (uniform within species) / **5.89 × 10⁻⁶** (family conditioned); native-unlisted left tails **1.08 × 10⁻⁹** / **5.57 × 10⁻⁹**. **These are diagnostic tails of a false-exchangeability source model and MUST NOT be reported as ecological p-values, selection, invasion effects, or independent evidence of rewiring.** Plant regional availability, literature reporting, phylogeny, butterfly distribution and host identity violate their independence assumptions.

There are 45 butterfly species with at least one documentable host pair in **both** native and introduced plant regions, 44 source-permutable species in the within-butterfly null and 36 species with permutation freedom even after host-family constraints. That exceeds the frozen >=20 informative-species and >=5 expected-overlap diagnostic gates, but **does not** remove observational confounding.

## Geographic status sensitivity, no numerical fishing

Reusing the independently original 80% whole-country origin-dominance sensitivity (without changing original source labels, WCVP taxonomy or host lists):

| Source paired overlap | Observed | Uniform source conditional mean (95% interval) | Host family conditioned mean (95% interval) |
|---|---:|---:|---:|
| Native ∩ introduced | **12** | 30.03 (22–38) | 28.84 (21–37) |
| Native ∩ WCVP-unlisted | **22** | 63.57 (53–75) | 59.29 (49–70) |

Introduced documented accepted pairs rise **156→168**; 12 remain co-documented native out of 168, and WCVP-unlisted documented pairs rise 282→450. Both control and introduced-source documentation still show strong geographic source partition rather than globally interchangeable host reportage. Their small ratio differences need a paired genus/plant-biased observational model before any biological specificity claim.

## Ecological implications and what is NOT established

**Positive source-evidence result:** HOSTS contains geographically source-qualified host-use reports on introduced plant areas across **62 butterflies**, and the indexed native and introduced host lists overlap in just **11 of 156 source pairs**, less than simple accepted-host-list and family-conditioned source label-exchangeability expects. However the *generic unlisted source* comparison shows qualitatively the same disproportionate source list segregation.

**Rejected conclusion:** this is not 145 novel host-plant adoptions. The absence of a mapped geographically classifiable native-region HOSTS report reflects documented source coverage, potentially taxonomy, geographic sampling, unknown butterfly nativity, and source collection behavior. This test has not observed selection, actual evolution, butterfly demographic invasion, larval completion, the date of adoption, or climate-niche shifts.

**Next falsifiable biological design:** choose independently confirmed *butterfly* native/introduced range pairs and targeted locality-specific host-use/immature-stage reports with effort-matched source native host records, then test if observed introduced-host-use shows extra host identity turnover beyond matched recording and plant-region availability. Existing PR #75 invasive-butterfly historical roster has only five original-panel butterfly species and sparse original imported-range geographic HOSTS annotations, so published feeding/host fitness or a second original interaction database may be necessary. The exact source code and publication-frozen original GEB result remain separate.

## Reproducible references

- [Pre-outcome fixed documentation-overlap protocol](../../provenance/reviewer_defenses/protocols/HOSTS_NATIVE_INTRODUCED_SOURCE_OVERLAP_NULL_20261027.md)
- [Fully independent conditional source null code](../../scripts/analyze_hosts_native_introduced_pair_overlap_null_20261027.py)
- [Synthetic exact-marginal probability controls](../../tests/test_hosts_native_introduced_pair_overlap_null_20261027.py)
- [Source-verified run and original derived numeric result](https://github.com/zuizui0223/chocho/actions/runs/40461539426)
