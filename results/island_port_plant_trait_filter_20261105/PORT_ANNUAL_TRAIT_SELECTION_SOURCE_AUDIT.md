# Global island port × introduced-plant annual/herb trait selection: empirical gate

**2026-11-05 | Hypothesis test, not manuscript claim.**  
Repo: \`zuizui0223/chocho\` branch \`analysis/island-port-trait-selection-gate-20261105\`  
Primary [WCVP13 × GIFT3.2 126-island workflow](https://github.com/zuizui0223/chocho/actions/runs/41946522894) — **SUCCESS**.  
Independent [offline coverage, omission and multiple-testing audit](https://github.com/zuizui0223/chocho/actions/runs/41947241042) — **SUCCESS**.

## Ecological question

**Do port-rich island regions contain introduced plant assemblages selectively enriched in annual or annual-herb life-history traits, beyond simply containing more introduced flora?**

A positive coefficient for the absolute introduced-plant trait share is necessary for the particular proposed second-stage mechanism. A *relative* introduced-minus-native trait difference is a separate outcome, and a change in native flora does not establish selection of the introduced plants. Neither annuality nor herbaceous growth form is equivalent to nutrient-site preference, leaf nutrient concentration or ruderal abundance.

## Frozen sources and trait construction

- WCVP v13 accepted names / native and introduced WGSRPD3 distributions, **exact pinned \`matildabrown/rWCVPdata\` commit \`65bed76bae9d644ccb6ad200c05f9f5071d89e05\`**. Both original RData files verified against git blob checksums. The accepted species table contains **360,064 species**, and the extracted set contains **224,021 accepted plant × one of 126 strict island botanical regions**.
- GIFT API stable \`index3.2.php\` officially accessible. **32 signed-by-content source response pages** collected for its entire species list, Lifecycle_1 (trait \`2.1.1\`), and Growth_form_1 (trait \`1.2.1\`). Source contents are SHA256-recorded in the primary artifact, and original raw trait data need not be replicated in Git.
- GIFT taxonomy harmonized species name matched **exactly** (canonical genus + binomial epithet) to the WCVP accepted taxon name. No fuzzy matching, genus imputation, or downgrading of missing values to “non-annual”.
- Trait acceptance: **agreement ≥0.66**, non-biased-source and derivation filtering \`biasref=0&biasderiv=0\`. Annual excludes biennial, perennial, mixed/variable or unknown. “Annual herb” requires both unambiguous annual lifecycle and \`Growth_form_1 == herb\`.
- Original WCVP region-specific *entire flora* denominators reproduced exactly against source-frozen island counts (not merely trait-matched subsample). A first implementation failed this integrity test because it dropped taxa whose names could not be reduced to exact binomials; it was fixed by keeping those species as trait-unknown, and the rerun passed.
- All covariates: **log₂(1 + 2023 indexed WPI port count), ln WGSRPD3 geographic area, ln(1+GBIF all-plant indexed records), ln(1+total introduced flora), ln(1+native flora), GIFT trait-coverage proportions in both groups, and Level 1 botanical geography.** 1,999 descriptive permutations shuffle port exposure only within Level 1 × within-Level 1 area-half strata. Each botanical region, not each plant/port, is the independent observational unit.

## Main result (predeclared minimum 20 known species and ≥20% trait coverage in both flora groups)

| Response vs indexed port richness | Island regions | Fitted β per log₂(1 + port count) | Exploratory two-sided permutation p |
|---|---:|---:|---:|
| **Annual fraction among introduced plants with known lifecycle** | 90 | **−0.00939** | **0.4085** |
| Annual fraction among native plants with known lifecycle | 90 | −0.02530 | 0.0285 |
| Annual introduced − annual native fraction | 90 | +0.01592 | 0.0645 |
| **Annual-herb fraction among introduced plants with both traits known** | 87 | **−0.01085** | **0.4595** |
| Annual-herb fraction among native plants with both traits known | 87 | −0.03490 | 0.0085 |
| Annual-herb introduced − annual-herb native fraction | 87 | +0.02405 | 0.0190 |

**Main decision: the absolute introduced-flora annual or annual-herb enrichment prediction is not supported.** The positive relative annual-herb difference is largely due to the *native* denominator declining in more port-rich regions. It does not demonstrate that port-adjacent land filters alien species into an annual-herb community. None of these models measures nitrogen deposition, wastewater application, leaf N, ruderal habitat or increased plant abundance.

## Multiple comparisons and geographic-coverage audit

The six exploratory nominal p-values were BH-adjusted **only as a descriptive multiplicity screen**, not a preregistered hypothesis test. Corresponding *q* values are 0.4595 (introduced annual), 0.057 (native annual), 0.0968 (relative annual), 0.4595 (introduced annual herb), **0.051** (native annual herb), and **0.057** (relative annual herb). None is below 0.05.

**GIFT matched-trait coverage is unequally distributed by origin.** Within the 98 islands with GBIF botanical record measures, median lifecycle coverage was **86.4% in introduced** vs **60.1% in native** plant taxa; median coverage requiring both lifecycle and growth form was **78.8% introduced** vs **56.3% native**. Native coverage had Spearman ρ≈**−0.32** with botanical-indexing effort in this sample, so source trait-known subsets are not representative random samples.

Under pre-specified source-flora set size and post-hoc coverage restrictions:

| Requirement | Eligible island regions | β introduced annual herb | β introduced − native annual herb | Descriptive p for relative |
|---|---:|---:|---:|---:|
| ≥20% trait coverage, ≥20 trait-known taxa | 87 | −0.01085 | +0.02405 | 0.015 (999 resamplings) |
| ≥30%, ≥20 taxa | 83 | −0.00413 | +0.02854 | 0.020 |
| ≥40%, ≥20 taxa | 66 | +0.00147 | +0.03658 | 0.004 |
| ≥20%, ≥40 taxa | 77 | −0.00592 | +0.02216 | 0.056 |
| ≥30%, ≥40 taxa | 73 | −0.00098 | +0.02367 | 0.102 |
| **≥50%, ≥30 taxa** | **45** | **+0.00254** | **−0.00325** | **0.609** |

The introduced share itself stays near zero across quality restrictions; the relative effect changes and is **not robust in the most demanding trait coverage subset**. The 50% subset also has a different biogeographical composition, so disappearance is a sensitivity limit rather than proof that the true relative effect is zero.

Leave-one-Level1-out sensitivity gave positive relative annual-herb β in all nine variants, 0.017–0.036, but this does not resolve the unequal trait coverage or convert the native-denominator contrast into absolute introduced annual-herb enrichment.

## Integration with prior positive island/port results

**Retain without down-grading:**
- Island latent larval host-resource opportunities +118.9% vs mainland +40.9%, source WCVP×HOSTS 239 butterflies.
- Modern port counts positively associated with **all introduced flora** (β+0.2445, botanical effort controlled, n=98) and **distinct introduced HOSTS-linked plant taxa** (β+0.1639, n=98).
- The positive port–host relation does not show a host-specific positive residual after all introduced flora are accounted for (β−0.0415, source descriptive p=0.093).
- Candidate port larval-host flora composition and butterfly assemblage similarity correspond across island pairs (ρ=+0.300; 51 landmasses, 1,275 dependent pairs).

**Do not promote the proposed mediating step:** The available complete global native/alien flora × GIFT records do **not** support a positive relationship of indexed ports with the percentage of annual or annual-herb species among *introduced plants*. “Port → greater absolute alien flora, a subset of which are larval hosts” remains a defensible cross-sectional interpretation. “Port → preferential annual/herb ruderal hosts → disturbance-loving butterflies” is **not** a demonstrated four-stage causal mechanism.

**Separate nutrient evidence:** Previously observed species-level Ellenberg N association (partial rank ρ≈+0.210) reflects European nutrient-site affinity and global port-recording frequency of a restricted candidate host set, not a worldwide nitrogen input or local trait enrichment. The national NIES and 3km modelled human-excretion versus cropland N tests in Draft PR #122 likewise do not identify wastewater nitrogen reaching host plants.

## Next identifiable ecological test — conditional gate only

A true ruderal/nitrophilous selection pathway requires an **independent species-level nutrient-site-affinity or field-measured ruderal performance trait** distributed beyond Europe, and an island/port-local native/introduced plant community observation panel with comparable detection. WCVP annuality alone is insufficient. It is appropriate to retain the paper's already supported island/port botanical-assembly story rather than assume a missing plant-trait mediator or repeatedly tune samples until a nominal positive result appears.

**Data handling:** official GIFT source APIs use externally licensed trait records and require source citations; the GitHub workflow preserves source SHA and derived WCVP accepted-species crosswalk as a 30-day artifact, not as invented universally complete WCVP life-form columns.
