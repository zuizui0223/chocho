# What the open AAI experiment actually rules out in the Kyoto butterfly case

> **Follow-up source-format correction (2026-10-10):** All seven original paired joint losses are in source groups with simple `l.id` records; 11 other source groups have comma-composite IDs. Full-cohort baseline initial leaf area is smaller in joint-loss groups (6.09 versus 9.17; source-strain conditional p=0.0057), but within the 19 *single-ID source groups* it is 6.09 versus 6.82 and the same source-strain conditional test is **p=0.2207**. This is an exploratory, posthoc necessary confounding check; **do not claim leaf size caused joint loss**, nor use the whole-cohort p-value alone. See [final stop decision](KYOTO_AAI_PRETREATMENT_ATTRITION_FINAL_STOP_20261010.md). The overlap signal is a paired observation pattern, **not** an identified plant-chemical mechanism.

**2026-10-10 — original-source reanalysis with NEW precision diagnostics, NOT a new biological result.** Branch `exploration/realized-host-window-v01`; GEB submission PR #38 untouched.

## Authenticated original material and author codebook

- Hashimoto & Ohgushi (2023) *Ecology and Evolution* DOI [10.1002/ece3.10164](https://doi.org/10.1002/ece3.10164).
- Complete author data/codebook: [Figshare 23170898](https://doi.org/10.6084/m9.figshare.23170898).
- `Hashimoto_and_Ohgushi_2023_bioassay_larvae.csv` MD5 `982d4931da306a7ff8e8d550cffe5246`; `Hashimoto_and_Ohgushi_2023_bioassay_leafarea.csv` MD5 `ebe8c2404d9a13b73f08c16ecb9ddc22`. Both original MD5s were reverified independently through GitHub Actions.
- Original design **120 assigned third-instar larvae**, split evenly between butterfly species and equally between added aristolochic acid I (`treatment=a`) and solvent control (`treatment=c`); 30 original `l.id` blocks per butterfly species, two treatment records per block (some original leaf IDs themselves use comma-separated multiple leaf segments).
- Original 24-hour growth: (`day1.mass-initial.mass`)/`initial.mass`; feeding: `consumed.leaf.area/initial.mass`. Leaf area is recorded in square mm; the mass unit in original data/metadata must be checked before converting absolute feeding values to mass-based ecology.
- Original authors explicitly label `loss=1` as larvae **haphazardly dead during the experiment**, excluded from their growth/consumption analysis. Death circumstances and whether exclusion is affected by treatment or handling cannot be inferred from that code.

## Original omissions, not a new treatment effect

| Species | Assigned AAI/control | `loss=1` AAI/control | With valid endpoints | Complete AAI-control source-leaf blocks |
| --- | --- | --- | ---: | ---: |
| *Atrophaneura alcinous* | 30/30 | **13/8** | **39/60** | **16/30** |
| *Sericinus montela* | 30/30 | 0/0 | 60/60 | 30/30 |

Important: the 21 lost *Atrophaneura* are not mechanically equivalent to `dead_immature` from the longer cage competition experiment. Their biological versus procedural causes are unobserved; do not infer mortality-toxicity or claim missing at random. The 16 complete-leaf contrast is **selected after treatment** because only leaves with two surviving/tested larvae remain.

**One original coding disagreement:** `h.id=a4`, composite source leaf ID `1.4,16.2`, treatment control `c`; `l.or.r` is `r,l` in larval file and `r,r` in feeding file. Original herbivore ID, leaf ID, species and treatment agree. The audit *flags but does not silently repair* this one half-label difference; merging by matching original `h.id` and concordant source leaf/treatment retains the recorded feeding outcome.

Leaf-area values also contain zero substitution, not just measured positive feeding: the authors' `note=*` means a negative reconstructed consumed area was set to zero and `note=**` means no bite mark so zero. **12 original zero recorded consumption outcomes**: *S. montela* **11** (five starred correction, six no-bite zero), *A. alcinous* **one** (starred correction). These transformations were authored in the published dataset, not applied in our analysis.

## Leaf-block bootstrap reanalysis of the existing 24-hour experiment

Original data were inspected for structure and loss counts, then an explicit [source-and-precision protocol](KYOTO_AAI_SPLIT_LEAF_ROBUSTNESS_PROTOCOL_V01.json) was frozen **before** re-estimating treatment contrasts. All analyses retain the original excluded/damaged records in denominator audits and calculate short-term growth/consumption only where outcomes were measured. Source-leaf resampling (9,999 draws, seed=20261010) prevents treating the two experimental leaf segments as two independent plant replicates.

| Source species | Complete-block AAI−control 24h growth ratio difference | Bootstrap 95% interval | Both-treatment-valid source-leaf groups |
| --- | ---: | ---: | ---: |
| *A. alcinous* | **+0.0744** | **[−0.0778, +0.2364]** | 16 |
| *S. montela* | **−0.0018** | **[−0.1948, +0.1927]** | 30 |

Growth difference is absolute difference in fractional 24h growth, so **+0.0744 corresponds to +7.44 percentage points in relative 24h weight change**, not +7.44% in biomass. Both bootstrap intervals include biologically relevant positive and negative values; this is **not** an equivalence result and not a proof of no effect.

A complementary all-observed individual analysis preserving whole-leaf bootstrap (but not requiring surviving matched halves) finds *A. alcinous* +0.0892 [−0.0800, +0.2600], and *S. montela* −0.0018 [−0.1955, +0.1876]. These are **post-assignment selected** contrasts; neither corrects potentially informative loss.

Original normalized leaf-consumption difference for paired A larvae was +6.64 (source leaf-area / initial wet mass units) with a very broad bootstrap interval [−12.19, +27.83]. For S larvae it was −11.06 [−66.33, +42.03]. Do not interpret these as standardized assimilated dry mass or compare absolute magnitude across developmental contexts.

The original paper already reported a non-significant AAI effect on immediate feeding/growth. This reanalysis adds **complete-leaf uncertainty and dropout transparency**, rather than presenting that published biological result as our discovery. Software and source: [GitHub Action 38038029233](https://github.com/zuizui0223/chocho/actions/runs/38038029233), three tests green, exact source SHA/MD5 and JSON receipt uploaded. Earlier source-only inventory: [Action 38037793737](https://github.com/zuizui0223/chocho/actions/runs/38037793737), two tests green.

## Independent source changes the biological interpretation

[Nishida & Fukami (1989), *Journal of Chemical Ecology*, DOI 10.1007/BF01014731](https://doi.org/10.1007/BF01014731) purified **seven aristolochic-acid analogues** and showed **synergistic feeding stimulation** in *A. alcinous* with water-soluble leaf extract constituents. This directly undercuts the inference that a one-molecule AAI painting assay could adjudicate all leaf-quality or inducible-chemistry effects. The 2023 assay is a **24-hour single-instar, single-compound addition test**; it neither measures chemistry of leaves after *S. montela* prior herbivory nor tracks recruitment or adult eclosion.

The [2017 host-regrowth paper](https://doi.org/10.1007/s10144-016-0568-8) found that **after regrowth**, previously herbivore-damaged *Aristolochia* foliage did not alter measured subsequent larval growth. That is a different temporal state than tissue examined shortly after feeding. It also does not supply the missing immediate post-damage chemistry or instar-specific food-access measurements in the 2023 competitive cage experiment.

**Narrow, publishable novel ecological mechanism still lacking evidence:** When *S. montela* suppresses native *A. alcinous* without lowering terminal resource inventory, test a **randomized previous-herbivore conditioning × mechanical damage** intervention with donors removed, *and* an independently randomized stage-specific food-access clamp. A plant chemical link is credible only after measured induced metabolites and mechanism-specific manipulation. The existing AAI test does NOT eliminate the plant-chemistry hypothesis; it does constrain claims to what was assayed.

**Stop rule:** Do not reparameterize AAI lags, subset deaths, leaf pair identities or growth endpoints post hoc to find a significant result. Do not append these exploratory diagnostics to the independent GEB manuscript about potential host-plant redistribution.
