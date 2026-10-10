# AAI split-leaf assay: joint mortality, baseline leaf size and source-ID confounding — final data decision

**Date:** 2026-10-10 | **EXPLORATORY REANALYSIS CLOSED — NO NEW CAUSAL BIOLOGICAL DISCOVERY**  
Branch `exploration/realized-host-window-v01`. The independent Global Ecology and Biogeography paper PR #38 is unchanged.

## Scientifically relevant question

The published experiment by Hashimoto & Ohgushi (2023), *Ecology and Evolution*, DOI [10.1002/ece3.10164](https://doi.org/10.1002/ece3.10164), added aristolochic acid I (AAI) to leaves and reported no detectable 24-hour larval-growth or consumption response. The original [Figshare dataset 23170898](https://doi.org/10.6084/m9.figshare.23170898) actually assigned **60 *Atrophaneura alcinous*** larvae in **30 recorded matched `l.id` source groups** (30 AAI / 30 control); only **39** had measured growth outcomes. Thus the published no-significant-growth-effect is **not equivalence** and is conditioned on postassignment loss.

This source-level reanalysis asks whether paired loss is already related to pretreatment conditions, rather than claiming AAI toxicity or a new plant-induced chemical mechanism.

## Exact data and source quality

- Author original AAI larva table [file 41146985](https://ndownloader.figshare.com/files/41146985), MD5 `982d4931da306a7ff8e8d550cffe5246`.
- Original half-leaf area table [file 41146988](https://ndownloader.figshare.com/files/41146988), MD5 `ebe8c2404d9a13b73f08c16ecb9ddc22`.
- All **120 originally assigned larvae** and leaf-area rows were retained, joined by `h.id` and independently confirmed treatment/`l.id`; one half-side label mismatch for `a4` was preserved and reported rather than repaired without explanation.
- The native butterfly has 13 author-coded `loss=1` in AAI, 8 in control, with 7 groups missing both caterpillars. Exact McNemar paired two-sided p=0.125 for differential loss by AAI assignment. The original publication describes haphazard larval deaths. No independently verified toxicological causes.
- The 30 native source groups comprise **19 simple `l.id` strings and 11 comma-composite `l.id` strings**. A comma may represent several original leaf fragments/IDs; the exact physical provenance is *not* fully reconstructable from this file. These must NOT be assumed identical experimental source units without qualifications. **All seven double-loss groups had single IDs**. Source-ID structure is a necessary, originally overlooked inferential restriction.

## Actual source-verified new analysis

The pretreatment covariate comparison was defined in [protocol v0.1](KYOTO_AAI_BASELINE_LEAF_DEATH_SENSITIVITY_PROTOCOL_V01.json) before inspecting original baseline areas and masses (but after the earlier discovery of seven shared losses), with three predeclared baseline descriptors and 9,999 randomized label permutations preserving authors' four native butterfly `strain` strata.

| Baseline feature, original matched source groups | Dual-loss−other contrast | Strain-restricted two-sided p | Holm-adjusted p |
| --- | ---: | ---: | ---: |
| Mean initial leaf-half area (original numeric source units) | **−3.0807** | **0.0057** | **0.0171** |
| Mean log initial caterpillar mass | +0.0708 | 0.4684 | 0.4684 |
| Fraction of comma-composite source IDs | −0.4783 | 0.0279 | 0.0558 |

**Crucial confounder revealed:** In all 30 source IDs the mean initial leaf-half area was **6.0918 in the seven both-loss groups** versus **9.1725 in the 23 others**. However, 11/23 of the "others" had comma-composite IDs, versus **0/7** of the both-loss groups. The composite ID source format is tied to how original foliage is recorded and to measured starting available area.

Therefore we ran **one explicitly posthoc, necessary source-format sensitivity**: limit to the **19 single-ID source groups**, compare the seven dual-loss versus 12 other groups, and again permute the loss labels within the original butterfly strain groups, 9,999 draws.

| Comparable source-format subset | Double-loss n | Other n | Mean initial area double-loss | Mean initial area other | Difference | Permutation p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| **Single-ID source groups only** | **7** | **12** | **6.0918** | **6.8231** | **−0.7313** | **0.2207** |

**Scientific judgment: baseline leaf area is NOT a robustly established predictor of double-loss independent of source-ID grouping.** The initial full-cohort area association is sensitive to observation-unit composition, and the p=0.0057 should NOT be reported alone or promoted as "small leaves caused mortality."

The previously observed double-loss clustering remains a feature of the recorded paired data. Its source-strain-fixed exact overlap check p=0.00711 tested **paired co-loss under source-strain margins**, not a randomized leaf-size, chemical or co-housed stress manipulation. It cannot assign the common cause of loss.

## What public sources answer and fail to answer

**Source-supported:** A large fraction of the *Atrophaneura* AAI assay was lost after assignment; loss codes are correlated within source-block pairs; source leaf ID format and initial leaf area are confounded; accounting for that confounding removes the apparent robust leaf-size relationship.

**Not source-identified:** AAI lethal dose response, whether AAI altered native butterfly risk of accidental death, causal effects of available leaf mass versus leaf quality, chemically induced leaf effects following *Sericinus montela* feeding, direct butterfly interference, viability to adult eclosion, or general ecological trap effects.

No available 2017 original individual-level *A. debilis* regrowth dataset was verified in this continuation: the [2017 paper](https://doi.org/10.1007/s10144-016-0568-8) and supporting-information PDF describe *already-published* plant compensatory responses; the original raw per-plant regrowth treatment outcomes were not retrieved and must not be silently reconstructed from printed summaries.

## Publication and analysis STOP rule

1. **Stop** additional significance-seeking covariate scans/subsets of this 120-larva AAI dataset. The already planned single-ID source confounding check has been run; the nonrobust outcome should be preserved.
2. For a new field/controlled trial, record actual plant accession **and leaf-pair physical identity** unambiguously, host fresh leaf surface and age, source maternal lineage, experimental handling, individual larval death reason and all original fates.
3. The independent mechanistic research priority remains randomized *A. debilis* edible access across 3–5th instars with *S. montela* present/absent, followed if needed by separate species-specific prior-feeding versus mechanically damaged plant trials.
4. **No positive new biological mechanism** and no main-GEB paper modification was produced in this exploratory source audit.

## Verification

Full baseline+format and single-ID robustness audit: [GitHub Actions run 38044801095](https://github.com/zuizui0223/chocho/actions/runs/38044801095), **5 tests passed**, original source MD5s matched, and archival JSON receipt uploaded. Source: `scripts/audit_kyoto_aai_pretreatment_joint_loss.py`. Reanalysis tree contains prereviewed protocol and posthoc labeled format sensitivity `KYOTO_AAI_SINGLE_ID_LEAF_AREA_SENSITIVITY_V01.json`.
