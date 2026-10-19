# Butterfly-family and receiver-region heterogeneity of climate-sorted larval-resource opportunity
**2026-10-19 | New source-verified original MACRO ECOLOGY diagnostic | Exploratory observational resource geography, not actual butterfly colonisation.**

## Specific question and held-constant model

Prior original global WCVP/HOSTS/BCE/CHELSA analyses found a tiny but positive climate-alignment residual after the introduced host PLANTS were rearranged under butterfly→host links, plant family, native host geography, continental and climate-bin margins. The new ecological question is whether this residual generalises across butterfly taxa and recipient geography, or is concentrated.

We reused the **unchanged BIO1/BIO12 3×3 continent-local plant swap model** on the **same frozen 239 butterflies × 365 WGSRPD3 cells**, plus exactly 1,027 accepted BCE/Clarke partial host links in 81 butterflies, with 679 additional accepted plant IDs. The augmented resource has 37,869 introduced plant×WGSRPD3 edges; 235 butterflies (14,961 novel resource region cells) qualify by original native-host climate support. A fixed-original-eligibility set of 234 butterflies (14,903 novel cells) is evaluated using the **same null graph draws**. Butterfly family labels come from frozen LepTraits source, six families exactly. Source/archive SHA, families and exact native and introduced resource cells are checked.

**[Source-verified run 39229792077](https://github.com/zuizui0223/chocho/actions/runs/39229792077)** is GREEN after fixing the final CI test's status-string mismatch; the original source-based computation itself had already run correctly. Two independently seeded 199-sample chains passed global source/plant family/3×3 climate/recipient-degree invariants, ≥10% burn movement, ≥100 swaps per spacing, global ESS≥30 and interchain median agreement. **Derived JSON artifact:** [12176772863](https://github.com/zuizui0223/chocho/actions/runs/39229792077/artifacts/12176772863) contains every family and Level1 recipient group for both climate windows and both eligible-species panels.

Our group estimand is **observed minus null median** percentage points in fraction O/A, where O = newly mapped butterfly×region host-resource cells within the original native HOST-PLANT BIO1/BIO12 climate envelope and A = all added resource cells, recomputed in EACH null graph and group. This is not an actual butterfly climate physiology or occurrence analysis. MCMC null 95-percentile ranges characterize algorithmic reallocation, NOT biological uncertainty/confidence intervals for causal effects.

## 1. Butterfly family heterogeneity: primary native-host climate 5–95%

| Butterfly family | Eligible butterflies with added cells | Observed added resource cells | Observed minus plant-network/climate-matched null median | Prespecified status |
|---|---:|---:|---:|---|
| Hesperiidae (セセリチョウ科) | 47 | 3,436 | **+0.824 to +0.845 pp** | QUALIFIED, observed higher than both chains' 97.5% null quantiles |
| Nymphalidae (タテハチョウ科) | 89 | 5,806 | **+0.512 to +0.523 pp** | QUALIFIED, exceeds both 97.5% null quantiles |
| Papilionidae (アゲハチョウ科) | 14 | 735 | **+0.999 to +1.043 pp** | QUALIFIED, exceeds both 97.5% null quantiles |
| Lycaenidae (シジミチョウ科) | 25 | 2,169 | −0.195 to −0.196 pp | QUALIFIED numeric estimation, but **observed WITHIN both null 95% ranges**, no supported negative deviation |
| Pieridae (シロチョウ科) | 36 | 2,763 | −0.105 to −0.160 pp | QUALIFIED numeric estimation, **WITHIN null 95% range**, no supported negative deviation |
| Riodinidae (シジミタテハ科) | 1 | 52 | NOT EVALUABLE | Prespecified low-data flag (fewer than five species and 150 resource cells) |

Original total: **235 climate-eligible butterflies**, 14,961 added resource cells. This study can assess broad patterns in only FIVE adequate families; one Riodinidae species must not be portrayed as a sixth independent confirmation or failure.

At the tighter 10–90% analogue, numeric residuals for all five adequate families turn positive: Hesperiidae **+1.198–1.221**, Nymphalidae **+0.939–0.955**, Papilionidae **+0.893–1.024**, Lycaenidae **+0.566–0.603**, and Pieridae **+0.108–0.183** percentage points. Pieridae remains inside the chain-specific null percentile envelopes. Thus positive across all five *point estimates* at the narrower climate definition does not prove five independent confirmatory effects. The broad-climate signs for Lycaenidae and Pieridae also remain small/nonextreme on the **matched 234-species** cohort (Lycaenidae −0.204–−0.232, Pieridae −0.105–−0.160 pp); the newly eligible *Satyrium pruni* is not the sole cause.

## 2. Receiver biogeographic region heterogeneity

Each introduced-host butterfly×region resource cell is assigned to its **recipient WGS RPD Level1** (not known real butterfly occupancy). At broad 5–95% analogue:

| WGSRPD Level1 region | Observed added host-resource cells | Observed minus null median (pp; two chains) | Readout |
|---|---:|---:|---|
| Northern America (7) | 2,998 | **+1.224 to +1.225** | Largest positive residual; exceeds each null 97.5% quantile |
| Pacific (6) | 1,695 | +0.491 to +0.519 | Positive point value; within null percentile ranges |
| Tropical Asia (4) | 1,116 | +0.354 to +0.395 | Positive point value; within null |
| Africa (2) | 2,665 | +0.324 to +0.388 | Positive point value; within null |
| Southern America (8) | 2,503 | +0.154 to +0.198 | Positive point value; within null |
| Australasia (5) | 1,281 | +0.018 to +0.028 | Very small; within null |
| Europe (1) | 971 | −0.034 to 0.000 | Essentially no residual at this climate definition |
| Temperate Asia (3) | 1,208 | −0.161 to −0.067 | Slightly negative point; within null |
| Antarctic/subantarctic (9) | 524 | 0.000 | Structurally non-informative, NOT confirmed butterfly habitat |

At narrow 10–90% analogue, notable positive residuals are **Southern America +1.672–1.729 pp**, **Africa +1.231–1.290**, **Northern America +0.956–0.960**, and **Tropical Asia +0.612–0.631**; each exceeds the chain-specific 97.5% null percentile. Others are smaller or within percentile intervals. This sharp change in regional pattern with analogue width is a major warning against describing one universal climatic affinity rule.

**Do not infer resource-independent spatial abundance**: all region groups share butterfly taxa and host networks, WGSRPD3 regions have greatly unequal area, the HOSTS/BCE source captures a geographically selected butterfly sample, and no true butterfly occupancy or annual range expansion is measured.

## 3. Ecological mechanism and decision

1. **A non-obvious global average conceals clade/realm heterogeneity.** The small +0.387–0.415 pp global 3×3 residual is not equally positive for every adequately sampled butterfly family and every recipient region. Hesperiidae/Nymphalidae/Papilionidae contribute positive resource alignment under wide climatic analogues; Northern America is the clearest positive recipient region. Weak negative family means are not distinguishable from source-null redistribution variation.
2. **Window sensitivity matters:** at 10–90%, multiple other clades/regions become positive. The ecological conclusion is not that some butterflies adapted evolutionarily to anthropogenic plant introductions, but that plant-resource opportunity climate sorting is geographically and taxonomically structured at this macro spatial scale.
3. **Plant-side geography remains dominant:** earlier plant-native climate and 1,500km source-conditioned nulls reduce the original large +8.28pp climate analogue statistic to only ~+0.34–0.42 pp globally. Families/receivers are not independent randomised replicates of plant invasions.
4. **The candidate journal storyline** should emphasize *uneven, plant-mediated ecological opportunity in a fixed recorded butterfly-host network* and the limits of apparent universality, not butterfly demographic population benefits. This is a new user-directed macro analysis on the same original source data, not a re-analysis of individual published caterpillar experiments.

**Stop/go**: The next defensible biological layer, if pursued, is a large observer-effort- and climate-matched GBIF occurrence comparison in native-only versus added-host resource locations. This has NOT been performed by this result. Do not elevate the source-biased family/realm contrasts to conservation prioritisation or species recovery predictions. Keep exploratory [Draft PR #46](https://github.com/zuizui0223/chocho/pull/46); frozen GEB PR #38 and main unmodified.
