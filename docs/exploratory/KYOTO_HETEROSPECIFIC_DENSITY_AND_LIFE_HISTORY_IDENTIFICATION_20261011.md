# Kyoto蝶の競争は「餌量×密度による生活史分岐」を経由するか：2020年先行実験との統合

**2026-10-11｜先行研究＋元CSV構造の再検証で、新規性の範囲と検証に必要な条件を確定**  
**Status: publication-verified ecological PRIOR ART / newly formalized falsifiable cross-species question, NO new butterfly behavioral or experimental outcome.**  
**Working branch:** `exploration/realized-host-window-v01` only; **GEB PR #38 unchanged**.

## What is already known — and MUST NOT be claimed as newly discovered

**Kozuki & Takeda (2004)** [DOI 10.11257/jjeez.15.157](https://doi.org/10.11257/jjeez.15.157) documented among *Atrophaneura alcinous* Japanese population strains how host-plant conditions, starvation, overcrowding, photoperiod, temperature and isolation contribute to split life cycles and diapause. Their three habitat populations differed in host use and diapause dynamics. **Crowding and food as a cause of diapause is not a new concept.**

**Takagi & Miyashita (2008)** [DOI 10.1603/0013-8746(2008)101[392:HPQIDI]2.0.CO;2](https://academic.oup.com/aesa/article/101/2/392/8474) directly showed that tougher *A. kaempferi* host leaves change native swallowtail early survival, development and pupal diapause under controlled light/temperature. Living pupae must not be labeled dead merely because they do not emerge promptly.

**Nakahara, Horita, Booton & Yamaguchi (2020)** [DOI 10.1111/ens.12399](https://doi.org/10.1111/ens.12399), *Entomological Science* 23:57–65, explicitly crossed **host food quantity** and **CONSPEFIC *A. alcinous* larval density** in four experimental groups, following individual developmental strategies. The authors already found:
- food restriction increases extra molting, larval development time, and reduces larval and pupal body weights;
- overcrowding combined with food shortage increases **conspecific cannibalism and disease**;
- under food shortage, pupal **diapause emerges more under low conspecific density than high density**;
- extra-molting, cannibalism risk and diapause interact as alternative life-history fates, with developmental asynchrony proposed as a cannibalism-risk escape route.

This is **published empirical evidence**, not a new outcome. The published abstract and references are public, but we have **not obtained individual-level original 2020 trial data**; do not invent sample sizes, colony IDs, exact event-time tables, conditional hazard coefficients or causal mediation estimates from the abstract.

**Hashimoto & Ohgushi (2023)** [DOI 10.1002/ece3.10164](https://doi.org/10.1002/ece3.10164) independently manipulated **S. montela and A. alcinous densities** in *A. debilis* cages at **15 combinations ×2 cages =30 cages** and recorded native caterpillar stages, survival, developmental time and final leaf scores. Their published results are themselves novel prior art: increasing introduced *Sericinus* suppresses native *A. alcinous* survival; an **intraspecific×interspecific larval density interaction** affects *A. alcinous* developmental time, especially when native larvae are more numerous. However, that paper did **not** independently manipulate available food quantity and did not report individual extra molting, conspecific cannibalism, disease, pupal diapause or adult flight outcomes.

### Source-verified feasibility: what can be reanalysed from the actual 2023 author CSV?

Data deposited by the authors at [Figshare 10.6084/m9.figshare.23170898](https://doi.org/10.6084/m9.figshare.23170898), source `larvaldata.csv` [file 41146991](https://ndownloader.figshare.com/files/41146991), exact MD5 `0e9ff924e4ae38cb78529ec540887c42`: **690 original cage×date rows, 30 independent cages, 23 repeated surveys**.

Actual native columns: `a1,a2,a3,a4,a5,app,ap,cumul.ap,total.a` (cage-level instar 1–5, prepupae/pupae, cumulative pupation) and the two original *species densities*. Missing: **unique larval identity**; **sixth instar / extra-molt event**; **cannibalism event/who attacked whom**; **disease cause**; **pupal diapause**; **adult emergence and flight-capability**; **synchronous edible leaf supply**. The life-history pathway cannot be resolved by simply fitting more GLMs to these counts. Lower juvenile counts could be mortality, dispersal among plants, pupation or missing detection. The basic survival/development density interaction is already the article's main result, not a novel discovery to claim in a second manuscript.

Original authentic-source replay plus explicit binary visibility audit is in `scripts/audit_kyoto_stage_demand_structure.py`, run by `.github/workflows/kyoto-original-stage-demand-audit.yml`; no original person-level histories or uncollected larval molt counts are fabricated.

## New falsifiable question worth testing (NOT already answered by either paper)

> **Does the identity of a second herbivore species redirect native *A. alcinous* from compensatory extra molting/diapause into direct mortality, relative to equivalently crowded native *A. alcinous* cohorts, and does this change across edible host supply and native cohort density?**

This is *not* the generic already-known claim “starvation/crowding induces diapause,” nor does observing an extra molt alone establish an adaptive evolutionary advantage. Mechanistic alternative predictions:

1. **Shared resource demand / total burden:** *S. montela* and additional native larvae primarily alter food consumption and instar-specific resource access. Accessible host-food replacement should reduce the harm under both species-addition treatments; effects at equal numbers may still differ because an *A. alcinous* larva eats ≈**2.6×** as much in the prior published rearing context.
2. **Species-identity plant-quality/encounter pathways:** At comparable stage-specific actual food access, *S. montela* addition changes *A. alcinous* development, diapause or adult output **differently** from conspecific addition. This would indicate an identity-associated bundled effect, **not uniquely plant chemistry** unless prior herbivore-conditioning/wounding and foliage nutritional chemistry are independently randomized or assayed.
3. **Density-mediated life-history redirection:** The same heterospecific addition has a different effect at **low versus high focal native cohort density**, with measurable differences in extra molts, confirmed cannibalism, disease, live diapausing pupae, and eventual flight-capable adults. In low density/poor food, extra molts or diapause may rise; in crowded/poor food, direct mortality/cannibalism may dominate. Actual directional responses are **HYPOTHESES**, not measured outcomes.

## CRITICAL design correction to earlier food-clamp proposal: control total larval abundance

The existing `KYOTO_HIDDEN_RESOURCE_BOTTLENECK_FACTORIAL_V01.json` contains a valid **2×2 effect of adding *Sericinus* with natural vs supported food**. However, by adding larvae it **also changes total larval head count**. It CANNOT separate an **interspecific identity-specific effect** from the effect of adding any larvae. The 2020 density experiment makes this missing comparison biologically important. Do NOT quietly rewrite the previous 2×2 as if it already proves heterospecific specificity.

A minimal **independent source-block randomized extension** to test identity is **2 resource access treatments × 3 additional-larva assignments**:

| Arm | Focal initial *A. alcinous* cohort | Added larvae | Support for edible *A. debilis* |
| --- | --- | --- | --- |
| 1 | Fixed and tagged | none | natural |
| 2 | Fixed and tagged | **additional native A larvae** | natural |
| 3 | Fixed and tagged | **same NUMBER of invasive S larvae** | natural |
| 4 | Fixed and tagged | none | supported |
| 5 | Fixed and tagged | additional native A larvae | supported |
| 6 | Fixed and tagged | same number of invasive S larvae | supported |

The additional-native arm and introduced-species arm are matched in **head count**, and the original **focal** tagged cohort is held fixed as the recipient. Do not count extra native additions as focal recipients in adult emergence denominators. Propose minimally invasive, independently validated individual IDs (camera/photographic identification, genetics or other animal-welfare-compatible approach) *before field use*, and ensure introduced larvae cannot escape or be accidentally released.

**Predeclared estimands** at each food level:
- `C_A = native focal adult proportion(additional A) - proportion(no additional larvae)`;
- `C_S = native focal adult proportion(additional S) - proportion(no additional larvae)`;
- **identity contrast** `I = C_S - C_A = proportion(additional S) - proportion(additional A)`;
- **food supply interaction of identity** `I_supported - I_natural`.

**Interpretational limit:** Match *head counts* does **not match food demand, species size or consumption**. The 2023 article found *A. alcinous* larval consumption about **2.6 times** *S. montela*. To infer direct species-specific chemistry or physical interference (rather than different resource extraction), either calibrate and prespecify a consumption-matched independent subdesign or rigorously measure food access/mass/quality as a **postrandomization diagnostic**, without naively conditioning the primary ITT estimand on observed food use. The food clamp itself can change handling/tissue age, so matched inaccessible sham leaves remain necessary.

If the MAIN aim remains “**conspecific density gates interspecific life-history redirection**,” the necessary independent intervention adds a **low/high starting focal cohort density** crossed with these six treatment types (**2×3×2 =12 combinations**). This increases required resources substantially; the minimal six treatments answer the **identity** question at just one predeclared native baseline density, **not** the density-moderation hypothesis.

There is currently **NO sample size estimate validated from pilot between-cage effect variance**. The historical 30 cages (two per combination) cannot simply be reclassified as replicated randomized 6- or 12-arm resource-clamp treatments; choose allocation via pilot/simulation only after collecting appropriate variance/censoring information. The primary endpoint must separate **adult flight**, **live diapause**, **extra molt confirmed**, **directly observed cannibalism**, **confirmed disease**, **other death** and **unresolved**. Keep mother line, whole-host plant age, leaf node/chemistry, photoperiod, plant accessibility and proper independent cage unit measured/blocked.

## Novelty audit / precise publication claim

Potential future paper question: **“Does heterospecific rather than conspecific crowding alter life-history strategies and eventual emergence of a native specialist butterfly under controlled host resource access?”** This is a meaningful test of how the **IDENTITY of resource competitors changes the DISTRIBUTION of native life-history outcomes**, rather than a tautological claim that more hosts increase feeding potential.

**Cannot be claimed now:** measured heterospecific-driven extra molts, chemically induced diapause, adaptive competition avoidance, subsequent generation abundance, or any adult emergence shift. The 2020 source may already have related density-mediated life-history findings; 2004 has habitat-dependent split cycles. Novelty exists only in the **species identity and source-controlled interaction**, not any individual life-history effect.

**Stop rule:** Do not further mine the same 690 source rows for nonexistent sixth-molt, disease, cannibalism or diapause labels. Do not make the main GEB introduced host-plant paper a different species-level competitive interaction study. **No 2020 raw CSV found or retrieved** in searches of publicly indexed publisher/citation resources; availability may require asking the authors, and absence from search results is not proof of nonexistence.
