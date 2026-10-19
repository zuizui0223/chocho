# Melitaea host effect robustness, maternal crossover, and experimental-unit boundary

**Date: 19 October 2026. Status: SOURCE-VERIFIED EXPLORATORY POST-HOC REANALYSIS OF PRIOR PUBLISHED BIOLOGY; NOT A NEW FITNESS OR INVASION EFFECT.** This addendum is isolated to `analysis/cross-butterfly-lifestage-raw-20261019` and does not alter frozen GEB PR #38, `main`, or the alternative paper draft PR #43. The stage outcome is *alive at diapause*, not adult emergence or recruitment.

## 1. Confirmed effect, source units and uncertainty

From Verspagen et al. (2020), CC0 Dryad doi:[10.5061/dryad.kwh70rz2m](https://doi.org/10.5061/dryad.kwh70rz2m), the original data are **1,680 larval fates in 14 maternal families, 18 egg clutches, 4 day temperatures (28,30,32,34°C, with 15°C nights), 2 host species, 15 larvae per mother×temperature×host cell (112 groups)**. The identical source bytes and source MD5 **373acf78433faddd218bc9d9c427fc9a** are enforced by code.

Using the **source `Dead` code** and equal-weighting each mother over 4 temperatures, `Veronica spicata − Plantago lanceolata` in diapause arrival is **+45/840 = +5.357 percentage points**; the exploratory 9,999-draw mother-cluster percentile bootstrap gives **+1.071 to +10.000 pp**. Mother leave-one-out point estimates range **+3.974 to +6.282 pp**. For 34°C versus 28°C, the host-effect interaction is only **+3.333 pp** with mother-cluster bootstrap 95% **−4.286 to +13.810 pp**: **not a resolved sign reversal or stronger heat advantage**. These are prior-art data, newly audited and reaggregated, not a novel butterfly observation.

## 2. Four contradictory survival/time endpoints do not erase the mean effect

Four source records say `Dead=0` (reached diapause alive) but `Time2=4` (where source README requires `Time2=NA` for surviving larvae). They occur at family×temperature×host:
- 6 / 28°C / Veronica
- 44 / 34°C / Veronica
- 13 / 32°C / Plantago
- 13 / 28°C / Veronica

The code holds initial group denominators at 15 and exhaustively evaluates **all 16 possible binary outcomes of only these four ambiguous cases**. The mean measured host advantage can range **+5.000 to +5.476 percentage points** (the originally coded case is +5.357 pp; forcing all four to failures yields +5.119 pp). These are **partial-identification bounds conditioned on only these four recorded ambiguities**; they do not include other unmeasured bias or dependence.

## 3. Maternal crossover findings are model-sensitive, not established biological reversals

In the original code, **7/14** mothers have both positive and negative `Veronica−Plantago` *sample* contrasts across the four temperatures; under 16 alternative classifications the count ranges **7–8/14**. Seven apparent sign switches must not be called seven independently observed gene×environment×host effects. Small per-cell counts, host choice not measured, temperature-group dependence and sample fluctuation remain.

A purely numerical exact **2^14 mother sign-symmetry diagnostic** for the average host contrast is **0.043945** two-sided using event codes, shifting across the 16 classifications to **0.0391–0.0620**. This is NOT an actual randomization-test p-value. A threshold-crossing change from only four ambiguous rows argues against labeling it a robust confirmatory result.

A further simulation fits **maternal-family-specific constant host log odds**, allows family×temperature baseline changes and an overall common host×temperature shift, and generates new *hypothetical* observations from the same 112 group structure. With 9,999 simulated worlds per assumed within-group correlation:

| Assumed dependency among larvae in the same group | Simulated mean mothers with both signs | Simulated fraction with ≥7 mothers showing both signs |
| --- | ---: | ---: |
| Independent Bernoulli (`rho=0`) | **3.752** | **3.67%** |
| Beta-binomial shared group correlation `rho=0.05` | **4.200** | **6.85%** |
| Hypothetical `rho=0.15` | **4.053** | **5.70%** |

**These group correlations are explicitly assumed, not measured.** Parameters are fitted using the same source observations; these tail fractions are **posterior-predictive-style model sensitivity, not confirmatory biological p-values**. Outcomes are near the ceiling (many larvae survive), so tail fractions need not be monotone in `rho`. No legitimate claim that the sample crossovers exceed an accurately measured biological null is made.

## 4. What this genuinely establishes and what it cannot

- Within this published experiment, **a mean difference in reaching diapause associated with assigned host plant is evident**, and that point contrast is not an artifact of four contradictory time records. No claim of adult reproduction.
- The appealing **maternal sign-crossover story is much more fragile** than the overall host mean. Both recording ambiguity and unknown group correlation affect its statistical interpretation. The authors' earlier published work already identified family-specific plasticity, so it is also not a new biological discovery.
- The full `chocho` 239-species native/current potential host envelopes **do not contain local egg allocation, tagged larval survival, host availability, recruitment or movement**. One cannot transport these diapause probabilities to the butterfly's global distribution cells. 
- Differences among the other available butterflies (lab survivor-only development in two *Danaus*, early field survival in `D. erippus`, unusable pupation outcome in `Vanessa cardui`) are **different biological endpoints**, and still cannot be pooled into a 4-species adult-recruitment law.
- Scientific next gate: acquire **multiple independent butterfly × site × season** systems with *locally coavailable native and introduced host alternatives*, egg allocation or experimentally assigned oviposition opportunity, consistently followed individual larvae to adults/offspring, full plot/mother IDs and habitat context. Predeclare one endpoint, include all eligible study outcomes regardless sign, and test whether quality×choice×phenology predicts held-out recruitment better than host-presence-only.

## Verifiability and branch separation

- Source download/checksum restoration: `scripts/download_exact_cross_butterfly_sources_20261019.py`.
- Four-source-study life-stage receipt: `scripts/analyze_cross_butterfly_lifestage_raw_20261019.py`.
- Exact 16 ambiguous-fate scenarios, 14-family sign symmetry and hypothetical beta-binomial dependence: `scripts/analyze_melitaea_ambiguity_and_crossovers_20261019.py`.
- Regression safeguards: `tests/test_melitaea_ambiguity_and_crossovers_20261019.py`.
- Green, byte-verified [GitHub Actions run 39203671158](https://github.com/zuizui0223/chocho/actions/runs/39203671158), which outputs both original-source stage and sensitivity JSON receipts as its separately labeled artifact.

**Author decision: no merge with `main` or PR #38 and no promotion of these results to a globally supported nonmonotonic butterfly recruitment effect.** A real new ecological advance must be demonstrated on independent butterfly cohorts with choice, season, local availability and adult outcome jointly measured.
