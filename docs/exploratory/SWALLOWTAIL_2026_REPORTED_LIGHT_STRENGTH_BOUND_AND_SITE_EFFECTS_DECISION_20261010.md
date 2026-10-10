# What the published light-intensity odds ratios do and do not explain about larval host co-use

**2026-10-10 | New mathematical conditional-sensitivity result; NO newly estimated ecological/competition effect.**  
Branch: `exploration/realized-host-window-v01`. **GEB PR #38 unchanged.**

## Source and original evidence

Jang et al. (2026), *Journal of Ecology and Environment* 50:15, DOI https://doi.org/10.5141/jee.26.017. The authors report microhabitat differentiation between *Sericinus montela* and *Atrophaneura alcinous* on Korean ***Aristolochia contorta***, and a 2×2 set of original HOST-RAMET observations:
- neither: **542**;
- *Sericinus* only: **204**;
- *Atrophaneura* only: **124**;
- both: **5**;
- total **875**, with *Sericinus* marginal **209** and *Atrophaneura* marginal **129**.

The article's reported adjusted odds ratios for relative light intensity, as transcribed in the earlier source audit, are **OR_S=1.009** and **OR_A=0.979** per +1 relative light percentage point (paper's denominators and regression level differ from a simple 875-ramet marginal model). During this continuation the journal web host blocked a fresh retrieval of full-text coefficients; these are **previously transcribed numbers, not independently re-verified against original PDF on this run**. The article's 255/215 discrepancy is in the **quadrat** sample, not this internally consistent 875-ramet contingency table. The site×month×light×ramet raw file is not currently available.

An earlier chocho mathematical stress-test [here](SWALLOWTAIL_2026_MICROSITE_OVERLAP_EXACT_NULL_AND_NO_CAUSAL_INFERENCE_20261010.md) proposed invented 400 'sunny' and 475 'shady' ramet classes in which within-class species occurrence could be independent while still reproducing all four printed pooled counts exactly. **But that demonstration had not calibrated its light responses to the authors' coefficients.**

## A. Quantitative correction: invented two-class sorting required much stronger light slopes

Take that earlier **invented** example's group marginals:

| Arbitrary class | N | Sericinus presence | Atrophaneura presence | Both present (a possible realization) |
| --- | ---: | ---: | ---: | ---: |
| 'Sunny' | 400 | 200 | 5 | 3 |
| 'Shady' | 475 | 9 | 124 | 2 |
| All | 875 | 209 | 129 | 5 |

If one **forces** 'sunny' and 'shady' to be exactly 100 and 0 percentage points relative light, and assumes no other habitat effect, their incidence probabilities imply:

- OR_S per one +1% light: **1.04026**, far steeper than reported **1.009**;
- OR_A per +1% light: **0.96726**, steeper in the negative direction than reported **0.979**.

Hence the old toy **must NOT be presented as quantitatively compatible with the fitted light coefficients in a light-only logistic probability model**, even though it establishes categorical environmental-sorting nonidentifiability from pooled counts.

## B. Rigorous lower bound on EXPECTED overlap within a restrictive light-only logistic model

Assume ONLY for this calculation:

1. Every ramet's relative light `x` lies in **0–100 percentage points**;
2. `p_S(x)=logistic(alpha_S+ln(1.009)*x)` with ONE intercept shared across all sites, months, plant ages and ramets; `p_A(x)=logistic(alpha_A+ln(0.979)*x)` with its own ONE shared intercept;
3. Each species' expected **marginal** count equals the printed 209 and 129 over 875 observations;
4. The species are independent conditional on `x`, **and** other host, month, source-site, detection and cohort differences have zero additional effects.

Since `p_S(x)` is monotone and its odds can increase by no more than `1.009^100` across the allowed light range, the mean `209/875` gives the rigorous pointwise bound

`p_S(x) >= logistic(logit(209/875) - 100*ln(1.009)) = 0.11355522`.

Multiplying by `sum p_A(x)=129`, the expected number of co-occurrences under conditional independence satisfies:

`E[both] = sum p_S(x)*p_A(x) >= 129*0.11355522 = **14.64862**`.

This is a **proven conservative mathematical bound** for any unknown distribution of `x` *within those exact restrictions*; it is **NOT** a general bound for butterfly competition, different site intercepts, nonlinear light effects, or **observed** co-occurrence. Seeing 5 remains stochastically possible even if expected values are larger.

As a more concrete (not globally optimized) stress test, concentrate relative light at its two extremes `x=0` and `x=100`, with unknown fraction sunny. Calibrate each intercept to print marginal expected S=209 and A=129. A grid over all fractions in increments 1/2000 finds the smallest expected overlap **22.93268**, at **sunny fraction 0.6205**. This minimum applies ONLY to the tested two-point light distributions, NOT all possible 0–100 light histograms.

Comparison:

| Calculation | Expected both | Logical status |
| --- | ---: | --- |
| Pooled naive independent incidence | **30.81257** | No environmental sorting |
| Model with reported ORs, two extreme light bins, grid-minimum | **22.93268** | Restrictive scenario, not a theorem for all RLI distributions |
| **Rigorous any-light-distribution lower bound** under a single intercept per species | **14.64862** | Theorem for this explicitly restrictive model |
| Published number of BOTH-positive ramet **records** | **5** | Observed count, NOT expectation |

It is correct to say: **the fitted-strength RLI coefficient alone in a common-intercept single-axis model is insufficient to make the mean co-occurrence as low as 5, given the pooled margins**. It is INCORRECT to say this rejects light heterogeneity in the real field or proves interspecific negative effects.

## C. Why this STILL does not prove competition: unrestricted source intercepts

Retain **exactly the authors' reported light slopes** `1.009` and `0.979`. Now permit species-specific intercepts for different unobserved environmental/site/host cohorts:

`logit(p_{S,g}(x)) = alpha_{S,g} + log(1.009)*x`, and analogously for `A`.

It is mathematically valid to assign the earlier 400- and 475-host groups to two otherwise distinct **hypothetical** source groups at the illustrative `x=50`, with each group's intercept calibrated to its printed group incidence probability. Species presence remains independent *within group*, expected total joint presence becomes **4.84947**, very close to the **5** published joint records, and the chosen integer group contingency realization reproduces the source's entire 542/204/124/5 total.

This is a second constructive **nonidentifiability example**, not a model fitted to two real Korean sites; group intercepts here are deliberately unconstrained and may represent host stem stage, microclimate, site, sampling or phenology. It proves mathematically why the single-intercept lower bound must never be transported to the published adjusted mixed models without original observations and random-effects structure.

## D. What this actually changes about the biology

The previously proposed contrast `host species available → immediate shared-host interaction` is insufficient. A more discriminating hierarchical question is:

**How much of the apparent butterfly co-use deficit is explained by measured site × date × relative light × plant size, how much by unobserved environmental structure, and only then what causal interaction follows from actual shared-plant encounter?**

A genuine empirical answer requires **real** ramet/quadrat survey records with stable unique stem IDs, nested quadrat/site, time/visit, larval presence for BOTH species including verified absences, relative light **in observed units**, plant leaf availability, plant size, host source, and survey effort. A joint occupancy or hierarchical encounter model must include measured environmental determinants BEFORE attributing residual segregation to interaction; a residual coefficient alone still cannot distinguish adaptive preference, spatial population dynamics or negative competition without manipulation.

In parallel, the Kyoto *Aristolochia debilis* species-pair cage work requires the independently specified randomized competitor presence/absence × stage-specific edible leaf access manipulation with adult viability, not host-cover proxies. Korea *A. contorta* host ecology is informative prior art, NOT independent causal validation of this Japanese experiment.

## File integrity, validation and STOP

- Frozen assumptions **before the calculations**: `SWALLOWTAIL_REPORTED_LIGHT_SLOPE_BOUND_PROTOCOL_V01.json`.
- Reproducible no-dependency script: `scripts/check_swallowtail_reported_light_slope_bound.py`.
- GitHub Action [light coefficient constraint audit](https://github.com/zuizui0223/chocho/actions/runs/38057988837): **6 unit tests** (mathematical bounds, source margins, old-toy implied ORs, counterexample with group intercepts and no competition assertion); JSON receipt uploaded. All calculations are on **published margins and conditional mathematical assumptions** only; no new field data, no newly identified butterfly ecological effect.
- The previously printed site/month data contain an unresolved **255 vs 215** quadrat denominator error; this continuation does not alter those counts or assert correction.
- **STOP** searching alternative invented light grids or site breakups to manufacture a desired result; the two mutually compatible explanatory structures already demonstrate the model identification boundary.
- **GEB submission PR #38 is unchanged**. The global introduced-HOST-PLANT availability paper is a separate estimand and cannot be made into a direct same-stem competition paper through these aggregated sources.
