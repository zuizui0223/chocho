# Temporal host-supply timing versus butterfly competition: an explicitly hypothetical sign-reversal model

**Date:** 2026-10-14. **Status:** ANALYTIC THEORY / SYNTHETIC EXAMPLE / NOT FITTED TO BUTTERFLY OUTCOMES.
**Branch:** theory/temporal-resource-competition-v01, derived from exploration/realized-host-window-v01.
**Submission boundary:** Do not merge into the GEB main manuscript or infer recruitment from HOSTS × WCVP host presence.

## Biological question and distinguishing result

At **identical total seasonal supply of actually edible host material**, can shifting the *calendar timing* of supply reverse which of two co-occurring butterfly species is more harmed by a heterospecific competitor?

The relevant field endpoint is the **number/proportion of originally assigned larvae that produce flight-capable adults**, not host occupancy, egg counts, overlapping mapped host ranges, adult observations, or a calculated union of available dates.

Seasonal resource pulses, stage-mediated priority effects and temporally modulated competition already have theoretical precedent. The toy result below is **not a claim of a new general theoretical mechanism**. The open empirical result would be a controlled, adult-fate-validated *butterfly-specific* heterospecific-competition sign reversal, with host tissue quantity and stage documented.

## Fully stated assumptions and equations

Two species A and B share an edible host resource across two **calendar** periods E (early) and L (late). Cohort food demands are $d_{i,t}>0$. Edible supply in period $t$ is $R_t$. Food cannot be stored or carried between periods. If total demand exceeds supply, food is allocated proportionally to contemporaneous demand:

$$
q_{i,t}=d_{i,t}\min\left(1,\frac{R_t}{d_{A,t}+d_{B,t}}\right).
$$

The fractional unmet demand for each species is

$$
u_i=\frac{\sum_t(d_{i,t}-q_{i,t})}{\sum_t d_{i,t}}.
$$

**Solely for the synthetic example**, assume an exponential link between fractional unmet demand and the probability of a flight-capable adult:

$$
p_i=p_{0i}\exp(-k_i u_i),\quad k_i>0.
$$

Neither $k_i$ nor $p_{0i}$ is estimated from actual butterflies. This is a demonstration that an ecological effect is logically possible under specified assumptions, not a fitted demographic function.

Let $T_E=d_{A,E}+d_{B,E}$, $T_L=d_{A,L}+d_{B,L}$, and let the total experimental resource deficit be $\delta$. Compare:

- Early deficit: $(R_E,R_L)=(T_E-\delta,T_L)$.
- Distributed deficit: $(R_E,R_L)=(T_E-\delta/2,T_L-\delta/2)$.
- Late deficit: $(R_E,R_L)=(T_E,T_L-\delta)$.

All arms supply exactly $T_E+T_L-\delta$ material. Require $0\le\delta<\min_{i,t}d_{i,t}$ so each species *alone* is adequately supplied under all three timing treatments; resource limitation appears in mixed-species cages.

Define the modeled competition penalty for species $i$ under regime $r$ as a **negative log risk ratio** relative to its own solo-cage control under that same timing:

$$
H_i(r)=-\log\left(\frac{p_i(\mathrm{mixed},r)}{p_i(\mathrm{solo},r)}\right), \qquad
D(r)=H_A(r)-H_B(r).
$$

A **relative competitive-disadvantage reversal** is $D(\mathrm{early})>0$ and $D(\mathrm{late})<0$. This is not equivalent to demonstrated population competitive exclusion or evolutionary fitness reversal.

Let $a=d_{A,E}/\sum_t d_{A,t}$ and $b=d_{B,E}/\sum_t d_{B,t}$. Then under these assumptions,

$$
D(\mathrm{early})=\frac{\delta(k_A a-k_B b)}{T_E},
\qquad
D(\mathrm{late})=\frac{\delta(k_A(1-a)-k_B(1-b))}{T_L}.
$$

For $a>b$, a reversal occurs exactly when

$$
\boxed{\frac{b}{a}<\frac{k_A}{k_B}<\frac{1-b}{1-a}.}
$$

Thus:
1. **Distinct relative demand timing is necessary** for this particular modeled direction of reversal. If $a=b$, no positive $k_A/k_B$ yields a sign reversal.
2. Even when temporal demand profiles differ, **strongly unequal scarcity-to-emergence sensitivities can prevent reversal**.
3. The prediction depends on contemporaneous allocation, no inter-period carryover, and the chosen dose-response rule. Plant induction, predation, enemy dynamics and development plasticity could change it.

## Transparent worked example (NOT field observations)

Choose the following invented demands: A=(0.75, 0.25), B=(0.25, 0.75). Choose $\delta=0.20$, $p_{0A}=p_{0B}=0.8$, and $k_A=k_B=3$. The seasonal food total is 1.8 in **every** regime.

| Treatment | Supply E / L | Unmet fraction A / B | Assumed adult-emergence probability A / B | D(r) |
| --- | --- | --- | --- | --- |
| Early deficit | 0.8 / 1.0 | 0.15 / 0.05 | 0.510 / 0.689 | +0.30 |
| Distributed deficit | 0.9 / 0.9 | 0.10 / 0.10 | 0.593 / 0.593 | 0.00 |
| Late deficit | 1.0 / 0.8 | 0.05 / 0.15 | 0.689 / 0.510 | -0.30 |

Here the admissible sensitivity-ratio interval is $1/3<k_A/k_B<3$. This is a **construction**, not a significant effect, estimated coefficient, prediction about *A. alcinous*, or evidence for a general butterfly ecological law.

Executable negative controls enforce:
- $\delta=0$ gives no competitor penalty;
- identical normalized demand timing gives no relative sign reversal, even if total food requirement differs;
- $k_A/k_B$ outside the derived interval prevents the illustrated reversal;
- changing total species demand (e.g. 2.6:1) alone does not guarantee or prohibit the result.

## What the Kyoto butterfly experiment actually establishes

Hashimoto & Ohgushi (2023), DOI https://doi.org/10.1002/ece3.10164, reported food demand of **1,024.89 mg** for *Atrophaneura alcinous* and **396.13 mg** for *Sericinus montela*. Both consumed most of this during instars 3–5 (1003.23/1024.89 and 394.25/396.13 respectively). Early/late **instars are not early/late calendar weeks**. Opposite calendar-time peaks in larval food demand are not established; actual adult emergence was not independently quantified as our endpoint.

The archive has **690 repeated larval stage observations from 30 cages**. Its four final defoliation variables refer to **four distinct plants within each cage assessed at one endpoint**, not time-series biomass or assigned regrowth timing. The published result already showed asymmetric density responses. The 2017 paper also investigated regrowth; neither previous result can be relabeled as a new timing manipulation.

A **calendar-time larval demand** pilot must first measure egg-placement dates, development rates, instar-specific consumption by date, overlap of the two species' cohorts, and the supply of genuinely edible new tissues. A source-supported 2.6-fold *season-total* food-demand difference is not evidence for the demand-profile offset in this toy model.

## Public-data role and originality threshold

- Hashimoto & Ohgushi 2023 archived stages and density design can check cage/instar and density feasibility, **not** the effect of a randomized supply-time intervention.
- Erickson, Schultz & Crone (2025), DOI https://doi.org/10.1002/ecs2.70259, independently counted monarch eggs/larvae across monthly plant surveys and can inform observational seasonality, **not** egg-to-adult cohort survival or the two-species competitive sign reversal. The chocho route-forward host-identity prediction did not validate robustly.
- The original chocho 239-butterfly GEB potential host geography remains a **separate estimand**. Moving a host on a map never implies local reproduction.

**Prior theory:** Holt (2008), resource pulses, DOI https://doi.org/10.1890/07-0348.1; Zou, Schreiber & Rudolf (2023), stage-mediated priority effects and season lengths, DOI https://doi.org/10.1098/rspb.2023.1217. The novelty cannot consist merely of discovering that timing matters.

## Forward, falsifiable experiment

**Design:** Randomize host enclosures/plants to 3 timing arms (early / distributed / late supply) crossed with focal species and competitor absent/present. Block by plant batch, parental family, experimental date, and ideally independent source populations. Hold original focal cohort size fixed. A separate **density-matched conspecific addition** arm is desirable if the claim specifically contrasts heterospecific competition with density/crowding effects.

**Measured manipulation:** Keep season-integrated fresh *edible* mass, plant species, initial biomass, handling and substrate equivalent; regularly check *achieved* time-resolved edible mass and foliage stage. Record new-leaf quality and induced effects. The total mass assigned is a design constraint, not evidence that actual accessible resources were identical in total or quality. Do not condition only on final defoliation or surviving larvae.

**Primary outcome:** Number of flight-capable adults per *originally assigned* focal larvae, with failures and zero-emergence cages retained. Primary unit of randomization and statistical replication is the **cage/plant**, not larvae inside the cage. Estimate and display $H_A(r)$, $H_B(r)$, and $D(r)$ with uncertainty and absolute emergence risks.

**Prediction gate:** Genuine sign reversal requires independently supported $D(\mathrm{early})>0$ AND $D(\mathrm{late})<0$ under achieved timing separation, not merely a small interaction p-value or opposite signs of noisy point estimates. An interaction without a sign flip supports *temporal modification*, not sign reversal.

**Critical stop rules:** If pilot cohorts do not show distinct calendar-time food-demand profiles, do not assert the sign-reversal hypothesis for this pair. If the intervention does not yield separated achieved leaf-access trajectories, classify the experiment as inconclusive. If suppression remains with continuously adequate edible supply, test food-quality induction, direct interference and enemies rather than forcing a biomass-only conclusion. No inference about introduced *plant* effects is possible from the Kyoto system: the shared plant is native, while one butterfly is introduced.

## Reproduce without external data

After installing the existing project in editable mode, run the demonstration with:

    python scripts/simulate_temporal_supply_competition.py

Run the synthetic mathematical checks:

    pytest -q tests/test_temporal_supply_competition.py

**Final status:** DERIVED and UNIT-TESTED analytic prediction only. NO causal butterfly result, validated demographic conversion, or evidence for host-redistribution fitness benefits. Do not modify GEB PR #38.
