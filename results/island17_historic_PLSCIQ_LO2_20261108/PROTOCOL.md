# Pre-fit protocol: historical UNCTAD 2006–2010 port connectivity and later butterfly communities

Frozen November 8, 2026. Exploratory observational predictiveness, not independently preregistered: previous GEB +54.9% potential larval-host resource expansion, positive later butterfly AIS signal, adverse earlier butterfly/general flora effects and PR137 negative local AIS-label swap results already known. No original manuscript modified. No journal submission.

## Research question

Can port-level historical UNCTAD PLSCIQ liner connectivity INDEX in 2006–2010 predict **later 2011–2022 first-indexed butterfly assemblage Jaccard similarities**, beyond geography, source port opportunity, historical botanical record effort, and *earlier 2001–2010 observed butterfly Jaccard*? This is a proxy for port connectedness, NOT verified bilateral vessel cargo journeys or butterfly colonization. The earlier 2001–2010 butterfly similarity is a descriptive, contemporaneous contrast with 2006–2010 PLSCIQ, NOT a historical forecast.

## Frozen outcome-independent source cohort

Original 194 WPI ports/51 island groups/1,275 pairs and butterfly source observations remain unchanged; original WPI geographical Git blob a4681c3e3c1406d3a37f51464a298286e9d0ee4f. From the 24 individually PortWatch-geoverified original WPI ports in Draft PR148 (15+9 distinct, 17 islands), keep only old UNCTAD PLSCIQ port series with at least **16 nonmissing of 20 quarters 2006Q1–2010Q4**. This selects **22 source original WPI ports** spanning **17 original island groups**, yielding **136 unordered original source island pairs**. The other two historically under-indexed ports are NOT biologically absent. This historical source version was exactly recovered in Draft PR146: three original API response SHA-256 values identical to prior 2026-10-30 source acquisition, 1300 port series, 81900 values and 57312 nonmissing. Never use a newer index edition silently. Original 2018 GLSN model edges/modern AIS OD are NOT source predictors in the primary model.

At each selected island: the prespecified primary index is median of the port-level historical 2006–2010 *mean PLSCIQ index* across those source-qualified verified WPI ports; alternative uses MAX over them (hub-opportunity sensitivity), with log1p transform. Original Borneo-like landmass group ID 9 includes multiple countries and ports; keep its original island geography and preregister a separate source-only excluded-group-9 sensitivity, do not split or drop after looking at butterfly outcomes. Save all original 51 island source eligibility, 22 port index units, and 136 original pair identities.

## Frozen model design

Both targets: later 2011–2022 butterfly Jaccard (PRIMARY) and earlier 2001–2010 Jaccard (contemporaneous context). Baseline: log1p original great-circle port distance and squared term; sum and absolute difference of log1p original WPI botanical port counts; original representative source latitude absolute separation; sum and absolute difference of log1p original pre2011 all-general-plant WPI-buffer observation record counts. The LATER target baseline additionally includes earlier butterfly Jaccard (prior period). NO 2001–2024 pooled candidate-host network, modern 2015–2021 AIS, or 2018 GLSN model edge in these time-correct baselines.

Exactly four model families per target and population:
BASE geography+previous butterfly+pre2011 botanical sampling effort;
MED_SUM = BASE + sum of the two island log1p historical MEDIAN PLSCIQ;
MED_SUM_DIFF = MED_SUM + absolute difference of the same island-index values;
MAX_SUM_DIFF = BASE + sum and abs difference of island MAX historical port-index values.

Five predeclared contrasts: BASE→MED_SUM, BASE→MED_SUM_DIFF, MED_SUM→MED_SUM_DIFF, BASE→MAX_SUM_DIFF, MED_SUM_DIFF→MAX_SUM_DIFF (nonnested descriptive difference, NOT causal increment). Positive relative MSE/MAE reduction means better out-of-sample pairwise prediction. Save all poor as well as good comparisons. Do not form naive p-values/CI from correlated island pairs.

## Strict double-island-unseen forecasts

For every eligible 136 island-pair heldout test, train ONLY on eligible island pairs with NEITHER test endpoint included; each train set has C(15,2)=105 source pairs, exactly ONE test pair. Fold-specific model column centering/scaling uses TRAIN pairs only, OLS with intercept, clipped predicted Jaccard to [0,1]. Repeat fixed source-only sensitivity omitting source group 9: 16 eligible islands / 120 test pairs / C(14,2)=91 training pairs per fold. Both targets × 4 model families × (136+120) heldout pairs = **2048 prediction rows**; 2 targets × 5 contrasts × 2 scopes = **20 comparison rows**. Primary 17-island model first; 16-island sensitivity retained in same output, no group-9 result-based selection.

Archive source island index ledger for all 51 islands, every original test pair ID, all 2048 model predictions with errors/target/fold and 20 comparison rows, plus exact input SHA receipts. Independently written standard-library audit rebuilds all 22 eligible source port IDs/means and 17 island medians/max, all 136/120 source pair indices, double-unseen train size and observed source outcomes, all 2048 prediction losses, comparison MSE/MAE, and no pair missingness. Independent offline refit must confirm the coefficients/predictions through fresh run; error-check audit alone is NOT independent proof of fitted model coefficients.

The sampling frame is NOT the global original 51 island population; 17 islands are selected by historical index/provider coverage and independent registry geography, with some countries combined into a single original landmass. Neither positive nor negative prediction gains identify actual vessel movements, introduced hosts, butterfly breeding, or directed colonization. Prior original biological studies remain unchanged.
