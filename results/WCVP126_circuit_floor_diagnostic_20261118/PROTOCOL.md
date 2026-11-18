# Frozen exploratory diagnostic — does the native near-zero similarity floor explain the electrical-circuit contrast?

**Date:** 2026-11-18, after observing native −0.044% / introduced +2.735% strict both-endpoints-unseen Jaccard MSE gains. **Do not represent subsequent variants as preregistered independent confirmations.** This evaluates one alternative explanation of the original island–port–flora story.

## Original population and invariant predictors

Reuse 126 WCVP accepted-name botanical regions, original 104 GSHHG-coast-anchored regions, 5,355 pairwise *both-origin-Jaccard-defined* shared native/introduced flora pairs, five original source-only SHA256 islands bins, strict test set requiring BOTH island endpoints unseen in training (original 1,089 within-fold tests), original current WPI pair-index and original GSHHG uniform finite-volume electrical resistance, unchanged ridge penalty alpha=20 and training-only standardization. Original geographic baseline B is log1p geodesic distance and square, source Level1 equality, latitude gap, min/max same-origin WCVP accepted species richness and area; P uses modern WPI pair port counts; C adds circuit resistance; PC both. **No island heldout endpoints may be seen during fitting.** No response-selected folds or penalty tuning. First reproduce exact source Jaccard B→C native −0.04430509181% and alien +2.73504068049%, plus the source matching 1,089 test pairs, as a fail-closed gate.

## Two non-equivalent ways in which near-zero Jaccard could matter

**Support/zero-floor (primary descriptive diagnostic).** For all 5,355 source pairs, and separately 1,089 strict test pairs, compute number and share of pairs having Jaccard exactly 0 (no shared accepted species), positive Jaccard proportion, descriptive medians for all and nonzero pairs, and nonzero interval/quantiles, native vs alien. A low native Jaccard *alone* does not algebraically enforce zero MSE percentage gain; scaling a whole target by a constant leaves MSE *ratio* unchanged. A high proportion of exact zero outcomes, however, can make an undirected continuous descriptor weak or useful in different ways.

**Presence-of-any-shared-taxon (primary model sensitivity).** Convert the native and introduced Jaccard target into binary `Y>0`. Fit the **same original ridge linear regression and clipping** on each outcome with the B/C/P/PC predictors and original five held-island folds, evaluate Brier (binary MSE) on exactly the original 1,089 strict pairs. Report equal-fold/pooled-pair Brier relative changes, each of five individual fold scores and #events; also compare fitted Brier to a same-training-prevalence-only null, never claim ROC without balanced event support. Any sign gain reports only how ocean structure predicts *sharing at least one accepted native or alien species*, not community ecological movement.

## Additional descriptive sensitivities (not new claims)

**Overlap magnitude among observed positive pairs:** Fit original B and C on positive Jaccard training pairs and evaluate on positives within held fold only. Before interpreting, require >=30 positive training pairs and >=5 positive held pairs per fold for BOTH origins; otherwise mark as NOT_ESTIMABLE. This is a post hoc outcome-conditioned population, not a forecast that predicts positivity. Also consider original whole-pair B/C scores restricted to same-positive test pairs for both origin classes; report the new denominator explicitly.

**Alternative community overlap denominators:** Exactly the source same-pair accepted-species Sørensen and Simpson overlap columns (not 2 extra independent biological datasets). Repeat the B→C strict source fold/test comparison on their unaltered original numerical values. A Sørensen transformation of Jaccard and Simpson overlap differ in denominator and zero prevalence. Report individual fold and descriptive common-fold intervals; no cherry-picking whichever transform is positive.

## Decision rule

Promote the native-vs-introduced *predictive* ocean circuit difference above the existing nonconfirmatory status only if multiple non-redundant response definitions are coherent across source five folds *and* paired fold uncertainty excludes zero, not merely point estimate sign. Even then no measured directional sea currents, cargo delivery, seed viability, historical introduction or butterfly dispersal is identifiable. The primary article remains a resource redistribution result; original positive +4.239% and strict +2.735% values must be retained even if floor tests are adverse.

## Source and reproducibility

- Original source data archives: original accepted WCVP pairs [run 44010307551], original source shortest-water geometry [44010673580], original full GSHHG *electrical* circuit [44013658955].
- Source model and five hashed folds: `scripts/audit_WCVP126_crossed_native_alien_geometry_port_predictions_20261117.py`, `scripts/audit_WCVP126_native_alien_double_unseen_islands_20261117.py`, `scripts/audit_WCVP126_uniform_circuit_double_unseen_20261118.py`.
- Original strict electrical source [run 44290887759], [source-locked receipt](../WCVP126_circuit_strict_geographic_transfer_20261118/VERIFIED_STRICT_CIRCUIT_SCORE.json).
