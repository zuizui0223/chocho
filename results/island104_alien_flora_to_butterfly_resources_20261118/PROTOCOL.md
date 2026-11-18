# Fixed direct bridge: alien floristic similarity → potential butterfly larval-resource similarity

**Date:** 2026-11-18, protocol committed before this *new 104-island bridge model* is fit. Exploratory study, designed after earlier WCVP/GloNAF/global findings were already known. Does **not** identify independent ecological cause or actual butterfly assembly.

## Frozen source and exact question

For **104** GSHHG-connected WCVP strict island botanical regions already used by the native/introduced flora ocean-circuit test, join **5,355** botanical region pairs with Jaccard defined in *both* native and introduced WCVP flora to the **original 239-butterfly potential-resource bit matrix** from the archived PR #118 commit `bb4d8f00fbdb97fa2d90f606a8c9e7e0620646f2` (39? no new insect-source data). Use the original species order and exact 26,530 native/14,553 gained/41,083 contemporary resource bits as immutable input checks. Unmatched island/region codes fail rather than being inferred or geo-nearest mapped.

Question: Does **accepted introduced WHOLE plant species-composition Jaccard** of island pairs help conditionally predict how similar their **contemporary sets of potentially resource-supported butterfly species** are, after controlling geographical and native host-resource structure and the number of opportunities? Plant taxa and butterfly species are DISTINCT response universes. The consumer resource Jaccard is **built from some of the same WCVP flora and HOSTS**: no truly independent cross-trophic observation claim.

## Fixed five hashed folds / test population

Use exact original source-only `chocho-WCVP126-double-unseen-v1|code` SHA256 fold labeling, five bins, 1,089 same-bin pairwise held-*both*-endpoints unseen test dyads of the original 5,355 pairs. No test island used as any training dyad endpoint. Use original source `src.ridge()` with training-only standardization, ridge alpha **20**, and clip predictions to [0,1]. No hyperparameter tuning, no response-selected dropped islands, no fold selection.

## Models frozen before seeing results

Response **Y_current** = Jaccard of 239 fixed butterfly species whose *currently recorded native or introduced known larval-host plants* are present in the island, conditional on the original WCVP/HOSTS mapping. Retain separate **Y_native** resource Jaccard as baseline explanatory feature (NOT as outcome for the alien-flora hypothesis). If both current sets empty, mark undefined and report, never impute as zero. Original 355-region support filtered only when explicitly describing native-only resource Jaccard; source-pair model includes source 104 plant-defined pairs and reports any undefined resource pairs.

For each island pair, fixed baseline E:
- log1p geodesic distance and its square, TDWG Level1 equality, absolute latitude difference, log1p(min/max original island area), log1p modern port-count product.
- native *whole vascular flora* Jaccard; native *potential butterfly resource* Jaccard; zero if native resource Jaccard undefined **only as an explicitly separate missing-source coding with missingness flag**, never disguised as exact observed zero. A source-native-resource-empty island is retained only if source-defined missingness is recorded and a sensitivity excluding these dyads is separately shown.
- log1p(min/max **native** butterfly potential-resource richness), log1p(min/max source **native** and **introduced** whole plant richness).

Baseline S adds log1p(min/max **contemporary butterfly potential-resource richness**) to E, and is **strictly retrospective** because it conditions on the target-derived current resource margins. E and S are different model targets; do not relabel S as a prospective prediction. Two primary comparisons: `E → E+alien_WCVP_plant_Jaccard` and `S → S+alien_WCVP_plant_Jaccard`. Fit and evaluate same dyads and fixed geographic folds for each nested pair. This measures incremental *geographic information* of shared alien whole flora, not a host-specific edge effect.

Prespecified negative/sensitivity:
- same E and S model pairs on **native-supported-islands only**, deleting all island pairs with source-empty native candidate-host resource to avoid missingness coding claims; new denominators reported; original 1,089 test pairs *not* silently reused.
- as an internal descriptive baseline, report source Kendall/Spearman rank co-association between introduced whole-flora Jaccard and contemporary resource Jaccard only in all pair cohort, noting shared WCVP source and dyad nonindependence. A raw correlation is not an independent test.

Score = held-island equal-fold MSE reduction % (positive better), pair-weighted sensitivity, 5-fold source-resampled descriptive 95% intervals (seed 20261118, 4,999 draws), per-fold signs, n-dyads, response-zero or absent rate and relative error. Do not convert descriptive 5-fold intervals into iid p-values.

## Scientific decision criterion / mandatory reporting

A positive alien-whole-flora incremental gain **is source-circular to some degree** because the potential-resource response was constructed from WCVP distributions and fixed HOSTS; only independent GBIF 224/155 indexed butterfly responses can test realized consumer community geography, and their host-specific out-of-geography gains are far smaller. If the 104-island incremental gain turns negative under S, report that without downranking earlier +54.9% and direct endemic larval use. If E is positive but S near zero, interpret as potentially shared **quantity of available resource**, not necessarily shared composition.

Do NOT equate WCVP introduced whole-flora similarity, particular HOSTS larvae, or physical sea-current/cargo paths. **No colonization, plant arrival timing, host quality, population increase, or competition estimand is fit.**

Original data receipts: WCVP source plant pairs [GitHub CI #44010307551], island coastline source [#44010673580], original 239×369 resource matrix [PR #118, immutable SHA], original native/current opportunities [GEB]. Full results and adverse folds in Supplement, not new headlines chosen post-result.
