# Fixed longitude-block double-unseen island-pair evaluation

## Purpose and prior knowledge
The earlier 51-island leave-one-island-out comparisons in PR128, PR130 and PR131 were seen. In PR131, host Jaccard added to native AIS/wind circuits + marginal host richness + early all-botanical recording effort improved later MSE +2.0211%, and circuits added to hosts + richness + effort improved +3.1576%. This new, stricter out-of-block design is specified after those findings. It is exploratory retrospective validation on the SAME source population, not an independent preregistered replicate.

## Source and covariates
Preserve exactly the original 51 landmass IDs, 1,275 unordered pair records, 225-butterfly original candidate-host panel, 194 biological ports, 182 marine-qualified ports, original 2001–2024 pooled host Jaccard, marginal host-count controls, archival 2001–2010 all-plant port-buffer record-count controls, and pre-existing Circuitscape uniform, native AIS and wind resistance fields. Do not add any new plant or butterfly taxon, geography or circuit; pin the audited atlas's island_geography.csv Git blob SHA.

## Fixed, target-free 5-block partition
Take original landmass representative-longitude coordinates from the fully audited atlas, sort by ascending longitude, breaking exact longitude ties by numeric landmass ID. Divide the sorted list into five ordered blocks of 11, 10, 10, 10, 10 islands. Store every island-to-block assignment and each block's longitude range. These blocks are purely longitude strata, NOT ecological regions, independent continents, shipping corridors or natural biogeographic boundaries. The first block may cover a wide longitude range; do not retune boundaries based on prediction scores.

## Primary evaluation: both endpoints absent from training
For each unordered block pair (g,h), g<=h (15 folds), remove ALL observations involving ANY island in either block from the training set. The test set contains precisely the original unordered island pairs with one endpoint in block g and one endpoint in h (both in the same block for g=h). Every original unordered pair is tested exactly ONCE, yielding 1,275 fixed primary predictions per model and target. Never include a test pair, either endpoint, or another incident pair in its training fold. Within-block test sets: 235 pairs across five folds; between-block sets: 1,040 pairs across ten folds. Single-block training support is 780 or 820 pairs; two-block training is 435 or 465 pairs. These train sizes and per-fold test sizes are fixed by the assignment. Do not pool held-out endpoints into original leave-one-island-out score or compare the resulting percentages as if they shared the same train/test split.

## Fixed models
Use PR131's three candidate-host-richness + pre-2011 botanical-recording-effort models, without a feature change:
- CRE = uniform + AIS + wind circuits + candidate-host richness controls + geographic/port/earlier-butterfly controls + early botanical effort;
- HRE = uniform + candidate-host Jaccard + candidate-host richness controls + geographic/port/earlier-butterfly controls + early botanical effort;
- HCRE = their combination.

Compare CRE→HCRE (candidate-host Jaccard beyond circuits) and HRE→HCRE (native AIS/wind circuits beyond hosts), for both later 2011–2022 first-indexed butterflies and earlier 2001–2010 context, with all signs retained. OLS with intercept, train-only predictor standardization and [0,1] clipped predictions, exactly as PR131. Primary metric: pair-weighted pooled held-out MSE over all 1,275 pairs. Secondary metric: MAE. Also report the separately conditioned 235 within-block and 1,040 between-block pairs. Report group-fold improvements as purely descriptive (folds share training islands) and individual pair improvement counts as also nonindependent. Do not estimate naive p-values/CIs.

## Requirements and inference boundary
Verify the 51-island/1,275-pair source and geographic identities, historical source GBIF archive hashes (as in PR131), frozen network/circuit source locks, exactly 15 disjoint test partitions, no endpoint leakage, and no pair reweighting. Save full predictions, all 90 model-fold scores, 60 contrast-fold scores, 12 prespecified comparison summaries and an independent standard-library auditor that recomputes each error, fold, comparison and split/source identity. Record both positive and adverse results. Retain parent predictions/source tables unchanged.

**Crucial limitation:** host Jaccard still pools 2001–2024, so post-butterfly host entries and detection-effort covariance remain possible. Missing original historical plant-species×port data prevents a valid pre-2011 host-identity rebuild. Geographic blocks are longitude strata of already selected islands, not independent regional replicates, and a two-unseen-endpoint comparison is not a forecast of colonization or ship-mediated transport.
