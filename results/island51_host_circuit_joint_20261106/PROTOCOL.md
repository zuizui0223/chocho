# Joint candidate-host and observed marine-circuit comparison

## Question, decision and prior evidence

Do candidate introduced-host similarity and independently constructed marine circuit measures each add held-out-island predictive information about later butterfly-record similarity when the other is already included?

This is the next analysis after commit `a32975a45d372d2521ec9968bf4e4999fb540f09` / PR127. The previous circuit results and the original positive host–butterfly correlations are already known. This specification is fixed before fitting the **new joint models**, not before all prior outcome exploration. The new comparisons are exploratory extensions on the same sample; no new independent validation population is claimed.

The original 51 eligible islands and 1,275 unordered pairs remain fixed. Candidate-host Jaccard represents the union of introduced larval-host candidates recorded at the original completed ports, with unanimous regional classification, pooled over2001–2024. Retain that exact archived predictor, including zero similarities and the original ports that do not qualify as marine terminals. Never redefine the host inventory from only the182 marine-access ports.

Primary response is `new_butterfly_jaccard`: species first indexed on each island in2011–2022 within the source record history. Secondary response is `earlier_butterfly_jaccard`, for2001–2010. The earlier response is temporal context, not a null required to have zero association with modern geography. Neither response measures verified arrival or breeding. Do not fit `host_jaccard` as a response in this joint analysis because it is a predictor.

## Approach selected before the new fits

1. **Selected: five fixed joint prediction models.** This directly connects the original resource question with the now-verified circuit matrices using existing data and the original island holdouts. It can establish conditional predictive information, with limited mechanistic interpretation.
2. Directional wind transport would require a directed transport measure and directional/time-ordered biological outcomes. The present undirected resistance and symmetric Jaccard targets do not identify that process.
3. Adding new climate or occurrence-effort data could address additional confounding but would be a separate acquisition and model specification. Existing port counts are only a partial effort proxy; the present analysis does not independently control climate, survey days or data-provider effort.

The analyst selected the first approach under the user's instruction to continue the existing island/port/resource work autonomously. No parameter or model is selected by the new scores.

## Inputs and geometry retained

Use the frozen original pair/island tables; the saved uniform and wind distributed-port resistance matrices; and the **direct-native AIS** resistance matrix from the completed source-quality check. Direct-native AIS is chosen because it averages every native raster pixel, rather than a nearest-subsampled overview. Both source variants' prior predictive results are already known; this is explicitly the current source-quality choice, not a new blind source comparison.

Preserve the original coastline graph,182 geometrically qualifying ports from the194-port ledger, distributed port weights, original port-count controls and all51 islands. The source graph and all matrices have already undergone the independent official Circuitscape checks in PR127. Reuse their byte-identical matrices; this extension does not rebuild circuits or choose new transport surfaces.

Verify all files listed in `input_lock.json`, the protocol hash and the prior frozen source bundle before loading data for the new fits. Matrix island IDs must equal the canonical source order; matrices must be finite, symmetric, zero-diagonal and positive off-diagonal. Port counts must agree between the archived island table and original source-preparation receipt. Retain all1,275 pairs with their original targets.

## Five models, unchanged estimation

Let `B` contain log1p(original great-circle distance), its square, the sum and absolute difference of log1p(original port counts), and—for the later response only—earlier butterfly Jaccard.

For each circuit matrix, retain the previous feature `r = log1p(R / median_positive_pair_R)`. The normalization uses response-independent geography over the fixed island set. Define `c_A = -r_A` and `c_W = -r_W` only for interaction construction: higher values mean greater circuit connectivity. Let `h` be the archived host Jaccard.

| Model | Predictors beyond B |
|---|---|
| G | uniform-ocean resistance feature |
| H | uniform + h |
| C | uniform + native AIS + wind resistance features |
| HC | uniform + h + native AIS + wind |
| I | all HC terms + h×c_A + h×c_W |

The two interaction terms enter together as **one secondary block**. Retain all constituent main effects; no stepwise selection, thresholds, outcome-selected transformations or alternative interaction search. Raw feature products are formed without using response values, then all columns are standardized using training rows only. Since all main effects and the intercept are included, centering the products is not a different prediction model.

Use the existing ordinary-least-squares solver, training-only feature standardization, intercept and predictions clipped to[0,1]. Do not change fitting method or loss to obtain improvement.

## Fixed comparisons and holdouts

The two primary conditional comparisons are **C→HC** (host information beyond circuits) and **H→HC** (AIS+wind jointly beyond hosts). Report **G→H** and **G→C** for context and exact prior-result reproduction. Report **HC→I** as the secondary interaction comparison. Evaluate and retain all five comparisons for both responses, irrespective of direction.

In every fold, withhold one island and all50 pairs containing it; fit on the1,225 pairs among the other50 islands. This predicts one unrepresented island paired with otherwise represented islands, not two entirely unseen islands. Use all51 folds. Every unique pair has two held-out predictions, one for each withheld endpoint. Features for the test island, including its host inventory and geography, are available; this is retrospective geographic transfer, not forecasting an island with unknown resources.

Fit each of the five models once per fold and response. Save all25,500 individual predictions (5models×2responses×51folds×50pairs), all510 model-fold losses, and all510 comparison-fold records. For each comparison report reference/augmented MSE and MAE, percentage MSE reduction, counts of improved/worse/tied islands, and all fold losses. Percentage reduction is100×(reference−augmented)/reference; never add percentages with different reference denominators. No pair-independent or independent-fold p-values, confidence intervals or standard-error bars.

Require G and C to reproduce the existing native-AIS comparison for both butterfly responses, including every held-out fold, within absolute and relative numerical tolerance1e-12. Stop before interpreting the new models if this identity fails.

## Interpretation and completion

Positive C→HC supports additional resource-composition information conditional on the represented marine fields. Positive H→HC supports additional circuit information conditional on the recorded candidate-host composition. Either can reflect residual climate, biogeographic, recording-effort or other geographic structure. Improvement in the secondary interaction model is descriptive predictive nonadditivity and does not establish a transport×resource mechanism.

Preserve all existing positive results and distinguish their estimands: original resource-opportunity expansion, host–butterfly correspondence, marine-circuit prediction, and the present conditional comparisons. An adverse or small new incremental score does not overturn a marginal correlation or the global resource-opportunity result. Host records can postdate butterfly records; no temporal mediation, causal plant carriage, verified butterfly movement or establishment claim is licensed.

Complete when the frozen comparison has run; all predictions and losses are independently reconstructable; meaningful leakage/source-alignment tests and the repository suite pass; a concise result report and exact-data figure are saved; and the changes are persisted in a continuation branch and draft PR with reproducible CI.
