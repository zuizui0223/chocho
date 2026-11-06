# Geographic accounting of fixed host/circuit prediction gains

## Purpose and scope

Identify where the saved conditional host/circuit prediction gains occur among the original 51 islands, and how positive and adverse contributions combine into the previously reported overall result. The parent is PR128 / commit `17a4768cbb6cddf249449bd67e3ea1124e0db2f5`.

This is a descriptive, post-hoc geographic interpretation of existing predictions. The parent outcomes and aggregate results are already known. The accounting and display rules below are fixed before this atlas is calculated; this is neither a new predictive validation nor an independently preregistered ecological test. No models, graph surfaces, parameters, targets, source dates or island membership are changed or refitted.

The analyst selected this bounded next step under the user's instruction to continue the island/port/resource work autonomously. It makes the existing result geographically interpretable. Directional transport models and new environmental regressions require separate data/specifications and are outside this continuation.

## Frozen inputs and geographic identity

Use the complete parent prediction, fold, summary and provenance files, the original eligible-island roster, source-reconstructed original port groups and representative points, the frozen NOAA WPI2019 layer, the 194-port attachment ledger, and the source-preparation receipt. Pin every direct input and this protocol in `input_lock.json`. Run the existing independent parent auditor in process with `audit(parent_results)` before summarizing; it does not fit models or write files. Do not run its CLI on the parent directory because the CLI would rewrite a parent receipt.

Retain exactly 51 islands, 1,275 unordered pairs, 194 original biological ports and 182 marine-qualified ports. The archived predictions contain 25,500 rows, including two endpoint-specific predictions for every model/target/pair. All 510 parent comparison-fold records (five contrasts × two targets × 51 islands) remain in the atlas, including adverse and tied values. Each fold describes 50 test pairs and 1,225 training pairs.

Join source ports by `str(int(OBJECTID))`, never by feature order or port name. Require a one-to-one port ledger join, exact canonical island coverage, source name/coordinate agreement, and original/marine counts consistent with the parent source receipt. Preserve the 12 excluded marine ports in the complete port table; their absent marine coordinates and displacements remain blank, not zero.

Map `verified_source_51_mean_port_points`: the spherical means of each original port group, not island polygon centroids or means of marine attachments. Independently reproduce those spherical means, including dateline handling, and retain the original coordinates. Country codes describe sampled WPI port records, not the complete political coverage of the island. The source lacks an authoritative complete island-name crosswalk. Use canonical landmass IDs with explicitly labelled example-port names and sampled-country-code sets. Choose the example port by smallest numeric OBJECTID before adding gain values. Preserve every port name in the companion table.

The coastline-grid background may be derived from the frozen `fine_cell_to_node` mask, displayed as 0.5-degree averaged land fraction over its recorded extent (180°W–180°E, 85°S–85°N). It is a location aid, not a new landmass classification or a new circuit graph. Do not add basin/region groupings or infer actual shipping routes.

## Signed gains and additive accounting

For each saved island-fold comparison, define positive gain as reference loss minus augmented loss, separately for MSE and MAE. This is the negative of the saved `delta_mse_augmented_minus_reference` field. Validate both definitions against the parent values.

Let `g_i = reference_mse_i - augmented_mse_i`, `N = 51`, and `S = sum(reference_mse_i)` for that comparison/target.

- Island relative gain: `100 * g_i / reference_mse_i`. If that island's reference loss is zero, retain the row and record the percentage as undefined, even for a zero-versus-zero tie.
- Contribution to the overall mean MSE gain: `g_i / N`.
- Contribution to the overall percentage reduction, in percentage points: `100 * g_i / S`. The 51 contributions sum to the parent's relative MSE reduction. If the global reference sum is zero, retain absolute quantities and make relative quantities undefined.
- Apply the same additive and relative definitions separately to MAE; never combine MSE and MAE units.

The ordinary mean of island-relative percentages is not the global percentage. Show reference errors alongside gains because large absolute gains can reflect greater initial prediction difficulty. Strict improvement/worsening/tie signs follow the parent's exact comparisons, without an outcome-selected tolerance.

For each of all ten target/contrast combinations, define positive mass `P = sum(max(g_i, 0))` and adverse mass `A = sum(max(-g_i, 0))`. Save their mean-loss quantities `P/51`, `A/51`, and net `(P-A)/51`, and their percentage-point equivalents using S. Require those identities and the parent summary to agree within absolute and relative tolerance 1e-12.

Concentration uses absolute MSE gain magnitudes. Fix k = 1, 3, 5: report the share of gross positive mass contributed by the largest k positive gains, and the share of gross adverse mass contributed by the largest k worsening magnitudes. Also record the smallest contributor count reaching 50% of each gross mass. Use ascending numeric island ID to resolve exact ranking ties. If a mass is zero, its shares and half-mass count are undefined; its mass and contributor count remain zero. Keep ranked positive and adverse IDs with the summary. Do not divide positive contributions by signed net gain and call the result a share.

No fold deletion/refitting claim is made. No new regressions, correlations against port count or island attributes, permutation tests, p-values or independent-fold confidence intervals are added.

## The two conditional comparisons and geographic display

Primary later response: `new_butterfly_jaccard` (first indexed 2011–2022). Earlier response: `earlier_butterfly_jaccard` (2001–2010), retained as temporal context, not a null expected to be unrelated to geography. The two conditional contrasts remain C→HC (hosts beyond circuits) and H→HC (AIS + wind beyond hosts).

Join those two MSE signs for every target/island into five display categories:

1. Both gains positive.
2. Host gain positive and circuit gain negative.
3. Circuit gain positive and host gain negative.
4. Both gains negative.
5. At least one exact tie.

Retain both underlying signs and all MAE values/signs; additionally record whether both MAE gains are positive. “Both improve” means that HC beats both C and H for that island's 50 test pairs. Those contrasts share HC predictions and are not independent confirmations.

Show all 51 original representative positions with equal-size markers for both responses. Use the same category colors in both panels. A companion ranking figure/table must retain all 51 islands and both primary contrasts. If map labels are selective, use the union of the five largest absolute contributions from each primary later contrast, with numeric-ID tie breaking; identify them as descriptive selections. Labels must not replace a port name with an asserted island name. Show the actual signs and adverse contributions, and do not rank map-marker size by island-relative percentage.

## Outputs, limits and completion

Required tables: 51-row `island_geography.csv`, 194-row `original_port_geography.csv`, 510-row `island_gain_decomposition.csv`, 102-row `island_conditional_status.csv`, and ten-row `gain_concentration_summary.csv`. Persist the fresh parent audit evidence, strict support/identity checks, all input/output/implementation hashes, source labels, and complete sign-category counts. A compact parent audit receipt may omit its reconstructed table bodies if the original audit/source hashes and successful row counts/checks are preserved.

A geographic gain summarizes a held-out island's 50 pair predictions; it is not local abundance, establishment, a port-specific effect or proof of transport. Observed ranks are descriptive and are not independently validated geographic drivers. Shared endpoints and overlapping training data remain. Hosts pool 2001–2024, can postdate later records, and remain candidate feeding resources. Port counts only partially represent inventory effort; host richness, climate, land use and detailed recording effort are not independently resolved by this atlas. Prior positive host associations, resource-opportunity expansion and circuit prediction results retain their own estimands.

Complete when the frozen accounting runs, all parent summary identities and geographic joins pass, an independent reviewer reconstructs the new summaries, figures are rendered and inspected, meaningful tests pass both in the complete repository and within the existing anonymous-package boundary, and the continuation is saved in a draft PR with successful reproducible CI. Parent scientific artifacts remain byte-identical.
