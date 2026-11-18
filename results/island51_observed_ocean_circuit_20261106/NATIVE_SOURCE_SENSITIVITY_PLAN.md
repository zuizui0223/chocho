# Native AIS source-quality sensitivity: fixed after the primary analysis

## Why this additional check is required

The primary overview-based analysis and its biological results have already been seen. Its frozen protocol, input manifest, matrices, scores and official Circuitscape checks remain unchanged.

An independent source audit subsequently found that all 159,744 sampled first-overview pixels match the top-left native 2 × 2 pixel. Across 35 sampled cells with positive native means, the median absolute relative difference between overview-based and direct-native 0.5° means is 19.39%. This motivates a single **post hoc source-quality sensitivity**, not conductivity tuning or selection by predictive performance. This plan is saved before fitting any models that use the new native AIS grid.

## Fixed change and unchanged design

Replace only the AIS input by an area-overlap mean of every original **0.005° native raster pixel** on the existing exact 0.5° grid. Read the checksum/CRC-verified native TIFF, not its external overview or internal pyramid. Retain published integer units with scale 1 and offset 0, exclude original NoData, and retain valid zero including land-zero pixels. Record original source identity, transformation, checksum and grid coverage in a separate receipt. Validate the new grid against the independent 39-cell native overlap means before fitting biological outcomes.

Keep the original 51 islands, all 1,275 pairs, 182 geometrically qualifying ports, 25 km attachment threshold, coastline graph, terminal weights, uniform and wind matrices, original port-count controls, all four model comparisons, all three outcomes, regression procedure, transformations and fold assignments unchanged. In particular, AIS conductivity remains `1 + ln(1 + x)`. Do not choose a different transformation, resolution, denominator, geographic filter or model after seeing the native-grid scores.

## Required outputs and interpretation

Persist the native grid separately, then compute its distributed-port resistance matrix on the identical graph and repeat the same fixed comparisons. Save all 51 folds for each comparison and all three targets, the full pair table, source and protocol hashes, and a direct overview-versus-native score comparison. Uniform and wind results must match the primary run; they serve as implementation invariants for changing only AIS.

Report the primary and native-source results side by side, including deterioration or reversal. State explicitly that this sensitivity was chosen after viewing the primary analysis because a source-quality audit identified subsampling, and is not an independent validation set. If the AIS or combined effect depends materially on source averaging, qualify that effect; the uniform and wind-only conclusions stand on their own source inputs.

Do not rewrite the original pre-fit protocol or claim that this native comparison was pre-registered before the primary outcome analysis. The original official Julia validation covers the original three networks. Any official check of the new native-AIS network must be recorded separately with its actual scope.
