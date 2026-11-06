# Observed marine circuit continuation: fixed analysis specification

This specification precedes fitting the new AIS/wind circuits to the biological outcomes. The existing 51 eligible islands and all 1,275 archived pairs remain fixed. No coastal resolution, port allocation, predictor transformation or model is selected using the new prediction results.

## Geometry and source recovery

The original 1° center-only sea mask failed a response-independent geometry check: Mediterranean and Baltic anchors occupy separate components; three Danish anchors coincide; the mean location of Borneo ports is 296 km inland of a sea center. These failures are retained in a separate audit.

Audit the194 original WPI2019 ports in the geographically reconstructed51 source groups. Verify the original WPI SHA and all1,275 archived centroid distances. Attach each actual port independently to its nearest GSHHG **marine** subcell within25km. Marine means the largest source-defined connected sea component; isolated inland waters do not count as ocean. No greedy unique-cell assignment or invented marine bridge is permitted. Every original port remains in the audit and every island remains in the biological analysis.

Resolve sea topology using original GSHHG2.3.7 intermediate coastlines, with four orthogonal neighbors and longitude wrap. Within each0.5° output block, keep distinct connected water patches as distinct graph nodes. Adjacent graph nodes exist only when fine sea cells share a face. Initial0.125° subcells left two distinct island port distributions aliased; the predefined finer0.0625° subcells resolve that geometry. The fixed primary therefore uses0.0625° coastline subcells and0.5° conductance blocks. Aggregation preserves the fine graph's connectivity; it is not an exact electrical/Kron reduction. Fine coast geometry does not create finer AIS or wind observations.

### Geometry-only amendment, before fitting any new biological model

The initial requirement that all194 ports directly enter an ocean-only graph was infeasible: the original WPI roster includes inland river ports, and coastline data contain no navigable river network. The final scope is therefore the independently defined coastal access points within25km of the largest marine component. The complete194-port audit is retained, including exclusions; original plant/butterfly outcomes and all51 islands remain unchanged. This is a marine connectivity analysis, not a complete river-plus-ocean shipping model.

At0.0625°,182 of194 ports qualify, all51 islands have at least one qualifying port, and all51 terminal distributions are distinct. The12 exclusions comprise7 source harbor-type RN ports,1 OR,1 RT,1 LC and2 CN. These source labels are recorded without claiming that every exclusion is a river: unresolved coastal geometry or source-location error can also matter. Maximum attachment distance is24.482km and median4.079km. The main component has189,550 nodes; no sea connection is added to accommodate a port. This amendment is based only on source geometry and harbor metadata, not on prediction scores.

## Independent conductance layers

AIS is averaged directly from the CRC-verified World Bank/IMF0.01° external overview onto an exactly aligned0.5° grid. The original native-to-overview resampling method is unspecified. AIS values are position-report density, including stationary vessels, not bilateral voyages or cargo. Original NoData is excluded, never replaced by zero. Valid land-zero pixels are part of the source cell average: this is not a water-area-normalized mean.

Wind is the published NOAA NCEP Reanalysis1 mean magnitude of monthly10m U/V vectors over Jan2015–Feb2021 (74 months). Repeat the already derived1° cell values over their0.5° children without claiming additional observational resolution. Wind speed loses transport direction and represents neither butterfly flight-height winds nor individual transport events.

Use the identical observed water domain for three conductance fields:

- Uniform water: conductivity1.
- AIS: conductivity `1 + ln(1 + observed AIS density)`.
- Wind: conductivity equal to the positive mean monthly-vector speed.

An edge has conductance equal to the harmonic mean of its endpoint conductivities, multiplied by the source-resolved shared water-face width in kilometers, divided by the spherical separation of the two0.5° block centers. Four-neighbor finite-volume geometry accounts for latitude and narrow face openings. Each water patch inherits its block's source field. Resistance within an aggregated patch is unresolved; narrow channels within a block are not assigned an independently measured conductance.

## Island terminals and numerical validation

Each island injects one total unit, spread equally over its qualifying original marine-access ports; weights of ports sharing a graph node are added. The182-port allocation is identical for uniform, AIS and wind models. There is no artificial short circuit between ports on the same island. The island-pair measure is electrical dissipation `(w_i-w_j)' L^+ (w_i-w_j)`. It equals ordinary effective resistance for single-node terminals but must be described as a distributed-port circuit measure for islands with several ports.

Require all182 qualifying port attachments in one source-defined marine component, all51 distinct source distributions, and finite, symmetric, nonnegative pair results. Require the full194-port eligibility ledger; never omit failed ports from the denominator. Unused isolated sea basins may be omitted from the numerical factorization. Verify series/parallel circuits and distributed-current fixtures. Independently compare geographically selected real pairs with the official Circuitscape implementation. Report actual versions, error and numerical residuals.

## Held-out comparison

Leave out every island together with its50 pairs; train on the other1,225 pairs. This gives51 dependent folds and2,550 predictions for1,275 unique pairs. Do not use pair-independent p-values or treat folds as independent replications.

The base model includes log1p(original great-circle distance), its square, the sum and absolute difference of log1p(original eligible port counts), and, for the later butterfly target only, earlier butterfly Jaccard. The same terms and support enter every comparison. Port counts address part of the relationship between terminal averaging, observation effort and the number of source ports; they do not eliminate all sampling bias.

Normalize each circuit matrix by its response-independent median positive pair value, then use log1p of the normalized value. Compare:

1. Base versus base+uniform-water circuit.
2. Base+uniform versus additionally AIS.
3. Base+uniform versus additionally wind.
4. Base+uniform versus additionally both AIS and wind.

Evaluate host Jaccard and later butterfly Jaccard; earlier butterfly composition is a temporal comparison, not a valid null expected to have zero association with modern geography. Report MSE, MAE and all51 fold losses. No conductivity transformation is tuned on these outcomes. This is retrospective geographic association and held-out-island prediction, not prospective establishment, causal plant carriage, or proof of historical shipping routes.
