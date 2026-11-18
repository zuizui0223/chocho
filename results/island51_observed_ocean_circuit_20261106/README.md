# Observed marine circuits for the original 51 islands

## Current result and scope

The observed marine analysis is executed for **all 51 original islands and all 1,275 original island pairs**. Coastline geometry, commercial AIS and monthly wind now enter a shared marine graph. The numerical calculation is independently checked against the actual Julia implementation of Circuitscape. This continues the earlier distance-only stepping-stone benchmark; its archived result remains available and is not overwritten.

The main predictive result is a **6.98% reduction in later-butterfly held-out MSE** when AIS and wind circuit measures are added to the uniform-ocean model, with improvement in **40 of 51** island holdouts. A subsequent direct-native-pixel AIS sensitivity gives **7.01%**, also improving **40 of 51** holdouts. Wind alone gives a **5.50% reduction**, whereas AIS alone gives **0.56%** in the primary run and **0.66%** with native averaging. For candidate introduced-host similarity, AIS alone gives **−0.15%** and **−0.11%**, respectively. Thus this run supports additional predictive information for later butterfly composition, especially from the wind-weighted marine field. It does not establish ship-carried plant introductions or directed wind dispersal.

**AIS source-value audit completed:** 638,976 native pixels at 39 fixed geographic locations confirm that the large published integer values exist in the original raster, with scale 1 and offset 0. Independent area-overlap means reproduce all 39 checked output cells within **4.68 × 10⁻⁸ relative error**. No source rescaling or outcome-dependent source correction was applied. The values are reported as the published AIS density index, not inferred voyages or cargo.

## Fixed held-out comparison

Values are percentage reductions in mean squared prediction error, `100 × (MSE_reference − MSE_augmented) / MSE_reference`; positive values are improvements. The reference differs between the first row and subsequent rows, as shown explicitly.

| Comparison | Candidate introduced-host similarity | Later butterfly similarity | Earlier butterfly similarity |
|---|---:|---:|---:|
| Base → base + uniform ocean | +4.091% | +5.082% | −0.032% |
| Base + uniform → additionally AIS | −0.152% | +0.558% | −0.039% |
| Base + uniform → additionally wind | +0.384% | +5.496% | −1.037% |
| Base + uniform → additionally AIS and wind | +0.321% | +6.980% | −1.177% |

The later-butterfly uniform model has MSE **0.00645231** and MAE **0.0571634**. Adding both fields gives MSE **0.00600192** and MAE **0.0541205**. The wind-only model improves 38/51 island holdouts, and the combined model improves 40/51. These counts are descriptive; the folds share training observations and are not independent statistical replicates.

Each fold excludes one island and **all 50 pairs involving that island**, fits on the other **1,225 pairs**, and predicts those 50 held-out pairs. All 51 folds are retained, giving 2,550 predictions for 1,275 unique pairs. No individual pair is in the training data of a fold in which it is tested. The base model contains log1p(original great-circle distance), its square, and the sum and absolute difference of log1p(original port counts). The later-butterfly target also controls earlier butterfly similarity. All models use identical targets, fold assignments and marine terminals. Models use train-only standardization, ordinary least squares, and predictions clipped to the Jaccard range [0, 1]. Circuit predictors are log1p transformed after normalization by the response-independent median positive resistance.

Earlier butterfly similarity is a **temporal comparison**, not a null outcome expected to lack association with modern geography. The difference between its prediction scores and later scores is not a formal interaction test or proof that a transport mechanism changed over time. The frozen comparisons were specified before fitting these new outcome models; there is no independent external validation set for choosing among the four comparisons after inspecting this table.

The later target is the original Jaccard similarity of species first indexed on each island in **2011–2022**, with **2001–2010** records as the earlier comparator. Plant candidates are pooled over **2001–2024**. First observation is not an established arrival or breeding date; these outcome definitions are retained from the [original source lock](../../provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_lock.json).

The numerical source is [model_comparison_summary.csv](model_comparison_summary.csv); [all_island_holdout_losses.csv](all_island_holdout_losses.csv) retains every fold, including adverse results.

### Direct-native AIS sensitivity: same interpretation

The source audit identified spatial subsampling in the published overview, so one additional comparison averaged every native 0.005° pixel directly onto the same 0.5° grid. This source-quality sensitivity was fixed **after viewing the primary results and before fitting the native-grid results**. It changes only AIS input values; it is not an independent validation set and was not used to select a different conductivity transformation or model. The [plan and timing](NATIVE_SOURCE_SENSITIVITY_PLAN.md), [source derivative](../../data/exploratory/island51_observed_ocean_20261106/ais_native_quality_check/README.md), and [complete sensitivity results](native_source_sensitivity/README.md) are preserved separately.

| Added predictors beyond uniform ocean | Host: overview | Host: direct native | Later butterflies: overview | Later butterflies: direct native |
|---|---:|---:|---:|---:|
| AIS | −0.152% | −0.106% | +0.558% | +0.659% |
| Wind | +0.384% | +0.384% | +5.496% | +5.496% |
| AIS + wind | +0.321% | +0.352% | +6.980% | +7.006% |

All geometry-only and wind-only control results, including individual fold losses, are **exactly unchanged**. The native combined later model has MSE **0.00600023**, compared with **0.00600192** in the primary run, and both improve 40/51 holdouts. Its MAE changes slightly in the adverse direction, **0.05412049 → 0.05413185**, while remaining below the uniform model's 0.05716343; the native variant is not superior on every loss measure. Earlier-butterfly combined performance remains adverse, changing from **−1.177%** to **−1.254%**. Native versus overview AIS resistances correlate at **0.999838**, but individual resistances change by up to **5.55%**; stability of the prediction conclusion does not imply identical source fields. All 12 comparisons and 612 fold records are retained in the sensitivity output, alongside the [direct comparison CSV](native_source_sensitivity/primary_vs_native_comparison.csv).

## What now runs

### Authentic original geography

The recovered WPI2019 source has 194 original focal ports assigned to the original 51 island groups. Reconstructing their spherical mean locations reproduces every archived pair distance, with maximum absolute difference **4.58 × 10⁻¹¹ km**. Original biological pair and island tables are checksum verified.

A 1° cell-center ocean mask failed independently of the biological outcomes: the Mediterranean and Baltic were disconnected, several Danish islands shared terminal cells, and Borneo's average port location lay inland. The replacement uses original **GSHHG 2.3.7 intermediate coastlines**, **0.0625°** water-cell centers and four-neighbor connections with longitude wrap. Separate connected water patches inside a **0.5°** conductance block remain separate graph nodes. Fine face contacts determine inter-block connections; no invented sea bridge is added to accommodate an island.

Ports attach independently to their nearest fine cell in the largest source-defined marine component, within the fixed **25 km** limit. Marine access is an operational definition based on this geographic proximity test, not verification of a navigable harbor-to-ocean path; the displacement can cross unresolved land or waterways. **182 of 194** qualify, retaining all 51 islands with distinct terminal distributions. The maximum included attachment distance is **24.482 km** and the median is **4.079 km**. The full marine graph contains **189,550 nodes and 370,996 edges**; island terminals occupy 157 distinct graph nodes.

The 12 excluded original port IDs are `283, 402, 408, 416, 417, 748, 656, 878, 1579, 31, 1152, 1153`. All remain in [original194_port_attachment_audit.csv](original194_port_attachment_audit.csv), with their source harbor types, coordinates and exclusion reasons. Distances are recorded for accepted attachments; excluded rows have no accepted attachment distance. They include inland or river ports and other locations without qualifying marine access at this resolution. This is an explicitly scoped **marine-access** analysis, not a complete navigable river-plus-ocean network. The geometry-only amendment and earlier failed audits are preserved in the [protocol](PROTOCOL.md) and [source bundle](../../data/exploratory/island51_observed_ocean_20261106/README.md).

### Three conductance fields, one graph

- **Uniform:** cell conductivity 1.
- **AIS:** `1 + ln(1 + published commercial AIS index)`, from the World Bank/IMF January 2015–February 2021 source. The CRC-verified 0.01° external overview is averaged directly onto the exactly aligned 0.5° grid. Original NoData is excluded; valid source zero is retained. The publisher's native-to-overview resampling method is unspecified.
- **Wind:** NOAA NCEP/NCAR Reanalysis 1 mean magnitude of 74 monthly 10 m U/V vectors over the same period. The native approximately 1.875° longitude × 1.9° Gaussian-latitude fields were interpolated to the existing derived 1° grid; those values are repeated over their 0.5° children, without claiming finer observational information.

The additional [native-source audit](ais_source_quality/README.md) found that all **159,744 sampled first-overview pixels** equal the top-left pixel of their native 2 × 2 block, with zero mask disagreements; 41,099 differ from the native block mean. This is evidence consistent with nearest-neighbor overview construction at the sampled locations. Across 35 sampled cells with positive native means, the median absolute relative difference in coarse averages is **19.39%**. Consequently the primary 0.5° product is an average of the published **0.01° overview**, not a reconstructed average of every native 0.005° observation. The pre-fit protocol already specified that distinction and is unchanged. A separate [native-source sensitivity plan](NATIVE_SOURCE_SENSITIVITY_PLAN.md) was fixed after seeing the primary results and before fitting native-grid results, to assess this concrete approximation. The published integer scale is retained exactly: the chosen logarithmic conductivity is specific to that scale and was not calibrated against actual transport. Numerical solver agreement cannot establish that physical calibration.

An edge's conductance is the harmonic mean of its endpoint conductivities multiplied by resolved shared water-face width, divided by spherical block-center distance. This accounts for latitude and the number of fine water contacts across a block boundary. The block graph preserves connectivity; it is **not an exact electrical reduction** of the fine water grid. Within-block bottlenecks and travel distances remain unresolved.

Each island injects a total of one unit, equally distributed over its qualifying original ports. Multiple ports at one graph node have their weights added. Ports belonging to the same island are **not electrically short-circuited**. For island pair `i, j`, the reported circuit measure is the energy

`R_ij = (w_i − w_j)' L⁺ (w_i − w_j)`.

This is ordinary effective resistance for single-node terminals, and a distributed-port circuit measure otherwise. The 182-port weights and marine graph support are identical for all three fields. Original 194-port counts, rather than only included counts, enter the statistical controls.

## Numerical and software verification

The complete 51-island matrices use SciPy sparse grounded LU with batched current right-hand sides. Their maximum linear-system relative residual is **2.54 × 10⁻¹³**.

Actual **Julia 1.13.3 / Circuitscape 6.1.0**, CHOLMOD, double precision independently solved two geography-selected real pairs on all three complete marine networks: `9–15` and `715–1577`, for **six checks**. The maximum relative resistance difference from SciPy is **7.99 × 10⁻¹⁰**; the maximum all-node Kirchhoff current residual is **8.81 × 10⁻⁹**. Seven known series, ring and distributed-chain assertions also passed. Full 1,275-pair matrices are SciPy outputs; the official Julia verification scope is these six cases, not every pair.

The [official validation receipt](official_circuitscape_validation.json) pins versions, source and matrix hashes, configuration, current-injection representation, numerical tolerances and actual differences. [Per-case checks](official_circuitscape_validation_cases.csv) give both solvers' values. The distributed current representation uses one negative-injection node as direct ground and retains the other signed injections; the ground supplies the omitted sink by current conservation. This avoids an official advanced-mode zero-total-source skip while preserving the same electrical problem. Full-node current balance and the edge-energy identity independently check this equivalence.

The [separate native-AIS official receipt](native_source_sensitivity/native_AIS_official_validation.json) records actual Julia checks for the same two pairs on the new native-AIS network. Maximum relative resistance difference is **3.97 × 10⁻¹¹**, and maximum Kirchhoff relative residual is **4.37 × 10⁻¹⁰**. The four unchanged uniform/wind controls also passed in that reexecution. This adds two new network/pair cases to the original six; it does not expand the official calculation to every island pair.

The repository suite including the native-source sensitivity reports **189 passed, 1 skipped**. The skipped file requires the absent optional `xarray` dependency for rebuilding the upstream NOAA wind raster; the frozen wind input, its hashes, alignment and its use in this circuit run were checked. The existing sealed repository layout whitelist was retained. No analysis requirements or tests were weakened to obtain a pass.

An [independent geometry and fold review](independent_geometry_and_fold_review.json) reconstructed all 3,093 fine-water components without false merges or splits, all 2,912,808 inter-node fine contacts, all 157 distributed node weights, and every one of the 612 fold records' summary contributions. This review found no blocking issue within that scope.

## Reproduce and inspect

Install the declared environment and execute from the repository root:

```bash
python -m pip install -e ".[test,analysis]"
python -m pytest -q
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python scripts/run_island51_observed_ocean_circuits.py \
  --input-root data/exploratory/island51_observed_ocean_20261106 \
  --coast-dir data/exploratory/island51_observed_ocean_20261106/coast_0p0625 \
  --marine-ports-only \
  --out-dir /tmp/chocho-observed-ocean-results \
  --network-dir /tmp/chocho-observed-ocean-networks
```

The runner verifies the frozen inventory and protocol before accessing biological outcomes. The persisted small input bundle is sufficient for this reproduction; downloading the original large raster archives is necessary only to rebuild the external source derivations. The [source README](../../data/exploratory/island51_observed_ocean_20261106/README.md) records original URLs, transformations and source-license notices. The dedicated [GitHub Actions workflow](../../.github/workflows/chocho-island51-observed-ocean-circuit.yml) executes the same frozen-source analysis and archives its outputs.

The separately executed [official Julia validation](official_validation/README.md) includes the pinned environment, portable six-case verification entry point, and original successful execution evidence. The optional figure can be regenerated after installing `.[figure]`:

```bash
python scripts/render_island51_observed_ocean_summary.py \
  --input-root data/exploratory/island51_observed_ocean_20261106 \
  --result-dir results/island51_observed_ocean_circuit_20261106
```

| Output | Purpose |
|---|---|
| [PROTOCOL.md](PROTOCOL.md) | Pre-fit specification and the source-only geometry amendment |
| [source_preparation.json](source_preparation.json) | Input identities, geometry, port allocation and source checks |
| [island51_observed_circuit_pairs.csv](island51_observed_circuit_pairs.csv) | Original outcomes and all three circuit predictors for all 1,275 pairs |
| [observed_ocean_circuit_holdout.json](observed_ocean_circuit_holdout.json) | Complete model definitions, diagnostics and summaries |
| `*_distributed_port_resistance.npz` | Three 51 × 51 distributed-current circuit matrices |
| [Summary figure](observed_ocean_circuit_summary.png) / [SVG](observed_ocean_circuit_summary.svg) | Source geography and incremental held-out prediction |

## Interpretation boundary and later work

This is retrospective community-composition prediction. AIS position density includes stationary vessels and is not observed bilateral voyages, cargo, plant propagule flow or historical introduction routes. Monthly 10 m wind speed is scalar, loses direction, and is neither flight-altitude wind nor a record of butterfly transport events. Observed later records do not establish true colonization, and pooled plant observations can postdate butterfly records. Geography, climate, recording effort and uneven port coverage remain possible explanations for associations after the declared controls. Coarse-grid land-zero AIS pixels also enter the source means; the means are not normalized by water area.

The plant target represents candidate introduced larval hosts from the original biological source, not confirmed feeding on those plants in every focal island. The present result is a completed source-backed feasibility analysis with actual out-of-island prediction and official numerical checks. A later mechanistic extension would need directional and temporally aligned transport evidence, finer electrical resolution or documented resolution sensitivity, and an explicit river-access treatment. Such extensions are separate from these fixed results. No claim that modern AIS independently explains introduced host composition is supported by this run, and the previous positive host–butterfly and island–mainland results are not invalidated by the modern AIS comparison.
