# Frozen exploratory AIS latitude-band vertex relabeling diagnostic

Date: 2026-11-07. Previous PR135 positive AIS prediction increments and adverse wind results were known: AIS-only +1.0020%, +1.3953%, +0.8609%, +1.2959% later MSE improvements in full and port-heavy training exclusions. This is post-hoc exploratory falsification, NOT an independent preregistration or new validation population.

Question: Is the REAL allocation of original native-AIS maritime circuit resistance to island pairs more predictive than full-graph vertex-relabelings among nearby same-hemisphere latitudes? Preserve every source pair, original butterfly target, baseline covariate, and test fold. Entire AIS vertex identities are permuted, NOT individual pairwise values, so the symmetric circuit graph's edge distribution and topological structure are preserved up to relabeling. The true island-to-connectivity link is disrupted. This is a structural placebo, not a mechanistic or causal null and not an exact distance-conditioned randomization test.

## Fixed source
Original 51 eligible landmass IDs, all 1,275 unordered source pairs, original 225 butterfly candidates, 194 botanical ports, 182 marine-qualified ports, original native AIS Circuitscape matrix, original host Jaccard pooled 2001–2024, candidate-host richness, original port counts and distances, pre-2011 all-botanical recording-effort controls, four PR134 latitude/hemisphere controls, uniform ocean and original mean wind resistance. The original target-independent 15 longitude-block double-unseen folds retain 235 within and 1,040 between longitude-block test pairs; BOTH endpoints are absent from their training set. The two butterfly community targets are later 2011–2022 (primary) and earlier 2001–2010 (context).

## Fixed signed-latitude permutation groups
Read the original 51 representative port latitudes from the source-locked island geography with Git blob f0c5b1139e54b4872d28682f0991f89f2368ab43. Verify exact 10 south, 41 north, 0 equatorial islands. Sort separately by signed latitude, breaking ties with numeric landmass ID. Partition in ten *contiguous* bands: south [5,5], north [6,5,5,5,5,5,5,5]. A vertex may only exchange its original AIS label with another source island within its band. Original signed hemisphere preserved. The bands are arbitrary, not natural bioregions; do not retune them based on scores.

## Fixed trial sequence
Seed 20261107; NumPy PCG64 default_rng; exactly 199 UNIQUE nonidentity AIS vertex relabelings, all bands shuffled independently for each candidate using fixed source band order. Reject only identity or duplicate maps BEFORE fitting any response. Replicate 0 is original identity, 1..199 null. Save all 200*51 source island-to-donor island assignments and exact band IDs. Use the SAME null maps across both target periods and both comparisons; no new AIS circuit simulation.

## Fixed two AIS-addition contrasts
1. HRE_LAT to HRE_AIS_LAT: AIS added on top of original host, richness, historic all-plant recording effort, geography, longitude and latitude controls, but with wind absent.
2. HRE_WIND_LAT to HCRE_LAT: AIS added conditional on the original wind-speed circuit.

For each null, replace ONLY the AIS circuit feature in the augmented model; the reference model is identical for all replicates. Refit the augmented OLS model separately within every original both-endpoints-unseen fold, train-only z-standardization, intercept, [0,1] prediction clipping. Original AIS replicate 0 must reproduce every parent PR135 FULL comparison among both targets and three scopes (12 comparisons, MSE and MAE to 1e-9), or FAIL.

Evaluate exactly 1,275 unique pairs per target/contrast/replicate, no outcome-selected subsets. Compare observed pair-weighted MSE and MAE improvements to the distribution across 199 perturbed matrices, in all-pair, within-longitude and between-longitude classes. Record the null 2.5%, median and 97.5% percentiles and count of null MSE/MAE gains >= observed, and descriptive rank fraction (1+exceedances)/200; this is NOT a formal independent-pair p-value.

Expected outputs: 199 null + 1 original assignment times 51 island rows = 10,200 maps; 2 targets × 2 contrasts × 1,275 baseline predictions = 5,100; 200 maps × 2 targets × 2 contrasts × 1,275 augmented held-out predictions = 1,020,000 records stored as compressed CSV; 200 × 2 × 2 × 3 scope comparison summaries = 2,400 rows; 12 observed-vs-null diagnostic rows. Save frozen seed, null maps and complete raw predictions, rather than only final averages.

An independent standard-library auditor must verify original source pair identities and both target observations, signed-latitude band assignments, all 200 vertex maps and their graph-label constraints, no heldout source island leakage, each raw error, all scope summary rows, and all observed-vs-null ranks. No naive independent-pair or fold confidence interval.

## Inference limits
The placebo destroys original per-island AIS strength and can create geographically implausible links even within broad latitude bands. It preserves the undirected AIS resistance graph only under labels; it is not a generative shipping-route null, not an exact conditional test, and does not establish actual vessel carriage or causal ecological mediation. AIS is indexed vessel-density proxy, not tracked bilateral cargo movements, while wind is undirected mean speed. Historical host species-by-port source is unavailable: the original pooled 2001–2024 host Jaccard can postdate later butterflies. First indexed butterfly observations are NOT confirmed arrivals or colonization. This reuses the SAME selected 51 islands and is exploratory.
