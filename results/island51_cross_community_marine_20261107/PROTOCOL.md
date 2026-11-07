# Frozen three-community marine connectivity validation on the same source-eligible island pairs

Date 2026-11-07. The previous PR138 public GBIF complete 194-port historical GENERAL Plantae speciesKey census and its positive +1.4740% later butterfly held-out MSE increment from AIS beyond recorded 2001–2010 general flora were known before this protocol. Earlier butterfly AIS effects and original local 750-km graph-label placebo limitations were also known. This is POST-HOC EXPLORATORY cross-community predictive comparison on the SAME original source, not independent preregistration or an independent source population.

## Question
When matching the original historical 2001–2010 indexed general plant communities, early (2001–2010) indexed butterfly communities, and later (2011–2022) first-indexed butterfly communities on **the exact same island pairs and frozen geographical train/test folds**, does the original marine AIS/wind information improve prediction similarly across taxa/time windows? Does AIS remain informative after adjusting contemporary/potential larval host information for later butterflies? This is comparison of historical source *indexed general flora*, not historical larval-host identity, and AIS is a modern 2015–2021 layer; NO causally ordered plant–marine pathway is inferred.

## Population and source-only inclusion
Use original 51 landmass IDs, all 194 original WPI botanical ports, 182 qualified marine ports, original butterfly source pair table with 1,275 unordered pairs and the source-pinned complete PR138 GBIF Plantae speciesKey snapshot:
- results/island51_pre2011_public_general_flora_20261107/source/GBIF_original194_port_2001_2010_speciesKey.csv.gz, sha256 42bca05d45e97dd38f1b879fc0934828fc9b7451a2b7848a29eac11c12e23061.
- its archived original194 WPI/GBIF count/source provenance and source receipt, locked by PR138.
- original fixed AIS, wind and uniform complete Circuitscape matrices and longitude-block-pair folds.

An island belongs to the main 2001–2010 general-flora comparison only if at least one Plantae speciesKey was indexed at any of its source WPI port-buffer 5-km searches, yielding exactly **45 source islands**, 6 source islands with zero speciesKey (NOT biological plant absence). The common population is all **45 choose 2 = 990 unordered island pairs** for which BOTH original endpoints have an indexed early plant speciesKey. This is a SOURCE-ONLY, TARGET-INDEPENDENT filter, never a butterfly-outcome-selected sample. Save an original-pair-index manifest and original 51-island coverage explicitly, including the six omitted no-index islands. Do not claim the 990 pair cohort is representative of all 51 islands or independent whole-island census.

## Baseline model, covariate rules
The **common geography + botanical recording opportunity** baseline for all comparable models includes:
- logged verified original WPI port great-circle distance AND its square;
- sum and absolute difference of log1p original per-island WPI biological port counts;
- original uniform ocean Circuitscape resistance (geographic barrier control);
- the four pinned source latitude/hemisphere pair controls (signed latitude sum, absolute latitude sum, absolute latitude difference, opposite hemisphere indicator);
- the two original 2001–2010 GBIF *general plant record-count* recording-opportunity pair controls from PR131;
- the two original PR138 2001–2010 indexed GENERAL plant speciesKey marginal richness pair controls (sum and abs difference log1p indexed speciesKey richness);
- two original GBIF requery 2001–2010 summed general plant-port-record count pair controls (sum and absolute difference log1p), to expose current GBIF index opportunity heterogeneity independently of old archive counts.

**Historical general flora** response: Jaccard of the two 2001–2010 Plantae speciesKey sets (defined because both sets are nonempty in the fixed 990-pair cohort), without plant Jaccard leaking into its model covariates, and without any butterfly-derived covariates.
**Earlier butterfly** response: original 2001–2010 butterfly Jaccard; no predictor derived from later butterflies or retrospective candidate HOSTS (2001–2024).
**Later butterfly** response: original 2011–2022 butterfly Jaccard, with original earlier 2001–2010 butterfly Jaccard included as a fixed baseline covariate (as in PR138). This later panel does not include the original pooled candidate-host (2001–2024) network.
**Later butterfly with retrospective HOSTS controls**, clearly secondary: same later response as above but additionally include the original candidate-host Jaccard and original marginal host richness pair controls, which can contain records that POSTDATE later butterfly indexing. It is a sensitivity panel, not time-correct independent validation.

No model contains historical general flora Jaccard as a covariate when predicting itself. All new outcome-independent original source predictors and explicit panel feature lists must be saved and audited.

## Fixed four marine models and comparisons per panel
Models BASE (no native AIS/wind), AIS (BASE + native AIS), WIND (BASE + original wind-speed resistance) and BOTH (BASE + both native AIS and wind), with no AIS × wind interaction terms. For each of the four prespecified panels, contrast BASE→AIS, BASE→WIND, AIS→BOTH (wind beyond AIS), WIND→BOTH (AIS beyond wind), BASE→BOTH. Evaluate in all **15 already frozen original longitude block-pair splits**: train and test only on the fixed source-only **990** eligible pairs, always with BOTH original longitude-block endpoints absent from training for their test split. The resulting 990 test pairs must each occur exactly once per model and panel, and no training pair contains any test endpoint. OLS with intercept, predictors standardized from that fold's training rows only, clip predicted Jaccard to [0,1]. Absolutely no outcome-driven bin, threshold, source filtering or feature tuning.

Expected: 4 panels × 4 model designs × 990 source-pair predictions = **15,840 model-specific predictions**, 4 × 4 × 15 = **240 model-fold MSE/MAE records**, 4 × 5 × 15 = **300 fold contrast records**, and 4 × 5 × 4 = **80 scope comparisons** (all indexed 990; within longitude block; between longitude blocks; and source-only ≥10 original 2001–2010 indexed general-plant speciesKeys at BOTH endpoints). Missing source keys remain excluded only from *this botanical-comparable cohort*; all original 51-island PR138 results remain unchanged.

## Verification and scientific boundaries
- No network requery; use only saved original PR138 source snapshot and original frozen GitHub files.
- Verify historical 194 port-source identities and speciesKey union from their raw snapshot; exactly 45 islands indexed; exactly 990 source pairs. Verify each plant Jaccard and original earlier/later butterfly observation from original source pairs.
- Preserve each source-pair original index, response, panel, fold, prediction and squared/absolute loss. An independent Python standard-library auditor reconstructs original 15 longitude blocks, all 990 plant-eligible original source pairs, all 15,840 saved prediction losses, 240 folds, 300 contrast folds and 80 summary rows; verify no training-test island endpoint leakage. Source-pinned ecological scope ensures no post-hoc selective results.
- Retain all favorable and adverse MSE/MAE signs. No naive pairwise independent p-values. The historical general flora Jaccard values are *indexed port-buffer records*, not real whole-island plant abundances, host acceptance, vegetation similarity in 2001–2010, or proof of biological introduction. AIS is a modern 2015–2021 vessel-density-derived undirected maritime resistance; earlier plants/butterflies predate it. The network's predictive information does not establish actual ship itineraries, cargo, butterflies visiting ports, effective dispersal or colonization. Indexed butterfly first records are not arrivals or breeding. Original cross-taxon differences can reflect nonuniform sampling coverage, temporal mismatch and calibration; no assumption of historical process change is allowed.
