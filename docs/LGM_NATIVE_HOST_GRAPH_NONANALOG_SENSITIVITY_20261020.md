# 21 ka WCVP-native plant climatic graph: source-replayed nonanalog uncertainty audit (2026-10-20)

## Evidence status

**Source-pinned exploratory macroecology.** This audit is **not** ancient plant occupancy reconstruction, an observed biogeographic dispersal corridor, an insect's genetic break, or a confirmed larval-host-specific historical mechanism. The TTF [v0.5 negative host-identity control](https://github.com/zuizui0223/TTF/pull/149) is preserved unchanged.

## Actual finished work

The original, trained-on-0BP, 15-term quadratic (4 CHELSA variable) logistic function was reconstructed **from its frozen 0BP 6,337-location predicted probabilities only**, by weighted inverse-logit linear regression. All original 416 native WCVP plant models were reproduced without refitting the underlying modern presence/background data. Saved original 21 ka probabilities were reserved exclusively for verifying the recovered function on a separate climate set.

- **416/416** original models passed.
- Maximum absolute replay error on original 0 BP scores: **4.70945×10⁻⁸**.
- Maximum absolute replay error on 21 ka validation scores: **7.99087×10⁻⁸**.
- Disagreements between recovered logistic model and original v0.2 supported 0.5° main masks at the source suitability decision threshold: **0 cells** at 0 BP and **0 cells** at 21 ka.
- Original main 0.5° 416-plant climate masks and thresholds are **unaltered**. New extrapolation scenario includes previously excluded *finite* climate values outside each modern plant's observed native-range environmental training domain. It does not impute missing CHELSA values.

Source-verification run: [Actions #39387916509](https://github.com/zuizui0223/chocho/actions/runs/39387916509) **SUCCESS**. Frozen source-input main packed mask SHA256 `169b9da084d0043abba5c89f008137d6867b554dfb7c7dbe283b750142d220af`; restored extrapolation mask SHA256 `67c0d99852b8b8896559dbc0da07504bcb07a6ad7c017afa0eb86e31d40acc2a`. Both provided with machine readable model-level tables and receipts.

The previously observed `Paullinia pinnata` (`2543867`) source-retraining discrepancy of **1.205×10⁻⁴** was not treated by relaxing the fixed fit-replay error, retuning model parameters or dropping the species. Rather, inverse-logit recovery of the **saved original 0BP response function** independently restored its original 21 ka prediction, and all 416 plants were verified.

## What remains of the directional ecological pattern?

| Conditional potential-climate graph metric | Original analog-only | Include extrapolated climate |
|---|---:|---:|
| WCVP native plant models with smaller 21 ka potentially suitable area | **289/416 (69.5%)** | **243/416 (58.4%)** |
| Median paired per-plant potentially suitable area delta, km² | **−2,903,383** | **−1,233,515** |
| Plants with lower largest-component *area fraction* | **271/416** | **253/416** |
| Insects with all known recorded host SDMs completed | **513** | **513** |
| Insects with smaller 21ka potentially suitable host area | **399/513** | **322/513** |
| Insects with fewer 4-neighbor components | **477/513** | **457/513** |
| Insects with lower largest-component *area fraction* (4 neighbors) | **371/513 (72.3%)** | **326/513 (63.5%)** |
| Insects with lower largest-component *area fraction* (8 neighbors) | **379/513** | **320/513** |
| Median paired within-insect 4-neighbor largest-component *fraction* delta | **−0.08582** | **−0.04128** |
| Distinct accepted-host combinations | **349** | **349** |
| Distinct accepted-host combinations with lower largest-component *area fraction* (4 neighbors) | **252/349** | **226/349** |

The registered source-sampled 0.5° climate study and PMIP4 21 ka terrain reconstruction were frozen before this sensitivity step. The component count is not by itself a measure of better connectivity: with a smaller suitable area, some disconnected components disappear entirely, even while the largest component becomes a *smaller share* of remaining modelled suitable habitat.

Among the 513 insects, the original analog-only 4-neighbor largest-share decline is 371. With hypothetical extrapolation, 305 of those still show decline, **66 switch to no decline**, and 21 of the original no-decline cases switch to decline, giving **326 declines**. For the *area* direction, 82 of the 399 analog-only declines change to nondecline, while five original no-decline taxa now decline, yielding 322 declines.

### Source and dependence

The original source taxonomy panel contained 641 insect species and 1,037 insect–recorded host links (1–3 hosts per insect). Only the 513 taxa with model predictions for **every recorded host** are scored, and 349 **distinct host combinations** exist among those 513 insects (363 distinct source plant host models). Hence the modelled resource landscape replicates only **349 different host combinations**, not 513 independent host landscapes, and even the 349 share plant models. Genus/host-combination dependence needs explicit modelling before organismal inference.

See full [Actions #39388770609](https://github.com/zuizui0223/chocho/actions/runs/39388770609) **SUCCESS**, source-verifying all original 513 species 4- and 8-neighbor area/component/fraction metrics while recomputing matching inferred signals for the same IDs. Full per-species results and archived SHA receipt: `results/lgm_513_insect_nonanalog_climate_graph_sensitivity_v04/`.

## Inferential boundary / next mechanistic research dependency

The reported landscapes are **0.5° potential climatic suitability on PMIP4 reconstructed ice-free land**, created with classifiers trained on 2010–2026 native-region screened plant observations. These are not populations or presence/absence maps in the LGM; no ancestral occupancy, host dispersal speed, true insect resource supply or insect genetic breaks are known. The source modern climate envelope excludes about **49%** of the physical 21 ka grid for a median plant model, and the extrapolated predictions in those novel climates are not empirically calibrated. Thus neither analog censoring nor unrestricted extrapolation can be called the one true historical landscape.

**Defensible ecological finding:** A majority of the *modelled* taxon-specific climatic resource landscapes show less available area and reduced large-component dominance in 21 ka under both explicit alternative climate extrapolation conventions, but the magnitude and participating taxa are substantially uncertainty-dependent.

**Still unproven:** The central species-specific historical larval-host *occupancy and true connectivity* hypothesis or any genetic spatial boundary. A next biology design needs independent paleobotanical occupancy constraints where available, regionally appropriate environmental calibration, taxon-specific dispersal uncertainty, narrow-coastline / climate-resolution sensitivity, and a truly new source-disjoint genetic response test. No reopening of the negative TTF v0.5 result.
