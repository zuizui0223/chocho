# CIRAD 2022 local port N gradient: independent offline falsification (2026-11-17)

**Status: exploratory only; preserve positive global resource and butterfly results in PR #157.** This is a new analysis of the source-verified inputs already used by nutrient pilot PR #122. No changes to the canonical manuscript and no claim that port sewage N reaches terrestrial host plants.

## Question and source-locked methods

Compare modeled 2022 agricultural N within **0–10 km** of each indexed WPI2019 port with the **10–25 km ring** around that same port. The 3-km CIRAD human-excretion and cropland-input rasters (original units 10 kg N/ha/yr) are *MODELED INPUT POTENTIAL*, not measured local soil N. The outer ring mean was independently reconstructed from the 10/25 km valid-pixel **counts and means**; masked sea/no-data pixels were never filled with zero.

Same WPI port-ID crosswalk among three original independent Actions artifacts, each 1,596 ports:
- PANGAEA Jones 2015 wastewater/source port CSV: [run 41563940689](https://github.com/zuizui0223/chocho/actions/runs/41563940689), exact SHA256 `9e3d35516927741f1f176a029b9d41bb03c70c850c066d15e5fe9730c2036482`
- ESA WorldCover2021 source port CSV: [run 41562997326](https://github.com/zuizui0223/chocho/actions/runs/41562997326), SHA256 `1cb891076e13aa9d93bced1b40e6471d447771af2b8b028f134ed052d13bf071`
- CIRAD 2022 25/10 km grids: [run 41923811745](https://github.com/zuizui0223/chocho/actions/runs/41923811745), SHA256 `7d934dec9bf3d62ceb1b7583ca5d84a54c57b0440147913ff27d9895026c3c0c`

All 1,596 port IDs, fields shared by independent Jones/WorldCover plant/WPI sources, and locked raw CSV SHA256 were checked. Deliberate byte tampering aborts.

Responses are (A) candidate-host indexed log((alien+0.5)/(native+0.5)) with indexed-host count effort proxy, and (B) recorded alien share alien/(alien+native) **only where denominator >0**, fractional-logit one equal-weight port. Not actual plant abundance or butterfly feeding. Country FE and country-cluster 95% intervals. Covariates: log modeled 2015 wastewater, log WWTP25 km count, 2021 500-m built/crop fractions, latitude and longitude harmonics, log both 2022 human/crop N 25-km mean. Exposure: `ΔcropN = log1p(mean0–10km) - log1p(mean10–25km)`; human N equivalent separately. Each raster must have >=8 valid pixels and >=10% valid coverage **both** near and ring; countries >=5 qualifying ports and predictor variation.

## Findings: positive catalogue-level association, not demonstrated N fertilization

| Definition | ΔcropN beta (95% country cluster) | Matched sample |
|---|---|---|
| Smoothed alien/native **indexed ratio**, country FE OLS joint crop/human gradient | **+0.0862 [+0.0238,+0.1486]** | 1,013 ports / 40 countries |
| **Recorded alien fraction**, one-port fractional logit, same joint gradients | **+0.1014 [+0.0377,+0.1652]** | 843 ports / 36 countries |
| Alien indexed **richness** (log1p), with both gradients | −0.0158 [−0.0947,+0.0630] | 1,013 / 40 |
| Native indexed **richness** (log1p), with both gradients | −0.0897 [−0.1849,+0.0056] | 1,013 / 40 |

**Hard source-grid falsification:**
1. Requiring >=50% valid pixel coverage in **each** near/ring for human and agricultural N reduces crop-gradient alien fraction to **+0.0491 [−0.0285,+0.1267]** (579 ports/28 countries) and log-ratio to +0.0814 [−0.0016,+0.1645] (654 ports/31 countries).
2. Reconstructing the difference in **proportion of nonzero modeled-N grid pixels** between 0–10 and 10–25 km from the original count-weighted zero fractions: crop-gradient alien fraction becomes **+0.1205 [−0.0041,+0.2450]** in the full 843-port sample and **+0.0173 [−0.0685,+0.1032]** with >=50% pixels. For smoothed ratio on strict 654 ports, concentration-gradient beta +0.0199 [−0.0675,+0.1074], but nonzero area gradient beta **+0.4660 [+0.0343,+0.8977]**. This could reflect crop landuse/source-model occupancy rather than nutrient-driven biological filtering.
3. Limiting to >=5 **nonzero input raster pixels in each ring** and fitting nonzero-cell mean intensity difference plus nonzero-area terms: alien fraction intensity delta **+0.0980 [−0.0391,+0.2352]** (608 ports/29 countries); no clear independent gradient intensity signature.
4. Leave-one-country-out **alien indexed fraction** prediction (without unavailable country FE): add crop gradient -> country-equal MSE reduction **+2.7437%** (36 country-bootstrap descriptive 95% [−0.0015%,+5.1723%], 21/36 held countries improved; n=843). With >=50% grid pixels, improvement **+0.0782%** ([−2.6738%,+2.8748%], 16/28 countries improved; n=579). Human gradient alone is −0.7897% (full), and joint human+crop +1.9687% with interval spanning 0. These results are neither temporal forecasting nor butterfly colonization.

## Scientific decision

**Retain the exploratory observational positive `ΔcropN`/indexed relative candidate-host association**, which was previously untested at this geographic scale. **Do not upgrade it to local sewage nutrient enrichment, agriculture-induced alien-host proliferation, or butterfly establishment**. The signal fails strict land-coverage, nonzero raster-cell composition and held-country transfer reliability gates. The nonzero N-cell share may be an agrarian land-cover/recording opportunity proxy. Unmatched dates (WPI2019, wastewater2015, WorldCover2021, crop/human N2022, pooled GBIF2001–2024) and unmatched plant 5 km vs N 0–10 and 10–25 km windows preclude physical causal interpretation. Repeated exploratory tests make p-values descriptive only.

**Provenance:** Scripts `chocho_n_gradient_pilot.py`, `chocho_n_gradient_transfer.py`, `chocho_n_gradient_cell_control.py` and corresponding full JSON receipts were run from original Actions CSVs, without public API re-queries; the standalone source-locked code and JSON can be attached to coauthor-review / supplemental exploratory materials. **No third-party source CSV redistributed in Git.** Keep PR Draft; do not merge this into the canonical GEB manuscript or claim a new independent N mechanism.