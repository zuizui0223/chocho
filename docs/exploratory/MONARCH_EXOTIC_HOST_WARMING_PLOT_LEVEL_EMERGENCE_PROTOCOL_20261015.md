# An introduced milkweed can yield adults under one temperature and fewer under another — source-verified emergence calibration

**2026-10-15 | POST-HOC REANALYSIS OF ALREADY-PUBLISHED HOST×TEMPERATURE EXPERIMENT.** The paper's host-by-warming finding was known before choosing this analysis, and its eight treatment-cell adult counts were inspected before writing this file. It is NOT registered or a new independent experiment. The separate GEB manuscript PR #38 / main are unchanged.

## Scientific motivation and limits

The global 239-butterfly chocho reconstruction maps **known-host botanical presence**. It does not observe production of adults. An independently conducted published experiment with a plant introduced to the region demonstrates the ecological filter at the relevant fitness endpoint: **host plant species × temperature (with pathogen exposure held in the randomized design)**.

Original article: Ragonese et al. (2025), *Ecological Entomology*, DOI [10.1111/een.70010](https://doi.org/10.1111/een.70010).
Original source author repository pinned to commit `74e3e2cd8079401c9ccd4359f9e39e2d52d1f7d2`, path `MWwarming_comp_May16.csv`, Git blob SHA-1 `e7a625911d8dc9b871a52d7657314a5d05472216`.
Site: Athens, Georgia, original **2021 experiment**, 30 plot blocks / half ambient and half heated, **four plant units per plot**, **two monarch larvae per host plant**. Total 240 assigned original larvae, of which 239 have known adult status and one has an unknown fate. Fixed experimental factors: milkweed species `tropical` (*Asclepias curassavica*, cultivated and not native in Georgia) versus `swamp` (*A. incarnata*, native); temperature `ambient` or `elevated`; OE inoculation `infected` versus `control`.

The original paper already reported that host identity interacts with warming in some development/survival contrasts and parasite-tolerance results. A possible posthoc source verification **cannot rebrand this as a new discovered interaction**.

The study does not randomize native-vs-exotic origin as a replicated multi-species factor: origin and *species* are fully confounded. It neither measures host occurrence in the chocho's 82 introduced-only global monarch resource cells nor proves butterflies are wild-colonizing any of them. The study followed larvae to **adult eclosion**, not flight capability, mating or reproductive offspring.

## Primary descriptive question

**At identical intended stage/parasite composition, does the *plant species* contrast in adult emergence change between ambient and warmed experimental plots?**

Treat **30 experimental plots** as the independent higher-level temperature assignment units. Within each plot, each host has **two distinct plant units**, one infected and one control, **two originally assigned larvae per plant**. Thus the plot-level marginal adult probability per plant species is `success/4`, *not* individual larva independence. Compare hosts within each plot, then compare the 15 paired host differences in warmed plots with the 15 paired host differences in ambient plots.

For a plot `j`, `d_j = adult_fraction(tropical)_j − adult_fraction(swamp)_j`.

```text
Delta = mean(d_j for elevated plots) - mean(d_j for ambient plots)
```

**Do NOT** select plots/host plants based on diagnosed infection status (post-randomization) or on adult survival/plant condition. Parasite inoculation assignment remains equally represented within each plot host.

### Ambiguous assigned larval fate

Source originally assigned 240 larvae but adult outcome missing for one. **Primary lower-bound treatment coding** counts unknown as no *observed* adult; a **separate upper-bound sensitivity** codes it as potentially producing one adult. This single missing case is in `ambient × tropical × control`. Both are recorded; neither is a claim that the larva biologically died/survived. Do not silently drop one larva.

Check source pupation/eclosion `Surv_pupa` vs `Surv_adult` in all 239 known records, and retain any pupation-without-eclosion events rather than imputing adulthood. This describes **stage endpoint counts**, not natural mortality causes.

### Uncertainty and diagnostics

Prospectively *for this posthoc reanalysis* (after knowing the original article):
- **Stratified plot bootstrap**, fixed 15 plot units per temperature, replacement **of whole plots** (never larvae/individual plants). 9,999 draws, deterministic seed 20261015. Calculate percentile 95% interval for `Delta`, and paired host contrast at ambient/elevated separately. Bootstrap may be unstable with only 15 plot units in each temp and does not fix physical plot spatial dependence or treatment-temperature layout.
- **Missing outcome extreme sensitivity**: rerun for unknown as failed vs one possible successful adult; do not condition on a known-only subset as primary.
- **Leave-one-plot-out**: remove each source plot, recompute interaction, and report range; this is sensitivity, not an independent plot-varying treatment effect test.
- **Within-pupa-to-adult audit**: explicitly report pooled and eight-cell original endpoint counts, checking all larvae initially assigned.
- No per-larva independent-binomial significance, no source-published model claims as new, no automatic conclusion of an ecological trap, no extrapolation to native vs introduced host plants generally.

### Biological finding threshold and stop rule

A large point difference is informative **as the original paper's published context dependence**, but not enough for a novel butterfly ecology manuscript. If plot-cluster 95% interval includes zero, label the exact interaction uncertain despite source-published posthoc pairwise tests. If it excludes zero, still label it as secondary descriptive reproduction of prior art and avoid native-exotic generalized statements.

The value of this reanalysis to chocho is **calibrating why geography cannot be treated as adult-emergence-equivalent host resources**: one known and accessible host species can have differing larval-to-adult success under different heat/pathogen regimes. Native host `A. incarnata` and nonnative `A. curassavica` are specific *species* in a particular trial. Combined independent species and site replication, locally matched planted/wild occurrence, egg cohort fate and adult reproductive offspring would be needed for broader ecological causal inference.

All code should fail closed on source blob identity, exact 240 original individual records, plot/plant nesting, 8 balanced 30-larva cells, and non-finite or unexpected source success fields.
