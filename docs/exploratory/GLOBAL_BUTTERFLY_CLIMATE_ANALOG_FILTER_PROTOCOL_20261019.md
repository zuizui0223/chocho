# Global climatic analogue test of redistributed butterfly host resources
**Exploratory protocol fixed 19 October 2026 BEFORE climate-overlap response was evaluated.** This is a new original macro-geography test on the existing pinned 239-butterfly × 365-WGSRPD3 host-resource panel, not a review of individual papers.

## Why this is the ecological next test
The preceding **trait-level model** suggested that disturbance-associated butterflies have greater proportions of the globally unused host-region map filled by introduced hosts; however, the simplistic denominator `365−native_host_regions` includes regions that are climatically implausible for the focal butterfly. A better ecological test must distinguish added host geography that matches the **climatic geography of the butterfly's existing native hosts** from opportunities that are climatic outliers. Those envelopes reflect available host plants, NOT the insect's directly measured thermal physiology.

**Primary question:** For each butterfly, do newly supplied plant-host regions occur disproportionately within temperature-and-precipitation conditions analogous to its native-host region envelope, compared with the same number of opportunities reallocated among native-unsupported regions *within the same broad WGSRPD Level1 area*?

**Secondary hypothesis:** Within climate-analog added opportunity, do disturbance-affiliated butterflies retain a stronger residual opportunity than disturbance-avoidant butterflies after host richness, native resource breadth and butterfly family adjustment? Do not promote this subgroup analysis over the primary host-resource climate comparison.

## Fixed data and computational chain
- Read exact original source ZIP from already verified author-only [run 39006554609](https://github.com/zuizui0223/chocho/actions/runs/39006554609), artifact `butterfly-author-anonymous-offline-source-replication-v01`, contained inner archive SHA256 `408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323`. This includes the independently reconstructed WCVP/HOSTS plant × WGSRPD3 native/current ranges, all known 239 butterfly–host accepted-ID edges, 365 region mapping, and CHELSA v2.1 1981–2010 unit climate.
- **No public redistribution of source archive**: licence review for TDWG/Kew-derived geometry and WCVP/HOSTS sources remains pending. CI downloads and uses locally. Upload **derived aggregate JSON only**, without individual copyrighted raw plant-region records.
- Cross-check every species' native, contemporary and added units against committed source panel (26,530 / 41,083 / 14,553 aggregate; 239 butterflies).
- Use CHELSA `bio1` mean annual temperature and `bio12` annual precipitation only. No response-driven variable substitutions or added selection. All region climate values must match source archive manifest and original climate checks. Keep 365 region classes fixed.

## Climate-analog definition and estimand
For each butterfly with at least ten native-host WGSRPD3 regions with complete climate:
1. Compute 5th–95th empirical percentile interval for `bio1` and independently for `bio12` across native-host resource regions. A candidate region is climatically analogous only if **both** values lie in these intervals. The shape is a simple rectangle, not a fitted SDM or a realized butterfly climate niche.
2. Define candidate space as all source climate-covered WGSRPD3 regions **not already native-host-supported** by the butterfly. Let `K` be the count of candidate regions climate analogous, `M` the total candidate regions, `A` the number of introduced-host added regions.
3. Count `O`: added butterfly × region opportunities inside analog climate.
4. Build a **within-Level1 geographic-composition benchmark**: within each WGSRPD Level1 area, hold `a` newly added regions fixed and take the expected analog count `a × (eligible climate analogous not-native regions)/(all climate-covered not-native regions)`. Sum across Level1 areas to `E`. This is an analytically evaluated stratified hypergeometric mean, NOT an unobserved climate-tracking butterfly record.
5. Estimate per-butterfly `(O−E)/A` only for `A>0`; separately report the fraction of `O/A`, the expected fraction `E/A`, and counts of species with zero analog candidates. Report both species-equal-weighted and added-region-weighted effects. Bootstrap genera (not each region incidence) to quantify between-species pattern uncertainty; genus resampling does NOT correct dependence from shared host plants or spatial autocorrelation.
6. **Predeclared threshold sensitivity**: repeat with 10th–90th percentile bounds, not by searching cutoffs; show eligible n and direction. No fitted butterfly SDM, larval quality, or population effect.

## Biological confounding and stopping gates
- Conditioning on native host geography does not prove butterflies can survive within the climatic rectangle. Native-range climate percentiles can overgeneralize common plants and exclude biologically viable climates.
- The within-Level1 geometric reference fixes continent-level introduced opportunity counts, **not** climate filtering of globally introduced *host plant species*, introduction history, local land cover, natural barriers, observer effort, or the number and size of individual plant ranges. Do not present a strong `O−E` as a causal trait or as a plant-network-matched null.
- Area of WGSRPD3 cells is unequal; cell counts are not land area or abundance. Some unvisited/non-host regions are not colonizable even if climate analog; don't overinterpret the remaining candidate denominator.
- Preserve 239 species even if some are ineligible for climate analysis (report excluded and why; do NOT code A=0 for missing climate). For the disturbance comparison, use the previously frozen strong definitions and predeclared all-macro covariates (native size quadratic, host richness, host family richness, wing proxy, voltinism, taxonomy family FE, HOSTS adequacy), genus-cluster interval.
- A primary climate analogue alignment result only measures **potential redistributed larval resource geographic overlap with climatic analogs**, not butterfly range changes, breeding, resource quality, habitat selection or human benefits.
- If evidence collapses on 10th–90th sensitivity, incomplete climates or source integrity, report NOT ROBUST/NOT EVALUABLE; no alternative posthoc climate axis search.

## Reproduction checkpoint
Reconstruct `model` with original pinned source archive's `build_model`, restoring `region_codes` in the same ascending order as that code, and verify every species' source-native and contemporary counts against the original `anthropogenic_species_metrics.csv`. The source archive is ephemeral GitHub Actions review material; any durable independent publication must rebuild source distributions from pinned WCVP/HOSTS and seek rights clearance.
