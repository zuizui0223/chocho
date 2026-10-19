# Species × region fixed-margin null for climatic resource alignment

**Date: 19 October 2026. EXPLORATORY, explicitly POST-HOC to the first uniformly assigned within-continent benchmark, not a preregistered confirmatory ecological p-value.** No literature case studies; original 239 butterfly×365 region macro WCVP/HOSTS/CHELSA source only.

## Specific biological falsification
Original host-range reconstruction found that 9,451 of 14,481 (65.3%) new butterfly×WGSRPD3 host-opportunity cells in climate-covered species lie within the central 5–95% annual-temperature/precipitation rectangle of their own native host-resource regions, versus 57.0% under a uniform-within-Level1 geometric expectation. A more demanding question is:

> **Is this climate-analog alignment more than the outcome of which regions globally gain the most butterfly resources, how many resource regions each butterfly gains, and which continent supplies those gains?**

If not, demote climate filtering to geometric sorting by region and introduction intensity; do not claim a new butterfly-specific alignment.

## Frozen estimand and source
Use exact source ZIP from run 39006554609, SHA256 `408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323`, without redistributing raw sidecars. Reconstruct exactly 239 butterfly species, 365 regions, 14,553 original added opportunity cells and their native exclusions. All 239 species are included in the null graph even if one of the five lacks sufficient native climate; the test statistic includes only the source-evaluable 234.

For each butterfly, compute native-host climate-analogue set using fixed CHELSA BIO1/BIO12 percentile rectangle at **5–95%** (primary) and **10–90%** (sensitivity) as in the already audited source-relative climate test. Count `O_q`, the number of *added* butterfly×region incidences where region belongs to the focal butterfly's analogue set. Main observed `O_05=9451`, sensitivity `O_10=6948`. No climate cutoffs or taxa chosen after test outcome.

## Strict structural null: within-Level1 bipartite swaps
- Represent the added-resource incidence matrix `A[butterfly,WGSRPD3_region]` with 14,553 1s.
- For each swap proposal choose TWO *observed added-incidence edges* within the SAME WGSRPD Level1 biogeographic area and **different butterflies and regions**: `(b1,r1),(b2,r2)` become `(b1,r2),(b2,r1)` only if both crossed cells are currently unoccupied AND neither is part of that butterfly's native-host resource envelope.
- Preserve (i) each butterfly's added regions **by continent**, (ii) each region's added-butterfly total, (iii) the whole 239×365 original native resource exclusion matrix, and (iv) all non-resource regions and climatic vectors. This also fixes butterfly×region total added incidence and introduced-host-resource opportunity volume. It DOES NOT preserve identities of introduced host plants or local larval recruitment and therefore is a **butterfly×region structural** null, not a plant-family climate-matched causal null.
- Use exactly **two predeclared independent chains** with seeds `20261019` and `20261020`, each a burn-in of `10×14,553` pair proposals, followed by **99 draws** separated by `2×14,553` proposals. Accept/reject rules deterministic from fixed seed, with accepted proposals and fraction of edge identities changed recorded. Retain all original region/species margins at every audited draw.
- Compute for each 99-draw chain a median, observed-minus-null absolute and share, lag-1 autocorrelation, approximate ESS, and between-chain median agreement. Quality gates: post-burn incidence identity turnover ≥0.1, valid source-matched denominators, all native exclusions, ≥50 accepted swaps per sampling interval on average, approximate ESS≥30 for both temperature-band statistics, and between-chain medians within 2× larger individual chain SD. If gates fail, result remains NOT EVALUABLE and tail proportions must not be promoted.
- The tail fraction has resolution 1/100 for a single 99-draw correlated chain and is **not an exact independent permutation p-value** even when simulation diagnostics pass. Do not merge two chains into 198 independent test replicates.

## Interpretation rule
- If `O_05` does not exceed structurally matched null medians consistently or 10–90% sensitivity reverses, the earlier +8.3-point climate-analog excess is mainly region intensity/geographic sorting, not unique butterfly-host climate alignment.
- If an excess remains, report its **effect size relative to the full 14,481 evaluated added cells**, not only a nominal Monte Carlo tail. Avoid claims of butterfly abundance, actual range expansion, climate-adapted larvae, successful immigration, or introduced plant benefits.
- The historical WCVP `Introduced` indicator is not dated; no species-level independent true butterfly climate niche, plant biomass, local observation-effort adjustment, site proximity or shared-host dependence is measured.

## Reproducibility stop
Only derived aggregate source-conditioned JSON/figures should be distributed while TDWG/WCVP distribution permissions are pending. Keep all code and outputs exclusively on global macro exploratory branch and PR #46; neither main nor frozen submission PR #38 is to be changed.
