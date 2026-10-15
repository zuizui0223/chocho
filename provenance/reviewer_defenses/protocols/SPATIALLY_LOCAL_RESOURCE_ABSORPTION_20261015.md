# Spatially local plant-range null — post-hoc geographic sensitivity (2026-10-15)

**Status:** exploratory specification frozen before inspecting this new spatial null's outcomes. The climate-matched absorption excess from PR #43 was already observed when this design was chosen. The GEB submission and its main figures are not modified.

## Ecological question and primary contrast

Previously, the 239-butterfly reconstruction showed **46,387 of 75,346 introduced host × focal butterfly × WGSRPD3 incidences** fall where that butterfly already has a *native* larval host opportunity. With plant-family, continent, and coarse CHELSA climate-bin structure preserved, excess native-resource absorption remained +288 to +558 incidences (climate-null v0.1).

**Question:** Does absorption remain above constrained null expectation when plants can exchange introduced regions only with *geographically proximate places* resembling each other's broad climate? This is not a test of whether plants have moved or whether larvae successfully reproduce.

Primary statistic: absolute **native-absorbed candidate incidences** `A`; divide by the fixed total `H=75,346` for effect-size communication. Secondary: unique new butterfly×region opportunity `U`, applying `H=A+N`, `U=N-R`. Single-host share, exact shared-host geography and co-user counts are descriptive only.

## Fixed inputs and spatial representation

- Identical HOSTS/WCVP commits, S1 239-species manifest and 365 WGSRPD3 regions as in the completed climate null.
- Reuse the previously archived CHELSA V2.1 BIO1, BIO7, BIO12, BIO15 regional medians; use **BIO1 annual mean temperature and BIO12 annual precipitation 2×2 median bins within WGSRPD Level1** as the climate condition.
- WGSRPD Level3 representative points from the official `tdwg/wgsrpd` repository at commit `52da7828aba9d461dd133c27b3bd7a4407161f54`, file `geojson/level3.geojson`, `shapely.geometry.shape(...).representative_point()`. Use shortest great-circle distance. This represents polygons very coarsely, particularly multi-island territories or high-area regions.
- Fail on duplicate codes or discrepancies with the pinned `data/frozen/wgsrpd3_level1_lookup_v0.1.json`. Do not impute missing region coordinates.

## Three spatial restrictions, chosen before new response

Apply plant-family × WGSRPD Level1 × BIO1/BIO12 2×2 stratified double-edge swaps, *plus* both conditions:

1. The two regions exchanging introduced plants must have representative-point great-circle distance <= the specified caliper.
2. Each plant's new candidate introduced region must lie within the same caliper of **at least one of that plant's own original introduced regions**. This original-footprint tether prevents unrestricted diffusion by repeated local swaps.

- **Primary caliper: 1,500 km** (comparatively broad because WGSRPD3 units are coarse).
- Stricter sensitivity: **750 km**.
- Looser sensitivity: **3,000 km**.

These kilometer thresholds are a spatial-sensitivity *design choice*, not biologically identified butterfly flight or plant dispersal radii. The null holds fixed the original host identities, introduced host counts per plant within each family×Level1×climate group, introduced plant-family totals per destination, native geographic exclusions and every unchanged immobile edge.

A proposed swap is accepted only if it is valid in both directions and respects all geographic, climatic and native-exclusion constraints. One chain is not evidence of exchangeability; the algorithm's support can be disconnected under these constraints.

## Predeclared admissibility and sampling

- Geographic and climate data must cover **all 365** resource regions. No environmental/centroid imputation.
- At least **40%** of original introduced plant×region edges must be part of a potential same-family×continent×climate swap group.
- At least **5%** of all introduced incidences must belong to host-region edges with some *geographically valid* alternative two-edge swap at baseline (precomputed).
- Burn-in: 5× total E proposal attempts; **>=100** accepted swaps and **>=2%** of original introduced edge incidences must be changed. These are geometry/evaluability screens, not convergence guarantees.
- Two independently initialized seeds **20261015** and **20261016**; **199 sampled worlds** each, spaced by 1× total E proposals.
- Report accepted-swaps per sample, burn mobility, lag-one correlation, effective sample size approximation only if applicable, and empirical distribution quantiles. Never pool MCMC draws as exchangeable independent permutations without calibration.
- If any condition fails, write `NOT_EVALUABLE`; do **not** claim that the ecological residual disappeared.

## Decision mapping

A reproducible residual across evaluable local-caliper chains strengthens the descriptive finding of non-random **potential native-resource overlap** beyond the named structural controls. An attenuated or reversed residual indicates that geography explains some/all of the earlier result under this null; it cannot identify a causal dispersal process. Inconsistent/non-mixing chains require abstention, not tuned replacement calipers.

**Strict claim ceiling:** No local host use, larval survival, phenology, ecological competition, range colonization, native plant abundance, elevation, habitat continuity or introduction-era history is measured. Geographic position is proxied by one point per broad WGSRPD3 polygon. Even a stable result remains post-hoc and conditional on fixed HOSTS/WCVP inventory and unmodelled reporting bias.

Original climate-null run: https://github.com/zuizui0223/chocho/actions/runs/38696007445
