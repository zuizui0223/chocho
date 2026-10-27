# Pre-outcome contract — cultivation rather than plant origin and observed monarch larval season

2026-10-28. This is an **exploratory observational reanalysis** of the publicly shared original dataset used by Freedman (2026), DOI:10.64898/2026.04.06.716508. The author's 209/348 nonnative-only winter cells and 30,320 core annotations **are prior art**. No GEB manuscript/main edits, no population or individual survival claim.

Original public ZIP SHA256 `4e4def63669b34738bcbfb8ba32ccf6c61d02cc6992b8f1ebfa8bc350a3987db`; source `v1/v1_thinned_larval_host_records.csv` is exactly 12,942 author-thinned rows. External URL from https://www.freedman-lab.ca/resources. Independent source label audit (no season inspected): cultivated categories: yes 1,815, likely 2,835, maybe 2,602, unlikely 2,487, no 3,199, missing four; documented source native 10,494, nonnative 2,448; source rearing=no 12,765. Source license requires credit; do not publish or ingest author observer IDs and raw locality data into chocho Git.

## Frozen primary statistical question

Among *source-native host plant taxa*, do **observed monarch larvae** associated with clearly cultivated plants (cultivated=`yes`) occur disproportionately in **Nov–Feb** relative to larvae associated with clearly uncultivated plants (cultivated=`no`), *within the same exact milkweed host species, source state/province and calendar year*?

- Source unit: author's already spatial-week-host-thinned individual photo record, used only once. Do not treat this as sampled host plant absence / population replication.
- Filter: host confidence `certain/high/medium` (already author core), source `rearing=no`, cultivated strictly `yes/no`, source `native_status=native`, non-ambiguous E/W field (exclude `na`), source setting not `obscured`, host exact two-word scientific binomial, parse valid `month` and `year`, `hex_id`, `State`.
- Binary outcome: cold season Nov–Feb vs other months Mar–Oct. A Mantel–Haenszel common **recording composition** odds ratio is assembled from 2×2 cells in host×State×year strata with both cultivated and wild observations and both seasonal values. Zero cross-product denominator => insufficient information, never pseudo-count.
- **Information gate frozen before seasonal outcomes:** at least **20** jointly-informative host×State×year strata across **10** exact native host species and **10** state/province units; ≥15 winter-cultivated and ≥15 winter-uncultivated accepted photo records in informative strata; no one state dominating ≥50% of informative records. The reported effect stays **exploratory** even if this gate passes.
- Site-/state-level bootstrap (999 resamples, seed 20261028) produces a descriptive interval only when ≥950 finite draws; not an ecological causal CI.
- Precision-controlled sensitivity: exact host×100km hex×year stratification under identical record filter; do not widen this if it fails.
- Label robustness: `yes/likely` vs `no/unlikely` (exclude maybe), same stratum and season outcome.
- Observation robustness: high/certain host-confidence only, and four-month shoulder seasons (March–April, September–October) against Nov–Feb; no rewired DOY bin search.

## Separate native-status comparator (NOT equivalent to species-matched cultivation)

In source uncultivated (`no`) observations, compare nonnative vs native larval **host records** in winter vs rest within State×year strata. Do not claim this controls host species, because nativeness generally cannot vary within a fixed host×region. Do not interpret unobserved species as a nonhost, plant absence, developmental failure or observed rescue. Both effects may be unknown/collinear.

## Frozen claim classification

Only a *relative seasonal composition of host-labelled monarch larval photos* is estimable. Source annotation partly selected 'interesting' geographies, iNaturalist effort is unequal, cultivation labels do not prove green leaf persistence, and there are **no parallel surveyed host stems lacking larvae** or egg-to-adult survival records. No evidence can establish new climate tolerance, demographic winter growth, competitive displacement, or that cultivating native plants is beneficial. A result in the same direction as Freedman (2026) must be labelled **postpublication secondary contrast**, not original discovery. Missing information gate => report descriptive sample sizes and no mechanism estimate.
