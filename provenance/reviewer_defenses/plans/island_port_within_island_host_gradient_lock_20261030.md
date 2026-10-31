# Island-port ecological gateway: within-island harbor-gradient diagnostic (2026-10-30)

**Primary research question was superseded:** [maritime alien-host rescue and stepping-stone hypothesis](maritime_host_rescue_steppingstone_20261030.md). This archived urbanization-gradient specification is **secondary descriptive support only**, not a test of ship-borne seed transport, nonharbor enrichment, butterfly establishment or stepping-stone movement. The earlier fixed specification and any exploratory outputs remain intact.


**Stage:** focused exploratory test of the first gateway prerequisite. The full harbor-versus-nonharbor host-plant comparison is **not possible yet** with existing published snapshots, because only harbor-centered botanical inventories are available. Do not confuse this proxy test with harbor-specific plant establishment, transport or subsequent butterfly invasion.

## Source and pre-outcome feasibility (audited before the first new gradient coefficient)

Pinned PR #104 baseline commit: `5bc347e3160f2a376377e2098598e8b3bb040a2e`.

1. `provenance/reviewer_defenses/results/coastal_port_vegetation_20261030/coastal_control_vegetation_pairs.csv`, Git blob SHA `10360b39e20a720000c5ac8f9f156073c67e7c9b`. Contains 625 harbor × island pairs at 1 km (the larger 5 km ring is a **separate support set**), with ESA WorldCover 2021 built near/far fractions.
2. `provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv`, Git blob SHA `fdf73b11ddcb1644f46730123350ae52cc61cae6`. Includes complete 2001–2024 harbor-5-km botanical host catalogue count, by current *regional catalogue* classification.
3. The join is exact existing `port_id`; harbor island `landmass_id` comes only from the vegetation source. No spatial reallocation and no name-based matching.

**Support check, computed before scores:** 625 1-km vegetation ports; 210 harbor queries complete, 97 islands; 161 harbors / 77 islands with classified host richness >= 1. Requiring >= 2 classified ports per island and nonconstant built-near-minus-far gives **24 islands, 104 harbors**. These are neither all harbors nor a balanced global sample. Zero classified records are **missing fraction**, not alien-host absence.

## Locked variables and statistical quantities

Primary response, each port: `alien_fraction = contemporary_not_native_in_catalogue / (contemporary_not_native_in_catalogue + native_in_catalogue)`. These are numbers of host **species** in a compiled 2001–2024 search area, not introduced individuals; ambiguous/no-region/unknown statuses are excluded from the denominator but retained in auditing. An observed zero alien numerator among positive denominators is a **catalogue observation** only. These fractions do not establish historic introduction, true absence or local larval use.

Primary exposure: `built_difference = built_near - built_far` for the 1-km harbor vs same-island far-shore WorldCover sampling. No radius selection after results. It measures 2021 landscape contrast, not an exogenous experimental port intervention.

Primary coefficient: within-island linear regression with each island given total weight 1 (each harbor weight 1/number of eligible harbors in its island), demeaning response, exposure and covariate by each island's **weighted within-island** mean. Adjust `log(1 + native + alien classified host species)` as a proxy for botanical record/catalogue completeness. Report the slope as **percentage points of observed alien-catalogue fraction per +10 percentage points of built contrast**. This is a descriptive linear association; bounded response and duplicate host taxa across ports mean the species counts are not independent trials. Do NOT use binomial trial SE.

Secondary prespecified checks, explicitly not alternative primaries: (a) same regression without the completeness covariate, (b) use `built_near` instead of `built_difference`, (c) leave one island out (minimum and maximum estimates, sign consistency), (d) island-cluster bootstrap 2,000 resamples (seed 20261030); percentile 2.5–97.5% interval and number of successful full-rank resamples, (e) repeat primary on the intersection with the previously eligible 51-island host-network cohort **without changing the full primary sample**. Show the sample sizes for every subset. Also show raw within-island covariance and any pathological rank/denominator failures.

Fail closed if key uniqueness, pinned blob identities, sample-support counts, raw range [0,1], or reported vegetation differences fail. If predictor and effort proxy are collinear, report *unidentified* coefficient rather than a spurious value. This pre-specified finite cohort is observational and exploratory; formal `p` values are not needed, and no threshold-dependent claim is permitted.

## Explicitly untested gateway link

The 40 islands with joint harbor-vegetation and older harbor-derived food-plant-network coverage do **not** have the matched **nonharbor botanical inventories** required to identify harbor-specific food plant enrichment. Other missing inputs: contemporary botanical effort/providership matching, a nonharbor urban coast control, independently known plant introduction years, and dated port shipping traffic aligned with butterfly events. The present gradient check cannot assign port-related resource concentration relative to nonports, nor say anything causal about butterfly colonization. 

**Next true endpoint (not estimated here):** for fixed islands and matched harbor/nonharbor urban/natural coast points, obtain plant species-by-year lists and comparable effort in equal-area/equal-season windows; compare *host-specific enrichment versus nonhost plant controls*, then test if genuinely **earlier** resources help predict later butterfly detections after species, source-distance, geography and observation adjustments. Retain 225 butterfly candidate species and missing states.

**Publication boundary:** Keep exploratory continuation in draft PR #104; do not change the GEB manuscript or the PR #102 cohort.
