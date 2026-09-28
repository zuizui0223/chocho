# Reviewer-defense sensitivity protocol v0.1

**Status:** FROZEN_BEFORE_SAME_FAMILY_HOST_NULL_RESULT_OPENED  
**Frozen on:** 2026-09-28  
**Scope:** Post-hoc reviewer-defense analyses for the butterfly specialization manuscript. These analyses respond to explicit concerns about structural concentration metrics, finite geographic headroom, geographic accessibility of climate comparisons, occurrence validation, taxonomic non-independence, and geographic sampling bias. They do not replace the original independent climate-test contract.

## 1. Same-family, same-host-count host null

### Question
Do the observed larval-host portfolios generate anthropogenic geographic resource expansion or host-contribution concentration that differs from what would be expected from host richness and plant-family identity alone?

### Null construction
For each butterfly, preserve the exact number of resolved accepted host species drawn from each WCVP plant family. Within each family, replace the observed hosts without replacement by a deterministic pseudo-random sample from all accepted WCVP-resolved plant species represented in the fixed HOSTS mirror and that family. Recompute native and contemporary WGSRPD3 unions from the fixed WCVP v13 snapshot.

### Fixed inputs
- HOSTS commit: 808e0b869f9ec1adf8efff87cf6a395adda103e0
- rWCVPdata/WCVP commit: 65bed76bae9d644ccb6ad200c05f9f5071d89e05
- S1 resource descriptors: SHA-256 894f48dbca1760fc4fa75bfee8f663540ab4b9380f8b9daf2bc09440e8bb0cdc
- WGSRPD level 3 geography as used in the original reconstruction
- 999 deterministic iterations for the reviewer-defense run

### Endpoints
For each eligible butterfly:
1. native WGSRPD3 resource breadth;
2. contemporary WGSRPD3 resource breadth;
3. introduced-added WGSRPD3 units;
4. log proportional expansion;
5. effective contributor number among hosts producing introduced-added units;
6. maximum single-host fractional contribution.

Primary reviewer-defense summaries are observed-minus-null-median residuals and the fraction of species above the corresponding null median. Host-family-breadth associations in these residuals are secondary.

### Claim boundary
The null pool is all WCVP-resolved plant species represented in HOSTS, not all vascular plant species. It asks whether the documented butterfly host set is unusual relative to alternative documented HOSTS plants from the same plant families and with the same family-specific host counts. It does not model host evolution, availability, preference, or physiological feasibility.

## 2. Structural concentration diagnostic

For the host-taxonomy-adequate expanded subset, quantify:
- Spearman(host-family breadth, contributing-host species count);
- Spearman(contributing-host count, effective contributor number);
- Spearman(contributing-host count, maximum single-host share);
- partial Spearman(host-family breadth, each concentration metric | contributing-host count).

The raw architecture contrast is not treated as an independent ecological effect if it largely disappears after conditioning on contributor count.

## 3. Geographic-headroom sensitivity

Report both absolute introduced-added WGSRPD3 units and proportional expansion. Diagnose finite headroom using the official total number of WGSRPD3 units and:
- partial Spearman(host-family breadth, log proportional expansion | native resource breadth);
- fraction of available global headroom filled = added / (total WGSRPD3 units - native units).

The headroom fraction is a sensitivity quantity, not a replacement primary endpoint.

## 4. Contemporary-occurrence validation

Using the independently acquired 2010-2026 butterfly GBIF panel, count observed butterfly species x WGSRPD3 units outside the native host envelope and ask what fraction fall inside the contemporary envelope only after introduced host distributions are retained. Report this for:
1. the full 32-species independent panel;
2. the 24 species passing the frozen pre-climate quality gate.

This validates opportunity geographically; it does not establish that the introduced host was used at the occurrence locality.

## 5. Same-WGSRPD-level-1 climate sensitivity

Recompute the climate-filtering comparison only for held-out observed versus never-observed resource units belonging to the same WGSRPD level-1 region. A species is informative only if at least one eligible cross-set pair remains. Summarize the score distribution and rerun the frozen host-breadth test with contemporary resource breadth as control.

This is a coarse accessibility sensitivity, not a mechanistic dispersal model.

## 6. Butterfly-family blocking

Use the LepTraits butterfly Family field to remove between-family mean differences before rank associations. Report family-blocked sensitivities for:
- host-family breadth versus proportional expansion, controlling native resource breadth;
- host-family breadth versus architecture metrics, controlling contributing-host count.

This is a taxonomic-clustering sensitivity, not a substitute for a dated phylogenetic covariance model.

## 7. Coarse geographic robustness

Assign each resource-eligible species to the WGSRPD level-1 region containing the largest share of its native reconstructed host-resource units. Report host-family breadth versus proportional expansion separately in sufficiently populated regions. This diagnoses whether the pooled pattern is dominated by one broad biogeographic region; it is not a correction for HOSTS/LepTraits sampling completeness.

## 8. Independent climate-test precision

For the original n=24 independent partial Spearman result, report an approximate Fisher-z 95% interval and the approximate one-sided alpha=0.05 correlation magnitude required for 80% power. These are descriptive precision diagnostics; the original permutation p-value remains the inferential result.
