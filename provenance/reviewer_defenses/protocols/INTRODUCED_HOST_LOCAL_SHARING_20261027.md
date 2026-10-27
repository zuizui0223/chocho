# Frozen source-limited ecological test: Do relocated plants become multi-consumer butterfly host hubs?

**2026-10-27. Frozen BEFORE computing the focal butterfly×host plant×HOSTS Location shared-host results.** This is an independent *HOSTS publication-catalogue network* macroanalysis on the fixed 239 butterfly species from chocho, stacked in a draft branch beyond PR #79. Do not alter the GEB submission, main, or earlier stage-evidence records.

## Ecological question

The relocated larval-host resource in chocho is an opportunity for a butterfly; it might also become a shared resource for several butterfly species, allowing altered interaction networks. The Chile field paper documents introduced *Pieris brassicae* and an indigenous *Tatochila* caterpillar on the same *Tropaeolum majus*, but **single-paper co-occurrence is not measured competition**. We ask a global, data-backed *source-documentation* question first:

**Within the SAME explicit original HOSTS Location, are WCVP-introduced larval host plants reported for multiple distinct butterflies more often than WCVP-native host plants after controlling for recorded plant-family composition, locality, and host-report sampling intensity?**

This is not proof of sympatry/larval resource competition: HOSTS Location often is a country or large subregion; its records are from unknown dates and source literature.

## Fixed source and structural units

1. The original 239 butterflies in the SHA-verified author-only WCVP/HOSTS archive from original GEB project (inner source ZIP SHA256 `408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323`, 11 member hashes); globalbioticinteractions/HOSTS `resource.csv` pinned commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`, 140,485 rows, 86,826 rows with full species, 4,695 original 239-butterfly rows and 3,898 exact frozen accepted butterfly×host pair original rows.
2. Original accepted plant IDs and *host plant family* from `inputs/original/insect_host_accepted.csv`; no fuzzy synonym reconciliation. Exact botanical origin according to PR #76: plant distribution WCVP native, introduced, WCVP-unlisted and source ambiguous from original 969 country×WGSRPD3 disjoint atoms (commit `4104c784ddcf97318c4a1bb8bd949b3147628494`), original ≥95% mapped country area and ≥95% uniform plant-status dominance (80% status-only sensitivity). Explicit Hawaii→HAW and West Malaysia→MLY where source country intersection is verifiable. Keep all broad Nearctic/Oriental/unknown/locality-unsafe origins as excluded but counted in denominators.
3. Raw source row grain `butterfly species × exact accepted WCVP host plant ID × original Location × original HOSTS source row`. Deduplicate *edges* to `butterfly species × plant ID × Location`, expected **3,586** accepted edges. Main plant-node grain is `accepted plant ID × exact original Location`. Its recorded **local consumer degree** is number of distinct original 239 butterflies; shared host if ≥2 distinct butterfly species, **cross-genus shared** if consumers include ≥2 insect genera. Report **plant node counts** and **edge counts** separately. Raw HOSTS row count per plant-node is a publication/reporting intensity proxy and *not biological abundance*.
4. Primary WCVP statuses of source nodes: `NATIVE_CONFIRMED` and `INTRODUCED_CONFIRMED`, after source-guard. `WCVP_UNLISTED` serves as independent negative *documentation* control (not confirmed cultivated-only). Mixed/unscorable remain explicit.

## Predeclared source-matched references

**Primary descriptive estimand**: among resolved plant×Location nodes, proportion with ≥2 documented butterfly species for introduced versus native plant origins. Also source local *edge-weighted* fraction of unique butterfly×host×Location reports located on such shared nodes. Report source reporting intensity in both origins before biological interpretation. No test based on raw HOSTS records as independent butterflies.

**Source-label matched exact artificial null**: For each `original Location × accepted plant family` stratum with both *native* and *introduced* resolved plant nodes, condition on `N` total nodes, `K` nodes shared by multiple butterflies, `m` introduced plant nodes; random label exchange yields `X~Hypergeom(N,K,m)`. Convolve exact distributions and report observed introduced shared node count, fixed-margin mean/interval/tails, and **number of genuinely permutable strata and source locations**. Discard (but count) strata with only native or only introduced plants from the *matched contrast*; they cannot identify an introduced-specific effect.

**Publication-intensity sensitivity**: additionally match on original source `raw HOSTS rows per plant×Location` binned as **1**, **2**, **3–5**, **≥6**. Hold source families and original location fixed as above. If no exchangeable cells remain, information gate fails. The reporter count is an endogenous source detection proxy, not the true abundance of either host or butterfly.

**Plant-area-unlisted control**: repeat the identical within-location/family and reporting-bin references using `WCVP_UNLISTED` vs native, as a way to distinguish introduced status from generic source geographic documentation. This control is imperfect because unlisted may include cultivated/undocumented and plant geography conflicts.

**Sensitivity**: redo the entire procedure with original precommitted WCVP 80% geographic-area status dominance; source taxonomic cohort and original 257 location strings unchanged. No posthoc reweighting that promotes broad labels to countries.

## Information gate

Treat an artificial matched contrast as **not ecologically interpretable** if fewer than 10 distinct *permutable source locations*, fewer than 15 permutable location×family strata, fewer than 20 source introduced host nodes in matched strata, or fewer than 10 total shared host nodes in matched strata. The structural result can be descriptive even if the gate fails. A convolved mathematical source tail is NOT a biological p-value: host introduction is not randomly assigned and butterfly preferences, plant habitat and research effort co-vary.

**No assertion of actual competition, displacement, co-occurring larvae, wild completion or improved insect fitness** unless independently local time-matched larval observational data test it. Most published HOSTS localities are geopolitical areas, not surveyed plots, and report dates are not standardized. This panel has 239 selected butterflies, not a global census, and consumer native/alien distribution is not mapped.

## Output and publication firewall

Export only aggregate counts and location-level/source-family source QA **without accepted butterfly×plant pair identifiers**, original copyrighted HOSTS rows, observers or exact original localities of point samples. Region labels may be shared only as aggregated locality counts/summaries; no original interaction matrix. Preserve full original 3586 host triples and 179 introduced triples as checksum-controlled numeric margins. Source scientific claims remain documentation co-reporting, not direct ecological interactions.

Results are a separate candidate follow-up; original GEB paper and previously completed PR #75–#79 unaffected.
