# Monarch: cultivation versus botanical origin as competing hypotheses (source feasibility only)

Date: 2026-10-28. Exploratory continuation of chocho PR #88; frozen GEB main/manuscript are unchanged.

## Why this is new and what is already prior art

Micah Freedman's 2026 preprint (DOI:10.64898/2026.04.06.716508) already annotated 38,856 iNaturalist monarch-larval records and retained 30,320 host-assigned records, including host species, cultivation status, time and locality; its winter (Nov–Feb) finding of 209/348 cells with **only non-native host detections** is not a novel chocho result. The dataset is available through the author's public resources page (https://www.freedman-lab.ca/resources). Its observation label assignment is partly targeted to particular geographies, and the published 100-km × week × host thinning is not a complete standardized sampling effort design.

Our distinct test, **conditional on source access and fields**, will ask whether the *relative composition of observed monarch larvae across season* is more strongly associated with **plant cultivation status** than **botanical native/non-native status**, with local space, calendar year and observer effort matching, before entertaining an abiotic-niche mediation claim.

## Independent hypotheses and source gates

- H-native: recorded larvae on recipient-nonnative host plants disproportionately occur in cold-season periods, holding cultivation status and available observation effort constant.
- H-cultivation: recorded larvae on likely cultivated plants disproportionately occur in cold-season periods for both *native and alien plant species*, holding exact region/year/host taxonomic blocks where possible.
- Distinct downstream vital-rate gate: larval occurrence is **not** egg-to-adult emergence or population growth, and absence of records is **not** larval absence. Actual leaf persistence/senescence is **not directly measured** by an author's horticultural-cultivation annotation.

A valid test requires raw record identity, sufficiently precise spatial and monthly date fields, high/confidence host labels, cultivated yes/no/unlikely/unknown, and biologically grounded native status by receiving region. Public data retrieval must have the author's dataset and code terms respected. If only the paper's 209/348 aggregate exists, refuse to re-label it as new analysis.

Pre-specify an identification matrix before opening outcomes:
  1. at least 20 matched region×year blocks with both cultivated/noncultivated records AND native/nonnative botanical observations if testing their interaction;
  2. at least 10 distinct milkweed taxa in both cultivation categories across seasonal periods, or downgrade species-level cultivation effect to descriptive;
  3. spatial leave-region-out replication and observer/site/source-record thinning; no tiny p-values treating repeated photo records as independent sites;
  4. primary outcome is a conditional **larval host composition** contrast, not breeding success, productivity or range expansion;
  5. native status in tropical Mexico below 24N ambiguous for *Asclepias curassavica*: exclude or status-unknown, do not silently label alien.
  6. do not treat host taxonomy native status as geographically fixed globally; do not copy prior author's native class wholesale outside documented receiver origins.

## Current stage

First acquire and audit **metadata only** from the author's explicitly public Dropbox folder. No raw observation coordinates, observer identifiers or annotated image rows will be committed. Pilot may fail if Dropbox blocks GitHub runners; a failed transport is not proof of insufficient biological samples and is not a result.
