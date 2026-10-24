# Spatial redundancy of indexed candidate host plants in 26 GBIF butterfly country cases

**24 October 2026 — explicitly POSTHOC audit, defined after viewing the
existing 26-case source and its coarse GADM1 host-record overlap.**
Not preregistered and not a new experiment or abundance model.

## Fixed records and experimental unit

Use only these archived no-new-network source inputs:
- `gbif_all26_host_gadm_source_20261023.json` (24 strict WCVP
  candidate introduced host taxa; 44 country×plant facets; 26
  original butterfly×country cells, 16 butterflies, 11 countries).
- `gbif_positive48_gadm_lifestage_facets_20261020.json`
  (source-exact 48 positive butterfly-country cells, including the 26).
- `gbif_juvenile_gadm1_conditional_20261021.json` and
  `gbif_source_introduced_hosts_five_cases_20261021.json`
  (the original five stage-labelled added-host cases, 24 source
  immature-stage labels).
Do not refit or change WCVP/HOSTS/BCE source host identities or the
229 butterfly country model. Require exact 26/44/48/5 source
cohort checks before interpreting the derived counts.

The two biological quantities that must be distinguished are:
(1) how many distinct introduced host species appear on the source
catalogue for a butterfly in a country (plant **opportunity inventory**);
(2) how many of these species have >=1 **GBIF indexed observation**
in the same GADM level-1 administrative division as an indexed butterfly
record (local **recorded co-coverage**). Neither implies a wild,
accessible plant, observed feeding or successful butterfly recruitment.

## Fixed descriptive sensitivity design

Count for each original focal butterfly record in its GADM1 division
whether the original accepted candidate-host plant list has **0, 1,
or >=2 distinct indexed plant species** there (do not count repeated
observations of a single plant as separate resources). Repeat for
the 106 original butterfly-positive GADM1 units, and among
the **seven butterfly-country cases with >=2 source candidate plants**
separately from the 19 cases with exactly one listed host.

Vary only the minimum **per-candidate-plant, within-GADM1 GBIF record
count** across 1, 2, 3, 5; this is a source-detection robustness check,
not ecological absence or known plant abundance. Do not cherry-pick
the threshold. Record largest potential loss of co-covered *GBIF
butterfly records* upon deleting **one** candidate plant's indexed
administrative-unit footprint, for each of the seven multi-candidate
cases. This is an index-level single-source fragility scenario and
**not** a simulated botanical extinction.

For the previously fixed **five stage-labelled positive added-host cases**,
report the same 0/1/>=2 categorical plant-GBIF coverage for the **24
source immature-stage labels**, preserving their 11-cell/39-label parent
provenance. No host identification is inferred from stage labels.

Audit fixed original integer totals: 26 source cases, 453 mapped
butterfly records, 231 with >=1 candidate-plant GADM record, 106
butterfly-positive GADM units and 40 with >=1 candidate-plant
co-record; 44 plant×country source queries (3,595 country-level
plant occurrence records). Every result must reproduce original
country case counts and preserve the observed-source denominator.
No raw occurrence records or author-only geospatial geometries.

## Inference boundary

This is an exact source-index network **measurement** result, not a
sampling-based null or independent evidence of ecological
specialization, realized functional redundancy, host suitability,
larval survival, establishment, dispersal or plant population loss.
A 0 means `no plant GBIF record under fixed filters`, never
biological host absence. Country-area attribution and GADM1 unit
names reflect reporting systems that have different biases and may
include cultivation. Report local taxon source insufficiency and
unequal GBIF plant-record effort prominently. Published raw-host
availability and a precise same-site host identification remain
necessary for the conservation mechanism.

All new derived files are excluded from the **frozen GEB v0.2
anonymous review snapshot**; the original GEB manuscript is unchanged.
