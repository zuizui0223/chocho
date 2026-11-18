# Port-indexed butterfly host plants: nutrient affinity beyond worldwide occurrence prevalence

2026-11-02 | [source-locked GH Actions run](https://github.com/zuizui0223/chocho/actions/runs/41403516863) | **exploratory evidence, not manuscript promotion**

## Finding

**The positive plant nutrient-affinity ↔ port-recording association is not solely a byproduct of worldwide GBIF country occupancy or global record volume.** It remains positive after ranking both global reporting proxies out of the N–port association. This gives a quantitative step beyond the original +0.341 Spearman observation. It does not establish a direct nutrient or sewage mechanism.

| N score × number of indexed WPI ports | Spearman or partial Spearman | Descriptive genus-bootstrap 95% interval |
|---|---:|---:|
| Raw, 182 accepted GBIF species | **+0.338489** | **[+0.19684, +0.46920]** |
| Adjust for 2001–2024 global GBIF record count | **+0.334322** | **[+0.19593, +0.46440]** |
| Adjust for number of GBIF-indexed countries | **+0.213587** | **[+0.06771, +0.34860]** |
| **Adjust for both worldwide reporting proxies** | **+0.209710** | **[+0.06490, +0.34839]** |

Input population: 687 candidate non-native butterfly larval-host plant species recorded in 1,596 ports (5km GBIF botanical searches, 2001–2024). Exact German Ellenberg nutrient indicator N available for 183 (26.6%), then 182 strict exact accepted GBIF species-key matches, representing 126 genera. `Salix alba` was the one N-scored species not resolved under the strict accepted-species condition. Median global GBIF indexed records among the 182 species: 105,644; median GBIF-indexed countries: 56. Source IDs and code in the machine-readable [receipt](GLOBAL_GBIF_PREVALENCE_CONTROL_20261102.json).

## Reproducible operations

- Uses original port-plant species frequency table, original frozen code, and same pinned Germany Ellenberg plant indicator table, [commit `73b3239`](https://github.com/bathyscapher/EVA/commit/73b323983702a74cb9b6867bcebe00a307ea6566).
- Re-matched plant names to GBIF with exact accepted-species keys. For each taxon, queried GBIF occurrence/search in 2001–2024 with geographical coordinates and without geospatial issues; collected global `count` and `country` facets. Historical index date cannot be reconstructed from this modern query alone.
- Spearman uses midranks. Partial Spearman computes residualized ranks for N and port count given globally indexed occurrence count and/or number of countries. Genus-resampling `n=2999` preserves same-genus species grouping; it is a descriptive uncertainty screen rather than a standalone independent causal test.
- Executable script: `scripts/exploratory_port_nutrient_global_prevalence.py` and offline regression checks `tests/test_exploratory_port_nutrient_global_prevalence.py`. Workflow `.github/workflows/chocho-port-nutrient-global-prevalence.yml`. GH Actions source-receipt artifact includes 182 plant-level inputs and exact mapped GBIF species keys; it has a 30-day artifact-retention window, while this narrative and receipt are committed as durable summary.

## What can and cannot be concluded

A positive correlation after country coverage adjustment is not an independent experimental comparison. Global GBIF occurrence count includes a subset of the same data from ports (part–whole overlap). Countries differ in size, sampling effort, degree of urbanization, plant origin and botanical collection histories. Both global covariates are imperfect statistical controls and a nonzero partial rank is not proof of shipping or enrichment.

The hypothesis that ports locally accumulate these candidate host plants **through nutrient loading** has still not passed the separate independent botanical check. The [strict 78-port / 35-island non-overlapping 5km within-island study](STRICT_WPI78_SOURCE_AUDIT_20261102.md) found mean N port minus nonport +0.113 early (9 independent eligible islands; island-bootstrap CI −0.034 to +0.302), and −0.012 late (22 islands; CI −0.317 to +0.281). On the **same 30 pairs / 9 islands** the late-minus-early N gap was −0.162 (CI −0.591 to +0.172); no robust direction identified. Neither local exposure levels nor butterfly larval growth, true plant abundance, historical native/alien status on each island, or actual settlement was measured.

**Keep both findings in the ecological synthesis:** globally, nutrient-affinity correlates with port-record frequency even after controlling two worldwide reporting proxies; within an island, present data do not identify consistent port-specific host N enrichment. These can coexist, and different scales/testable mechanisms must not be conflated. Original GEB +54.9%, temporal prediction improvements and positive island community composition associations remain untouched, as does the frozen manuscript.
