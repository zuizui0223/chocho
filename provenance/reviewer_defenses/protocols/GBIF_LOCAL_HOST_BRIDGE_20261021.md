# Strict HOSTS/WCVP-to-GADM1 plant co-detection bridge: pre-outcome design

**21 October 2026, exploratory, source-bounded.** This bridges the
fixed butterfly-host plant network to GBIF first-level administrative-unit
facets without equating plant geography, butterfly detection, life-stage
labels, feeding or recruitment. The original GEB manuscript is unchanged.

## Frozen butterfly and plant cohorts

The exact current country analysis uses 239 original butterflies, a
BCE/Clarke-supplemented HOSTS host graph, fixed accepted WCVP plant IDs,
source-frozen native/current WGSRPD3 ranges and CHELSA data; 229
noncolliding butterfly taxa contribute country records. We retain all
48 existing butterfly×country-positive cells and derive source candidate
introduced host plants strictly from that **unchanged** network,
using the same country→WGSRPD3 material intersection links.
This is a **source-host candidate inventory**, not observation.

The initial focal *plant-observation* cohort is fixed by the stage evidence
**before inspecting any candidate host names or plant occurrence responses**:
five added-resource butterfly×country cells with previously recorded
Egg/Larva/Pupa/Caterpillar stage strings: `Anteos clorinde` (CU),
`Colobura dirce` (VE), `Papilio demoleus` (CU, DO, HT).
All six remaining juvenile-labeled country cells are retained as
comparison/context but are **not falsely coded as host-absent sites**.

## Host candidate mapping

Restore only the original author-only source artifact (run
39006554609, SHA256
`408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323`)
inside a permitted GitHub Actions runner. Reconstruct the exact
`source_augmented_model` with audited HOSTS/WCVP/BCE data.
For every butterfly×country in the 48 source-positive cells:
- Require the same `material_link=1` country-WGSRPD3 intersections;
  require the original country source to lack any native host in
  `guard_link=1` regions.
- Record each accepted WCVP host-plant **species name/ID** only when
  that host has a source introduced-range cell in at least one
  material-linked WGSRPD3 unit. Record the matching botanical units.
- Independently reproduce the **26 positive added-host country flags**
  and **22 non-added flags**; fail closed on mismatch.
- Do not alter original source host list, map introduced hosts into native
  country status, turn plants in cultivation into wild hosts, or
  silently reconcile ambiguous plant scientific names.
- Keep raw author-only HOSTS/WCVP/TDWG/CHELSA files outside the repo.
  The derived **five-case candidate taxon list** is the only necessary
  publicly reviewable host subset.

## Subcountry observational link (subsequent step)

After passing original source identity and candidate-name checks,
match *exact* botanical accepted binomials against current GBIF species
taxonomy, with rejected/ambiguous/unknown hosts separately counted.
For source-recognized host taxa in the five fixed cells, request GBIF
`limit=0&facet=gadmLevel1Gid` within the same country,
`2010–2025`, with coordinates and no geospatial issues, from accepted
specimen/observation records. Preserve each response hash, GBIF taxon key,
query URL, and rejected source. Plant observations do not verify introduced
status, naturalization, local plant accessibility, blooming/leafing season
or that larvae fed on the recorded plant.

Compare the original conditional immature-butterfly GADM1 units against
GADM1 units with *at least one indexed observation of an eligible host*.
The GADM1 intersection is administrative co-detection only. A region
with no GBIF plant record is **unrecorded**, not plant-absent.
Do not infer more butterflies, reproduction, survival or a positive
host mechanism. Report raw unmatched cases and effort/access constraints
without selectively dropping negative plant matches.

A meaningful future ecological test requires GADM1-wise non-detection
background, source-verified plant habitat, dated introduction and larval
host-use/juvenile survival. This bridge is a feasibility/measurement
result, not a new causal macroecological discovery.
