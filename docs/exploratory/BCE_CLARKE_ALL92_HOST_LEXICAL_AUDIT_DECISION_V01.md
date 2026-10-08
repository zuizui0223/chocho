# All-92 exact European butterfly source-referenced host-list comparison

**Date:** 2026-10-08. **Status:** FULL SOURCE AUDIT COMPLETE; TAXONOMY & EVIDENCE-PROVENANCE NOT YET ADJUDICATED. **No change to GEB manuscript / PR #38.**

## Scope and prior status

Original Clarke (2024) raw Dryad files remained inaccessible via the authorized automated routes (HTTP 401/403), and official publisher appendix returned a challenge page. Therefore the input here is expressly a **BCE DERIVED** open checklist, which the website states is based on Clarke (2024), DOI [10.1002/ece3.10834](https://doi.org/10.1002/ece3.10834), with host evidence restricted to levels 1–3 (breeding after wild use, wild larvae, or wild eggs/oviposition). BCE does **not** return which of those three classes or original paper supports any particular species-level interaction. These cannot be treated as independent direct field observations, nor as the complete original Clarke tables.

The independently compiled BCE checklist was assessed against the fixed `globalbioticinteractions/HOSTS` commit `808e0b869f9ec1adf8efff87cf6a395adda103e0` and the existing frozen *chocho* 239 butterfly names. The full source set and exact match method were frozen after a positive 12-species pilot but **before** the outcomes of all 92 were known, in `BCE_CLARKE_92_SPECIES_HOST_LINK_AUDIT_PROTOCOL_V01.json`. Outcome-based butterfly selection was not permitted.

## Reproduced exact-name result

[Successful GitHub Actions run 37777242991](https://github.com/zuizui0223/chocho/actions/runs/37777242991) and archived artifact `chocho-bce-clarke-all92-lexical-source-audit-v01` (ID `11551280843`) contain SHA256 provenance for BCE taxonomy page and each species page, source URLs, all exact and excluded host binomial strings and error statuses. Original 239 chocho butterfly names intersected the BCE 501-species European roster at **92 distinct EXACT binomials**. In the 92:

- **83** butterfly pages have at least one species-rank BCE plant record; **9** do not provide a species-rank table usable under the fixed rule (7 no species-rank table, 2 only genus/infraspecific names). These 9 are **not** biological hostless species.
- For the 83 pages, **1,390** BCE/Clarke-evidence-filtered distinct butterfly–plant exact binomial relations and **1,068** raw original HOSTS exact-binomial relations.
- Only **338** match both butterfly and plant exact spelling. **1,052** BCE names lack an exact textual counterpart in HOSTS, and **730** HOSTS names lack an exact counterpart in BCE.
- A further **240** BCE page rows with generic, infraspecific or otherwise unresolved plant identities were **excluded from the species-rank source comparison**. They are not assumed to be false observations.
- Differences occur in both directions: e.g., `Anthocharis cardamines` BCE=61 / HOSTS=11 / exact shared=7, whereas `Celastrina argiolus` BCE=41 / HOSTS=93 / exact shared=5. Neither direction is biological sensitivity/specificity without taxonomic and evidence provenance.

**Core finding:** Two broadly used host inventories differ substantially in species-level *reported source names*. **This is NOT the number or rate of confirmed biological host-link omissions.** Reasons include taxonomic synonymy/WCVP version, butterfly genus and subspecies concepts, European versus worldwide geography, literature coverage and distinct criteria for what constitutes a larval host.

## Decision and next gate

The primary confirmation needed is unique WCVP accepted plant ID matching for all 1,052 BCE-only exact spellings, tested against the **original fixed E. species × accepted plant ID** interaction sidecar. Exact same accepted ID despite different spelling means not a missing butterfly–plant pair. An exact new accepted ID with corroborated evidence would be a candidate additional link for a **geographic sensitivity**, not a proven local demographic effect.

The next protocol `BCE_CLARKE_TAXONOMIC_HOST_CROSSWALK_PROTOCOL_V01.json` was frozen before the 92-species result was known. Its implementation is versioned as `scripts/resolve_bce_clarke_plant_synonyms_wcvp.R` and `scripts/audit_bce_clarke_taxonomic_crosswalk.py`; the subsequent GitHub workflow is `.github/workflows/butterfly-bce-clarke-wcvp-host-crosswalk.yml`.

**Do not** update the GEB figures, 239-species result, association nulls or manuscript until source matches, denominators, reference overlap and geographically stratified sensitivity are independently audited. No colonization, host fitness, consumer abundance or causal introduction response is inferred here.
