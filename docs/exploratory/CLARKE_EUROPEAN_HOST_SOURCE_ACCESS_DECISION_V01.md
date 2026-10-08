# Clarke (2024) European butterfly host evidence source audit

**Date:** 2026-10-08. **Status: ORIGINAL DRYAD DATA UNAVAILABLE FROM CI / CONSERVATIVE DERIVED-PAGE PILOT ONLY.** Main GEB manuscript, maps and frozen 239 species unaffected.

## Verified original-source blockers, not biological findings

Clarke H. E. (2024), *Ecology and Evolution* 14:e10834, DOI [10.1002/ece3.10834](https://doi.org/10.1002/ece3.10834), reports a curated evidence-ranked host-plant inventory for 464 European butterflies, synthesizing 1,119 references and recording 5,589 supported butterfly–plant associations. The dataset is deposited at [Dryad 10.5061/dryad.1vhhmgr12](https://doi.org/10.5061/dryad.1vhhmgr12); original database tables have relational identifiers, and source evidence 1–6 classes. The paper distinguishes 1 wild use plus completed breeding, 2 wild larvae and 3 wild eggs/oviposition from 4 undocumented claims, 5 captive-only rearing and 6 mistakes.

Two attempts in GitHub Actions failed:
- [Run 37774669703](https://github.com/zuizui0223/chocho/actions/runs/37774669703): original Dryad file endpoint HTTP **401**, file_stream HTTP **403**, PMC Appendix S1 response not a valid XLSX.
- [Run 37775172709](https://github.com/zuizui0223/chocho/actions/runs/37775172709): same official endpoint failures; Europe PMC official `supplementaryFiles` ZIP endpoint timed out twice (35s each). Verified original source = **none**.

These are **data-access limitations, not evidence of absence of host records, failed ecological predictions, or zero missing links**. Do not claim an original Clarke-2024 per-link audit has been completed. No source-derived correction factors to global opportunity, competition or fitness can be calculated.

## Separately sourced and less granular alternative

Butterfly Conservation Europe (BCE) [European butterfly taxonomy](https://www.bc-europe.eu/taxonomy.php) and its species [foodplant listings](https://www.bc-europe.eu/butterfly.php?genus=Iphiclides&species=podalirius) explicitly state that the larval foodplants come from Clarke (2024) **only ranks 1–3**. BCE is a **DERIVED public checklist**, not new field observations, not the full relational Dryad schema, and not byte-identical to the archived tables. In particular it does not expose the per-link evidence tier, individual bibliography, country or original uncertain/mistaken-link records.

A constrained source-first pilot is defined in `BCE_CLARKE_HIGH_EVIDENCE_HOST_PILOT_PROTOCOL_V01.json`: exact species binomial overlap of the published BCE European list and frozen chocho 239 butterflies; alphabetically first 12 eligible species selected without looking at their hosts. The code `scripts/audit_bce_clarke_evidence_host_pilot.py` records source page SHA256s, original HOSTS pinned pair spelling, exact accepted simple binomials, excluded non-specific taxa, and lexical mismatches. A BCE-only name is **not** evidence that HOSTS lacks the biological association until synonymy is independently checked. Even if the pilot succeeds, it does not represent all 239 focal butterflies or verify colonization and population fitness.

## Biological decision

The relevant global inference test is whether botanical redistribution *when conditioned on genuinely used host relations* still produces resource geography expansion/homogenization in a **representatively matched sample**, not whether a published checklist contains more plant names. We cannot make that claim yet. Stop adding area/network null tests until original or verified derived link data and exact taxon reconciliation meet their source/coverage gates.

**No changes to GEB manuscript or PR #38.**
