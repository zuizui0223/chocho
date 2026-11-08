# Integrated JBI draft v0.2 — islands, introduced resources, ports and marine access

**8 November 2026. Draft for coauthor review, NOT submitted. No existing GEB manuscript, raw source, circuit matrix, current map or negative sensitivity result was rewritten.**

## Deliverables

1. [Complete integrated English manuscript](ISLAND_RESOURCES_MARINE_ACCESS_MANUSCRIPT_V02.md): structured abstract, Introduction, Materials and methods, Results, Discussion, Conclusions, figure programme and explicit inference limits.
2. [Machine-readable evidence ledger](EVIDENCE_LEDGER_V02.json): 15 source-linked results at their genuine independent population denominators and negative/null controls, plus claims withheld from publication.
3. [Five-figure architecture and draft captions](FIGURES_CAPTIONS_V02.md): one coherent resource/access story; panels are **not** all rendered.
4. [Source and number audit](../../scripts/verify_island_resources_marine_jbi_v02.py) with [dedicated GitHub Actions workflow](../../.github/workflows/chocho-island-resources-marine-jbi-v02.yml). This is **not** a statistical re-fit; it compares described numerical claims with saved source files.

## Distinct samples

| Study layer | Population and actual unit | Source status |
|---|---|---|
| Global resource contrast | 239 butterfly species, 126 strict insular + 237 mainland WGSRPD3 regions, six mixed excluded | Source-locked existing results |
| Port/whole-introduced-flora contrast | 98 GBIF botanical-effort-complete insular WGSRPD3 regions | Three existing source-adjusted coefficient estimates |
| Host/butterfly island assemblage | Original *different* 225-butterfly candidate network across 51 landmass groups / 1,275 dependent pairs | Existing source-locked index-based Jaccard and heldout models |
| Ocean current field | 182 verified marine-access ports among 194 original biological ports, 189,550 nodes / 370,996 edges | Existing verified geometry, exact numeric current run and image in PR #150 |
| Historical liner index | 17 source-covered landmass groups / 136 pairs (22 matching original ports) | Exploratory index comparison, not voyages |

## Interpretation contract

**Keep:** potential island resources +118.9% versus mainland +40.9%, matched regional/area excess +10.7 percentage points, port ↔ all introduced flora β+0.2445 and port ↔ distinct candidate hosts β+0.1639, absence of positive host-specific residual β−0.0415, host ↔ later butterfly ρ+0.300, exploratory conditional host/circuit predictive improvements.

**Also keep:** stronger earlier host–butterfly correspondence, temporal leakage of pooled host candidates, adverse latitude-adjusted wind increment, local AIS-placebo non-discrimination, full 17-island historical index failure, and the successful current-density map's inability to identify biological relay islands.

**Do not add:** unsupported “1.4% realised” number without a defined and source-verified numerator and denominator; full-area-adjustment “14 pp” in place of the actual matched 10.7 pp; direct shipping/colonisation/larval feeding causation; nutrient-driven host enrichment (owned by gesui).

## What is finished versus missing

- **Finished:** integrated manuscript text, captions/specs, evidence source crosswalk, independent existing current-map run, and manuscript-only number-audit logic.
- **Not finished:** composing the actual five multipanel final figures; a comprehensive literature/source-license audit; an independent fresh statistical reproduction of *all* previously saved analyses; journal-specific formatting and confirmed scientific overlap disclosure; coauthor approval and journal submission.
- **Coauthor decision:** the global GEB paper and this JBI draft share the 239-butterfly global resource reconstruction, so manuscript overlap must be disclosed and judged before initial submission.

## Run local manuscript-only audit

```bash
python scripts/verify_island_resources_marine_jbi_v02.py
```

The audit checks that referenced original committed source files exist and contain numerical values matching the report. The **current-map maximum numerical residual** is an exception: it lives in the successful source-locked workflow artifact, not the pre-analysis protocol committed to Git. That value must be checked using the [Figure 4 original workflow](https://github.com/zuizui0223/chocho/actions/runs/42439384640), not falsely represented as rechecked by this script.

Prior draft [v0.1](../island_port_flora_assembly_20261105/README.md) remains intact. The GEB manuscript remains intact. Do not merge into the default branch or submit without explicit approval.
