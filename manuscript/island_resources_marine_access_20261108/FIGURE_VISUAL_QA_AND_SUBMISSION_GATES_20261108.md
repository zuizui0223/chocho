# JBI v0.2 visual figure QA and source-inference gate — 2026-11-08

**Status:** scientific-draft visual audit of the actual rendered figures. No new model fit, new data selection, journal submission, or changes to the global GEB paper. All figure results are exploratory or source-reconstructed at their declared units.

## Outcome of the inspection and corrections

| Figure | Source units | Checked correction and decision | Source-verified runner |
|---|---|---|---|
| 1: framework | 239 butterflies × WGSRPD3 botanical regions, with separate 225-candidate/51-island source | **PASS after repair** — flora assembly points toward reconstructed host opportunity; dashed arrows to butterfly records state association/prediction, not observed migration | [Corrected Fig1 run](https://github.com/zuizui0223/chocho/actions/runs/42444942484) |
| 2: island resource contrast | 126 insular and 237 mainland regions; 15 common Level1 × area strata | **PASS after editorial repair** — 4,829→10,570 (+118.9%) versus 21,215→29,901 (+40.9%); panel C now displays all 15 matched contrasts and the pre-existing +10.7 pp weighted result. The distinct 13/920 GBIF indexed-overlap screen is Supplementary context, not a biological realization rate | [Matched-strata Fig2 run](https://github.com/zuizui0223/chocho/actions/runs/42445180860) |
| 3: port floras | 98 botanical-effort-complete insular regions | **PASS** — +0.2445 whole introduced flora, +0.1639 introduced distinct host plants, −0.0415 conditional host β; current rendered annotation no longer collides with long category labels, and no invented confidence intervals are shown | [Current Fig3 run](https://github.com/zuizui0223/chocho/actions/runs/42444942444) |
| 4: current, host concordance and heldout predictions | 51 original landmass groups, 1,275 dependent pairs, 182 qualifying marine terminals | **PASS after repair** — 4A is uniform-water electrical current, not shipping or real relay; 4B uses legible, log-count hexbin source-pair density with colourbar and preserves original rho values; 4C separates one-island and two-island-unseen reference baselines | [Revised full Fig4 composite](https://github.com/zuizui0223/chocho/actions/runs/42444942451) |
| 5: maritime specificity | 51-island modern proxy test and distinct 17-island historical source | **PASS** — +1.002% latitude-adjusted AIS coexists with 100/199 locally constrained surrogates at least as predictive; historical 17-island −2.923% is the primary result, while +8.146% is a training-dependent 16-island sensitivity | [Current Fig5 run](https://github.com/zuizui0223/chocho/actions/runs/42444942459) |

## Separate technical gates

- Existing verified marine source [PR #150](https://github.com/zuizui0223/chocho/pull/150) provides 189,550 connected marine nodes and 370,996 undirected edges, and all 1,275 uniform-pair current solutions reconstructed by a 50-potential superposition procedure. The measured maximum grounded linear residual is approximately 1.33e-13 in the original executable receipt. This is **not** a test of the biological reality of the conductance map.
- The [JBI v0.2 original-source auditor](../../scripts/verify_island_resources_marine_jbi_v02.py) checks ledger and recorded claim numbers against existing committed results. It is **not** an independent statistical refit.
- The [five-figure workflow index](README.md) and [source-safe figure captions](FIGURES_CAPTIONS_V02.md) maintain the distinct units and temporal limitations. The current rendered figures were visually inspected as PNGs after the recorded repairs, not merely inferred from successful Python exit codes.
- The evidence ledger retains stronger host correspondence to earlier butterflies (ρ+0.418) than later butterflies (ρ+0.300), unlike a prospective arrival claim. Candidate hosts are indexed over 2001–2024 and can postdate the 2011–2022 later-butterfly response. There are 1,275 dependent source pairs, not independent field replicates.
- Spatially local AIS graph-label surrogates are as predictive as the actual AIS assignment in 100/199 cases. The historical 17-island port-connectivity index degrades the full sample prediction. Neither finding erases source-locked positive botanical and conditional-prediction observations; both prevent promotion to **identified ship effects**.
- The global GEB source shares the 239-butterfly HOSTS × WCVP reconstruction. The global +54.9% result is shared, not an independent second result.

## Remaining pre-submission scientific and editorial gates

1. Verify bibliography and terminology against the original literature, not generated citations.
2. Verify all source licenses, raw data and figure availability; ensure the complete figure provenance (GitHub Actions assets expire) has durable journal-suitable storage.
3. Review and explicitly disclose scientific overlap with the GEB paper to coauthors and editors as appropriate.
4. Decide whether the exploratory, retrospectively selected 51-island predictive comparisons are strong enough as main figures and how to communicate the unavoidable spatial/time limitations.
5. Obtain coauthor approval and journal format/data statements before submission.

**Editorial judgement:** These figures support a coherent observational resource-and-access geography article. They do not yet support a claim that islands are verified butterfly-transit relays, that ports preferentially carry larval hosts, or that newly indexed butterflies colonised those exact added resource cells.
