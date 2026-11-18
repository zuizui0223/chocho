# chocho source-qualified results consolidated on main — 2026-11-18

## Purpose
Preserve original source-derived positive findings **and their counterchecks/limits** from the integrated island, port, flora, ocean and butterfly research branches without accepting all stacked draft PRs as submission-ready. These are read-only evidence assets; they do not alter the frozen canonical GEB v0.2 or require that the GEB v0.3 draft has received author/coauthor approval.

## Merge method and exact scope
- Base main commit before import: `13341c97d9c2c54889e2ad605282794a9d2c2849` (base tree `527316a60e2d66449f230bc56fd72e4446395a18`).
- Source branches (PR # and pinned commit): [#177](https://github.com/zuizui0223/chocho/pull/177) `417bbc77f2736d0bff84b9b73547b83525d3f2c4`; [#125](https://github.com/zuizui0223/chocho/pull/125) `ff642c878e67850dc6429cd9deea61cf1ff0cce7`; [#117](https://github.com/zuizui0223/chocho/pull/117) `bb9f7d08e09416635725359eb4c7c3c48b866a33`; [#122](https://github.com/zuizui0223/chocho/pull/122) `b4f331c4eacc2105187644ed7b84e178a333248f`; [#120](https://github.com/zuizui0223/chocho/pull/120) `dfd49cfcf9fb01a3db47cc611aae8cdd9c53cda0`.
- **811 new, distinct source blobs** selected from branch trees, with **zero overwrites of preexisting main files** and **zero differing-source-path collisions**. Per source branch: #177: 503, #125: 174, #117: 82, #122: 18, #120: 34.
- Added source files live at original paths under `results/`, `provenance/reviewer_defenses/results/`, selected `manuscript/integrated_butterfly_resource_access_20261109/` and source-labeled `figures/integrated_island_port_butterfly_v03/`. Only textual reports, machine-readable numerical summaries/tables (max 100KB per blob), source SVG/JSON captions were copied. Large raw inputs, images/binary intermediate arrays, scripts, CI workflows, and old competing manuscripts **were not wholesale merged**.
- Git file SHA for every new path is preserved **byte-exact** (without changing scientific values), recorded in [POSITIVE_RESULT_FILE_MANIFEST_20261118.tsv](POSITIVE_RESULT_FILE_MANIFEST_20261118.tsv).
- This is a **selective evidence/result merge**, NOT the GitHub UI action “Merge pull request” on the open stacked development branches; those contain overlapping and untested combinations. Source scripts and complete original PR revisions remain available at linked immutable heads. Main CI result on this new commit must be checked independently.

## Claim map with limitations
Read [full positive/countercheck register](POSITIVE_EVIDENCE_AND_ONE_STORY_20261118.md), [11–18 November priority results](LAST_WEEK_POSITIVE_EVIDENCE_20261118.md), and [all 176 PR audit](ALL_PR_AUDIT_20261118.md). All are needed together; no cherry-picked p-values without later falsifications.

### Strong primary observations
1. 239 originally HOSTS/WCVP-eligible butterflies: documented fixed-host plant geography +54.9% source-assembled **potential** butterfly×region units; excludes absolute actual butterfly population response.
2. Insular botanical region relative gains +118.9% vs mainland +40.9% are **initial native-resource baseline amplified**, not disproportionately higher island gap filling.
3. Port counts positively correlate with island-wide introduced *plant floras*; port-specific host enrichment independent of all introduced plants and causal historical shipping are not established.
4. Earlier known-host geography adds information to 155 known-species later GBIF index predictions (+4.12% MSE reduction); static contemporary plant predictor is not timed natural introduction.
5. Oceanic circuit/water geometry predicts plant and index-butterfly community structure on **different geographic populations**; claimed host×maritime synergy and mechanism-specific ship/flow paths do not survive identification gates.
6. Independently sourced `Pieris cheiranthi` on `Tropaeolum majus` provides direct **larval host use**, not colonial transport, new island establishment, survival, or population gain.
7. Independent checklist 41 living Macaronesian butterflies includes 18 extant endemics, **0/18** in original 239. 8 species with source-identified host species are not representative; global potential 4/8 positive source-list outcomes are not endemic mean prevalence.

### Explicitly preserved limiting tests
- Location-attributed HOSTS reports are island-alien biased, but family/whole-flora adjusted extra biological preference not established.
- Port×island-area/native-flora associations are exploratory and do not transfer reliably across botanical regions.
- NOAA currents and AIS historical links remain sparse and cannot demonstrate plant seed landing or port-caused butterfly range expansion.
- Globally independent GBIF-indexed assemblage beta model shows **small** additional known-host predictive contribution (+0.157% MSE) and record-effort dependence.
- Historical pre-2001 review finds **145/317** supposedly newly GBIF-indexed later cells had earlier records: first index ≠ first biological arrival.
- Distinct monarch, LGM, Kyoto greenhouse and UKBMS satellite projects have separate question/target; archive PR links in the all-PR ledger but do not install off-topic code in this one GEB paper.

## Publication/editorial state
This main evidence package is source-linked, **not an accepted peer-reviewed causal synthesis**. The separate integrated v0.3 paper remains a draft candidate and must be checked for data rights, DOI source attribution, figure readability, model-to-claim congruence and coauthor approval before submission. This archive was not an occasion to change any coefficient, reconstruct original raw restricted observations, or imply a full 4-stage ship→plant→larva→breeding test.

## Next housekeeping
Close *separately scoped/superseded* off-topic PRs **only after keeping their branch names and source citations**; leave core source-relevant stacked PRs open until a tested integrated branch handles full code/CI/manuscript semantics. Closing a PR does not delete its branch or observations.
