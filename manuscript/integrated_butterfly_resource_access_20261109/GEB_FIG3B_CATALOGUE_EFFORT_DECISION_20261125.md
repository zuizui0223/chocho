# GEB Figure 3B — catalogue record-depth audit and fixed decision (25 November 2026)

User decision: stop exploratory directional-current analysis, retain all outcomes in Appendix S3.5, and concentrate on Figure 3B. No further HYCOM/OSCAR threshold, geographic subset, release-date or model search. Paper #191 changes the primary text to an explicit lack of direction-specific support; main positive island/port/morphology/circuit results are unchanged.

## Source-bounded Figure 3B question

Does the original fixed-239 species association between pre-existing native global host-resource geography and added island resource counts persist after accounting for **the actual number of catalogue rows**, beyond the prior distinct plant species/families proxy?

Existing values: raw Spearman rho +0.578896407; partial rho +0.348264 controlling resolved documented host plants and families. Those distinct taxa are both ecological variables and documentation-richness proxies. They cannot establish independence from research effort.

The exact GloBI HOSTS mirror 808e0b869f9ec1adf8efff87cf6a395adda103e0 is an NHM query DOI 10.5519/qd.bsucrxdz with 140,485 rows. Its declared schema has _id, HOSTS ID, insect/plant taxonomic fields, location, lab rearing and damage, but no explicit publication/reference/DOI identifier. Neither identifier is assumed to represent a study. Do not invent bibliographic effort. Raw catalogue rows are not independent observations.

## Locked next test

[Workflow commit 244d6c5](https://github.com/zuizui0223/chocho/commit/244d6c5ed3b0234433049fd450f3039c2394c220); [run 45544732456](https://github.com/zuizui0223/chocho/actions/runs/45544732456). At writing the run is in progress, not a result.

- Exact original HOSTS CSV blob, WCVP v13 and original mapping code, fixed239 source metrics and 126-island additions. Require all 239 exact butterfly labels and reconstructed resolved-host counts to agree; mismatches trigger source diagnosis, not selective exclusion or an automatic synonym change.
- Primary additional control: rank of all catalogue rows for each original butterfly label (adjusted for resolved species and family counts). Report all predefined sensitivities using exact non-ID row-content deduplication, resolved row counts, deduplicated resolved rows and resolved rows per host. Distinct identifier strings and locations are coverage diagnostics only, not publication counts.
- Reproduce raw and prior adjusted associations first. Descriptive 2,000 species bootstrap and all leave-one-species checks; same fixed adequate-taxonomy and one-family subsets. Do not claim independent effort, phylogenetic correction, causal identification or realized butterfly establishment. All results, including weak or adverse sensitivities, will be retained.

If exact mapping fails, stop fitting and reconcile the original source definition. If the new proxy weakens the association, revise Figure 3B interpretation rather than search for an alternative favorable specification. If it persists, retain it as a bounded same-catalogue diagnostic, not proof of effort independence. Neither paper figures nor scientific numbers are changed by this workflow. Internal aggregate-only output, no raw source redistribution or submission.

## Completed result — 25 November

[Run 45544732456](https://github.com/zuizui0223/chocho/actions/runs/45544732456) PASS; [signed internal artifact #14990294609](https://github.com/zuizui0223/chocho/actions/runs/45544732456/artifacts/14990294609), SHA-256 `b65982da1e8c1051462f36f15347fb4e8acef5a14fb2996500836bb86a108e50`. All 239 frozen species and resolved host counts match. Raw 6672 rows, non-ID deduplicated 6649; complete binomial 4695, resolved 3898, resolved deduplicated 3890. No explicit independent publication/reference field. Prior +0.348264; primary row-adjusted +0.357252 (descriptive 95% +0.238345,+0.474346), four alternatives +0.356964,+0.348422,+0.348357,+0.347793; all 2000 bootstrap samples finite and all reported intervals positive. Raw rows correlate +0.911700 with resolved host counts, so no independent effort identification. Retain all original source-defined patterns; no further proxy searching.

Paper [commit 1568c06](https://github.com/zuizui0223/chocho/commit/1568c066b59e579cdc4d3504683af5221763f928) integrates a bounded main-text sentence and full S4.1 diagnostics with exact logged JSON, all previous figure sources unchanged. Export QA pending at checkpoint. This is not authorization of source redistribution or reviewer access. Internal 217-file archive predates this new S4.1 result and must be updated separately.
