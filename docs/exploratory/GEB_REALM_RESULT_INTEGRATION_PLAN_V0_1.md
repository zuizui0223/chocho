# Conditional GEB realm-result integration plan v0.1

Status: **frozen before opening the recovered realized-realm result**.

This plan controls manuscript restructuring only. It does not alter any ecological threshold, species panel, occurrence query, realm definition, statistical estimand, or claim-promotion rule.

## Current manuscript hierarchy

The current manuscript has a strong 239-species resource-side core and a much thinner 32-species independent occurrence panel / 24-species climate-informative extension.

The realm audit is designed to answer a different and more central GEB question: whether the already-observed portfolio architecture generalizes across the butterflies' own realized broad-scale zoogeographic contexts.

## Conditional integration

### If the realm audit is NOT_EVALUABLE

- Do not alter the current main-figure sequence.
- Do not add a realized-realm statement to the Abstract, title, or main conclusion.
- Preserve the NOT_EVALUABLE audit and transport provenance outside the main narrative.
- The climate extension retains its current role.

### If the claim-promotion level is `descriptive_main_text`

- Add one short Results paragraph after the family/generalization sensitivity.
- Put the realm map and realm-stratified coefficients in the Supplement.
- Describe the result only as a post-hoc realized-geography audit.
- Do not add realm wording to the Abstract.
- Keep the climate panel in the main manuscript.

### If the claim-promotion level is `robust_main_text`

Promote realized-realm replication into the main paper.

The preferred main-text sequence becomes:

1. anthropogenic host redistribution changes reconstructed butterfly resource geography;
2. proportional gain is weakly related to host-family breadth;
3. specialists and generalists assemble similar aggregate gain through different host portfolios;
4. the portfolio architecture repeats across butterfly families;
5. the architecture persists within realized zoogeographic realms and under spatial-thinning / leave-one-realm-out checks;
6. climate filtering is retained as a narrower independent boundary test.

In this state, the realm analysis is more central to GEB than the 24-species climate result.

### If the claim-promotion level is `abstract_eligible`

Use the robust-main-text structure above and allow one concise Abstract clause about realized zoogeographic robustness.

Do **not** put realm names in the title. The title should remain about butterfly specialization, resource geography, and host portfolios.

## Figure policy

If the realm audit reaches `robust_main_text`, replace the current climate-centered final main figure with a realm-generality figure and move the climate figure to the Supplement.

### Realm-generality main figure

**Panel A — Global biogeographic coverage**

- Holt et al. (2013) 11-realm world map;
- no raw GBIF point cloud;
- annotate the number of realm-informative butterflies assigned to each interpretable primary realm;
- distinguish non-interpretable small realm cells visually but do not pool or merge them post hoc.

**Panel B — Realm-specific portfolio architecture**

For every realm with n >= 10:

- Spearman(host-family breadth, effective contributor number);
- sign-oriented Spearman(host-family breadth, 1 - maximum single-host share), or equivalently display the negative of the maximum-share correlation;
- show n for every realm.

The panel must visually distinguish realm-specific descriptive estimates from the pooled permutation result.

**Panel C — Magnitude versus architecture**

On the identical interpretable-realm species population show:

- pooled within-realm association with total log resource expansion;
- pooled within-realm association with effective contributor number;
- dominance-oriented pooled within-realm association;
- stratified-bootstrap 95% intervals for the architecture-minus-magnitude contrasts.

This panel is the direct visual statement of the ecological result: specialization may say more about **how opportunity is assembled** than about **how much opportunity is gained**.

## Climate-panel policy

The 24-species climate result remains scientifically useful because it is independently frozen and because it establishes that resource opportunity is not equivalent to realized distribution.

However, if realm replication is robust, the climate result should no longer carry the paper's global-generalization burden.

Move the climate figure to the Supplement and keep a concise main-text statement:

- climatic filtering is common within reconstructed resource opportunity;
- broader host-family diets did not detectably weaken that filtering in the independent 24-species panel;
- this supports the dimensional-separation argument but is not the evidence for global biogeographic replication.

## Story after robust realm replication

The ecological story should be phrased around a hierarchy:

**taxonomic diet breadth -> species-level host portfolio -> geographic resource opportunity -> realized biogeographic context**

Anthropogenic redistribution acts on host geography. Similar proportional opportunity gains can therefore emerge at different taxonomic breadths, while the underlying interaction portfolio remains strongly structured by specialization. The family and realm audits ask whether that portfolio rule is local to particular lineages or regions. Only if the frozen promotion criteria pass should the manuscript say that the rule persists across both.

## Forbidden post-result restructuring

- Do not merge Holt realms to rescue small cells.
- Do not lower the n >= 10 realm threshold.
- Do not lower the 150-species evaluability threshold.
- Do not choose raw-record versus 1-degree-cell weighting based on which gives the stronger result.
- Do not replace rejected GBIF taxa.
- Do not move the realm result into the Abstract unless the frozen `abstract_eligible` rule passes.
- Do not claim equivalence of the magnitude association to zero from this audit.

Machine-readable promotion rule: `provenance/geb/geb_realm_claim_promotion_rule_v0.1.json`.
