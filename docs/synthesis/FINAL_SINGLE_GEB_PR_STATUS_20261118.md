# chocho — one active GEB manuscript PR (2026-11-18, post-cleanup)

## Verified GitHub repository state
- **178 existing PRs:** 47 merged, 130 closed-unmerged (reversible, no source branch deletion), **1 currently open**: [#179](https://github.com/zuizui0223/chocho/pull/179).
- #179 is **Draft** on `paper/geb-single-island-port-endemics-20261118`, based directly on `main` at `75f2c67c`, with original stacked manuscript #157 and independent endemic source #177 now closed as superseded. Their original branches and citations remain accessible.
- #127 and #128 are closed because clean-main source model replay succeeded, **not** because marine transport or demographic colonization was established.
- Original GEB v0.2 submission manuscript and source inputs remain frozen on `main`. The integrated GEB v0.3 candidate text is changed **in the draft PR only**, not yet merged into main.

## Actual source and CI gates completed
- [Main source 51-island ocean replay](https://github.com/zuizui0223/chocho/actions/runs/44251956538): original 12 sea-circuit source-model comparisons reproduced, including 5.082% held-island late butterfly Jaccard prediction gain. Retrospective geography only.
- [Main source 51-island joint host/marine replay](https://github.com/zuizui0223/chocho/actions/runs/44251956521): 45 original tests, 25,500 predictions, 510 model folds and locked original comparisons reproduced.
- [#179 manuscript source-denominator and field evidence test](https://github.com/zuizui0223/chocho/actions/runs/44255809108): PASS, including exact 0/18 original-endemic selection, eight source-documented host species, two named-island caterpillar records and matched whole-/family-flora limits.
- At PR #179 head `cfd2123f886fdcb5fa4f9c6caee2c452d23b6c73` all five requested checks succeeded: manuscript source audit, paper CI, blinded DOCX, v0.2 figure render, anonymous reviewer bundle. The earlier main anonymous-bundle failure from a v0.3 source verifier incorrectly included in v0.2 snapshot was corrected in main commit `75f2c67c`.

## One-paper ecological claim
Human redistribution of plants already documented as larval food reshapes **potential** butterfly resource geography, disproportionately on the sampled islands because initial native-host coverage is low (239-source global +54.9%; island +118.9%; mainland +40.9%). Port counts relate to introduced whole floras; independent wild/garden larval observations document actual use of one introduced plant by a Canary endemic; known-host geography adds 4.12% error reduction to later GBIF recording prediction for **known butterflies**, and 51-island ocean geography adds separate retrospective prediction information. These layers are complementary, not a single identified historical plant delivery → new butterfly population causal chain.

## Independent checklist / source-selection distinction
The externally enumerated Macaronesian population has 41 living butterflies and 18 independently classified living endemics. **0/18** appear in the original final 239, so those sample effects do not directly estimate an island-endemic mean. The eight named-host species are source-selected; 4/8 have additional *global botanical resource regions*, none has newly host-present original home-archipelago regions under partial catalogues. Observed Tenerife wild L4 and La Palma garden larvae of *Pieris cheiranthi* on introduced *Tropaeolum majus* demonstrate **larval plant use**, not adult recruitment, newly occupied islands or average endemic host dependence. A family-matched flora baseline leaves no identified surplus alien-host preference beyond plant supply.

## Remaining author/coauthor release gates — not executed in this task
1. Final author and coauthor approval for scientific narrative, figure selection and the submitted journal version.
2. Check exact source photo photographer attribution/reuse rights; observations are referenced by source URL only, no third-party photographs embedded.
3. Dataset DOI/version, original data permissions, license requirements and reference formatting; ensure any journal supplementary archive excludes author-only raw HOSTS edges and human identity.
4. Final GEB manuscript/figure production as approved; a successful v0.2 DOCX job is **not** approval of the new v0.3 draft.
5. Formal submission only upon corresponding author direction, and no claims of historically verified ship-mediated plant arrival, new butterfly colonization, widespread realized homogenization or population increase.

## Starting point
- [PR #179 clean GEB candidate](https://github.com/zuizui0223/chocho/pull/179)
- [Original 22-claim registered results and bounds](POSITIVE_EVIDENCE_AND_ONE_STORY_20261118.md)
- [Main integrated ocean and original network evidence](MARINE_RAW_REPLAY_SOURCE_20261118.md)
- [One-paper evidence bridge in the clean candidate branch](https://github.com/zuizui0223/chocho/blob/paper/geb-single-island-port-endemics-20261118/manuscript/integrated_butterfly_resource_access_20261109/ONE_PAPER_EVIDENCE_BRIDGE_20261118.md) (PR draft, not yet on main).
