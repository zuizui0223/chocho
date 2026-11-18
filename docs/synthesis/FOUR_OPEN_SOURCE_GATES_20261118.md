# chocho — final four active PRs and remaining gates (2026-11-18)

## Actual GitHub state at this audit
**176 original PRs: 46 merged, 126 closed-unmerged, 4 open.** All closures were reversible PR-state changes; no branch deletions, forced ref updates, mass raw-code merges or journal submissions.

### Primary goal
**One** Global Ecology and Biogeography butterfly paper: how redistribution of already-known larval host plants changes *potential* butterfly food-resource geography, why low native-host coverage amplifies island gains, and how actual host use and marine-record geography filter ecological expression. Positive findings are in [22-source claim register](POSITIVE_EVIDENCE_AND_ONE_STORY_20261118.md), with later counterchecks retained.

## 4 remaining PRs (no duplicate submitted papers)

| PR | Purpose and exact next gate | Status |
|---|---|---|
| [#157](https://github.com/zuizui0223/chocho/pull/157) | The **only** integrated GEB v0.3 candidate manuscript/figures. Confirm all Results, claims, source inputs, protected data licensing and attribution, fifth-figure captions and manuscript bundle. Because original PR is stacked on an exploratory base, build a **clean candidate from main**, do not merge the full stacked history. Complete actual coauthor review before submission. | Keep open, Draft; not approved/submitted |
| [#177](https://github.com/zuizui0223/chocho/pull/177) | Independent 41-living/18-endemic Macaronesian butterfly external-validity, 8/18 exactly documented host-species source panel, directly photographed endemic *Pieris cheiranthi* larvae using alien *Tropaeolum majus* on Tenerife and La Palma. Complete source-type caption and **add as SI plus qualified field-use example to the SAME GEB article**, not a population prevalence or new-island arrival estimate. PR base is closed #176 but its branch still exists. | Keep open, Draft; source audit to reconcile with v0.3 |
| [#127](https://github.com/zuizui0223/chocho/pull/127) | Original 51-island observed uniform-ocean circuit geometry input, GSHHG ocean mask, grid/source controls, AIS/native-source alternative, leave-whole-island-out output. Restore original large masks/required inputs from immutable source and perform full raw-data replay or clearly document external-reproduction-only boundary. | Keep open for original-source restoration |
| [#128](https://github.com/zuizui0223/chocho/pull/128) | Exact same 51-island butterfly/host network: recompute original host/circuit model comparison, 1,275 dyads and both endpoints unseen sensitivity. Source-geography, survey-effort, native-host and AIS/wind controls must remain visible. Fold/interaction MSE scores must not be treated as causal shipping/settlement. Could combine with #127 into **one** reproducible marine source workstream. | Keep open pending model-source replay |

## Work already finished in main
- [Initial 811-result file exact-blob archive](SELECTIVE_RESULT_MERGE_20261118.md) and [71 missing source-result/verification files](D46_SOURCE_RECOVERY_20261118.md), both no overwrite of original frozen GEB v0.2.
- [32 source-code analysis script archive](CORE_SOURCE_CODE_BRIDGE_20261118.md) plus untouched main underlying original frozen resources. Source of every new blob and audit retained.
- [18-input arithmetic/source-identity/claim boundary CI](INTEGRATED_EVIDENCE_VERIFICATION.md), Actions [44247388248](https://github.com/zuizui0223/chocho/actions/runs/44247388248) SUCCESS; flags 126-region port native-host-resource covariate and separate 123-native-flora model.
- [Original scientific-source script verification](https://github.com/zuizui0223/chocho/actions/runs/44249244575) SUCCESS in 5 steps: 155-known-species prediction losses, 239 species island/mainland source contrast, HOSTS reported larval host source, independently enumerated 42 butterfly origin/selection screen, marine Figure 4 source.
- [44 stacked source-PR closures with head-SHA audit](D46_CLOSURE_DECISION_20261118.md), [12 core source-PR closures](CORE12_SOURCE_RECONCILIATION_20261118.md), and older 67 closures (56 outside + 11 redundant).
- Main original v0.2 paper and blinded release package unchanged; v0.3 candidate remains a distinct archived Draft and needs coauthor approval.

## Remaining scientific limitation (not silently assumed solved)
Current verified outcomes have different observational support: 239 source-selected butterflies' botanical region potential, 155 temporally scored GBIF index species, 224 global indexed butterfly beta panel, 126/123/98 port flora areas, 104/51 sea-geographic panels and nonrepresentative 8/18 external endemic exact-host cases. **No unified observational case series proves a chain ship -> introduced host -> wild larval recruitment -> new island butterfly colonization**. Original source processing for selected models succeeded in archived Actions, but the latest main CI is a frozen receipt/source-code replay gate, not a raw refit of every historical experiment.

## End condition
After source masks and experimental CI are reconciled for #127/#128, close them (without deleting branches). Incorporate #177 as SI/field-use example in a *clean* #157 GEB candidate branch. Then a **single active manuscript PR** remains for article approval, or none after formally approved merge. Do not submit without explicit author/coauthor decision.
