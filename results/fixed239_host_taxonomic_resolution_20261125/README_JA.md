# Fixed239 × 126 hierarchy: successful source-aligned reproduction (25 November 2026)

**Run:** [45560130341](https://github.com/zuizui0223/chocho/actions/runs/45560130341), exact analysis head 77c76099bd492d6a9f5cf989196c84946e7d16aa, all steps passed. [Artifact 14997553306](https://github.com/zuizui0223/chocho/actions/runs/45560130341/artifacts/14997553306) downloaded and its SHA-256 verified locally against the Actions digest 4e665dd376606782a84fbb5e2517cc32cd54abcedf0a0a44f3591d09b001ceb5. The five ZIP members include the full report, source gate, per-island counts, mismatch table and 37 extra resolved host-ID taxonomy records. No raw third-party data published by this receipt.

## Exact source checks
All 239 individual butterfly global native/current opportunity counts, all 126 regional native/current counts, and original SHA/row/edge gates passed. Zero island mismatches; island-native 4,829, current 10,570 butterfly-region opportunities exactly. BOU retained with zero counts. 37 resolved host IDs outside the strict accepted-species background retained in the original Q1 host distribution; broader ranks use the background union locally present resolved hosts. Four scratch scripts inspected, but their input pickles and original output CSV are still unavailable: this reproduces the pattern with the source-aligned analysis, not the exact scratch files.

## Supported and bounded findings
| Outcome | Base beta | HC3 p | Level1 t8 p | MEM10 beta | MEM10 HC3 p | MEM10 + t8 p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Native family coverage count | -0.0126 | .7793 | .7869 | +0.0051 | .9236 | .9245 |
| Native genus coverage count | -0.0928 | .1087 | .2132 | -0.0628 | .3169 | .4038 |
| Native species coverage count | -0.1676 | .0000954 | .0140 | -0.1330 | .0204 | .1095 |
| Direct paired native species minus family z-slopes | -0.1550 | .000237 | .00809 | -0.1382 | .00879 | .0587 |
| Direct paired native species minus genus z-slopes | -0.0748 | .0666 | .0881 | -0.0702 | .1038 | .1241 |
| Native log((species+.5)/(family+.5)) | -0.2565 | .0000227 | .00892 | -0.2215 | .00721 | .0793 |
| Native plants belonging to documented host families (log odds, n=125) | -0.2639 | .000195 | .00967 | -0.1726 | .0157 | .0281 |
| Single native-region share (standardized logit, n=125) | +0.3278 | .00163 | .00988 | +0.2662 | .0187 | .0493 |

Nearest/farthest distance quartiles each n=32: species coverage medians 74.5/6.5, genus 158/54.5, family 227/189.5; not matched adjusted effects. Removing all butterflies with any recorded Fabaceae/Poaceae host leaves 120; native family beta -0.0130 (t8 p=.849), species -0.2178 (p=.00480), MEM10 + t8 species p=.0586. Mean documented-family-set coverage shows weaker negative association (standardized beta -.1104, HC3 p=.0740, t8 p=.1255). Presence saturation is not fully excluded. The proportion of distinct native families containing a focal documented host has uncertain cluster/MEM10 effects. In contrast the proportion of native **species** belonging to such families declines robustly even after MEM10 plus cluster inference. Distinguish family presence from within-family representation; do not declare classical disharmony absent.

98 common islands with botanical indexing-volume covariate: species beta -.2162 (t8 p=.0147), log species/family -.2976 (p=.0140), host-family species share -.2463 (p=.00640). Family coverage now also negative (HC3 p=.0321, t8 p=.1206); genus HC3 p=.00312, t8 p=.0738. No proof of exactly zero family effect. Plant indexing is not host-interaction documentation effort. For the same 125 nonempty regions the native log species/genus distance coefficient changes from -.1790 to -.1130 after standardized logit single-native-region share (t8 p=.156); the species/family ratio moves -.2439 to -.2157 (p=.067). This is partial attenuation, not causal mediation, and does not reproduce a literal halving or verified endemic identity.

Contemporary all three rank-count slopes are near zero/positive and not resolved (base species +.0192, genus +.0573, family +.0340). Paired current-minus-native slopes are positive but broad-cluster uncertain: species +.1868, t8 p=.0677; genus +.1501, p=.0919; family +.0466, p=.2452. No equality or complete erasure claim; no stronger evidence for a uniquely species-specific introduced-plant effect without a direct three-way contrast.

## Candidate manuscript wording (not inserted)
"Across the same 126 botanical island regions, isolation was more negatively associated with the coverage of documented larval host identities than with the presence of their host families (difference in standardized slopes, beta = -0.155; nine-Level1 clustered p = 0.008; MEM10 plus clustered p = 0.059). More isolated islands also had lower representation of native plant species in the documented host families. These patterns suggest taxonomic identity and floristic composition jointly structure the mapped resource gradient. They cannot distinguish unrecorded host associations from genuine incompatibility."

For Discussion: broad host-family representation can remain despite reduced coverage by familiar host identities; remote floras contain more narrowly distributed native plants, and global catalogue gaps may systematically undercount island-specific interactions. Neither available family members nor observed plant geographic occurrence imply successful larval use. Preserve the localized HOSTS result separately and do not ascribe the attenuation to ports or infer observed colonization. These are exploratory endpoints with shared inputs and limited broad geographic replication.

The frozen original GEB draft PR179 has not changed; this PR192 retains the source audit and candidate wording only. Raw scratch input recovery and verified independent host documentation are separate tasks. No merge, submission, reviewer access, or coauthor approval.
