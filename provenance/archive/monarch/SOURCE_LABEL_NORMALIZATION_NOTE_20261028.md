# Source vocabulary normalization audit — source correction, not hypothesis revision

Date: 2026-10-28. The pre-outcome scientific contract is unchanged: native hosts, cultivated **yes** vs **no**, November–February against March–October, fixed host×State×year matching, all original minimum information gates.

The first numerical workflow [#40541074571](https://github.com/zuizui0223/chocho/actions/runs/40541074571) returned 0 eligible rows because code expected full genus–species binomials in `host`. It was a parsing error, **not** biological absence. Independently run source format audit [#40541666884](https://github.com/zuizui0223/chocho/actions/runs/40541666884) validated original author vocabulary *without opening seasonal outcome*:

- All 12,942 author-thinned rows have `State` and `hex_id` strings. Original `host` has 12,149 identified labels and 793 missing; **0** exact fully spelled binomial labels.
- Original `host` labels are mostly **bare milkweed epithets**, e.g., `syriaca` (3,638), `curassavica` (2,231), `incarnata` (1,896), `tuberosa` (1,228); other labels `Calotropis sp.` (156), `Cynanchum sp.` (252), `Gomphocarpus sp.` (217) are unresolved to species and excluded by the *same original exact-species gate*.
- Correct normalization: a single lowercase valid epithet in the original source nomenclature becomes `Asclepias <epithet>` using the original study convention; explicit other-family binomials remain as recorded; `sp.`, hybrid/complex/multispecies, missing or genus-only labels remain excluded. This step uses **no outcome, location, season or native/alien effect labels**.
- A dedicated regression test now asserts `curassavica → Asclepias curassavica`, and `complex speciosa/syriaca` and `Gomphocarpus sp.` remain excluded.

Do not interpret the initial zero result as a failed sample gate. Only the subsequent source-aligned run may be used to evaluate the unchanged pre-outcome information threshold.
