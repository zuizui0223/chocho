# BCE/Clarke additional accepted hosts: geography sensitivity and scientific decision

**Date:** 2026-10-08. **Status:** SOURCE-REPRODUCED EXPLORATORY UPPER-BOUND EUROPEAN HOST AUGMENTATION / GLOBAL ORIGINAL MANUSCRIPT NOT CHANGED. **No novel realized population-fitness result.**

## Source-audited chain

1. Fixed *chocho* core: original global 239 butterflies, 26,530 native butterfly × WGSRPD3 resource units, 41,083 contemporary and 14,553 introduced-added (54.8549%, 206/239 expanded). These totals were **exactly reproduced** in the current original HOSTS x WCVP reconstruction before augmentation.
2. Independent bibliography-derived *candidate* input: BCE website states its butterfly foodplants derive from Clarke (2024, doi:[10.1002/ece3.10834](https://doi.org/10.1002/ece3.10834)) and retain host evidence ranks 1–3. The original 2024 Dryad relational files were not available in CI; BCE is a **secondary** evidence-filtered list with incomplete per-link provenance, not a new field dataset. The exact overlap with the frozen chocho 239 species is 92 BCE butterfly binomials, of which 83 have a species-level foodplant page; original source retrieval [37777242991](https://github.com/zuizui0223/chocho/actions/runs/37777242991).
3. Fixed WCVP accepted taxon ID reconciliation [37777881774](https://github.com/zuizui0223/chocho/actions/runs/37777881774): 1,052 BCE-only exact botanical spellings resolve to 25 already-existing accepted-ID host pairs and **1,027 candidate accepted-ID host links absent** from frozen HOSTS. The 1,027 accepted-ID edge candidates span **679 distinct accepted botanical species IDs**, and affect **81 butterflies**. Candidate associations are literature checklists, not 1,027 independent witnessed feeding events.
4. Follow-up sensitivity was **frozen before seeing the accepted-ID gap counts** in `BCE_CLARKE_ONE_DIRECTION_HOST_AUGMENTATION_PROTOCOL_V01.json`; source-first script `analyze_bce_clarke_added_hosts_geography.py` runs exact 239 baselines and only adds the accepted-ID-absent BCE candidate edges, without modifying other species or altering the original source data.
5. [Successful source-audited GitHub Actions run 37778961442](https://github.com/zuizui0223/chocho/actions/runs/37778961442), artifact `chocho-bce-clarke-accepted-host-europe-geography-v01` (ID 11551469909) has the full 239 species-by-species sensitivity, source identities, and rerun receipt. Its WCVP additional botanical range builder found native/contemporary records for **all 679** plant IDs. Source snapshot: HOSTS commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`; WCVP commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`.

## Results: different estimands change in different directions

| Conditional resource metric, across original 239-species panel | Frozen HOSTS | One-direction BCE accepted-host augmentation | Change |
| --- | ---: | ---: | ---: |
| Native butterfly × WGSRPD3 resource units | 26,530 | 29,852 | **+3,322** |
| Contemporary butterfly × WGSRPD3 resource units | 41,083 | 44,819 | **+3,736** |
| Contemporary minus native (added) resource units | 14,553 | 14,967 | **+414** |
| Introduced-added / native ratio | 54.85% | 50.14% | **−4.72 percentage points** |
| Butterflies with any additional resource region | 206 | 213 | **+7** |

The exact 83-butterfly source-verifiable European host-page subset (not an unbiased representative sample of the full 239) has original native `9,225` → augmented `12,547`, and original introduced-added `4,796` → augmented `5,210`. Thus the same +3,322 native and +414 added units come from the European exact source-selected subset; the remaining 156 globally analyzed butterflies are held unchanged in this sensitivity.

**Central ecological/data interpretation:** supplementary documented butterfly–host links do not merely add new extra-regional resource opportunities. They also *repair a larger baseline of native potential host geography*. Under this one-direction augmentation, geographic resource expansion persists, and its absolute increment grows; but its **relative increment decreases** because the native baseline grows even faster. This is a real measured **dataset-conditioning sensitivity**, not a newly discovered mechanism of species colonization or ecological host quality. The sign of a strictly nested current-minus-native added-resource envelope is structurally nonnegative; the interesting uncertainty concerns its size, source-dependent spatial arrangement, and how many butterflies actually use the hosts.

## Crucial scope, reasons NOT to revise the main paper yet

- **Upper-bound extrapolation:** BCE's European records were hypothetically used as species-wide host associations at all botanical WCVP3 regions worldwide. No local butterfly feeding observations or recruitment data validate those global extrapolations. Using both published host compendia can create different sampling/ascertainment biases.
- **Geographic source selection:** only the 83 BCE source-visible European butterflies received additional host links, of 239 frozen butterflies. This is a **partial European-species sensitivity**, not a corrected global effect and not a representative sample-based confidence interval.
- **Evidence strength:** BCE retains Clarke ranks 1–3, but 3 is oviposition/eggs, not proof of larval survival to adulthood. Some host records may be repeated across the two literature databases.
- **No false-negative rate:** 1,027 is a source-conditional missing-*accepted-ID pair* count; it is not an estimate of the number of currently viable locally accessible hosts, and cannot be expressed as global HOSTS sensitivity from the nonrepresentative 83 species.
- **Main community claim unresolved:** the secondary geography sensitivity above checks **resource breadth**, not the principal fixed-margin-null comparison for 355-region biotic homogenization or butterfly-pair shared-host overlap. These need a separate analysis with frozen candidate identity and counterfactual regions; do not infer their exact magnitudes from the breadth shift.
- **No demographic inference:** resource opportunity is not abundance, colonization, competition or survival. The original occurrence validation classifications also depend on the adopted native HOSTS list and require a dedicated host-inventory sensitivity if promoted.
- **Keep original:** preserve the original GEB 239-species `+54.9%` estimate as precisely *conditional on the fixed HOSTS/WCVP reconstruction*. Do not silently replace it with `50.1%` or treat this as corroborated truth. State the conditionality and investigate reproducibly.

## Scientific decision

The source mismatch is too large to dismiss as a spelling issue: **1,027 accepted-ID candidate host links**, with a measurable **4.72 percentage-point sensitivity of proportional expansion** in a one-direction European host augmentation. That sensitivity is not a refutation of the original reconstruction but strengthens the need for source-quality and geography-aware reviewer defense. The next decisive structural question is whether **the original fixed-margin-excess resource homogenization**, rather than just absolute resource expansion, persists when the accepted-ID matched BCE host candidates are added under the same source scope. The follow-up must hold its Europe-derived sample fixed, not select favorable butterflies, and report null/reference sensitivity openly.

**No change to PR #38, current main-paper submission package, or its successful CI. No independent ecology paper promoted from these host-knowledge sensitivity analyses.**
