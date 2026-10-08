# The `Euphydryas editha` four-host local-to-global crosswalk

**Date:** 2026-10-08. **Status:** SOURCE-VERIFIED INPUT AUDIT / POST-HOC ONE-SPECIES STRUCTURAL RESULT / NOT A NEW LOCAL DEMOGRAPHIC EFFECT. **Scope:** analysis branch only. The GEB submission manuscript and PR #38 are unchanged.

## Reproduced fixed-panel baseline

GitHub Actions [run 37767325222](https://github.com/zuizui0223/chocho/actions/runs/37767325222) completed successfully against the pinned original HOSTS (`808e0b869f9ec1adf8efff87cf6a395adda103e0`) and WCVP (`65bed76bae9d644ccb6ad200c05f9f5071d89e05`) identities. `Euphydryas editha` had 66 original HOSTS consumer rows mapping to **30** exact accepted plant species, with **131 native** and **245 contemporary** potential resource WGSRPD3 regions, and **114 introduced-added** regions. This reproduces the frozen 239-species figure source exactly.

## A source-backed taxonomic completeness gap

Four exact local hosts documented in Singer & Parmesan (2018, *Nature*, doi:10.1038/s41586-018-0074-6) or Haan, Bowers & Bakker (2021, *Scientific Reports*, doi:10.1038/s41598-020-80413-y) were checked *before* reading exact accepted global-link outcomes. Two are present and two absent in the fixed *E. editha* interaction inventory:

| Local documented host | Exact raw HOSTS *E. editha* rows | Resolved in frozen WCVP-linked host inventory | Native / contemporary botanical WGSRPD3 |
| --- | ---: | --- | ---: |
| *Plantago lanceolata* (introduced host) | 2 | Yes, accepted ID 2569834 | 79 / 218 |
| *Collinsia parviflora* (ancestral host in Nevada) | 2 | Yes, accepted ID 2731115 | 26 / 27 |
| *Castilleja hispida* (native host in Taylor's checkerspot study) | 0 | **No exact butterfly link** | not assessed by prior fixed interaction crosswalk |
| *Castilleja levisecta* (native host in Taylor's checkerspot study) | 0 | **No exact butterfly link** | not assessed by prior fixed interaction crosswalk |

The missing Castilleja links are independently documented local biological hosts but absent from the exact species-wide pair table. Their missingness is an **interaction-knowledge completeness problem**, not evidence of evolutionary loss of function or local botanical absence. Since these are four deliberately selected species from two papers, the fraction 2/4 is **not an estimator of global HOSTS coverage**. The Nevada and Washington studies sample different subspecies/populations, habitats and history.

## Contribution of mapped introduced `Plantago lanceolata` to this butterfly

Its WCVP accepted botanical geography is **79 native / 218 contemporary** WGSRPD3 regions, hence **139 introduced-only botanical regions**. Those 139 botanical regions are not automatically 139 new butterfly-resource regions, because other known hosts can occur in the same cells.

Counterfactually excluding only *Plantago*'s **introduced** regions and retaining all 30 original known host identities and botanical native distributions removes **52/114 = 45.6%** of *E. editha*'s reconstructed introduced-added resource regions. The other 87 *Plantago* introduced-only botanical regions were already covered by at least one of the other 29 mapped host plants (or by native resources), so they are not uniquely attributable to it.

**Ecological limit:** 52 is **structural introduced-range reliance within the current database**, not the number of habitats occupied by a Plantago-dependent butterfly population, not field evidence of exclusive larval use, and not the population loss after plant removal. There is no direct causal link to the local Nevada extinction or Washington host choice data.

## Independently specified follow-up

The post-hoc follow-up protocol `EUPHYDRYAS_MISSING_CASTILLEJA_AUGMENTATION_PROTOCOL_V01.json` adds **exactly two** independently published Castilleja hosts to the original *E. editha* host list after exact accepted WCVP taxon and botanical-range verification. It then recomputes native, contemporary, introduced-added and Plantago-only resource regions. See `scripts/build_euphydryas_missing_hosts_wcvp.R` and `scripts/analyze_euphydryas_missing_hosts_augmentation.py`, executed via `.github/workflows/butterfly-euphydryas-missing-host-sensitivity.yml`.

A subtle but crucial diagnostic is **reclassification**: adding an omitted native host can make a region previously coded “introduced-added” become a native-resource region. Thus the introduced-added count can **decrease even as the total contemporary number of mapped host-resource regions increases**. This is a baseline-completeness property, not evidence that introduced hosts have harmed a butterfly.

**Stop/decision rule:** treat the update as a one-species, source-aware completeness sensitivity. Do not change the primary GEB panel, claim causal population recovery/extinction, or present the four handpicked hosts as representative of all butterflies. No extra host species selected based on outcomes.
