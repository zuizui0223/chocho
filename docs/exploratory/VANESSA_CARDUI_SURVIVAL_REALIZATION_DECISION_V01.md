# Vanessa cardui original experimental survival calibration: verified source, unreconciled endpoint

**Date:** 2026-10-08. **Status:** SOURCE VERIFIED / ADULT-FITNESS CALIBRATION ON HOLD / NOT A NEW ECOLOGICAL DISCOVERY. Analysis branch only; GEB manuscript unchanged.

## Frozen before inspecting this raw file

`VANESSA_CARDUI_SURVIVAL_REALIZATION_PROTOCOL_V01.json` was committed before the separate mirrored experiment was downloaded. It specifies a five-natural-plant source audit, source-aware observations and a stop rule when the original survival-event definition cannot be reconciled. This follow-up decision was written **after** seeing the raw record counts and checking the original paper. It is not preregistration.

## Why this case was selected

*Vanessa cardui* occurs in `data/frozen/s1_species_manifest.json` for the chocho butterfly analysis, and its 2024 feeding experiment has genuinely followed larvae rather than static plant geography.

Primary published source: Saldivar & Wilson-Rankin (2024), *Ecosphere* DOI [10.1002/ecs2.4810](https://doi.org/10.1002/ecs2.4810), original source DOI [10.5061/dryad.wdbrv15v5](https://doi.org/10.5061/dryad.wdbrv15v5). The same original archive can be accessed via [Zenodo record 10689711](https://zenodo.org/records/10689711), not a second experiment. Source `larv.plants.csv` verified MD5 `0af2cbb249494e11e8dc81a68dcf7bae` and **1,154** original rows in successful [GitHub source audit 37764390612](https://github.com/zuizui0223/chocho/actions/runs/37764390612).

## Raw data gate

| Natural putative host | Plant code | Rows | Published Table 1 mortality |
| --- | --- | ---: | ---: |
| *Abutilon palmeri* | ABPA | 215 | 100% |
| *Malacothamnus clementinus* | MACL | 240 | 100% |
| *M. fasciculatus* | MAFA | 272 | 99% |
| *Nicotiana glauca* | NIGL | 235 | 100% |
| *Sphaeralcea ambigua* | SPAM | 192 | 99% |

The archived natural-plant CSV has **1,154/1,154 records coded `status=2`**, and the archived README defines status=2 as death, status=1 as alive. But original article Table 1 reports 99%, rather than 100%, mortality for MAFA and SPAM, and the article reports a small number of successful pupations. These sources **do not furnish a reconciled individual adult-emergence indicator** in the inspected natural-plant CSV. We cannot infer that all larvae failed to pupate, or build an adult-emergence quality weight from this file.

**Critical original-Methods correction:** raw `source=oviposition` does not mean field-collected or wild larvae. Those eggs were laid by adult females that were reared from commercially sourced larvae. The article's 27 separately wild-collected larvae, discussed in its Table 2, are a **different cohort and not in the 1,154-row natural-plant archive**. Raw origin counts `carolina=869`, `oviposition=285` represent commercial eggs versus commercial-derived offspring, not commercial versus wild.

The archive can support narrow descriptive time-to-reported-death summaries, subject to its exact coding. It cannot calibrate *true adult recruitment* across the five hosts without the pupation/adult outcomes and source reconciliation. Original paper already reported the diet effect, near-total mortality of commercially sourced larvae and differences with field-collected larvae; no new physiological effect is claimed.

## Decision

1. **Access problem solved:** original source byte-exact mirror is reachable and reproducible. This is a reproducibility improvement, not a biological result.
2. **Demographic calibration gate blocked:** every archived natural-plant row has `status=2`; no source-confirmed per-individual emergence outcome is available in this CSV. The corrected script explicitly returns `adult_eclosion_estimate = null` and withholds survival curves.
3. **Biological inference bounded:** even if the two rare pupations are reconciled, a single commercial-rearing experiment in southern California cannot estimate worldwide realization of geographically available butterfly larval resources. Nursery versus wild-sourced plants and plant-trial groupings introduce additional confounding.
4. **Stop rather than rescue:** do not turn these already-published survival differences into a second ecology paper, infer a global fitness penalty for introduced hosts, or modify the GEB resource-opportunity manuscript. Future quantitative cross-system quality weighting requires independently measured adult emergence for multiple butterfly–host links across field-origin populations and matched sites.

**Main GEB paper remains unchanged.**
