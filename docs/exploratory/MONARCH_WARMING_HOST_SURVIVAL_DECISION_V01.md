# Monarch host × warming × parasite adult-survival reality check

**Date:** 2026-10-08. **Status:** PUBLISHED SOURCE RE-EXTRACTION; NOT AN INDEPENDENT ECOLOGICAL DISCOVERY. Experimental outcome values are transcribed from a pinned original author CSV and grouped deterministically; independent GitHub Actions reproducer may be inspected at `.github/workflows/butterfly-monarch-warming-fate-source.yml`. **No modification to the GEB manuscript or PR #38.**

## Original provenance and scientific status

Ragonese et al. (2025), *Ecological Entomology*, doi:[10.1111/een.70010](https://doi.org/10.1111/een.70010); archived source doi:[10.5061/dryad.prr4xgxx7](https://doi.org/10.5061/dryad.prr4xgxx7), original public author GitHub [IRagonese/MilkweedWarming2024](https://github.com/IRagonese/MilkweedWarming2024), source file `MWwarming_comp_May16.csv` pinned commit `74e3e2cd8079401c9ccd4359f9e39e2d52d1f7d2`. Published work already tested host, warming, parasitism and their interaction. This is a **descriptive audit of originally published individual records**, not a new source or a newly randomized experiment.

We viewed the published study abstract and first source rows before freezing the following simple 2×2×2 extraction. This is explicitly exploratory and nonconfirmatory.

## Adult survival, with every treatment arm retained

The original CSV contains **240 individual larvae**, 30 randomly assigned to each listed host × temperature × inoculation group. Of 240, **239 `Surv_adult` fields** contain binary 0/1; the other is missing in ambient, tropical, parasite-free controls. The table below uses all measured adult survival outcomes without deleting observed deaths. Survival is adult eclosion during the experiment, not survival of a free-living population.

| OE inoculation | Temperature | Tropical milkweed survived/observed | Swamp milkweed survived/observed |
| --- | --- | ---: | ---: |
| Control | ambient | **26/29 (89.7%)** | **28/30 (93.3%)** |
| Control | elevated | **23/30 (76.7%)** | **29/30 (96.7%)** |
| Infected/exposed | ambient | **28/30 (93.3%)** | **27/30 (90.0%)** |
| Infected/exposed | elevated | **28/30 (93.3%)** | **30/30 (100.0%)** |

In the **no-inoculation control** the tropical-minus-swamp raw adult survival difference was approximately **−3.7 percentage points at ambient** and **−20.0 points at elevated** temperature; their contrast was **−16.3 percentage points**. In inoculated groups, the host difference was **+3.3 points at ambient** and **−6.7 points at elevated**, contrast **−10.0 points**. These are unadjusted within-study arithmetic summaries; no interaction p-value, uncertainty interval, or superiority assertion is warranted without the randomized plant/plot/lineage allocation and its repeated larvae accounted for. Exposure (`OE_treatment=infected`) must **not** be equated with verified successful infection; conditioning on `InfectionStatus` after treatment would bias the design. The source's `Plant_ID` and `PlotNum` labels support reconstructing clustering; two larvae can share a plant. The authors' original model and publication remain the primary evidence for treatment-effect inference.

## What this does and does not test in chocho

*Danaus plexippus* is present among the 239 chocho butterfly species with **40 accepted host species**, **181 native** and **263 contemporary** potential host-resource WGSRPD3 regions, and **82 added resource regions**. This 2021 field-cage factorial experiment measures actual survival but does **not** sample those 82 geographically added regions as a colonization cohort. Its observational unit, treatments, host species, parasite lineage and life-history stage are not exchangeable with the frozen global opportunity unit. It **cannot** supply a unique conversion factor from a mapped host region to demographic fitness.

The study already found that temperature influences infection tolerance and that host plant identity affects development and survival; re-extracting the published interaction is not novel. The most defensible cross-scale conclusion is a **measurement boundary:** botanical host range is necessary for potential local host supply in the model, but is insufficient to infer larval eclosion probability in the absence of thermal and infection context. Such a discrepancy can be important even when `known-host × region` is unchanged.

**Stop:** Do not append hypothetical fitness weights to the 239-species resource map, claim a new host-quality mechanism, replace the original published experimental model with unclustered row tests, or merge this extraction into the GEB submission.
