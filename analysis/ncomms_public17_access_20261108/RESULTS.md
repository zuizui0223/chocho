# Nature Communications exploratory gate: 17 original islands, 225 butterflies, ocean circuits versus prior-record source distance

**2026-11-08 | source-locked quantitative analysis COMPLETED | Not a colonisation, shipping, or prior-HOSTS interaction test.**

## Data, scope and what was fitted

This new analysis is a re-use of a previously *independently queried*, **unselected for future first records** public 2026 GBIF index panel (original [October 31 run](https://github.com/zuizui0223/chocho/actions/runs/41196133727)) and the existing full **51-island source-reconstructed marine effective resistance matrix**. A 2026 live GBIF search index is not an independent biological census or the original private **83,573-row** annual raw Parquet. The 17 island groups were themselves selected by an earlier *near-island* study design. They are not a sample of all 51, let alone all global islands.

The analysis includes **all 225 original butterfly candidate species** at each of these 17 original WPI-associated GSHHG landmass IDs (3,825 initial possibilities). We define nonindexing strictly as a **GBIF recording outcome**, never a true ecological absence. In each time window we evaluate only butterflies not previously indexed on the target island; source distances use only prior-period records, never future records. The original period cohorts are:

| Retrospective holdout | At-risk record-index cells | First-indexed response | Prior source-index available among eligible cells | New records among those with a prior indexed source |
|---|---:|---:|---:|---:|
| 2011–2017 (source up to 2010) | **3,584** | **80** | 864 | 57 |
| 2018–2022 (source up to 2017) | **3,504** | **18** | 1,039 | 17 |

**Four fixed models**, one standard 2026 GBIF recording risk rather than biological colonisation:
- **B**: baseline with pre-window island-wide GBIF butterfly indexing effort, this species' prior indexed island count and historical record count, presence of any indexed source, original GSHHG island area, modern WPI port count, latitude/longitude, and 2001–2010 *general* plant recording richness as botanical-effort proxy.
- **D**: B + `log(1 + nearest previously indexed source island's great-circle port-centroid distance)`.
- **C**: B + `log(1 + minimum source-indexed uniform-ocean effective resistance)`.
- **DC**: B + D + C. The **increment beyond direct distance** is specifically D → DC. C contains NO ship positions or time direction, and source islands are not verified breeding exporters.

Train-only `StandardScaler` + fixed L2 logistic regression `C=1`, and **three independent refoldings of the same outcomes**: 17 source islands held out individually, six original connected geographic near-pair components held out, or 5 fixed species-key folds. All models are fit without test-label data. There is **no new climate or historical host-plant accepted-name panel in this comparison**, so these are observational-recording predictions only.

## Main source-verified heldout results

Values below are **percentage reduction in out-of-fold Brier loss** compared with the indicated reference, using the complete same target rows. Positive means lower prediction error; negative means it worsened. All results report both Brier and heldout log loss in the generated receipt.

| Period | Geographic holdout | Cases | B→D: nearest recorded-source distance | B→C: ocean resistance | **D→DC: ocean resistance beyond distance** |
|---|---|---:|---:|---:|---:|
| 2011–2017 | One **island** unseen (17 folds) | 80 | **+21.38%** | **+12.14%** | **−1.88%** |
| 2011–2017 | One **near-island component** unseen (6 folds) | 80 | **−75.79%** | **−83.14%** | **−16.19%** |
| 2011–2017 | **Species** unseen (5 folds) | 80 | **+22.81%** | **+12.34%** | **−1.06%** |
| 2018–2022 | One **island** unseen (17 folds) | 18 | **+6.09%** | **+4.55%** | **+0.08%** |
| 2018–2022 | One **near-island component** unseen (6 folds) | 18 | **+0.13%** | **−2.46%** | **−0.82%** |
| 2018–2022 | **Species** unseen (5 folds) | 18 | **+4.18%** | **+5.29%** | **+0.71%** |

**Independent log-loss guard:** Adding the circuit term to a model that *already contains direct source distance* (D→DC) worsens log loss in **all six** comparisons: respectively **−4.01%, −17.16%, −1.23%** for 2011–2017, and **−1.10%, −2.16%, −0.32%** for 2018–2022. The two small positive later-period Brier increments are not corroborated by this proper scoring metric. No source-invariant additional circuit information beyond direct distance is established.

The positive **distance-only** increments in held-one-island and held-species comparisons are genuine source-backed retrospective associations at their fixed risk sets. **However**, 2011–2017 performance collapses when one of six original geographic near-island components is excluded: the distance-only Brier loss increases from **0.03491** to **0.06136** (−75.79% gain). This negative geographic transfer test is mandatory to report; the strongest adverse components are original components 1 and 2. This prevents promotion of within-archipelago predictive gains to a universal island-dispersal law.

### Numbers with original units

2011–2017 **held-island** Brier: B **0.02153626**, D **0.01693214**, C **0.01892228**, DC **0.01725124**. Held-**component** Brier: B **0.03490556**, D **0.06136047**, C **0.06392443**, DC **0.07129598**.

2018–2022 **held-island** Brier: B **0.00511419**, D **0.00480268**, C **0.00488139**, DC **0.00479868**. Held-**component** Brier: B **0.00565651**, D **0.00564917**, C **0.00579555**, DC **0.00569567**.

The tiny later event count (18 across 3,504 risk cells), near-neighbor selection, species/island dependence and survey effort create substantial uncertainty. No independent confidence intervals or causal source attribution are inferred from fold-level counts.

## Relationship to the existing ecological positive evidence

**Preserve** the separate demonstrated outputs at their own grains:
- 239-species WGSRPD3 introduced resource opportunity increases **+118.9% on islands vs +40.9% mainland**, area/biogeography-matched island excess **+10.7 pp**.
- 98 region modern WPI port count associates with introduced **whole flora β+0.2445** and distinct introduced hosts **β+0.1639**, with no positive residual host-enrichment coefficient given full introduced flora **β−0.0415**.
- Original 51-pair geography and 225 candidate hosts have retrospective host–later butterfly similarity **Spearman ρ+0.300**; previous island-held source similarity models report exploratory positive conditional host and marine increments. They predict a **different outcome** on a **different GBIF source**, so this new binary first-index risk test does **not** negate those scores.
- Separate global 155-butterfly time-held and spatial-held host-based recording prediction improved species-equal MSE **+4.12% / +3.00%**; static host covariates cannot document actual prior plant arrival.
- The independent 35-butterfly / six-host-plant / 17-island *pre-2011 named-host* subset was separately source-recovered from [October 31 independent run](https://github.com/zuizui0223/chocho/actions/runs/41196373882). It contains **42 true HOSTS edges** and 595 source cells but was **posthoc selected**. Its earlier source-audited botanical conditional null did not establish host-specific exceptional first indexing. Do not treat this as a representative 225-species HOSTS panel.

## Science and inference limits

1. A prior GBIF record on a second island is a *potential source index*, **not proof of a breeding colony**, immigration vector or true population movement. A first indexed record is not verified first arrival.
2. Uniform-ocean effective resistance is an undirected modelled over-water circuit that does not allow island-stopover movement; modern indexed WPI ports and 0.5° ocean geometry are not reconstructed historical voyages.
3. The source17 cohort was chosen for geographic neighbor pairs, the 2026 GBIF re-query differs from the original private fixed 2001–2024 human observation source, and species/regions do not form independent replications. There is no climate-matched source panel.
4. **Most importantly**, this test has **no dated authentic larval host×island source predictor**. It is deliberately a separate access-layer sanity check; it **does not test the NC host × accessibility interaction**.
5. The original public 14.6 million-record GBIF DOI [10.15468/dl.umqmsf](https://doi.org/10.15468/dl.umqmsf) and the original GSHHG2.3.7 high-res source archives are known and source-hashed (Oct31 provenance), but the author's 83,573-row derived Parquet remains raw-byte-inaccessible in the current runtime. Rebuilding all 51 source islands from pinned original public occurrences is a separate source-reproducibility task, not an outcome inferred here.

## Reproduction and CI

[Source-locked script](../../scripts/run_ncomms_public17_ocean_access_holdout.py) and [automated original-cohort tests](../../tests/test_ncomms_public17_ocean_access.py) use original Git blob pins for all eight inputs. Three original public source files (225×17 panel, original query receipt and original GSHHG area frame) were byte-exactly copied into `data/exploratory/ncomms_public17_access_20261108/`; they were already on the prior Oct31 source branch and are not a newly fabricated GBIF query. The original uniform circuit matrix remains in `results/island51_observed_ocean_circuit_20261106/`. Every one of the six model×holdout combinations includes complete out-of-fold per-cell predictions and a source-locked numerical summary.

[Successful source-locked CI](https://github.com/zuizui0223/chocho/actions/runs/42458266523). Its artifact `ncomms-public17-source-access-heldout-all225` contains **six OOF prediction CSVs and JSON source receipt**. The CI checks all 80/18 cases, all 3,584/3,504 risk rows, 17/6/5 independent holdout folds and all four prediction columns, plus frozen SHA1 of source observations.

**Scientific decision:** Promote prior-record-source proximity as an *exploratory geography association*, NOT as a proven or geographically transferable dispersal pathway; **do not claim uniform-ocean circuits offer additional independent prediction beyond geographic distance**. The main Nature Communications resource×access hypothesis remains untested pending historical accepted HOSTS identities and more geographically representative risk sets.
