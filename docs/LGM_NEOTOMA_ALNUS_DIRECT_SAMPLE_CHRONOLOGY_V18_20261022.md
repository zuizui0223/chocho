# LGM Neotoma: six full original Alnus glutinosa pollen datasets, direct sample-age triage (v18)

**2026-10-22; verified source metadata, NOT independently calibrated fossil presence or ancestral host occupation.**

## Why this new source route was needed

The official Neotoma numeric taxon-ID + age-filtered occurrence search produced 27 metadata rows from six datasets (26 distinct `occid`, PR #63). The source dataset and chronology/units audit (PR #64) identified four numeric age discrepancies with source dataset envelopes and mixed calibrated/non-calibrated years BP. In PR #66, all SIX original complete source dataset downloads were obtained with each body SHA pinned. In the ensuing v17 one-occurrence-per-dataset crosswalk, the official individual occurrence records did **not** expose `sampleid` or `analysisunitid` for an exact join, and a deliberately impossible occurrenceid was echoed in a successful response. Therefore **do not join the earlier occurrence metadata to sample ages by matching calendar age, depth or row order**.

For this separate exploratory source audit, six full Neotoma dataset source bodies were replayed with exact source SHA256s (v16 artifact #12427662101) and each original `site.collectionunit.dataset.samples[]` was inspected directly. Dataset IDs 4396, 14556, 22662, 22683, 24197, 41351 and site IDs 3392, 9834, 14474, 14485, 15750, 23855 stayed fixed. Taxon ID **1299** is the database label *Alnus glutinosa*. Count as a source-positive pollen DATUM only a finite numeric `samples[].datum[].value > 0` with `datum[].taxonid == 1299`. Keep every associated `samples[].ages[]` record with original `chronologyid`, `agetype`, numeric age and bounds, rather than selecting whichever of multiple models falls in LGM.

## Actual source result

- Full datasets contain **557 total source samples**, of which **437 distinct dataset/sampleid entries** have a positive numeric taxonid1299 pollen datum (not necessarily 437 independent local populations).
- Those 437 sample entities include **788 separate source age-model records**, often *two age models for the same sample*. This is NOT 788 independent fossil occurrences.
- Under the SOURCE-LABEL-ONLY age criterion `agetype='Calibrated radiocarbon years BP'` AND `age.age in [18000,24000]` AND explicit `chronologyid`, **21 distinct source sample identities / 21 age-model rows** appear across **five** of the six sites. This is stronger metadata linkage than the older generic occurrence search, but still does **not verify age calibration independently**.
- Candidate source sample counts:
  - Sheheree Bog, dataset 14556 / source chronologyid **7725**: **9 source samples**, 18,564–23,611 source BP. All NINE lack older/younger bound fields. Note: three of these dated values 22,349–23,611 exceed that dataset's original advertised upper 21,974 cal BP, as already flagged in v15; keep as chronology-review conflicts.
  - Rifugio Mondovi, dataset 22662 / chronologyid **13010**: **1** candidate at 18,241 BP, missing age bounds and outside published maximum 16,079 cal BP in source v15.
  - Laghi dell'Orgials, dataset 22683 / chronologyid **13034**: **1** at 20,368 BP with source `ageyounger=20,044`, `ageolder=20,717`. The only candidate with both bounds in the direct source data.
  - Dar Fatma, dataset 24197 / chronologyid **14552**: **3** candidates ages 19,227–23,260 BP, no source bounds.
  - Saham Mere, dataset 41351 / chronologyid **25462**: **7** candidates ages 18,777–21,044 BP, no source bounds.
  - Plesheevo Lake, dataset 4396 / chronologyid **30953**, has **86** taxon1299-positive sample datums and associated `Radiocarbon years BP` age rows but **ZERO candidates meeting explicitly source-calibrated-BP label**. Previously flagged chronology very provisional, no arbitrary conversion between 14C age and calibrated BP.

**20/21** sample-labelled 18–24ka age candidates lack BOTH numerical age-uncertainty bounds. A source agtype string does not by itself certify dating controls, calibration curve validity or a species-diagnostic pollen morphotype. Pollen can be deposited by transported grains rather than local plant populations.

## Robust status and code

[Source-bound v18 GitHub Actions #39785217854](https://github.com/zuizui0223/chocho/actions/runs/39785217854) **SUCCESS**, with regression tests, six official original dataset SHA verifications, all positive/zero/invalid source `datum` classifications, exact `datasetid+sampleid` uniqueness checks and chronological record source-grain retention.

Durable source results:
- `results/lgm_neotoma_Alnus_direct_sample_age_models_v18/summary.json`
- `results/lgm_neotoma_Alnus_direct_sample_age_models_v18/all_source_Alnus_sample_datum_v18.csv`
- `results/lgm_neotoma_Alnus_direct_sample_age_models_v18/all_source_Alnus_sample_age_models_v18.csv`
- `results/lgm_neotoma_Alnus_direct_sample_age_models_v18/six_source_dataset_pollen_age_inventory_v18.csv`
- `results/lgm_neotoma_Alnus_direct_sample_age_models_v18/receipt.json`.
- Original source v16 raw download ZIP remains separately source-pinned with file/body digests, and does not get silently rewritten.

**Scientific inference boundary:** confirmed 21ka occupied *Alnus glutinosa* plant-population sites = **zero**. No fossil source plant distribution reconstructed, no insect genetic response opened or hypothesis test reopened. TTF #149 and chocho #59 modern-climate-matched OTHER REAL PLANT **negative** results remain negative. These new source-reconcilable direct sample candidates are potential inputs for a future calibrated age-depth/pollen species-identification verification, not evidence of larval resource refugia yet.
