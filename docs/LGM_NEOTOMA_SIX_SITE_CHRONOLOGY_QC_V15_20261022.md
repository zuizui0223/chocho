# 6 Neotoma Alnus pollen datasets: 18–24ka source age-scale and chronology QC

**Date:** 22 October 2026. **Status:** verified source metadata contradictions/uncertainty, **zero confirmed independent LGM Alnus glutinosa occupied site records**.

This diagnostic concerns the **27 records (26 distinct occurrence IDs) from six Neotoma dataset/site IDs** returned for numerical Neotoma taxon ID 1299 (*Alnus glutinosa*), pollen `datasettype`, and numerical sample central age `18,000–24,000 BP`. It is not a complete fossil survey, not proof of pollen morphology diagnostic to Alnus species, and not an insect genetic result.

## Exact source chain

1. Source [PR #61](https://github.com/zuizui0223/chocho/pull/61), [Actions #39722583906](https://github.com/zuizui0223/chocho/actions/runs/39722583906): generic `/occurrences?taxonname=...` returned identical source IDs to **a deliberately nonexistent taxon**. Its old results must never be treated as species-specific fossil data. The official taxon registry independently resolves Neotoma ID1299.
2. Direct-ID source [PR #62](https://github.com/zuizui0223/chocho/pull/62): `/taxa/1299/occurrences` gives distinct IDs from taxon3461, but first pages don't establish LGM absence.
3. Source bounded numeric-ID + exact pollen + 18–24ka age-value query [PR #63](https://github.com/zuizui0223/chocho/pull/63), [Actions #39723924751](https://github.com/zuizui0223/chocho/actions/runs/39723924751): **27 candidate rows /26 unique occurrence IDs /6 datasets**; original modelled/calibrated chronology uncertainty and pollen identification not validated.
4. Official source six `/datasets/{id}`, `/chronologies`, `/publications` metadata [Actions #39724671797](https://github.com/zuizui0223/chocho/actions/runs/39724671797) all 6×3 endpoints available, with source-citation snippets. Their chronology metadata are **not necessarily the exact age model assigned to each occurrence**.
5. Latest **independent original-byte-SHA source-only audit** [Actions #39725710005](https://github.com/zuizui0223/chocho/actions/runs/39725710005) **SUCCESS**: `scripts/audit_neotoma_six_dataset_age_scale_v15.py`, tests, 27-row source-level chronology flags and six-dataset citation matrix committed. Exact original source candidate + metadata SHA guards passed. No fossil occupancy admitted.

## Dataset-level age source flags

| Neotoma dataset | Site | Candidate rows / unique IDs | Published dataset age scale | Numerical source age issue / documentation |
|---|---|---:|---|---|
| 14556 | Sheheree Bog | **9 /9** | Calibrated radiocarbon years BP (upper bound 21,974) | **3** source occurrence central ages 22,349, 22,980, 23,611 are *above* the advertised dataset oldest numeric age. Whether a different chronology was used is unknown |
| 22662 | Rifugio Mondovi | **1 /1** | Calibrated radiocarbon years BP (upper bound 16,079) | **1** candidate central age 18,241 is numerically beyond the dataset age envelope. Original age notes cite IntCal04 calibration, model assignment unverified |
| 4396 | Plesheevo Lake | **3 /3** | Radiocarbon years BP | Source chronology notes explicitly: **“Chronology is very provisional. Very old at the base”**. This is NOT a validated 21ka calibrated age |
| 22683 | Laghi dell’Orgials | **1 /1** | Calibrated radiocarbon years BP | Numeric candidate lies within published broad envelope, but source chronologies explicitly cite the older **IntCal04** calibration; current re-calibration/source confidence not reviewed |
| 24197 | Dar Fatma | **6 /5** | **Both** radiocarbon years BP and calibrated radiocarbon years BP listed | Mixed scales and **one duplicate occurrence ID**. Do not double count the occurrence or pool dates as though already calibrated |
| 41351 | Saham Mere | **7 /7** | **Both** radiocarbon years BP and calibrated radiocarbon years BP listed | Mixed scales. Original reference Bennett 1988; correct sample-assigned calibrated age model not yet reviewed |

Thus **four of the 27 numeric central candidate age rows** fall outside **all listed dataset age envelopes** (3 Sheheree Bog, 1 Rifugio Mondovi). **This is a metadata discordance FLAG, not proof the fossil is misdated or absent**: source occurrence `age`, the dataset's reported age extent, the first chronology metadata entry and calibrations need not be on the same timescale or using the same chronology.

In all six datasets, original source chronology and publication metadata exist, but **none of the 27 sample-level calibrated 18–24ka chronological assignments** has yet been independently verified against the original age-depth controls and source publication. Twenty-six of 27 candidate rows lacked the occurrence-level `ageolder/ageyounger` uncertainty bounds in earlier source retrieval. No accepted species-level pollen morphology verification has occurred.

## What this can—and cannot—say

**Confirmed methodological finding:** A broad fossil API query can simultaneously return an apparently plausible species+age match and hide an invalid name filter, mixed calibration units, data-set/source age inconsistency and duplicate occurrences. Source species+date metadata are not automatically ancient occupied plant populations.

**Not confirmed:** Actual *Alnus glutinosa* occupation at six glacial sites, glacial resource connectivity, host-specific insect refugia or insect genetic breaks. The previous TTF [#149](https://github.com/zuizui0223/TTF/pull/149) and chocho [#59](https://github.com/zuizui0223/chocho/pull/59) climate-matched alternative-plant host-identity **negative** results stay frozen negative.

**Only defensible next biological step:** for each original occurrence ID (starting with the four flagged ages), retrieve its actual sample/depth and attached chronology version and age-depth controls, resolve whether `age` is calibrated or raw radiocarbon BP, review source publications/pollen morphotype resolution, and ensure pollen transport/source site geology make local species-level ancient occurrence interpretable. Until then, record zero fossil sites admitted into the true 21ka host-range model.

## Original source-reproducible output

- `results/lgm_neotoma_six_Alnus_source_age_scale_qc_v15/all_27_candidate_source_age_scale_flags_v15.csv`
- `results/lgm_neotoma_six_Alnus_source_age_scale_qc_v15/six_dataset_source_chronology_risk_v15.csv`
- `results/lgm_neotoma_six_Alnus_source_age_scale_qc_v15/summary.json` and `receipt.json`

The source success means these QC tables accurately reproduce the original Neotoma supplied evidence metadata; it does **not** validate any ancient organism location.
