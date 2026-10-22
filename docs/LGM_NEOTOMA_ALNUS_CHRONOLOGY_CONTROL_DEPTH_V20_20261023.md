# Neotoma LGM *Alnus* pollen — source chronology reliability and dating-controls review (v19/v20)

**23 October 2026. Source metadata-only result.** Zero independently established local 21ka *Alnus glutinosa* occupied populations, zero independent species-diagnostic pollen confirmations and zero genetic responses used. The previous frozen TTF v0.5 and chocho v0.8 real-alternative-plant **negative** host-specific climatic-history results remain unchanged.

## Objective and original evidence

The original six Neotoma official full pollen dataset downloads contain 437 distinct samples with positive numeric source `taxonid=1299` (*Alnus glutinosa* **database label**) and 788 associated sample age-model records, often multiple chronologies per sample. Earlier v18 identified **21 distinct positive-pollen samples** with a *source-labelled* `agetype=Calibrated radiocarbon years BP`, `age=18,000–24,000` and original explicit `chronologyid`, across five datasets. This label alone does not establish a calibrated age, local occupied trees, or larval resources.

The v19 source-only audit fetched all **11** original chronology IDs and **6** dataset chronology lists from the official Neotoma v2 REST API, source ZIP SHA256 `f01683540e34ba722ba65ed484c88139a323b23d6427da58a50be63f2e526650`, run [#39862870638](https://github.com/zuizui0223/chocho/actions/runs/39862870638) **SUCCESS**. Every source chronology ID and pollen-dataset linkage was explicitly verified. Unlike a naive scanner, the true API schema is `data[].chronology.reliableagespan.{younger,older}`, `data[].chronology.controls[]`, `modelType`, `default`, `agetype` and `datasets`.

## Important Neotoma API fragmentation correction

**The official chronology endpoint can return the SAME chronology as multiple `data[]` fragments containing different age controls.** A preliminary v20 calculation reading only `data[0]` incorrectly claimed 18 of 21 candidates were deeper than the deepest dated-control depth. After all repeated records were merged by source `chroncontrolid + depth + age`, the accurate count is **11/21 beyond the deepest returned control** and **10/21 inside the returned control depth extent**. The incorrect preliminary claim is explicitly superseded. The source genome/genetic data are not involved.

## Original 21 candidate sample age-model diagnostics

| Site / dataset | Source-labelled cal-18–24ka samples | Source actual chronology / default? | Source model-choice and direct-control diagnosis |
|---|---:|---|---|
| Sheheree Bog / 14556 | **9** | 7725 / **No** | Sample ages 18,564–23,611 BP, all depths 512–544 cm **below final returned control at 450 cm (8,786 BP)**; same exact nine samples in **DEFAULT** calibrated model 7726 have ages **8,592–9,413 BP** instead, no longer LGM |
| Rifugio Mondovi / 22662 | **1** | 13010 / **No** | 18,241 BP at sample depth 420 cm **below final original source control at 390 cm**; alternative DEFAULT calibrated model 13011 gives **14,751 BP**. Source reliable-age bound older 18,241 BP is reached at its outer limit |
| Laghi dell'Orgials / 22683 | **1** | 13034 / Yes | Central 20,368 BP, with original sample younger/older bounds **20,044–20,717 BP**, and depth 392 cm **within returned control depth extent 0–406 cm**, but source model `reliableagespan.older` is only **12,540 BP**. Alternative non-default calibrated model 13033 gives **17,830 BP**; do not swap to salvage a positive finding |
| Dar Fatma / 24197 | **3** | 14552 / Yes | Sample age 19,227–23,260 BP, depths 205–215 cm **between controls at 180cm/9,144 BP and 230cm/29,310 BP**. The chronological gap is **20,166 years**; model notes report multiple sediment hiatuses including at 180cm. Original sample-age uncertainty bounds absent. The alternative chronology 14551 uses **uncalibrated radiocarbon BP**, which cannot be numerically compared as the same age scale |
| Saham Mere / 41351 | **7** | 25462 / Yes | Source calibrated ages 18,777–21,044 BP; **six depths are within the source control interval 2354–2908cm**, one at 2912cm is 4cm below the deepest returned control. Six are within controls, including the 2908cm control that shares source age 20,761 BP; one outside. Original individual pollen sample age bounds absent. Source alternative chronology 25461 is raw radiocarbon years BP, not another calibrated direct comparison |

In the **21 original sample IDs**, v20 exact source-verified audit finds:

- **10** rely on a **non-default** chronology (9 Sheheree, 1 Rifugio); the same exact sample IDs' default CALIBRATED chronologies move all ten *out* of the 18–24ka window.
- **1** (*Laghi dell'Orgials*) is older than the source chronology metadata's reported reliable older age bound.
- **11** are deeper than their oldest/deepest *returned* source control; the other **10** are within returned control depth extents. Of the ten, three Dar Fatma lie in a very broad 20,166-year bracket.
- **20/21** lack BOTH older and younger individual source sample age uncertainty limits; presence of a chronology-level reliable age span is not a substitute for that missing uncertainty.
- The intersection **9/21** with a default calibrated source chronology, central age inside the source reported reliable age span, and original sample depth *inside returned control depths* consists of **Dar Fatma 3 + Saham Mere 6**. **These nine are NOT independently dated, species-diagnostic occupied plant populations.** This is a source-metadata triage for chronological follow-up, not a confirmatory hypothesis test.
- Metadata chronology controls sometimes contain a 14C laboratory identifier, but the source API entry is not an independently reproduced modern calibration curve, radiocarbon laboratory verification or a direct species diagnostic organismal sample.

## Reproducibility and next biological discriminator

- [Source-bound v20 Actions #39864856237](https://github.com/zuizui0223/chocho/actions/runs/39864856237) **SUCCESS**. Original v18 age/sample SHA and all 11 gzipped original v19 chronology-source body SHAs verified. Dataset/chronology IDs, duplicate control removal, original depth range, model default/alternative, reliable age bounds and original age uncertainty status all preserved; 21 sample and 11 chronology source tables, chronology control table and receipt committed.
- Result: `results/lgm_neotoma_21_Alnus_sample_chronology_control_qc_v20/` (11 model metadata, all returned distinct source controls, all21 candidate sample QC, five dataset summaries, machine source receipt).
- GitHub workflow uses only the SHA-frozen original Neotoma source data, not live query replacement, and does not promote any source label to historical occupied habitat.
- Primary next step if the fossil occupancy study continues: **published, independently cross-validated chronology calibration and pollen morphology** for the source six sites, prioritizing the 9 source-metadata-eligible Dar Fatma/Saham Mere sample IDs *without declaring them LGM occupied trees*. Preserve the 12 non-supportive/method-limited candidates and all source alternative age models. For a genuine test of host-mediated insect genetic differentiation a separate new genetic sample panel is still necessary.

**Interpretive restriction:** Fossil pollen might reflect wind-transported grains rather than trees growing at the core, even after reliable dating. WCVP modern native plant regions and modern insect sample coordinate anchors projected to 21ka do not establish ancient larval host occupancy. TTF negative real-other-plant matched identity control remains closed.
