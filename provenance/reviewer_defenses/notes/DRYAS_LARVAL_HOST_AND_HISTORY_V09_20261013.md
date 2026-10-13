# v09: Wild Dryas larval use reveals a model-correct place but an incomplete host identity

**2026-10-13 JST.** Targeted post-hoc study of the fixed `chocho` resource reconstruction. **Not an independently validated global colonization process**, and does not change the frozen GEB submission PR #38.

## Ecological finding 1: correct introduced-resource region, wrong documented supporting plant

Burg et al. (2014, [PLOS ONE](https://doi.org/10.1371/journal.pone.0104076)) document that non-native *Dryas iulia* established in Thailand, with published wild larval feeding observations on introduced *Passiflora foetida* (reported by Les Day via personal communication). The fixed HOSTS/WCVP resource reconstruction already classifies Thailand (THA) as an introduced-only host region but supports that class **solely by *Passiflora edulis***; its *Dryas* HOSTS portfolio lacks the *P. foetida* relationship. Adding this documented host increases the THA potential host list from 1 to 2, without changing the THA region's classification.

In Malaya (MLY; TDWG level3 includes Singapore), [ButterflyCircle local life-history documentation](https://butterflycircle.blogspot.com/2021/12/life-history-of-julia-heliconian.html) portrays *Dryas iulia* oviposition, caterpillar feeding on locally named *Passiflora suberosa*, and stages through eclosion. The same species is **already** represented in frozen *Dryas* HOSTS under ID 2560211, and its pinned WCVP region is introduced. Thus MLY also stays introduced-only. Site-level host voucher and a continuous unmanipulated field cohort were not verified. 2015 and 2023 photo captions appended to a 2021 blog post must retain their actual dates.

The frozen HOSTS list also includes *Lantana camara*. This is documented as adult nectar in Burg et al.; larval use by the focal butterfly was not confirmed in the reviewed sources. Removing this link is **not** a verified error repair, only a question-mark life-stage sensitivity.

| One-species host portfolio | *Dryas* native cells | *Dryas* contemporary cells | *Dryas* introduced-only additions | Whole-panel 239 introduced-only total |
| --- | ---: | ---: | ---: | ---: |
| Frozen | 45 | 139 | 94 | 15072 |
| Add reported *P. foetida* | 47 | 170 | 123 | 15101 |
| Provisional removal of *Lantana* | 42 | 129 | 87 | 15065 |
| Both changes | 46 | 163 | 117 | 15095 |

All changes are document/host sensitivities, **not observed butterfly range expansion**. Do not silently edit submission tables, standardize local *P. suberosa* names globally, or claim that *P. edulis* cannot be eaten.

## Ecological finding 2: archived first butterfly occurrence is not first real observation

Burg et al. (2014) explicitly report wild *Dryas iulia* sightings in Thailand and Malaysia **starting in 2007**, including personal communications. The exact pinned 132 GBIF records preserved for this focused audit have **6 THA records, first dated 2022**, and **126 MLY records, first dated 2019**. The Thai literature-start vs archived-sample-start gap is therefore **15 calendar years**.

This is **not** a 15-year establishment lag, plant-before-butterfly lag, or evidence of an unsuccessful resource forecast. It is direct evidence that an occurrence-only retrospective prediction might falsely treat a known historical regional presence as a novel post-2017 detection. Any hindcast using the current sampled GBIF records should first anchor known historical presences from independent faunistic literature and evaluate sampling effort. No butterfly-specific original historical event date can be inferred from an incomplete occurrence snapshot alone.

## Updated evidence triage

The previous 64 candidate butterfly×region cells already include four Nacaduba units reclassified to native by the Acacia source sensitivity. The remaining 60 contained 37 unaudited. This targeted audit labels **THA and MLY** as literature-linked local larval/life-history evidence, reducing unaudited candidates **37 → 35** (572 → 440 occurrence records). The frozen GBIF rows do **not** retain lifeStage and are not evidence of larvae.

This targeted review does not newly validate an individual life-cycle survival cohort or a global prediction-success fraction. The 22 Pyrgus source-screened records and one Dryas/Florida taxonomic HOLD remain unchanged.

## Reproducibility

Source-locked inputs, literature provenance, complete 64-cell candidate table, four scenario records, independently checked 2019–2026 GBIF year counts, deterministic auditor, independent NumPy/Pandas verification and test suite are packaged in `chocho_dryas_stage_validation_verified_v09_20261013.zip` (SHA256 `6bdef58eac544d7a870c1aa7b720633aaf60ed9552cfbc9503568a1d0c38f9b1`).

- Original v06 full regression suite freshly re-run: **169 passed**.
- New v09: **6 passed**, and all 6 also passed after ZIP extraction.
- Independent Boolean-matrix analysis reproduced every *Dryas* portfolio and all 239-species totals; independently aggregated THA/MLY chronology and 64 candidate rows.
- No code was integrated into the remote analysis branch; only these provenance receipts/notes. The frozen GEB PR #38 remains the same SHA `88209b8a215ac1bc9c4c5e67622fa104179c8895`.

**Research decision:** Prefer stage-verified host data and historical natural-history sources before claiming realized butterfly colonization from potential resource maps. The new mechanism suggestion (introduction pathways/propagule pressure + actually usable host phenology) is still an empirically testable hypothesis, not a causal effect estimate.
