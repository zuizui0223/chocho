# Eight residual port-host × sea-access GBIF first-index cases: original sources and historical status

**15 November 2026 | Completed targeted 8/8 source audit, exploratory and non-causal**

## Main decision

After screening **all 225 butterflies × 51 selected original landmasses (11,475 species–island cells)** for 1700–2000 GBIF indices, the corrected 2001–2010 candidate-host/sea-access model still had only **3 host-positive indexed first events in 2011–2017** and **5 in 2018–2022** in the ≥10-historical-port-plant-record panel.

The eight cases are now individually named, mapped to the **actual original GSHHG landmass** and independently source-triaged. **Five have historical island reports pre-dating the later GBIF first index.** The other **three have specific later source records but no documented local breeding/establishment or observed feeding on introduced hosts**. Thus **zero of the eight can presently be promoted to independently verified post-2010 island establishment caused by external larval resources**. This says nothing about the true count of unobserved colonizations.

### All eight cases, no post hoc case substitution

| Period | Butterfly | Fixed island geography | Evidence classification | Independent record / diagnostic |
|---|---|---|---|---|
| 2011–2017 | Carcharodus alceae | Great Britain (GSHHG 15) | Published older island occurrence | Surrey **1923** two individuals, historically considered accidental introductions; [Hantsmoths](https://www.hantsmoths.org.uk/lep.php?code=57.0011) |
| 2011–2017 | Issoria lathonia | Ireland (GSHHG 26) | Older Irish literature report | Two Irish vagrant records noted in an [Irish national checklist](https://www.habitas.org.uk/invertebrateireland/species.asp?item=9044), citing Baynes (1964); exact old dates and specimens not independently inspected |
| 2011–2017 | Lampides boeticus | Maui, Hawaii (GSHHG 214) | Published older island occurrence | **1983** Lepidopterists' Society newsletter explicitly reports observations on Maui; [original newsletter](https://images.peabody.yale.edu/lepsoc/nls/1980s/1983/1983_v25_n2.pdf) |
| 2018–2022 | Leptotes pirithous | Great Britain (GSHHG 15) | Published older island occurrence | Dorset **13 June 1938**; [Hantsmoths](https://www.hantsmoths.org.uk/lep.php?code=61.009) |
| 2018–2022 | Nymphalis polychloros | Ireland (GSHHG 26) | Unverified status / possible vagrant | **One** exact-key 2019 observation at Tramore, GBIF [4462800997](https://www.gbif.org/occurrence/4462800997). Regional official map [source](https://maps.biodiversityireland.ie/Dataset/363/Species/77534) does not supply independent breeding confirmation |
| 2018–2022 | Nymphalis antiopa | Sakhalin (GSHHG 28) | Published older island occurrence | *Vanessa antiopa* listed as species 28, with South-Saghalien localities by [Matsumura (1925)](https://eprints.lib.hokudai.ac.jp/repo/huscap/all/12580/15%283%29_p83-196.pdf) |
| 2018–2022 | Lycaena dispar | Madeira (GSHHG 391) | Checklist conflict / source verification required | **Four** 2020 human observations at exactly **one shared coordinate**, all from a single GBIF source provider, GBIF [5019945211](https://www.gbif.org/occurrence/5019945211), [5024616829](https://www.gbif.org/occurrence/5024616829), [5025406767](https://www.gbif.org/occurrence/5025406767), [5025680523](https://www.gbif.org/occurrence/5025680523); not on [specialist Madeira checklist](https://www.madeirabirds.com/node/82) or [BC Europe Madeira account](https://www.bc-europe.eu/webpage.php?name=endemic-butterflies-of-madeira) |
| 2018–2022 | Pararge aegeria | Orkney Mainland (GSHHG 456) | Single occurrence, establishment unresolved | **One** 2018 human observation near Kirkwall, GBIF [4453667332](https://www.gbif.org/occurrence/4453667332), [original iNaturalist 21648338](https://www.inaturalist.org/observations/21648338); reported coordinate uncertainty **2,094 m**. [Highlands & Islands atlas](https://butterfly-conservation.org/sites/default/files/2024-04/Highland%20Butterfly%20Atlas%202022%20Online.pdf) reports no northern-isles records within its earlier atlas coverage |

**The five earlier reports do not disprove any potential later recolonization**. They do demonstrate that a 1700–2000 *GBIF* search with no record can still miss published butterfly island occurrence, even after the original time-window correction. Likewise the three later exact species-key records do not establish resident status.

## Precise new source-level finding for the three uncertain late cases

Original-query audit [GitHub Actions #43666845606](https://github.com/zuizui0223/chocho/actions/runs/43666845606) **SUCCESS**, [full original GBIF individual-record and WKT/response-hash artifact](https://github.com/zuizui0223/chocho/actions/runs/43666845606/artifacts/14141399877). The analysis used the same GSHHG 2.3.7 island geometry and the exact species keys from the corrected full-panel frozen model. No re-defined geography or new case selection was introduced.

- **Madeira — Lycaena dispar**: all four GBIF events are HUMAN_OBSERVATION records on **9–12 March 2020**, all at **32.748955° N, 16.867933° W** (identical coordinate), all from datasetKey **6ac3f774-d9fb-4796-b3e9-92bf6c81c084**, with stated coordinate uncertainty 250 m. Four GBIF rows are **not** independent island colonizations or proof of larval development. A Madeira specialist checklist omits the taxon while reporting *Lycaena phlaeas*. Correct identification, provenance, and the original source site need further checking; do not assume either a misidentification or true established *L. dispar*.
- **Ireland — Nymphalis polychloros**: one exact speciesKey human-indexed observation **6 August 2019** in **Tramore, Ireland**, datasetKey **5275c9b8-f9be-4f89-b0e0-5e43518f4051**, reported accuracy 100 m. No larva, egg, or established population is documented by this record. Irish butterfly-status sources provide no independent breeding confirmation.
- **Orkney — Pararge aegeria**: one iNaturalist-derived human-indexed observation **18 July 2018** near **Kirkwall**, datasetKey **50c9509d-22c7-4a22-a47d-8c48425ef4a7**. The reported coordinate uncertainty is **2.094 km**. Observation confirms an indexed species occurrence under its assigned name, but not a population or multi-year breeding.

All six indexed records are source exact-species-key observations. They have **no recorded immature life stage** in this pulled subset. The photographic identification and geographical attribution of unusual records require independent checking before accepting resident status.

## Scientific interpretation for the GEB paper

The **positive, separately quantified** results are unchanged:
- Known-host redistribution increased *potential butterfly resource geography* across 239 butterflies by **+54.9%**; island regions had **+118.9%** relative increase compared with **+40.9%** on the mainland. This is not observed butterfly occupancy.
- Port number correlates with total introduced flora and known candidate host richness, but a *host-selective* surplus after whole-flora adjustment is **not** established.
- The separate 51-island original community target retained **+5.082% MSE reduction** from geometry-only Circuitscape in held-one-island retrospective prediction, with stringent geographic and historical transfer limitations.
- The 1700–2010 historically refitted host × uniform-marine-access **interaction did not show an independent positive effect** (≥10-plant-record panel coefficients **−1.249 pp per SD** early, **−0.311 pp per SD** late; very sparse positive host-indexed events). None of these butterfly GBIF indices directly measures island settlement.

### Implication for next analysis

Stop treating GBIF first-index rows as a colonization rate without independently dated island faunal history. Two viable follow-ons are (i) occurrence-index modeling without any claim of newly occupied islands, and (ii) an **independent natural-history validation panel** focused on dated local host utilization plus immature-stage / repeated breeding evidence, comparing arrival-access conditions. The current eight source-restricted records are a **source integrity check**, not a validation of port-mediated establishment or shipborne dispersal.

**Data:** [machine-readable fixed-eight literary crosswalk](../port_host_ocean_long_history_20261114/fixed8_external_history_crosswalk_20261115.json) and [machine-readable three-case record identity summary](verified_record_identity_summary.json). No older GEB manuscript figure was replaced, and no source-positive host-resource result was erased.
