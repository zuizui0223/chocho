# v11 | Dated natural history changes the apparent butterfly colonization timeline

**2026-10-13 JST · Post-hoc source audit, not a prediction of colonization.** Frozen GEB PR #38 remains at `88209b8a215ac1bc9c4c5e67622fa104179c8895`.

## What is supported

Among the **26** butterfly×WGSRPD3 units first retained after 2017 (215 occurrence records) that all four host-source reconstructions label `introduced_only`, historical literature already documents the butterfly before 2018 in **8 units (182 records)**. The previous three cases (Dryas THA, Dryas MLY, Historis BZS) now have five additions:

| Focal pair | First retained GBIF | Conservative published-before bound | Historical source |
|---|---:|---:|---|
| *Troides aeacus* / West Himalaya WHM | 2018 | 2011 | Smetacek 2011 DOI 10.18473/lepi.v65i2.a5, wild larval host *Aristolochia dilatata* and rearing |
| *Melanitis leda* / Marianas MRN | 2020 | 1942 | Swezey 1942, Guam; historical report of rearing caterpillar on corn in 1911 |
| *Papilio polytes* / Marianas MRN | 2024 | 1989 | Nafus & Schreiner 1989, Guam; abundant, larval hosts/citrus damage/parasitoid introductions |
| *Pieris brassicae* / Chile Central CLC | 2026 | 2014 | Zúñiga-Reinoso & Mardones 2014; Valparaiso in known early Chilean range |
| *Pieris brassicae* / Chile South CLS | 2026 | 2014 | Same paper; Valdivia in known historical distribution |

**One** extra species×region is a taxonomic ID HOLD (*Pyrgus communis*, Tennessee, five GBIF records). **17 remaining units with 28 records have not received complete historical auditing**. Because the eight confirmations were selected *post hoc* from a biased snapshot, this does not estimate a false-first detection rate and does not prove that no biological colonization occurred in the remaining regions.

## Strong additional check: WHM 'introduced-only' is wrong when a documented native host link is restored

The exact pinned WCVP source commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05` was queried using [GitHub Actions run 38399833118](https://github.com/zuizui0223/chocho/actions/runs/38399833118). Its accepted species *Aristolochia dilatata* has plant ID `2651355` with **native, extant, nondoubtful West Himalaya (WHM)** distribution. The frozen HOSTS list for *Troides aeacus* does not include that larval host although Smetacek (2011) reported it and attempted rearing.

A one-link, targeted **sensitivity** recodes *Troides aeacus*×WHM as `native_host_present`, not `introduced_only`, in all four source versions. In marker+BCE, the 239-species native/added resource cell totals change **29,855/15,072 → 29,856/15,071**; contemporary stays **44,927**. Only one cell is reclassified; the qualitative global expansion result is unchanged. Do **not** replace journal counts without rebuilding all final outputs.

Smetacek's observation also points to a biological hypothesis: *T. aeacus* presence may depend on **maturity, seedpod production and tissue quality** of its native vine host, not botanical geographic presence alone. This is an observational/local rearing insight, not a causal globally estimated effect.

## Source links and evidence boundaries

- Smetacek (2011): https://doi.org/10.18473/lepi.v65i2.a5
- Swezey (1942): https://hbs.bishopmuseum.org/pubs-online/pdf/b172p30-38.pdf
- Nafus & Schreiner (1989): https://micronesica.org/volumes/22
- Zúñiga-Reinoso & Mardones (2014): https://doi.org/10.4067/S0718-34292014000300015
- Burg et al. (2014): https://doi.org/10.1371/journal.pone.0104076
- Santos et al. (2011): https://doi.org/10.1590/S1676-06032011000300022

The 1911 report and 1942/1989 publications are **not butterfly establishment dates**. Historical literature proves prior regional presence, not actual local colonization time or generational persistence. The frozen GBIF export has no retained life-stage variable.

The complete 28-file self-contained v11 audit ZIP (SHA256 `95def5e9c8b6eaa943efa29dcc4c3de7a019ba6901b60d54a5adb1594995e958`) includes the v10 26-case input, four source portfolios, official WGSRPD3 geometry, pinned WCVP source rows, original HOSTS focal rows, eight dated literature anchors, deterministic output code, five passing tests, and a *separate* independent verifier. All checks passed again **after ZIP extraction**. Source notes and result receipt are saved here, but local v11 processing code is not remotely merged. The frozen paper PR remains unchanged.

**Research recommendation:** Stop treating GBIF first dates as colonization events. Restrict the next ecological temporal hypothesis test to independent, dated plant introduction plus host-specific wild egg/larva/offspring records; test plant developmental stage and usable tissue, not just presence.
