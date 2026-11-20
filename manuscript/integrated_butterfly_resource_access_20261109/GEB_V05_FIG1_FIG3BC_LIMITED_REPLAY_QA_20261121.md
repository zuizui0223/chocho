# GEB v0.5 — Fig. 1 and Fig. 3B/C limited replay: source-tested, anonymous-identifier screened, NOT full reviewer archive

**2026-11-21, coauthor QA on Draft PR #189 only.** Paper Draft PR #179 at fixed head `99b876421af2b4f51e58569862994622a02717d7` was not changed. Do not confuse the source-screened limited ZIP with a licensed, durable, anonymous deposit.

## What actually ran

- Source replay script: [`scripts/replay_GEB_v05_fig1_fig3BC_aggregates_20261121.py`](../../../scripts/replay_GEB_v05_fig1_fig3BC_aggregates_20261121.py). **Standard library only**, quoting-aware Python CSV, exact Git blob guards on all three input tables.
- Pinned source files: `results/island-mainland-resource-q1/region_metrics.csv` (**369 WGSRPD3 botanical regions**), `results/island_resource_gap_filling/by_species.csv` (**239 butterflies × 2 geographic classes**), and `data/frozen/figure_sources/anthropogenic_species_metrics.csv` (**239 butterfly trait/resource rows**).
- The source-pinned main figure numeric receipt and original plant-contribution report were checked from source Git blob SHA256-equivalent identities, but those original reports **are not** included unredacted in the author-identifier-screened ZIP.
- [**Successful Actions #44837370642**](https://github.com/zuizui0223/chocho/actions/runs/44837370642), [limited archived artifact #14679092683](https://github.com/zuizui0223/chocho/actions/runs/44837370642/artifacts/14679092683). **Seven-file limited ZIP** `GEB_V05_LIMITED_FIG1_FIG3BC_REPLAY_FIG3A_RECEIPT_ONLY_NOT_FULL_REVIEW.zip`: SHA256 `5f35415b81d26f87538381c53376f4ab8e4069a85069e7d510ddff70fdd9fcec`. Exact three data tables, executable replay script, CLI wrapper, scope/license README and stripped Fig.3A source-only numerical receipt.
- Seven-file ZIP CRC, member names/modes/time stamps, and direct account/name/email/home-directory/private-source path scans passed. After unzipping in a separate directory, **all source-derived Fig.1 and Fig.3 B/C numbers were re-executed and passed**; the same independent extraction was also checked in the assistant's working runtime. These are direct-string scans only, not proof of indirect anonymity or redistribution rights.

## Verified from the three included source CSVs

| Figure / check | Independently replayed numbers | Scope |
|---|---|---|
| Fig.1 — island resource opportunities | **126 botanical island regions**, 239 candidate butterflies; native **4,829**, current **10,570**, introduced-associated additions **5,741**, relative gain **+118.8859%** | Potential known-host geography, **not** observed butterfly food intake |
| Fig.1 — mainland opportunity comparison | **237** botanical mainland regions; native **21,215**, current **29,901**, additions **8,686**, relative gain **+40.9427%** | Same fixed candidate butterfly panel; not a causal plant-introduction intervention |
| Fig.1 — reasons for island amplification | Initial native coverage **16.04% islands** vs **37.45% mainland**; originally missing combinations newly supplied **22.705% islands** vs **24.517% mainland** | Low initial native baseline explains relative growth, not higher island gap-filling efficiency |
| Fig.1 — full 369-region source total | Original native **26,530** → current **41,083**, additions **14,553** | 126 islands + 237 mainlands + six mixed regions; excludes neither source stratum |
| Fig.3B — consumer geographic breadth | **239 butterflies**, source-fixed Spearman **ρ = +0.5788964069** between original global native-host-region count and added island opportunities; **5,741** island additions | Source-mapped host opportunity, not butterfly dispersal or establishment |
| Fig.3B — source rank quartiles | Broadest **59: 2,271/5,741 (39.56%)**, narrowest **60: 644/5,741 (11.22%)** | Two butterflies tie at source native geographical breadth **169** at the 59-species threshold; reverse lexical tie choice changes the broadest total to **2,276** (only 5/5,741 additional units) |
| Fig.3C — one host-family group | **97** butterflies with one documented host plant family; **82** gained potential global host-resource units. All butterflies gaining globally **206/239**, island-gaining **198/239** | Does not show realized advantage or survival of specialist butterflies |

**Fig.3A is not independently reconstructed.** Its fixed source reports agree: **670** contributing plant species, **38** contributing half of 14,553 added resource units, Gini **0.7773943**, top ten share **0.2510105**. The original butterfly × region × contributing plant fractional-credit edges are **not included** in the three exposed aggregate CSVs, so neither the original 670-plant ranking nor Gini can be recomputed ab initio. The ZIP includes a stripped, explicitly **`NOT_REPLAYED`** summary receipt and never labels this panel as independently reproduced.

## Where paper reproducibility now stands

1. **Fig.1**: already-derived aggregate counts/percentages independently repeat, NOT upstream WCVP/HOSTS data ingestion.
2. **Fig.2D**: earlier eight-file limited source ZIP independently reruns 126-island fixed239 WPI/host-gradient and Moran/MEM10 spatial sensitivity. [Successful CI #44769070926](https://github.com/zuizui0223/chocho/actions/runs/44769070926), SHA256 `c4a9f05f593f59924bdd188c88b7eb036b5d5a89804d665369ab0f4621ba9783`.
3. **Fig.3B/C**: independently repeated from already-derived butterfly tables. **Fig.3A only frozen-source numerical receipt**, not original plant-credit computation.
4. **Fig.2A–C, Fig.4 and other SI**: numerical audit/source receipts and limited legacy CIs exist, but a complete separate-extraction, approved-rights, anonymous-scope original data+code replay has **not been delivered**. The 155/224 observed butterfly cohort cites the verified **shared original GBIF DOI `10.15468/dl.umqmsf`**, but its raw source-to-model transformations are not independently bundled here.

**This does not certify full-paper reproduction, publication rights, a journal-ready anonymous archive, stable reviewer-access URL, complete original GBIF publisher attribution for the historic 225×51 or 1,596-port facet sources, or coauthor approval.** Those are open independent release gates. The older v0.2 anonymous CI is not itself evidence that the v0.5 GEB source tables can be redistributed anonymously.
