# v10 — A GBIF 'first detection' is not a butterfly colonization date

**2026-10-13 JST; retrospective source audit, exploratory, no new field observation.** The frozen GEB PR #38 (`88209b8a215ac1bc9c4c5e67622fa104179c8895`) remains unchanged.

## Exact frozen cohort
53,434 GBIF records for 31 mapped of 32 nominated butterflies. Recovered 48,728 mapped occurrence rows and 1,081 butterfly×WGSRPD3 cells with independent GeoPandas and Shapely geometry joins. Among 1,081 cells, 674 have early 2010–2017 records while 407 are first retained in the 2018–2026 snapshot. **28 of 31 butterflies have exactly 1,800 records**, suggesting a strong cap in the archived retrieval. **129 of 407 late-first cells come from four species with zero pre-2018 archived records.**

The first year must retain every publishing provider, not select an arbitrary first GBIF key. Exactly **278** late-first cells have iNaturalist alone in their first year, **29** iNaturalist and another provider, and **100** no iNaturalist. The earlier 293 iNaturalist first record estimate, based on one arbitrary row within ties, is superseded. Official dataset: https://www.gbif.org/dataset/50c9509d-22c7-4a22-a47d-8c48425ef4a7.

Among the 407, the source-enriched marker-corrected fixed host map has **358 native**, **26 introduced-only**, **23 no known host**. Seventeen introduced-only cells have observations in only one year.

## Four stringent cases and historical falsification

Post-hoc diagnostic screen = introduced-only, ≥2 recorded years, ≥5 retained GBIF records, ≥10 pre-2018 OTHER focal butterfly observations in the same region.

1. *Dryas iulia* **THA** — first archived 2022, but wild historical records since 2007 (Burg et al. 2014; DOI 10.1371/journal.pone.0104076).
2. *D. iulia* **MLY** — first archived 2019, but historically documented by 2014 (same Burg source).
3. *Historis acheronta* **BZS** (Brazil South) — first archived 2020, but *H. acheronta acheronta* is listed in Santos et al. 2011 *Biota Neotropica* 11(3), table 1, for Subtropical Atlantic Forest in Rio Grande do Sul (DOI 10.1590/S1676-06032011000300022). Survey includes 2001–2002 and 2006–2008, exact year of first field occurrence not allocated to this taxon. Maquiné is verified inside pinned WGSRPD3 BZS geometry; 2011 is a conservative publication bound.
4. *Pyrgus communis* **TEN** — first archived 2021 but **TAXON IDENTIFICATION HOLD**, because *P. communis* and *P. albescens* (now *Burnsius*) overlap and may not be reliably identified visually, and saved dataset lacks identifying voucher/genitalia. Sources [Tennessee Valley NABA](https://www.tnvalleynaba.org/butterflies-of-tennessee/skipper-hesperiidae/common-checkered-skipper-pyrgus-communis/) and [BugGuide](https://www.bugguide.net/node/view/354354). Do not mark species-specific pre-2018 presence without specimen verification.

**Outcome:** three pre-2018 historical regions, one taxon hold, **zero verified NEW colonizations among this screened set**. This is not a validation failure rate for all 407 or 26 or an assertion that butterflies cannot colonize introduced host regions. The remaining 22 weakly evidenced introduced-only late-first cells were not checked with the same depth.

## Negative-control observation thinning

582 cells have BOTH early and late retained records. Holding biological distributions fixed, randomly thinning early GBIF records to 50% (999 Monte Carlo replicates) makes a median **85** formerly-known early cells falsely appear first after 2017 (simulation 95% range **70–99**). The count of late-first cells would increase from 407 to median **492** with no colonization. This illustrates recording-process sensitivity, not statistical prediction of the missing biological population.

## Research consequence and limits

**Do not interpret GBIF first detection after 2017 as an arrival, establishment, dispersal, breeding or a host-before-consumer lag.** The contemporary WCVP host distributions do not encode years of arrival, and some host relations were not known to the model but were already in field literature. The 2018 partition is a retrospective screen with post-hoc selection based partly on later observations; it is not an independent time-validated forecast. Record exposure and iNaturalist provider effects require explicit modeling; genuinely documented larval host use, local species identity, reproduction and pre/post dated plant introductions are necessary for an ecological temporal test.

Artifacts: Complete source-pinned v10 ZIP provides scripts, 7 passing regression tests, independent Shapely reconstruction, all 407 and 26 cells, case ledger, input checksums and thinned-event simulation. GitHub machine-readable summary at `provenance/reviewer_defenses/results/butterfly_temporal_source_window_v10_20261013.json`. No main manuscript or PR changed.