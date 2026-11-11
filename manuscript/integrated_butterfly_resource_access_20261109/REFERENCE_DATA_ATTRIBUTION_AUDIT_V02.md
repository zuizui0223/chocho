# Bibliography, data attribution and figure QA — one island–port–butterfly paper v0.2

9 November 2026 | Scientific coauthor review, not a claim of journal submission readiness

## What was actually checked

- The saved **five-page integrated PDF** was downloaded from GitHub Actions and rendered at 125 dpi. Each page was visually inspected. No clipped label or Fig. 5 explanatory text collision was visible in this final saved version.
- The five-page PDF uses vector drawing objects throughout; Figures 2–4 originate from original SVGs converted to vector paths (not searchable text), and Figure 4 additionally has one raster density image. This is a scientific review bundle, not an independently certified publisher-compliant or accessible final-figure package.
- The original **12 references** were checked by DOI and article/dataset identity against publisher, PubMed or museum records. Their listed titles, publication years and DOIs matched; the complete subtitle of Slove & Janz (2011) was missing and has now been restored.
- Added the omitted 2021 butterfly migration review, original GSHHG paper and exact v2.3.7 source dataset; National Geospatial-Intelligence Agency's World Port Index (2019, Publication 150) is now explicitly named in Methods and References. The bibliographic list totals **16 entries**, 15 with journal/data DOIs and one official government WPI URL.

## Verified bibliographic identifiers

| Source | Verified DOI / official portal | Identity checked |
|---|---|---|
| Braga (2023) | https://doi.org/10.1016/j.cois.2023.101074 | Exotic host plants review, Current Opinion in Insect Science 58, 101074 |
| Braschler & Hill (2007) | https://doi.org/10.1111/j.1365-2656.2007.01217.x | Polygonia c-album, Journal of Animal Ecology 76:415–423 |
| Chowdhury et al. (2021) | https://doi.org/10.1111/brv.12714 | Butterfly migration review, Biological Reviews 96:1462–1483 |
| Daru et al. (2021) | https://doi.org/10.1038/s41467-021-27186-8 | Plant homogenization, Nature Communications 12:6983 |
| Govaerts et al. (2021) | https://doi.org/10.1038/s41597-021-00997-6 | WCVP Scientific Data article |
| Kawahara et al. (2023) | https://doi.org/10.1038/s41559-023-02041-9 | Butterfly phylogeny, Nature Ecology & Evolution 7:903–913 |
| Lancaster (2020) | https://doi.org/10.1038/s41559-020-1199-1 | Host breadth, Nature Ecology & Evolution 4:963–969 |
| Neu et al. (2021) | https://doi.org/10.1111/jbi.14258 | Pieris mannii, Journal of Biogeography 48:3016–3031 |
| Robinson et al. (2023) | https://doi.org/10.5519/havt50xw | NHM HOSTS original data set, five original contributors listed |
| Shirey et al. (2022) | https://doi.org/10.1038/s41597-022-01473-5 | LepTraits 1.0, Scientific Data 9:382 |
| Slove & Janz (2011) | https://doi.org/10.1371/journal.pone.0016057 | Nymphalinae paper; full title with global-scale subtitle restored |
| van Kleunen et al. (2015) | https://doi.org/10.1038/nature14910 | Global alien plants, Nature 525:100–103 |
| Wessel & Smith (1996) | https://doi.org/10.1029/96JB00104 | Original global high-resolution shoreline methods |
| Wessel & Smith (2017) | https://doi.org/10.5281/zenodo.7007502 | GSHHG exact version 2.3.7 archived dataset |
| Yoon & Read (2016) | https://doi.org/10.1007/s00442-016-3560-2 | Exotic host/ecological traps, Oecologia 181:985–996 |
| NGA (2019) | https://msi.nga.mil/Publications/WPI | World Port Index Pub. 150; original source snapshot needs reconciliation |

## Source rights and availability — what is known vs still pending

**NHM HOSTS:** Official archived dataset cites CC0. Source: https://data.nhm.ac.uk/dataset/hosts .

**GSHHG:** NOAA/NCEI describes GSHHG under LGPL; use exact 2017 version 2.3.7 plus the 1996 methodological citation. Source: https://www.ngdc.noaa.gov/mgg/shorelines/shorelines.html .

**LepTraits:** Scientific Data article is CC BY 4.0; check the exact underlying dataset and consensus table distribution licence separately. Source: https://www.nature.com/articles/s41597-022-01473-5 .

**WPI:** Official US-government 2019 Pub. 150 and data portal confirmed. Reconcile the paper's actual local port snapshot, conversion and identifier version with it before the final data availability statement. Source: https://msi.nga.mil/Publications/WPI .

**WCVP:** Source journal paper is verified. The actual source downloaded distribution snapshot, version and redistribution terms remain an independent open gate; a journal article's CC licence is not automatically a licence for underlying data.

**GBIF new butterfly/port index:** This is the highest-priority open gate. GBIF's research citation guidelines (https://www.gbif.org/citation-guidelines) require appropriate occurrence-data DOI citation, including data reached through the search API; for derived queries they recommend a derived-dataset DOI or traceable original dataset-publisher citations. The existing exact-record download DOI **10.15468/dl.pp5nc9** is ONLY for the earlier 53,434-record secondary butterfly validation. It does NOT identify the newly re-queried 225-butterfly×51-landmass indexed source, nor the 1,596-port plant queries. Query receipts and SHA-256 prove what was requested, but do not replace data-origin attribution. Complete the underlying datasetKey attribution or obtain the appropriate derived-data DOI before submission.

**Restricted historical source:** The older 83,573-row species–island–year Parquet is not an openly accessible raw-source file. Do not claim full unauthenticated replication of every historically dependent analysis until approved licence and access routes are established.

## Scientific inference / submission review gates remaining

1. **GBIF data-origin citation and exact WPI/WCVP source/version attribution.** No fake DOI or claim of permissions. These must be verified before a complete submission.
2. **Cross-grain interpretation review.** Global 239 butterflies; botanical 98 islands; 155 butterfly predictions; and a study-selected 225×51 island source are separate units. Correlation of source-overlapping HOSTS/WCVP variables is not an independently demonstrated mechanism. The observed 208/1,163 compatibility rate exceeds only its family/degree null. The historically prior 34/113 cohort is a DIFFERENT test, and its geographic plant-prevalence null eliminates excess.
3. **Butterfly geography/ports check.** A geographic conditional ocean-resistance association coexists with poor geographic-out-of-region predictive transfer and nonrobust trait moderation. The 111 literature-linked migration proxy candidates are positive-only; 114 missing species are UNKNOWN rather than nonmigrants.
4. **Journal-size figures, citations, anonymization and coauthor signoff.** All five figures and a five-page review PDF exist. Double-anonymous Word generation and its CI have passed, but individual figure typography at publisher column size, SVG path accessibility, final caption formatting, all dataset-use credits, coauthor/funding details, and external submission approval remain incomplete.

## Finished reproducible artifacts

- Island–port butterfly canonical article: ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V02.md.
- Five-source-verified figures and original renders: figures/integrated_island_port_butterfly_v02/ and FIGURE_CAPTIONS_V02.md.
- Source-typed result population and 27 audited claims: EVIDENCE_LEDGER_V02.json.
- Review DOCX and scientific figure build workflows: successful CI is a computational quality check, not journal approval.

**Next single gate:** Obtain/verify actual GBIF data citations for the independent port plant index and 225×51 butterfly source, and recover precise WPI/WCVP data versions before submitting the already-unified paper. Do not spin up a new separate GEB/JBI ecological project.
