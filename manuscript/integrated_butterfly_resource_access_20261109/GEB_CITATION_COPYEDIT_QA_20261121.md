# GEB journal-export citation copyedit — 21 November 2026

**Status: citation/reference alignment verified; journal submission remains blocked by the separately documented source-access and author gates.** This revision applies to Draft PR #179's current island–port–sea paper. It does not certify the earlier v0.2 review bundle.

## Result

| Check | Verified result |
|---|---|
| Main body / structured abstract | 4,979 / 292 words under the current exporter count; within the 5,000 / 300 checks |
| Journal-facing references | 24 main, nine SI; every entry has an actual author–date callout in its document |
| Historical bibliography | All 39 current source entries retained in the unchanged article; the complete original reference/appendix archive is byte-identical; 17 unused entries omitted only from the journal export |
| Scientific content | All 381 non-citation numerical tokens in the main body retain their order and value; scientific source/figure validator passes |
| Source integrity | SHA-256 guards cover the article, source SI, captions, original reference/appendix archive and numerical figure receipt |
| Display pieces | Four main and two SI figures remain embedded; eight SI appendices remain separately numbered |
| Independent review | Source metadata, author order, citation coverage and export code reviewed; material findings resolved |

The official [GEB author guidelines](https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html) specify APA6 production style and permit any consistent reference style at initial Free Format submission. This work completes copyediting; it is not journal approval or authorization to submit.

## Verified corrections and source decisions

1. **Citation style and coverage.** Expand the first citations of Couto, Halsch and Hatfield; alphabetize the Daru/van Kleunen citation groups; cite actually used data and geographic sources. Three-to-five-author first citations name everyone. References print all authors through seven and the first six plus final author for eight or more. Journal/volume italics, sentence-case titles, author punctuation and hanging reference indents are consistent. The full author rosters and primary metadata remain in `GEB_REFERENCE_COPYEDIT_V01.json`.
2. **EuDiS version identity.** The [version-specific Figshare record](https://api.figshare.com/v2/articles/22251028/versions/5) gives version 5 a publication date of 16 June 2023. Correct the exported year from 2025 to 2023 and retain the exact versioned DOI `10.6084/m9.figshare.22251028.v5`. The unused accompanying journal paper is omitted from the export, so no same-author 2023 a/b collision remains.
3. **Actual GBIF and WGSRPD sources.** The [specific GBIF citation endpoint](https://api.gbif.org/v1/occurrence/download/0013264-261023074918325/citation) prescribes GBIF.org and 29 October 2026 for `10.15468/dl.umqmsf`. This applies to the original 224-species source only. Cite Brummitt's 2001 second edition following [TDWG's recommended document citation](https://www.tdwg.org/standards/wgsrpd/), while retaining the actual polygon commit `52da7828aba9d461dd133c27b3bd7a4407161f54` in source provenance.
4. **Port/current provenance boundaries.** The [official 2019 NGA publication](https://msi.nga.mil/api/publications/download?key=16694622%2FSFH00000%2FPub150bk.pdf&lang=en) is the 27th edition, Publication 150; 35 denotes Notice to Mariners. Matching coordinates do not establish the study GeoJSON's missing download/conversion/edition chain. OSCAR analyses accessed the [NOAA legacy mirror](https://oceanwatch.pifsc.noaa.gov/erddap/info/yearly_336c_0b32_9cd3/index.html) of approximately five-day, one-third-degree gridded velocity estimates, distinct from the NASA v2 metadata-only access test and the undirected GSHHG resistance model. No byte-equivalence claim is added.
5. **Other bibliographic details.** Follow the article byline's `Read, Q.` for Yoon and Read; retain the 1967 Princeton imprint for MacArthur and Wilson; cite Govaerts as WCVP's registered dataset creator and record his additional editor role in metadata. Add verified journal issue locators and preserve each version-specific dataset identifier. No source dataset, model output or interpretation was changed.

## Verification and reproduction

Run the existing immutable scientific/editorial validator with the source snapshot `3dbf038dab64e8bbc9ba6c0f5b0e19ce1cb2b250` and the frozen v0.3 article. Then run:

```bash
python scripts/build_GEB_v05_blinded_main_and_SI.py --repo-root . --outdir dist/geb_v05
```

The build checks source hashes, exact replacement counts, author counts, complete first-use forms, reference/citation alignment, alphabetical citation groups, numerical preservation, word limits, anonymity and embedded figure counts. It emits `GEB_v05_CITATION_AUDIT.json` as well as the source-gated document preflight. The document workflow now triggers on the citation helper and manifest and includes these checks in its artifact.

The local final PDFs contain 21 main pages and 11 SI pages. All 32 pages were visually reviewed at 130 dpi, with independent review of main pages 1–11 and all SI pages; the primary reviewer inspected main pages 12–21. The SI Figure S5.1 lead-in was kept with its following image after review identified an orphaned lead-in. The final changed-page review and file hashes are recorded in `GEB_CITATION_COPYEDIT_CHECKS_20261121.json`. Local pagination can differ from a hosted LibreOffice environment; hosted checks do not replace this direct page review.

## Remaining Stage A limits

- Original contributor attribution for the historical first-index and 1,596-port GBIF acquisitions is incomplete. The verified 224-species DOI cannot be borrowed for either.
- The 83,573-row historical Parquet, 842,964-row original 155-species model panel and 12 associated model inputs remain unavailable. Stored-output validation is not an independent 155-species refit.
- A complete, anonymous current-version reviewer data/code archive and a stable publication deposit remain outstanding; the older v0.2 bundle does not satisfy this requirement.
- WPI conversion provenance and source/media rights still require resolution. Identifying author details, declarations, coauthor approval and cover-letter signature remain author decisions.

The prior Figure 2 head `cdefceca781e22b6a8864ad06a9e31a751c81460` completed all 14 Actions runs, including the [104-island HOSTS-masked negative-control audit](https://github.com/zuizui0223/chocho/actions/runs/44905322790). A newer branch head must be checked on its own run results. This task does not merge the Draft PR, contact coauthors or submit to a journal.
