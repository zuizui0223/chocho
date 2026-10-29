# UK butterfly historical prediction: public visit-count source gate

2026-10-28. Exploratory independent continuation, stacked on PR #96, original GEB main unchanged.

## Pivot justified by original data checks

PR #96 established 37/39 historically recordable butterflies but UKBMS GBIF DOI 10.15468/gmqvmk contains PRESENCE-only occurrence records and provides no explicit visit/absence denominator. Do not treat unrecorded species at reported site dates as ABSENT without a verified common all-species survey.

**New independently discovered source**: Natural England Long Term Monitoring Network (LTMN006), https://publications.naturalengland.org.uk/publication/5692468885979136, 2023 downloadable records:
- Original UKBMS **visit** file: https://publications.naturalengland.org.uk/file/6020843017338880 (XLS, nominal 2.6MB).
- Original UKBMS **weekly count** file: https://publications.naturalengland.org.uk/file/6589256320155648 (XLS, nominal 20.0MB).
- Original UKBMS **transect crosswalk**: https://publications.naturalengland.org.uk/file/6033934228652032 (XLSX, nominal 30.6kB).
- Historical UK Environmental Change Network butterfly observations 1993–2015, DOI 10.5285/5aeda581-b4f2-4e51-b1a6-890b6b3403a3, URL https://data-package.ceh.ac.uk/data/5aeda581-b4f2-4e51-b1a6-890b6b3403a3 (published comma-separated values).
- UKBMS published site index data 1976–2024, DOI 10.5285/53a1c705-4437-48bb-9bcd-b68977f8edd9, URL https://data-package.ceh.ac.uk/data/53a1c705-4437-48bb-9bcd-b68977f8edd9 (published comma-separated values, site index is **model-derived**, not raw zero).

Natural England notes 32/37 LTMN sites have associated UKBMS transects. Natural England 2023 original UKBMS LTMN counts and visit time may constitute real survey-level denominator, but only if **original visit site/week keys and species count site/week keys are demonstrated to match**. A literal zero in a species count may only be a true species-level nondetection if the source table includes complete checklists and its semantics prove zero was recorded at a visit. Absence of a row is NEVER converted to zero by default. If an original source's data licence limits sharing, retain only source checksums, schema, count summaries and bibliographic attribution, not original rows.

## Source-first stage (prior to inspecting data rows)
Attempt bounded public download from these exact publishers, no user credentials, no alternate scraped mirrors. Verify content format, bytes, SHA256, ZIP integrity, CSV header and row total for public EIDC data; XLS source binary OLE header and dataset declared size only at this first source gate. Any 403/HTML/login or format mismatch is *not* evidence of zero survey data. No original transect identifiers/locations or rows in Git or artifacts. Source retrieval receipts explicitly state blocked/available.

## If source permits (next scientific action)
Predefine historical training 1993–2003 and 2007–2015 test for ECN, or 1976–1995 and 2005–2022 for LTMN if its actual sample spans that range. Require ≥3 fully monitored transect sites with both time windows, ≥3 target butterfly species and >=100 verified visits per period, plus survey-wide rather than focal-targeted protocols. Compare heldout **actual species count** or detections, not GBIF index record counts. Minimum unbiased M0 climate/effort/site+habitat and M1 pre-epoch wild-host distribution; human-cultivated host only if a separate *independently dated* plant resource record is accessible (current host catalogue alone is inadequate).

All subsequent statistical choices must be frozen before opening out-of-time butterfly outcomes. If no valid visit+count join, retain feasibility/blocked status and stop without artificial hindcast. First-stage discovery is not proof of host mediation.
