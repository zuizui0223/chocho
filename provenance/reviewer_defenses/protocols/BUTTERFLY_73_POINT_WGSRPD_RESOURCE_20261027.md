# PR72 butterfly niche-entry: fixed point-to-botanical-resource-region test

**2026-10-27. Fixed before new point-to-WGSRPD polygon outcome analysis.** The original 73 species×country positive units (32 novel-dominant and 41 analogue-dominant) were source-frozen in successful GitHub workflow #40431411694. Their 6,076 underlying GBIF indexed butterfly records were independently re-acquired in complete workflow #40432323278, without any per-point botanical map or climate result.

## Aim and estimand

For all source-positive observed butterfly points, ask whether the actual observation coordinate falls inside a **WCVP-introduced host-resource supported WGSRPD Level-3** botanical distribution polygon, not merely inside a country where some fraction of area is supported.

Using the **same original 365-region CHELSA BIO1/BIO12 source** and exact native HOST-PLANT 5–95% and 10–90% temperature/precipitation rectangles, split added-source regions at the observation point into (1) within a documented introduced-host region *outside* native-host climate proxy (climate-novel resource-region co-index), (2) inside the analogue native-host climate proxy, (3) no introduced-added host resource in that region, (4) documented native-host region despite exclusion, or (5) unmappable/ambiguous/missing climate/invalid source.

This is a botanical regional **point-in-polygon recording proxy**, **not the actual BIO1/BIO12 value at the point**. Do not call observed butterfly fundamental or realized climatic niche expansion, actual larval host use, establishment, spread, or demographic effect.

## Data authority and retrieval contracts

- Reuse exact author-only archived input ZIP, SHA256 `408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323`, from the already audited GH Actions #39006554609. Check original manifest and every constituent SHA. Do not redistribute original HOSTS/WCVP/CHELSA/TDWG raw files.
- Original 239 butterflies, partial BCE supplement (1027 accepted host links in 81 source butterflies), 235 eligible host-climate natives, 365 CHELSA unit climate regions and exact `inputs/level3.geojson` source geometry must reproduce the existing source-bounded model and region mask.
- Exact 73 source cohort from PR #74 `pr72_exact_73_point_followup_roster_20261027.json`, no positive outcome selection; exact GBIF taxonKey/country/year 2010–2025, five original basisOfRecord filters and valid coordinates; per-pair cap 5,000 records / total 40,000, no sampling or inferred absences on cap/failure. Original vs live source count drift must be explicit.
- Point-map geometry `LEVEL3_COD` must be validated, multipart features consolidated and source region code identity checked. A coordinate matching zero valid source polygons stays unassigned; two or more is spatial ambiguity, not arbitrarily assigned. Never turn a point in a large botanical region into confirmed leaf/host availability at site.

## Fixed response and negative checks

1. For both country *host-resource* classes, report **equal-case** and **raw record weighted** fraction of records within (a) mapped introduced-supported analogous resource region, (b) mapped introduced-supported climate-novel resource region, (c) documented no-added host region, (d) other/ambiguous/unmapped.
2. Where all 73 original source positive queries are complete, report all 6,076 original indexed record identities validated with the re-query country totals; if drift, mark current snapshot and do not silently replace.
3. Report source documented <=1km coordinate uncertainty subset, duplication by rounded 0.01-degree coordinate cell, and two top original indexed recording case removals (Taiwan) as diagnostic, not cherry-picked ecological conclusion.
4. Compare original country-level fraction of introduced-host mapped area to fraction of indexed records mapped to an introduced-host botanical region, **without pretending** that observed GBIF point locations were uniform random draws of national land area. Provide no causal p-value.
5. Check 5–95% and 10–90% source climate analogues separately. No cherry-picked climate thresholds, taxon exclusions or juvenile labels.
6. Only claim **regional climatic resource-opportunity co-recording** if full source integrity and geometry checks pass. A new, genuinely fine-resolution same-product BIO1/BIO12 extraction for **both observed points and a comparable native baseline** would be necessary to call observed point climatic novelty; new wild-reproduction and host-use data are additionally necessary to claim colonization and interaction changes.

## Publication firewall

This work remains stacked on PR #73/#72, outside frozen GEB PR #38/#70 and main. Export aggregate outcomes and verified SHA receipts only, without exact coordinates, individual GBIF IDs, observer identifiers, or author-only plant/TDWG geometries.
