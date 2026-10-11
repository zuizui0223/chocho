# Completion evidence and a known-host mapping gap — v02

Date: 2026-10-11. Targeted retrospective audit; not a systematic review or a new field discovery. The frozen submission at 88209b8a215ac1bc9c4c5e67622fa104179c8895 remains untouched.

## Evidence gained

The official Fernandez-Hernandez & Minno (2015) PDF, printed p.10, was obtained and visually checked. Papilio demoleus laid eggs on a small citrus tree at Sierra del Chorrillo, Cuba, on 20 November 2013. One field-origin individual subsequently emerged as a female after rearing on young Citrus × limon leaves. Approximately 28 days is from hatching to adulthood, not oviposition to adulthood. The nine collected swallowtail eggs were not all confirmed to belong to the same species; 1/9 is therefore not an identifiable P. demoleus survival estimate. The collection tree was identified only to Citrus, and rearing food must not be assigned to it retrospectively.

The targeted ledger now contains three species-region cells supporting local breeding initiation, compared with two before full-text retrieval: Danaus plexippus/AZO, Pieris brassicae/NZS and Papilio demoleus/CUB. Strict frozen-host-specific field-larval correspondence remains one cell. One field-origin reared-adult event is supported, but no individual-linked field cohort survival rate or next-generation reproduction is estimated. The nearby Triphasia nondetection is not evidence of unusability; Las Tunas immatures identified as Heraclides andraemon are explicitly excluded.

## Mapping gap, not a new host switch

Pinned raw HOSTS contains four P. demoleus–Citrus limon rows, HOSTS IDs 92900, 92901, 92946 and 92962. In pinned WCVP, the accepted species is Citrus × limon, plant_name_id 2724178, and CUB is introduced, extant and nondoubtful. Literal taxon-name matching drops this known association because of the hybrid symbol. The species-rank filter itself passes. The existing six CUB host taxa would become seven under this single curated alias, without creating a new CUB resource cell. No global recalculation or change to +54.9% has been performed; broader name reconciliation may affect native as well as introduced envelopes.

## Seasonal interpretation

Phillips et al. (2020), Methods / Search timing and frequency, reports pupal overwintering, partial aestivation and 2–4 generations per year in Nelson. Delayed emergence is not synonymous with mortality. Its S1 Data has six sheets and 40 monthly detection rows, but no individual cohort IDs, host-species-specific fates or separate stage transitions. The verified annual totals are 263,999 inspections, 111,159 host-plant detections and 2,884 inspections detecting P. brassicae. Those are not individual survival denominators. The November 2000 Azores table-only nonemergence candidate remains provisional: original table image, denominator and complete follow-up were not verified.

## Source independence and verification

The acquired pinned HOSTS export has 140,485 rows and no bibliographic-reference or observation-date fields. Source independence cannot be certified from this export; this is not proof that the papers overlap it. No independent predictive-validation rate is claimed.

The full updated downloadable code package passed 112 tests locally, including exclusion of transported/wrong-taxon records, mixed-cohort denominator safeguards, and retention of genuine zero-survival cohorts when mortality and full follow-up are explicitly documented. An independent sparse-matrix implementation again reproduced 26,530 native, 41,083 contemporary and 14,553 added cells. The incremental patch was applied to the exact preceding v01 code and every changed file matched the v02 bytes. The local classifier and ledger remain in the downloadable patch/package; these receipt and note commits do not imply they have been merged into the GitHub code branch.

Detailed receipt: ../results/butterfly_completion_source_audit_v02_20261011.json (records the earlier 111-test checkpoint; the subsequent zero-survival regression adds one passing test without changing the real-data conclusions).

Primary sources: https://journals.flvc.org/troplep/article/view/89747 ; https://journals.flvc.org/troplep/article/download/89747/86082/116584 ; https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0236791 . Pinned source-acquisition runs: 38126096883 and 38126320915.
