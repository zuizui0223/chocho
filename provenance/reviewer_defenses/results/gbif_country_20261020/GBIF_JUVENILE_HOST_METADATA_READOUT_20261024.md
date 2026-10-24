# Public GBIF stage-labelled butterfly records contain no explicit verified host-use field in this source cohort

**23 October 2026 data-audit execution; archived 24 October 2026.**
All **15** exact fixed life-stage × butterfly taxon × country filters were
rechecked against the original **39** GBIF indexed records from
the first 48 source-positive butterfly-country cells. A source-locked
GitHub Actions runner inspected selected public record fields transiently
and retained only aggregate counts. No observer, exact locality,
record ID, photograph, exact date or free-text comment is persisted.

## Verified metadata availability

| Field coverage metric | Indexed GBIF records |
| --- | ---: |
| Frozen stage-labelled source records | **39** |
| Any one of six selected annotation/association fields nonempty | **12/39** |
| Any selected field containing source parser host/feeding action keyword | **0/39** |
| Exact matching botanical candidate plant binomial mentioned | **0/39** |
| Both accepted candidate plant binomial and host/feeding keyword | **0/39** |

Original fields inspected: `associatedTaxa`,
`associatedOccurrences`, `occurrenceRemarks`,
`organismRemarks`, `habitat` and `dynamicProperties`.
No exact candidate name appears among those selected fields in this
frozen original 39-source cohort.

**This is missing observational metadata, not evidence that larvae
had no hosts or did not recruit.** Furthermore, only five
source-added-host-country cases had the complete strict accepted
plant binomial candidate list; unlisted plants and historical
synonyms were not invented or automatically attributed.

Reconstruct the analysis with
[audit script](../../../scripts/audit_gbif_juvenile_host_annotations_20261023.py)
and the [completed source-locked workflow](https://github.com/zuizui0223/chocho/actions/runs/39919887018).
See the [numeric, privacy-safe country-cell aggregate summary](gbif_juvenile_host_annotation_coverage_20261023.json)
for exact source SHA-256 and all 11 source country-cell totals.
