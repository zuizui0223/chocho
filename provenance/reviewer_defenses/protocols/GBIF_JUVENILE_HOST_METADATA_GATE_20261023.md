# Source-locked GBIF juvenile-associated host metadata gate, 23 October 2026

**Posthoc feasibility audit, not a confirmatory test of butterfly recruitment.**
Freeze the query frame and output fields before accessing any new record bodies.
This is a source-audited effort to find whether **existing public indexed
early-stage butterfly records contain host-plant evidence**. It is not a
new global ecological sample or a search for a statistically favorable effect.

## Fixed sample and source identity

Read the **same 15 exact positive species×country×lifeStage categories** from
the archived, source-complete
`gbif_positive48_gadm_lifestage_facets_20261020.json` and
`gbif_juvenile_gadm1_conditional_20261021.json`.
Source contains exactly **39 GBIF indexed Egg/Larva/Caterpillar/Pupa labels**
across **11 species×country cells and three butterfly species**, each previously
matched to its original country GADM1 subset. Do NOT select only cases
whose source host opportunity or earlier spatial correspondence was positive.

Query GBIF `/v1/occurrence/search` with the exact fixed taxonKey, country,
lifeStage original raw spelling, year=2010,2025, hasCoordinate=true,
hasGeospatialIssue=false, occurrenceStatus=PRESENT and the original
five accepted basisOfRecord types. Use `limit=100`; all 15 source groups
have at most 17 records. Reject count drift, missing results, incomplete
pagination and inconsistent stage/taxon/country identity; never record such
a response as no ecological evidence.

## Privacy-aware outcome extraction

The runner may transiently inspect public GBIF record bodies, but **must not
save or upload raw records, observer identities, observer IDs, locations,
exact dates, references, original occurrence IDs, photos or free-text
remarks**. Persist only query URLs and response hashes; per-group counts
and fields-coverage summaries. Source text used for keyword checks is
ephemeral. Retain exact source candidate **plant binomials** from the
five previously fixed added-host juvenile country cases solely to check
literal matches; do not invent synonyms or fuzzy accepted links.

For all 39 records report counts with nonempty:
`associatedTaxa`, `associatedOccurrences`, `occurrenceRemarks`,
`organismRemarks`, `habitat`, and `dynamicProperties`. Independently
report records containing an **exact candidate accepted binomial** in
one of those fields and records with both that name and text explicitly
indicating feeding/oviposition/host use. String matches are
**unverified host mentions**, not larval feeding validation.
Report missingness separately from negative biological host evidence.
Distinct GBIF record keys are checked for duplication **but never saved**.

If provider text contains non-public personal information, reject it
from saved outputs. A source's `associatedTaxa` field could denote
predators, parasites or habitat associates rather than the larval host;
do not reclassify it into host-use evidence without manual record
adjudication. No dependence on output p values.

## Ecological boundary

This data describes only the accessibility of field annotations for the
preselected GBIF sample; it neither verifies wild breeding nor a new host
introduction, survival or population establishment. Only records with
verified plant association and developed offspring could address
potential host opportunity becoming realized butterfly recruitment.

Keep new files out of frozen GEB manuscript and anonymous reviewer pack.
