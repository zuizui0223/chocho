# Posthoc dataset registry metadata lookup for seven country-cell GBIF source facets

The preceding 7/7 source-facet audit was completed before this metadata lookup.
It found six distinct GBIF dataset keys. This follow-up is explicitly **posthoc**
and limited to reading their registry metadata, not further selecting or
excluding butterfly-country outcome cells. No new source-dependent butterfly
effect sizes are estimated.

For each of the exactly six distinct frozen GBIF dataset keys, read public
`GET https://api.gbif.org/v1/dataset/{key}` once, with bounded retries,
preserving a source URL, UTC timestamp, HTTP status, SHA-256 of the response,
title, publishingOrganizationKey, type and DOI (when provided). The original
dataset keys must appear exactly once. Do not infer independence from a
different dataset title, DOI or publisher; the same field observation may be
republished. Do not fabricate a replacement when a dataset is missing.

Retain source errors separately as unavailable, not ecological zeros.
Use the resulting metadata **only** for distinguishing metadata provenance,
especially the two 42/43-record publishing sources of `Papilio demoleus`
in Israel. Continue to restrict the primary ecological inference to
potential resource geography rather than actual larval recruitment.
