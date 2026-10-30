# Early-plant facet pagination repair

The first fixed28 acquisition finished with25 complete islands and3 explicit missing queries. Three source requests returned duplicated species keys across100-bin pages despite stable total counts; the duplicated keys sit at tied record-count boundaries. No partial result was classified complete or used as absence. The entire initial result, source archive and code are preserved.

This is a source-transport repair, fixed before host-record joins. It changes no target selection, host roster, early window, taxonomy namespace, biological predicate, polygon, butterfly outcome or scientific threshold. Only the3 failed original target queries are eligible; no complete island is reacquired or replaced. Repair results and original failures remain separate, with an integrated manifest explicitly identifying each row's source.

## Fixed repair contract

1. Request the same exact source predicate with all four facets, offset0 and facetLimit10000 (the same maximum10000 bins implied by the original100 pages of100). Require every field to contain fewer than10000 positive unique bins. Source response must remain count-only.
2. For each of SPECIES_KEY, BASIS_OF_RECORD, DATASET_KEY and YEAR, issue a facet-specific terminal probe at offset equal to that field's actual returned bin count, facetLimit1. Require the same total count and no positive terminal bin. This detects a silent server-side page-size cap; an oversized response cannot be misclassified complete.
3. Repeat the complete offset0 request. Require identical total and all four key-to-count maps; ordering among tied-count entries is irrelevant, because all keys are returned in one page. Require unchanged total through every terminal probe. All species bins must be unique and positive; YEAR and DATASET_KEY totals must reconcile, and every year must be2010–2017.
4. Keep the original1.5-second start interval,45-second timeout,5MB response cap and3-attempt retry bound. On any error, retain an explicit missing status. Do not adaptively enlarge bounds or repeatedly rerun failed repairs after inspecting host associations.
5. Preserve actual repair request/response/receipt triples, deterministic archive, first-run lock, source hashes, script hash and this plan hash. Create a separate merged result manifest with25 original complete rows and up to3 repaired rows. Preserve original25 hash-identical source results. No failed or partial query becomes an ecological zero.

The official occurrence OpenAPI (`https://techdocs.gbif.org/openapi/occurrence.json`, SHA256 `ffde21586ee99495708d1fcb37cb88f792c50605fbb7a263217b3936f07fce39`) documents facetLimit/Offset but no facetSort or maximum facet limit. Therefore acceptance of10000 and completeness must be verified by actual responses and terminal probes, not assumed from the schema.

The repair is an adaptation to observed pagination failure, not an independent ecological replication. Data were already queried in the first acquisition; no stronger claim about pre-registration is made.
