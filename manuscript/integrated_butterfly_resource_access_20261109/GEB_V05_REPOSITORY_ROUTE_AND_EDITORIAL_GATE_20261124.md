# GEB v0.5 — repository route and unresolved editorial/author gates

**Draft internal decision note, 24 November 2026. Nothing uploaded, sent or approved.** Current scientific manuscript unchanged. The verified 192-file staging package is retained; repository allocation below is a proposal, not a decision to release any file or to omit necessary supporting data.

## Verified policy boundary

- [Global Ecology and Biogeography author guidelines](https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html), Data and Code Provision, Storage and Documentation: data and code must be accessible during peer review, with a stable public archive for publication. GitHub alone is insufficient. The double-anonymous reviewer link must not identify authors. CC0 and CC BY are accepted; CC BY-NC requires EiC approval; other CC licences are not supported by the stated data provision policy. Ethical/legal exceptions are discretionary and must be justified, not presumed.
- [Dryad file requirements](https://datadryad.org/help/requirements/files) and [publication policies](https://datadryad.org/publication_policy): data published by Dryad must be compatible with CC0. Dryad supports linked Zenodo publication of software/supplementary files with separate applicable licensing options. Do not assume an NC-SA data file is accepted as supplementary information; confirm actual repository and journal conditions.
- [GBIF citation guidelines](https://www.gbif.org/citation-guidelines): historical API users are responsible for original source publisher attribution. Derived dataset registration requires actual source dataset DOI/key list, ideally counts. The separate original 155/224 DOI has verified 1,716 contributors. The original historical and port facets do not have an equivalent recovered provider list.

## Proposed route (pending rights and author approval)

1. **Primary research-data archive**: only author-controlled numerical outputs and transformed tables that are demonstrably CC0-compatible if using Dryad. Attribution remains scholarly practice; do not waive another rights holder's rights. Alternative stable archive with file-specific licence terms is a possible option subject to GEB policy.
2. **Code**: original author code and replay wrappers with an author-approved licence and third-party code notices, using a stable archive (e.g. linked Zenodo). No automatic MIT/CC0 assignment; repository is currently public with no declared repository-level licence.
3. **Public third-party inputs**: cite and retrieve exact HOSTS, WCVP, EuDiS, LepTraits, original GBIF DOI and other source versions from their upstream repositories, with checksums, transformations and explicit source terms. Do not duplicate raw third-party data if link-and-rebuild suffices. Source licence identification does not by itself establish compatibility of the derived material with CC0.
4. **Restricted/mixed derivatives**: retain internally until specific rights, applicable licence, suitable permitted repository and editorially accepted review access are established. CIRAD modeled N/P V1 is CC BY-NC-SA 4.0; raw grids and original source extracts are not in the current staging package. Source-level NC-SA does not necessarily attach to a statistical coefficient. Assess the particular material. The original GBIF download includes a mixed licence population; do not treat all aggregate tables or scores as licensed raw records, but also do not assert a blanket right to relicense them.
5. **Historical facets**: keep the original result tables/query scope/dates and incomplete provider attribution visible. Do not remove inputs merely to create an apparently runnable archive, invent an original DOI, or substitute the newer 1,596-port datasetKey query as if it were the original snapshot. An acceptable restricted access/attribution solution requires editorial and source-provider guidance.

## Author decisions required (not sent)

- Confirm authority to archive the authors' code and numerical outputs, and select explicit code and data terms separately. Confirm any third-party code or supplied data permissions and the appropriate release route.
- Confirm the repository route after checking actual file-specific rights. A private review link is not by itself a copyright permission. Confirm treatment of the public manuscript/source repository in double-anonymous review; do not change repository visibility or history automatically.
- Confirm all author identity, funding, conflict-of-interest, contribution and submission statements separately after technical preparation.

## Draft editorial clarification (for review, not sent)

Dear Editorial Office,

We are preparing a Research Article on butterfly larval-host resource geography. Before submission, we would appreciate clarification of the data and code provision requirements for a study combining multiple third-party sources with different reuse conditions.

We have verified the original GBIF download used for our butterfly-record analyses (DOI 10.15468/dl.umqmsf), including its contributing datasets and archived record-level licences. Separate historical plant/port and island facet queries are preserved as aggregate snapshots, but their original publisher-level identifiers were not retained. A later provider query is available, but cannot establish the original provider attribution. Our intended approach is to distinguish these sources explicitly rather than attach an unrelated or retrospective DOI.

Some upstream data, including a modeled nutrient product, have reuse conditions that are not necessarily compatible with a CC0 deposit. We intend to cite the exact public versions and provide acquisition and transformation code, and to archive the numerical supporting material only under the permissions applicable to it. We will not assume that every upstream licence applies to all statistical summaries, or that an archive link itself constitutes reuse permission.

Could you advise whether (i) a split, appropriately licensed stable archive with an anonymized review route, (ii) exact external-source links and reconstruction code for material not suitable for redistribution, and (iii) explicit disclosure of the original historical facet attribution limitation would satisfy the journal's requirements, or whether a specific exception or additional documentation is needed? We recognize that the supporting data and code must be accessible for assessment and that any exception requires editorial approval.

Kind regards,

[Corresponding author]

## Draft GBIF clarification (for review, not sent)

For original historical occurrence-search facets that preserved taxon-level counts, query parameters, geography and date windows but did not retain occurrence IDs or datasetKey, what attribution procedure is recommended when a later provider census cannot be shown to match the original acquisition? We have a separate original GBIF download DOI with complete contributor counts; it does not cover these historical facets. We seek an accurate attribution route, not retroactive assignment of a DOI to a different source.

## Current data availability wording boundary

No stable anonymized reviewer archive or permanent deposit has yet been created. Do not insert a fabricated review URL, DOI or statement that all source data/code are available. A final Data and Code Availability Statement must be updated from actual approved deposits and access tests, preserve the original 155 numerical reproduction limitations and missing historical attribution, and distinguish source reconstruction from model replay.

No manuscript edit, no licence grant, no editorial/contact transmission, no repository visibility change.
