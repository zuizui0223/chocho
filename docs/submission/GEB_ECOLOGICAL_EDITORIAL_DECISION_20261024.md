# GEB editorial choice and initial-submission handoff — 24 October 2026

**Decision:** advance the ecological resource-homogenization version as the *single* GEB Research Article candidate, without merging it into `main` or submitting it on behalf of the authors. This branch adds a bounded absolute-null-effect wording correction and a separately tested initial-submission gate. Author identities, declarations, acknowledgement wording and final submission remain author-controlled.

## Selected article and ecological question

**Working title:** *Plant globalization expands and homogenizes butterfly larval-resource geography*.

**Question:** When humans redistribute plants already documented as butterfly larval hosts, does the spatial organization of potential larval resources change for consumers and between regions, beyond a larger host map?

Three supported levels should be distinguished:

1. **Resource geography gain**: fixed 239 butterflies, 206 gains, +14,553 butterfly×WGSRPD3 mapped candidate resources (+54.9% versus native 26,530). Extent alone follows union logic and is not the primary surprise.
2. **Composition and geographic co-use**: across 355 native-active geographical regions, mean potential butterfly-resource assemblage Jaccard increased from 0.27698 to 0.46208 (+66.8%). The same constrained row/column butterfly×region margin null expects about 0.45348; only **+0.00860 absolute Jaccard** is left beyond that null (nominal 499-state Monte Carlo tail 0.002). Geography overlap among all 28,441 butterfly pairs increased from 0.22237 to 0.32788 but exceeded a null median 0.32145 by only **+0.00643**. These small residuals must accompany the raw relative effect whenever used as novelty evidence. The null does not isolate causal introductions, botanical host quality or actual consumer communities.
3. **Exact host-sharing space**: already host-sharing butterflies increased their modelled pair×region shared-resource units from 62,473 to 141,885 (+127.1%), including between-family opportunities. These are *potential co-use configurations*, not competition, demonstrated co-occurrence or new larval host acceptance. Of the 14,553 added units, 8,574 (58.9%) had only one documented introduced host support; no actual population resilience follows.

The most defensible butterfly-ecology contribution is **the geographic and compositional reorganization of candidate consumer resources when interaction partners move but fixed consumer-host links do not**. The theoretical generalization is about resource opportunity and scaling, not demographics. An alternative main narrative led by geographic novelty of candidate host records (draft PR #43) is comparatively source-inventory-sensitive and not demonstrably a more mechanistic ecological discovery.

## What is intentionally not elevated
- The 206-species negative association between log resource gain and the share of singly supported added regions is **nonexceptional** under a new qualified 2×199-draw botanical plant-family×continent graph null (observed adjusted beta −0.22072; null median −0.22215/−0.22226), including the 191-species source-adequate subset. It is a source-compositional statistic, not exceptional functional buffering. Refer to [separate null verdict](https://github.com/zuizui0223/chocho/blob/analysis/gain-support-plantnull-20261024/docs/exploratory/GAIN_SUPPORT_STRUCTURAL_NULL_VERDICT_20261024.md) but do not insert the exploratory branch into blinded reviewer materials.
- Disturbance-affinity gain advantage: positive in 183 butterflies but reduced in host-taxonomy adequate subset. No new general rule.
- Recent full-index GBIF country-level butterfly co-recording: +3.91 percentage points, genus bootstrap CI −6.42 to +14.94 with only 21 informative strata. Different grain from the current manuscript's restricted climate-stratified 66/115 presence result; neither confirms human-introduced hosts caused actual butterfly colonization.
- GloBI pilot: only 5/20 selected butterflies produced stage+plant+coordinate+trophic candidate source rows across two source IDs, and some coordinates derive from historical regional centroids. No general field larval feeding or recruitment validation.
- The 38 most-contributing host removal exercise is intentionally target-ranked. Its large modeled loss is a network stress, not measured threat from invasive plant removal.
- No new systematic posthoc predictor hunt, null retuning or separate butterfly ecological paper is needed.

## Verified editorial changes in this candidate
The selected PR #38 version originally advertised raw +66.8% regional and +47.4% consumer-envelope similarity increases alongside nominal fixed-margin tails. This audit adds **+0.0086** and **+0.0064 absolute Jaccard residuals** directly to the structured Abstract and corresponding cover letter so the effect is not represented as mostly unexplained. Exact source metrics, figures, results tables, botanical and butterfly identities, and journal title remain fixed.

Current [GEB author guidelines](https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html) require a double-anonymous separate title page; Research Articles typically have about 5,000 main-body words, structured Abstract <=300 words, one corresponding author and a separate cover-letter paragraph of <250 words on readership relevance. This candidate's independent preflight measured:
- Abstract **284** words.
- Main body **4,884** words.
- Cover-letter interest paragraph **187** words.
- Running title **39** characters; **7** keywords.
- Separate editable Supporting Information DOCX, anonymous data-and-code ZIP, and BCE/Clarke Table S9 accepted-host source ZIP paths all mapped.

Exact test run: [GEB ecological manuscript preflight](https://github.com/zuizui0223/chocho/actions/runs/40083643559) returned **PASS_SCIENTIFIC_FORMAT_PENDING_ONLY_AUTHOR_IDENTITIES_DECLARATIONS**, with *exactly six author-controlled placeholders* and no other blockers under the repository's frozen format and archive checks.

## Outstanding author-controlled items before actual submission
1. Final author order, names, institutions, postal/contact details and ORCID for the submitting author and intended coauthors, with **one** corresponding author.
2. CRediT roles approved for each author.
3. Acknowledgements approved by contributors, and complete funding/grant statements.
4. Conflict of interest declaration confirmed with all coauthors.
5. Corresponding-author signature/contact details in the cover letter.
6. Formal author review/approval of the final text and its submission, plus actual portal upload and file inspection.

Do not infer these from a GitHub author name, a student profile or an email address, and do not state that any coauthor has agreed until they have.

## Exact upload outputs after approval
- The double-anonymous, line-numbered **main DOCX** and separate native-table **Supporting Information DOCX** from the same successful [GEB blinded-DOCX workflow](https://github.com/zuizui0223/chocho/actions/workflows/build-blinded-review-docx.yml).
- Anonymous data-and-code ZIP from the successful [anonymous review bundle workflow](https://github.com/zuizui0223/chocho/actions/workflows/build-anonymous-review-bundle.yml).
- Separate BCE/Clarke **S9 source and reproducibility ZIP** from the same anonymous review bundle action as mapped in `manuscript/geb_initial_submission_manifest_v0.1.json`. Archive/download both reproducibility ZIPs before transient Actions retention expires.
- A separate identifying title-page document and signed cover-letter PDF supplied by authors.
- GBIF occurrence download DOI `10.15468/dl.pp5nc9` is already recorded. A public study-data/code DOI and public release tag are deferred until publication according to current repository manifest.

## Release boundary
This review does **not** submit to Wiley, merge PR #38, change `chocho/main`, invent author approval, or publish an identifying reproducibility link in the blinded files. Treat a green CI as technical readiness, not completed submission.
