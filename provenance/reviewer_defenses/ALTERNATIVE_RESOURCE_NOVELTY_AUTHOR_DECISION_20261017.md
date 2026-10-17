# Editorial decision gate — one GEB butterfly ecological paper (17 October 2026)

**Status:** explicit author-level recommendation after completion of the alternative full GEB research article draft, the three-source-verified candidate figures and the detached author-review bundle. This is not a second submission or a preregistration. Do not silently update or merge current GEB submission PR #38.

## Current decision

**Outcome: author-review candidate, not automatic replacement of the frozen v0.2 paper.** The alternative is coherent enough to read as a Research Article (~4,300 main words, ~250-word seven-heading structured abstract, 25 references, 3 proposed main figures), but the question of conceptual/ecological novelty still needs deliberate editorial judgment.

Working biological question: *When does an introduced larval host produce a genuinely new geographically available resource for a butterfly, rather than another potential host where a native resource already exists?*

The source scenario supplies:
- Original HOSTS panel: H=75,346 introduced host×butterfly×WGSRPD3 candidate incidences; A=46,387 already native-supported; R=14,406 extra introduced supports in new cells; U=14,553 first supports of distinct new butterfly×region cells (H=A+R+U).
- Independent Europe-based BCE/Clarke candidate evidence (81 butterflies, 1,027 additional links): H=103,349; A=67,740; R=20,642; U=14,967. Native-overlap fraction **61.6%→65.5%** and distinct-first-support fraction **19.3%→14.5%**: different incident denominators and source definitions.
- The 81 host-augmented butterflies' initial and augmented sets of introduced-added resource cells agree on only 3,225 of 6,445 union cells (**Jaccard 0.5004**), because 1,403 old-added cells are now native-supported and 1,817 other cells are newly recognised. These are **source label changes, not temporal colonisation**.
- Numerically qualified but post-hoc original-HOSTS climate-and-1,500-km matched null: the observed native overlap exceeds null median by only **271–276 host incidences, about 0.36 percentage points** of the fixed original H. This analysis was **not rerun** under BCE-augmented hosts.
- Independent occurrence check: **66 of 115** outside-native butterfly×region cells from a separately chosen 32-butterfly GBIF panel overlap added resources; continent-matched structural null median **49**, but no verified local larval host use or demographic impact.

## What is and is not novel

**Potential ecological contribution:** a fixed-link, consumer-centred geographic distinction between added host plant occurrence, already available native-host geography and newly opened butterfly regional access; linked to a transparent estimate of how incomplete host records affect that distinction.

**Not new independently:** that more plant ranges can enlarge a union, that non-native plants have ecological effects, or that plant introductions can homogenise plant and animal interaction networks. The equality H=A+R+U is **bookkeeping**, not a tested biological mechanism. The **50% host-source identity discordance** is an empirical sensitivity of one Europe-biased catalogue, not a global ecological phenomenon or causal range shift. The **0.36-point** constrained structural-null deviation is numerically reproducible but small and source-conditional.

The current strongest paper argument is therefore an **ecologically motivated inference-limit and opportunity-geography comparison**. Avoid adjectives such as 'dramatic ecological filtering', 'unexpected butterfly spread', 'universal host-network fragility', 'increased competition' or 'confirmed anthropogenic butterflies'.

## Reviewer objections to address *before* changing the frozen manuscript

1. **Mechanism:** The measured response is potential regional resource *access*, not survival or recruitment. Does the Introduction and Discussion make the implied consumer mechanism a *future field hypothesis* rather than a result? The full draft explicitly says yes.
2. **Representativeness:** Original 239 butterflies were selected by documented native resources from 339 descriptors. BCE augmentation covers 81 with European foodplant evidence; conclusions are **not probability estimates for the world's butterflies**.
3. **Source extrapolation:** Clarke/BCE evidence is regionally documented; linking each candidate host globally to WCVP accepted plant ranges assumes the link may be transferable. It does **not** certify butterfly host acceptance in every region.
4. **Metric unit mixing:** A/H and U/H partition host×butterfly×region *incidences*, whereas resource expansion counts unique butterfly×region cells; captions and legends must state both denominators.
5. **Counterfactual robustness:** The small original-HOSTS local-null residual has no BCE-aware matched counterpart; the +3.98-percentage-point source descriptive change and +0.36-percentage-point null deviation are **different estimands**.
6. **Spatial and temporal structure:** WCVP introduced vs native status does not carry introduction time; WGSRPD3 is geographically coarse; CHELSA 1981–2010 is static; raw occurrence records and botanical reference coverage are uneven.
7. **Significance:** Serial-swap MCMC worlds are correlated and the spatial protocol was decided post-hoc. Numerical mixing passed a finite gate but is neither exact inference nor a predeclared ecological hypothesis test.
8. **Journal fit:** GEB requires broad ecological/biogeographic conceptual relevance, not simply a more reliable checklist audit. The manuscript must lead with **what resource novelty means for consumers** and avoid making data-quality differences sound like observed ecological change.

## Author's choice (only one paper)

**Option A — retain frozen homogenisation manuscript (PR #38)**, and use the source-turnover analysis in Supporting Information or as a targeted limitation. Advantages: original suite and current 3 figures already mature. Risks: +66.8% raw homogenisation is largely fixed-margin-predictable; exact shared-host exposure alone does not demonstrate an ecological process.

**Option B — replace its storyline after deliberate review with the new full alternative (PR #43)**, which leads with regional novelty and source-inventory robustness, supports the framing with the small qualified geographic-null contrast and butterfly occurrence context. Advantages: a clearer, less tautological consumer-resource question and explicit uncertainty. Risks: the strongest unexpected numerical result (Jaccard=0.500) is a one-way, Europe-biased measurement sensitivity rather than a general causal ecological effect.

**Current recommendation:** keep PR #43 Draft as a fully readable alternative and do not merge into PR #38 without explicit author choice. If ecological mechanism becomes essential for GEB, independent geographically verified larval utilisation, host quality, seasonal leaf availability or adult emergence would be more decisive than further post-hoc geographic null tuning.

## Executed deliverables

- **Blinded complete Research Article draft:** `docs/exploratory/BUTTERFLY_RESOURCE_NOVELTY_FULL_DRAFT_20261017.md`.
- **Evidence consistency and journal-length tests:** `tests/test_butterfly_resource_novelty_full_draft.py` (latest paper CI success).
- **Main alternative figures:** candidate geographic resource-incidence partition (Fig. 1), 81-species host-inventory geographic identity turnover (Fig. 2), and specialist/native resource versus independent occurrence panel (Fig. 3). The SMALL original-HOSTS climate/distance null effects occupy a SEPARATE exploratory supporting Fig. S2 rather than being bundled with the main accounting figure.
- **Supporting source budget figure:** original HOSTS and BCE/Clarke incidence partitions with explicit different denominators.
- **One ZIP author-review packet:** `scripts/build_butterfly_resource_novelty_review_bundle.py` generated from frozen original source receipts and pinned successful source-rendered figure Actions. Includes 3 candidate main images, supporting composition chart, SVGs, full species CSV, raw evidence receipts, a numerical source manifest and readme spelling out non-submission status.

No official journal submission, public archive DOI, new species occurrence inference or local larval-use validation occurred during this step.


## Journal-compliance and rendering audit (17 October 2026, subsequent author-review revision)

The completed separate figure artwork [GitHub Actions run 38995025937](https://github.com/zuizui0223/chocho/actions/runs/38995025937) and [final source-checksum author-review ZIP run 38995600935](https://github.com/zuizui0223/chocho/actions/runs/38995600935) have been verified. The inner ZIP contains the complete alternative article, **3 individually assigned main figures, 2 separate supporting figures, 81-species CSV and numeric/figure-source evidence receipts**; all 26 source-member checksums matched.

The alternative full English article has a seven-heading structured abstract (~260 words), main text ~4,268 words, 25 references, eight alphabetised keywords, and an anonymisable running title <40 characters. These are within the Research Article guidance of *Global Ecology and Biogeography* for a <=300-word structured abstract, ~5,000 main words and approximately <=50 references. The journal also requires double-anonymous review, a separate author title page, continuous line numbering, and Supporting Information uploaded separately; see its [official author guidelines](https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html).

A locally rendered **21-page editable, anonymous review DOCX** was created from the latest source-verified article with the three numbered main figures. Its source-specific GitHub account/PR links and Word core property identifiers were masked; page images confirmed continuous line numbering, no clipping and readable figures. Supporting Figures S1–S2 were separated into a second, rendered 2-page DOCX, rather than remaining embedded with the main figure set. These Word outputs are **author-review files available in the conversation** and do not constitute publication, a public preprint or a commit to the original GEB submission branch.

**Submission blocker, not silently resolved:** a separately anonymised, COMPLETE executable code and data archive and an author-identified title page still need to be supplied to Wiley, along with final author approval of the alternate storyline. The GitHub-linked SOURCE review ZIP is an *internal author provenance archive* and should NOT be uploaded unmodified as double-anonymous supplementary material, because it contains source repository paths and PR URLs. The narrative's main inferential limitation (Europe-biased host evidence and no BCE-aware spatial null) is also unresolved. No reviewer/author should interpret source-Jaccard=0.5004 as temporal ecological change or the small original-HOSTS spatial residual as new globally causal biology.

## Source redistribution and blinded supplementary-package gate (17 October 2026)

A self-contained computational source package is being prepared separately from the already rendered 21-page double-anonymous article and its figures. The source ZIP includes only source-derived accepted-host sidecars, the original 239-butterfly metrics and S1 descriptors, the partial BCE/Clarke accepted-host crosswalk, archived climate and geographic constraints, verified executable scripts, and the original HOSTS-only calibrated MCMC series. Running the ZIP's standard-library entrypoint must reproduce both original and BCE-augmented `H=A+R+U` budget counts and source-specific 81-species spatial Jaccard without internet access. The complete serial MCMC draw history is provided separately because it requires substantially longer numerical work. The ZIP is forbidden from including the original author account handle, personal emails, or any external workflow-run IDs embedded in scripts. These technical checks do not guarantee anonymity without an editorial manual audit.

**Third-party data rights require independent review before this ZIP is uploaded or publicly deposited:**

- [Natural History Museum HOSTS data portal](https://data.nhm.ac.uk/dataset/hosts) identifies the hosted HOSTS dataset's licence as **CC0**, with citation DOI `10.5519/havt50xw`.
- The [WCVP Scientific Data publication](https://doi.org/10.1038/s41597-021-00997-6) describes its GBIF data download as **CC BY 4.0**. Confirm how this licence applies to the *exact pinned WCVP v13 distribution snapshot* used in the analysis, and include any required attribution.
- The particular [TDWG WGSRPD GitHub repository](https://github.com/tdwg/wgsrpd) does not establish an unambiguous licence for its converted Kew-source GeoJSON polygon files; [upstream issue #12](https://github.com/tdwg/wgsrpd/issues/12) explicitly describes unresolved provenance and geometry-quality concerns. Do **not** assume that a file being publicly downloadable gives permission for arbitrary redistribution. Substituting a differently licensed independent boundary product would require a new comparison because it can change the representative points used in the geographic null.
- BCE/Clarke candidate accepted-host associations were extracted through publication-derived Butterfly Conservation Europe pages, not through an independently globally licensed original host-occurrence census. Attribution, source-site terms and the permitted reuse of this compiled crosswalk also need author/publisher confirmation.

**Editorial position:** Keep the computational bundle as **author-only, licence-review pending**, even if code hashes and local reproduction tests all pass. Do not merge the exploratory PR #43 into PR #38, label the artifact “publication deposited,” or send identifiable GitHub workflow links to a double-anonymous reviewer. If any source cannot be redistributed under acceptable terms, prepare a public derivative-free reproduction script and cite the authoritative download path rather than silently substituting new geographic data. The observed resource results and the manuscript's inferential limitations do not change because of these administrative checks.
