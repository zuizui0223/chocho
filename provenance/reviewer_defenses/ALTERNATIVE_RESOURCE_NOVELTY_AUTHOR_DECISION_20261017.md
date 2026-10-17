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
- **Main alternative figures:** candidate resource-incidence partition + conditional null, 81-species host-inventory geographic identity turnover, and specialist/native resource versus independent occurrence panel.
- **Supporting source budget figure:** original HOSTS and BCE/Clarke incidence partitions with explicit different denominators.
- **One ZIP author-review packet:** `scripts/build_butterfly_resource_novelty_review_bundle.py` generated from frozen original source receipts and pinned successful source-rendered figure Actions. Includes 3 candidate main images, supporting composition chart, SVGs, full species CSV, raw evidence receipts, a numerical source manifest and readme spelling out non-submission status.

No official journal submission, public archive DOI, new species occurrence inference or local larval-use validation occurred during this step.
