# GEB original155 partial source-panel replay — 24 November 2026

**Outcome: exactly 842,964 keyed model rows and 93 comparable columns verified with zero numerical discrepancies, using the recovered original intermediate sources. This is not full raw-GBIF-to-panel reconstruction.**

## Source and operational provenance

[Source-bound GitHub Actions run 45316564720](https://github.com/zuizui0223/chocho/actions/runs/45316564720) completed successfully. The executed script `scripts/replay_GEB_v05_original155_partial_panel_20261124.py` first verifies the complete draft original archive SHA-256 `9e39b0068c3617741f9b275a1edfddd574c93d620853e956bd706ad067a53f9f` and the original file-manifest entries before extracting exactly ten required original sources, including the archived panel solely as the **comparison target**. Original panel SHA-256 is `5ae5b48b752fad05067537de5fb95d4ee74f17f3ba1114443d337bd7fd3695a1`.

The source snapshot's frozen public `forecast_all_225_disposition.csv` (SHA-256 `ae0a874b1d743657fdd95a1e653d7ca0939fc41`) supplies **stage labels only** to reproduce the original 155 eligible species out of 225. This is an independent alternate input for eligibility and is **not** the unavailable original `forecast_candidate_disposition.csv`. The original family-support selection was recomputed from the original `host_region_family_composition.parquet`; both all-listed and positive-share-only rules gave the same **18** plant families as the archived target.

## Verified comparisons

| Check | Result |
|---|---|
| Source-selected cohort | **155** butterfly species |
| Reconstructed species×cell domain | **280,988** pairs |
| Reconstructed original panel keys | **842,964 / 842,964** species×cell×period rows, sorted and aligned with no duplicate keys |
| Original↔recomputed common columns | **93**, all have zero cells outside tolerance |
| Distance from past occurrences | **842,964/842,964** within tolerance; maximum absolute error **1.4495071809506044 × 10⁻¹² km** |
| Observation proxy | Days, success, four season quarters, provider dominance and hypergeometric 10-provider-day recording probability all match **842,964/842,964**, max absolute error zero |
| Host opportunity and 18 family shares | Rejoined original upstream host-composition tables exactly for all 842,964 rows |
| Environment and geometry | Original source intermediate joins and spatial block fields exactly match |

The report uses row-key alignment, finite/NaN checks and preset numeric tolerance (relative `1e-6`, absolute `1e-8`) and records **all** compared columns, not samples. Values from the saved final panel were **never substituted as a generator input**.

## Three missing ORIGINAL input paths and defensible scope

The 3,032-entry source manifest verifies **10 of 13** exact original build-script input paths; **59** raw GBIF occurrence parquet partitions are registered, but were **not** reprocessed in this run.

- `outputs/forecast_candidate_disposition.csv` — original file missing; source-locked 225-species **published disposition** substituted **only for the predeclared stage-based eligibility selection**, with 155-species membership verified.
- `outputs/temporal_host_family_support.csv` — original file missing; **18 eligible plant families** reselected from the saved original upstream family-sharing table, then matched against the frozen panel.
- `work/remote_sensing/host_region_counts.csv` — original file absent, including zero matching `host_region_counts.csv` basename candidates in the source manifest. The **three origin-credit predictor columns** `native_host_count`, `introduced_host_count` and `introduced_share` were **excluded** from this reconstruction, not filled from the comparison panel.

Hence the result verifies the original **upstream-derived** panel components, **not** complete regeneration of the independent WCVP/HOSTS host counts, record-day source construction from the original occurrence partitions, or the full raw observation→host+environment→model workflow. It is a substantive partial reproducibility pass, not a new GBIF biological inference.

## Implications for manuscript and remaining gates

The prior wording that the 842,964-row panel was **not recovered** or **no primary refits were executed** is now historically obsolete. The full panel is recovered privately; independent temporal and spatial fresh model fits passed, the held-species fit is numerically unmatched, and this partial upstream replay verifies 93 columns. Two additional temporal sensitivity re-fits passed and two remain numerically discordant. None of those issues is erased by the panel comparison.

Raw-to-panel reconstruction of the three host-origin predictors, complete raw occurrence preprocessing, independent third-party rights and indirect double-anonymous access remain open. Do **not** export the original source rows, change frozen paper effect sizes, publish a draft release, create a fictitious DOI, mark the paper Ready, or send it for coauthor review without authorization.

[Aggregate-only Actions receipt](https://github.com/zuizui0223/chocho/actions/runs/45316564720/artifacts/14885815414) · [Original source manifest inventory](https://github.com/zuizui0223/chocho/actions/runs/45316151416)

