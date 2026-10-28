# Prospective climate × resource-season opportunity scenarios, 2026-10-28

**Phase**: CMIP6 source-backed exploratory thermal projection, never a claim about realized population growth, evolved climate niche, host switching or number of actual monarch generations. It reuses **Freedman (2026)** public annotated monarch larval host records for geographic sampling anchors. Freedman's winter larval geographic findings are prior art; this study does not claim novelty from their 209/348 winter nonnative-only grids. The exploratory branch must not alter the GEB manuscript.

## Source and immutable filters

- Original public author ZIP SHA256 `4e4def63669b34738bcbfb8ba32ccf6c61d02cc6992b8f1ebfa8bc350a3987db` (Freedman 2026, DOI `10.64898/2026.04.06.716508`). Original exact 30,320 `core` and 12,942 `thinned` photo records. Match `thinned.uri` to `core.uri` only in transient runner memory. Reject ambiguous URIs, absent coordinates, known coordinate uncertainty >20 km, unresolved plant binomial, reared photos, obscured source setting, and ambiguous host origin/cultivation and E/W. No new taxon, region, time window or record inclusion choice based on a climate outcome.
- All non-identifying output summaries weight **distinct original 100-km study hex locations equally** (no per-photo pseudo-independence). Each hex is represented by the median of original photo coordinates for extraction at 10 arc-min resolution. Hex keys, individual photo, coordinates, map references, observer IDs and source plant observation rows are never published or committed. Multiple 100-km hexes sharing one 10m climate pixel do not become distinct independent climate realizations; additionally report number of distinct climate pixels.
- Climate source: [WorldClim v2.1](https://www.worldclim.org/), baseline 1970–2000 and downscaled CMIP6 10-minute **same-product** monthly tmin/tmax. Frozen 2 GCMs: **MIROC6, IPSL-CM6A-LR**, two pathways **SSP1-2.6, SSP5-8.5**, two periods **2041–2060, 2081–2100** = 8 future GCM×SSP×period outcomes plus baseline. No substitution or dropping source files after observing results; fail-closed on missing/invalid downloads, grid offsets, source TIFF units or uncertain source joins. Keep original archive names and SHA256 for all raster ZIPs in derived receipts, but never export climate pixels.

## Biological thermal thresholds and limits

Developmental thermal degree-days are based on **Zalucki (1982)** and **Neupane et al. (2022), Ecosphere DOI: 10.1002/ecs2.3874**: lower `11.5 °C`, capped at `33 °C`, approximate **323 °C-days** from egg to adult. Calculated from monthly climatological midpoint `(tmin+tmax)/2` as `sum(days_in_month × max(0, min(monthly_mean,33)−11.5))` over November–February (30+31+31+28=120 days). Monthly means erase daily variability and cannot predict actual survival or full developmental completion; **323 threshold is a potential thermal-input screen, not an estimated realized annual extra generation**.

Three independent host leaf persistence scenarios apply to each identical sampled physical climate location:

1. Leaf present all four winter months = **maximal thermal-only opportunity** (323 GDD screen).
2. Leaf present for **half** winter temperature accumulation; a transparent hypothetical leaf-fraction intervention, not assigned to native/non-native/cultivated labels.
3. Leaves restricted to **November–December** or **January–February**. The absence-of-leaves case has 0 supported development even when thermal input is sufficient. Report `323/winter_GDD` as the **minimum host-effective winter-thermal fraction required** (undefined when thermally impossible).

Separate warm-season caution: record whether **mean monthly summer Tmax** exceeds `38°C` in any June–August month. This is **NOT** a model of extreme hourly heat exposure, mortality or known lethal event; constant-summer mean Tmax is much coarser than the physiological exposure thresholds described in Yang et al. (2022) Ecology and Evolution DOI: 10.1002/ece3.9039.

## Primary estimands, conditions and falsifiers

- Main panel: fraction of *geographically sampled, previously recorded monarch-larval* 100-km hexes with winter GDD screen ≥323 at baseline vs each GCM×SSP×future, paired at **fixed** source hex coordinates. `n hex`, `n distinct climate pixels`, mean/median winter GDD, required leaf persistence, optional summer monthly-heat screen. At least **30** linked valid original source hexes and 8 complete future scenarios are required; otherwise no numeric forecast. All group summaries (native vs non-native host records, cultivated yes/no, observed winter larval photos) are **descriptive geography-conditioned views**, not treatment effects, wild host occupancy or causal prediction.
- This study measures thermal opportunities at **current recorded larval localities only**. It does not model host plant geographic migration, host abundance, egg-laying, fecundity, parasite and disease impacts, diapause, migration or adult arrival; therefore **no future monarch distribution map, population count, colonization probability or predicted winter generations** is authorized.
- Climate baseline 1970–2000 differs from source recorded larval photo years (roughly 2000–2025); no inference about observed year-specific climate or phenological fitness can be made.
- The GCM/SSP ensemble is only a small scenario set and does **not** define probabilistic prediction intervals. No ad hoc change to 323, 11.5/33, years, spatial sampling, or host leaf scenarios in response to favorable numerical output.

**Biological follow-up gate:** independent **site×year** records of milkweed leaf availability and overwinter egg→adult stage completion on native and introduced/cultivated hosts with effort-controlled absent plants are required before a real climate-mediated population claim. This analysis only quantifies the potential intersection of climatic thermal energy and **hypothetical** host resource persistence.
