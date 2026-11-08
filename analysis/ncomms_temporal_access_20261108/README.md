# Nature Communications — time-prior resource × marine-access feasibility result (2026-11-08)

**Current decision: Nature Communications is the FIRST-SHOT TARGET; the chronological host×access mechanism is NOT YET TESTABLE using the committed source snapshot.** The pre-existing JBI/GEB positive observations are preserved without revision. This result is a completed reproducibility/source-coverage diagnostic rather than a new experimental or statistical ecological result.

## Successful frozen-source execution

- [Source-locked CI, 9 regression tests and full audit: SUCCESS](https://github.com/zuizui0223/chocho/actions/runs/42446879247) (GitHub Actions artifact contains exact JSON, README, and SHA-256 pinned source identities; artifact retention 30 days).
- [Frozen, explicit methodological protocol](PROTOCOL.md); [source auditor](../../scripts/audit_ncomms_temporal_gate.py); [tests](../../tests/test_ncomms_temporal_gate.py); [dedicated workflow](../../.github/workflows/chocho-ncomms-temporal-host-access.yml).
- **Status of original source and null:** all original 194 early-port botanical inventories are marked completed, 51 source landmass groups are source-indexed, and an independent original-matched comparison and family-structured null were verified. This is *not* the same as full historical host detectability.

## Additional independently archived positive prediction evidence (different estimand)

The original **155-butterfly** prediction panel is separate from both the **225-candidate/51-island** index dataset and **239-butterfly/WGSRPD3** region reconstruction. On the original **208,482** test rows, adding static host identity to the M3 background+climate+land-cover model (M5) lowers **species-equal MSE** by **4.1202%** in the temporal evaluation (0.032063883→0.030742795) and **2.9952%** in the spatial evaluation (0.032326905→0.031358648). The results remain reproducible against the original pinned `host_island_synthesis_20261030/forecast_temporal_results.json` and `forecast_spatial_results.json` and are a **positive result we do not discard**.

**Critical distinction:** this forecasts a *rarefied/standardized butterfly recording probability*, not confirmed colonisation; the host distribution variables do **not** document a plant introduction before a butterfly record. A temporal split of the butterfly observations is not a proof of historical host arrival. The separate earlier-historical-host test below is the mandatory NC mechanism gate.

### Quantitative bottleneck — same 1,163 later-recorded butterfly×island pairs

| Indicator | Fixed source observation | Meaning |
|---|---:|---|
| Original candidate butterfly panel | 225 species | Not the separate 239-species WGSRPD3 global reconstruction |
| Species with later butterfly record pairs | 82 | Earlier/later are archive indexing windows, not colonisation |
| 2011–2022 later butterfly records in common source cohort | 1,163 pairs | Positive indexed records only, **not a full risk set** |
| Matching host plant recorded during 2001–2024 | **208 / 1,163 (17.88%)** | Retrospective host inventory; contains possible future information |
| Matching host plant already indexed in 2001–2010 | **5 / 1,163 (0.43%)** | Only five early *documented* host matches, **not absence in the other 1,158** |
| Butterfly species accounting for five early matches | **5** | *Colias croceus*, *C. hyale*, *Nymphalis antiopa*, *Pieris rapae*, *Vanessa cardui* |
| Family-preserving surrogate draws ≥ observed historical host matches | **44/199** | Descriptive structured comparison, not a causal p-value |
| Newly reconstructed early **general** Plantae records | **149/194 ports; 45/51 island groups with at least one indexed plant speciesKey** | This verifies historical general flora coverage, **not** particular larval-host identities |

A **minimum support-screening** benchmark of 50 confirmed host-prior cases across ≥10 butterflies and ≥10 island groups is not met. The benchmark is not a formal statistical power calculation. Because we have only five documented prior-host matches *among positive later butterfly records*, the source does **not** support fitting a definitive interaction model for record-based colonisation, nor treating all remaining pairs as host absences. The source also lacks all nonrecord/effort-stratified at-risk butterfly×island combinations needed for event prediction.

## Preserve original positive evidence, but keep estimands separate

- **239-species WGSRPD3 resource context:** potential island host opportunities +118.9%, mainland +40.9%, area/biogeography-matched +10.7 percentage points.
- **98-region WPI/flora association:** introduced total plant β +0.2445; distinct introduced host β +0.1639; host-specific conditional β −0.0415.
- **51-island original candidate records:** host/later butterfly Jaccard Spearman +0.300; LOIO exploratory host gain +1.7312% MSE, marine-circuit gain +6.4880%. Stricter two-unseen-endpoint support and adverse AIS local placebo are included in [PR #151](https://github.com/zuizui0223/chocho/pull/151).
- These results describe a real *data-defined* resource pattern, but cannot be promoted into a demonstrated source-host-before-colonist mechanism through writing alone.

## Why existing data cannot be silently recycled for the NC main claim

1. The **2001–2024 candidate introduced-host inventory** overlaps and can postdate the 2011–2022 butterfly-index response. Its 208 corresponding pairs are not independent time-prior resources.
2. The **2001–2010 general GBIF Plantae speciesKey source** (27,598 port × speciesKey rows) does not identify HOSTS larval plants by itself. Adding it as an early host would introduce a taxonomy and treatment-classification error.
3. A model trained/evaluated on **1,163 observed later butterfly cases alone** cannot estimate species×island event probabilities or counterfactual non-events without the missing at-risk population and observation effort.
4. **2015–2021 AIS position density** is not a pre-2011 historical ship-link predictor and does not identify bilateral transport; scalar wind-speed circuits also lack direction.
5. The published 2026 Nature Communications study by Gross, Kawahara and Daru (*Climate and regional plant richness drive diet specialization in butterfly caterpillars*, DOI 10.1038/s41467-026-73236-4) already establishes a global plant richness/climate/diet specialization account. Our distinguishing test must be a **dynamic redistribution and accessibility finding**, not repackaging food-plant richness.

## Work that actually unlocks the Nature Communications question

**Main target**: construct an independently reproducible `butterfly species × island × time` risk-set table with observed pre-baseline index status, observational effort, and verified **named** historical host-plant identities on each source island before each butterfly outcome window. Dated herbivore resource must be measured *before* the butterfly target, with observation missingness separated from confirmed non-records. Source biological records and hosting must be subject to public data or review-ready controlled-access release.

Compare a fixed set of models—baseline geography/climate/effort, time-prior hosts alone, source-only uniform-ocean access alone, additive, and **prior host × access**—on the same held-out samples, including butterflies and islands absent from training. Compare conditional log-loss and calibration under suitable support, blocked by both geographic and taxonomic structure. Challenge the interaction with degree/family-constrained host rewiring, non-host plants, survey-index placebos and distance/AIS source controls. This is an **observational advance**, not proof of breeding/arrival unless independently supported life-stage records are acquired.

**Publication checkpoint:** Do not replace the JBI version or submit NC until the temporally prior host support and risk set are materially improved, the interaction shows independent predictive value rather than an assumption-driven sign, all original-source data/code are review-accessible, and manuscript overlap with the GEB 239-species study is explicitly agreed with coauthors.

## Journal constraints (verified against current official pages)

Nature Communications expects a genuine new disciplinary advance supported by strong evidence and does **not** accept pre-submission enquiries. The paper must not significantly overlap an under-review or in-press companion paper. The minimum data and core code necessary to inspect the claims must be available with a transparent Data/Code Availability statement. See [guide to authors](https://www.nature.com/ncomms/submit/guide-to-authors), [peer review](https://www.nature.com/ncomms/editorial-policies/peer-review), [reporting standards](https://www.nature.com/ncomms/editorial-policies/reporting-standards).

**No journal submission or PR merge was made.**
