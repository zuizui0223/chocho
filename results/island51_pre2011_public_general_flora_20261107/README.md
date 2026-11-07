# Reacquiring historical general plant identity at all 194 original ports

**Status: complete, source-locked exploratory macroanalysis. 2026-11-07.** Stacked on PR #137. Unlike earlier analyses using 2001–2010 GBIF **general botanical record counts** as sampling-effort controls, this analysis directly queries 2001–2010 **Plantae speciesKey identities** at every original WPI port and tests their independent predictive contribution across the original 51 landmass groups. The early speciesKey panel is **GENERAL PLANTS, not an authenticated butterfly HOSTS network**.

## 1. Completed external data acquisition

| Source audit | Checked result |
|---|---:|
| Original WPI2019 ports queried | **194/194** |
| Original islands with at least one source port | 51/51 |
| Ports with at least one live 2001–2010 indexed plant occurrence | **149/194** |
| Ports with zero indexed live plant occurrences | 45/194 |
| Islands with one or more indexed general plant speciesKeys | **45/51** |
| Islands with no indexed general plant speciesKey | 6/51 |
| Distinct general Plantae GBIF speciesKeys, entire 194-port cohort | **4,876** |
| Unique port × GBIF speciesKey rows | **27,598** |
| Live summed 5-km port-buffer occurrence indices (overlaps possible) | 173,113 |
| Archived original summed port-buffer occurrence indices (overlaps possible) | 155,002 |
| Port-level original versus live indexed-record counts not equal | **100/194** |

Original island IDs, all 194 original botanical port IDs and exact archived WPI coordinates are source-pinned in the [protocol](PROTOCOL.md). The new 2001–2010 query uses public GBIF occurrence/search with taxonKey=6 (Plantae), PRESENT occurrences, geographic coordinates without indexed geospatial issues, a 5-km radius around each original port and paginated speciesKey facets. It was applied **without selecting ports based on butterfly outcomes**.

The live public GBIF re-query is **not automatically equivalent to the historical author's original query**, because filters/search shapes and GBIF indexing have changed; counts differ at 100 ports. The original 194-port WPI source roster and old indexed-count roster were checked by fixed Git blob hashes. Individual 5-km port search circles can overlap; the 173,113 total is a sum of searches, **not 173,113 independent occurrence records**. A GBIF speciesKey nonrecord does **not** imply the species is ecologically absent. The six islands without indexing were retained with explicit pairwise indexing-missing indicators rather than removed or described as plant-free.

### Archived source evidence

- [Full 194-port plant speciesKey × indexed count snapshot, compressed](source/GBIF_original194_port_2001_2010_speciesKey.csv.gz), with all **27,598** distinct port × taxon keys; no restricted HOSTS interaction edges.
- [Port-source query and original vs live index ledger](source/GBIF_original194_port_2001_2010_coverage.csv), all 194 ports, query radius, original IDs, response-page count and per-port raw response receipt digest.
- [Original 51-island general flora indexing and union richness](source/GBIF_original51_early_general_flora_islands.csv) and [acquisition source receipt](source/source_receipt.json), including source Git blob pins and content SHA256.
- [Independent acquisition script](../../scripts/acquire_island51_pre2011_public_flora.py) with 194-source completeness, GBIF page count/facet checks, explicit rate limit/HTTP failures and 5km source geography.

These are independently re-fetched **indexed plants observed in 2001–2010**, not actual 2001–2010 introduction events or verified host use. The source snapshot is preserved in Git rather than solely a 30-day GitHub Actions artifact. For formal publication, the public GBIF source should also be assigned an appropriate reproducible GBIF occurrence-download DOI; this script's public API result is not itself a GBIF-cited DOI.

## 2. Fixed butterfly-prediction scope and models

The existing source original **51 islands, 1,275 unique unordered island pairs**, 225 butterfly candidate panel, 194 port groups (182 marine-qualified), candidate larval-host network Jaccard pooled over **2001–2024**, candidate-host marginal richness controls, geographic distance/original port counts, four latitude/hemisphere controls, archived early all-botanical record-count controls, and original uniform/modern native AIS Circuitscape resistance are unchanged. Later target: first-indexed butterfly community similarity (2011–2022) including earlier butterfly Jaccard as a baseline predictor. Earlier target: 2001–2010 butterfly community similarity, which is **contemporaneous** with newly queried general flora and thus only context.

The fixed 15 original longitude-block pair folds omit **both endpoints** of every test island pair from training; all **1,275** pairs are predicted once per model and target, with 235 within- and 1,040 between-longitude-block test pairs. Model fitting uses train-only-standardized OLS with intercept and predictions clipped to [0,1].

Exactly six models compare the original HRE_LAT baseline (host identity, host richness, geography, historic all-plant RECORD COUNTS), original HRE_AIS_LAT (+ AIS), a new HB baseline (+ historical general flora marginal richness and explicit zero-indexing indicators), HBA (+ AIS), HBF (+ historical general flora speciesKey Jaccard), and HBFA (+ both AIS and indexed historical general flora Jaccard). The zero-index island-pair Jaccard is set to 0 as a **modelling imputation only**, flagged in all new models; the raw source zero is never called ecological plant absence.

## 3. Main later-butterfly out-of-block findings, 2011–2022

| Predictive term added to the specified baseline | All original 1,275 pairs: MSE reduction | Between longitude blocks, 1,040 pairs | Within longitude blocks, 235 pairs |
|---|---:|---:|---:|
| Original AIS (before NEW taxa) | **+1.0020%** | +1.3347% | +0.3812% |
| New earlier general-flora richness plus indexing flags beyond original baseline | **−1.2990%** | **+5.2496%** | −13.5150% |
| Earlier general-flora composition beyond indexed species richness, AIS absent | **−1.1769%** | −3.9448% | +3.1330% |
| Earlier general-flora composition beyond indexed species richness, AIS present | **−0.9068%** | −3.4915% | +3.0240% |
| **AIS beyond new historical general-flora richness and index coverage** | **+1.2102%** | **+2.1264%** | −0.2159% |
| **AIS beyond historical general-flora composition, richness and index coverage** | **+1.4740%** | **+2.5533%** | −0.3299% |

All percentages are *out-of-block MSE reductions against the reference model listed in that row*, not additive percentages and not causal effect sizes. Positive signs reduce the held-out MSE; negative values worsen it. The [complete 48-row comparison table](comparison_scope_summary.csv) also preserves MAE and the source-only **990-pair subset** where both endpoint islands have indexed historical general plant taxa. In that 990-pair subset, AIS beyond general flora adds **+1.9391% MSE**, while the general-flora identity Jaccard itself reduces prediction accuracy under the fixed controls (**−4.1565% MSE when AIS absent**).

**Positive result:** the conditional AIS increment is not absorbed by directly measured **historically indexed general-plant species identity**. It remains +1.4740% MSE overall, including +2.5533% for between-longitude-block source pairs, under this fixed model family.

**Important negative/heterogeneity results:** historical general plant identity does not itself improve all-pair MSE and worsens it between longitude blocks. The indexed general flora marginal richness/coverage helps between-block prediction but worsens within-block and all-pair MSE. The simpler original AIS model still has a lower **absolute all-pair MSE 0.00719853** than the new richer AIS-plus-early-flora identity model **0.00734268**. This analysis therefore does **not** demonstrate that adding early flora makes the globally best-performing model better; it demonstrates a positive *conditional AIS addition* to a prespecified, richer botanical baseline.

The source-to-community units are novel direct GBIF source index taxa, but none of these added taxa is classified as a historically **introduced butterfly host** here. The original 2001–2024 host model remains separately retrospective.

## 4. Earlier first-indexed butterflies, 2001–2010

Earlier butterflies share the **same 2001–2010 interval** as the directly queried botanical speciesKeys, so this is not a plant-first-before-butterfly test. Negative results remain explicit: original AIS **−8.1571%** all-pair MSE; early botanical richness/recording flag terms **−20.5301%**; early floristic Jaccard beyond new richness **−3.1343%**; AIS beyond full early general flora **−4.4447%**. These source-time comparisons cannot identify temporal turnover of actual ship movements, native vegetation, insect larval feeding, or colonization.

## 5. Reproducibility and boundaries

The full live extraction + original target fit completed via [GitHub Actions](https://github.com/zuizui0223/chocho/actions/runs/42361257543). Source checks require all **194 WPI ports**, complete public GBIF facets and original port/island identity. No outcome-selected filtering was performed. The original AIS-only comparison reproduces the prior PR134 MSE/MAE scores to a maximum absolute discrepancy **4.79×10^−13**.

An independent Python-standard-library [auditor](verify_early_flora.py) reconstructs all **27,598 archived port×plant speciesKey rows**, 51 island taxon unions, all **1,275 historically indexed general-flora pair Jaccards/richness/missing flags**, all 15 held-out geographic fold memberships and original butterfly source pairs, all **15,300 held-out predictions**, 180 model-fold errors, 180 contrast-fold errors and **48 scope results**. Full prediction rows and [model/contrast fold losses](heldout_predictions.csv.gz) are archived in this branch. A new [network-free offline workflow](../../.github/workflows/chocho-island51-early-flora-offline.yml) regenerates all source-derived predictor and model/contrast rows from the saved public GBIF source snapshot, without querying GBIF again.

**Unresolved inference limits:** AIS is a modern (2015–2021) vessel-activity-derived undirected maritime resistance, **not observed shipping voyages or cargo carriage**. The original candidate-host Jaccard pools 2001–2024, so some candidate host records can POSTDATE the later butterfly outcome; conditioning on earlier *general* plants does not repair this taxon-specific time-ordering gap. The new 2001–2010 plant indices are GBIF historical observation dates, not verified introduction dates, independent whole-island plant censuses or classified larval hosts; the six indexed-zero islands are NOT biologically empty. First butterfly indexing in GBIF is not butterfly breeding, feeding, true arrival or colonization. Longitude-block folds are all from the SAME 51-island cohort; there are no naïve independent-pair p-values/CI or causal ship-mediated dispersal claims.

The earlier positive global GEB +54.9% potential larval-host resource portfolio expansion, prior island-host associations, earlier 17-island historical host-name matching, negative butterfly re-recording associations, and the PR137 local placebo result remain in their original distinct source populations and are neither removed nor overwritten.
