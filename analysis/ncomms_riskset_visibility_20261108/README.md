# Nature Communications: source-locked 225 × 51 butterfly temporal visibility

**8 November 2026. Result: exact source visibility, NOT a model-ready event/non-event risk set.** This branch is an observational evidence audit. It changes no original GEB/JBI results and does not submit to a journal.

## Original source population, independently reconciled

The 225 candidate butterflies and 51 original WPI-linked island/landmass groups represent 11,475 possible species×island cells. Their historical original-source margins agree independently by species and by island.

| Source observation (original GBIF *index*, not true occupancy) | Count |
|---|---:|
| All source species×island combinations | **11,475** |
| Source butterfly indexed in 2001–2010 (aggregate margins) | **524** |
| Source butterfly first indexed in 2011–2022 (aggregate margins) | **430** |
| First-index 2011–2017 cells with recovered source species+island+first-year identity | **302** |
| First-index 2018–2022: aggregate margin but identity not in this public data | **128** |
| Not previously indexed in 2001–2010, by subtraction ONLY | **10,951** |
| No index in either original window, by subtraction ONLY | **10,521** |

The 302 positive identities occur in **106 species and 48 original islands**. The residual 128 later indexed identities have 74 nonzero species margins and 41 nonzero island margins. A pair of marginal frequency tables is not a record of which species occurred on which island.

**Critical rule:** The 11,173 species×island rows WITHOUT a source-verified 2011–2017 event are not negatives. They mix 524 earlier records, 128 later records, and 10,521 not-indexed-through-2022 combinations, but the deposited source cannot identify these categories row by row. Non-indexing is not ecological absence or an effort-calibrated negative. The complete 225×51 CSV intentionally labels such states NOT_RECOVERED_IN_PUBLIC_COMMITTED_ROW_LEVEL_SOURCE. A regression fit to these rows as 0/1 would be invalid.

All 194 WPI port botanical inventories were indexed or explicitly zero-covered in the original earlier audit. A later re-query shows **45/51** source islands have an indexed 2001–2010 general plant speciesKey. General Flora IDs do not constitute pre-2011 *larval host* evidence; the pooled 2001–2024 modern host groups on the 302 positive events cannot be backdated.

## Distinct earlier 17-island / 60-butterfly risk analysis also recovered

The author previously built an explicitly at-risk source using the nonpublic original species×island×year data. Its original aggregate scientific result is retained at provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_port_island_realization_risk_20261031.json.

| Window and selected source group | Earlier host indexed | Earlier host NOT indexed | Raw recording RR | Original geographic-margin conditional null |
|---|---|---|---|---|
| 2011–2017: 905 source at-risk cells, 113 events | 34/207 = **16.43%** | 79/698 = **11.32%** | **1.45** | Observed host matches 34; expected **33.92**; 1047/1497 correlated draws ≥ observed |
| 2018–2022: 792 source at-risk cells, 33 events | 20/173 = **11.56%** | 13/619 = **2.10%** | **5.50** | Observed matches 20; expected **19.50**; 770/1497 correlated draws ≥ observed |

These cohorts are source-selected and not two independent experiments. The apparent raw positive association is almost identical to conditional reference matches preserving butterfly×regional component and island outcome margins. This is **not a supported host-driven colonisation effect** or a calibrated independent-pair p-value. Historic plant non-index is not proven host absence. This 17-island at-risk population is not the separate original 1,163 positive-record comparisons, in which only five prior-host matches were identified.

## Reproducible and fail-closed outputs

Execute the standard-library-only code:

    python -m unittest discover -s tests -p test_ncomms_island_risk_visibility.py -v
    python scripts/build_ncomms_island_risk_visibility.py --output-dir /tmp/chocho-ncomms-visibility

The script checks six original blob SHA1 identities before loading them, original 51-group/225-species keys, all 302 source positive identities, marginal totals 524/430 and residual 128, current botanical-source coverage, and the separate 17-island source conditional null. The complete CI artifact contains:

- species_island_225x51_visibility_NOT_EVENT_RISK.csv — 11,475 explicit source-visibility rows with no invented non-events
- original51_source_index_margins.csv — 51 source island margins and historical general-plant effort indicators
- original225_source_index_margins.csv — 225 source butterfly margins and missing positive identity counts
- risk_visibility_audit.json — pinned source and output file hashes, complete residuals and original 17-island results

Successful [source-locked GitHub Actions run](https://github.com/zuizui0223/chocho/actions/runs/42454155002). The original 83,573-row annual Parquet and full original named HOSTS×port×year source are not present in publicly committed data. The script must not fabricate either dataset.

## Nature Communications inference and next source gate

Existing positives remain intact: 239-species WGSRPD3 potential island resource expansion (+118.9% island vs +40.9% mainland), modern ports/entire introduced flora (+0.2445), original 51-island host/later butterfly community correspondence (rho +0.300), and the separate 155-species host-based temporal/spatial predictive gains (+4.12% / +3.00%). These are **distinct original populations and estimands**, not a joint causal chain.

The NC test — whether **time-prior accepted larval-host identity × marine accessibility** improves later butterfly record prediction under geography, climate, observation effort, and structured non-host/placebo tests — remains blocked at the raw source level. Recover the original 83,573-row global_original_species_island_year.parquet (SHA256 8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f) plus accepted species-specific pre-2011 HOSTS×island records and the original effort definition. Do not interpret NC source visibility or raw RR as butterfly colonisation, breeding, or ship-mediated plant introduction.

No merge and no submission was performed.
