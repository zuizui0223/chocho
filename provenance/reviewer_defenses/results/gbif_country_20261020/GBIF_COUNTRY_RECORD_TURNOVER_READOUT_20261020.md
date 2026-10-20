# Two eight-year GBIF windows do not resolve host-enabled butterfly establishment

**20 October 2026 | source-frozen exploratory country-grain recording-turnover analysis.**
This is not a colonization, local host use, recruitment, survival, or host plant
introduction-date test. The 2010–2017 and 2018–2025 windows have equal duration,
but GBIF observer participation, reporting, taxonomy, and verification vary.

## Exact comparison

Reuse the same full-index 229-taxon country facets (33,817,891 accepted-basis
records) and exactly the fixed 255 species × botanical Level-1 resource-contrast
strata (3,188 candidate cells). No new organism records or resource labels were
generated. Filter only on the number of **other** butterfly taxa recorded in
each country in **both** time windows. The focal species never qualifies its
own reporting effort. Compare newly *indexed* country detections, not arrivals
of butterflies or their hosts.

Primary: accepted bases, at least three other taxa per country in *each* time
window, focal observation threshold of at least one record:

| Paired-window cell state | Added-host country | Control country | Total |
| --- | ---: | ---: | ---: |
| 00: not recorded in either period | 1,323 | 1,315 | 2,638 |
| 01: recorded in 2018–2025 only | 9 | 7 | **16** |
| 10: recorded in 2010–2017 only | 5 | 2 | **7** |
| 11: recorded in both periods | 11 | 12 | **23** |
| All cells passing paired reporting gate | 1,348 | 1,336 | **2,684** |

There were 105 countries with adequate other-butterfly recording and at
least one qualifying candidate cell. The earlier pooled analysis contained
48 positive candidate cells; after applying paired-window recording gates,
46 are recorded in at least one window. This does not add a biological absence
class to GBIF records.

Within the 231 retained resource-contrast strata after early-unrecorded
conditioning, 01 counts are **K=16**, added-host overlap **O=9** and conditional
expected added-host overlap **8.096060475**. Hence the 01 excess is **+5.65
percentage points** (genus-bootstrap 95% interval **−15.02 to +22.04 pp**,
1,999 draws); exact conditional geographical upper tail **0.3943**. Only
**10 informative strata** and **10 positive species** remain, below the
50-stratum/25-species/25-positive-cell information gate.

The reverse 10 comparison has **K=7**, **O=5** and conditional expectation
**3.856327006**: **+16.34 pp**, genus interval **−21.82 to +57.84**, only
**6 informative strata**. Directional differences are not demonstrably
consistent with colonization, and the forward and reverse denominators are
different. All 12 pre-specified settings fail the information gate.

## All 12 source and reporting sensitivity settings

The two partially overlapping basis scopes show **identical country-level
presence states, conditional point estimates and confidence intervals** at
every paired effort/focal threshold when using common resampling streams.
They are a robustness check of the binary record response, not independent
source replication or additional data.

| Setting | Late-only K | Excess pp | Informative | Early-only K | Excess pp | Informative |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `all_accepted_bases/other1/focal1` | 17 | +8.92 | 11 | 7 | +17.27 | 6 |
| `all_accepted_bases/other1/focal5` | 7 | +17.52 | 5 | 1 | -16.67 | 1 |
| `all_accepted_bases/other3/focal1` | 16 | +5.65 | 10 | 7 | +16.34 | 6 |
| `all_accepted_bases/other3/focal5` | 7 | +17.01 | 5 | 1 | -16.67 | 1 |
| `all_accepted_bases/other5/focal1` | 16 | +5.66 | 10 | 7 | +15.43 | 6 |
| `all_accepted_bases/other5/focal5` | 7 | +17.01 | 5 | 1 | -16.67 | 1 |
| `core/other1/focal1` | 17 | +8.92 | 11 | 7 | +17.27 | 6 |
| `core/other1/focal5` | 7 | +17.52 | 5 | 1 | -16.67 | 1 |
| `core/other3/focal1` | 16 | +5.65 | 10 | 7 | +16.34 | 6 |
| `core/other3/focal5` | 7 | +17.01 | 5 | 1 | -16.67 | 1 |
| `core/other5/focal1` | 16 | +5.66 | 10 | 7 | +15.43 | 6 |
| `core/other5/focal5` | 7 | +17.01 | 5 | 1 | -16.67 | 1 |

Intervals for every estimable outcome overlap zero; where only one positive
cell remains, a genus-bootstrap interval cannot reliably be estimated.
The sign and size of an individual point estimate therefore do not establish
a temporal shift.

## Posthoc geographic dependency (NOT a new biological test)

| Region removed | Late-only K | Excess pp | Informative strata | Early-only K | Excess pp |
| --- | ---: | ---: | ---: | ---: | ---: |
| exclude_caribbean | 12 | -6.36 | 8 | 7 | +23.53 |
| exclude_central_america | 11 | +15.03 | 7 | 5 | +4.80 |
| exclude_caribbean_and_central_america | 7 | -6.98 | 5 | 5 | +13.94 |
| exclude_level1_8 | 6 | -5.37 | 4 | 2 | +24.85 |

Caribbean removal reverses the late-only point estimate, from **+5.65** to
**−6.36 pp**, but both intervals overlap zero and information drops from
10 to 8 strata. This mirrors the country-level pooled-result geographic
dependency; it is *not* proof that host-plant effects differ among continents
or that a specific invasion pathway produces more butterflies.

Late-only positive country counts: DE: 1, DO: 1, EC: 1, GB: 1, GT: 2, HN: 2, HT: 3, IL: 1, LB: 1, NI: 1, RW: 1, SE: 1.
Early-only positive country counts: DE: 1, GF: 1, GT: 1, GY: 1, HN: 1, NA: 1, VE: 1.
These are indexed country detections, not wild demographic events.

## Posthoc minimum later-indexed-record strength (distinct from prespecified 5/5 binning)

The fixed primary requires **zero** earlier 2010–2017 records. Holding that
early-zero criterion and all matched country candidates unchanged, we raised
only the minimum number of 2018–2025 indexed records. This does *not* assume
independence of GBIF records, wild reproduction or plant arrival dates.

| Later record threshold, earlier still zero | Qualified late-only cells | Added-host positive cells | Conditional expected | Excess pp | Informative strata | Genus 95% interval pp |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| >=1 (primary) | 16 | 9 | 8.096 | +5.65 | 10 | −15.02 to +22.04 |
| >=2 | 7 | 3 | 2.364 | +9.08 | 4 | −21.43 to +36.67 |
| >=5 | **2** | **1** | 0.467 | +26.67 | **2** | not estimable |

Exactly **9/16 late-only cells have only one record**, and **6/9 positive
added-host cells** have only one record. The increase in point estimates at
stricter thresholds is a consequence of the shrinking, highly selected
denominator, not stronger evidence. Every setting remains well below the
prior 50-informative-strata gate.

The only two species/country cells with zero earlier records and >=5 later
records are `Aguna asander` in Dominican Republic (9 indexed records;
**added-host-resource country** under the coarse fixed catalogue) and
`Papilio demoleus` in Israel (86 indexed records; **not** an added-host
country in the same catalogue). The complete 7-cell >=2 list is retained in
[the archived summary](gbif_country_record_turnover_summary_20261020.json).
These are **targets for finer-grain occurrence and record-provenance checks**,
not evidence of recent colonization, breeding, or population support.

## Scope, verification and replay

- Frozen protocol: [GBIF_COUNTRY_RECORD_TURNOVER_20261020.md](../../protocols/GBIF_COUNTRY_RECORD_TURNOVER_20261020.md)
- Summary (12 full setting summaries, country breakdowns and posthoc regional diagnostics): [JSON](gbif_country_record_turnover_summary_20261020.json)
- Exact full JSON artifact from [successful independent GitHub Actions run](https://github.com/zuizui0223/chocho/actions/runs/39386941953), artifact ID 12248345183. It expires after 30 days; this summary and reconstruction code remain in the repository.
- Source hashes: country CSV `d9f5207cd611e737279b7b579fe8f2347132d68bab36fc78a56e727f498c0359`; original 239-name source `18f84371cf2a9a71ca06202fac2246d9aa3df0465bd973d06b4e9e8bfa58e8d2`; locked primary strata `6eaa91001b282a8bb8bb3f78847386ef5c4fbefcf273f4b09d173d6a691100ec`.
- Code: [record-turnover analysis](../../../scripts/analyze_gbif_country_record_turnover_20261020.py); local tests and new workflow succeeded on commit `48b92fdb15dd9819f51311fb3f4ca8b7de67438e`.
- Reproducibility:

```bash
python scripts/analyze_gbif_country_record_turnover_20261020.py \
  --country-counts data/exploratory/gbif_country_20261020/gbif_country_counts.csv \
  --original-manifest data/exploratory/gbif_country_20261020/original_239_taxon_manifest.json \
  --frozen-strata provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_primary_strata.json \
  --bootstrap-draws 1999 \
  --output /tmp/gbif_country_record_turnover.json
```

**Ecological boundary:** human relocation of host resources has changed
*potential larval-resource geography*, but this country-level detection
exercise does not show that host redistribution preceded butterfly appearance,
improved reproductive success, or turned host availability into viable
populations. Source size cannot compensate for only 16 positive late-only
cells. The next empirical discriminator would require subcountry locations,
dated host establishment and evidence of larval recruitment.
