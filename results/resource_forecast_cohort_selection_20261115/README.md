# Selection between butterfly host-resource geography and actual recording forecasts

**15 November 2026 | Frozen source-complete cohort crosswalk; not a per-species prediction-effect test**

## Question and decision

Does the **155-butterfly recording-forecast evaluation cohort**, which supports the independent positive +4.12% later-period and +3.00% held-spatial improvement from known-host geography, resemble the **239-butterfly global host-resource reconstruction**?

**Not in every biologically relevant dimension.** The recording-predictable butterflies come from a source-selected original 225 candidate list and have more documented larval hosts, larger global potential host-resource expansions and a markedly more mainland-favoured **species-equal gap-filling fraction** than 70 original candidate butterflies not admitted to forecast fitting. The latter have insufficient recording or input support, **not measured failed forecasts**.

This finding does **not** diminish the strong original resource results, nor does it refute the observed positive MSE gains **within the eligible 155 species**. It quantifies why the published 155-species forecast gain cannot be automatically generalized to all 239 global resource butterflies.

## Original source sizes and independent evidence

- The global resource/trait table contains **239** butterfly species.
- The frozen forecast candidate table contains **225** of these, **155** fitted and evaluated under the prefit eligibility contract, **70** excluded. **14** of the 239 resource butterflies were never in the 225 forecasting candidates; they must not be treated as forecast failures.
- Same HOSTS×WCVP source catalogue and butterfly names were matched across (a) original 239 resource metrics, (b) the full 225 candidate disposition ledger and (c) the same-butterfly **126 island / 237 mainland** source resource gap-filling table.
- The 155 fitted species account for exactly **208,482** held test species×cell rows in the previously confirmed forecast comparison, not novel source data or a refit.
- Original per-butterfly M3/M5 prediction losses and the underlying prediction Parquet were **not archived in the source repository's forecast snapshot**; committed JSON/CSV supports only aggregate cross-validation losses and per-species selection support. Attempts to search the current repository and the source PR #101 did not reveal those prediction files. **No invented per-butterfly loss–resource correlation is reported here.**

## Cohort differences, identical source metrics and species weighting

| Original frozen species attribute | Fitted recording forecasts | Candidate but not fitted | 239 resource species outside forecast candidates |
|---|---:|---:|---:|
| Number of butterflies | **155** | **70** | **14** |
| Known larval-host plant species per butterfly, mean | **13.04** | **7.10** | **5.57** |
| Global introduced-host-added potential botanical regions per butterfly, mean | **65.92** | **55.07** | **34.36** |
| Species-equal island minus mainland native-host-gap filling | **−7.91 pp** | **+0.34 pp** | **−6.97 pp** |
| Butterflies with greater island gap-fill fraction | **34** | **30** | **1** |
| Butterflies with greater mainland gap-fill fraction | **110** | **23** | **8** |
| Native-based pooled island host-resource increase within the cohort | **+116.40%** | **+128.39%** | **+96.88%** |

Within the original 225 candidate panel, the modeled-minus-unmodeled **species-equal gap-filling difference is −8.247 pp**; its **4,999-draw, within-cohort-species bootstrap** 95% interval is **−11.644 to −4.813 pp**. This is **descriptive**; shared host networks and phylogeny make a simple independent-species inferential p-value unwarranted.

The source selection is also strongly associated with known-host catalogue breadth and butterfly recording support, so raw group contrasts are not a causal effect of eligibility or host-resource change.

## Do dietary breadth and native-resource differences fully account for it?

Two conditional relabeling diagnostics preserve the observed counts of forecast-included butterflies in source-covariate strata:

1. **Frozen host-family breadth stratum maintained**, all **225** candidate butterflies represented: observed model-included minus excluded gap filling = **−8.247 pp**, matched-strata null expectation **−0.476 pp**, descriptive two-sided permutation **p=0.0002**.
2. **Host-family breadth plus quartile of global native-resource region count maintained**, **208** butterflies in common-support strata: observed difference **−7.290 pp**, null expectation **−1.444 pp**, descriptive two-sided permutation **p=0.0008**.

Thus restricting to known major diet and original resource differences **does not fully erase the cohort contrast**. It remains an observational selection difference: record density, accessible domain, survey effort, geography, climate and shared phylogeny are not exchangeable or comprehensively controlled by these post-hoc source permutations.

The second comparison excludes 17 out-of-support butterflies, and therefore differs in estimand from the full 225 contrast. The 2 permutations are **diagnostics on the same data, not independent replicated experiments**.

## Relation to the GEB main story

**Retain the ecological positive:** Human redistribution of existing known host plants expands potential butterfly resource geography **+54.9%** globally among **239 species**, and proportionally on islands **+118.9%** compared with **+40.9%** on the mainland, largely because of lower native-host baseline coverage. The historical full source supports those exact resource accounting estimands.

**Retain the observational positive, with a bound:** For **155 eligible butterflies**, host resources and host-family composition added to the climate+land-cover background improved known-species recording MSE **+4.12%** in later periods and **+3.00%** in held spatial blocks. A species-held-out comparison gave **−0.21%**, i.e. no transfer evidence to previously unseen butterflies.

**New limit:** Forecast eligibility preferentially selects butterflies with richer known-host catalogues and a distinctive mainland-vs-island host-gap-filling pattern. The MSE improvement is **an accurate group-average improvement for that selected recording population**, not a general numerical gain for all 239 reconstructed resource species. Per-species MSE losses must be independently recovered before estimating the relationship between a butterfly's global host redistribution and its predictive benefit.

**Port / marine / establishment inference remains separate:** Ports correlate with broader alien flora, but not proven host-selective enrichment after whole-flora controls. Circuitscape's original 51-island retrospective community prediction retained positive gains subject to geographically strict transfer failures. Neither the forecast selection comparison nor resource redistribution proves a butterfly's long-distance arrival or successful larval use on introduced hosts.

## Reproducibility and governance

- [Executed source-locked GitHub Actions run #43673335738](https://github.com/zuizui0223/chocho/actions/runs/43673335738) **SUCCESS**, including 239/225/155/70/14 population asserts, exact 208,482 held evaluation rows, 4,999 species bootstraps and both conditional permutations.
- [Machine-readable 239-species source crosswalk, source SHA-256 hashes, original and conditional contrast results](https://github.com/zuizui0223/chocho/actions/runs/43673335738/artifacts/14144923549).
- Code: `scripts/audit_resource_forecast_cohort_selection_20261115.py` and `.github/workflows/chocho-resource-forecast-cohort-bridge.yml`. These consume only four original frozen source CSVs; they do not refit any ML model, use the inaccessible local private workspace or claim new GBIF observations.
- Source workflow generated an artifact that was independently checked for the original cohort counts and exact manuscript risk population, not a separate biological replication.

**Editorial placement:** Supplementary Results S.A (forecast transfer boundary) and a short section 5 Limitations statement. The primary GEB Results on resource opportunity and positive known-species butterfly-record forecasting remain unchanged. The PR stays Draft pending coauthor and final archive checks.
