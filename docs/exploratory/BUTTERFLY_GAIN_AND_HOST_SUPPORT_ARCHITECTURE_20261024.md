# Potential resource expansion and host-support multiplicity across butterflies

**Exploratory ecology result — 24 October 2026.** Independent exact-source GitHub Actions runs completed. **Not a demonstrated butterfly demographic buffer or a second standalone paper.** The frozen GEB submission manuscript and PR #38 are unchanged.

## What was actually tested

An initially *positive fragility* directional prediction was recorded before estimating its coefficient: larger proportional new resource geography should depend more frequently on a single introduced known host within each newly opened region. The prediction was contradicted by the source-locked species-level test.

Define each butterfly's relative potential resource gain as
`G = log(contemporary known-host regions / native known-host regions)`.
For any butterfly with at least one newly supported region, define
`F_single = 1 − fraction of added regions with >=2 documented introduced hosts`.
Both are deterministic projections from the same original HOSTS/WCVP plant distribution data. `F_single` is not wild plant occupancy, direct larval feeding, host quality or local extinction vulnerability.

Three SHA-pinned frozen CSVs provide resource counts, host-support metrics and source-descriptor metadata. Exactly **239** butterfly species join by scientific binomial; **206** have new opportunities, representing **124 genera**. Their native and contemporary resource-cell arithmetic reproduces the **14,553 added butterfly×WGSRPD3 combinations**.

## Main estimate: source-conditional negative gain–single-host association

Species-weighted OLS response `F_single`, predictor `G`, with `log(native_resource_units)`, `log1p(resolved_host_species)`, and `host_family_count` controls. Genus-cluster resampling 5,000 times, fixed seed; those intervals are not phylogenetically independent confidence guarantees or causal p-values.

| Fixed analysis | n butterflies | n genera | Adjusted slope of single-host support share versus log resource gain | 95% genus-bootstrap percentile interval |
|---|---:|---:|---:|---:|
| **All gain-positive species** | **206** | 124 | **−0.2207** | **−0.2984 to −0.1647** |
| Host taxonomy adequate | **191** | 119 | **−0.2711** | **−0.3690 to −0.1967** |
| Secondary alternative outcome, maximum single contributing plant fractional credit | 206 | 124 | **−0.1501** | **−0.2154 to −0.0935** |

The raw species-rank correlation between resource gain and `F_single` is **−0.200**, while the covariate-partial residual Pearson correlation is **−0.409**. All **5,000/5,000** genus-bootstrap fits were successful in each model. Leave-one-genus-out slopes remain negative (full model **−0.2406 to −0.2114**).

**Interpretation of β:** Within this chosen regression, a doubling of the contemporary/native resource-region ratio (`log(2)` change) corresponds to an adjusted **15.3-percentage-point lower recorded single-host share** (bootstrap-transformed interval ≈11.4–20.7 percentage points lower). This is a cross-species slope interpretation, **not** the predicted outcome of experimentally doubling host supply for one species.

Descriptively, the low versus high third of the gain-positive species have mean `F_single` **76.9% versus 68.1%** (posthoc descriptive quantile strata only). Species-median `F_single` is **68.8%**. Independently, the original reconstructed added butterfly×region cells are **8,574/14,553 = 58.9% singly supported**; pooled-cell and species-equal metrics have different denominators and must not be interchanged.

**Result:** The prespecified claim that geographic resource gain accompanies higher reconstructed single-host dependence is not supported; a negative within-source association is reproducible. This does **not** imply that the currently expanding butterflies enjoy robust functional host substitution.

## Independent host-catalogue extension: a related but DIFFERENT response

A source-pinned 2026-10-17 BCE/Clarke supplementary archive (original [Actions run](https://github.com/zuizui0223/chocho/actions/runs/38988757959), artifact 12084636746) added **1,027 candidate accepted plant–butterfly host links** to 81 previously fixed butterflies. This extension is heavily European in ascertainment. Full geographic source identities remain unchanged.

The BCE archived per-butterfly `H=A+R+U` results identify
`support intensity = (H−A)/U = 1 + R/U`
for species with U>0. This is *mean introduced host incidence count per new region*, **not** the fraction of regions with ≥2 supports and **not** a direct confirmation of the main `F_single` slope.

Within the same **70** butterflies with a positive U under **both** host catalogues (**53 genera**), we compared the separate regressions `log(support intensity) ~ log(relative gain) + log(native regions) + log1p(source-specific recorded hosts) + host-family count`:

| Paired 70 species | Adjusted slope | 95% genus-bootstrap percentile interval | Median supports per added resource region |
|---|---:|---:|---:|
| Original HOSTS inventory | **+0.2975** | **+0.1766 to +0.5719** | **1.44** |
| HOSTS plus BCE candidate links | **+0.5773** | **+0.2948 to +0.9756** | **2.21** |

Both support-intensity slopes are positive. In **84.3%** of the paired 70 species, the source-derived support-intensity measure increased after catalogue enrichment. Yet the correlation between the *within-species change* in log resource gain and the *within-species change* in log support intensity was **−0.019**, approximately zero. Those changes do not substantiate a causal within-species complementarity mechanism.

The 81-species group totals show original average supports per new region **2.13**, versus BCE-augmented **3.27** (source-incidence weighted). The original resource-geographic novelty count changed from 4,628 to 5,042 due to simultaneous additions and native-baseline reclassifications. Never treat the increased support count as field-confirmed plant abundance or viable larval fallback. Do not call this an independent global reanalysis of plant colonization.

## Ecological contribution, boundaries and reviewer guidance

The strongest **source-conditional** insight is that *geographic extent of potential resource gain* and *multiplicity of known-host support* need not trade off across butterflies. In this catalogue, stronger expansion coincides with higher mapped redundancy, even while many individual added regions have only one recorded host. The two perspectives coexist at different aggregation levels.

**What remains unidentified:** the extent and support multiplicity are both built from the same host-resource unions. Shared geography, introduced incidence budget, correlated plant distributions, sampling completeness and taxonomic clustering can produce such slopes without ecological selection. The genus bootstrap does not construct a plant-network preserving null or recover actual caterpillar dynamics. Adding BCE hosts also changes the native-resource reference, so shifts in `U` are source reclassifications rather than time series. Any conservation prescription involving invasive-host removal requires verified local native and introduced plant supply, oviposition, survival and suitable replacement hosts.

**Recommended manuscript use:** Add this as a tightly bounded *supporting* resource-architecture analysis to the **single** GEB resource geography paper, not a new stand-alone ecological discovery. The primary paper's resource expansion across specialization and nonuniform host contributions remain the primary GEB evidence. Explicitly state why a direct claim of increased population stability, competition or demographic benefit is not justified. Do not change `main` or PR #38 unless a deliberate manuscript revision is approved.

### Reproducible sources
- [Focal analysis protocol](../protocols/BUTTERFLY_GEOGRAPHIC_GAIN_FRAGILITY_20261024.md)
- [BCE posthoc comparison protocol](../protocols/BUTTERFLY_BCE_GAIN_SUPPORT_SOURCE_SENSITIVITY_20261024.md)
- [Primary machine-readable results](butterfly_geographic_gain_fragility_20261024.json)
- [BCE source-specific results](butterfly_bce_paired_gain_support_20261024.json)
- [239-species analysis workflow](https://github.com/zuizui0223/chocho/actions/runs/40079582355) — **SUCCESS**
- [BCE archive comparison workflow](https://github.com/zuizui0223/chocho/actions/runs/40079912815) — **SUCCESS**

Source interpretation boundary: data measure modeled butterfly×geographic-region *potential larval host resources*, not observed adults, larval resource productivity, historical colonizations, protected populations, or plant-management effects.


## Audited figure and interpretation

[View two-panel PNG](figures/butterfly_gain_host_support_20261024.png) · [View editable SVG](figures/butterfly_gain_host_support_20261024.svg). [Successful source-locked rendering workflow](https://github.com/zuizui0223/chocho/actions/runs/40080207205).

Panel A presents each of the 206 gain-positive butterflies without smoothing, with the host-catalogue adequacy labels exposed rather than hidden. Panel B residualizes both the resource-gain predictor and the single-host-support outcome against native resource count, resolved host count and host-family breadth. The adjusted coefficient and genus-bootstrap slope envelope were checked against the immutable numerical receipt before rendering. The wedge is a **slope-uncertainty display**, not a prediction interval for local butterfly populations or host availability.

The companion [candidate Results, Discussion and legend](../../manuscript/butterfly_gain_support_exploratory_insert_v0.1.md) is an isolated draft addendum. Publication claims should remain explicitly about reconstructed potential host-resource architecture.
