# Cultivation vs plant origin in wild monarch larval photo records: information gate readout

**Date: 2026-10-28. Exploratory source-backed conclusion: NOT IDENTIFIABLE under the pre-fixed geographically replicated test.** This negative information-gate decision is not a hypothesis about biologically absent winter larvae. Original GEB submission and main unchanged.

## Original independent source and identification

- Original author: Micah Freedman (2026), preprint DOI https://doi.org/10.64898/2026.04.06.716508, public annotated host-plant dataset linked from https://www.freedman-lab.ca/resources. Original results already established that **209/348** winter-season grid cells had nonnative-only observed host records; **do not claim this as our discovery**.
- Exact downloaded public source ZIP SHA256 `4e4def63669b34738bcbfb8ba32ccf6c61d02cc6992b8f1ebfa8bc350a3987db`; author 12,942 source thinned monarch larval records (100-km/week/host thinning). After the source-format, wild-rearing, geography and accepted host-label filters: **10,237** eligible photo records. Source quality: **{'common_source_eligible': 10237, 'source_cultivation_missing': 3, 'source_geography_ambiguous_or_obscured': 1123, 'source_reared_or_rearing_ambiguous': 161, 'source_rows': 12942, 'unresolved_original_host_taxon_label': 1418}**.
- First original run returned zero due to erroneous expectation that the source stored full binomials. A separate pre-outcome source-format audit showed host labels are bare *Asclepias* epithets (e.g. `syriaca`), with other genus-only labels such as `Cynanchum sp.` unresolved. Source-only name normalization fixed parsing **without changing original biological strata, filters, thresholds or testing direction**. The initial zero was software, not biology. [Independent non-outcome audit #40541666884](https://github.com/zuizui0223/chocho/actions/runs/40541666884) and [source normalization note](../SOURCE_LABEL_NORMALIZATION_NOTE_20261028.md).
- Original source experiment [GitHub Actions #40542008389](https://github.com/zuizui0223/chocho/actions/runs/40542008389) **SUCCESS**; 8 scientific/source synthetic tests passed (later strict taxon correction test added); original output JSON SHA256 `a24bd23ee1faff9ed26f7232ba8011ef61b899a4a85693b25a119120ab1accaf` and GitHub archived aggregate only artifact ID `12750193213`. No original sighting IDs, observer names, coordinates, locality records or photos appear in this repository. 

## Frozen primary within-species biological test and failure

Frozen in `provenance/archive/monarch/SEASONAL_CULTIVATION_CONDITIONAL_TEST_20261028.md` before opening the seasonal outcomes: within host **source-native** *Asclepias* species × receiving State × year, compare observed larval-photo winter Nov–Feb composition in `cultivated=yes` vs `cultivated=no`, using source `rearing=no` records only. **Not a sampled-plant denominator nor larval production**. Original pre-fixed information gate ≥20 informative strata, ≥10 hosts, ≥10 States, ≥15 winter cultivated and ≥15 winter uncultivated records, no State ≥50% of matched record counts. Whole-State resampling 999 draws must give ≥950 finite draws to report even a descriptive bootstrap interval.

| Pre-fixed diagnostic | Actual source result | Threshold / interpretation |
|---|---:|---|
| Informative host × State × year strata | **17** | ≥20 required; **FAIL** |
| Native *Asclepias* taxa in informative strata | **7** | ≥10 required; **FAIL** |
| Distinct States with informative matched records | **4** | ≥10 required; **FAIL** |
| Total original photo observations in informative strata | **200** | conditional source counts, not independent plants |
| Largest State's proportion of matched records | **96.5%** | <50% required; **FAIL** |
| Cultivated winter / other-season photos | **31 / 86** | individual photos, not host availability |
| Uncultivated winter / other-season photos | **5 / 78** | individual photos, not ecological absence |
| Posthoc descriptive Mantel–Haenszel OR | **5.313** | NOT a qualified causal or general population effect |
| Finite State-cluster bootstrap resamples | **941/999** | ≥950 required, **FAIL**, no interval |

**Scientifically correct verdict:** the apparently large OR is dominated by one State and sparse within-host winter controls, so the stated cultivation-mediated seasonal release is **not demonstrated**. The original study is itself nonrandom sampling of *already observed larvae*, so no survival, leaf duration, real host use probability, adapted physiology, population growth or seasonal range expansion is measured.

## Frozen sensitivity and specificity diagnostics

| Contrast | Informative strata | Host taxa | States | Largest State | Photo records | Descriptive OR | Source gate |
|---|---:|---:|---:|---:|---:|---:|---|
| Add `likely` vs `unlikely` labels | 41 | 10 | 10 | 73.1% | 568 | 2.041 | **FAIL** |
| Host × 100-km hex × State × year | 7 | 4 | 2 | 91.7% | 24 | 1.829 | **FAIL** |
| High-confidence hosts only | 10 | 5 | 2 | 98.5% | 136 | 5.260 | **FAIL** |
| Equal-length shoulder months Mar–Apr/Sep–Oct vs winter | 12 | 6 | 3 | 95.5% | 88 | 3.270 | **FAIL** |
| Nonnative host only: cultivated vs wild | 4 | 1 | 2 | 93.8% | 32 | 0.094 | **FAIL** |
| Wild only: nonnative vs native hosts, State × year | 3 | not applicable | 1 | 100.0% | 20 | 3.900 | **FAIL** |

Including less-certain cultivation labels yields a point OR of 2.04, but the **descriptive** State-bootstrap 95% interval [0.31, 4.65] includes 1 and still fails the 50% State-dominance gate. Tight 100-km matching provides only 24 matched photos. Do not optimize the original thresholds after observing this.

## Consequence and decision

A species-specific host-labelled larval dataset is a major methodological improvement over GBIF adult occurrences, yet **positive-only photographic observations cannot identify a cultivation-induced biological seasonal niche shift**, and the source lacks geographically independent matched wild/cultivated same-species winter controls. Stop this source as a **global mechanism test**; retain as a failed geographic- and information-gate result. A future primary macro hypothesis requires replicated native-host wild versus planted plots with effort/number-of-host-stems/green-leaf availability and wild larval development or overwintering outcomes before claiming interaction-mediated niche release. Original preprint authors' winter host geography stays acknowledged as prior art.
