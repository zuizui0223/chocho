# From potential butterfly host resources to observed caterpillars: independent field-data cross-check

**2026-10-27 — completed, source-locked external empirical validation in draft chocho PR #81.** This is exploratory cross-system macro ecology, not the original GEB result and not a claim of experimentally identified host choice. GEB main/manuscript is unchanged.

## 1. Two different ecological observation stages must not be conflated

The original chocho GEB paper reconstructs geography of known larval-host plants with the source accepted butterfly–plant pairs fixed. Original aggregate potential host-resource area increased **+54.9%** across the original focal butterfly panel. This is a host-resource **opportunity calculation**, not observed wild larval occupancy, abundance or fitness.

A **separate published independent field project**, Whitener & Hurlbert (2026), [*Ecosphere* 17:e70615](https://doi.org/10.1002/ecs2.70615), used Caterpillars Count! standardized June–July summer foliage branch surveys in eastern North America. Their original [public data and R source](https://github.com/hurlbertlab/caterpillars-on-plants) are read at **commit `2c5690ac7c6dba7eca506da886c49781c37617bd`**, with original public Table S2 Git blob `bc0d66a2fa2e44c4c1d959fa4700f4256cbb238e`. Source analyzes **all Lepidoptera caterpillars (butterflies AND moths)**, not identifiable original 239 butterfly host-specific detections. Plant native vs alien origin is the paper's North American **USDA classification**, not WCVP host origin in the previous global country model.

The authors already reported **8.1% of surveys on native plants**, **2.2% on alien plants** containing a caterpillar among approximately 68,392 field surveys. These are **already published results**, not independent discoveries by the new cross-reconstruction.

## 2. A source-verified direct connection to original chocho's 239 accepted butterfly host portfolios

[Original public Table S2×HOSTS source CI #40476370333](https://github.com/zuizui0223/chocho/actions/runs/40476370333): SUCCESS, six synthetic checks passed. Exact original author-only WCVP/HOSTS 11-member SHA checks plus external public original Git blob. Original 239 butterflies and **2,596** accepted original butterfly×plant source host pairs restored before cross-joining.

The Table S2 published filter contains 100 plant taxa with ≥40 independent published survey events each: 86 USDA native, 14 USDA alien. **One of the 86 native taxa is only the genus *Zanthoxylum***; no species-level WCVP host ID was invented. The **99** unambiguous original field species names produce:

- **30/99** exact accepted original-239 host-plant matches (27 native, 3 alien), 69 original field species without a source-exact original-239 link. Those 69 are **NOT** true biological non-hosts.
- Matched 30 have 11,621 published branch survey events, 1,070 caterpillar-positive branch events. **19/30** plants have exactly ONE original-239 butterfly consumer source link; source portfolio size is poorly resolved on this small set.
- Native/alien source descriptive published 100-taxon panel: 52,826 native branch surveys (4,565 caterpillar positive, **8.64%**), 1,916 alien surveys (45 positive, **2.35%**). This differs from the paper's larger 68k-survey analysis because Table S2 prespecifies ≥40 surveys per plant taxon and includes only 100 taxa. The counts are reconstructed from Table S2 percent values rounded to two decimals and passed 0.0051-percentage-point exact-rounding tests.

### Source host-portfolio count does not generalize to *all caterpillar* field detections

On the 30 exact 239-butterfly accepted host plant species, a grouped binomial source model with USDA native/alien origin and log(1+documented original239 butterfly consumers) yields **odds ratio 0.92 per 1 standard deviation of source log host count**, plant-GENUS bootstrap 499 draws 95% **0.77–1.08**; genus-weighted/descriptive Spearman correlation of host count vs observed caterpillar fraction is **−0.036**. Leaving out each **plant genus** when training, the complete source model has a **−0.202 log-score gain per 1,000 held-out branch surveys** versus USDA native/alien alone (negative = worse prediction).

Two important source-outcome alignment checks:
- Treating original239 unmatched plant names as **0 documented original links**, not zero true hosts, yields 99 field species, **OR 1.06** (genus 95% interval **0.92–1.21**), out-of-genus predictive gain **−0.163/1,000 surveys**.
- Using documented HOSTS source consumers of **all original Lepidoptera** (including moths and butterflies), more relevant taxonomically to the field response, over 99 survey species yields **OR 1.14**, genus 95% **1.00–1.37** (lower numerical bound 0.9988), held-out genus prediction still does **not improve** (**−0.075/1,000 surveys**). This is a source recording coverage comparator, not the 239-only butterfly hypothesis.

The first 30-plant result includes just THREE alien plant species; its statistical information gate passed exactly its predeclared minimum, but this is **weak origin/consumer coverage**. No evidence that the original 239 HOSTS host density specifically predicts larvae in this taxonomically broader field study.

## 3. Rebuilt original raw survey IDs and matched same sites AND years

A second independent public-source run, [raw field audit CI #40477558557](https://github.com/zuizui0223/chocho/actions/runs/40477558557), recovered actual raw field source:
- **274,530** long-format arthropod/report rows;
- **68,561** unique surveyed branch-event IDs after prespecified June–July (Julian 152–212), longitude ≥−100, no wet leaves, excluding Coweeta, source botanical-name/USDA origin validity and matching distinct ID quality. The original authors report **68,741** in their somewhat different field plant-name pipeline; the difference **180** is documented as an incomplete reproduction, not adjusted by hand to match the paper.
- **212** distinct sites and **14** field survey years; 65,111 native and 3,450 alien source surveys; 5,326 survey IDs with source caterpillar detections.
- For **126 distinct site×year contexts at 71 sites** with both native and alien branch surveys (≥5 each), caterpillars were recorded in **2,189 / 28,215 = 7.76%** native plant surveys and **76 / 3,371 = 2.25%** alien plant surveys. Descriptive site-year Mantel–Haenszel native-vs-alien **OR 2.21**. In tighter same-site/same-7day-week comparisons (567 site×week contexts at 72 sites, minimum ≥2 each), native 1,927/24,923 (7.73%) vs alien 70/3,225 (2.17%), **OR 2.25**. These are still observational plant-origin associations, not an experimental intervention nor original biological novelty: the 2026 field study already established a substantial native-vs-alien caterpillar difference.

## 4. Does original chocho host-portfolio richness improve prediction WITHIN the same survey environments?

The critical new bridge is a source-exact 239-butterfly HOSTS link to public raw field-site survey IDs. In the wider raw dataset, **66** distinct surveyed plant names link exactly to an original 239 WCVP accepted host (51 native and 15 alien) at least once. Among the **38 site×year matched contexts at 22 field sites** with both native and alien linked plants:
- **37** distinct field plant species: **28 native and 9 alien**; 20 plant genera.
- 2,033 original branch surveys: 1,209 native (88 caterpillar detections, **7.28%**), 824 alien (34 positive, **4.13%**).
- Fixed site×year binomial-logit source model (38 stratum intercepts, original plant native/alien and log1p original239 documented butterfly consumer count) gives **original239 host count OR 1.128 per one original included-plant SD**, native-vs-alien conditional **OR 2.525**. Design condition number 54.0 and full rank are acceptable for this *source index model*.
- **Posthoc** uncertainty check after the initial source model results: resample entire original field sites (199 bootstrap replicates); host count source OR central percentile 95% **0.530–1.510**, encompassing no association; native-plant conditional OR **1.377–3.893**. Bootstrapping sites respects some survey dependence but does not correct plant species, phenology or botanical habit confounding and is NOT a natural experiment.
- Crucial precommitted **leave-one-whole-site-out predictive check** across all 22 sites: source HOSTS consumer count offers **−0.0873 log-score units per 1,000 held-out branch surveys** beyond USDA plant native status (negative = slightly worse, very small magnitude). Thus the imported global HOSTS host consumer richness has no demonstrated incremental **out-of-site** predictive skill in this observation setup.

[Final all-source site-year workflow CI #40478643768](https://github.com/zuizui0223/chocho/actions/runs/40478643768): all 13 offline source-independent tests PASS; 68,561 unique original field ID, original author-only 239 accepted host graph fixed, site-year gate and 199 full site-block resamples completed.

## 5. Conclusion and scientifically correct scope

**Supported empirical conclusion:** geographic reconstruction of butterfly potential larval host resources and realized field caterpillar detections are separate empirical axes. The published independent field data show a strong plant-origin-associated difference in wild caterpillar *encounters* that persists in site/year-matched subsets. The NEW direct cross-source analysis shows global documented butterfly consumer source richness **does not improve predictions of wild all-Lepidoptera larval detection** after source plant origin when evaluating unseen plant genera or sites.

**Not supported:** evidence that original 239 butterflies fail to use introduced hosts, that plant introductions reduce original butterfly population abundance, that host switching has a global rate, that introduced plants specifically generate interspecific competition, that caterpillars complete their life cycle on matched host plants, or that butterfly climate niches expanded. Original source mismatch: 239 BUTTERFLIES are potential-consumer predictors, whereas outcome includes all moth/butterfly larvae; 99 source species in the published table, only 30 exact 239 matches, and raw matched 37 species (9 alien). The source HOSTS catalogue is reporting biased, "no accepted HOSTS link" is NOT biological absence, and native woody plant taxonomy and survey season strongly affect detection.

**Research bridge:** original chocho potential host geography (`resource opportunity`) → independent observed branch-level Lepidoptera larvae (`field encounter`) → independent butterfly larval survival/fitness (`recruitment filter`, still not measured in this macro dataset) → potential demographic expansion/competition (still unmeasured).

A different published field/lab experiment by Ferreira & Rodrigues (2021), [*Ecology and Evolution* DOI 10.1002/ece3.7821](https://doi.org/10.1002/ece3.7821), demonstrates a potential survival/performance bottleneck in two *Danaus* butterfly taxa on exotic versus native milkweeds, BUT the native and exotic field hosts were sampled in **different sites 365 km apart**, so the apparent in-field plant-origin effect is confounded with environmental locality and cannot be used here as causal support for an alien host fitness penalty. It does show why mere botanical availability or larval observation is not a full-life-cycle demographic metric.

## Source files and executable evidence

- [Original source-blind cross-data protocol](../../provenance/reviewer_defenses/protocols/FIELD_CATERPILLAR_REALIZATION_FILTER_20261027.md) and [source matched site-year protocol](../../provenance/reviewer_defenses/protocols/FIELD_RAW_SITE_YEAR_CATERPILLAR_REALIZATION_20261027.md).
- [Original accepted HOSTS×public 100-plant source comparison](../../scripts/analyze_field_caterpillar_host_realization_filter_20261027.py), [original raw survey ID/site/year and USDA plant origin audit](../../scripts/audit_field_raw_site_year_caterpillar_20261027.py), [site-year fixed effect/out-of-site model](../../scripts/analyze_field_site_year_HOSTS_portfolio_predictiveness_20261027.py).
- [First source CI #40476370333](https://github.com/zuizui0223/chocho/actions/runs/40476370333), [final all source site-year CI #40478643768](https://github.com/zuizui0223/chocho/actions/runs/40478643768), numerical source-only artifacts.
- All raw original private WCVP/HOSTS interaction link identities, field survey IDs/site coordinates, and potentially identifying individual source records are withheld from public artifacts. Only field plant species-level aggregates, botanical name-match QC and source-field site/year aggregate figures were committed or archived.

No changes to original GEB main manuscript or figures. Draft stacked PR #81 remains an **exploratory source-linked ecology** continuation.
