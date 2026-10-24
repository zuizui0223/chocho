# Butterfly resource niches, actual abundance and competition — evidence hierarchy (2026-10-25)

## Biological question
When human plant redistribution increases the potential geography of shared larval hosts, does that produce (a) positive local co-occurrence or population growth, (b) heterospecific larval resource competition, or (c) coexistence through niche partitioning? **Each level demands different observational grain, different denominators and different evidence.** A rise in shared potential resource geography cannot identify demographic effects.

## Existing chocho biological-scope facts
- The 239 butterfly HOSTS/WCVP reconstruction adds 14,553 butterfly × region potential resources (+54.9%) across 206 butterflies; 670 host plants contribute, with 38 accounting for ~half of the added fractional credit.
- Potential exact-host shared butterfly pair × region instances increase 62,473 → 141,885 (+127.1%). This is a **set-union resource opportunity** statistic, not contemporary butterfly competition, increased recruitment or even site-level co-occurrence.
- Potential region-resource Jaccard increases 0.277 → 0.462, mostly explained by marginal addition; the excess over fixed-margin expected ~0.453 is small (~+0.0086 absolute Jaccard).
- The proposed macro relation "greater geographic gain means exceptional multiple-host demographic buffering" is **falsified at the structural level**: observed species-level β=−0.2207 sits inside two source-locked family/continent fixed-degree botanical nulls (β≈−0.222). Avoid reusing it as an ecological advantage.
- New individual abundance data from **Henry et al. (2025), *Scientific Data*, 10.1038/s41597-025-05513-8** acquired and authenticated: Figshare archive doi:10.6084/m9.figshare.27934602.v1, MD5 3161b1c8802a0d67771b2044e18da928; **1,321,882 original reported species×survey records from 19 programs, explicit zero records 425,802, 121/239 accepted butterfly binomials represented**, 2000–2020. These records are actual adults counted (not absolute estimated breeding populations), and analyses must not pool inconsistent survey types. [Source quality readout](../../provenance/reviewer_defenses/results/sbus_original_survey_count_source_gate_20261025.json) and [program effort/repeatability report](../../provenance/reviewer_defenses/results/sbus_original_program_repeatability_20261025.json).

## Niche/competition evidence already published (not claimed as independent chocho discovery)

### Regional host-resource filtering
[Nakadai et al. (2018), Oecologia](https://doi.org/10.1007/s00442-018-4062-1): 232 Japanese butterflies at 1 km and 10 km grid scales; sharing larval hosts predicts less geographical exclusivity, compatible with niche filtering rather than competition at the **distributional scale**. This does not rule out larval competition after species are filtered into a locality.

### Experimentally manipulated asymmetric butterfly interaction
[Hashimoto & Ohgushi (2023), Ecology and Evolution](https://doi.org/10.1002/ece3.10164): experimental density manipulation on a shared *Aristolochia debilis* resource, with invasive *Sericinus montela* and native *Atrophaneura alcinous*. Raising *S. montela* larval density reduced survival and prolonged development of *A. alcinous*; reverse interaction lacked detectable effect. *A. alcinous* larvae consumed 2.6x more plant mass in the lab, but the *S. montela* effect was not attributable simply to reduced remaining plant quantity, implying unmeasured plant quality or other indirect responses. **One legitimate case of asymmetric interaction; cannot multiply one case into a global universal competition coefficient.**

The repo's archived matched starting density experiment reanalysis (16 of original 30 cages) shows an early difference in focal aggregate counts before any recorded pupation, but individual death causes/leaf mass are not available; this is **prior-art calibration**, not new original ecological effect.

### Microhabitat niche partition on the same host
[Jang et al. (2026), *Niche partitioning of Sericinus montela and Atrophaneura alcinous in small scale riverside habitats*](https://snu.elsevierpure.com/en/publications/niche-partitioning-of-sericinus-montela-and-atrophaneura-alcinous/): Both larvae use *Aristolochia contorta* but in different riverbank/light microhabitats, with opposite light-intensity associations and evidence of different leaf resources/preferences. This illustrates why sharing a botanical *species* need not imply sharing **suitable leaf biomass × light × instar × time**.

### Density dependence is not necessarily interspecific competition
[2009 Maculinea study](https://pubmed.ncbi.nlm.nih.gov/19484269/): density dependence explained considerable annual butterfly variance, whereas potential heterospecific competition effects were weaker (study-specific). [2017 monarch density study](https://www.nature.com/articles/s41598-017-14510-w) similarly shows adult/egg density dependence does not tell which local mechanisms limit growth. Repeated adult transect series confound larval crowding with weather/flight-season detectability.

### Population dynamics need not be positive after resource globalization
[Edwards et al. (2025), *Science*](https://pubmed.ncbi.nlm.nih.gov/40048533/): 35 U.S. monitoring programs suggest a **22% butterfly abundance decline across the United States between 2000 and 2020**, despite widespread observed plant globalization elsewhere. This is **not** an experiment testing host introductions and not evidence that plant globalization caused the decline; it motivates avoiding extrapolation from potential maps.

## Proposed nontrivial ecological mechanism (untested)
**Opportunity–density transition with microhabitat-dependent competition:** an introduced host can initially relax environmental filtering and make multiple larvae's potential geographic resources co-available. Competition becomes important only when simultaneously active larvae's total demand approaches or exceeds the *effective* host biomass produced and suitable during their larval development window. Same nominal plant species may segregate in leaf age, sunlight/shade, host tissue, time and natural-enemy context. With low effective demand/production ratio, shared hosts can promote regional co-occurrence without competition; above a density/quality threshold, local heterospecific density effects can become asymmetric. Neither a map of resource union nor adult counts alone measures that threshold.

**Test design:** 
- Regional (resource filter): butterfly host overlap versus spatial co-detection controlling abiotic similarity and survey effort.
- Local (growth/interaction proxy): same-site consecutive-year butterfly abundance transitions and *previous year* host-sharing peer abundance, controlling own abundance, comparator non-host-sharing peer, effort, species–site and year, and testing degree-preserving butterfly-host graph nulls.
- Mechanistic (competition): follow eggs to adults while manipulating competitors and edible host leaves; measure host biomass, phenology, and microhabitat light. This is required to separate food exploitative competition from plant-quality-mediated indirect interactions.
- Temporal plant-release mechanism: require time-stamped *local* introduced-host emergence/expansion, not WCVP region-level contemporary naturalization flags.

**Main paper decision:** the current GEB `main` and draft PR #70 remain separated. Until a numerically qualified independent monitoring source panel passes a guarded ecological comparison, do **not** upgrade competitive exclusion or demographic increase to a result.

## Source and licence
Actual U.S. butterfly source raw dataset: Henry et al. (2025), doi:10.1038/s41597-025-05513-8 and Figshare doi:10.6084/m9.figshare.27934602.v1, CC BY 4.0, source archive MD5 pinned; only aggregate counts/test results are stored in `chocho`. U.K. UKBMS site-index source doi:10.5285/53a1c705-4437-48bb-9bcd-b68977f8edd9 did not provide direct public API bytes without programmatic EIDC access and was not used to substitute GBIF occurrence counts for butterfly population size.
