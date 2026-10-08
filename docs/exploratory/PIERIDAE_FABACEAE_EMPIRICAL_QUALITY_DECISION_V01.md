# Pieridae–Fabaceae: from potential resources to realized ecological quality

**Date:** 2026-10-08. **Status:** Exploratory evidence audit, not a new biological result. Work confined to `analysis/resource-homogenization-v01`.

## Question

Do introduced Fabaceae that extend the mapped larval host-resource geography of Pieridae also yield viable butterflies, and are they safe from parasitoids?

## Observational evidence gate (predefined, failed)

The frozen 24-pierid GloBI pilot found 58 geographically mapped, dated and explicitly larval feeding-relation records. Twenty-two match an exact known accepted Fabaceae HOSTS link. Seven records (five butterflies, four WGSRPD3 regions) mention a known Fabaceae host classified as introduced in the recording region; **none** occurs in an independently validated introduced-only resource cell. The predefined minimum of 10 records across three butterfly species and three regions therefore **failed**.

Four of the seven observations had directly retrievable iNaturalist source IDs; the source audit found photo and larval annotation for all four, but only one had an explicit Eating field. One record marked the plant as cultivated. The other three records were literature-indexed observations and lacked direct occurrence IDs. Do not reinterpret missing data as true local absence, or annotation as photographic verification of ingestion.

Sources: [GloBI pilot run 37715427511](https://github.com/zuizui0223/chocho/actions/runs/37715427511) and [original-source audit run 37715634355](https://github.com/zuizui0223/chocho/actions/runs/37715634355).

## Independent published rearing evidence

Koptur S, Salas Primoli A, Ferreira Paulino-Neto H, Whitfield J (2024). *Pierid Butterflies, Legume Hostplants, and Parasitoids in Urban Areas of Southern Florida*. Insects 15:123. DOI [10.3390/insects15020123](https://doi.org/10.3390/insects15020123).

The study monitored two native (`Senna chapmanii`, `S. ligustrina`) and two non-native (`S. polyphylla`, `S. surattensis`) plants in three urban Miami locations. It monitored immature `Phoebis` butterflies, rearing them to either adult emergence or parasitoid emergence. The authors' printed all-site summary reports 17 parasitoid emergences among 235 native-host resolved outcomes and 49/283 among exotic-host outcomes. **Important: the native values are internally inconsistent:** the three site totals are 3+6+7=16 (not 17). In the `S. ligustrina` row, the site parasite counts 2+0+7=9 also differ from that row's printed total of 7. The source discrepancy is retained rather than silently fixed. The site-stratified Mantel–Haenszel OR uses the published site cells (16/235 native versus 49/283 introduced), whereas the published all-sites headline uses 17/235 native. The reported host-origin contrast goes *against* a simple enemy-escape hypothesis.

A reproducible reanalysis of the four-species × three-site Table 1 is in `scripts/audit_pieridae_fabaceae_empirical_quality_anchor.py`. It uses a site-stratified Mantel–Haenszel odds ratio, not a causal model. The paper itself documents different adult-emergence preferences between the two `Phoebis` species and includes parasitoid taxa on both native and introduced hosts. This establishes that **local larval success and enemy exposure vary by actual host identity**, but these differences are already published.

**Essential limitations:** Native/exotic status is confounded with which `Senna` species were sampled; repeated weekly counts on the same plants are not independent; 176 individuals found had an unresolved fate and are omitted from Counted; species-specific parasitoid attack risk is not available. The article's abstract has an inconsistent butterfly species epithet; its body/figures use `Phoebis sennae` and `P. philea`. Do not treat these numbers as evidence of regional population growth or new colonization.

## Exact HOSTS–WCVP link coverage check

The frozen HOSTS accepted sidecar omits `Phoebis philea × Senna surattensis`, `P. sennae × Senna surattensis`, and `Phoebis philea/sennae × S. ligustrina`, although the cited study documents usage/rearing on those plants in south Florida. The WCVP sidecar includes `S. surattensis` and `S. ligustrina`, with Florida introduced and native respectively.

A fixed-source sensitivity adds those four independently documented exact host links *without altering any plant native/contemporary map*. Critically, this is still a **global extrapolation of a Florida observation**, not a species-level global validation. `S. polyphylla`, another recorded exotic food plant, was not in the frozen accepted-host plant sidecar, and is not assigned an invented ID.

This is a positive biological link-source check, but it can change the resource-map sums in either direction because newly documented native host ranges affect the baseline as well as contemporary resource space. The whole 24-pierid panel's Fabaceae ecological benefit is not estimated.

## Scientific decision

1. **Structural geography:** within-family host-identity effects exist but are modest (+0.0105 Jaccard vs strict null, p≈0.03) and dependent on Fabaceae, Pieridae and a highly incomplete source-to-reproduction bridge.
2. **Observed larval use:** independent GloBI feasibility gate failed; handfuls of field records are not evidence of native-butterfly colonization in new resource-only regions.
3. **Fitness and enemies:** a published Florida rearing study demonstrates adults can be obtained from nonnative `Senna`, while parasitoid emergence among resolved outcomes is **higher** for exotic hosts in the same system. A potential host-resource gain need not equal a demographic gain.
4. **Publication:** do not promote a new ecological mechanism paper without an independently sourced, multi-species site-level success/abundance endpoint including control opportunities or a replicated experimentally observed fitness outcome; do not update main GEB manuscript.

The next substantive test must separate **host exploitation**, **larval-to-adult fitness**, and **enemy mediation** across introduced and native resources at matched sites. More static host-list randomizations will not accomplish this.

## A distinct anthropogenic route: cultivated host plants outside naturalization maps

The 2024 south-Florida study documents `Senna polyphylla` growing as an introduced ornamental host plant used by the monitored Phoebis butterflies, with adults successfully reared from larvae collected on these plants. In the study's Table 1, `S. polyphylla` accounts for **44 juveniles with resolved fates, including 8 parasitoid emergences** (18.2% among resolved outcomes); this is pooled over the Phoebis butterflies and sites.

Yet [Kew POWO's contemporary geographic distribution](https://powo.science.kew.org/taxon/234621-2) lists `S. polyphylla` as native to tropical America and the Caribbean and introduced in certain parts of South America/Africa, **not Florida**. It is absent from the accepted host-plant universe of the frozen HOSTS-WCVP sidecar. Therefore the observed **cultivated** host opportunity in Florida is invisible to the reconstructed `S. polyphylla` introduced-range resource map. Neither source proves naturalization in Florida.

This identifies **two distinct exposure-generating pathways**:
- **Naturalization pathway:** an introduced plant becomes established or otherwise listed in the regional botanical distribution, making it countable in WCVP-based reconstructed resources.
- **Horticultural pathway:** deliberate plantings of an alien host plant may support larval development despite no corresponding naturalization footprint in the WCVP regional distribution.

**Crucial limitation:** A single published cultivation case is an existence proof of a *measurement gap*, not evidence that urban ornamental resources generally rescue butterflies, increase reproduction, or systematically broaden distributions. Such claims would require replicated site-level plant presence and effort, feeding, adult emergence, and parasitoid data. Do not treat Kew's absence of an introduced-status listing as evidence that the plant is biologically absent from local gardens.

This is potentially a more interesting mechanistic question than additional reshuffling of static host lists, because an identical host-species identity can yield different spatial ecological opportunities depending on whether populations of that plant are cultivated, naturalized, or absent.
