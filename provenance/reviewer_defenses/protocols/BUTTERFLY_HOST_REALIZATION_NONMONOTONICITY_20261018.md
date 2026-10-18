# Host redistribution is monotonically available but not necessarily reproductively beneficial (2026-10-18)

**Status:** prospective *new ecological hypothesis and feasibility gate*, developed after the HOSTS/WCVP paper's descriptive, original-HOSTS spatial, and BCE/Clarke cataloguing sensitivity results. These are NOT new empirical butterfly-fitness results, a preregistered confirmatory study, or permission to replace frozen GEB submission PR #38. Keep this as an independent conceptual development within *one* alternative butterfly resource manuscript (draft PR #43) until independent performance evidence is acquired.

## Core biological question

**When does a newly distributed larval-host plant *produce* butterflies, rather than merely become another mapped potential host or attract eggs away from higher-quality resources?** Human plant redistribution alters the *location of apparent interaction opportunity*, but local egg-laying, juvenile survival, seasonal resource availability and demographic source/sink dynamics govern whether it increases **real recruitment**. The latter effect **need not be monotone** in the amount of host occurrence or the number of host taxa.

This work is NOT a claim that all introduced plants are ecological traps, nor that the potential-range source Jaccard of 0.5004 measures actual recruitment. It must explicitly challenge the monotonicity of *effective* consumer support, rather than replay the elementary result that the union of planted-host regions expands as hosts are added.

## Structural monotonicity and its precise ecological break

Let `G_b(S) = union_{h in S} C_h` be the regional geographic envelope of a butterfly's recorded host plant set and distributions. If additional hosts are included without removing any original host range, `G_b(S) subset G_b(S union S_extra)` is guaranteed by set theory. The monotone union map is **not new biology**.

For an already natively host-supported local region and a fixed seasonal egg budget `E_t >= 0`, define:
- `q_t in [0,1]`: the **fraction of eggs actually diverted** from native to introduced host plants. This is a realized allocation fraction, not independent host-plant choice probability in every region.
- `w_N(t), w_I(t) in [0,1]`: expected recruit survival per egg, potentially including leaf phenology, local microclimate, parasites, predators and development (must be estimated, not assumed).
- `B_N = sum_t E_t*w_N(t)`, `B_mix = sum_t E_t[(1-q_t)w_N(t)+q_t*w_I(t)]`.

Then `B_mix - B_N = sum_t E_t*q_t*[w_I(t)-w_N(t)]`. This identity yields three **conditional diagnostic predictions**:
1. **No diversion, feedback or native-host removal**: adding an independent non-negative reproductive opportunity CANNOT reduce the pre-existing recruitment (`q_t=0`; additional recruits >=0). Poor-quality host presence alone is insufficient for an immediate negative effect.
2. **Substitution trap condition:** if `q_t>0` in active oviposition windows and the egg-to-recruit yield of the introduced host is lower than the native alternative sufficiently often, total recruitment can *fall* even though mapped resource coverage never falls. Host counts alone cannot establish diversion or yields.
3. **Seasonal role reversal:** `w_I(t)-w_N(t)` and `q_t` can change sign/magnitude across seasons. A host may contribute positively during one period and negatively during another. A temporally resolved data design must retain adult flight periods, host tissue phenology and egg-to-late-instar stage transitions instead of collapsing them into one annual regional occurrence.

For genuinely *new* host-supported region cells with no recorded native host (`G_N=0`), a newly accessible plant can add non-negative local recruitment **if actually used**, but if the plant attracts breeding females away from pre-existing demographic sources, whole-landscape recruitment can still decline. Write `B_total=(1-p)B_source+p*(1-m)*B_new`, where `p` is the redistributed breeder fraction and `m` movement mortality. Increasing `p` depresses output when `(1-m)*B_new < B_source`, even if `B_new>0`. This is a habitat/source–sink decision, not proof from current butterfly occurrences.

With year-to-year variability, even a non-declining arithmetic mean of annual growth need not imply non-declining long-term `exp(mean(log λ_year))`. Temporal failure, density dependence, parasite amplification and year-round foliage can change risk as well as mean reproduction. Quantify this **only with genuine year/site replication**, never from a self-chosen illustrative simulation.

## Tie precisely to the existing source-conditional macro-resource accounting

The original HOSTS-only 239-species reconstruction partitions 75,346 candidate host×butterfly×WGSRPD3 incidences:
- **A = 46,387**: introduced hosts in regions *already natively supported* in a coarse union envelope. They are candidates for substitution/oviposition redistribution questions, NOT evidence that females actually switch or that the native host is accessible in the same field patch.
- **R = 14,406**: extra introduced host incidences inside already newly supported regions. These are **structural redundancy only**. Whether they buffer failure depends on independent quality and asynchronous availability; a less suitable attractive host can *reduce* actual quality even if plant count grows.
- **U = 14,553**: uniquely new butterfly×region host opportunity cells. These are candidate geography, NOT source populations. Demographic sink risk is conceivable if migrating adults settle there and recruitment is low; local breeding data are required.

After European BCE/Clarke host inventory supplementation, A/H changes 61.6%→65.5%, U/H 19.3%→14.5%, and 81 species' region-identity Jaccard=0.5004. These changes reflect **catalogue completeness**, not field fitness. None can be used as a numerical proxy for `q_t`, `w_I/w_N`, or `λ`.

Thus the bridge from existing work is the **joint distribution of structural opportunity class, host quality, choice and phenology**. We do not have a global test of these demographic variables today.

## Independence and novelty audit: what previous literature already established

- Singer & Parmesan (2018), *Nature*, DOI 10.1038/s41586-018-0074-6: `Euphydryas editha` adopted non-native `Plantago lanceolata`, initially with superior survival; continued preferential host use created a management-sensitive dependency, and grazing withdrawal preceded local extinction. This **already demonstrates positive early fitness and later eco-evolutionary trap**.
- Yoon & Read (2016), *Oecologia*, DOI 10.1007/s00442-016-3560-2: meta-analysis of **76 studies** reported often lower juvenile success on exotic hosts and reduced lepidopteran assemblage abundance/richness in invaded vegetation. A general assertion that alien hosts can harm Lepidoptera is **not novel**.
- Thiemann & Cipollini (2025), *Insects* 16:331, DOI 10.3390/insects16040331: native `Anthocharis midea` showed strong late-season egg-laying on invasive `Alliaria petiolata` while the tested larvae failed to complete development, a **season-dependent preference–performance trap**. Seasonal mismatch as one counterexample is **not novel**.
- Greenstein, Steele & Taylor (2022), *PLOS ONE*, DOI 10.1371/journal.pone.0269701; [CC0 Dryad 10.5061/dryad.7wm37pvvp](https://doi.org/10.5061/dryad.7wm37pvvp): 127 historically reported *Danaus plexippus* candidate host plants were classified **34 high-performance, 42 low-performance, 33 non-host, 18 unsubstantiated**. Host catalogues do not measure functionally equal resources; this is a pre-existing, one-butterfly evidence set.
- Haan, Bowers & Bakker (2021), *Scientific Reports*, DOI 10.1038/s41598-020-80413-y: non-native *Plantago* sometimes supported high performance, but ancestral host use remained feasible. Do **not** simplify native=good/introduced=bad.
- Recent [CC0 Dryad IMMP California 2021–2022 site dataset](https://doi.org/10.5061/dryad.mkkwh718x), published September 2026: open milkweed abundance, repeated visits and egg/larva stages. Its 2026-09-02 Version 6 contains site/visit CSV and an `Act2` worksheet in the study XLSX; host-quality response variables, sampling denominators, and native versus non-native host comparability require an **observed-data schema audit before an empirical claim**. Other national IMMP site data require request and must not be assumed freely available.

**Novel angle worth pursuing** is a *comparative sign-change boundary* across butterflies and regions: given fixed mapped host abundance, does seasonal attraction to introduced versus native hosts predict a **reversal from more structural opportunity to less late-stage juvenile production**, after real sampling effort and habitat differences are accounted for? Its independent test remains future work. Numerical toy examples merely establish the model's logical possibility.

## A direct bridge inside the existing 239-butterfly panel: monarch host-quality sensitivity

The original, frozen `data/frozen/figure_sources/anthropogenic_species_metrics.csv` explicitly contains **`Danaus plexippus`**, with **40** resolved recorded host species, **181** native-supported butterfly×region resource cells, **263** contemporary cells and **82** distinct introduced-added cells. Thus independent life-stage evidence can challenge an **actual existing chocho macro-resource estimate**, not a hypothetical new example.

Greenstein et al. (2022), [PLOS article and supplementary Appendix S1](https://doi.org/10.1371/journal.pone.0269701), DOI [Dryad 10.5061/dryad.7wm37pvvp](https://doi.org/10.5061/dryad.7wm37pvvp), report **127** candidate monarch host plants classified as **34 high-performance / 42 low-performance / 33 non-host / 18 unsubstantiated**. A stronger, independently sourced stress test would reconcile the 40 recorded HOSTS hosts to Greenstein's name-resolved performance classes **before geographic outcomes are inspected**, then reconstruct *native and contemporary* envelopes for (i) all documented hosts, (ii) only high-performance hosts, and (iii) all high + low excluding confirmed non-hosts. Match accepted WCVP IDs and expose unmatched taxa rather than misclassifying or treating no evidence as low performance.

**Target estimand:** how many of the original **82** mapped introduced-added `Danaus plexippus` cells remain **new** with at least one high-performance host, and how many of 181 native cells remain after the quality filter? Use unique butterfly×region units, not the 127-candidate-classification denominator. This is **quality-filtered potential resource geography**, *not observed recruitment*, since per-plant performance measured elsewhere need not transfer across global climate, infection, season or native/introduced range. Conditioning the *native* baseline on host quality may both remove and add apparent new cells: report sets and source turnover, not only new-cell count. A meaningful result could be that source-defined **resource gain is not robust to quality**, but is not established until the crosswalk and analysis complete.

Potential data sources: the published PLOS [CC BY supplement](https://doi.org/10.1371/journal.pone.0269701.s001) and Dryad CC0 underlying classification; [the public Zenodo record](https://zenodo.org/records/6761802) describes the same study, **not an independent replication**. Avoid treating the 2022 paper's sample as global butterfly-level fitness; only `D. plexippus` can be evaluated in this specific crosswalk.

For later macro generality, the scope gate requires independently measured quality on **multiple butterfly species**, otherwise the manuscript can present monarch only as a prespecified mechanistic case study, not as a second global population-effect estimator.

**Access status as of 18 October 2026:** Dryad v6 2026 monarch IMMP California `Act2` workbook public metadata are available, but the anonymous API file endpoint returned **HTTP 401**, consistent with Dryad's documented authentication requirement for file downloads. The [source-only schema gate](../../../../scripts/audit_external_immp_lifestage_schema.py) correctly records **NOT_EVALUABLE_DRYAD_AUTHENTICATION_REQUIRED** rather than reporting juvenile survival, ecological trap evidence or a result-driven data selection. No credentials or unauthorised bypass have been attempted.

## Analysis and disconfirmation design: don't search for desirable significance

### Phase 0: public-source eligibility gate (before any result)

For each proposed public dataset, audit actual rows, species IDs, dates/site-year/visit denominators, native/introduced host identities, egg and juvenile instar counts, adult counts/effort, missingness, duplication, sampling protocol and licence. Record **PASS for specified estimand** or **NOT IDENTIFIABLE**, including no valid native/exotic paired comparison. Both the 2026 IMMP data and 2022 monarch host-performance dataset are **single-species systems**; they cannot on their own demonstrate a 239-butterfly macroecological law.

### Phase 1: independent plant/larval outcomes

Within matched site × season × butterfly species × plant taxon, estimate:
- `egg exposure`: eggs per monitored host plant × verified effort and host cover, not adult GBIF record count.
- `juvenile stage index`: late instar counts relative to lagged earlier eggs with interval and revisit modeling, after accounting for detectability and cohort dynamics. **A raw larva/egg ratio on the same date is NOT survival**, because eggs and larvae belong to overlapping, differently detected cohorts.
- `survival`: only where tagged egg clutches or explicit rearing/survivorship records follow the same cohort to a specified endpoint.
- `functional quality`: experimental larval performance and biologically matched leaf phenology; adult oviposition preference is distinct from actual eggs laid and leaf availability.

Standardize for effort, plant cover, geographic patch, month, weather, year and observer. Use proper overdispersed/hierarchical/zero-inflated models where identified, with held-out site-years and measurement uncertainty. Without independent measures of egg diversion or matched cohorts, do NOT claim a causal trap.

### Phase 2: cross-source linked macro hypothesis

Link only geographically validated host performance records to the 239-butterfly HOSTS/BCE network at compatible resolution. Classify site-level interactions as potential `A`, `R`, or `U` **only if local native alternatives are independently documented**; a WGSRPD3 source-only label does not prove patch-level alternatives. Fit interaction terms for host performance × season × native alternative access × introduced-host abundance; compare independent, response-held-out site-years against a monotonic plant-presence-only benchmark. A genuine extension would estimate the **fraction of independent butterfly–site–season observations with opposite direction of mapped opportunity and juvenile recruitment**, complete with sampling uncertainty and range restriction. If data are insufficient, report this non-identifiable, rather than synthesizing 239 imaginary reproduction outcomes.

### Phase 3: conservation counterfactual, conditional on new evidence

Compare **native host restoration**, **conditional exotic-host management**, and **phenological buffering** under observed site-specific recruitment, switching and seasonal use. A good conservation strategy retains or restores reliably productive *native* larval resources and tests native-host accessibility before exotic removal. Do **not** suggest preserving invasive plants automatically merely because a butterfly uses them; do not suggest indiscriminate removal where consumers have management-sensitive dependence without a local fallback assessment. *Functional rescue* is demonstrated by survival/reproduction, not occupancy or presence.

## Hard negative controls / stopping rules

1. Reproduce abundance-only predictions; if the new model improves neither held-out juvenile outcomes nor verified independent source occupancy, it is **not** a published new principle.
2. Fit native and exotic host performance after controlling plant species, site, season and effort. If no within-site/season native/exotic contrast, confounding prevents a trap conclusion.
3. Adult occurrence signals alone can validate potential geography but cannot estimate preference, recruitment or selection.
4. Do not treat all 34 high and 42 low *monarch* host categories as global butterfly fitness data.
5. Do not overwrite frozen GEB PR #38 or call a deterministic mathematical illustration an observed empirical discovery.
