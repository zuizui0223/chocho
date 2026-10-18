# When does moving a larval host plant actually produce more butterflies?

**Biogeographic sign-reversal perspective — 18 October 2026.** This is an exploratory, falsifiable new view of the *existing* global host-resource reconstruction in chocho, NOT a demonstrated global population effect, not a new confirmatory statistical result and not a second GEB paper. Original submitted GEB PR #38 is frozen. Keep the discussion within draft PR #43 pending new independent life-stage observations.

## One ecological question, three very different kinds of monotonicity

Human relocation of a plant can change at least three fundamentally non-equivalent objects.

| Ecological object | Quantity | Is it monotone with adding known larval-host plants? | Supported here? |
|---|---|---|---|
| **Host-resource presence** | Contemporary butterfly envelope `G_C(S)=union_{h in S}C_h` | **Yes**, provided original resources are not removed. A set-union theorem; not biological discovery. | Original 239-butterfly descriptive database |
| **Geographic novelty** | `U(S)=|G_C(S)\G_N(S)|` | **Not necessarily**. Adding a newly recognized native host can remove previously "new" cells, while adding an introduced-range host can create other cells. | **Yes, catalogue-sensitive**: 81 BCE-augmented butterflies have original/augmented source Jaccard 0.5004. Direct independent host quality example below. |
| **Real butterfly recruitment** | `B=eggs × (host-specific and season-specific survival to recruit)` | **Not necessarily** if host switching reallocates oviposition or adults, changes parasites, or changes habitat. Without such diversion/feedback, independent positive additions cannot reduce recruitment. | Mathematically conditional and supported by individual published natural-history examples, **not tested as a general global pattern** in chocho. |

**The actual scientific problem is the gap between these three levels**, not the set-theoretic expansion of a host-map union. Existing butterfly occurrence validation (66/115 outside-native observed butterfly×region cells associated with added host geography, region-matched null median 49) supports some potential relevance of the plant maps, **not an estimate of breeding on those hosts or causal butterfly range expansion**.

## A new, independent example within the global butterfly panel

The 239-species source includes **Danaus plexippus**, with 40 recorded accepted HOSTS species, native-resource 181 WGSRPD3 units, contemporary 263, and **82 distinct regions added by introduced plants**. To test how sensitive this mapping is to plant quality rather than merely the number of known host taxa, an independent, published source was joined on strict accepted/input taxon binomials:

- Greenstein, Steele & Taylor (2022), *PLOS ONE*, DOI [10.1371/journal.pone.0269701](https://doi.org/10.1371/journal.pone.0269701); open S1 source DOI [10.1371/journal.pone.0269701.s001](https://doi.org/10.1371/journal.pone.0269701.s001), SHA-256 `bc19e33823a3e860443b0086af129e39f4218663233fdceeccb6b08a15b3af79`. The independently downloaded and parsed appendix contained exactly **127** butterfly–plant source class entries: 34 H1–H3 high, 42 L1–L3 low, 33 N non-host, 18 U unsubstantiated; source publication is CC BY.
- On the *fixed original* 40 HOSTS accepted species, **23** were classified H1–H3, **8** strict experimentally established high H3, **8** low L, **2** source-classified N, **2** U, and **5** had no strict binomial match to the published 127. The last five are NOT "non-hosts" or low fitness. Fuzzy/synonym matches were not invented or auto accepted.
- From the **82 original new resource regions**, **78** contain at least one introduced HOSTS plant independently assigned an H1–H3 high class; **37** contain an introduced H3-class plant. Both statements are **source-level evidence coverage**, not observed larval performance at every global region or a direct estimate of survival.
- When **both native and contemporary maps** are rebuilt using ONLY H1–H3 species, the reconstructed native envelope **shrinks 181→152** and contemporary **263→253**, yet `U` **increases 82→101**. This counterintuitive **+19** follows because the native baseline shrinks proportionally more; **23 of the 101 high-filtered new cells were already native-resource-supported under full HOSTS** and are newly labelled 'added' by the narrower quality-defined native baseline. With only H3 species, native **120**, contemporary **188**, U **68**, but only **37 of the original 82** show introduced H3-class supports. **Different denominators and source-quality thresholds must not be compared as if they were one estimate of butterfly expansion.**
- The original 40-HOSTS source includes two names classed N in Greenstein's compilation. This is an evidence discordance between interaction sources, *not proof that the 40-plant HOSTS graph is wrong globally* or that actual larvae cannot utilize those plants under any conditions. It warrants taxonomic and regional evidence review before conservation claims.

The verified [quality-source workflow](https://github.com/zuizui0223/chocho/actions/runs/39096026037) fixed both external publication source and original HOSTS/WCVP source ZIP by byte SHA-256, reproduced 40/181/263/82 exactly, and checked publication class counts before evaluating geography. The [frozen result receipt](../../provenance/reviewer_defenses/results/butterfly_monarch_independent_host_quality_resource_novelty_v0.1.json) gives all class-specific results and names requiring manual review. This is a **genuine independent source-data sensitivity for one species**, not a new global fitness or population-size estimate.

## Under what conditions can added hosts lower production?

The source-conditioned, explicit model is in [`host_realization_nonmonotonicity.py`](../../src/butterfly_specialization_ecology/host_realization_nonmonotonicity.py), covered by [deterministic tests](../../tests/test_host_realization_nonmonotonicity.py). For seasons `t`, with a fixed egg budget `E_t`, `q_t` fraction of eggs diverted toward an introduced host, and per-egg recruit yields `w_N(t)` and `w_I(t)`:

`ΔB = sum_t E_t q_t [w_I(t)-w_N(t)]`.

If `q_t=0`, original native-host recruitment is unchanged; a bad new host that attracts **no eggs** cannot *by itself* lower total recruitment. When `q_t>0` and `w_I(t) < w_N(t)`, it can. The model labels its examples **illustrative-only**: in early season `E=100,q=.6,w_N=.5,w_I=.8` gives `ΔB=+18`, while in late season `E=100,q=.6,w_N=.7,w_I=.2` gives `ΔB=−30` (annual `−12`). No actual monarch population is asserted to have those parameter values.

**New source–sink pathway:** Even where the destination is a genuinely new U cell with no recognized native host, `B_new>0` can coexist with `ΔB_total<0` **if a fraction of breeders would otherwise breed in a more productive native-source population and become attracted to the new habitat**. A local gain is not automatically a metapopulation gain; explicit movement, survivorship and breeding site reallocation must be shown, not inferred from adult GBIF occupancy.

**Long-term pathway:** Disturbance-dependent introduced plants can give some unusually productive years and some failure years. Arithmetic-mean production and long-run geometric-mean growth need not share the sign of their differences; examples in the code establish this mathematically but are **not measured temporal fluctuations**. This links the benefits of a new resource to the conservation value of a reliable native alternative and phenological buffering.

Already-known empirical precedents, not unclaimed discoveries:
- Singer & Parmesan (2018), *Nature*, DOI [10.1038/s41586-018-0074-6](https://doi.org/10.1038/s41586-018-0074-6), documented an *Euphydryas editha* population that initially benefited from exotic `Plantago` but later lost native-host reliance and crashed after land-use change.
- Thiemann & Cipollini (2025), *Insects*, DOI [10.3390/insects16040331](https://doi.org/10.3390/insects16040331), documented late-season oviposition preference for invasive `Alliaria` paired with larval failure in `Anthocharis midea`.
- Yoon & Read (2016), *Oecologia*, DOI [10.1007/s00442-016-3560-2](https://doi.org/10.1007/s00442-016-3560-2), reviewed **76 studies** and already demonstrated frequent juvenile-performance reduction on exotics.
- Haan et al. (2021), *Scientific Reports*, DOI [10.1038/s41598-020-80413-y](https://doi.org/10.1038/s41598-020-80413-y), showed that an exotic host can perform well and ancestral host capacity need not be lost. **Do not equate exotic with universally bad.**

## Conservation decision rule: put *source function* before host presence

| Host map and field evidence | Main risk | What conservation should test before intervention |
|---|---|---|
| **A**: introduced host overlaps mapped native-host access, **and** local butterflies lay eggs on both | **Oviposition substitution or trap**, if the introduced option is attractive but lower-yield | Verify native host is accessible in the same patch; compare eggs per host AND independently followed larval cohorts and phenological windows |
| **R**: several introduced hosts support one mapped newly accessible butterfly region | **False security from structural redundancy** if all hosts disappear together or one diverts eggs from a better host | Test which plant species are actually successful local hosts and whether their leaf-quality windows are temporally complementary |
| **U**: introduced host apparently supplies the only new mapped regional access | **Anthropogenic dependence** or a source–sink draw from older populations | Confirm successful local development, demographic source contribution, adult movement and genuinely absent local native fallback before choosing invasive-host management |
| Apparent source-dependent novelty in **A/R/U** due to catalogue incompleteness | **Conservation misclassification** | Collect region-specific host evidence (eggs, larvae, plant taxonomy, successful development) and rerun the native-versus-introduced baseline |

Neither "introductions are good because they add butterflies" nor "remove all introduced hosts because they are traps" follows. The conservative management sequence is **assess existing native host quality, restore usable native alternatives, verify current dependence and stage-specific development, then manage introduced plants to protect both native communities and herbivore recruitment**. Some invasive-plant removal decisions may require transitional native-host provision; maintaining an invasive solely to help a native butterfly is not a general policy.

## Realistic pathway to a new high-impact macroecological study

The test that would go beyond the previous butterfly literature is **not** another simulation showing that traps are possible. It is a *replicated, response-blinded comparative empirical test* across butterfly species and sites: whether the sign of **juvenile production** in the presence of introduced hosts changes predictably with (1) native host alternatives, (2) realized egg diversion, (3) independent host survival classes, (4) host leaf seasonality, and (5) source–sink connectivity, against a host-presence-only benchmark. New scientific information requires multiple independently sampled butterfly systems and cohort/effort-aware outcomes.

Available primary inputs with different levels of identifiability:
- **Global mapping:** chocho HOSTS/WCVP 239 butterflies, 365 regions; fixed reconstruction, but **no observed reproduction**.
- **Independent performance contrast:** Greenstein et al. 2022 PLOS appendix, **127 candidates for one butterfly**; 35 of 40 original HOSTS host species have strict name matches. This produces a defensible single-species *quality-conditioned mapped resource* sensitivity, not across-butterfly demography.
- **Independent stage monitoring:** 2026-09-02 CC0 published IMMP California data DOI [10.5061/dryad.mkkwh718x](https://doi.org/10.5061/dryad.mkkwh718x) list egg/larva survey records in an `Act2` sheet. However its anonymous file-download API returned HTTP **401**; [machine-readable eligibility gate](https://github.com/zuizui0223/chocho/actions/runs/39095277102) records **NOT_EVALUABLE_DRYAD_AUTHENTICATION_REQUIRED**, not a negative result. Moreover same-day stage ratios would not equal cohort survival even after legitimate access.
- **Haan et al. 2021 field experiment:** published separate Dryad DOI [10.5061/dryad.612jm642h](https://doi.org/10.5061/dryad.612jm642h), with evidence about the quality of different ancestral/new host options, but one *Euphydryas* system only.

**Publication decision:** This is promising as a *next ecological question* and already sharpens what the current single GEB paper must not claim; it has not yet passed the minimum independent multi-species reproductive-data gate for a new globally supported mechanism. Do NOT add modeled negative butterfly fitness to the Results of the frozen submission PR #38, imply that 78/82 regions support actual successful reproduction, or make an automatic second paper. The source-specific monarch quality result may be valuable as a clearly limited supplementary exploratory bridge within PR #43 if the author approves.
