# HOSTS report geography × WCVP accepted-plant origin: first source-verified macro result

**27 October 2026. Source-replayed, exploratory, not proof of new butterfly host adoption or invasion.** Independent branch stacked on source Location inventory PR #75. Original GEB v0.2 and main unchanged.

## Source and analysis status

- [Complete GitHub Actions source run #40439954723](https://github.com/zuizui0223/chocho/actions/runs/40439954723): **SUCCESS**, 7 source-independent synthetic tests passed, original 140,485 HOSTS source rows and 86,826 species-level HOSTS rows reproduced, exact 11-member author-only original WCVP/HOSTS archive checks passed.
- Original **239 focal butterflies**; 4,695 original HOSTS source rows matched insect species, 3,898 original raw rows joined by exact pair name to accepted WCVP host IDs; **3,586 unique butterfly×accepted-plant×raw-Location triples** (not independent populations or citation-level observations). The other 797 original focal insect rows lacked a source-fixed accepted pair.
- Original country–WGSRPD3 geography from SHA-verified 969 **disjoint area atoms** and 247 country map units, source commit `4104c784ddcf97318c4a1bb8bd949b3147628494`. Plant origin status requires ≥95% of **entire original country land area** to share native, introduced or WCVP-unlisted source status; country mapped coverage must be ≥95% and ISO geometry must be source-qualified. Original Hawaii→HAW and West Malaysia→MLY exact island/territory source codes separately checked; unrelated East Malaysia/Java/Ryukyu labels NOT guessed. Threshold ≥80% tested separately.

## 1. Global source evidence by REPORT Location, main 95%-country area threshold

| Exact accepted 239-butterfly × host × original Location unique triples | n |
|---|---:|
| Native-source plant region | **399** |
| Introduced-source plant region | **179** |
| Extant nondoubtful WCVP plant distribution **not listed** in source location | **303** |
| Mixed or partial WCVP botanical status within country | 521 |
| Original country ISO geography not source safe | 309 |
| Original subcountry location lacks unambiguous botanical region | 42 |
| Broad biogeographic realm, provisional, interception, unreconciled label | 1,833 |
| **Total, no silent exclusions** | **3,586** |

Thus 881/3586 = **24.6%** are classified directly as native / introduced / WCVP-unlisted, with the rest retained as explicitly unresolved. Excluding unlisted from the supported extant-plant-native/introduced categories leaves 578 source triples, **179/578 = 31.0%** reported against introduced plant range. This fraction is a **filtered reporting sample**, not a prevalence estimate for worldwide butterflies, all larval hosts, or invaded butterfly populations.

## 2. Apparent source-native versus introduced documentation is highly asymmetric

Introduced-zone host-use catalogue records comprise **179 distinct report-locality triples**, for **156 unique butterfly–accepted host plant pairs** spanning **62 butterfly species**. Only **12/179 introduction-locality triples** and **11/156 distinct pairs (7.1%)** have the *same exact accepted pair* documented at another original HOSTS Location that was independently confirmed as the plant's source-native region under the same 95% rule. The rest **167 triples and 145 pairs** have no source-native-Location record.

**NOT evidence of 145 proven new host switches.** These could have been used in native areas but never indexed to a precise report location. WCVP accepted host origin, original HOSTS report effort and the native ranges of the butterflies are three different constructs. In fact 62 source species have some introduced-area host report but only **45** have any source-native plant-area host report at all, highlighting source asymmetry.

## 3. Frozen geography sensitivity and heterogeneity

If the required fraction of a country's area dominated by a native/introduced/unlisted WCVP status is relaxed from **95% to 80%**, while the original source identities and geographic safety gates stay fixed:

- Native 399→424; introduced **179→197**; unlisted 303→475; mixed 521→306; all source-unsafe/unmapped numbers unchanged.
- Introduced-location triples sharing the same accepted butterfly-host link with a source-native-location report **12→14**; no documented native-area report **167→183**.
- This indicates that country size, partial coverage and WCVP range generalization matter; never select a more favorable cutoff after observing the answer.

Top families among the **156 introduced-area accepted species×host pairs** (descriptive, not enrichment against available or reported host families): Fabaceae 50, Poaceae 22, Brassicaceae 21. Those 93 pairs are 59.6% of introduced-area pairs; source host composition and crop status must be controlled before interpreting these as a crop/anthropogenic ecological effect.

Original high-volume Location reporting examples at 95% (WCVP origin matched to an exact accepted focal butterfly/host pair):

| Original Location | Introduced triples | Native triples | Unlisted triples | Main qualification |
|---|---:|---:|---:|---|
| Taiwan | 20 | 74 | 48 | Country unit source matched |
| West Malaysia | 29 | 17 | 31 | Exact verified MLY only, not all MY |
| Japan | 9 | 36 | 19 | Country area status |
| Brazil | 8 | 26 | 30 | 50 more mixed botanical status |
| USA | 14 | 6 | 70 | 369 mixed; US cannot be one homogeneous botanical unit |
| Hawaii | 12 | 4 | 8 | Exact HAW botanical region, not US-wide |
| India | 0 | 0 | 0 | 130 classified COUNTRY_MAPPING_UNSAFE due original disputed/missing ISO source geometry, not a plant-use absence |

The historical preliminary '776 introduced-area records, 49 matching source-native' used a different source population (18 selected geography labels and all potentially eligible insect species) and different geography rules. It should **not** be called replicated by this strict original-239 179/12 result, despite superficial numerical resemblance.

## 4. What this unlocks, without overclaiming

The prior GEB paper estimated geography of *potential known-host resources* under plant introductions: 206/239 butterflies expand reconstructed resource geography and aggregate area +54.9%. The new HOSTS Location × WCVP evidence demonstrates real **geographically reported host links within introduced plant areas**, but source-native reporting evidence for most of those links is missing. This supports a new testable distinction: **geographic continuation of established consumer–resource edges** versus **post-introduction source recording of host relationships without a native-area counterpart**. It does not establish that a butterfly adopted that host for the first time, that the insect is invasive in the recipient area, that larvae complete development, or that climate niche expanded.

For a biological test, link verified source-native **consumer** distribution and dated receiving establishment to each candidate host-use site, and compare host-record discovery relative to source sampling effort / source botanical geography. The currently validated invasive-butterfly historical roster from PR #75 has 5 original focal species and 13 initial regional events (3 additional butterfly species outside 239), but their original HOSTS Location evidence is too sparse for population-level host-adoption inference.

**Next statistical experiment:** separate the numerator (newly documented introduced-area host reports) from recording exposure (per-species HOSTS source intensity and mapped Location availability), native-area reference completeness and original host breadth. Do not promote a simple '145/156' to a real ecological-host-adoption rate.

## Reproducibility / data distribution boundaries

- [Exact source-locked analysis code](../../scripts/analyze_hosts_location_wcvp_origin_20261027.py), [synthetic source controls](../../tests/test_hosts_location_wcvp_origin_20261027.py), [frozen protocol](../../provenance/reviewer_defenses/protocols/HOSTS_LOCATION_WCVP_ORIGIN_20261027.md).
- [Successful complete source replay](https://github.com/zuizui0223/chocho/actions/runs/40439954723) with [aggregate-only JSON artifact](https://github.com/zuizui0223/chocho/actions/runs/40439954723/artifacts/12706166594).
- The original HOSTS/WCVP raw interaction names, exact accepted butterfly×host×location triples and original botanical polygons are intentionally not redistributed.
