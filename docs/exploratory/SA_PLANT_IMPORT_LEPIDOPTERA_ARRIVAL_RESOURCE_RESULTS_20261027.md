# Arrival channels versus larval-host resources: first actual commodity-level Lepidoptera interception analysis

**2026-10-27 — independently source-verified South African public plant-import laboratory data; exploratory descriptive, not causal arrival→establishment inference.** Draft PR #83 on PR #82. Main GEB manuscript/figures unchanged.

## 1. Real original data acquired, reproducibility and sample-grain correction

The ORIGINAL public article is Saccaggi et al. (2021), *Scientific Data* 8:83, DOI **10.1038/s41597-021-00869-z**. Its original Figshare DOI is the **collection** **10.6084/m9.figshare.c.5180681**, with:
- Figshare article **13118177**, original `Saccaggi et al 1 list of contaminants on SA plant imports.csv` (file ID 25166393, SHA256 `c12d4af0c96a247d520f71ba6c2dfa67eacb33b53dd5aa329265c850d3fe15da`), **26,295 original raw CSV lines × 92 columns**.
- Figshare article **13668839**, original `Saccaggi et al Metadata of contaminants on SA plant imports.csv` (file ID 26235020, SHA256 `0f7ad1b9b44dad51a53d55a1f6bad193d623b0b8beb7853eb0b6f025dfcf640c`), **878 original taxon rows × 14 columns**. The exact Lepidoptera result field is `Insecta: Lepidoptera`.

[Full source-only acquisition with five independent tests](https://github.com/zuizui0223/chocho/actions/runs/40489503139) SUCCESS. Earlier Figshare ARTICLE ID 14393723 is a 2022 *Ecological Applications* follow-up dataset, **NOT** the source of the 2021 laboratory record matrix; it was rejected from this study rather than silently substituted. The two Dryad data sources DOI 10.5061/dryad.dv41ns243 and 10.5061/dryad.31zcrjdnw were independently identified, but their public ZIP API requests returned HTTP401 in GitHub Actions. No species×commodity data were invented from those Dryad files.

**Critical source denominator fix:** The raw Figshare samples sheet has 26,295 CSV rows. Exactly **4** rows have neither original `sample.source` nor `crop.commodity`: 3 have a nonzero literal in `Insecta: Lepidoptera` and 1 has no test label. These four cannot be treated as verified inspected consignments. The remaining **26,291 source-regime-tagged original sample rows** equal the paper's *initial* 26,291 lab records; its final cleaned publication panel uses **25,279**, so a further **1,012 eligibility/exclusion decisions** have not been reproduced. All commodity-specific denominators are explicitly **pre-paper-final source-regime-tagged**; no claim of rerunning the authors' exact final cleaned population.

Among the 26,291 source-regime-tagged sample rows:
- **15,848** have a Lepidoptera test status, **15,673** have literal `0` (tested negative), and **175** have nonzero Lepidoptera-positive taxonomic labels;
- **10,443** have `NA` (not tested), which must NEVER enter the denominator of a detection probability;
- the original four unidentified source rows account for the additional 3 apparent positives (raw all-line count 178), not source-verified shipment detections.

In other words, the source-eligible conditional detection frequency is **175/15,848 = 1.104%**, NOT 178/26,295 or 175/26,291. This is still only a **laboratory selected inspection fraction**, not a population-wide arrival rate.

## 2. Original source biological-material classes: source detection, not origin effect

The predeclared nonexclusive biological *material class* is constructed from the exact original 30 commodity labels (plus one missing source label). Original detail remains in the machine-readable output. One source-vocabulary correction assigned `Ex vitro plantlets` to LIVE_PROPAGATION along with `In vitro plantlets`, seedlings and cuttings; this changes live tested denominator by +57 but no positive cases.

| Original imported plant material, source-regime-tagged | Lepidoptera-positive / Lepidoptera-tested | Conditional tested-detection % |
|---|---:|---:|
| LIVE PROPAGATION: rooted plants, cuttings/budwood, bulbs, tubers, live plantlets | **41 / 10,337** | **0.397%** |
| SEED/GRAIN: seeds, beans, grain, nuts | **30 / 389** | **7.712%** |
| FRUIT/LEAF/FLORAL PRODUCE and dried goods | **104 / 5,039** | **2.064%** |
| OTHER/UNKNOWN source-regime-tagged | 0 / 83 | 0% |

**The visually dramatic seed/grain advantage is dominated by grain and NOT by the type of seed originally planted to make host forage.** Source original:
- **Grain:** **22/160** tested positive (**13.75%**); of 22 positive specimen labels, **6 explicitly mention larvae**, 16 do not specify stage.
- **Beans:** **4/121** (3.31%), 1 explicit larval label, 1 pupa, 2 unspecified.
- **Nuts:** 3/15 (20%; tiny sample, high uncertainty).
- **Seeds:** only **1/93** (**1.08%**), with life stage unspecified.

Thus a **grain shipment** (product potentially carrying moth pests) is not the same exposure as a **seed consignment that established a relocated forage host** in the chocho WCVP map. Seed-propagated plants can create postarrival larval food, but seed itself need not be the vector for insects that feed on mature leaves. Conversely a grain shipment can carry food/feed insect stages without producing a planting, distinct from a live ornamental plant that serves simultaneously as commodity and host.

Stage information in the **175** valid taxon-positive source samples: 86 source literal explicitly labelled larval, 4 explicitly pupal, the remainder had no definite stage annotation. Source split before excluding the four unknown rows: LIVE material 28 explicit larva/41 positive; FRUIT/LEAF etc 51 explicit larva +3 pupae/104; SEED/GRAIN 7 explicit larva +1 pupa/30. This is **not** rearing to adult or proof that each specimen was using the cargo as its actual host.

### Source sampling regime is overwhelmingly confounded with material

The original manuscript distinguishes audit shipments with no obvious contamination, intervention after a concern is raised, and PEQ propagation plant inspections. The source conditional rates illustrate their nonequivalence:

| Source regime | LIVE samples detected / tested | SEED+GRAIN samples detected / tested |
|---|---:|---:|
| **audit** | **11/544 = 2.02%** | **26/256 = 10.16%** |
| **intervention** | 3/166 = 1.81% | 4/112 = 3.57% |
| **PEQ** | 27/9,627 = 0.280% | 0/21 = 0% |

Audit (n=10,629 raw samples) is *ad hoc*, not randomized over imports. PEQ contributes almost all live-plant tested specimens, but few comparable seed specimens. This is a **source-selected laboratory detection comparison**, not a tested claim that living plants generally have fewer vector-transported larvae than grain. The original years span 1994–2019 with substantial changes in recording before/after ~2005, and the source lacks the total volume and size of imports by commodity.

## 3. For original chocho 239 butterflies, actual commodity/taxon-level co-detection exists but is rare

Only **two original-239 butterfly species** can be exactly identified from these source `Insecta: Lepidoptera` positive labels among 43 original Lepidoptera taxon metadata rows:

1. ***Lampides boeticus*** (Lycaenidae), **2018 audit**, on cargo identified as ***Phaseolus vulgaris***, material **Beans**, one positive sample. **That same accepted plant species is an original documented larval host in the original SHA-verified chocho HOSTS source.** This is actual import *butterfly taxon × listed original host plant commodity* co-detection, not proof that the intercepted insect was a larva or feeding on the beans. No established South African butterfly population or founding shipment is inferred.
2. ***Pieris brassicae*** (Pieridae), **1997 intervention**, on commodity **Fresh foliage**, botanical crop recorded ***Brassica sp.***, one positive sample. Genus only: cannot equate this with an exact accepted *Brassica* species larval-host match. Original life stage unspecified. Demonstrates butterflies can be intercepted on cut foliage without a whole rooted live plant, not proof of larval host use or establishment.

The original source has many moths/higher labels: in the all-raw-file 178 positive cells, **86 explicit source larval labels** plus four pupa labels and 88 unspecified or species labels. The 2 original chocho butterfly taxa have **no source-explicit immature stage**, so do not infer egg/larva attachment even where cargo is a documented larval host.

**Historical establishment warning counterexample:** Original taxon metadata classifies ***Spodoptera frugiperda*** (a MOTH, not one of the chocho butterflies) as historical South African absent and subsequently present **(2017)**. There are **three exact-species positive import samples** on ***Zea mays***: two **Plants** and one **Grain** in the **audit** sampling regime. All were from **2018–2019**, *after* the recorded establishment year. Thus the recorded shipments cannot be used as an empirically documented pre-establishment propagule arrival that caused that particular invasion. This is a powerful reason to keep source interception chronology and established-population chronology separate.

## 4. Which ecological two-channel mechanism survives?

**Source-supported:**
- insects can be detected on live plants, cuttings, fruit/foliage produce and seed/grain goods; these are **distinct carrier modes** with different likelihood of transporting an immature stage;
- the same known larval-host botanical species may appear as an intercepted commodity (*Lampides* × bean), but other consumers may arrive on unrelated plant products (hay hypothesis, some moths on stored grains);
- *plant resource redistribution* and *insect arrival* represent **two exposure processes** and need separate, species/destination/year provenance;
- this dataset carries a rare historical taxon known to later establish in SA, but its identified positive interception samples postdate establishment, so a bad species-wide join would falsely infer predicting invasion.

**Not supported/testable from these data:**
- no source-verified `butterfly × first international shipment × stage × establishment` complete panel beyond two butterfly-positive commodity co-detections;
- no comparable global original WCVP potential host plant opportunity map vs actual dated arrival *for each butterfly* in the same receiving country with effort controls;
- no confidence that seed-propagated host resource opportunities have a quantifiable 'orphan' fraction;
- no causal claim that live plants carry more or fewer butterflies, because moths dominate source detections, inspection effort is selected, and plant-product classes vary.

This is an empirical **arrival-stage identifiability and mechanism qualification**, not a new invasion outcome model. The next decisive data collection target is species-level intercepted *butterfly* identity with recorded original host plant/cargo and stage, **before** the species is reported established in receiving region, combined with independent wild larval reproduction time series and host resource mapping. Without that, a conditional arrival×host interaction coefficient is not estimable.

## Reproducibility and manuscripts firewall
- Frozen [source plan](../../provenance/reviewer_defenses/protocols/SOUTH_AFRICA_PLANT_IMPORT_INTERCEPTION_20261027.md).
- [Original full Figshare two original member verification](https://github.com/zuizui0223/chocho/actions/runs/40489503139) and [source-only schema artifact](https://github.com/zuizui0223/chocho/actions/runs/40489503139/artifacts/12728335315).
- [Deterministic order Lepidoptera by original commodity and inspection regime](../../scripts/analyze_sa_import_lepidoptera_detection_by_commodity_20261027.py), [CI source original run](https://github.com/zuizui0223/chocho/actions/runs/40490377309), [derived source numerical artifact](https://github.com/zuizui0223/chocho/actions/runs/40490377309/artifacts/12727506641) (the later rerun corrects original ex-vitro plantlets class).
- [Source-verified two original-239 butterfly imports and one Spodoptera temporal check](../../scripts/audit_SA_intercepted_butterfly_hosts_and_establishment_20261027.py), [seven unit tests and genuine original source CI](https://github.com/zuizui0223/chocho/actions/runs/40491506044), [source-only imported plant/taxon aggregate](https://github.com/zuizui0223/chocho/actions/runs/40491506044/artifacts/12729003017).
- The source dataset authors Saccaggi et al. 2021 already contributed the basic interception dataset and biosecurity framing. The source-positive rates, host-resource conceptual decomposition and links to original chocho 239 host identities here are an independent exploratory reanalysis; do **not** claim original discovery of interceptions or prior taxonomic fact.

All original shipment IDs, anonymous importer/exporter IDs, original raw consignment rows, and unpublished source WCVP/HOSTS edges are withheld from derived PR artifacts. Existing GEB manuscript, main and other chocho PRs remain unchanged.
