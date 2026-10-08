# Fabaceae-feeding Pieridae: independent larval-use feasibility decision v0.1

**8 October 2026.** Exploratory, not a confirmatory test. Scope: `analysis/resource-homogenization-v01`; GEB main manuscript unchanged.

## Why this question

The previous strict butterfly Family × botanical Family link-identity null was modest (original excess +0.010538, p=0.030; second seed excess +0.009536, p=0.026). Removing the **Fabaceae** host family reduces the excess to +0.005236 (p=0.180) and +0.004641 (p=0.196), respectively; removing all **Pieridae** also weakens it (p=0.116 and p=0.156). Thus the broad consumer-specific claim is not robust to some taxonomic exclusions.

The frozen 239 butterfly cohort contains 43 Pieridae. Of these, 24 species have 238 exact accepted Fabaceae HOSTS associations. The current resource reconstruction gives 2,798 newly added Pieridae × WGSRPD3 region opportunities, of which 1,377 are uniquely supported by the presence of a contemporary Fabaceae known host. These are opportunities in an inventory, not fitness measurements or realized larval populations.

## Outcome-blind species panel and original sources

Panel protocol: `PIERIDAE_FABACEAE_LARVAL_USE_PILOT_PROTOCOL_V01.json`. All 24 Fabaceae-feeding Pieridae were selected before consulting GloBI source outcomes, based on HOSTS and butterfly taxonomy only.

Two reproducible GitHub runs:

- GloBI API extraction, [37715427511](https://github.com/zuizui0223/chocho/actions/runs/37715427511), artifact `butterfly-pieridae-fabaceae-larval-use-v01`.
- Original source verification, [37715634355](https://github.com/zuizui0223/chocho/actions/runs/37715634355), artifact `butterfly-pieridae-fabaceae-original-source-audit-v01`.

Feasibility counts:

| Stage | Records / coverage |
|---|---:|
| Selected butterfly species | 24 |
| Source exact butterfly matches in API responses | 4,788 |
| Dated/georeferenced plant-feeding candidates | 2,613 |
| Explicitly annotated larval records | 63 |
| Larval records from non-HOSTS provenance mapped to WGSRPD3 | 58 |
| Exact known Fabaceae host is WCVP-introduced in the region | **7** |
| Number of butterfly species represented in those 7 records | **5** |
| Number of WGSRPD3 regions | **4** |
| Strict introduced-only resource regions with larval evidence | **0** |

The frozen pilot gate was **>=10 introduced-host observation records from >=3 butterfly species and >=3 regions**. With 7 records, it **FAILED**. `Catopsilia pomona` reached the 2-page API query cap, so zero subsequent records for that species are **unknown**, not biological absence.

## Why the seven are not seven confirmed field feeding events

Three older GloBI rows (*Colias croceus* 1960 and 1975, *Colias hyale* 1975 on *Medicago sativa* in Germany) came from a **literature-based 2025 BEXIS interaction compilation**, with no unique original observation IDs; one event predates the title's stated 1975–2005 coverage. They are separately indexed literature associations, **not independently verified field photographs**.

The remaining four GloBI rows have original iNaturalist observation identifiers, all found in the current iNaturalist API, with matching current butterfly taxa, years, photos, larval annotations, research grade, and `captive=false`. The original observation-field evidence differs:

| Butterfly | Introduced Fabaceae plant | Original iNaturalist observation | Source field interpretation |
|---|---|---|---|
| *Eurema hecabe* | *Caesalpinia pulcherrima* | [244915185](https://www.inaturalist.org/observations/244915185) | `Host Plant ID`; plant also marked `Cultivated`, without an explicit Eating/Feeding field |
| *Phoebis philea* | *Senna alata* | [53608634](https://www.inaturalist.org/observations/53608634) | Explicit `Eating` field and larval annotation |
| *Phoebis sennae* | *Senna obtusifolia* | [180437332](https://www.inaturalist.org/observations/180437332) | `Insect Host Plant` and `Interaction->Herbivore of` assertion |
| *Phoebis sennae* | *Senna obtusifolia* | [234133454](https://www.inaturalist.org/observations/234133454) | Same two host/herbivory assertions, with an associated plant observation link |

**One of the four originals has an explicit Eating field**, two have an observer-declared Herbivore-of field, and one only an observer-declared host field. Source metadata have been checked; photographs have NOT been independently adjudicated for actual ingestion, nor has larval survival been tested. No original-observation proof of new butterfly-range colonization follows from these records.

## Claim boundary and decision

1. The overall plant-family-conditioned network result is a real **conditional structural signal**, but weak/lineage-sensitive. The most defensible focal class is Fabaceae-related Pieridae; it is **not** established as a causal Pieridae × Fabaceae interaction.
2. The source-traced larval-use pilot fails its predeclared coverage requirement, with **zero strict introduced-only resource cells**. Do **not** promote realized range filling, ecological rescue, or a new ecology manuscript based on these records.
3. Even verified larval eating is not larval performance, completed development or lifetime fitness. Yoon and Read (2016; doi:10.1007/s00442-016-3560-2) already showed mean costs of exotic hosts in a meta-analysis of 76 studies. Braga (2023; doi:10.1016/j.cois.2023.101074) reviews cases of traps and beneficial alien-host use. The question of which imported resources contribute useful demographic opportunity remains **unanswered by the reconstruction**.
4. The earlier 11-species independent host-use pilot exposed at least one introduced-only false positive caused by a missing native butterfly–host edge in HOSTS; do not interpret inventory absence as ecological absence.

**Stop rule:** Do not expand the same presence-only GloBI panel simply to find enough positive cases. A genuinely new outcome requires independent larval development/survival, verified field reproduction, or unbiased occupancy with effort and historical baselines. Keep all exploratory results off main GEB until those data are available.
