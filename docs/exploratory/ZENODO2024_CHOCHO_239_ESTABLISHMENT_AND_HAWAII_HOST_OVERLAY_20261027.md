# chocho: actual alien butterfly establishment records intersect reconstructed larval-host resource geography

**2026-10-27. Source-locked exploratory original analysis; NOT an invasion success-rate study.** Draft PR #84 stacked on PR #83, original GEB manuscript and main untouched. The question is whether established nonnative butterflies have documented relocated host-resource opportunities in the actual recipient region, and whether that source overlap is more discriminatory than ordinary host resource prevalence.

## 1. Species-level independent establishment-source accession verified

Turner, Blake, Mally & Liebhold, *International Non-native Insect Establishment Data*, **Zenodo 2024 v0.2.1**, [record 11519262](https://zenodo.org/records/11519262). The actual original source members, not wrongly expected filenames from the 2021 version:

- `occurrence_table_20230317.csv`: **13,047** alien insect establishment *positive* taxon×receiving-region rows, SHA256 `05483a2a9a7c926fb5105c44058ee52a4598039dc82c5bb10d1db284b3b66fe7`; original fields include `genus_species`, `family`, `order`, `rank`, `year`, `region`, `intentional_release`, `established_indoors_or_outdoors`, and `eradicated`.
- `attribute_table_20230317.csv`: **8,514** taxon attributes rows, SHA256 `d184315959746f867fef7f377d5532c383cc7d7fc1b3e135f6bf967c586275c5`.
- `source_table_20230317.csv`: **13,210** original primary bibliographic/support rows, SHA256 `ad2258907d07e6216b3a708d8bdee8780c6ba66a0b8046952b9342575e273ec4`.

[Source-only audit #40498308116](https://github.com/zuizui0223/chocho/actions/runs/40498308116) SUCCESS, five synthetic source tests. All original source public member checksums are validated. Original chocho resource denominators independently verified: **239 accepted butterfly species**, **206 with source WCVP resources expanded**, **14,553 newly represented butterfly×WGSRPD3 resource units**.

The strict butterfly binomial intersection of the original 239 with Zenodo 2024 alien-established positive rows gave **18 species** and **35 species×receiving-region positive records** across **11 named source receiving areas**. The larger original dataset has **838** unique source Lepidoptera species in its positive establishment roster (but most are moths). Species positives by receiving area:
Hawaii **10**, North America **5**, New Zealand **4**, Japan **3**, Okinawa **3**, South Africa **2**, Korea **2**, Ogasawara **2**, Europe **2**, Chile **1**, Galapagos **1**.

[Exact original source species join #40498879270](https://github.com/zuizui0223/chocho/actions/runs/40498879270) SUCCESS, 6 synthetic tests; [derived source positive-only case receipt](https://github.com/zuizui0223/chocho/actions/runs/40498879270/artifacts/12730973939).

**NO biologically negative outcome:** the other 221 original chocho butterflies are NOT DOCUMENTED in this version's eleven-region alien-established-positive source list. These are not failed introduction attempts or region-wide butterfly absences; source zones and inspection exposure are incomplete and uneven. In the 18 documented alien-established source butterflies only, the descriptive Spearman rho of number of documented alien recipient regions with original source **global** WCVP host-range newly gained units is **+0.378**, versus original documented host-family breadth **+0.531**. These are reporting/positive-list correlations without a measured arrival denominator, so no causal statement or p-value for establishment.

## 2. Hawaiian source-local overlap: 10 actually documented alien butterflies, all have WCVP *introduced wild* hosts

For the **ten Hawaii-positive original chocho butterflies**, overlay the **same source original accepted HOSTS plant-ID × WCVP Hawaii accepted botanical region code HAW**, separately native/introduced/wild unlisted. The [original 11-member author-only source SHA ZIP](https://github.com/zuizui0223/chocho/actions/runs/39006554609) and the exact HAW polygon code are independently validated, not an island substituted with all mainland USA.

[Reproducible HAW co-geography result #40499612707](https://github.com/zuizui0223/chocho/actions/runs/40499612707): **10/10** positive alien butterflies have ≥1 WCVP-introduced wild Hawaii source plant accepted as one of their HOSTS worldwide larval hosts; **7** have only introduced accepted WCVP wild-range hosts in HAW, **3** have BOTH native and introduced source hosts; none is native-only or wild-unlisted-only. The exact ten butterfly × accepted larval-host plant-source tally is **82 introduced-Hawaii host IDs, 5 native-Hawaii host IDs, 251 original known hosts NOT recorded as extant wild in Hawaii**. The latter 251 are NOT evidence plants/larval resources are absent from gardens, nurseries or agricultural plots.

Source examples of intro-only WCVP host geography include *Danaus plexippus*, *Agraulis vanillae*, *Hylephila phyleus*, *Vanessa cardui* and *Vanessa virginiensis*. The ten original butterfly×accepted-host sets were documented **worldwide**, and geographic plant availability does not demonstrate Hawaii larval feeding, oviposition, survival or carrying the insect on a plant shipment.

## 3. CRUCIAL original source pool control: ten of ten is structurally unsurprising

[Original-all-239 host saturation #40500199080](https://github.com/zuizui0223/chocho/actions/runs/40500199080) SUCCESS, three new source tests and all original 2596 accepted butterfly-host links independently verified.

Of **ALL** original 239 butterflies (regardless of alien establishment listing in Hawaii), WCVP original wild source plants already provide:

| Original 239 source Hawaii WCVP host geography | Number of butterflies |
|---|---:|
| only WCVP introduced wild hosts | **148** |
| both native and introduced wild hosts | **8** |
| WCVP native wild host only | **0** |
| no accepted plant host in HAW wild WCVP map | **83** |
| **total** | **239** |

Therefore **156/239 = 65.27%** of the entire original 239 have ≥1 Hawaii-introduced source host, even without labeling them as invaders or testing their arrival. Most butterflies with ≥21 known accepted hosts have some Hawaii-introduced plant in the reconstruction. To assess whether ten of ten adds discriminating power beyond the known host-list documentation scale, we FROZE source accepted HOSTS breadth strata **1–5, 6–20, 21–50, ≥51**, fixed number of Hawaii source positives in each original stratum (1,2,5,2) and calculated exact artificial hypergeom label exchange among original 239 **without calling unlisted species failed arrivals**. The artificial expected number of ten documented Hawaii aliens with any WCVP-introduced host is **8.998**, source central 95% interval **7–10**. Observed 10 is not extraordinary under this host-richness-matched artificial source reference (mathematical right tail 0.2695, **NOT an ecological p-value**, since butterfly natural dispersal, propagule arrival and biogeographic source pools were not randomized).

Original host richness bin source counts:
- 1–5 accepted hosts: **54/119** have Hawaii-introduced wild hosts; **1/1** source Hawaii aliens does.
- 6–20 hosts: **61/79** source support; **2/2** positives.
- 21–50 hosts: **35/35** source support; **5/5** positives.
- ≥51 hosts: **6/6** source support; **2/2** positives.

Thus the apparent dramatic “all established aliens have introduced plant hosts” is **partly a generic consequence of globally popular introduced host plants and larger documented host portfolios**. A positive-only registry plus a 239-species arbitrary global species pool does not identify a real ecological establishment effect, and no invasion success probability, orphan-opportunity percentage or fraction attributable to seed vs live plant can be estimated.

## 4. Independent transport-stage overlap and chronology highlight why the two channels must remain separate

The original SA Figshare 1994–2019 source had exact intercepted original 239 butterfly labels **two species**, both already marked `present` in historical/current SA presence metadata: *Lampides boeticus* on bean cargo (*Phaseolus vulgaris*; accepted worldwide HOSTS larval host) in **2018**, *Pieris brassicae* on *Brassica* sp. fresh foliage in **1997**. Source first documented *P. brassicae* alien establishment in South Africa in Zenodo **1994**, BEFORE its detected 1997 interception, so the later intercepted butterfly cannot be used as pre-establishment arrival causing original founding. *Lampides boeticus* is not listed as new alien established SA in the 2024 source, even though original SA specimen metadata marks its presence. Exact source contradiction may reflect species natural status or registry scope, NOT colony failure. [South Africa original source case #40492231178](https://github.com/zuizui0223/chocho/actions/runs/40492231178).

The Zenodo 2024 year field also contains **1758 for *Vanessa atalanta* in Hawaii**, a suspicious ancient record that should NOT be treated as verified arrival year. One source-positive species×recipient entry, *Pieris brassicae* in New Zealand (2010), has original `eradicated=Yes`: the positive alien historical occurrence should NOT be silently interpreted as a presently persistent wild population. Historical first-detection years, interception dates, natural arrival, subsequent establishment and local extirpation are separate endpoints.

## 5. Mechanistic status and decision for the paper

What is independently verified now:
1. Plants redistributed by humans create reconstructed larval-host food geography for a large original butterfly host-network universe (GEB main).
2. A separately published **actual alien-established Lepidoptera positive roster** can be joined at exact species to the original 239 butterflies (18 taxa across 35 positive source recipient units).
3. Hawaii WCVP introduced host plants are found for **all ten** source-documented alien-established butterflies (seven without any Hawaii-native wild plant source host), but that is **not surprising after conditioning on original host-list breadth**: expected ~9/10 among source taxon-exchanged candidates.
4. Real South African border detections prove insects can be observed on known-host goods, but for the two original butterfly species the dated specimens do not supply first pre-establishment import events; the recorded changed-status *Spodoptera frugiperda* interception years 2018–19 also follow recorded establishment 2017.

What remains unknown: the original 239 study's butterfly propagule pressure/arrival frequency by species×destination×year, transport commodity/historical insect stage, dated recipient cultivated-only host resource, actual repeated local larvae on hosts, and standardized true failed-arrival outcomes. A meaningful two-process arrival × postarrival-resource interaction coefficient is STILL not identifiable. Treat this as an important new **positive-only geographically aligned resource-feasibility analysis plus falsification of naive establishment-risk conclusions**, not a validated invasion mechanisms paper.

**Next biological design:** measure actual recipient-specific inspection/interception positives BEFORE first occurrence, carefully distinguish recorded interception from arrival and persistent wild reproduction, and assemble comparable inspected/shipment effort in the SAME species/destination/time units; do not simply regress presence on introduced-host availability in all butterfly species, where source labels contain no true negatives.

## Sources and reproducibility
- [Zenodo 2024 original 13,047 source schema CI](https://github.com/zuizui0223/chocho/actions/runs/40498308116), [source positive 18 study taxa result](https://github.com/zuizui0223/chocho/actions/runs/40498879270).
- [Ten Hawaii positive source HAW plant native/introduced footprints and CI](https://github.com/zuizui0223/chocho/actions/runs/40499612707), [source-only derived HAW JSON](https://github.com/zuizui0223/chocho/actions/runs/40499612707/artifacts/12732390563).
- [All 239 source Hawaiian host opportunity and exact host-richness matched source-reference test](https://github.com/zuizui0223/chocho/actions/runs/40500199080), [derived source-only artifact](https://github.com/zuizui0223/chocho/actions/runs/40500199080/artifacts/12731634421).
- [Original source-preserving 2024 Zenodo exact species parser](../../scripts/analyze_zenodo_2024_established_butterflies_original239_20261027.py), [HAW ten established host-status source run code](../../scripts/audit_Hawaii_established_2024_butterfly_WCVP_host_coverage_20261027.py), [all source239 background](../../scripts/audit_Hawaii_all239_WCVP_host_opportunity_baseline_20261027.py).

Everything remains DRAFT in PR #84, stacked beyond #83 and outside the GEB manuscript/main. No author-only accepted plant×butterfly pair identity lists or original specimen-level customs data are published.
