# Does anthropogenic host redistribution create local butterfly multi-consumer hubs?

**Date:** 2026-10-27. **Status: completed, reproducible original-source HOSTS/WCVP macro documentation audit; hypothesis NOT established as an ecological process.** Draft PR #80 stacked on PR #79; original GEB submission and main unchanged.

## Main positive source finding: plant introductions can be reported as resources for several butterfly species

This is a different estimand from chocho's GEB main reconstructed +54.9% potential host-resource geography. Instead of a host-presence envelope for one butterfly, the new analysis counts **different original 239 butterfly taxa recorded by HOSTS as users of the same accepted WCVP plant taxon in the SAME original HOSTS Location**.

Verified original source [GitHub Actions #40470959824](https://github.com/zuizui0223/chocho/actions/runs/40470959824): **SUCCESS, 13 synthetic source tests**. Exact original SHA-verified 11-member author-only WCVP/HOSTS source ZIP, original HOSTS commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`, original accepted 239-butterfly manifest, and original country/TWDG geography 969 disjoint area atoms. Original 3,586 butterfly×WCVP accepted plant×original source Location edges reconstructed as **3,035 distinct plant×Location nodes** before botanical status exclusions.

The original primary botanical status source gate uses ≥95% *country mapped area* and ≥95% source WCVP plant-origin area dominance. Unknown national, mixed, broad source labels are NOT imputed. Hawaii and West Malaysia retain explicit accepted WGSRPD3 source units. Excluded source records remain in the full 3,586-edge accounting.

### Documented resource sharing, original source ≥95% WCVP botanical classification

| Original plant-region status | Original accepted plant×Location host nodes | Nodes with ≥2 documented butterflies | Percentage of documented nodes | Butterfly×plant×Location source edges | Edges on shared host nodes | Cross-butterfly-genus shared nodes |
|---|---:|---:|---:|---:|---:|---:|
| **WCVP introduced** | **141** | **31** | **21.99%** | **179** | 69 | **22** |
| **WCVP native** | 329 | 54 | 16.41% | 399 | 124 | 39 |
| **WCVP-unlisted** | 247 | 39 | 15.79% | 303 | 95 | 16 |

**Supported descriptive claim:** actual HOSTS records document that introduced-range plant hosts can be used by multiple butterfly species within an original named source region, with some cross-genus co-use. This is *network co-reporting*, not measured contemporaneous co-occurrence or proven interspecific competition; the source named region may be a whole country.

The observed difference of **+5.58 percentage points** in shared-plant node frequency (introduced vs native) is an **unadjusted catalogue comparison**. The number of ORIGINAL HOSTS source raw records also differs: 196 raw records supporting 141 introduced plant nodes (1.39 per plant×Location) and 430 raw records for 329 native nodes (1.31 per plant×Location). Repeated reports and coarseness can increase the chance that two butterflies are indexed against the same plant.

## Source-matched falsification: why the raw +5.6 pp is not a demonstrated invasion mechanism

Before inspecting these results, the [protocol](../../provenance/reviewer_defenses/protocols/INTRODUCED_HOST_LOCAL_SHARING_20261027.md) froze exact plant-status label-exchange references. We only compare source-accepted plant nodes that have **both WCVP introduced and native plants reported under the same Location and the same accepted plant family**; all other source source strata are excluded from the MATCHED contrast but counted in the source descriptive sample.

| Original ≥95% WCVP test / reporting reference | Target introduced plant nodes in matched strata | Target observed shared nodes | Exact source conditional expected | Central 95% source reference interval | Regions where genuine permutation possible | Source information gate |
|---|---:|---:|---:|---:|---:|---|
| Same HOSTS Location + plant family | 69 | **17** | **14.05** | 10–18 | **7** | FAIL (need ≥10 locations and ≥15 informative strata; only 11 strata) |
| Plus matched original HOSTS raw-record count bins (1,2,3–5,≥6) | 61 | **14** | **13.30** | 12–15 | **2** | FAIL |
| Native vs WCVP-unlisted *generic documentation control*: same Location + family | 124 | **26** | **21.15** | 16–26 | 11 | PASS a source-level documentation gate, NOT biological specificity |
| Same + raw record count bins | 100 | **16** | **16.00** | 15–17 | 1 | FAIL |

The introduced-source observed 17 is **inside** its artificial 95% reference interval (10–18), with mathematical fixed-source right tail 0.120. After controlling explicit original HOSTS source reporting counts, observed 14 vs expected 13.30 (right tail 0.400) is even less distinctive, yet only **two** locations permit a genuine label exchange and the result cannot be generalized.

The source unlisted control exhibits a comparable descriptive excess (26 vs 21.15 in the location+family reference) which disappears after raw report matching (16 vs 16.00). This pattern is compatible with source publication concentration, host taxon ecology, or plant geographic opportunity rather than special new ecological co-use following introduction. Fixed-label hypergeometric mathematical tails are **not valid ecological p values**, because plant status is not randomly assigned and butterfly distributions and source reports are not exchangeable.

### Sensitivity at the originally fixed 80% WCVP country plant-origin threshold

At 80% area dominance, introduced plant nodes rise **141→152**, shared nodes **31→37** (24.34%). Native nodes rise 329→348, shared nodes 54→59 (16.95%). Location×family-matched introduced nodes 86 have **25 observed shared vs 21.73 expected** (reference 17–26; mathematical right tail 0.119), but just **9** permutable locations. Adding source raw reporting bins gives **17 observed vs 16.13 expected**, reference 14–18, just **3** permutable locations (information gate FAIL). Thus the direction of unadjusted documentation difference is reasonably stable, **while the ecological introduced-host hub mechanism remains statistically unidentifiable**.

## Concentration in a handful of reported geographies

In the original primary 95% source status, **17 original Location names** contribute introduced-plant records; locations with the most recorded introduced-host butterfly links are:

| HOSTS recorded source locality | Source butterfly×introduced-host links | Introduced plant×Location nodes | Shared host nodes |
|---|---:|---:|---:|
| West Malaysia | 29 | 23 | 4 |
| Finland | 25 | 17 | **7** |
| Taiwan | 20 | 18 | 2 |
| Sri Lanka | 15 | 11 | 4 |
| Thailand | 15 | 10 | 4 |
| USA | 14 | 8 | 4 |
| Cuba | 13 | 11 | 2 |
| Hawaii | 12 | 12 | 0 |

Finland alone provides **7/31** source introduced shared plant nodes; these are records for original broadly classified insect taxa, not a proof of real larval interactions at the same plant/season. Finland plus West Malaysia, Sri Lanka, Thailand and USA provide **23/31** shared nodes, demonstrating source-geographic concentration. The source `Location` granularity and publication effort are large downstream limitations.

## Biological case: reported native/alien consumer sharing in Chile, not experimental competition

An original Chile historical source [INIA retrospective account](https://books.instituto-idema.org/sites/default/files/Control%20Biologico%20de%20Plagas%20en%20Chile.pdf) reports larval feeding by invasive *Pieris brassicae* on *Tropaeolum majus* near Viña del Mar–Concón in 1971. Olivares & Hormazábal (1998), [original Revista Chilena de Entomología](https://www.insectachile.cl/rchen/pdfs/1998v25/Olivares_Hormazabal_1998.pdf), explicitly report larval co-use of *T. majus* by indigenous *Tatochila autodice blanchardi* and *P. brassicae*, without identifying direct competitive suppression, individual survival effects or the exact shared plant observation time. Kew lists *T. majus* as Peru-origin and introduced in Chile Central, but the exact shared-field source WCVP botanical subregion status remains not independently verified.

This is a **hypothesis-generating direct multi-consumer ecology example**. The new global HOSTS source macro is structurally analogous, not independent field replication, and includes no consumer-native/native classification. Independent butterfly work [Japanese butterfly co-occurrence/host sharing](https://pubmed.ncbi.nlm.nih.gov/29392443/) has found patterns more consistent with shared-resource filtering than resource competition, demonstrating that these are competing interpretations rather than synonyms.

## Main defensible conclusions and next biological gate

1. The original known-host redistribution greatly expands **potential resource geography** (+54.9%, original GEB). This is not equivalent to local actual adoption.
2. Documented **locality-level host-network sharing** occurs for **31 source introduced plant×Location nodes (22 cross-butterfly-genus)** across the original 239 butterflies, and its unadjusted recorded frequency exceeds WCVP-native plants, but only in a geographically narrow, literature-selected sample.
3. Under locality/family and especially original raw publication/reporting-intensity controls, observed overlap is close to conditional documentation expectations and there are too few informative localities to claim introduced hosts specifically generate greater shared butterfly interaction networks.
4. To test real competition or consumer coexistence, require **same-site, same-period wild juvenile records**, butterfly-native/invader identifiers, host patch occupancy, replicate shared vs unshared plants, larval success/abundance, parasitoid and habitat controls. Country-scale HOSTS co-reports cannot substitute.

**The current viable general research question is not simply whether introduced hosts boost shared network counts, but WHEN a geographically redistributed host is converted into a realized multi-consumer resource, and how that changes realized coexistence and possible interaction strengths.** At present the source supports only geographic report-level co-use, not the strength of these ecological pathways.

## Reproducibility and firewall

- [Froze source and structural matched control before observation](../../provenance/reviewer_defenses/protocols/INTRODUCED_HOST_LOCAL_SHARING_20261027.md).
- [Source-pinned exact 239 butterfly host-node reanalysis](../../scripts/analyze_introduced_host_local_multi_consumer_network_20261027.py).
- [No-network exact native vs introduced plant-node null tests](../../tests/test_introduced_host_local_multi_consumer_network_20261027.py).
- [Successful original source CI and numerical results artifact](https://github.com/zuizui0223/chocho/actions/runs/40470959824), [aggregate-only ZIP](https://github.com/zuizui0223/chocho/actions/runs/40470959824/artifacts/12718308477).

No original copyrighted HOSTS×WCVP individual interaction matrix, host IDs, butterfly pairing or individual article text was redistributed. The original GEB manuscript and previous PR #75–#79 remain unchanged.
