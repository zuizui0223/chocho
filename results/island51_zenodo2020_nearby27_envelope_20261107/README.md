# Can all 27 geographically close but UNVERIFIED external AIS port names rescue original-island direct route coverage?

**Completed 2026-11-07; exploratory source-only sensitivity. NO additional verified voyages, ports or butterfly models.** This report adds a deliberately overinclusive *hypothetical identity assumption* to the authoritative source crosswalk in [PR #141](https://github.com/zuizui0223/chocho/pull/141), without changing the original positive or adverse island/port and GEB results.

## Source and fixed comparison

Duvallet et al. (2024), public AIS-derived 2020 origin/destination triplet source, [Zenodo DOI 10.5281/zenodo.10783698](https://doi.org/10.5281/zenodo.10783698), CC BY 4.0. The source-original 2020 ZIP MD5 is af27e1bef1f706e2e29ed53aee47ff10, and 645,581 source voyage rows (departure dates 2020 Jan 1–Oct 7) were independently re-read with complete original 194-WPI-port geographic context.

The frozen parent PR141 linked **46 external endpoint nodes** uniquely to 46 original WPI ports by exact canonical name plus ≤25 km original source geography, representing 26 original island groups. The 46 verified endpoint nodes cannot be changed by this scenario. Among the previously non-matched names, exactly **27 external source nodes** lie within **10 km** of a geographically nearest original WPI port, despite failing the exact name requirement; 18 are within 5 km. The [27-row candidate list](../island51_zenodo2020_portvoyage_audit_20261107/source_only_10km_unverified_port_candidates.csv) explicitly marks **all 27 identities UNVERIFIED**. Felixstowe versus HARWICH, London versus GRAVESEND and Koh Si Chang Transhipment Area versus MANILA exemplify why proximity-only matching is not valid port identity.

The [source-frozen protocol](PROTOCOL.md) specifies two scenarios: (A) the immutable original 46 **VERIFIED ONLY** WPI-matched endpoints, and (B) a deliberate *what-if* where all 27 unverified source nodes are temporarily assigned their nearest original WPI port solely to assess whether resolving all close geography names could plausibly attain the **80 original island-pair link minimum**. This does not validate the 27 links, does not make those source OD counts real verified voyages, and must never be substituted for the genuine port-to-port OD evidence.

## Decisive source-coverage result

| Original 51-island source coverage | **VERIFIED (observed matched)** | **HYPOTHETICAL: all 27 unverified names assigned to nearest WPI port** |
|---|---:|---:|
| External source nodes assigned a WPI identifier | **46** | 73, including **27 unverified** |
| Distinct WPI ports reached by the assignments | **46** | 64, including unverified identity assumptions |
| Original island groups reached | **26** | 32, including unverified identity assumptions |
| 2020 source rows with both WPI endpoint labels under scenario | **442** | 1,573, with **1,131 candidate-dependent** |
| Cross-original-island source OD rows | **133** | 298, with **165 candidate-dependent** |
| Different original unordered island pairs with a source OD episode | **12** | 28, of which **16 newly require unverified name identities** |
| Meets prespecified minimum 80 VERIFIED island-pair links | **No** | **No**, even if all 27 false-positive possibilities were granted |

**Conclusion:** Even if all 27 geographically close but name-discordant source nodes were truly the same WPI ports, the original 51-island source would cover only **28 of the 80 island-pair links** required by the source-only feasibility gate. It is still short by **52 pair types**. Thus this particular proximity-only name-correction opportunity **cannot** rescue the globally comparable direct-voyage model. The source actual *verified* network remains **12 original island-pair links and 133 cross-island source records**, unchanged.

The 165 extra source triplet rows are NOT confirmed ship voyages at the assumed original WPI ports. They are reference calculations that depend on incorrectly accepting all name-discordant nodes. Of those 165, **99 augment original already observed pair 15–26**, 22 would create assumed pair 15–206, and 13 would create assumed pair 15–152. These numbers demonstrate concentration of the counterfactual and why simply counting additional event rows overstates independent network coverage. All 1,275 original possible island pairs, including the 1,247 pairs without hypothetical source links, retain their indexing status; a zero remains **not evidence of no ship traffic**.

## Independent direct-source audit

The original public Zenodo source ZIP was downloaded and its publisher MD5 checked. The [source-only scenario runner](../../scripts/audit_island51_zenodo2020_nearby27_whatif_coverage.py) rereads the full original 645,581 source triplet rows and MUST reproduce exact parent PR141 verified totals (**442 both-endpoint source rows, 133 between-island rows, 12 unordered island pairs, 94 WPI directed port links**) before computing the hypothetical scenario. Original WPI IDs, names and island groups were not rewritten.

A separately written [raw-source independent auditor](verify_27_envelope_from_raw.py) independently reparsed the full 645,581-row ZIP, validated the original 46 exact WPI source labels and 27 fixed geography-only labels, and reconstructed every one of the [1,275 original verified-versus-hypothetical source-pair rows](original51_verified_vs_UNVERIFIED_geonear_alias_2020_OD_pair_scenarios.csv), including original real counts and extra source rows that rely on false-identity assumptions. The [machine-readable source receipt](source_only_geonear27_whatif_receipt.json) and [independent raw-source audit receipt](independent_raw_source_audit.json) preserve all enumerated source-coverage totals, source publisher MD5, the early research stopping criterion and explicit NOT-VERIFIED identity flags. No vessel MMSI data, entire vendor ZIP or private butterfly-host network is republished.

## Inference limits and decisive next data need

This is **not** an additional positive ecological shipping result: it is a source-availability test that explains why the current public 2020 AIS origin/destination extract cannot demonstrate a globally usable 51-island shipping network even after optimistic local-name harmonization. It does not establish actual ship voyage absence, 2001–2010 maritime routes, plant/larval-host transport or butterfly establishment. The 27 candidates may be different ports, and source 2020 Jan–Oct cannot be equated to earlier botanical or butterfly records. We do NOT fit a butterfly response model on the 28 pseudo-links, because they are unverified and far below the prespecified support minimum.

**Next empirical need:** a broader independently identified port-call/OD archive with verified port identifiers and historical timestamps, rather than relaxing WPI names or zero-imputing source voyage counts. Earlier PR137 geographic AIS null-source limitation, PR138 194-port early GENERAL plant data, PR140 cross-community conditional AIS results and GEB +54.9% *potential* larval-host resource portfolio expansion remain separate observations, not negated or newly validated by this source-only feasibility bound.
