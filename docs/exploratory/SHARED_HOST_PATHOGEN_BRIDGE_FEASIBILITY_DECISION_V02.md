# Shared introduced host as pathogen bridge: source audit and hard decision

**Date:** 2026-10-10 | **Status:** HYPOTHESIS PLAUSIBLE; FOCAL SPECIES PAIR FAILS LOCAL CO-USE EVIDENCE GATE | **No new ecological cross-species effect**

## Why this differs from the current butterfly paper

The GEB submission reconstructs potential larval resource geography for 239 butterflies and potential exact-host sharing. It does **not** measure a plant-mediated cross-species pathogen transfer, infection establishment or viable adult offspring. A novel *ecological* result would require such independently observed or experimentally manipulated outcomes; the map alone cannot imply competition or disease.

## External prior art: what was already published

1. *Euphydryas phaeton* on introduced *Plantago lanceolata* shows lower cellular immune response and, in 2016 field samples, higher JcDV burden than on native *Chelone glabra*, but no clear between-host survival effect (Muchoney et al. 2022, DOI 10.1002/ece3.8723).
2. On the **same introduced plant and virus**, *Anartia jatrophae* showed experimental evidence of lower pathogen burden and higher survival than on native *Bacopa monnieri* (Muchoney et al. 2023, DOI 10.1111/ele.14162).
3. Host × dose contrasts comparing the above two butterfly species were **already published** (Muchoney et al. 2024, PMID 39159850). Do not call species differences new.
4. Within-*E. phaeton* plant-mediated horizontal/vertical transmission and infected-plant host differences are **already published** (Christensen et al. 2024, DOI 10.1002/ecy.4282).
5. Viral hotspots across *E. phaeton* landscapes tied to local host plant configuration and phytochemistry have also been studied (2025, PMCID PMC11858745). Again, no heterospecific donor result implied.
6. Per-capita parasite infection can **decline** under higher larval crowding even on shared foodplant substrates (2025, PMCID PMC12419877); more resource co-use does not mechanically mean more disease.

The hypothesis requiring new evidence is **heterospecific, source-verified infection transmitted through a shared introduced larval host and altering recipient adult production**, not mere pathogen susceptibility across species.

## Actual original-source audit conducted

The focal candidate `Junonia coenia × Anartia jatrophae × Plantago lanceolata` was selected before screening GloBI feeding observations. Both butterflies are in frozen chocho metrics. *Euphydryas phaeton* is **not** in the 239-species frozen species-metrics panel, despite being a useful prior-art focal taxon.

Workflow: https://github.com/zuizui0223/chocho/actions/runs/38015898425 (3/3 tests passed; source receipts saved). Prespecified `interactionType=eats`, `includeObservations=true`, up to two pages of 256 records per butterfly, original taxon, exact plant, explicit larval annotation, date, coordinates and independent source provenance. The API returned shorter first pages, so no 256-row pagination cap was triggered in this test.

| Species | GloBI eats rows retrieved | Exact species binomial records | Exact `Plantago lanceolata` records | Strict dated, georeferenced original larval feeding |
| --- | ---: | ---: | ---: | ---: |
| *Junonia coenia* | **217** | 216 | **0** | **0** |
| *Anartia jatrophae* | **33** | 32 | **0** | **0** |

Raw API page SHA-256s and search URLs are retained in the artifact. **Zero is zero in this specific indexed GloBI query, NOT evidence that either butterfly fails to use the plant.** Independent original experiments demonstrate both butterflies are physiologically capable of using *Plantago lanceolata*. The GloBI observation API is neither an exhaustive host list nor a controlled effort-to-detect dataset. The query therefore does not establish local co-use or rule it out.

**Consequence:** the focal `J. coenia × A. jatrophae` *natural cross-species pathogen bridge* cannot be promoted as an established study system from current local original observations. Stop at the predeclared eligibility gate. Do not silently choose different butterflies after seeing this result or claim an ecological disease effect.

## Independent original monarch cohort as a calibration of endpoints

An original, immutable public source exists for **a different system**, Ragonese et al. 2025 (DOI 10.1111/een.70010), field factorial of *Danaus plexippus* on tropical/swamp milkweed with ambient/elevated temperature and OE exposure.

Audited code: `scripts/audit_original_monarch_survival_bridge_calibration.py` and [successful workflow 38015801797](https://github.com/zuizui0223/chocho/actions/runs/38015801797).

- Author repository commit `74e3e2cd8079401c9ccd4359f9e39e2d52d1f7d2`; exact data blob `e7a625911d8dc9b871a52d7657314a5d05472216` was Git-blob-hash verified.
- **240 larvae initially assigned**, 30 per eight cells, **120 host plants** with two larvae per plant: plant is a shared exposure unit.
- One original accidental-injury fate (`ID=22b`) is coded NA and **was not treated as a natural mortality event**.
- Of **239 recorded adult-fate outcomes**, **219 produced adults** (91.6% of known outcomes). Within the 8 randomized study cells, adult eclosion ranged from **23/30** to **30/30**, with one 26/29 known outcome cell. This is a descriptive source audit, not a hypothesis test or new factorial discovery.
- Existing 2025 article already analyzed host/temperature/infection effects, including transmission/protozoan infection outcomes. This source does **not** test two consumer species sharing one introduced host or cross-species pathogen transfer.

**Critical distinction:** parasite inoculation assignment, infection assay status, adult eclosion, adult longevity and flight-capable adult reproduction are different outcomes. The original study includes adult emergence, not a demonstrated post-removal wild reproductive replacement rate. Never interpret egg counts, viral DNA or source presence as viable butterfly demographic benefit.

## Final research decision

The appropriate next scientific question remains:

> Can contemporaneous use of a *single, locally established introduced larval host* transfer a shared pathogen between butterfly species, and can the direction of the disease-mediated demographic effect be positive, zero or negative depending on plant physiology?

But a publishable answer requires a genuinely source-verified **local co-use / heterospecific transmission** system and offspring-to-adult endpoint, or a controlled experiment explicitly testing those outcomes with required approvals. The presently audited literature establishes within-species pathways and differences in host quality, but **not** this new ecological interaction.

**Stop further post-hoc rearrangements of the same HOSTS/WCVP geography or GloBI index to chase a positive source overlap.** Direct field observations and independently adjudicated biological source evidence would need to be collected. Keep the GEB source paper separate.

**Current work is complete as a feasibility rejection, not an empirical ecological mechanism discovery.** No changes to GEB PR #38.
