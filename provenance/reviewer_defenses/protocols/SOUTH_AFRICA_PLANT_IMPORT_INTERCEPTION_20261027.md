# Independent intercepted-insect commodity exposure and subsequent establishment from actual inspected consignments

**2026-10-27 source addendum, frozen before first reading original South African Figshare CSV rows.** This source complements, but cannot be imputed into, the 239-butterfly GEB resource map. Keep original manuscript/main untouched; draft PR #83.

## Dataset identity (source not outcome)
Saccaggi, Arendse, Wilson & Terblanche (2021), *Scientific Data* 8:83, DOI `10.1038/s41597-021-00869-z`, public Figshare `10.6084/m9.figshare.14393723`. Original 1994–2019 South African import-laboratory records of plant-product interceptions: 26,291 initially recorded consignments, **25,279** finalized sample records, **13,731** contaminant observations, and two public CSV source products:
- `List of contaminants on SA plant imports 1994–2019.csv` (raw sample×commodity×origin×time and 79 typed taxon result fields). Source result `NA` is **not tested**, source `0` is **tested but not detected**. Never turn `NA` into tested negative.
- `Metadata of contaminants on SA plant imports 1994–2019.csv` (taxon species and higher group, whether recorded in South Africa at interception and since established). Publication reports **7 species across all taxa** initially not established in South Africa that later established, but does NOT state how many are Lepidoptera. Cannot assume this count for butterflies.

Original commodity distinguishes `propagation material` (plants, live cuttings, rooted stock) from `fresh produce`, `seed`, `grain`, `soil` and other categories; original three sampling regimes `intervention`, `audit`, `postentry quarantine` are NOT randomly equivalent. Source reported interception is a sampled consignment observation, not true arrival frequency of all incoming trade. A tagged contaminant is not necessarily a larval host user of the named commodity; incidental contamination and adult movement are possible.

## Stage A: raw source access/schemas (first)
Obtain original versioned Figshare API metadata, file download links, byte md5/SHA, original schema names, record grain, positive/zero/NA semantics, sample regime and raw status vocabulary. Never redistribute sample-level importer/exporter codes. If an ID or file unavailable, stop and report HTTP/error/source, not invented counts. Do not assume the two CSVs are linked by sample ID: validate taxon lookup mechanism (79 taxon cells in sample file may store a taxon name directly or code).

## Stage B: source-backed arrival commodity contrast, if schemas support it
Compare positive **Lepidoptera taxon detection** across **distinct INSPECTED** sample IDs, grouped a priori into:
1. `LIVE_PROPAGATION`: live plants, rooted/leafy cuttings, bulbs/plantlets (distinct source original types).
2. `SEEDS_GRAIN`: seeds, grains, beans, nuts.
3. `FRESH_PRODUCE`: fresh fruits/vegetables.
4. `OTHER_OR_UNKNOWN`: other raw classes, separately marked.

Use only samples with Lepidoptera actually tested, and explicitly report **tested**, **positive**, **not tested** separately by commodity and sample regime (intervention, audit, PEQ). No mixing raw taxon-negative assumptions. Estimate detection rate difference within standardized sampling regime and overlapping years (if identification adequate), otherwise descriptive only. Factor record multiplicity is sample-level: 13,731 organism records do NOT represent 13,731 separate positive samples. Source group 'insects' is broader than Lepidoptera, and most recorded insects are not original study butterflies.

## Stage C: establishment observation
Source species-level metadata says whether each detected taxon was native/established/absent THEN, and later status. Only if original columns permit an unambiguous 1:1 join to sample-level Lepidoptera IDs and status chronology, count species intercepted by commodity that later established in South Africa. Never infer causal establishment probabilities because follow-up time, regulations, inspection selection and exposure unknown. Seven later established across all taxa is a source-wide reference, not seven butterflies.

## Scientific decision gate
Present observed commodity-specific Lepidoptera detection only if at least **30 tested samples in two or more commodity classes**, at least **5 Lepidoptera positives total**, and **10 years with valid comparable reporting**. If sample source `NA` is not interpretable, rate estimates not issued; first assess the exact source column values. Native/alien butterfly panels must remain separate and exact taxonomic identity confirmed before any chocho species intersection.

Link to WCVP host geography only at species×receiving-location if original source identifies exact butterfly host and match, which may be NONE. Figshare commodity/source richness alone cannot establish that a given larval plant genus carried a given invasive butterfly species, and original 239 resource figures must remain unchanged.

This first analysis is an **independent empirical test of whether early-stage Lepidoptera detections occur in multiple plant-product route types**, not yet the novel causal test of coupling between arrival route and postarrival host resources.
