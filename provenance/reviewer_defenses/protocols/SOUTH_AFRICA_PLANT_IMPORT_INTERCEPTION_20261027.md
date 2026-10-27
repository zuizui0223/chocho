# Independent intercepted-insect commodity exposure and subsequent establishment from actual inspected consignments

**2026-10-27 source addendum, frozen before first reading original South African Figshare CSV rows.** This source complements, but cannot be imputed into, the 239-butterfly GEB resource map. Keep original manuscript/main untouched; draft PR #83.

## Dataset identity (source not outcome)
Saccaggi, Arendse, Wilson & Terblanche (2021), *Scientific Data* 8:83, DOI `10.1038/s41597-021-00869-z`, public Figshare `10.6084/m9.figshare.c.5180681`. Original 1994–2019 South African import-laboratory records of plant-product interceptions: 26,291 initially recorded consignments, **25,279** finalized sample records, **13,731** contaminant observations, and two public CSV source products:
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


## Original source grain addendum — after successful schema-only retrieval, before commodity tests
Original Figshare **collection** c.5180681 member articles #13118177 (sample file ID #25166393) and #13668839 (metadata file ID #26235020) are exact. Source file SHA256:
- samples (26,295 ORIGINAL UNFILTERED CSV rows, **92 columns**): c12d4af0c96a247d520f71ba6c2dfa67eacb33b53dd5aa329265c850d3fe15da
- taxa (878 rows, **14 columns**): 0f7ad1b9b44dad51a53d55a1f6bad193d623b0b6f025dfcf640c

The paper's 25,279 **FINAL FILTERED** sample count is not the raw file's 26,295 rows. No sample subset or 25,279 estimate will be fabricated without replicating the authors' exact exclusions. The original sample has exact columns sample.source, crop.commodity, year, and Insecta: Lepidoptera. Under source-only catalogue parsing, Insecta: Lepidoptera is 15,673 literal tested negative '0', 10,443 literal not tested 'NA', 178 nonzero/non-NA potentially positive, and one blank. The 178 may contain taxonomically broad larval categories, rather than butterfly species. Thus the observed source denominator among nonmissing is 15,851 tested, **not** all 26,295 raw rows.

For a source-transparent primary analysis, use:
- detailed raw commodity \`crop.commodity\`, preserved exactly in output, and source \`sample.source\` in PEQ/audit/intervention. First require all 30 original commodity vocabulary values are accounted for in original source metadata; unknown commodity types retain original label and are NOT merged into a preferred transport type after outcome analysis.
- Predeclared **three** biologically interpretable *prespecified* broad routes, with a fourth explicit mixed/other pool: **LIVE_PROPAGATION** plants, cuttings/budwood, bulbs/rhizomes/tubers and in-vitro propagules; **SEED_GRAIN** seeds, grains, beans, nuts; **FRESH_OR_DRY_PRODUCE** fresh/dried fruits, vegetables and foliage; **OTHER_UNKNOWN** all remaining commodity labels or missing. Live foliage is NOT assumed a rooted seedling or proof of oviposited larvae on that specific leaf. Output original granular commodity row counts always, so group mapping is auditably non-silent.
- Count **any positive insect from order Lepidoptera** for the exact taxon result field as the primary stage-1 interception outcome (many moths). Maintain special interpretation for 'Lepidoptera larva' as juvenile category and any binomial-specific names separately; DO NOT claim all positives are larval butterflies. Report positives, TESTED negatives, NOT_TESTED, and MISSING values per source regime and commodity.
- Source sampling regimes **audit**, **intervention**, **PEQ** are preexisting nonexchangeable strata. Report crude detection rates by product type and within each regime; treat PEQ x live plants vs audit x seed comparisons as **not comparable in selection**. If possible compare routes only in common regime and year bins (1994–2004; 2005–2019). No causal confidence interval or biological import-frequency estimator from preferential inspection.
- For exact species-level Lepidoptera names, compare to original chocho frozen 239 butterflies only by exact accepted binomial; lack of match is NOT proof no butterflies are intercepted. Source taxon metadata historical/current South Africa status can only classify taxon (no linkage to an individual founding event). A raw source sample positive is not a successful established insect.
