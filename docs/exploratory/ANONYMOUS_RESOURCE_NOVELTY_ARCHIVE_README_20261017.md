# Anonymous butterfly host-resource archive — offline reviewer reproduction

**Scope:** the same 239-butterfly, single-paper potential larval-resource reconstruction. This is an *author-anonymised review copy*, not a deposited, legally cleared dataset or a new independent experiment. The full upstream bulk HOSTS/WCVP files are not redistributed.

## Reproduce numerically without internet access

Extract the ZIP and, using Python 3.10 or later with no third-party modules, run from its root:

```bash
python RUN_ME.py
```

This validates SHA-256 checksums of every included script and input, reconstructs the original 239-butterfly native/newly added resource accounting, adds all 1,027 curated BCE/Clarke candidate accepted-plant host relations for 81 butterflies, and recalculates the source-based geographic identity turnover. Recalculated summary data and a full 81-species CSV appear in `reproduced/`.

The expected original candidate host×butterfly×WGSRPD3 region incidences are H=75,346 (A=46,387 native-overlapping; R=14,406 extra introduced support in new cells; U=14,553 unique new butterfly×region opportunity). After the partial Europe-derived host addition the respective numbers are H=103,349; A=67,740; R=20,642; U=14,967. Exactly 3,225 of the 6,445 butterfly×region cells called *new* under either host inventory overlap, for Jaccard=0.5004 in the 81 updated butterflies.

## What is included

The `inputs/` directory contains the reconstructed accepted host–plant association table and extant/non-doubtful WCVP native and contemporary plant×WGSRPD3 sidecars used by the original analyses; 239-butterfly descriptors, observed species metrics and WGSRPD3 continent lookup; the authenticated BCE candidate accepted-host-ID evidence and its additional host range sidecars; and frozen CHELSA region climate values and official TDWG WGSRPD3 geometry. `code/` contains the original analytic scripts without source-linked account metadata, and `evidence/` contains the original HOSTS-specific full MCMC draw series for the numerically qualified 1,500-km geographic sensitivity.

The main offline command reruns all descriptive host-source results, **not** the costly post-hoc MCMC simulation. The full source-specific serial output and numerical mixing diagnostics are included for audit. For an independent rerun of the *original HOSTS*, geographically restricted null, first install Shapely and use:

```bash
python code/audit_butterfly_spatial_null_longchain.py --descriptors-csv inputs/s1_resource_descriptors.csv --insect-host-csv inputs/original/insect_host_accepted.csv --native-distribution-csv inputs/original/native_extant_nondoubtful_wgsrpd3.csv --contemporary-distribution-csv inputs/original/extant_nondoubtful_wgsrpd3.csv --level1-json inputs/level1_lookup.json --unit-climate-csv inputs/unit_climate_crossfit.csv --level3-geojson inputs/level3.geojson --output-json reproduced/qualifying_original_HOSTS_spatial_mcmc.json
```

This optional command entails tens of millions of swap proposals, is post-hoc and **has not been repeated with BCE/Clarke added host links**. The included R scripts document how the original source-specific WCVP/HOSTS sidecars were prepared, but the full raw provider datasets must be separately retrieved and source permissions respected to repeat that upstream taxonomic resolution step.

## Inference, blinding and licensing ceilings

Each host×butterfly×region incidence represents **potential mapped larval-resource opportunity**, not plant biomass, actual host acceptance, larval survival, local butterfly occupancy, or a colonisation event. Host-list source Jaccard is **catalogue agreement**, not temporal change. The BCE evidence is European and partial, not a representative independent worldwide host inventory. The WCVP introduced/native labels are not dated introduction histories, and the original HOSTS climate/geography-matched null excess of approximately 0.36 percentage points was not source-validated under BCE augmentation.

All members are checked for direct author/source account markers and explicit personal data; this scanner cannot prove perfect anonymity. Full non-repository metadata (document properties, citations, supplementary materials, filenames) still require human editorial review. **Source data licences, restrictions on sharing derived plant-distribution sidecars, and journal supplementary/embargo policies must be confirmed before this bundle can be publicly deposited or submitted.** Do not treat this ZIP as authorization to redistribute raw third-party databases.

Primary data/source citations: HOSTS Lepidopteran Hostplants (DOI 10.5519/havt50xw); WCVP (*Scientific Data*, DOI 10.1038/s41597-021-00997-6); LepTraits (*Scientific Data*, DOI 10.1038/s41597-022-01473-5); Clarke (2024), *Ecology and Evolution* DOI 10.1002/ece3.10834; TDWG World Geographical Scheme and CHELSA V2.1 (DOI 10.16904/envidat.228).
