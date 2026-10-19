# Four-butterfly life-stage experiment: source provenance and exact-byte restoration (19 October 2026)

This directory and its scripts are **exploratory only**, separate from the frozen Global Ecology and Biogeography v0.2 manuscript PR #38. They reanalyze *published prior-art* observations, not new animal measurements. Sources are public CC0 Zenodo-hosted copies of Dryad deposits. All survival inferences are bound to original source-defined life stages.

## Data sources and permissions

| Butterfly/study | Original Dryad DOI | Original frozen Zenodo mirror | Required files |
| --- | --- | --- | --- |
| *Vanessa cardui*, Saldivar & Wilson-Rankin | https://doi.org/10.5061/dryad.wdbrv15v5 | https://zenodo.org/records/10689711 | `larv.plants.csv` (69,772 original bytes); `larv.switched.csv` (33,479 original bytes) |
| *Danaus erippus* and *D. gilippus*, Ferreira & Rodrigues | https://doi.org/10.5061/dryad.3bk3j9kk2 | https://zenodo.org/records/7465143 | `field_performance.txt` (2,764 original bytes); `lab_performance.txt` (1,911 original bytes) |
| *Melitaea cinxia*, Verspagen et al. | https://doi.org/10.5061/dryad.kwh70rz2m | https://zenodo.org/records/4279915 | `Survival_data.csv` (48,291 original bytes) |

All Zenodo records declare CC0, and their original MD5 hashes and exact original file lengths are preserved in `scripts/download_exact_cross_butterfly_sources_20261019.py` and checked independently in `scripts/analyze_cross_butterfly_lifestage_raw_20261019.py`.

**Integrity correction:** The connector's browser-readable **previews** normalize some original byte sequences (notably text BOM/newlines/spacing). They preserve parsed rows but FAILED original source MD5. Four such browser-previews were deliberately **removed** from this branch. **Do not treat a missing downloaded source file as inaccessible or silently use previews.** Only Melitaea `Survival_data.csv` is source-byte-identical in git. The exact other four original files are acquired at workflow/runtime from their frozen Zenodo URLs, verified byte count and MD5, then used for statistical calculations. Neither the GitHub preview nor downloaded generated JSON should be called an original source archive. Do not disable the fail-closed checksum checks.

To reproduce with internet for original bytes:
```sh
python scripts/download_exact_cross_butterfly_sources_20261019.py
python scripts/analyze_cross_butterfly_lifestage_raw_20261019.py --output /tmp/original_lifestage_receipt.json
```
Once exact originals have been fetched, the second command is entirely offline and standard-library-only. The associated GitHub Actions workflow regenerates a checksum-verified JSON *audit receipt* and uploads the receipt alone, not the four original external data files.

## Endpoints: no fabricated adult production

- *Vanessa*: 1,684 rows (1,154 plant-trial; 530 diet-switch) and **all have `status=2`**, while original publication describes four individuals reaching the pupation phase and two completing the chrysalis. The CSV contains no independently identifiable individual pupation successes; analyze event times only, not a fitted 0% or positive pupation rate. Five species codes in first dataset, three in switch dataset; different larval sources.
- *Danaus*: 40 lab assignments split two butterfly taxa×two hosts×10, and 54 early-instar field assignments with four unknown outcomes due to cattle. Report mean development time *survivors only*; native Asclepias and exotic Calotropis are **different field sites** so host identity and location are confounded in the field and bagging may alter microclimate.
- *Melitaea*: `Dead=3` = larva died; `Dead=0` = survived **to diapause**, not adult. Source design is 14 families × 4 **day** temperatures × 2 hosts × 15 larvae; 15 °C nights were used in source experiment. Four `Dead=0` rows have unexpected numeric `Time2=4` and are flagged without silently recoding.
- Do not pool early-instar field presence, lab survival, larval death time and diapause arrival as the same demographic endpoint.
- No complete same-cohort **oviposition diversion × host-specific survival × adult offspring** panel across these studies. Neither nonmonotone butterfly population output nor causal effects of introduced host geography has been demonstrated.
- The `chocho` 239-species potential host-resource geography panel has exact named links for *Danaus gilippus* (16 introduced-added units) and *Vanessa cardui* (78), not the other two species. These global counts **cannot** be scaled by the source studies' lab survival proportions to obtain realized recruitment.
