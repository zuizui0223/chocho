# Independent CC0 butterfly larval observations — exploratory only (19 October 2026)

Raw files from published Dryad datasets mirrored by Zenodo, fetched from each archived downloadable text/CSV link on 19 October 2026. Original filenames and raw ASCII bytes are retained *within taxon subdirectories*. All five file sizes and their Zenodo original MD5 are checked by `scripts/analyze_cross_butterfly_lifestage_raw_20261019.py`, which rejects any byte drift.

| Subfolder / files | Source/DOI | Zenodo version | Checksum in script |
| --- | --- | --- | --- |
| `vanessa/larv.plants.csv`; `vanessa/larv.switched.csv` | Saldivar & Wilson-Rankin (2024), https://doi.org/10.5061/dryad.wdbrv15v5 ; associated article https://doi.org/10.1002/ecs2.4810 | https://zenodo.org/records/10689711 | CC0 |
| `danaus/field_performance.txt`; `danaus/lab_performance.txt` | Ferreira & Rodrigues (2022), https://doi.org/10.5061/dryad.3bk3j9kk2 ; associated article https://doi.org/10.1002/ece3.7821 | https://zenodo.org/records/7465143 | CC0 |
| `melitaea/Survival_data.csv` | Verspagen et al. (2020), https://doi.org/10.5061/dryad.kwh70rz2m ; associated article https://doi.org/10.1098/rspb.2020.2577 | https://zenodo.org/records/4279915 | CC0 |

Original READMEs are accessible at their published records and were manually cross-checked before defining stages. Copying files here does not imply the authors were part of this analysis or that their original experiments are newly observed results.

## Critical endpoint boundaries

- *Vanessa cardui*: 1,684 individual records, all `status=2` (README calls this death), but original article reports 4 reaching pupation / 2 forming pupae. The raw fate code does not independently identify pupation/adult emergence. Use recorded longevity only, and never estimate pupation from the CSV's all-dead code.
- *Danaus erippus*, *D. gilippus*: 40 lab larval assignments with survival and developmental time; developmental means use **survivors only**. The 54 *D. erippus* field entries are early-instar presence/survival under bagging vs open exposure. Native Asclepias and exotic Calotropis were in separate sites, so no causal cross-host field effect; 4 grazing-related NA must be kept.
- *Melitaea cinxia*: 1,680 larval assignments in exactly 14 maternal families × 4 day temperatures (28,30,32,34°C) × 2 host plants × 15 larvae, with death (`Dead=3`) versus *reaching diapause* (`Dead=0`). Not adult emergence. Four `Dead=0` rows have non-NA `Time2=4`, contrary to the README's interval-time convention. The original `Dead` column remains unchanged and exceptions are exposed in the result.
- The ~239-butterfly geographic HOSTS/WCVP estimate is independent of these original experimental outcomes; the latter are **not** a global adult recruitment or land-cover/plant-introduction time series.

Run `python scripts/analyze_cross_butterfly_lifestage_raw_20261019.py --output /tmp/larval_audit.json` offline, without libraries other than Python standard library. The exploration-only workflow checks source MD5/row counts, non-fabrication, fatal value ranges and writes a machine-readable artifact.

This directory is ONLY for the isolated exploratory branch; it is not the frozen GEB manuscript data-and-code supplementary review bundle. If the branch is ever merged into the v0.2 manuscript branch, the anonymous review staging EXCLUDE_RELATIVE list must first be updated and CI must confirm these exploratory files cannot leak into the original double-anonymous submission.
