# Source-first host-leaf phenology bridge — 28 October 2026

## Main question
Can independent USA National Phenology Network (USA-NPN) repeated observations of living milkweed **Leaves** provide real information on the seasonal host-leaf access that was explicitly *assumed* in chocho's climate-only CMIP6 monarch winter-thermal scenario PR #90?

## Frozen, pre-observation source acquisition
- Native service, as documented by USA-NPN `rnpn` (usa-npn/rnpn master, `R/base-reqs.R`, `R/npn_get_species.R`, `R/npn_phenophases.R`, `R/npn_data_download.R`): `https://services.usanpn.org/npn_portal/` with GET `species/getSpecies.json`, GET `phenophases/getPhenophasesForSpecies.json`, POST `observations/getObservations.json`.
- Exactly six preselected milkweeds (do not substitute after inspecting records): *Asclepias syriaca*, *A. incarnata*, *A. tuberosa*, *A. speciosa*, *A. fascicularis*, *A. curassavica*. If a name is missing from the NPN registry, report it as taxonomic coverage missing, not field absence.
- Completed calendar years **2018–2025**. Retrieve each year in a separately bounded query, using ONLY verified registry species IDs and that species' explicitly identified `Leaves` phenophase IDs. Do not silently substitute `initial growth`, `colored leaves`, flowering, GBIF plant records or phenophase first/last observation.
- A `Leaves` **yes** denotes at least one live, unfolded leaf; **no** denotes explicitly checked no live unfolded leaf. Unknown, absent visit, incompatible protocol date and conflicting same-date reports are not negative.
- Retain site/individual ID only in process memory to de-duplicate; no coordinates, observer identity, source row body or exact sites in committed artifacts. Output only aggregate species × year × month summaries, distinct independent sites and independent site×year counts, and source request provenance.

## Scientific gate before opening observations
1. At least **six calendar years** successfully accessed with complete, time-stamped status data; fail closed on partial years.
2. To estimate a species-specific winter leaf-access *descriptive envelope*, require ≥**15 distinct site×years** with at least **two separate winter observation dates**, across ≥**3 states**, and both explicit YES and NO winter observations. Otherwise winter leaf persistence remains *UNIDENTIFIED* for that species.
3. Winter (Nov–Feb) site×year crosses calendar years; never summarize a single January–December without reassigning Nov–Dec to the next winter label. Do not interpolate across missing months. A site×winter without directly sampled months is **unknown** for all unvisited months.
4. Count recorded YES/NO as visit-conditional leaf presence, not probability a host is available at any random site or an edible biomass measure. Site-year and observer effort are nonrandom.
5. Before coupling to 840 monarch photo 100-km hexes in PR #90, demand species, climate pixel, month and receiving-region overlap. Missing overlap means no global parameter replacement. Historical NPN 2018–2025 leaves and WorldClim 1970–2000 thermal baseline mismatch in epoch; no direct raw historic model calibration.
6. An observational plant leaf is not evidence of sufficient quality, oviposition, larval survival, completed winter generation, or benefit/harm from invasive milkweeds. Native status and cultivated status are not in Nature's Notebook leaf protocol: do NOT use these observations to identify a cultivation *effect*.

## Interpretation
Successful source retrieval can support historical **host phenology constraints** and possible geographically conditional scenario masks, not forecast leaf senescence under future CMIP6. An explicit leaf-phenology response model under warming requires annual climate/site covariates, held-out-site validation, freezing/drought covariates and a transferability gate. Keep original GEB manuscript/main and established CMIP6 climate-only PR #90 unchanged.

## Attribution
Use USA-NPN required acknowledgment in any research use: “Data were provided by the USA National Phenology Network and the many participants who contribute to its Nature's Notebook program.” See https://www.usanpn.org/about/terms and https://www.usanpn.org/data/code.
