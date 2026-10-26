# Prespecified Caterpillars Count seasonal larval-use test (2026-10-26)

**Exploratory source-first test, fixed before reading season outcomes.** The source directly observes *caterpillars*, mostly unverified moth/butterfly larval Lepidoptera, **NOT identified butterfly abundance**. The frozen GEB manuscript and submission PR #70 are unaffected.

## Primary data and prior art
Source: hurlbertlab/caterpillars-on-plants, exact commit 2c5690ac7c6dba7eca506da886c49781c37617bd. Original raw survey file data/fullDataset_2025-04-18.csv, Git blob SHA 5e9f8c8c98b7e941c53a3ba32baf9be9f6d0692d (77,610,567 bytes). Mapping files data/plant_origin_status.csv (SHA fd14396dd9b83fee12355de115df20caa501ca3d), inferredPlantNames_2024-08-20.csv and officialPlantList2024-10-04.csv are pinned at the same commit.

Whitener et al. (2026), *Ecosphere* DOI 10.1002/ecs2.70615, already found nearly fourfold higher occurrence of caterpillars on native vs alien woody plants from ~68,741 original surveys in June–July, controlling botanical/geographic structure. That published native–alien contrast is **NOT a new finding**.

## Distinct, testable seasonal question
Does the *relative deficit of observed larvae on alien versus native woody plant branches* change across survey seasons after holding botanical plant family, study site, year and survey method constant? This is an observational test of *temporal realized resource use*, NOT a test of individual survival, species-specific butterfly production or insect competition.

## Fixed unit and temporal definitions
- Observational replicate is original unique field branch survey ID, not arthropod-sighting row. Collapse all rows with the same survey ID to caterpillar-positive if Group exactly "caterpillar" and Quantity > 0; otherwise an eligible surveyed branch is caterpillar-nondetection (not site absence). Abort on inconsistent plant identification/site/year/observation method or conflicting native-vs-alien classification within a survey.
- Reproduce the authors' geographical and data-quality inclusion: exact plant scientific name (with source inferred plant name only if Species=N/A and NameConfidence >=2), source native/alien origin join; exclude WetLeaves, sites "Coweeta - BB", "Coweeta - BS", "Coweeta - RK"; source Longitude >= -100; source verified plants and plant family, valid survey year 2010–2024. Do not select new sites or dates by results.
- Three predefined seasonal DOY windows: **early 91–151**, **peak 152–212** (original published June–July), **late 213–273**. Exclude outside periods rather than silently moving boundaries.
- Exact stratum: field site Name × survey Year × woody plant Family × ObservationMethod × one fixed season. Retain only strata with >=1 eligible native branch survey and >=1 eligible alien branch survey; all count denominators are unique survey IDs.
- Within each season estimate Mantel–Haenszel pooled **alien-versus-native odds ratio of any observed caterpillar** from stratum 2×2 tables, OR = sum_i(a_i*d_i/n_i) / sum_i(b_i*c_i/n_i), where a/b are alien +/−, c/d native +/−. Never add arbitrary pseudo-counts; if MH terms vanish mark NOT_IDENTIFIABLE.
- Primary between-season change = log(OR_late)−log(OR_peak); **positive** hypothesized temporal relative compensation on alien hosts. Secondary early−peak. This was fixed without looking at outcomes.
- 2,000 **site-cluster** bootstrap replicates, seed 20261026, each time resampling site names with replacement and retaining all that site's strata. Report 2.5/50/97.5% percentiles, skipped degenerates; do not treat branches as independent biological site replicates. Visual-only source-method subpanel is predeclared sensitivity if all three seasons meet source gates.
- Prespecified source gate **each season**: >=30 matched site×year×family×method strata, >=10 sites, >=200 surveyed branches for native AND alien, and >=10 caterpillar-positive surveys total. All must pass before interpreting late vs peak and early vs peak. Confirm original authors' 68,741 unique June–July surveys under their exact quality filters or explicitly explain source-version discrepancy.

## Interpretation/stop
Even a positive late−peak contrast cannot establish that alien hosts are phenological refugia, causal environmental rescue, or that any butterfly species increased. Plant species identity changes across regions and time, plant age/leaf biomass are not measured, and most caterpillars are moths. No individual egg-to-adult recruitment, competition or change after plant introduction can be inferred. The source code already analyses native–alien differences and latitude, so this test must be described as **secondary exploratory seasonal structure** rather than globally novel conclusions. Failure of source or bootstrap gates is not ecological evidence of no use.

Source and derived aggregate only on analysis/larval-host-seasonal-realization-20261026; main GEB paper isolated.
