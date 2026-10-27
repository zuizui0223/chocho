# Where do invasive butterflies actually recruit? Evidence-first source cohort

**2026-10-27. Status: first verified historical event roster and original HOSTS Location feasibility audit complete; mechanism tests not yet established.** This is independent of the frozen GEB paper in `main`; no main claims, original 54.9% reconstructed potential host-range increase, or manuscript figures have been changed.

## 1. Source-proven candidate roster, by butterfly species × recipient region

The [13-entry evidence register](../../data/exploratory/invasive_butterfly_verified_entry_events_20261027.csv) stores original range description, recipient country/territory, earliest documented calendar year or date window, potential established-by year separately, accepted scientific name, host-resource context, source URL and confidence. The original 239 source cohort is checked against the exact SHA-verified author-only source in CI, NOT inferred from the 339 species S1 manifest.

| Butterfly | Original 239? | Source-native macro region | Recipient and earliest record | Date/arrival caution |
|---|:---:|---|---|---|
| Pieris rapae | Yes | Europe; prehistoric Eurasian expansion | Quebec **1860**, New Zealand **1930**, Australia/Melbourne **1929** | Australia **~1937** is frequently the time of establishment rather than first observation; native Asia origin of some subspecies separate |
| Thymelicus lineola | Yes | Palearctic | London, Ontario **1910** | Introduction occurred by/before 1910, not proven during 1910 |
| Pieris brassicae | Yes | Europe/Eurasia | Chile, Valparaiso **1970–1971** | Reported spring/summer arrival interval, subsequent southern and northern dispersal |
| Papilio demoleus | Yes | Asia/Australasia | Dominican Republic **2004**, Puerto Rico **2006** | Asia-to-Caribbean source genetically investigated; Puerto Rico may be secondary spread, not a second independent introduction |
| Danaus plexippus | Yes | North/Central America ancestral population | Hawaii **1841**, Australia **1870/71**, New Zealand **1873** | Likely or possible natural island colonization; human transport of adult/larvae NOT proven |
| Papilio xuthus | **No** | East Asia (Japan, Korea, China) | Hawaii/Oahu **1971** | Historical recorded larval Citrus use, but not in the actual 239 study dataset |
| Cacyreus marshalli | **No** | Southern Africa | Mallorca **1987–1989** | Earliest retrospectively recognized damage versus scientifically verified established population **1989** |
| Chilades pandava / Luthrodes pandava | **No** | South/Southeast Asia; exact Taiwanese source boundaries require adjudication | Guam **2005** | Real larval use of endemic **Cycas micronesica** reported: critical example of an alien butterfly colonizing a NATIVE plant |

**Source evidence (representative):**
- Pieris rapae first North American record [University of Florida](https://ask.ifas.ufl.edu/publication/IN283), NZ [Te Ara](https://teara.govt.nz/en/1966/butterfly-white), Australian first 1929 [Museums Victoria](https://collections.museumsvictoria.com.au/species/11463) versus c. 1937 establishment [South Australian Butterflies](https://sabutterflies.org.au/pier/rapae.html).
- Thymelicus lineola original Canadian spread [Canadian Entomologist](https://www.cambridge.org/core/journals/canadian-entomologist/article/abs/skipper-thymelicus-lineola-ochs-lepidoptera-hesperiidae-and-its-parasites-in-ontario/E208E8CA601627682D33C1E5C455A46C).
- Pieris brassicae introduced Chile [Ministerio del Medio Ambiente](https://simbio.mma.gob.cl/Especies/Details/1729).
- Papilio demoleus 2004 genetic source [Harvard](https://farrell.oeb.harvard.edu/publications/provenance-old-world-swallowtail-butterflies-papilio-demoleus-lepidoptera), Puerto Rico [EPPO](https://gd.eppo.int/reporting/article-1282).
- Danaus plexippus Australia [Clarke & Zalucki 2004](https://doi.org/10.1023/B:BINV.0000010120.29634.db), NZ first [Te Ara](https://teara.govt.nz/en/photograph/12431/monarch-butterflies), Hawaii [Freedman et al.](https://academic.oup.com/evolut/article/74/2/377/6726937).
- Papilio xuthus introduced to Hawaii [BAMONA](https://www.butterfliesandmoths.org/species/Papilio-xuthus), historical 1971 original report must be checked in the primary scientific note before publication.
- Cacyreus marshalli [Biological Invasions](https://link.springer.com/article/10.1007/s10530-022-02895-3), first report conflicts [EPPO](https://gd.eppo.int/reporting/article-3868).
- Chilades pandava invasion and endemic cycad use [experimental/field study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5022404/), [Cycas recruitment conservation study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3306358/).

## 2. The critical feasibility test: HOSTS Location is informative, but sparse for actual butterfly invasion sites

The original raw [HOSTS Location catalogue scan](https://github.com/zuizui0223/chocho/actions/runs/40436090416) yielded **140,485 original records; 86,826 raw rows with insect- and host-species fields; 257 nonempty Location strings**. The original `chocho` builder discarded those location strings, and the current frozen source join can only recover locations by exact source binomial name + WCVP accepted host ID. In [source run 40438202193](https://github.com/zuizui0223/chocho/actions/runs/40438202193), the original five-butterfly source-audited results are:

| Original 239 butterfly | Raw species-level original HOSTS rows | Exact accepted-pair/location rows | Location labels | Raw host records with exact documented recipient labels |
|---|---:|---:|---:|---|
| Pieris rapae | 122 | 97 | 8 | USA 5; Hawaii 3; Canada/Australia/NZ **0** |
| Thymelicus lineola | 14 | 10 | 4 | Canada/USA **0**, broad Nearctic 3 |
| Pieris brassicae | 56 | 43 | 6 | Chile **0** |
| Danaus plexippus | 126 | 81 | 11 | Hawaii 3; Australia/NZ **0** |
| Papilio demoleus | 85 | 47 | 7 | Dominican Republic/Puerto Rico **0** |

**These are recording-coverage zeroes, NOT biological interaction or host absence!** An original North American or broad taxon label may contain invasive-region links yet be too coarsely georeferenced to prove them. The 239-species butterfly study cannot straightforwardly assume HOSTS has introduced-region records for all historical invasion targets.

For the five focal butterfly species together: **403** original species rows, **278** exact accepted-pair rows. The selected geographically explicit recipient strings furnish only **11 raw records in USA/Hawaii**, with no explicitly labelled original first-invasion countries Canada, Chile, Australia/NZ, or Dominican Republic for their respective focal invaders. HOSTS `Location` should therefore be used as a positive-evidence atlas and recording-effort denominator, supplemented by published larval feeding histories, NOT as a comprehensive host-switching reconstruction.

## 3. The ecological mechanism decomposition (necessary conditions, not assumed evidence)

An observed butterfly point in an alien receiving area can be tagged by three orthogonal ecological axes:

**A. Adult/larval butterfly establishment:** explicitly independently native vs nonnative consumer; recipient first record, recurrent wild presence, reproductive evidence. A wind-blown island arrival and a transported, established butterfly are different claims.

**B. Plant-host interaction history:** documented larval host pair observed in original butterfly native range vs introduced recipient; host plant itself WCVP native vs introduced vs WCVP-unlisted in recipient; direct new larval feeding rather than inference from host range; horticultural-only resource as its own category. Missing native-area records create **NO_NATIVE_AREA_RECORD_FOUND**, never `novel_host_adopted` by deduction.

**C. Consumer climate:** observed butterfly point climate vs an observation-effort-adjusted butterfly NATIVE-area climatic baseline from the same bioclimatic product, independent of HOST-PLANT envelope. Source `chocho` climate opportunity is currently plant resource geography, not consumer physiology.

These combine into documented host-continuity/novel-host/focal-host-native status × climate-analog/novel × wild-recruitment confirmed/unknown. In particular, a *Cacyreus marshalli* on introduced *Pelargonium* and a *Chilades pandava* larva on native Guam *Cycas micronesica* instantiate sharply contrasting **qualitative mechanisms**, but their numerical contributions cannot be estimated from this five-species HOSTS Location subset alone.

## 4. Full original 257-label geography scope gate (completed)

[Source-pinned run #40438692867](https://github.com/zuizui0223/chocho/actions/runs/40438692867) passed with **257/257** source labels and **86,826/86,826** species-level HOSTS report rows conserved, and regression tests for country/subnational/realm, historical and interception names. Original label scope is:

| First-pass label scope | Source rows | Interpretation |
|---|---:|---|
| Exact named country or coded territory | **43,507** | Can be linked to a country polygon, but WCVP status may vary within it |
| Biogeographic region or multicountry | **37,550** | Cannot assign a unique national WCVP host origin |
| Exact subnational or island locality | **4,192** | Requires original botanical polygon intersection (e.g. Hawaii vs mainland USA; West Malaysia vs Sabah) |
| Historical/multinational names | 487 | Requires historical-boundary adjudication |
| Explicit interception | 276 | Not treated as evidence of wild larval use |
| Provisional location labels | 144 | Excluded from primary status |
| Unreviewed labels | 670 | Need manual resolution; never force to a country |

The fraction with an initially **single-country/territory label is 50.1%**, not the fraction safely classifiable as WCVP-native or introduced. Another **43.2%** is explicitly region/multicountry. All records remain in the source denominator; no location label was silently discarded or classified as a botanical origin. [Full source-only 257-label scope receipt](https://github.com/zuizui0223/chocho/actions/runs/40438692867/artifacts/12704809380).

## 5. Next empirical gate: WCVP plant-origin status (not yet completed)

The correct follow-up reconstructs all 257 Location labels (country vs specific botanical L3 vs historical/multiple countries vs broad biogeographic regions vs interceptions). Against WCVP native/current extant accepted plant codes, label each confidently georeferenced `butterfly × accepted host × report region` as native, introduced, unlisted, mixed or uncertain **without recoding the 0 location records as zero ecological use**.

Then build a bias-aware source-host adoption model with species×locality as the unit, accepted host-family distance, diet breadth, crop status and matching reporting effort. For historical invasive butterflies, native consumer polygons/arrival chronology must be constructed from independent evidence; use separate case studies for the 3 panel-excluded but mechanistically valuable species.
