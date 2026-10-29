# Pre-source test: BSBI Plant Atlas historical dated host geography for UKBMS butterfly hindcast

2026-10-29. Exploratory second data source for [chocho PR #97](https://github.com/zuizui0223/chocho/pull/97). This original plant-atlas work is **independent** from the source-confirmed UKBMS butterfly visit and NASA POWER meteorology. No original GEB manuscript/main changes.

## Scientific identification logic and prerequisite

The original historical 19-transect, 27-butterfly NASA POWER (1981–1995 train, 2005–2021 test) backtest failed its climate-only predictive increment gate: out-of-time site×species-equal log-loss baseline 0.410004 vs weather 0.418668 (weather WORSE by 0.008664), 1/5 held-out transect folds favorable. Thus **no valid future butterfly range forecast** can be claimed. Importantly, a failed M0+weather comparison does not prove host availability changes matter. A source-dated plant geography is needed to evaluate any host mechanism.

## Independent original dated plant source

BSBI/UKCEH **Plant Atlas 2020**, exact published 10 km hectad×date-class source [Zenodo record 10998003](https://zenodo.org/records/10998003) DOI:10.5281/zenodo.10998003, data file `PA2020_10by10kmSquareData.zip`, published 2024, original MD5 `66a8bb4305df484839587b5fe8acf1d2`. Approx 3,497 plant taxa, documented presence-only 10km geography and multiple date classes (pre1930, 1930–69, 1970–86, 1987–99, 2000–19), based on historical volunteer botanical recording.

**Do not claim plant absence from no record**; GBIF and atlas plant survey effort changes. **Do not claim human maintenance/cultivation** from introduced status or a grassland/crop landscape. Atlas presence and historical source discovery do not prove larval host use at a UKBMS route.

## Freeze source schema and temporal honesty before reading rows

1. First retrieve exact original Zenodo ZIP, verify original MD5, ZIP CRC and safe member names. Record only filenames, original CSV file header fields, row counts when bounded, and original archives' byte SHA256. Never commit original plant×hectad rows or UKBMS transect coordinates.
2. Require date-class definitions to be inspectable; confirm whether each matrix is a location×species×date-class row, wide spatial raster or other shape. No optimistic default column assumptions.
3. The UKBMS NASA historical train years **1981–1995** overlap Plant Atlas `1987–1999` class. Treating that class as a *training-era* source at 1981 would be **future information leakage**. For a strict pre-1987 plant source, use only classes ending by 1986 (but such source can be sparse and may include records later than training visits in 1981–85). A truly past-only predictive test should therefore restrict its origin to **1987 or later**, with site/visit sample recalculated from source efforts before butterfly outcomes, or use a different historical cutoff with adequate era coverage.
4. Plant Atlas `2000–19` overlaps UKBMS test years 2005–2021. It is a **post-outcome descriptive audit**, not a predictor allowed into the historical forecast. Do not use 2000–19 occurrence to select source plants, host taxa or old training sites.
5. Establish exact butterfly→larval-host botanical taxon mappings from independent pre-existing HOSTS/LepTraits source; never substitute all plants in hectad or assume each observed grass/crop is an accepted host. Author-only restricted source artefacts must remain transient and not be redistributed.
6. Before attempting M0 vs M0+hosts, require original ≥10 historically paired fully monitored UKBMS transects with georeferencing and ≥10 independently dated, botanical-taxonomically resolved host taxa yielding spatial contrast in pre-outcome historical presence. Host presence at 10km is a **landscape proxy**, not co-located plant stems or local plant abundance. Report species×site coverage + survey-effort matched **nonhost botanical placebo**; if no nonhost control, no credible interaction interpretation.
7. Any standard model built after viewing the 27-species butterfly changes is **exploratory**, even if protocol is frozen for a given source analysis. Out-of-time predictions alone do not establish causal temperature response or whether humans maintain the host.
8. Original UKBMS 19-route source may have weak spatial independence at 10km, and NASA POWER coarse reanalysis may share meteorological pixels among route replicates: audit unique plant/hectad and climate pixels, never call the resulting 19 routes 19 independent climatic environments.

## Stage decision

This stage only audits historical BSBI plant database **availability, correctness and host-relevant temporal support**. The original climate-only hindcast has failed its predictive gate; no climate tracking, host rescue, future CMIP6 or human-managed host causal effect is asserted. Keep analysis drafts stacked; GEB manuscript/main remain untouched.
