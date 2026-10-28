# UKBMS: historical presence records are abundant but the GBIF version lacks a verified non-detection denominator

**28 October 2026 — public-source-only survey feasibility. Not an occupancy/range shift or host-interaction study.** Exact fixed GBIF dataset: [UKBMS occurrence registry](https://www.gbif.org/dataset/1e266c3d-92ef-4d5a-8e4a-c04742c772c3), registered as UK Butterfly Monitoring Scheme (UKBMS), CC BY 4.0. These UKBMS source-total figures are **not the 38-butterfly panel** from the separate 114-query source test; they represent all taxa included by this particular GBIF publisher dataset.

## Verified source windows

| GBIF UKBMS dataset, Great Britain, coordinates | 1970–1990 source window | 2010–2024 source window |
|---|---:|---:|
| Indexed `PRESENT` occurrence records | **595,688** | **6,293,165** |
| Indexed `ABSENT` occurrence records | **0** | **0** |
| Years represented in facet | **18** (1973–1990) | **15** (2010–2024) |
| Fixed first 12 source records with `locationID` | **12/12** | **12/12** |
| Fixed first 12 source records with `eventDate` | **12/12** | **12/12** |
| Fixed first 12 source records with `eventID` | **0/12** | **0/12** |
| Fixed first 12 source records with `samplingProtocol` | **0/12** | **0/12** |
| Fixed first 12 source records with `individualCount` | **0/12** | **0/12** |

The **entire original dataset's explicit GBIF `ABSENT` status index is zero under these filters**. The field presence checks are only 12 arbitrary first-returned records per era, not a random schema sample. `locationID` might refer to a location grid instead of a distinct transect station. The presence of `eventDate` cannot reconstruct survey visits with zero target butterflies, zero any species, standard route effort, or detection probability. The original UKBMS monitoring scheme likely collected more effort data than this particular GBIF occurrence projection exposes; this is an **information transport limitation**, not evidence that no observers ever recorded absences or an absence of the butterflies in Britain.

## Research decision

**Historical time coverage PASS; verified repeat-survey / explicit-negative occupancy information FAIL in the GBIF occurrence representation.** This is strong evidence to prioritize a license-authorized UKBMS transect visit×site×date panel instead of generating false survey absences from GBIF. The original chocho 39-species GBIF source gate also passed, but cannot by itself supply biological colonization or climate tracking. More rows of presence-only GBIF historical records are not a substitute for site×visit detection histories.

A valid historical host-mediated climate tracking test needs to obtain: (1) repeated original UKBMS survey visits at comparable sites in 1976–1990 and 2010–2024, including focal species **zero detections**, known effort/protocol and taxonomic reporting; (2) historical 0.1° [E-OBS](https://doi.org/10.24381/cds.151d3ec6) annual/weather trajectories rather than 1970–2000 single climatology; and (3) **historically dated larval host localities** and independent 1970–2000 crop/ornamental maintenance records, plus nonhost plant and land-cover controls. Whether an introduced wild host is human-planted or continuously maintained is **not** resolved by WCVP botanical range alone. Only compare climate-only against host-conditioned models using *pre-cutoff predictors* and independently held-out later surveyed cells.

[UKBMS source access and data request](https://ukbms.org/request-data) explicitly covers more detailed annual and weekly site counts; a formal request/license may be needed. [eBMS data policy](https://butterfly-monitoring.net/ebms-data%20access) likewise requires a scientific data request for transect histories. No access was requested or granted in this analysis, and no protected dataset was assumed to be available.

## Provenance

- Fixed `UKBMS` GBIF dataset registration and query protocol, saved before source outcomes: `provenance/archive/host_hindcast/UKBMS_SAMPLING_EVENT_SOURCE_GATE_20261028.md`.
- Exact original-source Actions [#40631986374](https://github.com/zuizui0223/chocho/actions/runs/40631986374) **SUCCESS**; 3 source-independent synthetic field/negative-response tests. [Derived source-only original artifact](https://github.com/zuizui0223/chocho/actions/runs/40631986374/artifacts/12792106237). ZIP SHA256 `1b2315745b36ac784f3a05406c54d06da8d4e309e19080ced7ad0b12d75a89e7`.
- Derived source-only JSON SHA256 `e6972f8c3bbca1f19cc9f7e9f8ab1fea57354de6909391a4f21d9826ce766daf`; includes only field counts and years, source request URLs and SHA256, not a single original GBIF occurrence body, site ID, observer, exact coordinate or event ID.

Keep GEB `main` and original frozen manuscript unchanged. **No climate-host predictive benefit claimed.**
