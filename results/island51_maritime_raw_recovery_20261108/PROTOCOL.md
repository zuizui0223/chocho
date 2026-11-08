# Archived GLSN2018 network-edge source and UNCTAD PLSCIQ time-series recovery

2026-11-08. Exploratory **data recovery and independent port-identity audit only**. Original 51 island/landmass groups, 194 original WPI2019 botanical ports, 182 marine-qualified ports, both butterfly observation periods, prior GEB larval-host resource results and PR137/140/141 marine-source limitations are unchanged. No butterfly prediction is fitted; no manuscript is submitted.

## 2018 GLOBAL LINER SHIPPING NETWORK (GLSN)

Original open MIT-licensed research repo: https://github.com/Network-Maritime-Complexity/GLSN-efficiency-vulnerability . Read-only pinned two **original binary XLSX** files:
- `data/port 2018.xlsx` — Git blob SHA-1 `5486d8c204a7c49b375a14762fe5cb36c7877d96`; size 34675 bytes.
- `data/Edge 2018_with great circle distance (GC) and weight.xlsx` — Git blob SHA-1 `18098cab432ddbcd7e0a7091fc1c4cedbda40fb7`; size 1018504 bytes.

Retrieve from the publisher repository in GitHub Actions (not the connector's Unicode-broken XLSX parser), and reject bytes unless Git blob-hash checks match. Use stdlib ZIP/XML for workbook extraction, record exact sheet names, dimensions and column headers, sheet row counts, numerical edge weights and source WPI identity metadata. Earlier paper/public summary cites 887 network ports and 14052 undirected weighted edges: verify rather than assume these counts. A link indicates co-membership in a liner service network and is NOT an independently confirmed nonstop physical ship voyage. **No imaginary OD traffic or biological colonization effect.**

Join original WPI194 to GLSN port data only where source port label, country (or confidently mapped source country), and original centroid agree and a unique candidate survives. If GLSN original port sheet has no authoritative coordinates/ISO country fields, freeze all WPI joins as NAME-ONLY UNVERIFIED; do not silently use published summary-only names to claim geographic identity. Do not interpret edge absence or missing WPI port as evidence of no liner service. Preserve all source-only matched-edge outcomes separately from uncertain candidates. Do not republish raw third-party network XLSX in the public chocho Git tree; source blob SHA, sheet counts, permitted derived WPI identity/edge metadata and original source URL are sufficient. Original source may remain in a GitHub Actions artifact (time-limited) where licensed.

## UNCTAD/DBnomics port-level PLSCIQ 2006Q1–2021Q3

The project already records an original public DBnomics retrieval from 2026-10-30, provider edition updated 2021-09-15 and indexed 2022-04-04. Archived endpoint: `https://api.db.nomics.world/v22/series/UNCTAD/PLSCIQ`, at `offset=0,500,1000`, `limit=500,observations=1`, with original response SHA256:
- offset0 `29a343bb6b1113b1699a9ca1821952cce1cb8f150831921da68e248e78fb3e1b`
- offset500 `c82e2d2119b0f4b03e4afb18c04f9c4442f19a2bfce96efbc8d30a87aed42933`
- offset1000 `7eb1eaf4349008b784d9e6d2bd834cadf29c51fa493c2cdd2f83e2d767978255`.

Old audit: 1300 series, 81900 port×quarter rows, 57312 nonmissing, 2006-Q1–2021-Q3, base `Index (Maximum Q1 2006=100)`, directory hash `5cfa0ae00d0c42645dae30be6a5fad068dc35c5f`. Rerequest the same public API responses with bounded timeout, retries and explicit checks. If SHA matches exactly, mark historical raw source RECOVERED. If hashes differ, never call the new file old-version recovered: save new upstream hashes/version and drift summary instead. Validate all 3 pages, each series code's uniqueness, source-label country+port and per-period values, missing != zero. Source raw downloaded files should not be redistributed publicly until rights are verified. Archive source-only summary of 2006–2010 index availability and mean/trend among **original 194 WPI port-name/country candidates**, not copy all 1300 full provider time series into the public PR.

The **PLSCIQ is a port-level liner connection intensity index**, not count of ship voyages, cargo throughput or butterfly arrivals; it contains no bilateral route direction and is not a replacement for the 2020 Zenodo physical OD source.

## Scientific stop and independent check

Publish exact Git blob hashes, old/new DBnomics page SHA256, source vintage, source/API status, source source header/sheets, counts, 194 original source-port join ledger with match confidence, and a standalone consistency auditor. Old source scripts/reports at `provenance/reviewer_defenses/results/port_connectivity_archive_20261030/` are not overwritten. User authorized GitHub commits and PR edits, **not journal submission**.

Any WPI match requiring extra validation stays unverified, and no butterfly models are fitted until usable source coverage is established. Distinguish successful source acquisition from evidence of colonization.
