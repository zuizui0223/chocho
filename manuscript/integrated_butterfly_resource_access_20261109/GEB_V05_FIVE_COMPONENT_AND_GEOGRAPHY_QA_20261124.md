# GEB v0.5 — five-component internal archive and geography audit, 24 November 2026

**Both PRs remain Draft.** [Run 45355093236](https://github.com/zuizui0223/chocho/actions/runs/45355093236), workflow head `412d71220b071feab939ba71176f060ac8a5f018`. Both jobs completed successfully. This is a limited internal replay, not a reviewer-ready archive.

## Integrated replay

Recovered and hash-verified the existing 39-member archive and the separate S3 EuDiS artifact. Integrated 7 EuDiS component members and 2 source receipts. Preserved 36/39 original members byte-for-byte; only the three wrappers changed (manifest, README, runner). Built and independently extracted a new 48-member ZIP. All five component replay checks passed, including exact EuDiS paired CSV and JSON output comparison. The existing archived scores and biological numbers remain unchanged.

| Component | Files | Scope |
|---|---:|---|
| Fig1/Fig3B–C | 7 | Aggregate resource/geographic replay; Fig3A separately verified |
| Fig2A–C | 11 | Derived WCVP/WPI and marine regression inputs |
| Fig2D | 8 | Fixed239/MEM10 derived-input model |
| Fig4 | 5 | Archived scores only, not fresh 155/224 refits |
| S3 EuDiS | 7 | Archived regional counts and paired-region bootstrap |
| Source receipts | 7 | Includes exact v5 metadata and historical 155 receipt |
| Wrappers | 3 | Internal scope, runner and manifest |
| Total | 48 | All five extracted replays passed |

Inner ZIP: 591,276 bytes; SHA-256 `1f0bd3cbb439b2828c3b1c9415f2ab10c55f3aec1fefb336396a98a4b6472566`.
Outer Actions artifact [14901913291](https://github.com/zuizui0223/chocho/actions/runs/45355093236/artifacts/14901913291): 593,921 bytes; SHA-256 `fe0d8782ceea27396e13eeb201e5cb6a3c21a9aff675e1beb4cfecf544446b46`.
Direct identity pattern scan passed, but indirect anonymity, archive-wide redistribution rights, actual reviewer access, permanent DOI, upstream raw-data and whole-SI replay remain open. Do not send this internal artifact to the journal. The old 155 missing-panel receipt is explicitly historical and superseded; original source recovery and mode-specific reproduction outcomes are separate.

## Geography provenance

The frozen NOAA-recovered 3,630-point WPI GeoJSON matches the original marine source-bundle manifest. For all 1,596 complete plant-query ports, normalized ID/coordinates match the original private source hashes at 6 **and 7 decimal places**. The whole NOAA recovery GeoJSON (SHA-256 `62159ee5ce4a31aff90ee46efdb3cdf11a9b368472a147392eb742148f340d83`) is not byte-identical to the original private GeoJSON. Current NOAA metadata explicitly attributes NGA 2019 World Port Index. Historical original conversion provenance and original-file reuse remain open; no private coordinates newly redistributed.

The preserved GSHHG source notice, LGPL text and README matched the frozen manifest; original licence notice states LGPL version 3 or later. Apply notice conditions and source-specific review for redistribution. Exact source notices and checksums in [rights matrix](GEB_V05_SOURCE_RIGHTS_AND_REUSE_MATRIX_20261121.md). [Geography receipts artifact 14902500472](https://github.com/zuizui0223/chocho/actions/runs/45355093236/artifacts/14902500472).

Next: complete current SI replay inventory and permitted data/code delivery; resolve remaining historic GBIF attribution or explicitly documented limits; do not fabricate historical DOI. Original manuscript PR #179 remains unchanged. No journal submission, no coauthor request, no public deposit.
