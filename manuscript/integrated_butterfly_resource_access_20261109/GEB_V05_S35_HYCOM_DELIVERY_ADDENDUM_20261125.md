# GEB v0.5 — S3.5 HYCOM reviewer-package delta (25 November 2026)

Internal Stage A, **not a reviewer-ready package**. The 24-Nov 192-file staging archive and exact delivery register remain unchanged. PR #191 paper head `6d44689362379136fc58808c295e391e39227daf` adds S3.5 and bounded directional-current text to Results 3.2 and Discussion 4.2; PR #179 and its original export are retained. PR #190 analysis source head `d8ee0228acdf0d66612a760aef76947bde5649730` is frozen for the new source-bound replay.

## Source, attribution and interpretation

- HYCOM+NCODA GOFS 3.0 GLBu0.08/expt_19.1, Naval Research Laboratory, daily 00Z 2005, 0.08-degree, 15 m, eastward/northward velocity; six dated launch cohorts, 306 fixed starts. [Official product](https://www.hycom.org/data/glbu0pt08/expt-19pt1), [NCSS catalog](https://ncss.hycom.org/thredds/catalogs/GLBu0.08/expt_19.1.html), [specific recommended acknowledgement](https://www.hycom.org/publications/acknowledgements/ocean-reanalysis-data). Official site states public release and unlimited distribution; catalog says freely available; no blanket CC0 claim. Include recommended funding/compute attribution in approved export. The original physical code uses the daily NCSS 2005 product, not the 3-hourly product or newer GOFS.
- WCVP v13 and EuDiS exact source file and hash retain original pins; GSHHG, TDWG and WPI source-derived geography inherit existing source-specific review decisions. Mixed derived tables require file-specific rights assessment; no raw source archive, ocean velocity grids or new permission grant in component.
- Preserve all six coefficients, global paired +0.397 pp/100 km, full-design interval +0.222 to +0.572; 52/306 completed endpoints, 12 sources, 104 destinations, 1,236 primary dyads; Japan exclusion +0.006 and restricted distances/region intervals span zero; orientation diagnostic and original fixed-anchor bounds retained. No completed newly within-50-km endpoint at another fixed recipient; one censored trajectory approached Corsica. No observed seed landfall, causal transport or butterfly establishment.

## Exact additional source/replay gate

Workflow commit [59912c6](https://github.com/zuizui0223/chocho/commit/59912c62b8a09c27e21cdcb1429f98444f074c3a), [CI #45542370685](https://github.com/zuizui0223/chocho/actions/runs/45542370685): **in progress at checkpoint, not PASS**.

1. Recover hash-verified original launch roster, six date-shard receipts, 306 physical endpoints, original trait, robustness and anchor results; reconcile 306 unique starts and all archived track objects.
2. Rebuild exact family trait inputs from pinned WCVP/EuDiS and archived geography/coverage inputs. Stage original scripts and derived inputs, exact file SHA-256/byte register, with all delivery statuses held internally.
3. Re-extract the component and execute original robustness and fixed-anchor programs; compare complete original result trees with fixed rtol 1e-7 and atol 1e-9 plus exact nonnumeric values. No fresh velocity integration, full upstream coastline reconstruction, new ecological hypothesis or manual document review claimed. Record failure if comparison does not pass.
4. Only after successful execution and verification, reconcile component file hashes and source families against 192-file staging package. Neutralize identifying internal provenance in a distinct reviewer candidate, preserve numerical evidence, and decide permission/attribution per file. Do not silently overwrite existing register or declare 192-file candidate covers S3.5.

Current delivery still HOLD: author-controlled code permissions, mixed-source derivative rights and historic GBIF attribution, GEB/Dryad-compatible archive route, verified anonymous reviewer link, and exact final manuscript/appendix QA remain open. No external contact, deposit, licence grant, PR merge, ready-for-review action or submission.
