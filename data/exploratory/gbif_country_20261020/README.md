# Full-index GBIF country facets, 20 October 2026

This is an independent country-grain dataset, separate from the original
bounded WGSRPD3 occurrence sample. Read the [result and limitations](../../../docs/exploratory/GBIF_COUNTRY_MACRO_READOUT_20261020.md).

## Frozen inputs

| File | Meaning |
| --- | --- |
| `gbif_country_counts.csv` | 23,433 country/species/window/basis rows from all 924 completed queries |
| `gbif_country_facet_receipts.json.gz` | Byte-preserving compressed query receipts, complete-window decisions, exact keys, URLs, timestamps, raw-body hashes and retry history |
| `gbif_country_raw_facets.tar.gz` | All 924 original compact JSON responses; no individual occurrence records |
| `original_239_taxon_manifest.json` | Unchanged original 239-name receipt; 135 completed samples, 96 incomplete samples, eight unresolved taxa |
| `original_135_wgsrpd_counts.csv` | Unchanged original sampled WGSRPD3 aggregate for the archived baseline replay |
| `country_wgsrpd_links.csv` / `country_wgsrpd_summary.json` | Source-pinned political-country / botanical-region intersection weights, ISO decisions and exclusions |
| `country_area_atoms.csv` / `country_area_atoms_summary.json` | Exact disjoint geographic membership areas for arbitrary host/climate subset unions |

The two basis scopes **overlap**. Do not sum them. All 231 resolved source names
completed both 2010–2017 and 2018–2025 under both scopes. They represent 230
accepted keys because Oeneis bore and Oeneis polixenes share 9804158. Both names
are excluded from the ecological comparison, leaving 229 taxa. The retained
country-assigned record totals are 33,817,891 (accepted bases) and 33,621,039
(human observation/preserved specimen core). Records are neither independent
individuals nor verified wild breeding populations.

Retrieval: 2026-10-20 04:36:45–04:53:07 UTC. All queries used `limit=0`,
`facet=country`, `facetLimit=300`, coordinates present, no GBIF geospatial issue,
and occurrence status PRESENT. All raw responses and aggregates were reparsed
and hash-checked. Missing country assignments remain explicit in the receipts.
The archive is a live-index retrieval, not a GBIF DOI occurrence download.

## Sources and boundaries

GBIF occurrence API: https://api.gbif.org/v1/occurrence/search

Natural Earth 1:50m admin-0 map units are pinned to commit
`117488dc884bad03366ff727eca013e434615127` of
https://github.com/nvkelso/natural-earth-vector . The country crosswalk records
all ISO exceptions and the source-file SHA256. Source-frozen TDWG boundaries
are from the existing author archive. Only derived areas are included here.
No raw author-only HOSTS, WCVP, CHELSA or TDWG geometry is redistributed.

## Exact analysis replay

The author must separately supply the previously verified source-root directory
containing `inputs/` and `code/`. No new network fetch is needed for the frozen
country comparison.

```bash
python -m pip install -e ".[analysis,figure]"
gzip -dc data/exploratory/gbif_country_20261020/gbif_country_facet_receipts.json.gz > /tmp/chocho_country_receipts.json
python scripts/analyze_gbif_country_host_opportunity_20261020.py \
  --source-root /path/to/author-source-root \
  --original-manifest data/exploratory/gbif_country_20261020/original_239_taxon_manifest.json \
  --gbif-counts data/exploratory/gbif_country_20261020/gbif_country_counts.csv \
  --gbif-manifest /tmp/chocho_country_receipts.json \
  --crosswalk-summary data/exploratory/gbif_country_20261020/country_wgsrpd_summary.json \
  --crosswalk-links data/exploratory/gbif_country_20261020/country_wgsrpd_links.csv \
  --crosswalk-atoms data/exploratory/gbif_country_20261020/country_area_atoms.csv \
  --atoms-summary data/exploratory/gbif_country_20261020/country_area_atoms_summary.json \
  --output-dir /tmp/chocho_country_replay
```

To acquire a **new live snapshot**, use
`scripts/acquire_gbif_country_facets_20261020.py` with this original manifest,
a fresh output directory and `--requests-per-second 1`. A later live snapshot
is not byte-equivalent to this archive; existing valid checkpoints can be
resumed. The original 96-species record-body recovery remains a separate,
transport-blocked acquisition route.
