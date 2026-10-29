# Baseline-only census of120 previously unseen islands

All120 original-source islands in the fixed near-continent/far-WPI deficient cell were queried for2010–2017 only. The original96-island roster is excluded. Selection requires100<=area<20,000km², representative latitude within±55°, no Antarctic ambiguity, continental-coast-gap upper<100km and WPI-distance lower>100km.

The acquisition completed120/120 queries, with114 complete zeros and6 positive islands containing44 indexed records. The6 counts are1,1,2,10,13 and17. **Zero islands pass the unchanged baseline recording screen** (at least20 species-assigned records, at least90% assignment and at least3 recorded years). The predeclared stop rule was reached, so no2018–2025 source query was made for this roster. Indexed zero is not butterfly absence or a surveyed nondetection.

The plan, source/selection identities and all120 exact geometry hashes were frozen before queries. Selection SHA256: `d96226d1ad1ee153f8ccd4a935026aad3f54d0929c86abcce8cedde50c46d08e`. Full acquisition ran2026-10-29T13:52:41.608287Z–13:58:43.990708Z, after selection lock at13:51:47.449149Z. There were240 actual HTTP calls: a complete facet response and separate first-page check for every island.

- `unseen_island_baseline_selection.json` and `.csv`: complete geographic roster and source/geometry identities.
- `unseen_island_baseline_manifest.json`: query contract, complete facet responses, source-byte receipt hashes, statuses and stopping decision.
- `unseen_island_baseline_counts.csv`: all120 outcomes and eligibility/exclusion reasons.
- Four `*_facets.csv` files: species, year, dataset and basis margins.
- `continuation_verification.json`: cross-task independent review, source/numerical checks, tests, figure review and original-bundle identity.

The full interpretation is in `provenance/archive/host_hindcast/ISLAND_BASELINE_HOST_RESULTS_20261029.md`. Preserve the first acquisition. A replay with a fresh source cache is a new live-source snapshot:

```bash
python scripts/audit_unseen_island_baseline_20261029.py \
  --source-cache /tmp/chocho-island-sources \
  --cache-dir /tmp/chocho-baseline-replay-cache \
  --output-dir /tmp/chocho-baseline-replay
```
