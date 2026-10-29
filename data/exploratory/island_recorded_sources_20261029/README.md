# Prior recorded-source proximity: post-pilot diagnostic

The fixed 96-island early/late outcomes were already observed. This is an exploratory conditional recording diagnostic, not a temporal holdout or a colonization test.

The frozen plan is `provenance/archive/host_hindcast/ISLAND_BASELINE_HOST_PLAN_20261029.md`. The original richness/port plan and failed gates are unchanged.

| Analysis | Candidate cells | Later-record positives | Observed T | Null mean T | Geometric-distance ratio | Qualified approximate tail |
|---|---:|---:|---:|---:|---:|---:|
| primary | 2385 | 66 | 7.43258 | 8.14948 | 0.488267 | not qualified |
| source_records_ge3 | 2197 | 62 | 7.72152 | 8.32173 | 0.548696 | not qualified |
| late_records_ge2 | 2385 | 35 | 7.81798 | 8.35105 | 0.586798 | not qualified |

T is mean log(1 + nearest sampled-source coast gap in km) among later-record positive cells. Source candidates had earlier focal records on other members of the sampled96. Upper coast-gap estimates supply T; the JSON retains sampling intervals and lower-bound observed T.

The null fixes species-by-longitude-block and target-island positive margins. Four fixed chains each use100,000 attempted burn-in steps and4,000 samples every50attempts. Rejections and0.5lazy holds count as attempts. Chains are never extended after observing a result. Full diagnostics, qualification failures and finite tail counts are in the JSON; samples, margins and the explicit allowed risk cells are supplied for independent review.

Three other recorded taxa do not establish visits or adequate detectability. Later-record outcomes are not arrival or establishment. The source pool omits mainland and many island populations. Fixed island margins prevent estimating standalone port, area or mainland-distance effects. The distance ratio is descriptive, not a causal effect or an odds ratio.

Replay from repository root:

```bash
/tmp/chocho-geo-venv/bin/python scripts/analyze_island_recorded_sources_20261029.py --source-cache /tmp/chocho-island-sources --output-dir /tmp/chocho-recorded-source-replay
```

Use a fresh replay directory to preserve the first run lock and exact original samples. The rendered PNG/SVG summarize those preserved samples without rerunning the sampler; the main continuation report and verification receipt cover these later presentation artifacts.
