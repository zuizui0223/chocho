# Task 3 independent review — latitude-zone conditional source-distance null

## Decision

**PASS within the reviewed scope; no blocking defect found.** The new partition implements the fixed specification, preserves the old eligible cells/outcomes/distances, and produces correct exact counts, cell occupancies, conditional moments and finite Monte Carlo tail bounds. The result is suitable for delivery with its existing ecological limitations. No production code, test, source data or measured output was edited; this review report is the only written deliverable. No new random samples, remote operations, review comments or subagents were used.

Review date: 2026-10-30 UTC. Repository reviewed: `/workspace/scratch/10d2ad986b7b/chocho`.

## Scope and evidence

Read the Task 3 brief, implementation report and review diff, then inspected the actual new script and tests plus the reused exact engine and first-run products. Principal files:

- `scripts/analyze_island_latitude_null_20261030.py`
- `tests/test_island_latitude_null_20261030.py`
- `scripts/analyze_island_source_exact_null_20261030.py`
- `data/exploratory/island_recorded_sources_20261029/recorded_source_risk_cells.csv`
- `data/exploratory/island_recorded_sources_20261029/recorded_source_target_support.csv`
- `data/exploratory/island_expansion_20261029/expanded_island_panel_selection.csv`
- `data/exploratory/island_latitude_null_20261030/latitude_null_diagnostic.json`
- The fixed plan, first-run lock, exact occupancy/block CSVs, all three saved sample CSVs and post-run validation in the same output directory.

The implementer's focused test evidence was reused rather than rerunning the suite: the retained new-suite transcript reports 15 passing tests; the implementation report and post-run receipt report 34 passing old-plus-new exact tests. Independent numerical work below did not import either production analysis module or use their counting, unranking, occupancy, moment or binomial helpers.

## Specification and source preservation

| Requirement | Independent evidence | Assessment |
| --- | --- | --- |
| Fixed geographic zones | Reconstructed directly from archived representative latitude: south `< -23.5`, tropical `-23.5 <= latitude <= 23.5`, north `> 23.5`; reviewed boundary and hemisphere tests | Matches specification |
| Frozen panel and targets | Direct source reads identify exactly 96 unique panel islands and the same 28 supported targets; zones contain 4 south, 17 tropical and 7 north targets | Preserved |
| Original risk cells and labels | Joined new occupancy records to old risk rows by scenario, species key and target island; keys are unique and sets agree for all 2,385 / 2,197 / 2,385 cells | Preserved |
| Original outcomes and scores | Checked every old/new observed label, column index, original row index, longitude and binary64 score bytes | Preserved exactly |
| Correct refined row identity | Every new row index equals `3 * original_row_index + fixed_zone_index`; the target zone and latitude match the archived panel | Correct |
| Refined and original margins | Exact integer occupancy totals agree with observed refined row, collapsed old row and original island margins multiplied by the whole-fiber count | All reconcile |
| Independent uniform draws | Reviewed unchanged integer product-rank sampler; replayed every saved rank from the declared seed and checked index, scenario and seed fields | Correct fixed stream, with replacement |
| Fixed budgets and limits | First-run lock specifies 100,000 draws per scenario, seeds 2026103011/12/13 and original 100,000-state/2,000,000-transition block caps; every saved file has exactly 100,000 rows | Matches specification |
| Incomplete and fully forced outcomes | Reviewed inherited count failures and new exception/report path: incomplete execution receives `INCOMPLETE_NO_QUALIFIED_BOUND`; fully forced fibers receive explicit no-contrast status and unqualified bounds | Correct; actual three fibers are variable |
| Source and run identity | Independently recomputed every pinned input hash and every measured output hash listed by the diagnostic, plus script, plan and lock hashes | All match |

The production script hash remains `048cc2b194f1e69b539f8c412357cdc548c355ed67c13bf885f8c1ecbff777d4`; the reused exact-engine hash remains `ae37f7902910c86d3ec9173824920f110636bcd0989dc9c41249d40c37f349f0`. The fixed plan hash is `3acbcac76c4a4d554426da2d8d0a652b20fb19db841fb4e6ac50f9fcd136e277`.

## Independent exhaustive enumeration

A separate standard-library DFS enumerator was written in the review command, directly from the original risk CSV and archived latitudes. It builds each longitude-by-zone block from allowed species–target pairs, chooses each active row's allowed column subsets, updates residual island margins and explicitly visits every completed table. It does not use the production dynamic-programming recurrence or its counts. Zero-margin rows/columns contribute unique zero states. At every completed table, it accumulates exact integer occupancies and its score using the original cells.

All 36 scenario blocks were completely enumerated. The block cardinalities are independently reproduced, in longitude 0–3 and south/tropical/north order:

- Primary: `1, 1, 1, 1, 52542, 1, 1, 10, 1, 1, 10, 28`.
- Source records >=3: `1, 1, 1, 1, 52542, 1, 1, 6, 1, 1, 10, 21`.
- Later records >=2: `1, 1, 1, 1, 459, 1, 1, 1, 1, 1, 2, 15`.

Exact block occupancies were multiplied by the counts of the other blocks to obtain whole-fiber numerators. **Every one of the 6,967 published cell numerators and denominators matches independently enumerated results**, as do statuses and displayed occupancy probabilities. Direct means and population variances of the enumerated block score lists were combined using block independence; none uses the production moment helper.

| Scenario | Whole fiber | Independently recomputed observed T | Null mean T | Null variance T | Ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| Primary | 147,117,600 | 7.4325844968663395 | 7.684739024113289 | 0.005372013155875627 | 0.7771246418586739 |
| Source records >=3 | 66,202,920 | 7.721522419804050 | 7.934990710715287 | 0.006112050077094749 | 0.8077777736415962 |
| Later records >=2 | 13,770 | 7.8179783834080325 | 8.157308104167132 | 0.011026488350979736 | 0.7122475675600188 |

Differences from stored summaries are at floating rounding level and below 1e-12. All original forced positive contributions remain in the full statistic.

| Scenario | Full positive denominator | Forced positives | Variable positives | Variable allowed cells | Variable species | Variable islands |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Primary | 66 | 30 | 36 | 89 | 26 | 11 |
| Source records >=3 | 62 | 28 | 34 | 85 | 24 | 11 |
| Later records >=2 | 35 | 17 | 18 | 44 | 15 | 8 |

These denominators and support counts independently match the saved outputs. The new constraints restrict the original table space; the narrower support must not be counted as increased ecological replication.

### Complete enumeration of the smallest global fiber

The later-record >=2 scenario additionally permits complete enumeration of the Cartesian product of all blocks. All **13,770 full tables** were scored without drawing any new samples. Exactly **3 tables** satisfy the specified inclusive closer-tail condition, giving **3 / 13,770 = 0.00021786492374727668**. Their scores are at the observed minimum within the specified numerical tolerance; the minimum is 7.8179783834080325. The full-enumeration mean, 8.157308104167134, matches the saved moment result.

This is supplementary post-run verification, not an alteration or replacement of the prescribed first-run Monte Carlo output. It confirms that repeated saved ranks in this small fiber are appropriate and that the reported Monte Carlo upper bounds are compatible with its true finite-table tail probability.

## Saved samples and binomial bounds

Read all 300,000 saved rows, without running the production sampler. Every rank is in bounds and matches the locked Python random seed stream exactly; all sample indices, scenario labels, seed fields and file hashes agree. Recomputed tail counts and sample means agree exactly. There are 99,964 / 99,932 / 13,758 distinct ranks in the three 100,000-row files: repetitions are retained as required by sampling with replacement.

The one-sided binomial limits were independently solved by bisection on the binomial CDF, using log binomial coefficients and log-sum-exp summation rather than the production beta-quantile helper. The solutions agree within 1e-12.

| Scenario | Closer hits / draws | Individual 95% upper | Three-scenario simultaneous 95% upper |
| --- | ---: | ---: | ---: |
| Primary | 33 / 100,000 | 0.0004412262763574 | 0.0004755362745642 |
| Source records >=3 | 416 / 100,000 | 0.0045106833940118 | 0.0046151932885148 |
| Later records >=2 | 24 / 100,000 | 0.0003375075740103 | 0.0003678138781162 |

These bounds concern the specified fixed uniform-table probability, conditional on the observed margins and risk mask. They are not confidence limits for an ecological effect. Exact integer counts and exact uniform sampling do not make the Monte Carlo tail estimate fully enumerated; the supplementary 13,770-table audit above is explicitly separate.

## Scientific interpretation and remaining limitations

The code and results support the following calibrated account:

> 固定した熱帯・南北温帯相当の緯度区分を追加すると、記録上の距離比は主解析で0.488から0.777、感度解析で0.550から0.808、0.588から0.712へ1に近づいた。新しい条件付き比較でも近い側への偏りは残るが、比較可能な変動陽性は36・34・18セルまで減り、主解析の残差の約95.9%は90°W–0°の熱帯群に集中する。したがって、広い経度区分内の緯度の違いが元の比較に関わっていたことを示す探索結果であり、分散・定着・宿主・港の機構を確定する結果ではない。

The independently recomputed signed longitude-1 tropical contributions are:

| Scenario | Longitude-1 tropical contribution | Global observed-minus-null difference | Signed share of net difference |
| --- | ---: | ---: | ---: |
| Primary | -0.2417915626398877 | -0.2521545272469492 | 95.8902% |
| Source records >=3 | -0.2278174747186151 | -0.2134682909112362 | 106.7219% |
| Later records >=2 | -0.3361294669305218 | -0.3393297207590997 | 99.0569% |

The source>=3 share exceeds 100% because another block offsets part of the signed contribution. These values are arithmetic decompositions, not probabilities or causal proportions explained.

The analysis is exploratory after the old outcomes and concentration were known. A representative latitude zone is a coarse geographic proxy, not measured climate or habitat. Species ranges, source sampling, observer processes, mainland sources and most islands remain incompletely represented. The target roster already conditions on known later all-butterfly recording richness. Early nondetection is not biological absence, and later recording is not first arrival or establishment. Many species–island cells remain dependent through shared geography and sources; the independent units of the Monte Carlo run are complete tables, not independent ecological replicates. Island margins still absorb island-constant port/area/mainland main effects, so this result supplies no independent port effect.

The implementation report already states these limitations adequately. Preserve them beside the numerical result in the integrated handoff. No new favorable zone, cutoff or sensitivity should be selected from this decomposition.

## Coverage limits and disposition

No source geometry was remeasured, no raw occurrence record was acquired, and no old ecological model was refitted. Full global table enumeration was performed only for the 13,770-table scenario; exact block enumeration verifies all scenario counts, occupancies and moments, while the existing tested rank/table bijection and retained execution evidence support the sampling implementation. The entire project suite was not repeated; root owns integration testing. There are **no required fixes or unresolved blocking findings in Task 3**.
