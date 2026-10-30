# Exact conditional recorded-source comparison

This run completes the independent numerical follow-up specified in the frozen [implementation plan](../../../provenance/archive/host_hindcast/ISLAND_EXACT_NULL_HOST_ROSTER_PLAN_20261030.md). It evaluates the **same conditional recording comparison** as the 20261029 analysis. Every earlier scientific input, risk mask, observed table, distance, margin, chain and failed MCMC qualification remains unchanged.

The finite table counts and sampling law are exact. Distance-score moments use floating arithmetic, and the reported tail bounds retain finite Monte Carlo uncertainty. This is post-pilot numerical validation, not an untouched ecological test or evidence of first arrival, dispersal or establishment.

## Results

The statistic is the mean stored `log(1 + upper nearest-source coast gap in km)` over every positive cell. The denominator remains 66, 62 or 35, respectively, including all forced positives. The closer-source tail includes `T_null <= original observed T + 1e-12`.

| Frozen scenario | Candidate cells / positives | Observed T | Conditional null mean T | Conditional null SD | Descriptive geometric-distance ratio | Tail hits / independent draws | One-sided 95% upper bound | Simultaneous three-scenario 95% upper bound |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Primary | 2,385 / 66 | 7.432584497 | 8.149904124 | 0.106967799 | 0.488058684 | **0 / 100,000** | 0.0000299569 | 0.0000409426 |
| Baseline source records ≥3 | 2,197 / 62 | 7.721522420 | 8.319693155 | 0.104100483 | 0.549816477 | **0 / 100,000** | 0.0000299569 | 0.0000409426 |
| Later target records ≥2 | 2,385 / 35 | 7.817978383 | 8.349500000 | 0.123391767 | 0.587710020 | **1 / 100,000** | 0.0000474378 | 0.0000604679 |

The ratio is `exp(observed T − conditional null mean T)`. It is a descriptive contrast on the geometric scale of `1 + gap`, not a causal effect, odds ratio, colonization probability or confidence interval for an ecological effect.

Zero tail hits do not mean zero tail probability. The one-sided Clopper–Pearson upper bound is `-expm1(log(.05)/100000)` when no hits occur. The simultaneous family uses per-scenario alpha `.05/3`; the reported three bounds therefore have at least 95% joint coverage under the independent uniform-table sampling law. The later-record ≥2 scenario had **one actual tail hit**; its upper bound uses the corresponding binomial quantile. No exact enumerated tail probability is claimed.

## Finite table space and exact support

| Scenario | Longitude 0 count | Longitude 1 count | Longitude 2 count | Longitude 3 count | Whole labeled-table count |
|---|---:|---:|---:|---:|---:|
| Primary | 1 | 6,870,383,478 | 91,827 | 5,541,566,504 | 3,496,100,624,678,845,194,886,224 |
| Source records ≥3 | 1 | 1,403,180,478 | 18,727 | 2,678,695,685 | 70,389,053,018,969,220,551,610 |
| Later records ≥2 | 1 | 4,692 | 1 | 952,792 | 4,470,500,064 |

Fixed zero margins sharply reduce the part of the allowed mask that can vary. Exact cell occupancy counts identify 259 variable allowed cells and 65 variable positives in the primary comparison, plus one forced positive; the source ≥3 comparison has 241 variable allowed cells and 61 variable positives, plus one forced positive. The later ≥2 comparison has 110 variable allowed cells and 26 variable positives, plus nine forced positives. All nine positives in its longitude-2 block are fixed; the remaining seven allowed cells among that block's active rows and columns can never be occupied.

The primary observed-minus-null difference in T is −0.717319627. Longitude 1 contributes −0.522070679, longitude 2 contributes −0.090170923 and longitude 3 contributes −0.105078024; longitude 0 is deterministic. These are descriptive contributions to the original statistic, not additional hypothesis tests. They show that most of this contrast comes from one coarse longitude block.

## Why the draws are uniform and independent

After retaining the original labels and removing only zero-margin rows and columns, each longitude block uses a suffix recurrence `C(i, residual_columns)`. For row `i`, the recurrence sums exact Python-integer completion counts over all lexicographic subsets of allowed columns having that row's required positive count. It keeps the full labeled residual-column vector. Negative residuals, incorrect residual totals and residuals exceeding the suffix's remaining allowed capacity have zero completions.

Each row-pattern branch receives a consecutive integer interval of length equal to its child's completion count. Induction from the unique empty suffix gives a bijection between ranks `0,...,C−1` and all feasible labeled block tables. The longitude blocks form a Cartesian product. One arbitrary-integer global rank is drawn uniformly from the whole product, split into block ranks by integer division, and unranked. New ranks are drawn independently **with replacement**; repeated ranks are retained. No floating branch probabilities, MCMC transitions or configuration-model approximation are used.

The frozen budget is 100,000 complete tables per scenario, with seeds 2026103001, 2026103002 and 2026103003. Every one of the 300,000 reconstructed binary tables was checked against its full original mask and both original margins. All budgets completed. No budget, threshold or model was adapted after viewing a new result.

The fixed per-block caps were 100,000 memoized states and 2,000,000 attempted row-pattern transitions. The largest actual block required 6,401 states and 67,074 attempted transitions. Exceeding either cap raises an error and creates an incomplete receipt without qualified bounds.

Conditional sum-score moments use a centered mixture recursion weighted by exact completion counts. Cell occupancy numerators are independently computed from forward prefix multiplicities times suffix completion counts. All occupancy row and column sums reconcile exactly to the fixed margins multiplied by the fiber count. Their weighted score mean agrees with the centered-mixture mean. Block moments add, and all fixed positive contributions remain in the full statistic.

## Numerical validation

The new test file contains 19 passing tests. Exhaustive testing covers every 3×3 binary structural mask and every feasible margin fiber, including all 19,683 compatible mask/table pairs: exact completion counts match independent enumeration, and every rank produces exactly one distinct feasible table. Other checks cover forced and zero-margin fibers, infeasible constraints, both resource caps, structural-zero cycles, labeled unequal margins, product ranks, integers above floating-point precision, sampling with replacement, score-independent table identity, conditional moments, deterministic gzip, frozen source hashes, full-table reconstruction and the original forced-cell denominator.

Tests were written and observed failing before implementation. An independent reviewer reran all 19 tests and checked the actual state-cap and transition-cap failure paths. Both produced `INCOMPLETE_NO_QUALIFIED_BOUND`, with no qualified probability fields or sample files. The original script and immutable run lock were checked again at successful completion. The execution receipt preserves the command and the six stdout progress lines returned by the original run.

The independent post-run review passed every input/output/script/lock hash, all 300,000 saved sample indexes, seeds, rank bounds and tail counts, and a deterministic replay of every saved rank from its declared seed. It independently recomputed all 6,967 cell occupancy numerators and the conditional mean and variance with forward reward moments. Sixteen selected saved tables were independently unranked and checked against the stored statistic, mask and margins. They include the later-record ≥2 tail hit: sample 8,704, global rank 2,072,598,585, T = 7.817510900516905. The K = 1 bounds were separately solved from the binomial CDF by bisection and matched the reported values.

The old MCMC run's failed ESS/tail-precision gates remain historical evidence. The new bounds apply because these are independent uniform draws: the tail-hit count has a binomial distribution at the fixed sample size. An MCMC ESS estimate or a requirement for 20 tail hits is unnecessary for an exact one-sided binomial upper confidence bound. The observed tail count still limits precision; it is not converted into an exact tail probability.

## Files

- `exact_null_diagnostic.json`: full source identities, cardinalities, block contributions, moments, sample summaries, bounds and interpretation.
- `exact_null_run_lock.json`: immutable pre-run plan/input/script hashes, seeds, budgets and caps.
- `primary_rank_T_samples.csv.gz`, `source_records_ge3_rank_T_samples.csv.gz`, `late_records_ge2_rank_T_samples.csv.gz`: every independent global rank and full statistic, in draw order.
- `exact_cell_occupancy_counts.csv.gz`: every allowed cell's exact whole-fiber numerator and denominator, occupancy status, probability and observed label.
- `exact_block_counts_and_moments.csv.gz`: block cardinalities, computational work, support and descriptive score contributions.
- `exact_null_execution_receipt.json`: successful original command, exit status, captured progress output and test/review evidence.
- `exact_null_delivery_manifest.json`: file sizes and SHA256 identities for the delivered files; the manifest does not hash itself.

All CSV files use deterministic gzip with no embedded filename or timestamp. Arbitrary-size ranks and exact counts are preserved as decimal integers. The production script SHA256 is `ae37f7902910c86d3ec9173824920f110636bcd0989dc9c41249d40c37f349f0`; the plan SHA256 is `ce04b052a37ad4f07d2c81e87916dbf3addac3083a7117ebecea47545a4f99e9`; the immutable run-lock SHA256 is `1d2e7ccc251abeece13ab3123ec6af7cfc392465e1045e4ebd0d66ee7a4468f5`.

## Replay

Run tests from the repository root:

```bash
/tmp/chocho-followup-venv/bin/python -m pytest tests/test_island_source_exact_null_20261030.py -q
```

The original execution command was:

```bash
/tmp/chocho-followup-venv/bin/python scripts/analyze_island_source_exact_null_20261030.py
```

Use a new, nonexistent output directory for a replay. Existing output directories are refused, so the first-run artifacts cannot be silently replaced:

```bash
/tmp/chocho-followup-venv/bin/python scripts/analyze_island_source_exact_null_20261030.py --output-dir /tmp/chocho-island-exact-null-replay
```

## Scientific scope

The 96-island panel and its early/late outcomes were already observed before this numerical follow-up. The uniform null conditions on the existing species-by-longitude-block and target-island positive margins, structural zeros, thresholds and target roster. Early nondetection remains different from biological absence, and a later record remains different from first arrival or establishment. Mainland and most island sources are missing; habitat, climate, observation effort, taxonomy and source-sampling biases are not fully controlled by these margins. Fixed island margins also prevent estimating standalone port, area or mainland-distance effects.

The result resolves the earlier numerical mixing limitation for this particular conditional recording comparison. It does not repair the original failed ecological validation, justify a new favorable specification, or establish a dispersal mechanism.
