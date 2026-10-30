# Fixed exploratory comparison: recorded hosts and later butterfly occurrences (2026-10-30)

**Status:** exploratory diagnostic, not preregistered, not a test of biological host absence. Earlier work already revealed the 20/66 recorded-host overlap and 11/36 variable-positive overlap. This specification is fixed **before computing the exact recorded-host table-score distribution**. No change to the already observed outcomes is allowed.

## Scientific question

Given the already selected and documented 96-island recorded-source comparison (28 targets, 93 early-recorded butterflies, 2,385 eligible butterfly–island cells and 66 later-positive cells), are positive 2018–2025 butterfly-recording cells unusually enriched for **documented prior known-host plant recording** in the same island?

The exposure is strictly `confirmed_host_record == True`: one or more conservatively matched host taxa from the frozen original literature-derived host list has a 2010–2017 indexed species occurrence in the exact target-island polygon. All other cells are **not documented early host co-recording**, which *does not* mean absence of biological hosts or unusable resources. Explicit unresolved taxonomic states remain intact.

## Locked finite-law score and one-sided direction

Let `Y_si` be the already observed binary later-positive butterfly recording, and `R_si` the frozen recorded-known-host flag. Define `H = sum(Y_si * R_si)` over all 2,385 originally allowed cells, including forced-positive cells. Larger H means stronger recorded-host co-location, not increased population abundance.

Primary mathematical reference law: uniform over allowed 0/1 tables preserving (a) original species-by-longitude-by-latitude-zone positive totals; (b) original target-island positive totals; (c) all allowed structural cells from the prior frozen source-distance analysis. This is **exactly the previous primary latitude-conditioned fiber**, with reported cardinality 147,117,600. Under this law, island-only recording-effort totals are fixed by column margins, species-by-zone occurrence totals by row margins; neither is estimated as a coefficient.

Primary report: H_obs, exact E[H], exact support [min,max], exact count of all allowed full tables with H >= H_obs, and inclusive upper-tail probability. Decompose the count into forced H and variable H. Use integer-count enumeration by fixed blocks, convolving blockwise score-count distributions; no Monte Carlo. Check independent sum of multiplicities equals 147,117,600 and per-block counts agree with existing diagnostic. Use all cells, *not* a new selected 89-cell-only dataset (the latter is a display of existing conditional support).

## Interpretation guardrail

A small tail probability would be evidence only of nonrandom **catalogue-compatible botanical and butterfly co-recording** relative to this fixed selected recording panel. It would not show actual larval utilization, establishment, planted/native host status, introduction mechanism, a port effect, causation, or a test with genuinely host-free cells. The baseline panels have strong, biased geographic and effort coverage; 72/89 variable no-match cells retain unresolved host catalogue entries. Those are not valid host-absent controls.

The separately documented source-distance association is strong; a first-pass H association may reflect source distance. A prespecified descriptive companion will calculate exact finite-fiber covariance/correlation between H and the old `sum(log1p(source_gap_km) * Y_si)` score **only if it can be computed with the same integer-count block enumeration**, without adjusting the primary probability. If the two indicators covary, explicitly mark biological independence unestablished. No new port coefficient, time-window changes, taxonomy promotion, null reruns or newly queried islands.

## Verification conditions

- No modified old source or model inputs, no relaxed taxonomy, no new p-value test selection.
- Verify all 2,385 butterfly–island keys and exactly 66 later positives; every target island present in the panel.
- Reproduce 12 existing latitude-zone blocks (including 2 empty) and each block's prior exact table count; full product 147,117,600.
- Verify the forced/variable support counts: 30/36 observed positives, and recorded-host among them 9/11, as previously reported.
- Enumerate score frequencies using independent per-row subset recursion and column capacities, with BigInt/integer accumulation; compare to a tiny brute-force synthetic enumerator in tests.
- Report the actual output as an additional exploratory follow-up in draft PR #102; do not promote this result into the independent GEB manuscript.

## Completion criteria

A reproducible source script, fixed stdout/JSON results, test, result memo, and an explicit label of what remains scientifically unidentifiable. If the exact count is incompatible with the previous fiber, stop and report a failed verification rather than inventing a p-value.
