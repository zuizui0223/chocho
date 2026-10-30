# Latitude-zone continuation of the island recorded-source comparison

This is exploratory work on the same frozen 96-island panel, 28 target islands and previously observed butterfly outcomes. It does not replace the original longitude-only result or establish dispersal, arrival or establishment.

The primary observed/reference ratio on the geometric scale of `(1 + upper coast gap km)` changes from 0.488059 to 0.777125 after adding three fixed latitude zones to the species-by-longitude margins. The two fixed recording-threshold scenarios give 0.807778 and 0.712248. The primary remaining signed distance contrast is 95.8902% concentrated in the tropical 90°W–0° block; only 36 positive cells across 26 species and 11 islands can vary under the new margins. Latitude zones are geographic proxies, not a complete climate adjustment.

## Evidence

- `latitude_null_diagnostic.json`: complete results, support, source/code/output hashes and interpretation limits.
- `latitude_null_fixed_plan.md` and `latitude_null_run_lock.json`: first-run specification, seeds and caps.
- `independent_input_assertions.json`: exact preservation of original cells, outcomes, scores and target-island margins.
- `latitude_exact_block_counts_and_moments.csv.gz` and `latitude_exact_cell_occupancy_counts.csv.gz`: exact integer counts, moments and forced/variable support.
- The three `*_rank_T_samples.csv.gz` files: all 100,000 independent uniform product ranks and statistics per scenario. Repeated ranks are retained.
- `post_run_validation.json` and `independent_latitude_review.md`: implementation verification and independent exhaustive block enumeration, including all 6,967 occupancy numerators.
- `latitude_conditioned_distance_comparison_20261030.png` / `.svg`: reviewed comparison figure; ratios have no ecological-effect confidence intervals.

The independent review also fully enumerates the smallest global scenario: 3 of 13,770 tables lie in its inclusive closer tail. This supplementary audit does not overwrite the specified Monte Carlo result.

## Reproduction

The analysis entry point is `scripts/analyze_island_latitude_null_20261030.py`. Read its fixed-plan CLI and use a new output directory; existing first-run output is protected against overwrite. The renderer `scripts/render_island_latitude_comparison_20261030.py` only reads stored results and does not sample again. A fresh execution is a reproduction, not a new independent ecological data set.

See the [integrated result](../../../provenance/archive/host_hindcast/ISLAND_DATED_HOST_LATITUDE_RESULTS_20261030.md) for the joint island, port and host status.
