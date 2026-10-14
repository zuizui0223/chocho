#!/usr/bin/env python3
"""Conservative precision SCREEN for a not-yet-run Kyoto butterfly experiment.

Uses the *assumed* between-INDEPENDENT-BLOCK standard deviation of a FOUR-ARM
invader-by-food contrast. Both four- and six-arm balanced designs estimate
this same contrast. For B >= 10 independent source blocks, t_{.975,9}=2.2622
upper-bounds the Student critical value for higher B. The resulting number of
blocks achieves the requested approximate two-sided 95% half-width assuming
stationary, independent, fully observed block estimates.

This is NOT sample size recommendation or power; it ignores uncertainty in the
hypothetical standard deviation, finite block heterogeneity, selection, maternal
dependence, censoring, resource noncompliance and food quality. Use pilot data
and design-specific simulation prior to animal allocation. 10 is a display
floor for this illustration, not an adequacy threshold.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

T_DF9_UPPER_BOUND = 2.2621571627409915
MIN_DISPLAY_BLOCKS = 10


def scenario(sd_block_contrast: float, target_halfwidth: float) -> dict:
    if (not math.isfinite(sd_block_contrast) or sd_block_contrast <= 0
            or not math.isfinite(target_halfwidth) or target_halfwidth <= 0):
        raise ValueError("Block SD and half-width must be positive finite numbers")
    raw = (T_DF9_UPPER_BOUND * sd_block_contrast / target_halfwidth) ** 2
    if not math.isfinite(raw) or raw > 100_000:
        raise ValueError("Requested precision is unsupported or implausibly high")
    blocks = max(MIN_DISPLAY_BLOCKS, math.ceil(raw - 1e-12))
    achieved_conservative_halfwidth = (
        T_DF9_UPPER_BOUND * sd_block_contrast / math.sqrt(blocks)
    )
    if achieved_conservative_halfwidth > target_halfwidth + 1e-12:
        raise RuntimeError("Illustrative precision inequality violated")
    return {
        "assumed_independent_block_SD_risk_difference": sd_block_contrast,
        "requested_two_sided_halfwidth_risk_difference": target_halfwidth,
        "independent_complete_blocks_illustrative": blocks,
        "minimum_whole_cages_for_primary_FOUR_arm": blocks * 4,
        "minimum_whole_cages_for_extended_SIX_arm": blocks * 6,
        "additional_cages_to_estimate_identity_arm_comparisons": blocks * 2,
        "conservative_critical_value_t_df9_for_B_gte10": T_DF9_UPPER_BOUND,
        "calculated_conservative_halfwidth": achieved_conservative_halfwidth,
        "not_power_or_empirical_variance": True,
    }


def grid(sd_values: list[float], halfwidth_values: list[float]) -> dict:
    if not sd_values or not halfwidth_values:
        raise ValueError("The hypothetical input grid cannot be empty")
    return {
        "schema": "kyoto_six_arm_experimental_precision_sensitivity_v01",
        "status": "HYPOTHETICAL_PILOT_PLANNING_ONLY_ZERO_OBSERVED_VARIANCE",
        "primary_estimand": "S_competitor_x_food",
        "assumptions": [
            "Source blocks truly independent; each block contains one whole cage per assigned arm",
            "Approximate Student t interval for means of independent primary block risk differences",
            "Conservative t(9) critical value 2.2622 applied for block counts >=10",
            "No clustering across blocks due to shared maternal families/donor plants or site",
            "Six arms add two cages per block to primary four-arm test for identity contrasts",
            "No unobserved initial larvae; followup complete and treatment fidelity achieved",
            "This is NOT inferential power and cannot justify animal sample sizes by itself",
        ],
        "scenarios": [
            scenario(sd, target) for sd in sd_values for target in halfwidth_values
        ],
        "authorizations_or_animal_work_completed": False,
        "main_GEB_unchanged": True,
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--block-sds", nargs="+", type=float, default=[0.10, 0.20, 0.30, 0.40])
    p.add_argument("--halfwidths", nargs="+", type=float, default=[0.10, 0.15, 0.20])
    p.add_argument("--receipt", required=True, type=Path)
    args = p.parse_args()
    out = grid(args.block_sds, args.halfwidths)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
