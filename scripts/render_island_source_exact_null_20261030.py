#!/usr/bin/env python3
"""Render the saved exact-null result without changing or rerunning the analysis.

The plotted contrasts are descriptive. Monte Carlo tail-probability bounds are
deliberately not drawn as uncertainty intervals around these distance contrasts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter


REPO = Path(__file__).resolve().parents[1]
DEFAULT_DIRECTORY = REPO / "data/exploratory/island_source_exact_null_20261030"
DIAGNOSTIC_SHA256 = "9b2b21478c6a73dce547328995628f58e7d3c0f351649a037c7fee9fa74c939c"
SCENARIOS = ("primary", "source_records_ge3", "late_records_ge2")
EXPECTED_RATIOS = (0.4880586837411645, 0.5498164765760946, 0.5877100195080515)
BLOCK_LABELS = ("−180° to −90°", "−90° to 0°", "0° to 90°", "90° to 180°")


def reviewed_values(path: Path) -> tuple[list[float], list[float], float, float]:
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == DIAGNOSTIC_SHA256, "Frozen diagnostic identity differs"
    diagnostic = json.loads(raw)
    assert diagnostic["qualified_bounds"] is True
    analyses = diagnostic["analyses"]
    ratios = []
    for name, expected, positives in zip(SCENARIOS, EXPECTED_RATIOS, (66, 62, 35)):
        analysis = analyses[name]
        ratio = analysis["observed_to_null_geometric_mean_one_plus_gap_ratio"]
        assert math.isclose(ratio, expected, rel_tol=0, abs_tol=5e-13)
        assert math.isclose(ratio, math.exp(analysis["observed_T"] - analysis["null_mean_T"]),
                            rel_tol=0, abs_tol=5e-13)
        assert analysis["support"]["positive_cells"] == positives
        ratios.append(ratio)
    primary = analyses["primary"]
    blocks = primary["blocks"]
    assert [block["longitude_block"] for block in blocks] == [f"longitude_{i}" for i in range(4)]
    contributions = [block["observed_minus_null_contribution_to_full_T"] for block in blocks]
    total_difference = primary["observed_T"] - primary["null_mean_T"]
    assert math.isclose(math.fsum(contributions), total_difference, rel_tol=0, abs_tol=1e-12)
    assert contributions[0] == 0.0
    western_share = contributions[1] / total_difference
    assert math.isclose(western_share, 0.7278076048233534, rel_tol=0, abs_tol=1e-12)
    return ratios, contributions, total_difference, western_share


def render(source: Path, output_stem: Path) -> dict:
    ratios, contributions, total_difference, western_share = reviewed_values(source)
    navy, muted_blue, ink, muted, grid = "#225E8A", "#9BAFBE", "#172B3A", "#526674", "#E3E9ED"
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "text.color": ink, "axes.labelcolor": ink,
        "xtick.color": muted, "ytick.color": ink,
        "axes.edgecolor": "#8796A1", "axes.linewidth": 0.8,
        "axes.unicode_minus": True, "savefig.facecolor": "white",
        "svg.fonttype": "none", "svg.hashsalt": "chocho-exact-null-20261030",
    })
    fig = plt.figure(figsize=(13.2, 7.2), facecolor="white")
    left = fig.add_axes((0.185, 0.345, 0.295, 0.405))
    right = fig.add_axes((0.675, 0.345, 0.285, 0.405))

    fig.text(0.045, 0.935, "Exact conditional recorded-source comparison",
             fontsize=21, weight="bold", va="top")
    fig.text(0.045, 0.881, "Frozen candidates, distances and margins · three fixed scenarios",
             fontsize=11.5, color=muted, va="top")
    fig.text(0.045, 0.812, "A  Distance contrast", fontsize=14, weight="bold")
    fig.text(0.565, 0.812, "B  Primary contrast by longitude", fontsize=14, weight="bold")

    left.set_axisbelow(True)
    left.grid(axis="x", color=grid, linewidth=0.8)
    left.axvline(1.0, color=muted, linewidth=1.2, linestyle=(0, (4, 3)))
    positions = (2, 1, 0)
    left.scatter(ratios, positions, s=90, color=navy, zorder=3)
    for x, y in zip(ratios, positions):
        left.annotate(f"{x:.3f}", (x, y), xytext=(10, 0), textcoords="offset points",
                      ha="left", va="center", fontsize=12, weight="bold", color=navy)
    left.set_yticks(positions, ("Primary", "Source records ≥3", "Later records ≥2"))
    left.tick_params(axis="y", length=0, pad=10)
    left.set_xlim(0.0, 1.12)
    left.set_ylim(-0.55, 2.75)
    left.set_xticks((0, 0.25, 0.5, 0.75, 1.0), ("0", "0.25", "0.50", "0.75", "1.00"))
    left.text(1.0, 2.57, "Null = 1", ha="center", va="center", fontsize=10, color=muted,
              bbox={"facecolor": "white", "edgecolor": "none", "pad": 2})
    left.set_xlabel(r"$\exp(T_{\mathrm{obs}}-\mathrm{E}[T_{\mathrm{null}}])$" +
                    "\nRatio on geometric scale of (1 + gap)", labelpad=13)

    right.set_axisbelow(True)
    right.grid(axis="x", color=grid, linewidth=0.8)
    right.axvline(0, color=muted, linewidth=1.2)
    right_positions = (3, 2, 1, 0)
    right.barh(right_positions, contributions, height=0.47,
               color=(muted_blue, navy, muted_blue, muted_blue), edgecolor="none", zorder=3)
    right.plot(0, 3, marker="o", markersize=4.5, color=muted, zorder=4)
    for x, y in zip(contributions, right_positions):
        right.annotate(f"{x:.3f}".replace("-", "−"), (x, y),
                       xytext=(-7, 0), textcoords="offset points", ha="right", va="center",
                       fontsize=11.5, color=ink)
    right.set_yticks(right_positions, BLOCK_LABELS)
    right.tick_params(axis="y", length=0, pad=10)
    right.set_xlim(-0.655, 0.04)
    right.set_ylim(-0.55, 3.55)
    right.xaxis.set_major_locator(FixedLocator((-0.6, -0.4, -0.2, 0)))
    right.xaxis.set_major_formatter(FuncFormatter(lambda x, _: "0" if x == 0 else f"{x:.1f}".replace("-", "−")))
    right.set_xlabel("Contribution to observed − null T\n(mean log units)", labelpad=13)

    for ax in (left, right):
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="x", length=4, pad=6)

    fig.text(0.045, 0.205,
             "Descriptive ratios; no effect confidence intervals.\n"
             "Positive pairs: primary 66; source ≥3 62; later ≥2 35.",
             fontsize=10.3, color=muted, va="top", linespacing=1.6)
    fig.text(0.565, 0.205,
             f"{western_share:.1%} of the primary difference comes from −90° to 0°.\n"
             f"The four contributions sum to {total_difference:.3f}.".replace("-", "−"),
             fontsize=10.3, color=navy, va="top", linespacing=1.6)
    fig.text(0.045, 0.106,
             "T = mean log(1 + nearest recorded-source coast-gap upper bound in km).",
             fontsize=10.2, color=ink)
    fig.text(0.045, 0.071,
             "Fixed recorded-source conditional comparison; no causal or first-arrival claim.",
             fontsize=10.2, color=ink)
    fig.text(0.045, 0.032, "Source: exact_null_diagnostic.json · fixed analysis dated 2026-10-30",
             fontsize=9, color=muted)

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    png, svg = output_stem.with_suffix(".png"), output_stem.with_suffix(".svg")
    fig.savefig(png, dpi=240, metadata={"Software": "chocho exact conditional result renderer"})
    fig.savefig(svg, metadata={"Date": None, "Creator": "chocho exact conditional result renderer"})
    plt.close(fig)
    return {
        "source": str(source), "source_sha256": DIAGNOSTIC_SHA256,
        "outputs": [str(png), str(svg)], "ratios": ratios,
        "primary_block_contributions": contributions,
        "primary_total_difference": total_difference,
        "longitude_1_share": western_share,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_DIRECTORY / "exact_null_diagnostic.json")
    parser.add_argument("--output-stem", type=Path,
                        default=DEFAULT_DIRECTORY / "recorded_source_exact_null_summary_20261030")
    args = parser.parse_args()
    print(json.dumps(render(args.source, args.output_stem), indent=2))


if __name__ == "__main__":
    main()
