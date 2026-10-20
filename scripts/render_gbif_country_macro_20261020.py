#!/usr/bin/env python3
"""Render the frozen country-grain result and its prespecified sensitivities."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = json.loads(args.result.read_text())
    grid, diag = result["all_72_prespecified_settings"], result["diagnostics"]
    rows = [
        ("Primary country comparison", result["primary"]),
        ("Human observations + specimens", grid["core|pooled|E3|q0.05|records1"]),
        ("Narrower climate analogue (10–90%)", grid["all_accepted_bases|pooled|E3|q0.10|records1"]),
        ("At least 5 focal records per country", grid["all_accepted_bases|pooled|E3|q0.05|records5"]),
        ("Additional country-area strata", diag["level1_area"]),
        ("Additional reporting-intensity strata", diag["level1_record_intensity"]),
        ("Added resources cover at least 90%", diag["added90"]),
        ("One botanical region per country", diag["single_region"]),
    ]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(12, 7.5))
    fig.subplots_adjust(left=.35, right=.83, top=.75, bottom=.23)
    fig.text(.05, .94, "Country records leave the host-opportunity association unresolved",
             fontsize=17, weight="bold", color="#16324f")
    fig.text(.05, .885, "229 taxa  |  140 eligible countries  |  33.82 million indexed country records",
             fontsize=12, color="#46576b")
    fig.text(.05, .835, "Primary: +3.91 pp (95% genus-bootstrap interval −6.42 to +14.94)",
             fontsize=13, weight="bold", color="#16324f")
    ax.axvline(0, color="#738294", linestyle="--", linewidth=1.2, zorder=1)
    for i, (label, row) in enumerate(rows):
        y = len(rows) - 1 - i
        point = row["observed_minus_expected_pp"]
        low, high = row["genus_bootstrap_excess_pp_95ci"]
        color = "#176b80" if i == 0 else "#7b96ac"
        ax.errorbar(point, y, xerr=[[point-low], [high-point]], fmt="o", color=color,
                    elinewidth=2, capsize=4, markersize=8 if i == 0 else 6, zorder=2)
        n = row["positive_country_species_n"]
        info = row["information"]["nonzero_conditional_variance_strata"]
        ax.text(1.03, y, f"{n:2d}   /   {info:2d}", transform=ax.get_yaxis_transform(),
                va="center", fontsize=10, color="#344558")
    ax.text(1.03, 1.08, "Positive cells /\ninformative strata", transform=ax.transAxes,
            fontsize=10, color="#344558", linespacing=1.4)
    ax.set_yticks(range(len(rows)), [label for label, _ in reversed(rows)])
    ax.set_ylim(-.6, len(rows)-.4)
    ax.set_xlim(-15, 30)
    ax.set_xticks([-10, 0, 10, 20, 30])
    ax.grid(axis="x", color="#e6ebf0", linewidth=.8)
    ax.tick_params(axis="y", length=0, pad=12)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Observed minus conditional expected overlap (percentage points)", labelpad=13)
    fig.text(.05, .13, "Information gate: 21 of 50 required strata in the primary comparison.",
             fontsize=12, weight="bold", color="#8a4d1c")
    fig.text(.05, .087, "All 72 prespecified settings and all 8 diagnostics have intervals spanning zero.\n"
             "Country co-coverage does not establish local host use, establishment, or a causal effect.",
             fontsize=10.5, color="#46576b", linespacing=1.5)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        fig.savefig(args.output_dir / f"gbif_country_macro_summary.{extension}",
                    dpi=180, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
