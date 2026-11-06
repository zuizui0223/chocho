#!/usr/bin/env python3
"""Render the observed marine source domain and fixed island-heldout results."""
from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def render(input_root, result_dir):
    root, out = Path(input_root), Path(result_dir)
    summaries = pd.read_csv(out / "model_comparison_summary.csv")
    ports = pd.read_csv(out / "original194_port_attachment_audit.csv")
    included = ports["included"].astype(str).str.lower().eq("true")
    with np.load(root / "ais/commercial_AIS_mean_0p5deg_2015_2021.npz") as source:
        ais = source["AIS_mean"].astype(float)
        valid = source["source_valid_mask"]
    with np.load(root / "coast_0p0625/island51_coastal_ocean_graph.npz") as coast:
        components = coast["node_component"]
        sizes = np.bincount(components, weights=coast["node_fine_water_cells"])
        marine = components == np.argmax(sizes)
        domain = np.zeros_like(valid)
        domain[coast["node_block_row"][marine], coast["node_block_col"][marine]] = True
    field = np.ma.array(np.log10(1+ais), mask=~(domain & valid))
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none"})
    fig = plt.figure(figsize=(12, 10.2), layout="constrained", facecolor="white")
    grid = fig.add_gridspec(2, 1, height_ratios=[1.15, 1.])
    ax = fig.add_subplot(grid[0])
    ax.set_facecolor("#e5e7eb")
    image = ax.imshow(field, extent=(-180, 180, -85, 85), origin="upper",
                      cmap="viridis", vmin=0,
                      vmax=float(np.quantile(field.compressed(), .995)),
                      interpolation="nearest", aspect="auto", rasterized=True)
    ax.scatter(ports.loc[included, "longitude"], ports.loc[included, "latitude"],
               s=20, c="#fb923c", edgecolor="white", linewidth=.4,
               label=f"Marine-access ports ({int(included.sum())})", zorder=3)
    ax.scatter(ports.loc[~included, "longitude"], ports.loc[~included, "latitude"],
               s=26, marker="x", c="#ef4444", linewidth=1,
               label=f"Outside marine-access criterion ({int((~included).sum())})", zorder=4)
    ax.set(xlim=(-180, 180), ylim=(-70, 80), xticks=np.arange(-180, 181, 60),
           yticks=[-60, -30, 0, 30, 60], xlabel="Longitude", ylabel="Latitude")
    ax.set_title("A   Observed shipping field and the original island ports", loc="left", weight="bold", pad=12)
    ax.legend(loc="lower left", frameon=True, facecolor="white", framealpha=.94, fontsize=9)
    bar = fig.colorbar(image, ax=ax, shrink=.92, pad=.015, fraction=.028)
    bar.set_label("log10(1 + published AIS density index)", fontsize=10)
    ax.text(.99, .025, "World Bank / IMF, Jan 2015–Feb 2021\n0.5° source overview averages; color clipped at the 99.5th percentile",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            bbox={"facecolor": "white", "alpha": .85, "edgecolor": "none", "pad": 4})

    bx = fig.add_subplot(grid[1])
    comparisons = ["uniform_vs_ais", "uniform_vs_wind", "uniform_vs_ais_and_wind"]
    names = ["AIS shipping density", "Mean wind speed", "AIS + mean wind speed"]
    targets = [("host_jaccard", "Candidate introduced-host similarity", "#0f766e"),
               ("new_butterfly_jaccard", "Later butterfly similarity", "#7c3aed")]
    yy = np.arange(3)
    values = []
    for k, (target, label, color) in enumerate(targets):
        rows = summaries.set_index(["comparison", "target"])
        gain = np.array([float(rows.loc[(c, target), "relative_MSE_reduction_percent"]) for c in comparisons])
        improved = [int(rows.loc[(c, target), "islands_improved"]) for c in comparisons]
        positions = yy + (k-.5)*.30
        bx.barh(positions, gain, height=.26, color=color, label=label, zorder=2)
        values.extend(gain.tolist())
        for y, x, count in zip(positions, gain, improved):
            bx.annotate(f"{x:+.2f}%   ({count}/51 islands)", xy=(x, y),
                        xytext=(5 if x >= 0 else -5, 0), textcoords="offset points",
                        ha="left" if x >= 0 else "right", va="center", fontsize=10, color=color)
    span = max(max(values)-min(values), 2.)
    bx.set_xlim(min(0, min(values))-span*.65, max(0, max(values))+span*.65)
    bx.axvline(0, color="#334155", linewidth=1)
    bx.set(yticks=yy, yticklabels=names, xlabel="Reduction in held-out prediction MSE (%)")
    bx.invert_yaxis()
    bx.grid(axis="x", color="#e5e7eb", linewidth=.7, zorder=0)
    bx.set_title("B   Information added beyond geography and the uniform-ocean circuit", loc="left", weight="bold", pad=14)
    bx.legend(loc="upper center", bbox_to_anchor=(.5, -.17), ncol=2, frameon=False)
    fig.suptitle("51 islands · 1,275 island pairs · 182 marine-access ports", fontsize=18, weight="bold")
    fig.supxlabel("Each fold holds out all 50 pairs involving one island. Controls: distance, original port counts, and earlier butterfly similarity for the later target.\nPositive values indicate lower MSE. The 51 folds share data; these are descriptive prediction results, without independent-fold confidence intervals.",
                  fontsize=9, y=-.025)
    for extension in ("png", "svg"):
        fig.savefig(out / f"observed_ocean_circuit_summary.{extension}", dpi=180,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path, required=True)
    args = parser.parse_args()
    render(args.input_root, args.result_dir)
