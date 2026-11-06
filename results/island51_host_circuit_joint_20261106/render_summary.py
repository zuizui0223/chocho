#!/usr/bin/env python3
"""Render the fixed model scores and conditional increments from saved CSVs."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


def render(results_dir: Path) -> None:
    scores_path = results_dir / "model_summary.csv"
    comparisons_path = results_dir / "comparison_summary.csv"
    scores = pd.read_csv(scores_path)
    comparisons = pd.read_csv(comparisons_path)
    targets = ("new_butterfly_jaccard", "earlier_butterfly_jaccard")
    models = ("G", "H", "C", "HC", "I")
    labels = (
        "G   Geography + uniform ocean",
        "H   G + candidate-host similarity",
        "C   G + native AIS + wind",
        "HC   Hosts + circuits together",
        "I   HC + two interaction terms",
    )
    selected = scores.loc[scores.target.eq(targets[0])].set_index("model").loc[list(models)]
    if len(selected) != 5 or selected.MSE.isna().any():
        raise ValueError("Expected one finite later-response MSE for each fixed model")
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#aab4b8", "axes.labelcolor": "#263a42",
        "text.color": "#1b3039", "xtick.color": "#42565f",
        "ytick.color": "#263a42", "svg.fonttype": "none",
        "svg.hashsalt": "chocho-island51-host-circuit-joint-20261106",
    })
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.3), gridspec_kw={"width_ratios": [1.1, 1]})
    fig.subplots_adjust(left=.255, right=.965, top=.76, bottom=.22, wspace=.75)
    fig.suptitle("Candidate hosts and marine circuits: conditional prediction", x=.04, ha="left", y=.97, fontsize=18, weight="bold")
    fig.text(.04, .9, "51 fixed islands · 1,275 pairs · one island withheld at a time · retrospective records", fontsize=11, color="#53666d")

    ax = axes[0]
    values = selected.MSE.to_numpy() * 1000
    y = np.arange(len(models))
    colors = ["#7c8990", "#598d7f", "#476d96", "#075e68", "#8b7297"]
    ax.scatter(values, y, s=90, c=colors, zorder=4, edgecolors="white", linewidths=.8)
    ax.axvline(values[0], color="#b4bec2", ls="--", lw=1, zorder=1)
    span = max(np.ptp(values), float(values.max()) * .1)
    for v, yy in zip(values, y):
        ax.text(v + span * .045, yy, f"{v:.4f}", va="center", fontsize=10)
    ax.set_xlim(max(0, float(values.min()) - span * .2), float(values.max()) + span * .6)
    ax.set_ylim(4.55, -.55)
    ax.set_yticks(y, labels)
    ax.tick_params(axis="y", length=0, pad=10)
    ax.grid(axis="x", color="#edf1f2", zorder=0)
    ax.set_xlabel("Held-island MSE (×10⁻³); lower is better", labelpad=12)
    ax.set_title("A  Later butterfly records (2011–2022)", loc="left", fontsize=12, weight="bold", pad=20)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))

    ax = axes[1]
    contrasts = ("C_to_HC", "H_to_HC")
    contrast_labels = ("Add hosts\nbeyond circuits", "Add AIS + wind\nbeyond hosts")
    target_labels = ("Later (primary)", "Earlier (temporal context)")
    target_colors = ("#075e68", "#ad6d3c")
    all_values = []
    points = []
    for ti, (target, color, label) in enumerate(zip(targets, target_colors, target_labels)):
        rows = comparisons.loc[comparisons.target.eq(target)].set_index("comparison")
        for ci, contrast in enumerate(contrasts):
            row = rows.loc[contrast]
            value = float(row.relative_MSE_reduction_percent)
            yy = ci + (-.14 if ti == 0 else .14)
            ax.scatter(value, yy, s=85, color=color, marker="o" if ti == 0 else "D", label=label if ci == 0 else None, zorder=4)
            all_values.append(value)
            points.append((value, yy, int(row.islands_improved), color))
    span = max(10., max(all_values + [0]) - min(all_values + [0]))
    for value, yy, improved, color in points:
        ax.text(value + .045 * span, yy, f"{value:+.2f}%  ({improved}/51)", va="center", fontsize=10, color=color)
    ax.axvline(0, color="#9caeb4", lw=1)
    ax.set_xlim(min(all_values + [0]) - .14 * span, max(all_values + [0]) + .85 * span)
    ax.set_ylim(1.5, -.6)
    ax.set_yticks([0, 1], contrast_labels)
    ax.tick_params(axis="y", length=0, pad=10)
    ax.grid(axis="x", color="#edf1f2", zorder=0)
    ax.set_xlabel("MSE reduction from that comparison's reference (%)", labelpad=12)
    ax.set_title("B  The two fixed conditional comparisons", loc="left", fontsize=12, weight="bold", pad=20)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(-.03, -.35), fontsize=10)

    fig.text(.04, .11, "Labels in B: relative MSE reduction; parentheses: islands with lower MSE. Positive values indicate improvement.", fontsize=10, color="#53666d")
    fig.text(.04, .065, "Folds share training data and pairs share islands. No independent-fold uncertainty intervals. Host records span 2001–2024.", fontsize=10, color="#53666d")
    fig.text(.04, .025, "All models retain distance and original port-count controls; later-response models also include earlier butterfly similarity.", fontsize=10, color="#53666d")
    png = results_dir / "joint_host_circuit_summary.png"
    svg = results_dir / "joint_host_circuit_summary.svg"
    fig.savefig(png, dpi=220, facecolor="white", metadata={"Software": "matplotlib; render_summary.py"})
    fig.savefig(svg, facecolor="white", metadata={"Date": None})
    plt.close(fig)
    receipt = {
        "description": "Exact saved model MSEs and the two conditional comparisons; no independent-fold intervals",
        "inputs_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (scores_path, comparisons_path)},
        "outputs": [png.name, svg.name],
        "outputs_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (png, svg)},
        "renderer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "values_in_left_panel": selected.reset_index()[["model", "MSE", "MAE"]].to_dict(orient="records"),
        "values_in_right_panel": comparisons.loc[comparisons.comparison.isin(contrasts)].to_dict(orient="records"),
    }
    (results_dir / "figure_receipt.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(f"Saved {png} and {svg}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    render(args.results_dir)
