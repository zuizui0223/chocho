#!/usr/bin/env python3
"""Render the preserved exploratory chains without refitting or resampling."""
from pathlib import Path
import csv
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/exploratory/island_recorded_sources_20261029"


def main():
    report = json.loads((DATA / "recorded_source_diagnostic.json").read_text())
    with (DATA / "recorded_source_chain_samples.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    scenarios = [("primary", "Primary: later records ≥1"),
                 ("source_records_ge3", "Earlier source records ≥3"),
                 ("late_records_ge2", "Later records ≥2")]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "bold", "svg.hashsalt": "island-recorded-sources-20261029"})
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 5.9), sharex=True, sharey=True)
    fig.patch.set_facecolor("#f8fafc")
    for axis, (key, title) in zip(axes, scenarios):
        result = report["analyses"][key]
        diagnostic = result["diagnostic"]
        values = np.array([float(r["mean_log1p_source_gap_km"]) for r in rows if r["analysis"] == key])
        if len(values) != diagnostic["stored_samples_total"] or not np.isclose(values.mean(), diagnostic["null_mean_T"]):
            raise ValueError("Published sample count or mean differs from diagnostic")
        if diagnostic["qualified_tail"]:
            raise ValueError("This figure's explicit unqualified-tail caption needs review")
        axis.hist(values, bins=np.linspace(7.2, 8.85, 56), color="#6798af", alpha=.9)
        axis.axvline(diagnostic["observed_T"], color="#bb4b28", lw=2.5, label="Observed")
        axis.axvline(diagnostic["null_mean_T"], color="#24465a", lw=1.7,
                     linestyle="--", label="Sampled reference mean")
        axis.set_title(title, loc="left", fontsize=11.5, pad=13)
        axis.set_xlim(7.2, 8.85)
        axis.set_xticks([7.3, 7.6, 7.9, 8.2, 8.5, 8.8])
        axis.set_xlabel("Mean log(1 + source-coast gap in km)", labelpad=10)
        axis.grid(axis="y", alpha=.17)
        axis.set_axisbelow(True)
        support = result["support"]
        ratio = diagnostic["observed_to_null_geometric_mean_one_plus_gap_ratio"]
        annotation = (f"{support['positive_cells']} positive species–island pairs\n"
                      f"Descriptive geometric ratio: {ratio:.3f}\n"
                      f"R̂ {diagnostic['split_rhat']:.4f}  |  ESS {diagnostic['pooled_ess']:,.0f}\n"
                      f"Lower-tail draws: 0 / {len(values):,}")
        axis.text(.025, .97, annotation, transform=axis.transAxes, va="top", fontsize=9.5,
                  bbox={"boxstyle": "round,pad=.5", "fc": "white", "ec": "#dce3e8", "alpha": .96})
    axes[0].set_ylabel("Stored MCMC draws (dependent)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.5, .19),
               ncol=2, frameon=False, fontsize=10)
    fig.suptitle("Earlier recorded sources are closer in the observed configurations", x=.055,
                 y=.98, ha="left", fontsize=17, fontweight="bold", color="#17374a")
    fig.text(.055, .915, "Exploratory fixed-96 panel  •  Species × longitude-block and island totals preserved",
             fontsize=11, color="#445d6b")
    fig.text(.055, .145, "No qualified tail probability in any scenario.", fontsize=12,
             fontweight="bold", color="#a34225")
    fig.text(.055, .115, "Primary and source ≥3 also fail ESS ≥1,000. Zero tail draws leave tail precision unresolved.\n"
             "Ratios describe sampled configurations; they are not dispersal, colonization or causal port effects.\n"
             "Recorded sources include selected islands only; mainland and unsampled islands are omitted.",
             fontsize=9.5, linespacing=1.5, color="#445d6b", va="top")
    fig.subplots_adjust(left=.06, right=.99, top=.83, bottom=.32, wspace=.12)
    fig.savefig(DATA / "recorded_source_diagnostic.png", dpi=180, facecolor=fig.get_facecolor())
    fig.savefig(DATA / "recorded_source_diagnostic.svg", facecolor=fig.get_facecolor(),
                metadata={"Date": None})
    plt.close(fig)


if __name__ == "__main__":
    main()
