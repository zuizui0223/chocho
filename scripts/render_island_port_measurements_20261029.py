#!/usr/bin/env python3
"""Render the measured geographic support, without a butterfly-effect claim."""
from pathlib import Path
import os
import json

import numpy as np
import pandas as pd
os.environ.setdefault("MPLCONFIGDIR", "/tmp/chocho-matplotlib-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/exploratory/island_port_20261029"


def run():
    manifest = json.loads((OUT / "island_geography_manifest.json").read_text())
    frame = pd.read_csv(OUT / "island_geography_exposures.csv")
    frame = frame[frame.antarctica_could_be_nearer == 0]
    remote = frame[frame.mainland_coast_gap_lower_km >= 500]
    assert len(frame) == manifest["summary"]["mainland_comparison_unambiguous_islands"]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "bold", "svg.fonttype": "none"})
    fig, (ax, bar) = plt.subplots(1, 2, figsize=(12.8, 6), gridspec_kw={"width_ratios": [1.3, 1]})
    color_near, color_far = "#286b9e", "#d6dce1"
    near = frame.wpi_entries_within_50km_certain > 0
    for mask, color, name in [(~near, "#929ca5", "No WPI entry certainly within 50 km"),
                               (near, color_near, "WPI entry certainly within 50 km")]:
        subset = frame[mask]
        ax.scatter(subset.mainland_coast_gap_upper_km,
                   subset.nearest_wpi_coast_distance_upper_km,
                   s=7, alpha=.35, c=color, linewidths=0, rasterized=True, label=name)
    ax.set_xscale("symlog", linthresh=1)
    ax.set_yscale("symlog", linthresh=1)
    ticks = [0, 1, 10, 100, 1000, 10000]
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(FixedLocator(ticks))
        axis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.axvline(500, ls="--", lw=1, c="#333333")
    ax.axhline(50, ls="--", lw=1, c="#333333")
    ax.set_xlim(0, 8000)
    ax.set_ylim(0, 8000)
    ax.set_xlabel("Distance to a mapped continental coast (km)")
    ax.set_ylabel("Distance from island coast to nearest WPI entry (km)")
    ax.set_title(f"A  Geographic exposures across {len(frame):,} islands", loc="left", pad=15)
    ax.text(530, 1.35, f"{manifest['summary']['remote500_WPI_certainly_within50km']:,} islands\nremote + near a WPI entry",
            fontsize=9, color="#164665", ha="left", va="bottom",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": .9, "pad": 4})
    ax.legend(loc="upper left", fontsize=8, frameon=False, markerscale=2)

    bounds = [(5, 20), (20, 100), (100, 1000), (1000, np.inf)]
    labels = ["5–20", "20–100", "100–1,000", "≥1,000"]
    counts = []
    for lo, hi in bounds:
        sub = remote[(remote.area_km2 >= lo) & (remote.area_km2 < hi)]
        n = int((sub.wpi_entries_within_50km_certain > 0).sum())
        f = int((sub.wpi_entries_within_50km_possible == 0).sum())
        counts.append((n, f, len(sub)-n-f, len(sub)))
    counts = np.array(counts)
    y = np.arange(4)
    nearpct = counts[:, 0]/counts[:, 3]*100
    farpct = counts[:, 1]/counts[:, 3]*100
    ambpct = counts[:, 2]/counts[:, 3]*100
    bar.barh(y, nearpct, color=color_near, height=.6, label="Within 50 km")
    bar.barh(y, farpct, left=nearpct, color=color_far, height=.6, label="Beyond 50 km")
    bar.barh(y, ambpct, left=nearpct+farpct, color="#cc9b42", height=.6, label="Sampling ambiguity")
    for pos, (n, f, ambiguous, total) in enumerate(counts):
        bar.text(nearpct[pos]/2, pos, f"{n:,}", color="white", ha="center", va="center", weight="bold")
        bar.text(nearpct[pos]+farpct[pos]/2, pos, f"{f:,}", ha="center", va="center", color="#26323a")
        bar.text(103, pos, f"n={total:,}", va="center", fontsize=9)
    bar.set_yticks(y, labels)
    bar.invert_yaxis()
    bar.set_xlim(0, 122)
    bar.set_xticks([0, 25, 50, 75, 100], ["0%", "25%", "50%", "75%", "100%"])
    bar.set_ylabel("Island area (km²)")
    bar.set_xlabel("Share of islands; labels give island counts")
    bar.set_title(f"B  Remote islands (≥500 km; n={len(remote):,})", loc="left", pad=15)
    bar.legend(loc="lower left", bbox_to_anchor=(0, -.24), ncols=1, frameon=False, fontsize=8)
    fig.suptitle("Natural isolation and mapped port proximity provide different comparisons", x=.065,
                 y=.98, ha="left", fontsize=15, weight="bold")
    fig.subplots_adjust(left=.075, right=.97, top=.85, bottom=.25, wspace=.34)
    fig.text(.075, .085,
             "GSHHG 2.3.7 L1 islands ≥5 km²; WPI 2025: 3,688 usable entries. 120 islands with unresolved Antarctic comparison excluded.\n"
             "Spherical coast samples ≤1 km apart; axes linear below 1 km, logarithmic above. Three remote islands straddle the 50 km port threshold.\n"
             "WPI lists major ports, not all harbors; nearby entries may be on other landmasses. Geographic support only, no butterfly or shipping effect.",
             fontsize=8.5, ha="left", va="top", linespacing=1.5, color="#45515c")
    for extension in ("png", "svg"):
        fig.savefig(OUT / f"island_port_geographic_support.{extension}", dpi=180, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    run()
