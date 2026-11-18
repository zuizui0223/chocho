#!/usr/bin/env python3
"""Render all-island sign maps and additive gains from the frozen atlas tables."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SPEC = Path(__file__).resolve().parent
BACKGROUND = ROOT / "data/exploratory/island51_observed_ocean_20261106/coast_0p0625/island51_coastal_ocean_graph.npz"
TARGETS = ("new_butterfly_jaccard", "earlier_butterfly_jaccard")
PRIMARY = ("C_to_HC", "H_to_HC")
CATEGORIES = {
    "both_positive": ("Both gains positive", "#0072b2"),
    "host_positive_circuit_negative": ("Hosts positive; circuits negative", "#009e73"),
    "host_negative_circuit_positive": ("Circuits positive; hosts negative", "#e69f00"),
    "both_negative": ("Both gains negative", "#777777"),
    "at_least_one_tie": ("At least one exact tie", "#b392c4"),
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def numeric_id(row):
    return int(row["island_id"])


def background():
    lock = json.loads((SPEC / "input_lock.json").read_text())
    pin = next(x for x in lock["files"] if x["path"] == str(BACKGROUND.relative_to(ROOT)))
    if digest(BACKGROUND) != pin["sha256"]:
        raise ValueError("The frozen source-grid background changed")
    with np.load(BACKGROUND, allow_pickle=False) as data:
        mask = data["fine_cell_to_node"] < 0
        lat = data["fine_latitude"]
        lon = data["fine_longitude"]
    if mask.shape != (2720, 5760):
        raise ValueError("Unexpected frozen coastline-grid geometry")
    fraction = mask.reshape(340, 8, 720, 8).mean(axis=(1, 3))
    extent = [float(lon.min()) - .03125, float(lon.max()) + .03125,
              float(lat.min()) - .03125, float(lat.max()) + .03125]
    return fraction, extent, "upper" if lat[0] > lat[-1] else "lower"


def place_labels(fig, ax, anchors):
    """Move text only; keep source-coordinate markers and leader endpoints fixed."""
    offsets = [(6, 8), (6, -16), (-30, 8), (-30, -16), (12, 27),
               (-36, 27), (14, -30), (-40, -30), (30, 5), (-55, 5),
               (12, 40), (-50, 40), (12, -44), (-50, -44)]
    occupied = []
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for island, x, y in anchors:
        annotation = ax.annotate(island, (x, y), xytext=offsets[0],
                                 textcoords="offset points", fontsize=8,
                                 color="#243b43", weight="bold",
                                 arrowprops={"arrowstyle": "-", "lw": .5, "color": "#87999f"},
                                 bbox={"facecolor": "white", "edgecolor": "none", "alpha": .88, "pad": .8},
                                 zorder=6)
        for offset in offsets:
            annotation.set_position(offset)
            annotation.update_positions(renderer)
            text_box = annotation.get_bbox_patch()
            annotation.update_bbox_position_size(renderer)
            box = text_box.get_window_extent(renderer).expanded(1.12, 1.18)
            inside = ax.bbox.contains(box.x0, box.y0) and ax.bbox.contains(box.x1, box.y1)
            if inside and not any(box.overlaps(other) for other in occupied):
                break
        occupied.append(box)


def save(fig, stem, out):
    png, svg = out / f"{stem}.png", out / f"{stem}.svg"
    fig.savefig(png, dpi=220, facecolor="white", metadata={"Software": "matplotlib; render_atlas.py"})
    fig.savefig(svg, facecolor="white", metadata={"Date": None})
    plt.close(fig)
    return [png, svg]


def render(out):
    geography_path = out / "island_geography.csv"
    status_path = out / "island_conditional_status.csv"
    gains_path = out / "island_gain_decomposition.csv"
    geography = rows(geography_path)
    statuses = rows(status_path)
    gains = rows(gains_path)
    geo = {r["island_id"]: r for r in geography}
    if len(geography) != 51 or len(geo) != 51 or len(statuses) != 102 or len(gains) != 510:
        raise ValueError("Atlas figure inputs have incomplete or duplicate support")
    for target in TARGETS:
        selected = [r for r in statuses if r["target"] == target]
        if len(selected) != 51 or {r["island_id"] for r in selected} != set(geo):
            raise ValueError("A map lacks one of the original island identities")
        if any(r["mse_category"] not in CATEGORIES for r in selected):
            raise ValueError("Unknown joint-sign category")
    later = {comparison: [r for r in gains if r["target"] == TARGETS[0] and r["comparison"] == comparison]
             for comparison in PRIMARY}
    for selected in later.values():
        if len(selected) != 51 or {r["island_id"] for r in selected} != set(geo):
            raise ValueError("A contribution panel lacks complete island support")
    label_ids = set()
    for selected in later.values():
        ranked = sorted(selected, key=lambda r: (-abs(float(r["mse_global_gain_contribution_pp"])), numeric_id(r)))
        label_ids.update(r["island_id"] for r in ranked[:5])

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#9dabb0", "text.color": "#233a43",
                         "axes.labelcolor": "#344b54", "xtick.color": "#53656d",
                         "ytick.color": "#344b54", "svg.fonttype": "none",
                         "svg.hashsalt": "chocho-island51-fixed-gain-atlas-20261106"})
    land, extent, origin = background()
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))
    fig.subplots_adjust(left=.045, right=.985, top=.76, bottom=.28, wspace=.09)
    fig.suptitle("Where the joint host–circuit model improves prediction", x=.045, ha="left", y=.965,
                 fontsize=19, weight="bold")
    fig.text(.045, .885, "All 51 original port groups · each point summarizes 50 held-island pair predictions", fontsize=11)
    counts = {}
    for ax, target, heading in zip(axes, TARGETS, ("A  Later records: 2011–2022 (primary)", "B  Earlier records: 2001–2010 (context)")):
        selected = [r for r in statuses if r["target"] == target]
        counts[target] = {category: sum(r["mse_category"] == category for r in selected) for category in CATEGORIES}
        ax.imshow(land, extent=extent, origin=origin,
                  cmap=LinearSegmentedColormap.from_list("source_land", ["#f5f9fb", "#d9e1e4"]),
                  vmin=0, vmax=1, interpolation="nearest", zorder=0)
        for category, (_, color) in CATEGORIES.items():
            group = [geo[r["island_id"]] for r in selected if r["mse_category"] == category]
            ax.scatter([float(r["representative_longitude"]) for r in group],
                       [float(r["representative_latitude"]) for r in group],
                       c=color, s=42, edgecolor="white", linewidth=.65, zorder=4)
        ax.set_xlim(-180, 180)
        ax.set_ylim(-85, 85)
        ax.set_aspect("auto")
        ax.set_xticks([-120, -60, 0, 60, 120])
        ax.set_yticks([-60, 0, 60])
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude" if target == TARGETS[0] else "")
        ax.grid(color="white", lw=.6, zorder=1)
        ax.set_title(heading + f"\nBoth gains positive: {counts[target]['both_positive']}/51", loc="left", fontsize=12, pad=10)
    anchors = [(island, float(geo[island]["representative_longitude"]), float(geo[island]["representative_latitude"]))
               for island in sorted(label_ids, key=int)]
    place_labels(fig, axes[0], anchors)
    handles = [Line2D([], [], marker="o", ls="", color=color, markersize=7, label=label)
               for label, color in CATEGORIES.values()]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.51, .115), ncol=3,
               frameon=False, fontsize=10)
    fig.text(.045, .075, "Hosts: C → HC. Circuits: H → HC. Both contrasts share HC predictions; their signs are not independent evidence.", fontsize=9, color="#53656d")
    fig.text(.045, .033, "Locations are spherical means of original ports. Numbered labels mark the fixed top-five absolute contributors per later contrast.", fontsize=9, color="#53656d")
    outputs = save(fig, "island_conditional_gain_maps", out)

    lookup = {comparison: {r["island_id"]: r for r in selected} for comparison, selected in later.items()}
    ordered = sorted(geo, key=lambda island: (-float(lookup["H_to_HC"][island]["mse_global_gain_contribution_pp"]), int(island)))
    labels = []
    for island in ordered:
        name = geo[island]["example_port_name"]
        if len(name) > 30:
            name = name[:29] + "…"
        codes = "/".join(json.loads(geo[island]["sampled_country_codes"]))
        labels.append(f"{island} · {name} [{codes}]")
    fig, axes = plt.subplots(1, 2, figsize=(14, 16), sharey=True)
    fig.subplots_adjust(left=.335, right=.96, top=.875, bottom=.08, wspace=.26)
    fig.suptitle("All 51 islands: contributions to the later-record prediction gain", x=.04, ha="left", y=.98,
                 fontsize=17, weight="bold")
    fig.text(.04, .95, "Each bar contributes to the overall MSE reduction. Positive bars improve the score; negative bars worsen it.", fontsize=10.5)
    fig.text(.04, .925, "Labels: landmass ID · smallest-ID example port [sampled-port country codes]. Port names are not island names.", fontsize=10)
    maxima = max(abs(float(r["mse_global_gain_contribution_pp"])) for selected in later.values() for r in selected)
    span = max(maxima, 1e-8)
    totals = {}
    for ax, contrast, title in zip(axes, PRIMARY, ("Hosts beyond circuits", "AIS + wind beyond hosts")):
        values = np.array([float(lookup[contrast][island]["mse_global_gain_contribution_pp"]) for island in ordered])
        totals[contrast] = float(values.sum())
        colors = ["#176b74" if value > 0 else "#b8663d" if value < 0 else "#777777" for value in values]
        y = np.arange(51)
        ax.barh(y, values, height=.68, color=colors, zorder=3)
        ax.axvline(0, color="#7f939b", lw=.8, zorder=4)
        ax.set_xlim(-span * 1.35, span * 1.35)
        ax.grid(axis="x", color="#e9eff1", zorder=0)
        ax.set_title(f"{title}\nSum = {totals[contrast]:+.4f} percentage points", fontsize=11, pad=12)
        ax.set_xlabel("Contribution (percentage points)", labelpad=9)
        ax.set_yticks(y)
        ax.tick_params(axis="y", length=0, pad=8)
        for yy, value in zip(y, values):
            ax.annotate(f"{value:+.3g}", (value, yy), xytext=(4 if value >= 0 else -4, 0),
                        textcoords="offset points", ha="left" if value >= 0 else "right", va="center", fontsize=7)
    axes[0].set_yticklabels(labels, fontsize=8)
    axes[0].set_ylim(50.8, -.8)
    axes[1].tick_params(axis="y", labelleft=False)
    fig.text(.04, .035, "Rows are ordered by the saved circuit increment, with numeric-ID tie breaking. All adverse contributions remain visible.", fontsize=10)
    fig.text(.04, .015, "These are descriptive contributions from pair predictions, not port-specific effects, local abundance or independently validated geographic drivers.", fontsize=9, color="#53656d")
    outputs += save(fig, "all_island_later_gain_contributions", out)
    receipt = {"schema": "chocho_island51_host_circuit_atlas_figures_v1",
               "input_sha256": {p.name: digest(p) for p in (geography_path, status_path, gains_path)},
               "background_path": str(BACKGROUND.relative_to(ROOT)), "background_sha256": digest(BACKGROUND),
               "renderer_sha256": digest(Path(__file__)),
               "output_sha256": {p.name: digest(p) for p in outputs},
               "mapped_island_ids": sorted(geo, key=int), "mapped_points_per_target": 51,
               "selected_later_label_ids": sorted(label_ids, key=int), "category_counts": counts,
               "all_island_contribution_order": ordered, "contribution_sums_pp": totals,
               "coordinate_meaning": "Original194-port group spherical means, not island polygon centroids or marine-terminal means.",
               "background_display": "0.5-degree averaged land fraction from the frozen0.0625-degree source-grid mask.",
               "no_independent_fold_intervals": True}
    (out / "figure_receipt.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"figures": [p.name for p in outputs], "all_islands": len(geo)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=SPEC)
    args = parser.parse_args()
    render(args.results_dir)
