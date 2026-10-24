#!/usr/bin/env python3
"""Two-panel source-locked illustration of potential resource gain and single-host support.

This is a visual description of reconstructed known-host geography,
NOT independent butterfly occupancy or resilience evidence.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "data/frozen/figure_sources/anthropogenic_species_metrics.csv"
H = ROOT / "data/frozen/figure_sources/host_contribution_metrics.csv"
RESULT = ROOT / "provenance/reviewer_defenses/results/butterfly_geographic_gain_fragility_20261024.json"
BLOBS = {
    A: "9f83d2d9ac8e46f1cfd4ded572c7a023273503dd",
    H: "4cf4851499ca2233ff2adadf1c4ee7f0e863f23d",
}


def source(path):
    data = path.read_bytes()
    sha = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\x00" + data).hexdigest()
    if sha != BLOBS[path]:
        raise ValueError("Source file changed: " + str(path))
    with path.open(newline="", encoding="utf-8") as fd:
        return {z["species"]: z for z in csv.DictReader(fd)}


def points():
    names, supporters = source(A), source(H)
    if len(names) != 239 or len(supporters) != 239 or set(names) != set(supporters):
        raise ValueError("Unexpected fixed species join")
    rows = []
    for name in sorted(names):
        a, h = names[name], supporters[name]
        old = float(a["native_resource_units"])
        added = float(a["introduced_added_units"])
        if old <= 0 or added < 0:
            raise ValueError("Nonpositive native or negative new resource unit")
        if added == 0:
            continue
        rows.append({
            "name": name,
            "gain": math.log(float(a["contemporary_resource_units"]) / old),
            "single": 1-float(h["fraction_added_units_with_multiple_contributing_hosts"]),
            "native": math.log(old),
            "hosts": math.log1p(float(a["resolved_host_species"])),
            "family": float(a["host_family_count"]),
            "adequate": a["host_taxonomy_lower_bound_adequate"] == "True",
        })
    if len(rows) != 206:
        raise ValueError("New resources cohort differs")
    return rows


def draw(out_svg: Path, out_png: Path):
    rows = points()
    r = json.loads(RESULT.read_text(encoding="utf-8"))
    target = r["analyses"]["primary_all"]
    if target["species"] != 206 or target["genera"] != 124:
        raise ValueError("Frozen numeric receipt not matching plotted sample")
    x = np.array([s["gain"] for s in rows])
    y = np.array([s["single"] for s in rows])
    w = np.array([s["adequate"] for s in rows], dtype=bool)
    Z = np.array([[1,s["native"],s["hosts"],s["family"]] for s in rows])
    xp = x-Z @ np.linalg.lstsq(Z,x,rcond=None)[0]
    yp = y-Z @ np.linalg.lstsq(Z,y,rcond=None)[0]
    coef = float(np.linalg.lstsq(
        np.column_stack([Z, x]), y, rcond=None)[0][-1])
    expected = target["coefficient_log_gain"]
    if not math.isclose(coef,expected,abs_tol=1e-9):
        raise ValueError(f"Plot adjusted slope differs from result: {coef} != {expected}")
    low, high = target["genus_bootstrap_95_percentile"]

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.85), constrained_layout=True)

    ax=axes[0]
    ax.scatter(x[w], 100*y[w], s=21, alpha=.50, label=f"Host information adequate (n={sum(w)})")
    ax.scatter(x[~w], 100*y[~w], marker="x", s=34, linewidths=1.1,
               label=f"Less complete catalogue (n={sum(~w)})")
    ax.set_xlabel("Log contemporary / native resource-region ratio")
    ax.set_ylabel("New regions with one recorded host (%)")
    ax.set_title("(A) Species-level mapped opportunity")
    ax.text(.04,.06,"Spearman ρ = −0.200\n206 expanded butterflies",
            transform=ax.transAxes, va="bottom", ha="left",
            bbox={"boxstyle":"round", "facecolor":"white","alpha":.85,"edgecolor":"none"})
    ax.set_ylim(-5,105)
    ax.grid(alpha=.18)
    ax.legend(frameon=False,loc="upper right",fontsize=7.5)

    ax=axes[1]
    ax.scatter(xp[w],yp[w]*100,s=21,alpha=.48)
    ax.scatter(xp[~w],yp[~w]*100,s=34,marker="x",linewidths=1.1)
    line = np.linspace(float(xp.min()),float(xp.max()),160)
    centre=100*coef*line
    bounds=np.array([100*low*line,100*high*line])
    ax.fill_between(line,bounds.min(axis=0),bounds.max(axis=0),alpha=.13,
                    label="Genus-bootstrap slope interval")
    ax.plot(line,centre,lw=2,label="Adjusted species-level association")
    ax.axhline(0,lw=.8,ls="--",alpha=.6)
    ax.set_xlabel("Resource gain residual after controls")
    ax.set_ylabel("Single-host share residual (percentage points)")
    ax.set_title("(B) After native area & host-breadth controls")
    ax.text(.04,.06,f"Adjusted β = {coef:+.3f}\n95% genus CI [{low:+.3f}, {high:+.3f}]",
            transform=ax.transAxes,va="bottom",ha="left",
            bbox={"boxstyle":"round","facecolor":"white","alpha":.85,"edgecolor":"none"})
    ax.grid(alpha=.18)
    ax.legend(frameon=False,loc="upper right",fontsize=7.5)
    fig.suptitle("Butterfly resource gain and catalogue-recorded host support",fontsize=14)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_svg,format="svg")
    fig.savefig(out_png,format="png",dpi=160)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--svg",type=Path,required=True)
    ap.add_argument("--png",type=Path,required=True)
    z=ap.parse_args()
    draw(z.svg,z.png)
    print("PASS: both panels source-verified against exact 239-species files and bootstrap receipt")


if __name__ == "__main__":
    main()
