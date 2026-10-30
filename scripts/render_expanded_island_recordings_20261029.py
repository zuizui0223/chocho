#!/usr/bin/env python3
"""Render the frozen panel's recording support; no ungated fitted trend lines."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/chocho-expanded-matplotlib-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/"data/exploratory/island_expansion_20261029"


def run(folder=DATA):
    with (folder/"expanded_island_recording_metrics.csv").open() as f:
        rows = list(csv.DictReader(f))
    late = [r for r in rows if r["window"] == "2018_2025"]
    eligible = [r for r in late if r["primary_eligible"] == "1"]
    result = json.loads((folder/"expanded_island_recording_analysis.json").read_text())
    groups = [(g,p) for g in ("under100", "atleast500") for p in ("within50", "beyond100")]
    all_n = [sum((r["gap_class"],r["port_class"]) == gp for r in late) for gp in groups]
    eligible_n = [sum((r["gap_class"],r["port_class"]) == gp for r in eligible) for gp in groups]
    colors = {"within50": "#156b77", "beyond100": "#c86b36"}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelcolor": "#263442", "text.color": "#263442"})
    fig, axes = plt.subplots(1,2,figsize=(12,5.4), gridspec_kw={"width_ratios": [1,1.25]})
    fig.subplots_adjust(left=.065, right=.98, bottom=.29, top=.79, wspace=.28)
    ax = axes[0]
    x = np.arange(4)
    ax.bar(x,all_n,width=.63,color="#e4e9ed",label="Fixed selected roster")
    ax.bar(x,eligible_n,width=.63,color=[colors[p] for g,p in groups],label="Meets recording screen")
    for j,(n,k) in enumerate(zip(all_n,eligible_n)):
        ax.text(j,n+.65,str(n),ha="center",fontsize=10,color="#657580")
        ax.text(j,max(k,.5)+.45,str(k),ha="center",fontweight="bold",fontsize=12)
    ax.set_xticks(x,["Near continent\nNear WPI", "Near continent\nFar WPI", "Remote\nNear WPI", "Remote\nFar WPI"],fontsize=9)
    ax.set_ylabel("Physical islands")
    ax.set_ylim(0,max(all_n,default=1)+5)
    ax.set_title("A  Selected islands vs usable recording pools",loc="left",fontsize=11,pad=13)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax = axes[1]
    for port,label,marker in (("within50","WPI within 50 km","o"),("beyond100","WPI beyond 100 km","^")):
        subset = [r for r in eligible if r["port_class"] == port]
        ax.scatter([float(r["mainland_coast_gap_upper_km"]) for r in subset],
                   [float(r["rarefied_species_20"]) for r in subset],
                   s=45,c=colors[port],marker=marker,alpha=.82,linewidth=.6,edgecolors="white",label=label)
    ax.set_xscale("symlog",linthresh=1)
    ax.set_xlabel("Continental coast gap (km; log scale above 1 km)")
    ax.set_ylabel("Expected species bins in 20 records")
    ax.set_ylim(.5,20.5)
    ax.set_title("B  Record-count-standardized recorded richness",loc="left",fontsize=11,pad=13)
    ax.legend(frameon=False,fontsize=9,loc="best")
    ax.grid(axis="y",alpha=.15)
    status = result["models"]["primary_late_m20"]["status"]
    qualifier = "Prespecified model support gates not met" if status != "FIT_EXPLORATORY_ASSOCIATION" else "Exploratory panel association; residual source and geographic confounding"
    fig.suptitle("Butterfly records across the fixed 96-island panel",x=.065,y=.97,ha="left",fontsize=18,fontweight="bold")
    fig.text(.065,.88,f"2018–2025  •  {len(eligible)} / {len(late)} islands meet the m=20 recording screen  •  {qualifier}",fontsize=10)
    fig.text(.065,.15,"Screen: ≥20 species-assigned records, ≥90% species assignment, ≥3 recorded years. Gray bars retain every selected island.",fontsize=9)
    fig.text(.065,.105,"Near / remote continent: <100 / ≥500 km. WPI classes refer to registered location proximity, including other landmasses.",fontsize=9)
    fig.text(.065,.06,"Finite-record richness is not true community richness. Empty records are not absences; no causal port or arrival effect is shown.",fontsize=9)
    for suffix in ("png","svg"):
        fig.savefig(folder/f"expanded_island_recording_support.{suffix}",dpi=180,facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    run()
