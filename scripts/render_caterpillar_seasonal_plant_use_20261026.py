#!/usr/bin/env python3
"""Original Caterpillars Count seasonal field use OR with site-bootstrap intervals.

Do not conflate larval Lepidoptera with identified butterfly-specific population
changes, larval survival, competition or benefit of alien plants.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/"provenance/reviewer_defenses/results/caterpillars_count_native_alien_seasonal_larval_MH_20261026.json"


def render(svg,png):
    d=json.loads(RESULT.read_text(encoding="utf-8"))
    data=d["main_all_methods"]
    if data["status"]!="ORIGINAL_SOURCE_MATCHED_SEASONAL_LARVAL_OCCUPANCY_ESTIMATE_QUALIFIED":
        raise ValueError("3-season original-source matched site bootstrap not qualified")
    seasons=["early","peak","late"]
    title=["Early (DOY 91–151)","Summer (DOY 152–212)","Late (DOY 213–273)"]
    y=np.arange(3)
    val=np.array([data["observed"][s]["OR_alien_vs_native"] for s in seasons])
    ci=np.array([data["observed"][s]["95pct_site_cluster_bootstrap_OR"] for s in seasons])
    if not np.isfinite(val).all() or not np.isfinite(ci).all() or np.any(ci<=0):
        raise ValueError("Invalid original data OR/boot source receipt")
    fig,ax=plt.subplots(1,2,figsize=(10.7,4.4),constrained_layout=True)
    a,b=ax
    a.errorbar(val,y,xerr=[val-ci[:,0],ci[:,1]-val],fmt="o",capsize=4,lw=1.7,markersize=7)
    a.axvline(1,ls="--",lw=1.2)
    a.set_yticks(y,title)
    a.set_xscale("log")
    a.set_xlim(.09,4.2)
    a.set_xlabel("Alien / native caterpillar survey odds ratio (MH)")
    a.set_title("(A) Site–year–family–method matched",fontsize=10.8,fontweight="bold")
    a.invert_yaxis()
    a.grid(axis="x",alpha=.17)
    difference=data["late_minus_peak_logOR"]
    margin=data["late_minus_peak_site_cluster_95_logOR"]
    secondary=data["early_minus_peak_logOR"]
    margin_secondary=data["early_minus_peak_site_cluster_95_logOR"]
    ys=np.array([0,1])
    ds=np.array([difference,secondary])
    rang=np.array([[margin[0],margin[2]],[margin_secondary[0],margin_secondary[2]]])
    b.errorbar(ds,ys,xerr=[ds-rang[:,0],rang[:,1]-ds],fmt="o",
               lw=1.7,capsize=4,markersize=7)
    b.axvline(0,ls="--",lw=1.2)
    b.set_yticks(ys,["Late minus summer","Early minus summer"])
    b.set_xlim(-2.0,2.7)
    b.set_xlabel("Difference in log(OR alien/native)")
    b.set_title("(B) No clear seasonal improvement",fontsize=10.8,fontweight="bold")
    b.grid(axis="x",alpha=.17)
    b.invert_yaxis()
    fig.suptitle("Observed caterpillar larvae on woody plants: seasonal alien–native contrast\n"
                 "Original Caterpillars Count! standardized branch surveys, 2010–2024",
                 fontsize=12.5)
    fig.supxlabel("2,000 site-cluster resamples (1,984 valid) · Caterpillars include moths: not butterfly-specific breeding success",
                  fontsize=8.5)
    svg.parent.mkdir(parents=True,exist_ok=True)
    png.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(svg)
    fig.savefig(png,dpi=185)
    plt.close(fig)
    print("Source-verified original caterpillar survey bootstrap seasonal contrast plotted")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--svg",type=Path,required=True)
    p.add_argument("--png",type=Path,required=True)
    a=p.parse_args()
    render(a.svg,a.png)
