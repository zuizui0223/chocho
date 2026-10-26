#!/usr/bin/env python3
"""Visual QA: monarch four-block pupation source-label contradiction.

Source-locked to archived output, not an effect image asserting biological
fitness or adult recruitment. The original treatment-label discrepancy is
front-and-centre and neither scenario is selected as 'truth'.
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
RECEIPT=ROOT/"provenance/reviewer_defenses/results/monarch_leaf_age_pupation_stage_conversion_20261026.json"


def draw(svg,png):
    d=json.loads(RECEIPT.read_text(encoding="utf-8"))
    if (d["status"]!="SOURCE_TREATMENT_LABEL_INDETERMINATE_TWO_SCENARIOS_RETAINED"
        or d["source_assigned_larvae"]!=480 or
        len(d["original_trial_label_disagreements"])!=2):
        raise ValueError("Pupation source label contradiction audit unqualified")
    cases=list(d["experimental_survival_scenarios"].values())
    if len(cases)!=2:raise ValueError("Both original source-trust scenarios are required")
    labels=["2023 trial 1","2023 trial 2","2024 trial 1 *","2024 trial 2"]
    keys=("20231","20232","20241","20242")
    y=np.arange(4)
    orig=np.array([cases[0]["within_4_trial_blocks"][k]["Top_minus_Bottom_probability"]*100
                   for k in keys])
    alt=np.array([cases[1]["within_4_trial_blocks"][k]["Top_minus_Bottom_probability"]*100
                  for k in keys])
    fig,axes=plt.subplots(1,2,figsize=(11.5,5.3),constrained_layout=True)
    ax=axes[0]
    ax.plot(orig,y,marker="o",ls="-",lw=1.35,label="S0 original survival summary")
    ax.plot(alt,y,marker="s",ls="--",lw=1.35,label="S1 hypothetical pupal record counts")
    ax.axvline(0,color="gray",ls=":",lw=1)
    ax.set_yticks(y,labels)
    ax.invert_yaxis()
    ax.set_xlim(-15,22)
    ax.set_xlabel("Top minus Bottom pupation (% of assigned larvae)")
    ax.set_title("(A) One disputed trial changes treatment direction",fontweight="bold",fontsize=10.5)
    ax.grid(axis="x",alpha=.15)
    ax.legend(loc="lower right",fontsize=7.7,frameon=False)
    ax.text(.04,.95,"* 2024 trial 1: summary 30/36, both pupal files 36/30",
            transform=ax.transAxes,ha="left",va="top",
            fontsize=8.2,bbox={"facecolor":"white","alpha":.87,"edgecolor":"none"})
    ax=axes[1]
    point=np.array([c["MH_common_OR_top_vs_bottom"] for c in cases])
    ci=np.array([c["trial_cluster_bootstrap"]["common_odds_ratio_bootstrap_95"]
                 for c in cases],dtype=float)
    ax.errorbar(point,[0,1],xerr=[point-ci[:,0],ci[:,1]-point],fmt="o",
                lw=1.8,markersize=7,capsize=5)
    ax.axvline(1,color="gray",ls="--",lw=1)
    ax.set_yticks([0,1],["S0: trust survival summary","S1: trust both pupal files"])
    ax.invert_yaxis()
    ax.set_xlim(.5,2.65)
    ax.set_xlabel("Trial-stratified pupation odds ratio, Top / Bottom")
    ax.set_title("(B) Source scenario determines uncertainty",fontweight="bold",fontsize=10.5)
    ax.grid(axis="x",alpha=.15)
    p0=cases[0]["conditional_fixed_margin_exact"]["two_sided_probability_order_p"]
    p1=cases[1]["conditional_fixed_margin_exact"]["two_sided_probability_order_p"]
    ax.text(.05,.53,
            f"Fixed-margin conditional tests\nS0 p={p0:.4f} · S1 p={p1:.4f} (hypothetical)",
            transform=ax.transAxes,ha="left",va="center",fontsize=8.0,
            bbox={"facecolor":"white","edgecolor":"none","alpha":.93})
    fig.suptitle("Observed monarch development on young versus old milkweed leaves\n"
                 "Source treatment-label inconsistency prevents a unique survival conclusion",
                 fontsize=13.0)
    fig.supxlabel("Original Fisher et al. (2026) experiment; 480 larvae in 4 YearTrial blocks. "
                  "Pupation ≠ adult emergence or butterfly population growth.",
                  fontsize=8.3)
    svg.parent.mkdir(parents=True,exist_ok=True)
    png.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(svg)
    fig.savefig(png,dpi=175)
    plt.close(fig)


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--svg",type=Path,required=True)
    ap.add_argument("--png",type=Path,required=True)
    a=ap.parse_args()
    draw(a.svg,a.png)
    print("PASS: 4 original monarch blocks and 2 contradictory pupation source scenarios verified")
