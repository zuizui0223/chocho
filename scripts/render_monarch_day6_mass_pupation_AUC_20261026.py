#!/usr/bin/env python3
"""Visualize exact monarch Day6 growth-to-pupation source-linked prognosis.

Figure is a within-original-experiment prognostic result, NOT a causal effect
of making a caterpillar heavier; Top/Bottom treatment 20241 is unresolved.
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
RECEIPT=ROOT/"provenance/reviewer_defenses/results/fisher_monarch_Day6_mass_pupation_ID_AUC_20261026.json"


def draw(svg,png):
    d=json.loads(RECEIPT.read_text(encoding="utf-8"))
    if (
        d["status"]!="ORIGINAL_INDIVIDUAL_STAGE_LINK_CONFIRMED_SOURCE_CONDITIONAL_PROGNOSIS_ONLY"
        or d["individual_stage_link_quality"]["Day6_individuals_measured"]!=270
        or d["individual_stage_link_quality"]["Day6_and_later_pupation_ID_matched"]!=190
        or not d["yearTrial20241_entire_recorded_TopBottom_group_rename_AUC_invariant"]
    ):
        raise ValueError("Original individual stage-prognostic source gate not qualified")
    rows=d["observed_six_original_strata"]
    st=d["fixed_prognostic_analysis"]
    if len(rows)!=6 or st["permutation_draws"]!=20000:
        raise ValueError("Fixed six source strata or 20k null history changed")
    vals=np.array([r["within_stratum_prognostic_mass_AUC"] for r in rows])
    overall=float(st["observed_equal_6_strata_Mann_Whitney_AUC"])
    null=np.asarray(st["null_q025_median_q975"],dtype=float)
    if not math.isclose(float(vals.mean()),overall,abs_tol=1e-10):
        raise ValueError("Actual source stratum ranks do not match archived primary AUC")
    labs=[
        "2023 T2: Top (19 / 34)",
        "2023 T2: Bottom (15 / 35)",
        "2024 T1*: Top (36 / 42)",
        "2024 T1*: Bottom (30 / 46)",
        "2024 T2: Top (50 / 59)",
        "2024 T2: Bottom (40 / 54)",
    ]
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9})
    fig,axs=plt.subplots(1,2,figsize=(11.8,5.2),constrained_layout=True)
    y=np.arange(6)
    a,b=axs
    a.scatter(vals,y,s=70,marker="o",label="Within-source recorded leaf group")
    a.axvline(.5,lw=1,ls="--",alpha=.65,label="No ranking ability")
    a.set_yticks(y,labs)
    a.invert_yaxis()
    a.set_xlim(.4,.9)
    a.set_xlabel("AUC: Day-6 mass ranking versus later pupation")
    a.set_title("(A) Three experimental trials, six fixed strata",fontweight="bold",fontsize=11)
    a.grid(axis="x",alpha=.16)
    a.text(.025,.04,
           "* Source YearTrial 20241 Top/Bottom identity unresolved;\n"
           "  within-group mass ranking is invariant to renaming labels",
           transform=a.transAxes,ha="left",va="bottom",fontsize=7.5,
           bbox={"facecolor":"white","edgecolor":"none","alpha":.9})
    b.axvspan(null[0],null[2],alpha=.17,label="20k fixed-margin null central 95%")
    b.axvline(.5,lw=1,ls="--",alpha=.75)
    xpoint=[st["trial_mean_AUC"][q] for q in ("20232","20241","20242")]
    b.scatter(xpoint,[0,1,2],s=63,label="Equal leaf groups within trial")
    b.scatter([overall],[3],s=88,marker="D",label="Equal six-strata AUC")
    b.set_yticks([0,1,2,3],["2023 T2","2024 T1*","2024 T2","All six strata"])
    b.invert_yaxis()
    b.set_xlim(.4,.9)
    b.set_xlabel("Concordance AUC")
    b.set_title("(B) Variability across only three Day-6 trials",
                fontweight="bold",fontsize=11)
    b.grid(axis="x",alpha=.16)
    b.text(.025,.91,f"Primary AUC {overall:.3f}\n"
           f"Permutation-null median {null[1]:.3f}\n"
           f"20,000 permutations; two-sided tail "
           f"{st['two_sided_deviation_from_half_empirical_fraction']:.5f}",
           transform=b.transAxes,ha="left",va="top",fontsize=8.3,
           bbox={"facecolor":"white","edgecolor":"none","alpha":.9})
    b.legend(loc="lower left",fontsize=7.4,frameon=False)
    fig.suptitle("Observed monarch day-six mass precedes later recorded pupation\n"
                 "Fisher et al. (2026) original individual source; n=270 measured, 190 pupated",
                 fontsize=13.0)
    fig.supxlabel("Original-study prognosis, not a causal mass intervention. 20241 treatment labels unresolved. Pupation is not adult emergence.",
                  fontsize=8.4)
    svg.parent.mkdir(parents=True,exist_ok=True)
    png.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(svg)
    fig.savefig(png,dpi=175)
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--svg",type=Path,required=True)
    p.add_argument("--png",type=Path,required=True)
    a=p.parse_args()
    draw(a.svg,a.png)
    print("PASS: source-checked monarch individual day-six→pupation AUC and null envelope visualized")
