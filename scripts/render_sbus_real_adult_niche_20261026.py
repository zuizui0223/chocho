#!/usr/bin/env python3
"""Visual evidence ladder: monitored butterflies versus host network nulls.

All data are source-locked aggregate scalar null histories; no site/census rows.
Figure does NOT show or imply an estimated interspecific competition effect.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
R=ROOT/"provenance/reviewer_defenses/results"
A=R/"sbus_illinois_adult_codetection_hostnull_20261026.json"
B=R/"sbus_illinois_temporal_host_cofluctuation_fragility_20261026.json"
C=R/"sbus_illinois_host_bipartite_degree_null_20261025.json"


def load(p):
    return json.loads(p.read_text(encoding="utf-8"))


def histogram(ax,data,obs,title,desc,bins=23):
    x=np.asarray(data,dtype=float)
    if len(x)!=199 or not np.isfinite(x).all():
        raise ValueError("Only exact 199 original null histories are reportable")
    lim=np.quantile(x,[.025,.5,.975])
    ax.hist(x,bins=bins,alpha=.6,density=True,label="Fixed-degree host-network null")
    ax.axvline(obs,ls="--",color="black",lw=1.8,label="Observed adult field association")
    ax.axvline(lim[1],ls=":",color="black",lw=1.2,label="Null median")
    ax.set_title(title,fontweight="bold",fontsize=10.3)
    ax.set_xlabel(desc,fontsize=9)
    ax.set_ylabel("Source-null density")
    ax.text(.03,.95,
            f"Observed {obs:+.4f}\nNull median {lim[1]:+.4f}\n95% null [{lim[0]:+.4f}, {lim[2]:+.4f}]",
            ha="left",va="top",transform=ax.transAxes,fontsize=8.2,
            bbox={"boxstyle":"round","facecolor":"white","edgecolor":"none","alpha":.88})
    ax.grid(alpha=.12)
    return lim


def draw(svg,png):
    a,b,c=load(A),load(B),load(C)
    if(a.get("null_phi",{}).get("status")!="EXPLORATORY_SOURCE_CONDITIONED_BIPARTITE_NULL_QUALIFIED"
       or b.get("independent_seed_null_equal_pairs",{}).get("status")!="SOURCE_CONDITIONED_GRAPH_NULL_QUALIFIED"
       or c.get("result",{}).get("status")!="POSTHOC_SOURCE_CONDITIONED_NETWORK_NULL_QUALIFIED"):
        raise ValueError("A prerequisite original exact-source null did not qualify")
    co=a["null_phi"]
    syn0=a["within_site_temporal"]
    syn1=b["independent_seed_null_equal_pairs"]
    competitive=c["result"]
    if not (co["null_draws_valid"]==199 and syn0["null_draws_valid"]==199 and
            syn1["n_graphs"]==199 and competitive["complete_numerically_valid"]==199):
        raise ValueError("Source annual/bipartite randomization history missing")
    fig,axes=plt.subplots(2,2,figsize=(12,7.8),constrained_layout=True)
    hist_co=histogram(
        axes[0,0],co["null_difference_history"],co["observed"]["difference"],
        "(A) Geographic co-detection",r"Shared minus nonshared adult species-pair φ"
    )
    ax=axes[0,1]
    h0=np.asarray(syn0["null_difference_history"])
    h1=np.asarray(syn1["history"])
    obs=float(syn1["observed"])
    ax.hist(h0,bins=21,density=True,alpha=.4,label="Source null 1")
    ax.hist(h1,bins=21,density=True,alpha=.4,label="Independent source null 2")
    ax.axvline(obs,color="black",ls="--",lw=2,label="Observed annual synchrony")
    ax.set_title("(B) Within-site annual co-fluctuation",fontweight="bold",fontsize=10.3)
    ax.set_xlabel("Shared minus nonshared correlation of adult annual rate anomalies")
    ax.set_ylabel("Source-null density")
    ax.text(.03,.95,
            f"Observed {obs:+.4f}\nNull medians: "
            f"{syn0['null_difference_q025_median_q975'][1]:+.4f} / "
            f"{syn1['null_q025_median_q975'][1]:+.4f}\n"
            f"Both right-tail fractions 0.030\n"
            f"Secondary: 3 endpoints evaluated",
            ha="left",va="top",transform=ax.transAxes,fontsize=8.2,
            bbox={"boxstyle":"round","facecolor":"white","edgecolor":"none","alpha":.88})
    ax.grid(alpha=.12)
    histogram(
        axes[1,0],competitive["null_beta_history"],
        competitive["observed_beta"],
        "(C) Lagged peer abundance association",
        "Shared minus nonshared lag coefficient"
    )
    ax=axes[1,1]
    ax.axis("off")
    ax.set_title("(D) Biological interpretation",fontweight="bold",fontsize=10.3)
    notes=(
        "Potential host geography: expanded (reconstruction)\n\n"
        "Field co-detection: descriptive overlap, within host null\n\n"
        "Adult annual synchrony: positive, borderline posthoc\n"
        "                              and not multiplicity-robust\n\n"
        "Lagged competition proxy: within host null; not established\n\n"
        "Not measured: local host biomass / larval feeding / survival,\n"
        "                  introduced host establishment dates"
    )
    ax.text(.025,.94,notes,transform=ax.transAxes,va="top",ha="left",
            fontsize=10,linespacing=1.35)
    axes[0,0].legend(loc="lower right",frameon=False,fontsize=7.2)
    axes[0,1].legend(loc="lower right",frameon=False,fontsize=7.2)
    axes[1,0].legend(loc="upper right",frameon=False,fontsize=7.2)
    fig.suptitle("Documented host sharing and adult butterfly community data:\n"
                 "distinct spatial, temporal and lagged tests",
                 fontsize=13.5)
    fig.supxlabel("Exploratory source-conditioned randomizations, not biological competition experiments",
                  fontsize=9.1)
    svg.parent.mkdir(parents=True,exist_ok=True)
    png.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(svg)
    fig.savefig(png,dpi=180)
    plt.close(fig)
    print("Source-locked null histograms verified, co-detection and temporal interpretation bounded")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--svg",type=Path,required=True)
    p.add_argument("--png",type=Path,required=True)
    args=p.parse_args()
    draw(args.svg,args.png)
