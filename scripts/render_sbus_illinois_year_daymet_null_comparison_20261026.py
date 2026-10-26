#!/usr/bin/env python3
"""Source-backed three-panel uncertainty comparison: host-sharing adult synchrony.

Not a 1-km site weather estimate, host-plant phenology, competition experiment,
or new GEB submission figure. Original field data and fixed 640 pair cohorts
remain unchanged in the three posthoc source diagnostics.
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
R=ROOT/"provenance/reviewer_defenses/results"


def read(file):
    return json.loads((R/file).read_text(encoding="utf-8"))


def panel(ax, name, obs, source, normal_note):
    if source["status"]!="EXPLORATORY_SOURCE_CONDITIONED_BIPARTITE_NULL_QUALIFIED":
        raise ValueError("Incomplete 199-graph source null cannot be shown as qualified")
    hist=np.asarray(source["null_difference_history"],dtype=float)
    if hist.shape!=(199,) or not np.isfinite(hist).all():
        raise ValueError("Missing or nonfinite fixed-source 199 graph null values")
    predicted=np.quantile(hist,[.025,.5,.975],method="linear")
    recorded=source["null_difference_q025_median_q975"]
    if not np.allclose(predicted,recorded,atol=1e-10,rtol=0):
        raise ValueError("Source null interval does not equal archived primary history")
    if not math.isclose(obs,source["observed"]["difference"],abs_tol=1e-10):
        raise ValueError("Observed field values differ from archived condition")
    ax.hist(hist,bins=20,density=True,alpha=.7,label="Family-conditioned host-graph null")
    ax.axvline(obs,ls="--",lw=2,color="black",label="Observed adult synchrony")
    ax.axvline(predicted[1],ls=":",lw=1.1,color="black",label="Null median")
    ax.grid(alpha=.12)
    ax.set_ylabel("Source network-null density",fontsize=8.7)
    ax.set_xlabel("Mean annual rate correlation: shared minus nonshared hosts",fontsize=8.7)
    ax.set_title(name,fontsize=10.5,fontweight="bold")
    q=source.get("right_tail_fraction_exploratory")
    if q is None:
        q=source.get("right_tail_fraction_exploratory_only")
    text=(f"Observed {obs:+.4f}\n"
          f"Null median {predicted[1]:+.4f}\n"
          f"95% [{predicted[0]:+.4f},{predicted[2]:+.4f}]\n"
          f"Exploratory upper fraction {q:.3f}")
    ax.text(.025,.965,text,transform=ax.transAxes,ha="left",va="top",
            fontsize=8.2,bbox={"facecolor":"white","edgecolor":"none","alpha":.9,"boxstyle":"round"})
    ax.text(.98,.035,normal_note,transform=ax.transAxes,ha="right",va="bottom",fontsize=7.5)


def draw(svg,png):
    a=read("sbus_illinois_10h_source_cofluctuation_20261026.json")
    b=read("sbus_illinois_year_shock_niche_synchrony_20261026.json")
    c=read("sbus_illinois_coarse_daymet_weather_host_niche_20261026.json")
    if (a["status"]!="SOURCE_10H_CORRECTION_FULL_REPLAY_QUALIFIED"
        or b["status"]!="EXPLORATORY_YEAR_SHOCK_SENSITIVITY_SOURCE_QUALIFIED"
        or c["status"]!="COARSE_14REGION_WEATHER_PROXY_SENSITIVITY_QUALIFIED_EXPLORATORY_ONLY"
        or b["calendar_year_correction"]["original_fixed_temporal_pairs"]!=640
        or c["source_weather_quality"]["fixed_temporal_species_pairs"]!=640):
        raise ValueError("At least one original or environmental sensitivity is source-unqualified")
    fig,ax=plt.subplots(2,2,figsize=(11.8,8.0),constrained_layout=True)
    panel(ax[0,0],"(A) Original source: five >10h surveys removed",
          a["corrected_observed"]["difference"],
          a["corrected_source_network_null"],"Adult source monitoring, 640 fixed pairs")
    panel(ax[0,1],"(B) Remove species-specific common-year shocks",
          b["corrected_observed"]["difference"],
          b["source_degree_graph_null"],"Other sites used to estimate year effect")
    panel(ax[1,0],"(C) Remove regional modeled weather proxies",
          c["coarse_seasonally_weather_adjusted_observed_pair_contrast"]["difference"],
          c["weather_adjusted_bipartite_host_null"],
          "Spring GDD5 + summer rain; 14 rounded bins")
    ax[1,1].axis("off")
    ax[1,1].set_title("(D) Inference boundary",fontsize=11.0,fontweight="bold")
    info=(
        "All tests: exact same 640 adult butterfly pairs\n"
        "62 shared-host, 578 nonsharing, 199 null graphs\n\n"
        "Original: marginal source-network association\n"
        "Calendar-year adjustment: inside null envelope\n"
        "Coarse independent Daymet: inside null envelope\n\n"
        "121 monitoring sites but only 14 rounded\n"
        "public climate-location classes\n\n"
        "No localized host-plant phenology, leaf biomass,\n"
        "larval survival or plant-arrival dates measured.\n"
        "Source nulls are NOT competition experiments.\n\n"
        "Three related discovery endpoints were examined;\n"
        "these exploratory tails are not confirmatory p-values."
    )
    ax[1,1].text(.035,.96,info,transform=ax[1,1].transAxes,
                 ha="left",va="top",fontsize=10.1,linespacing=1.25)
    ax[0,0].legend(frameon=False,fontsize=7.4,loc="lower right")
    fig.suptitle("Butterfly host-sharing annual synchrony: environment-falsification checks\n"
                 "Illinois 2000–2020 source-monitored adult encounter rates",
                 fontsize=13.3)
    fig.supxlabel("Values are observational fixed-source correlations; all weather-site positions are coarse")
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
    print("PASS: independently archived 199 graph nulls match observed adult temporal environment sensitivities")
