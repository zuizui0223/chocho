#!/usr/bin/env python3
"""Source-audited alternative Figure 3: taxonomic breadth and butterfly occurrence.

Only genuine frozen 239-species source CSV, the independent 32-species
occurrence CSV, and their original null receipts are plotted. Scatter jitter
is deterministic display-only and does not enter Spearman calculations.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

SOURCE_239_SHA256 = "7b2a20387d656dbfcd7b3c38ec6a7fe2505474ad51f783dd1f201eb6b7eabd30"


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8",newline="") as h:
        return list(csv.DictReader(h))


def midranks(values: list[float]) -> np.ndarray:
    x = np.asarray(values,dtype=float)
    ind = np.argsort(x,kind="stable")
    out = np.zeros(len(x),dtype=float)
    start = 0
    while start < len(ind):
        stop = start+1
        while stop<len(ind) and x[ind[stop]]==x[ind[start]]:
            stop += 1
        out[ind[start:stop]] = (start+stop+1)/2
        start = stop
    return out


def spearman(x: list[float], y: list[float]) -> float:
    if len(x)<3 or len(y)!=len(x):
        raise ValueError("malformed paired metrics")
    a, b = midranks(x), midranks(y)
    if np.std(a)==0 or np.std(b)==0:
        raise ValueError("noninformative ranks")
    return float(np.corrcoef(a,b)[0,1])


def infer_source_values(
    species_csv: Path,
    occurrence_csv: Path,
    climate_null_json: Path,
    occurrence_null_json: Path,
) -> dict:
    if hashlib.sha256(species_csv.read_bytes()).hexdigest()!=SOURCE_239_SHA256:
        raise ValueError("frozen species figure-source bytes have changed")
    sp = read_csv(species_csv)
    occ = read_csv(occurrence_csv)
    null = json.loads(occurrence_null_json.read_text(encoding="utf-8"))
    precision = json.loads(climate_null_json.read_text(encoding="utf-8"))
    if len(sp)!=239 or len(occ)!=32:
        raise ValueError("panel size drift")
    if len({r["species"] for r in sp})!=239 or len({r["species"] for r in occ})!=32:
        raise ValueError("species duplicates or missing identities")
    family=[float(r["host_family_count"]) for r in sp]
    native=[int(r["native_resource_units"]) for r in sp]
    log_gain=[float(r["log_resource_expansion"]) for r in sp]
    ratio=spearman(family,log_gain)
    rho_native=spearman(family,native)
    if not math.isclose(ratio,0.00798541214928876,abs_tol=1e-9):
        raise ValueError("original specialization gain correlation changed")
    if not math.isclose(rho_native,0.27625647485507027,abs_tol=1e-9):
        raise ValueError("native breadth specialization source correlation changed")
    if not math.isclose(ratio,float(precision["rho"]),abs_tol=1e-9):
        raise ValueError("original expansion effect receipt changed")
    adequate=[
        r for r in sp if r["host_taxonomy_lower_bound_adequate"].strip()=="True"
    ]
    if len(adequate)!=215:
        raise ValueError("adequately resolved host inventory species changed")
    q75=float(np.percentile([int(r["native_resource_units"]) for r in adequate],75))
    upper_one=[
        r for r in adequate
        if float(r["host_family_count"])==1
        and int(r["native_resource_units"])>=q75
    ]
    if len(upper_one)!=18:
        raise ValueError("widely distributed one-family butterfly count changed")
    informative=[r for r in occ if int(r["observed_outside_native"])>0]
    total=sum(int(r["observed_outside_native"]) for r in informative)
    found=sum(int(r["outside_native_units_explained_by_introduced_host_ranges"]) for r in informative)
    if len(informative)!=23 or total!=115 or found!=66:
        raise ValueError("occurrence panel overlap or effort support changed")
    if (null["species_with_outside_native_occurrences"]!=23 or
        null["observed_recovered_species_x_units"]!=66 or
        null["level1_composition_matched_null"]["median"]!=49 or
        null["level1_composition_matched_null"]["q025"]!=43 or
        null["level1_composition_matched_null"]["q975"]!=56):
        raise ValueError("region-matched observed-vs-null evidence drift")
    return {
        "schema":"chocho_butterfly_resource_novelty_breadth_occurrence_figure_v0.1",
        "239_species_source_SHA256":SOURCE_239_SHA256,
        "32_occurrence_source_SHA256":hashlib.sha256(occurrence_csv.read_bytes()).hexdigest(),
        "focal_species":239,"adequate_host_species":215,
        "source":sp,
        "host_family_vs_native_resource_spearman":rho_native,
        "host_family_vs_log_gain_spearman":ratio,
        "host_family_vs_log_gain_bootstrap_95ci":precision["bootstrap"]["ci95"],
        "adequate_native_q75":q75,
        "wide_one_family_species_n":len(upper_one),
        "wide_one_family_names":sorted(r["species"] for r in upper_one),
        "occurrence_panel_species":32,
        "informative_species":len(informative),
        "outside_native_species_region_units":total,
        "observed_recovered":found,
        "region_matched_null": {
            "median":49,
            "q025":43,
            "q975":56,
        },
        "scientific_boundary":"Host range and butterfly occurrence are distinct evidence types; the latter does not establish local larval use or dispersal causality.",
    }


def xjitter(species: str) -> float:
    h=int(hashlib.sha256(("butterfly-figure3-jitter|"+species).encode()).hexdigest()[:8],16)
    return ((h/0xffffffff)*2-1)*0.15


def render(data:dict,png:Path,svg:Path)->None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(15.2,6.4))
    gs=fig.add_gridspec(1,3,left=.073,right=.965,top=.78,bottom=.22,
                       width_ratios=[1.05,1.05,.80],wspace=.31)
    axs=[fig.add_subplot(gs[0,i]) for i in range(3)]
    rows=data["source"]
    fam=np.asarray([float(r["host_family_count"]) for r in rows])
    xr=np.asarray([float(r["host_family_count"])+xjitter(r["species"]) for r in rows])
    native=np.asarray([int(r["native_resource_units"]) for r in rows])
    gain=np.asarray([float(r["log_resource_expansion"]) for r in rows])
    eligible=np.asarray([r["host_taxonomy_lower_bound_adequate"]=="True" for r in rows])
    wide=set(data["wide_one_family_names"])
    highlighted=np.asarray([r["species"] in wide for r in rows])
    fig.suptitle("Resource geography is not interchangeable with diet breadth or butterfly occurrence",
                 fontsize=15.2,x=.073,y=.95,ha="left",weight="semibold")
    fig.text(.073,.887,
        "239-species potential host-resource panel; independent 32-species butterfly occurrence panel",
        fontsize=10.2)
    a,b,c=axs
    a.scatter(xr[~eligible],native[~eligible],s=20,alpha=.25,marker="x")
    a.scatter(xr[eligible],native[eligible],s=18,alpha=.46)
    a.scatter(xr[highlighted],native[highlighted],s=45,alpha=.95,facecolors="none",linewidths=1.2)
    a.set_yscale("log")
    a.set_xlim(0.2,min(18,max(fam)+1))
    a.set_title("A  |  Native resource breadth",loc="left",fontsize=11,pad=12)
    a.set_xlabel("Number of recorded host-plant families")
    a.set_ylabel("Native host-resource WGSRPD3 units (log scale)")
    a.text(.98,.98,
        f"Spearman ρ = {data['host_family_vs_native_resource_spearman']:.3f}\n"
        "18 one-family high-resource species\nin the adequately resolved subset (n=215)",
        transform=a.transAxes,ha="right",va="top",fontsize=9.1,
        bbox={"facecolor":"white","alpha":.8,"edgecolor":"none"})
    b.scatter(xr[~eligible],gain[~eligible],s=20,alpha=.25,marker="x")
    b.scatter(xr[eligible],gain[eligible],s=20,alpha=.49)
    b.scatter(xr[highlighted],gain[highlighted],s=45,facecolors="none",linewidths=1.2)
    b.set_title("B  |  Proportional resource increase",loc="left",fontsize=11,pad=12)
    b.set_xlabel("Number of recorded host-plant families")
    b.set_ylabel("log(contemporary / native resource units)")
    b.set_xlim(a.get_xlim())
    b.text(.035,.98,
        f"Spearman ρ = {data['host_family_vs_log_gain_spearman']:.3f}\n"
        "95% bootstrap CI: [−0.111, 0.128]",
        transform=b.transAxes,ha="left",va="top",fontsize=9.2,
        bbox={"facecolor":"white","alpha":.8,"edgecolor":"none"})
    n=data["outside_native_species_region_units"]
    null=data["region_matched_null"]
    values=[data["observed_recovered"],null["median"]]
    c.bar([0,1],values,width=.56,alpha=.72)
    c.errorbar([1],[null["median"]],
               yerr=[[null["median"]-null["q025"]],[null["q975"]-null["median"]]],
               fmt="none",capsize=6,ecolor="black",linewidth=1.2)
    for i,v in enumerate(values):
        c.text(i,v+4,f"{v}/{n}\n({100*v/n:.1f}%)",ha="center",fontsize=10.3)
    c.set_xticks([0,1],["Observed\nresource overlap","Region-matched\nnull median"])
    c.set_ylim(0,85)
    c.set_ylabel("Outside-native butterfly × region occurrence units")
    c.set_title("C  |  Presence is not larval host use",loc="left",fontsize=11,pad=12)
    for ax in axs:
        ax.spines[["right","top"]].set_visible(False)
        ax.grid(axis="y",alpha=.12)
        ax.set_axisbelow(True)
    fig.text(.073,.119,
        "A–B: visual jitter only; correlations use actual numeric counts. Hollow points mark geographically broad one-family species.",
        fontsize=9.2)
    fig.text(.073,.075,
        "C: 66 of 115 outside-native occurrence units (23 informative of 32 butterflies); matched null 49 (95% interval 43–56).",
        fontsize=9.2)
    fig.text(.073,.031,
        "Potential resource geography ≠ observed butterfly density, local larval feeding, successful recruitment or causal range expansion.",
        fontsize=8.8)
    png.parent.mkdir(parents=True,exist_ok=True)
    svg.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(png,dpi=180,bbox_inches="tight")
    fig.savefig(svg,bbox_inches="tight")
    plt.close(fig)


def main()->None:
    p=argparse.ArgumentParser()
    for key in ("species_csv","occurrence_csv","precision_json","overlap_null_json","output_png","output_svg","output_data_json"):
        p.add_argument("--"+key.replace("_","-"),required=True,type=Path)
    args=p.parse_args()
    data=infer_source_values(args.species_csv,args.occurrence_csv,
                             args.precision_json,args.overlap_null_json)
    render(data,args.output_png,args.output_svg)
    result={k:v for k,v in data.items() if k!="source"}
    args.output_data_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_data_json.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":"RENDERED_BUTTERFLY_RESOURCE_ECOLOGICAL_CONTEXT","sources":result,
        "png":str(args.output_png),"svg":str(args.output_svg)},indent=2))


if __name__=="__main__":
    main()
