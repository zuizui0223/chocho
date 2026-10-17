#!/usr/bin/env python3
"""Supplementary figure: butterfly-specific new-resource region identity turnover.

The Jaccard measures agreement BETWEEN TWO HOST-INVENTORY CLASSIFICATIONS,
not observed temporal redistribution, diet evolution, or larval establishment.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path


def render(data:dict,png:Path,svg:Path)->None:
    if data.get("schema")!="chocho_butterfly_bce_added_region_identity_turnover_v0.1":
        raise ValueError("unsupported source Jaccard summary")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sp=data["species"]
    agg=data["affected_81"]
    if len(sp)!=81 or not all(key in agg for key in (
        "intersection_over_union_jaccard","retained_in_both",
        "reclassified_as_native","newly_added_by_augmentation")):
        raise ValueError("incomplete paired source data")
    fig=plt.figure(figsize=(13.0,6.8))
    fig.subplots_adjust(left=.075,right=.96,bottom=.19,top=.80,wspace=.29)
    axs=fig.subplots(1,2,gridspec_kw={"width_ratios":[1.1,.9]})
    ax,ab=axs

    x=[v["base_added"] for v in sp]
    y=[v["augmented_added"] for v in sp]
    m=max(x+y)
    ax.scatter(x,y,s=21,alpha=.63)
    ax.plot([0,m+8],[0,m+8],linestyle="--",linewidth=1.1)
    ax.set_xlim(-6,m+12)
    ax.set_ylim(-6,m+12)
    ax.set_aspect("equal",adjustable="box")
    ax.set_title("A  |  New regional opportunity changes by species",loc="left",pad=12)
    ax.set_xlabel("Initially counted new resource regions (HOSTS)")
    ax.set_ylabel("Reclassified new regions (HOSTS + BCE)")
    ax.grid(alpha=.12)
    ax.set_axisbelow(True)
    for name in ("Lycaena dispar","Minois dryas","Thymelicus lineola"):
        rows=[v for v in sp if v["species"]==name]
        if len(rows)!=1:
            raise ValueError(f"example host-inventory case absent: {name}")
        v=rows[0]
        shift={
            "Lycaena dispar":(10,5),
            "Minois dryas":(-18,8),
            "Thymelicus lineola":(7,-12),
        }[name]
        ax.scatter([v["base_added"]],[v["augmented_added"]],s=64,
                   facecolors="none",edgecolors="black",linewidths=1.2,zorder=5)
        ax.annotate(name,(v["base_added"],v["augmented_added"]),
                    textcoords="offset points",xytext=shift,fontsize=9,
                    ha="left" if shift[0]>0 else "right")

    sizes=[
        agg["retained_in_both"],agg["reclassified_as_native"],
        agg["newly_added_by_augmentation"]
    ]
    text=[
        "New in both source inventories",
        "Original new, now native-supported",
        "Previously unknown new resources"
    ]
    cursor=0.
    colors=plt.rcParams["axes.prop_cycle"].by_key()["color"]
    legend=[]
    for i,n in enumerate(sizes):
        width=100*n/agg["union_of_either_source"]
        bar=ab.barh(0,width,left=cursor,height=.37,color=colors[i])
        ab.text(cursor+width/2,0,f"{n:,}\n{width:.1f}%",
                ha="center",va="center",weight="medium",fontsize=10)
        legend.append(bar[0])
        cursor+=width
    if abs(cursor-100)>1e-9:
        raise ValueError("source resource identities do not exhaust the union")
    ab.set_xlim(0,100)
    ab.set_ylim(-.57,.64)
    ab.set_yticks([])
    ab.set_xticks([0,20,40,60,80,100])
    ab.set_xlabel("Proportion of the 6,445 butterfly × region union")
    ab.set_title("B  |  50% disagreement in the location of new resources",
                 loc="left",pad=12)
    ab.grid(axis="x",alpha=.12)
    ab.set_axisbelow(True)
    ab.legend(legend,text,bbox_to_anchor=(0,-.20),loc="upper left",
              frameon=False,fontsize=9)

    fig.suptitle("The map of 'new' butterfly resources changes with host records",
        x=.075,ha="left",y=.96,fontsize=16,weight="semibold")
    fig.text(.075,.888,
        "81 butterflies gained evidence-filtered European candidate host links; "
        "counts reflect potential butterfly × WGSRPD3 resources, not colonization",
        fontsize=10.5)
    fig.text(.075,.07,
        f"Jaccard across {agg['union_of_either_source']:,} original-or-augmented new cells = "
        f"{agg['intersection_over_union_jaccard']:.3f}; "
        f"{agg['species_increase_net']} species gain net regions, "
        f"{agg['species_decrease_net']} lose, "
        f"{agg['species_equal_net']} retain the same count. "
        "Net change does not identify place-level stability.",fontsize=9.5)
    fig.text(.075,.035,
        "The BCE/Clarke host sample is Europe-biased; 'new' is a source-dependent "
        "map category, not confirmed larval use or chronological colonization.",fontsize=9)
    png.parent.mkdir(parents=True,exist_ok=True)
    svg.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(png,dpi=190,bbox_inches="tight")
    fig.savefig(svg,bbox_inches="tight")
    plt.close(fig)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--summary-json",type=Path,required=True)
    p.add_argument("--output-png",type=Path,required=True)
    p.add_argument("--output-svg",type=Path,required=True)
    a=p.parse_args()
    x=json.loads(a.summary_json.read_text(encoding="utf-8"))
    render(x,a.output_png,a.output_svg)
    print(json.dumps({"status":"RENDERED_SOURCE_CONDITIONAL_BUTTERFLY_REGION_IDENTITY",
       "png":str(a.output_png),"svg":str(a.output_svg),
       "affected_butterflies":x["affected_81"]["original_new_resource_cells"],
       "jaccard":x["affected_81"]["intersection_over_union_jaccard"]},indent=2))


if __name__=="__main__":
    main()
