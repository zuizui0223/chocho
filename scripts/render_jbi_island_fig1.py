#!/usr/bin/env python3
"""Render source-anchored Fig. 1, preserving three genuinely different analysis grains.

Illustrative dashed arrows represent corresponding statistical information,
NEVER observed shipping, propagule movement or colonization.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch

ROOT=Path(__file__).resolve().parents[1]
LEDGER=ROOT/"manuscript/island_resources_marine_access_20261108/EVIDENCE_LEDGER_V02.json"


def audit():
    l=json.loads(LEDGER.read_text(encoding="utf-8"))
    assert l["schema"]=="chocho_JBI_island_resource_marine_v02_evidence_ledger"
    ids={x["id"] for x in l["evidence"]}
    assert len(ids)==16
    for key in ["island_resource","matched_island_excess","port_flora",
                "port_hosts","pair_correspondence","marine_graph",
                "both_unseen","AIS_local_placebo","independent_gbif_record_overlap"]:
        assert key in ids
    get=lambda x:next(y for y in l["evidence"] if y["id"]==x)["values"]
    assert get("island_resource")["island_native"]==4829
    assert get("island_resource")["island_current"]==10570
    assert get("island_resource")["mainland_native"]==21215
    assert get("island_resource")["mainland_current"]==29901
    assert get("marine_graph")["original_ports"]==194
    assert get("marine_graph")["marine_ports"]==182
    assert get("pair_correspondence")["dependent_pairs"]==1275
    assert get("independent_gbif_record_overlap")["recorded_butterfly_in_new_opportunity_cells"]==13
    return l


def rectangle(ax,xy,w,h,face,title,detail,edge="#263a50",title_size=13,body_size=10.2):
    x,y=xy
    patch=FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.013,rounding_size=0.022",
                         facecolor=face,edgecolor=edge,linewidth=1.5)
    ax.add_patch(patch)
    ax.text(x+.021,y+h-.035,title,fontsize=title_size,weight="bold",va="top",ha="left")
    ax.text(x+.021,y+h-.102,detail,fontsize=body_size,va="top",ha="left",linespacing=1.42)
    return patch


def arrow(ax,from_xy,to_xy,dash=False,label=None,offset=0):
    v=FancyArrowPatch(from_xy,to_xy,arrowstyle="-|>",mutation_scale=13,
                      linewidth=1.55,color="#4d6173",linestyle="--" if dash else "-")
    ax.add_patch(v)
    if label:
        ax.text((from_xy[0]+to_xy[0])/2,(from_xy[1]+to_xy[1])/2+offset,label,
                fontsize=9,ha="center",va="center",color="#415569")
    return v


def render(ledger,out):
    out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":10.5,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig,ax=plt.subplots(figsize=(16.7,8.4))
    ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
    fig.suptitle("Two layers of island butterfly resource geography",
                 fontsize=18,weight="bold",y=.975)
    ax.text(.028,.93,"RESOURCES  ·  separate botanical-region population",fontsize=12,
            color="#325877",weight="bold")
    ax.text(.028,.565,"ACCESS  ·  separate physical island-group population",fontsize=12,
            color="#2d6550",weight="bold")

    # Plant redistribution changes candidate host availability, not the reverse.
    # The left-to-right conceptual arrow is NOT an observed causal transport
    # estimate: it only describes which constructed resource layer precedes
    # the separate recorded-butterfly correspondence.
    rectangle(ax,(.03,.69),.245,.20,"#edf3fa","Introduced plant assembly",
              "Introduced flora includes larval hosts\nMatched island excess: +10.7 pp\nPort models: 98 insular regions",
              title_size=12.6)
    rectangle(ax,(.37,.69),.245,.20,"#edf3fa","Regional host opportunity",
              "239 fixed butterflies\n126 insular / 237 mainland WGSRPD3\n+118.9% islands; +40.9% mainland",
              title_size=12.6)
    arrow(ax,(.29,.79),(.356,.79))
    rectangle(ax,(.03,.33),.245,.20,"#eef7f0","Over-water accessibility",
              "51 island/landmass groups\n194 indexed ports; 182 marine terminals\n189,550 nodes; 370,996 edges",
              title_size=12.6)
    rectangle(ax,(.37,.33),.245,.20,"#eef7f0","Marine circuit predictors",
              "1,275 non-independent island pairs\nUniform sea, AIS density, wind speed\nMapped current ≠ actual island relay",
              title_size=12.6)
    arrow(ax,(.29,.43),(.356,.43))
    rectangle(ax,(.716,.46),.255,.275,"#f9f5ea","Butterfly record geography",
              "Original separate 225-candidate panel\n2001–2010 / 2011–2022 indices\nHost vs later butterfly ρ = +0.300\nTwo-layer conditional predictions",
              title_size=12.6)
    arrow(ax,(.63,.78),(.70,.62),dash=True)
    arrow(ax,(.63,.43),(.70,.56),dash=True)
    ax.text(.81,.77,"Shared statistical information",ha="center",fontsize=10,color="#495765")
    ax.text(.50,.238,
            "Dashed links = retrospective association and conditional prediction, "
            "not observed plant shipments or butterfly dispersal.",
            ha="center",fontsize=11,fontweight="bold",color="#384d5e")
    ax.text(.50,.170,
            "WGSRPD3 botanical regions and WPI island groups are NOT the same ecological units. "
            "239 and 225 are not interchangeable sample denominators.",
            ha="center",fontsize=10.2)
    ax.text(.50,.111,
            "Potential resources ≠ feeding or population establishment  |  "
            "Current concentration ≠ a stopover island  |  Indexed first record ≠ arrival date",
            ha="center",fontsize=10.1)
    ax.text(.50,.055,
            "The manuscript keeps uncertainty from observation effort, chronology, island geography and maritime proxies explicit.",
            ha="center",fontsize=9.7)
    fig.subplots_adjust(left=.015,right=.99,top=.96,bottom=.035)
    checks={}
    for ext in ("png","svg","pdf"):
        p=out/f"fig1_two_layer_framework.{ext}"
        fig.savefig(p,bbox_inches="tight",dpi=180)
        checks[ext]=hashlib.sha256(p.read_bytes()).hexdigest()
    plt.close(fig)
    return checks


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out-dir",type=Path,default=ROOT/"_jbi_fig1")
    a=ap.parse_args()
    ledger=audit()
    checks=render(ledger,a.out_dir)
    receipt={"schema":"chocho_JBI_two_layer_source_anchored_figure1",
             "status":"SOURCE_ANCHORED_DIAGRAM_RENDERED",
             "evidence_entries":len(ledger["evidence"]),
             "ledger_sha256":hashlib.sha256(LEDGER.read_bytes()).hexdigest(),
             "figure_sha256":checks,
             "scope":"conceptual distinct-population diagram, not a causal directed-path model"}
    (a.out_dir/"fig1_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":
    main()
