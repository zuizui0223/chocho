"""Render the two main figures from the archived regional and port result tables."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures/integrated_island_port_butterfly_v03"
GAPS = ROOT / "results/island_resource_gap_filling"
PORT = ROOT / "provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_top10_pre2011_port_plant_overlap.json"
GREEN, SLATE = "#126a59", "#546d81"


def save(fig, stem):
    for ext, meta in [("png", {}), ("pdf", {"CreationDate": None}), ("svg", {"Date": None})]:
        path = OUT / f"{stem}.{ext}"
        fig.savefig(path, dpi=225, metadata=meta)
        if ext == "svg":
            lines = path.read_text(encoding="utf8").splitlines()
            path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf8")
    plt.close(fig)


def regional():
    with (GAPS / "pooled.csv").open() as stream:
        data = {r["geography"]: r for r in csv.DictReader(stream)}
    matched = json.loads((GAPS / "summary.json").read_text())
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.95))
    fig.subplots_adjust(left=.058, right=.98, top=.68, bottom=.23, wspace=.34)
    fig.text(.04,.92,"Why did potential larval-host coverage grow more on islands?",size=17,weight="bold")
    fig.text(.04,.81,"239 butterfly species   ·   126 island and 237 mainland botanical regions",size=11)
    for ax, key, title, divisor, limit in zip(axes, ["coverage","fill","gain"],
            ["A   Native host coverage","B   Empty regions filled","C   Gain over native coverage"],
            ["All species–region pairs","Pairs without a native host","Pairs with a native host"],[50,50,145]):
        vals = [100*float(data[g][key]) for g in ("island","mainland")]
        ax.bar([0,1],vals,width=.56,color=[GREEN,SLATE])
        ax.set(xticks=[0,1],xticklabels=["Island","Mainland"],ylim=(0,limit),ylabel="%")
        ax.set_title(title,loc="left",fontsize=11,weight="bold")
        ax.set_xlabel(divisor,fontsize=9)
        for i,v in enumerate(vals):
            ax.text(i,v+limit*.035,f"{v:.1f}%",ha="center",fontsize=11,weight="bold")
        ax.spines[["top","right"]].set_visible(False)
    fig.text(.04,.04,f"After matching area and broad geography, introduced hosts filled {matched['matched_island_fill_pct']:.1f}% of island gaps and {matched['matched_mainland_fill_pct']:.1f}% of mainland gaps.",size=10)
    save(fig,"Fig1_resource_geography_scales_v03")


def ports():
    src = json.loads(PORT.read_text())
    taxa = sorted(src["source_contributor_plants"],key=lambda r: -r["GEB_fraction_of_added_species_region_opportunity_attributed_to_this_plant"])
    fig,(a,b) = plt.subplots(1,2,figsize=(9,5.5),sharey=True)
    fig.subplots_adjust(left=.25,right=.95,top=.75,bottom=.18,wspace=.15)
    fig.text(.04,.93,"Which globally important host plants were recorded near island ports?",size=15,weight="bold")
    fig.text(.04,.85,f"{src['top10_species_with_plant_index_at_study_island_ports']} of the {len(taxa)} leading contributors occurred in the 2001–2010 inventory of {src['study_port_proxies']} ports on {src['study_islands']} landmasses.",size=10)
    for ax,vals,title,xlabel,xlim,color in [
        (a,[100*r["GEB_fraction_of_added_species_region_opportunity_attributed_to_this_plant"] for r in taxa],"A   Global resource gain","Share of added species–region pairs (%)",5.4,SLATE),
        (b,[r["n_port_buffers_with_pre2011_record"] for r in taxa],"B   Island port inventories","Port buffers with a record (of 38)",max(r["n_port_buffers_with_pre2011_record"] for r in taxa)+2,GREEN)]:
        ax.barh(range(len(taxa)),vals,height=.48,color=color)
        ax.set_title(title,loc="left",fontsize=11,weight="bold")
        ax.set_xlabel(xlabel,fontsize=9)
        ax.set_xlim(0,xlim)
        ax.spines[["top","right","left"]].set_visible(False)
        ax.tick_params(axis="y",length=0)
        for n,v in enumerate(vals):
            ax.text(v+xlim*.025,n,f"{v:.2f}" if ax is a else str(v),va="center",fontsize=9)
    a.set_yticks(range(len(taxa)),[r["plant_name"] for r in taxa],fontstyle="italic")
    a.invert_yaxis()
    fig.text(.04,.055,"Rows link accepted plant names across datasets. Zero denotes no record in the sampled ports.",size=10)
    save(fig,"Fig2_global_contributors_at_ports")


if __name__ == "__main__":
    plt.rcParams.update({"svg.fonttype":"none","svg.hashsalt":"chocho-main","pdf.fonttype":42})
    regional()
    ports()
