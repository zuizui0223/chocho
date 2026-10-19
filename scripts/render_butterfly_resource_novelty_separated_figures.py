#!/usr/bin/env python3
"""Generate manuscript-aligned, separate resource-accounting Fig1 and null FigS2.

The frozen, source-checked numerical data come from the already-audited
butterfly_resource_novelty_figure.extract_results. No new ecological response
or simulated null is estimated by this figure-only program.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from render_butterfly_resource_novelty_figure import extract_results, load


def split_figure_data(data: dict) -> tuple[dict, dict]:
    if (
        data["focal_butterflies"] != 239
        or data["resource_regions"] != 365
        or data["total_candidate_incidences"] != 75346
        or [v["n"] for v in data["accounting"]] != [46387, 14406, 14553]
        or len(data["controls"]) != 3
        or [m["draws_per_chain"] for m in data["controls"]] != [499, 299, 99]
    ):
        raise ValueError("source verified accounting/structural sensitivity input drift")
    if sum(v["n"] for v in data["accounting"]) != data["total_candidate_incidences"]:
        raise ValueError("potential resource-incidence accounting not conserved")
    control = data["controls"]
    if any(
        len(x[k]) != 2 or any(v <= 0 for v in x[k])
        for x in control for k in ("a_excess", "u_deficit")
    ):
        raise ValueError("missing or invalid dual-chain null outcome")
    main = {
        "schema": "chocho_butterfly_fig1_geographic_novelty_accounting_v0.1",
        "focal_butterflies":239,
        "wgsrpd3_resource_regions":365,
        "native_butterfly_region_cells":26530,
        "contemporary_butterfly_region_cells":41083,
        "novel_butterfly_region_cells":14553,
        "incidence_budget":75346,
        "accounting":data["accounting"],
        "identity": "H = A + R + U; each novel butterfly x region cell receives one first support incidence",
        "claim_boundary": "A, R, U are host x butterfly x region candidate-incidence units, NOT relative butterfly abundance or probability of local host feeding.",
    }
    supplementary = {
        "schema": "chocho_butterfly_figure_S2_matched_null_v0.1",
        "incidence_budget":75346,
        "observed_native_absorbed_incidence_count":46387,
        "observed_distinct_new_butterfly_region_cells":14553,
        "controls":control,
        "claim_boundary": "Post-hoc, source-conditional randomization results; bars are mean of two chain medians, NOT confidence intervals, nor causal variance partition, nor independent exact-permutation p-values.",
    }
    return main, supplementary


def _save(fig, png: Path, svg: Path):
    png.parent.mkdir(parents=True, exist_ok=True)
    svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")


def render_main1(data:dict, png:Path, svg:Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n=data["native_butterfly_region_cells"]
    x=data["novel_butterfly_region_cells"]
    total=data["contemporary_butterfly_region_cells"]
    H=data["incidence_budget"]
    if n+x != total:
        raise ValueError("native + introduced-added geographic union identity drift")

    fig=plt.figure(figsize=(12.5, 6.25), facecolor="white")
    fig.subplots_adjust(top=.77,bottom=.17,left=.105,right=.93)
    fig.suptitle(
        "Most mapped introduced host supply does not open a distinct resource region",
        fontsize=16, weight="semibold", x=.10, ha="left", y=.96
    )
    fig.text(
        .10,.89,
        "239 butterfly species  •  365 WGSRPD3 regions  •  fixed butterfly–host associations",
        fontsize=10.5
    )
    fig.text(
        .10,.82,
        "Potential butterfly × region resource cells:   "
        f"{n:,} native → {total:,} contemporary (+{x:,}; +{100*x/n:.1f}%)",
        fontsize=11,weight="medium"
    )

    ax=fig.add_axes([.105,.42,.82,.27])
    color=plt.rcParams["axes.prop_cycle"].by_key()["color"][:3]
    labels=[
        "A: Already native-supported",
        "R: Extra introduced host support",
        "U: First support for a new cell",
    ]
    left=0
    handles=[]
    for i,item in enumerate(data["accounting"]):
        value=item["n"]
        share=100*value/H
        bar=ax.barh(0,share,left=left,color=color[i],height=.62,linewidth=0)
        handles.append(bar[0])
        ax.text(left+share/2, 0, f"{value:,}\n{share:.1f}%",
                ha="center",va="center",fontsize=12,
                weight="semibold", color="white")
        left+=share
    ax.set_xlim(0,100);ax.set_ylim(-.59,.59)
    ax.set_yticks([])
    ax.set_xticks([0,20,40,60,80,100])
    ax.set_xlabel(
        "Share of 75,346 introduced host × butterfly × region candidate incidences",
        fontsize=10.5
    )
    ax.spines[["top","left","right"]].set_visible(False)

    fig.legend(handles,labels,bbox_to_anchor=(.095,.325),
        loc="upper left",ncol=3,frameon=False,fontsize=9)
    fig.text(.105,.215,
        "A + R + U = H  |  Each added butterfly × region cell contributes "
        "ONE first-support incidence (U); other introduced hosts in that cell are R.",
        fontsize=10)
    fig.text(.105,.120,
        "These are POTENTIAL larval host-resource maps. Incidence fractions are "
        "not butterfly demographic gains, host biomass, local feeding, or "
        "colonisation rates.",fontsize=9)
    _save(fig,png,svg)
    plt.close(fig)


def render_s2(data:dict,png:Path,svg:Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    controls=data["controls"]
    fig,axs=plt.subplots(1,2,figsize=(12.4,6.0))
    fig.subplots_adjust(top=.72,bottom=.23,left=.10,right=.96,wspace=.31)
    fig.suptitle("Small structural deviations persist under geographic constraints",
                 fontsize=15.5,weight="semibold",x=.10,ha="left",y=.96)
    fig.text(.10,.87,
        "The exact original-HOSTS incidence budget (H=75,346) is fixed "
        "within each plant-range null; models differ in their permitted swaps.",
        fontsize=10)
    xy=range(3)
    xt=["Family +\ncontinent", "+ Coarse climate\nBIO1 × BIO12",
        "+ Local geography\n≤1,500 km"]
    color=plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for ax,key,title,ylab in (
        (axs[0],"a_excess","A  |  More native overlap than expected",
         "Observed minus null median\n(host×butterfly×region incidences)"),
        (axs[1],"u_deficit","B  |  Fewer distinct new cells than expected",
         "Null median minus observed\n(distinct butterfly×region cells)"),
    ):
        med=[sum(v[key])/len(v[key]) for v in controls]
        ax.bar(xy,med,color=color[0],alpha=.45,width=.54)
        for k,c in enumerate(controls):
            xx=[k-.08,k+.08]
            ax.scatter(xx,c[key],s=35,color=color[1],zorder=3)
            ax.annotate(
                f"{min(c[key]):.0f}–{max(c[key]):.0f}",
                (k,max(c[key])),ha="center",textcoords="offset points",
                xytext=(0,7),fontsize=9
            )
        ax.set_ylim(0,max(med)*1.22)
        ax.set_xticks(list(xy),xt)
        ax.set_ylabel(ylab,fontsize=9.2)
        ax.set_title(title,fontsize=10.5,pad=17)
        ax.grid(axis="y",alpha=.13)
        ax.set_axisbelow(True)
        ax.spines[["top","right"]].set_visible(False)

    fig.text(
        .10,.12,
        "Bars = mean of two chain medians; dots = each chain's null median. "
        "499, 299 and 99 simulated worlds per chain, respectively.",fontsize=9
    )
    fig.text(
        .10,.065,
        "This POST-HOC comparison is source-conditional and not a causal "
        "fraction 'explained by' climate or distance; the local residual is "
        "≈0.36 percentage points of H.",fontsize=8.9
    )
    _save(fig,png,svg)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--plant-null-json",type=Path,required=True)
    p.add_argument("--climate-null-json",type=Path,required=True)
    p.add_argument("--local-null-json",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    src=extract_results(load(a.plant_null_json),load(a.climate_null_json),load(a.local_null_json))
    main_fig,supp=split_figure_data(src)
    a.output_dir.mkdir(parents=True,exist_ok=True)
    render_main1(main_fig,
        a.output_dir/"figure_main1_geographic_accounting.png",
        a.output_dir/"figure_main1_geographic_accounting.svg")
    render_s2(supp,
        a.output_dir/"figure_S2_matched_null_sensitivity.png",
        a.output_dir/"figure_S2_matched_null_sensitivity.svg")
    for name,obj in (
        ("figure_main1_geographic_accounting_values.json",main_fig),
        ("figure_S2_matched_null_sensitivity_values.json",supp),
    ):
        (a.output_dir/name).write_text(
            json.dumps(obj,indent=2,ensure_ascii=False)+"\n",encoding="utf-8"
        )
    print(json.dumps({
        "status":"RENDERED_SEPARATE_EXPLORATORY_FIG1_AND_FIGS2",
        "main_figure":main_fig["schema"],
        "supplementary_figure":supp["schema"],
        "files":sorted(p.name for p in a.output_dir.iterdir()),
    },indent=2))


if __name__=="__main__":
    main()
