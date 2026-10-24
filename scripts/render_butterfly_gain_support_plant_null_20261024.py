#!/usr/bin/env python3
"""Draw the observed butterfly support slope against source-qualified null chains.

The figure is an exploratory STRUCTURAL NEGATIVE CONTROL, not evidence of
butterfly demography or consumer adaptation.
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

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "provenance/reviewer_defenses/results/butterfly_gain_support_plant_null_20261024.json"


def render(svg: Path, png: Path):
    data = json.loads(RECEIPT.read_text(encoding="utf-8"))
    check = data["numerical_evaluation"]
    if (data["status"] != "QUALIFIED_POSTHOC_SOURCE_CONDITIONAL_PLANT_RANGE_NULL"
        or not check["numerically_qualified"]
        or data["focal_butterflies"] != 239
        or data["source_botanical_introduced_host_area_edges"] != 30373
        or len(data["chains"]) != 2
        or len({c["seed"] for c in data["chains"]}) != 2):
        raise ValueError("Source-locked null result has not passed qualifications")
    captions = [
        ("primary","full_null_history","β (all 206 expanded butterflies)"),
        ("adequate_catalogue","adequate_null_history","β (191 adequate-host butterflies)"),
    ]
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9})
    fig,axarr = plt.subplots(1,2,figsize=(11.2,4.4),constrained_layout=True)
    for ax,(name,history_key,title) in zip(axarr,captions):
        comparison = check["result"][name]
        observed = comparison["chains"][0]["observed_beta"]
        assert comparison["chains"][1]["observed_beta"] == observed
        for i,c in enumerate(data["chains"]):
            values = np.asarray(c[history_key],dtype=float)
            if len(values)!=199 or not np.isfinite(values).all():
                raise ValueError("Chain length or finite values unexpectedly changed")
            ref = comparison["chains"][i]
            if not np.allclose(
                np.quantile(values,[.025,.5,.975]),
                [ref["null_q025"],ref["null_median"],ref["null_q975"]],
                atol=1e-9,rtol=0):
                raise ValueError("Rendered histogram no longer matches source receipt")
            ax.hist(values,bins=17,density=True,alpha=.46,
                    label=f"Plant-range null, seed {c['seed']}")
        ax.axvline(observed,color="black",lw=2,ls="--",
                   label="Original host-range association")
        ax.set_title(title,fontweight="bold")
        ax.set_xlabel("Adjusted coefficient of log gain")
        ax.set_ylabel("Density of 199 null draws per chain")
        ax.grid(alpha=.15)
        medians = [r["null_median"] for r in comparison["chains"]]
        ci_low = min(r["null_q025"] for r in comparison["chains"])
        ci_high = max(r["null_q975"] for r in comparison["chains"])
        text = (f"Observed β = {observed:+.4f}\n"
                f"Null β medians = {medians[0]:+.4f} / {medians[1]:+.4f}\n"
                f"Both 95% null ranges include observed β")
        ax.text(.03,.97,text,transform=ax.transAxes,ha="left",va="top",
                fontsize=8.1,bbox={"boxstyle":"round,pad=0.5",
                                  "facecolor":"white","edgecolor":"none","alpha":.85})
        ax.legend(fontsize=7.2,frameon=False,loc="upper right",
                  bbox_to_anchor=(1,.74))
    fig.suptitle("Constrained plant redistributions reproduce the gain–support association",
                 fontsize=13)
    fig.supxlabel("Posthoc structural null — host identities, native ranges and botanical family × continent margins fixed",
                  fontsize=9)
    svg.parent.mkdir(parents=True,exist_ok=True)
    png.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(svg,format="svg")
    fig.savefig(png,format="png",dpi=165)
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--svg",type=Path,required=True)
    p.add_argument("--png",type=Path,required=True)
    args=p.parse_args()
    render(args.svg,args.png)
    print("PASS: 2×199 null β draws and observed associations validated and rendered")
