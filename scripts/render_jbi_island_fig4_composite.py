#!/usr/bin/env python3
"""Build the combined JBI Fig. 4 using the verified current map and source pair data.

A: PR150's fresh, verified uniform current-density image (modelled water flows).
B: all 1,275 dependent host/butterfly source-pair Jaccard values.
C: conditional prediction improvements from separately defined heldout models.

This does NOT calculate butterfly movements or stepping-stone passage.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

ROOT=Path(__file__).resolve().parents[1]
PINNED={
    "pairs":("provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_pairs.csv",
             "5663d2ffe42419f5c4496b92c9f3bbf76175907b"),
    "pair_audit":("provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_audit.json",
                  "16da0dc0f7dda124b9cce88ac81b1a63ce8c7c94"),
    "LOIO":("results/island51_host_circuit_joint_20261106/model_summary.csv",
            "1f1ac86ab086d5b1db14b330c7c65836a19c05ae"),
    "double_unseen":("results/island51_double_unseen_blocks_20261107/contrast_scope_summary.csv",
                     "528427cddda70ddc6340bc3758c0080cd3096941")
}


def blob_sha(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()


def audit_and_load():
    data, hashes={},{}
    for name,(rel,sha) in PINNED.items():
        raw=(ROOT/rel).read_bytes()
        assert blob_sha(raw)==sha,f"Source changed: {rel}"
        if rel.endswith(".json"):
            data[name]=json.loads(raw)
        else:
            data[name]=list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
        hashes[name]={"path":rel,"sha":sha}
    pairs=data["pairs"]; original=data["pair_audit"]
    assert len(pairs)==1275
    assert int(original["eligible_islands"])==51
    assert int(original["all_pairs"])==1275
    ids={x["island_a"] for x in pairs}|{x["island_b"] for x in pairs}
    assert len(ids)==51
    assert len({tuple(sorted((z["island_a"],z["island_b"]))) for z in pairs})==1275
    host=np.array([float(x["host_jaccard"]) for x in pairs])
    later=np.array([float(x["new_butterfly_jaccard"]) for x in pairs])
    earlier=np.array([float(x["earlier_butterfly_jaccard"]) for x in pairs])
    assert np.all(np.isfinite([host,later,earlier]))
    assert np.all((host>=0)&(host<=1))
    assert np.all((later>=0)&(later<=1))
    rho=np.corrcoef(rankdata(host),rankdata(later))[0,1]
    rho_early=np.corrcoef(rankdata(host),rankdata(earlier))[0,1]
    assert np.isclose(rho,original["spearman_host_new"],atol=1e-10)
    assert np.isclose(rho_early,original["spearman_host_earlier"],atol=1e-10)
    return data,hashes,host,later,earlier


def exact_heldout_scores(data):
    m={}
    for row in data["LOIO"]:
        if row["target"]=="new_butterfly_jaccard":
            m[row["model"]]=float(row["MSE"])
            assert int(row["n_folds"])==51
    assert set(["H","C","HC"]).issubset(m)
    host_score=100*(m["C"]-m["HC"])/m["C"]
    circuit_score=100*(m["H"]-m["HC"])/m["H"]
    assert np.isclose(host_score,1.7312,atol=0.0001)
    assert np.isclose(circuit_score,6.4880,atol=0.0001)

    strict={}
    for row in data["double_unseen"]:
        if row["target"]!="new_butterfly_jaccard" or row["pair_class"]!="all_pairs":
            continue
        if row["role"] not in ["host_composition_beyond_circuits","circuits_beyond_host_composition"]:
            continue
        strict[row["role"]]=float(row["MSE_reduction_percent"])
        assert int(row["unique_evaluated_pairs"])==1275
    assert len(strict)==2
    assert np.isclose(strict["host_composition_beyond_circuits"],0.8808,atol=.0001)
    assert np.isclose(strict["circuits_beyond_host_composition"],4.4744,atol=.0001)
    return [host_score,circuit_score], [
        strict["host_composition_beyond_circuits"],
        strict["circuits_beyond_host_composition"]
    ]


def verify_map(current_map_dir:Path):
    image=current_map_dir/"figure4_uniform_marine_current.png"
    receipt=current_map_dir/"execution_receipt.json"
    if not image.is_file() or not receipt.is_file():
        raise FileNotFoundError("Run original source-verified current map generator first")
    report=json.loads(receipt.read_text(encoding="utf-8"))
    assert report["schema"]=="chocho_island51_overwater_pair_current_map_v1"
    assert report["conductivity"]=="uniform"
    assert (report["islands"],report["island_pairs"])==(51,1275)
    assert (report["marine_nodes"],report["marine_edges"])==(189550,370996)
    assert report["max_linear_residual"] < 1e-7
    assert len(report["archived_resistance_energy_checks"])>=3
    return image,report


def render(data,host,later,earlier,mse_loio,mse_both,map_image,out):
    fig=plt.figure(figsize=(17.4,12.2),constrained_layout=False)
    grid=fig.add_gridspec(2,2,height_ratios=[1.30,0.78],width_ratios=[1.18,0.82],
                         left=.065,right=.965,top=.925,bottom=.115,
                         hspace=.24,wspace=.40)
    top=fig.add_subplot(grid[0,:])
    sc=fig.add_subplot(grid[1,0])
    bx=fig.add_subplot(grid[1,1])
    fig.suptitle("Resources and marine accessibility contain distinct geographic information",
                 fontsize=17,weight="bold",y=.983)

    image=plt.imread(map_image)
    # Reuse the inspected coastline/current panel but trim its own title and
    # bottom narrative before it is embedded under this figure's panel title.
    # Keep the original axes, centroids and colourbar.
    image=image[75:-45,:,:]
    top.imshow(image,aspect="equal")
    top.set_axis_off()
    top.set_title("A  Modelled over-water current concentration among 51 island groups",
                  loc="left",fontsize=13,pad=6)

    # The original pair table is strongly concentrated at zero Jaccard.
    # Ordinary linear Blues maps most nonempty source bins nearly white.
    # A logarithmic *count* colour scale and a non-white low end preserve all
    # 1,275 original observations without changing correlations or models.
    count_cmap=LinearSegmentedColormap.from_list(
        "island_pair_counts",["#96bedb","#3d7da9","#073a6e"])
    hb=sc.hexbin(host,later,gridsize=32,mincnt=1,cmap=count_cmap,
                 norm=LogNorm(vmin=1),extent=(0,1,0,1),linewidths=.08)
    colorbar=fig.colorbar(hb,ax=sc,fraction=.040,pad=.024,shrink=.83)
    colorbar.set_label("Island pairs / hexagon (log count)",fontsize=8.2)
    colorbar.ax.tick_params(labelsize=7.5)
    sc.set_xlim(-.02,1.02);sc.set_ylim(-.02,1.02)
    sc.set_xlabel("Jaccard similarity of introduced candidate hosts")
    sc.set_ylabel("Jaccard similarity of later-indexed butterflies")
    sc.set_title("B  Host–butterfly community correspondence\n51 islands, 1,275 dependent pairs",
                 loc="left",fontsize=12)
    sc.text(.98,.98, f"Later Spearman ρ = {data['pair_audit']['spearman_host_new']:+.3f}\n"
            f"Earlier ρ = {data['pair_audit']['spearman_host_earlier']:+.3f}\n"
            f"Partial later ρ = {data['pair_audit']['partial_rank_host_new']:+.3f}",
            transform=sc.transAxes,ha="right",va="top",fontsize=9.4,
            bbox=dict(facecolor="white",alpha=.9,edgecolor="none",pad=6))

    xs=np.arange(2)
    yy=["Hosts beyond marine terms","Marine terms beyond hosts"]
    width=.32
    bx.barh(xs+.19,mse_loio,height=width,label="One island unseen",
            color="#577d9d")
    bx.barh(xs-.19,mse_both,height=width,label="Both islands unseen",
            color="#d59d55")
    for y,val in zip(xs+.19,mse_loio):
        bx.text(val+.11,y,f"+{val:.2f}%",va="center",fontsize=9.6)
    for y,val in zip(xs-.19,mse_both):
        bx.text(val+.11,y,f"+{val:.2f}%",va="center",fontsize=9.6)
    bx.set_yticks(xs,yy)
    bx.tick_params(axis="y",labelsize=9.5)
    bx.set_xlim(0,8.4);bx.set_ylim(-.62,1.8)
    bx.set_xlabel("Decrease in held-out MSE (%)")
    bx.set_title("C  Conditional predictive information\nLater butterfly index similarity",
                 loc="left",fontsize=12)
    bx.legend(loc="lower right",frameon=False,fontsize=8.5)
    bx.grid(axis="x",alpha=.15)
    bx.set_axisbelow(True)
    fig.text(.5,.063,
             "Different held-out splits use different covariates and MSE baselines. "
             "Hosts indexed in 2001–2024 can postdate 2011–2022 butterfly records.",
             ha="center",fontsize=9)
    fig.text(.5,.037,
             "Corridors cannot identify stepping-stone islands (no modeled through-island transfer); "
             "prediction gains do not demonstrate transport, feeding or colonisation.",
             ha="center",fontsize=9)

    out.mkdir(parents=True,exist_ok=True)
    results={}
    for ext in ("png","svg","pdf"):
        path=out/f"fig4_resource_and_marine_circuit.{ext}"
        fig.savefig(path,bbox_inches="tight",dpi=185)
        results[ext]=hashlib.sha256(path.read_bytes()).hexdigest()
    plt.close(fig)
    return results


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--current-map-dir",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    args=ap.parse_args()
    data,sources,h,late,early=audit_and_load()
    scores,strict=exact_heldout_scores(data)
    image,current=verify_map(args.current_map_dir)
    exports=render(data,h,late,early,scores,strict,image,args.out_dir)
    receipt={
        "schema":"chocho_JBI_Fig4_full_two_layer_composite_v1",
        "status":"SOURCE_AUDITED_COMPOSITE_RENDERED",
        "source_git_blobs":sources,
        "source_current_run":{"islands":current["islands"],"pairs":current["island_pairs"],
                              "max_linear_residual":current["max_linear_residual"]},
        "later_spearman":data["pair_audit"]["spearman_host_new"],
        "earlier_spearman":data["pair_audit"]["spearman_host_earlier"],
        "partial_later":data["pair_audit"]["partial_rank_host_new"],
        "MSE_gain_LOIO":scores,
        "MSE_gain_both_unseen":strict,
        "exports_sha256":exports,
        "not_identified":["ship effects","island relay transport","host-before-butterfly chronological order"]
    }
    (args.out_dir/"fig4_composite_receipt.json").write_text(
        json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":receipt["status"],"sources":len(sources),"exports":exports},indent=2))


if __name__=="__main__":
    main()
