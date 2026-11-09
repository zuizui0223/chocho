#!/usr/bin/env python3
"""Source-linked new Figure 5 for ONE islands/ports/butterflies GEB article v0.3.

Uses existing original source-locked CSV/JSON results. Never fits a model,
inventing transport tracks or interpreting GBIF first index as colonisation.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

R=Path(__file__).resolve().parents[1]
OUT=R/"figures/integrated_island_port_butterfly_v03"
PATHS={
    "joint":"results/island51_host_circuit_joint_20261106/comparison_summary.csv",
    "historical_cross":"results/island51_cross_community_marine_20261107/comparison_scope_summary.csv",
    "early_host":"provenance/reviewer_defenses/results/original230_dated_host_join_20261031/original230_prior_host_recurrence_audit.json",
    "host_null":"provenance/reviewer_defenses/results/host_geographic_prevalence_null_20261031/host_geo_frequency_null_results.json",
    "public51":"data/exploratory/public51_butterfly_index_20261109/full51_source_validation.json",
    "firstidx":"results/public51_butterfly_trait_access_20261109/public51_butterfly_conditional_trait_marine.json",
    "firsttransfer":"results/public51_butterfly_trait_access_20261109/public51_trait_summary.json",
}
def read(k):
    p=R/PATHS[k]
    if p.suffix==".json":return json.loads(p.read_text(encoding="utf8"))
    with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def numbers():
    joint=read("joint");cross=read("historical_cross");early=read("early_host")
    null=read("host_null");g=read("public51");idx=read("firstidx");t=read("firsttransfer")
    def j(c):
        x=[float(z["relative_MSE_reduction_percent"]) for z in joint
           if z["target"]=="new_butterfly_jaccard" and z["comparison"]==c]
        if len(x)!=1:raise ValueError("joint source comparator absent "+c)
        return x[0]
    def c(panel,contrast):
        a=[z for z in cross if z["panel"]==panel and z["contrast"]==contrast
             and z["pair_class"]=="all_both_plants_indexed"]
        if len(a)!=1 or int(a[0]["unique_source_pairs"])!=990:raise ValueError("cross community source mismatched")
        return float(a[0]["MSE_gain_percent"])
    vals={
       "host_given_marine":j("C_to_HC"),
       "marine_given_host":j("H_to_HC"),
       "base_to_hosts":j("G_to_H"),
       "base_to_marine":j("G_to_C"),
       "flora_ais":c("early_general_flora","BASE_to_AIS"),
       "early_butterflies_ais":c("earlier_butterflies","BASE_to_AIS"),
       "later_butterflies_ais":c("later_butterflies","BASE_to_AIS"),
       "later_butterflies_ais_with_host":c("later_butterflies_retro_hosts","BASE_to_AIS"),
       "later_butterflies_ais_wind_with_host":c("later_butterflies_retro_hosts","BASE_to_BOTH"),
       "prior_host_share":early["host_network_structural_null"]["observed_prior_host_index_event_fraction"]*100,
       "family_null":early["host_network_structural_null"]["null_fraction_mean"]*100,
       "geo_null":null["more_constrained_structural_nulls"]["family_plus_region_signature_and_prevalence_band"]["null_mean_fraction"]*100,
       "first_pre":g["totals"]["first_indexed_in_2011_2017"],
       "first_late":g["totals"]["first_indexed_in_2018_2022"],
       "sea_first_transfer":t["periods"]["2018_2022"]["heldout"]["geographic_cluster"]["sea_over_geodesic"]["brier"],
    }
    # Use source's own nearest conspecific and residual-ocean fixed-margin
    tests=idx["periods"]["2018_2022"]["within_geographic_cluster_two_margin"]["tests"]
    vals["source_proximity_p"]=float(tests["geo_access_mean"]["p_two_sided_monte_carlo"])
    vals["residual_sea_p"]=float(tests["sea_residual_access_mean"]["p_two_sided_monte_carlo"])
    check={"host_given_marine":1.731162479353141,
           "marine_given_host":6.4880250544198574,
           "later_butterflies_ais_with_host":.67755574674502794,
           "later_butterflies_ais_wind_with_host":-4.6985954047498986,
           "sea_first_transfer":-6.479351063178629}
    for k,want in check.items():
        if abs(vals[k]-want)>1e-8:raise ValueError(f"Unfrozen numeric source {k}: {vals[k]} != {want}")
    if vals["first_pre"]!=221 or vals["first_late"]!=96:raise ValueError("GBIF source cohort changed")
    if abs(vals["prior_host_share"]-34/113*100)>1e-9:raise ValueError("historical host cohort changed")
    return vals

def build(v,output):
    fig,ax=plt.subplots(2,2,figsize=(14.8,9.3),layout="constrained")
    fig.suptitle("Figure 5  |  Host and marine geography: complementary records, uncertain pathways",
                 x=.02,ha="left",fontsize=15.2,fontweight="bold")
    # A: conditional same-fold held-island prediction.
    a=ax[0,0]; y=np.arange(3)
    data=[v["host_given_marine"],v["marine_given_host"],8.6164]
    labels=["Add hosts to AIS/wind\n(C → HC)","Add AIS/wind to hosts\n(H → HC)","Both vs base geography\n(G → HC; contextual)"]
    bars=a.barh(y,data,color=["#529b83","#2c6c8e","#8b929a"])
    a.invert_yaxis();a.set_yticks(y,labels,fontsize=8.6)
    a.set_xlim(0,10.3);a.set_xlabel("Relative held-island MSE reduction (%)")
    for q,x in zip(bars,data):
        a.text(x+.15,q.get_y()+q.get_height()/2,f"+{x:.2f}%",va="center",fontsize=10,fontweight="bold")
    a.set_title("A  51 islands: complementary conditional predictors",loc="left",fontsize=11)
    a.text(.01,.02,"Same leave-one-island-out folds; host lists visible on the test island.\nComparisons have different denominators—do not sum gains.",
           transform=a.transAxes,va="bottom",fontsize=8.0,color="#536070")

    # B: stricter independent double-unseen geography.
    b=ax[0,1];z=[v["flora_ais"],v["early_butterflies_ais"],v["later_butterflies_ais"],v["later_butterflies_ais_with_host"],v["later_butterflies_ais_wind_with_host"]]
    ys=np.arange(len(z));labs=["Earlier general plants: +AIS","Earlier butterflies: +AIS","Later butterflies: +AIS","Later butterflies, hosts held: +AIS","Later butterflies, hosts held: +AIS + wind"]
    b.barh(ys,z,color=["#a58c8b","#a58c8b","#598f9b","#598f9b","#a58c8b"])
    b.axvline(0,color="#1f2937",lw=.9);b.set_yticks(ys,labs,fontsize=8.0);b.invert_yaxis()
    b.set_xlim(-10.5,3);b.set_xlabel("MSE reduction (%) — negative means worse")
    for i,vv in enumerate(z):b.text(vv+(.23 if vv>=0 else -.23),i,f"{vv:+.2f}%",ha="left" if vv>=0 else "right",fontsize=9,va="center")
    b.set_title("B  990 pairs: both endpoints geographically unseen",loc="left",fontsize=11)
    b.text(.01,.02,"Modern AIS/wind do not identify transport paths or historical exposures.",
           transform=b.transAxes,va="bottom",fontsize=8,color="#536070")

    # C: independently sourced species x island first-index response
    c=ax[1,0]
    cc=c.bar([0,1],[v["first_pre"],v["first_late"]],width=.56,color=["#5c94ae","#8ca8b9"])
    c.set_xticks([0,1],["2011–2017","2018–2022"])
    c.set_ylabel("First-GBIF-index butterfly × island cells")
    c.set_ylim(0,265)
    for q,yy in zip(cc,[v["first_pre"],v["first_late"]]):
        c.text(q.get_x()+q.get_width()/2,yy+5,str(yy),ha="center",fontsize=11)
    c.set_title("C  New indexes and limits of marine transfer",loc="left",fontsize=11)
    c.text(.025,.89,"Late conditional source proximity p=%.5f\nLate residual ocean resistance p=%.5f\nBUT unseen-region uniform-sea model: %.2f%% Brier MSE"%
           (v["source_proximity_p"],v["residual_sea_p"],v["sea_first_transfer"]),
           transform=c.transAxes,va="top",fontsize=9.4,
           bbox=dict(boxstyle="round",facecolor="white",alpha=.90,edgecolor="#d5dce4"))
    c.text(.025,.04,"p: constrained retrospective 83-event subset, not flight paths.",transform=c.transAxes,
           va="bottom",fontsize=8,color="#536070")

    # D: historical prior host null explicitly retained.
    d=ax[1,1]
    xs=np.arange(3);vv=[v["prior_host_share"],v["family_null"],v["geo_null"]]
    bs=d.bar(xs,vv,width=.62,color=["#3d8b83","#9ab4ac","#a3adb5"])
    d.set_xticks(xs,["Observed\n34/113","Plant family/\ndegree null","Plant-geography/\nprevalence null"],fontsize=8.6)
    d.set_ylabel("First-index butterfly cells with prior host index (%)")
    d.set_ylim(0,40)
    for q,yy in zip(bs,vv):d.text(q.get_x()+q.get_width()/2,yy+.7,f"{yy:.1f}%",ha="center",fontsize=10)
    d.set_title("D  Prior host indexing is not exceptional under geography",loc="left",fontsize=11)
    d.text(.02,.96,"2001–2010 plant index → 2011–2017 butterfly index\n17 original landmasses; local use and establishment untested",
           transform=d.transAxes,va="top",fontsize=8.6,
           bbox=dict(boxstyle="round",facecolor="white",alpha=.92,edgecolor="#d5dce4"))
    for a in ax.flat:
        a.spines[["top","right"]].set_visible(False)
        a.grid(axis="y",alpha=.10)
    fig.text(.02,.018,"All panels concern retrospective GBIF indexing. Distinct sampling units and holdouts are NOT pooled; electrical resistance, AIS density and wind do not show realized migration.",fontsize=8,color="#526071")
    output.mkdir(parents=True,exist_ok=True)
    prefix="Fig5_hosts_marine_conditional_vs_transfer_v03"
    for ext in ("png","svg"):
        fig.savefig(output/(prefix+"."+ext),bbox_inches="tight",dpi=260 if ext=="png" else None)
    plt.close(fig)
    return prefix

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",type=Path,default=OUT)
    args=p.parse_args()
    vals=numbers()
    name=build(vals,args.out_dir)
    digest={k:{
      "sha256":hashlib.sha256((R/v).read_bytes()).hexdigest(),
      "git_blob":hashlib.sha1(b"blob "+str((R/v).stat().st_size).encode()+b"\0"+(R/v).read_bytes()).hexdigest()
     } for k,v in PATHS.items()}
    artifacts={}
    for ext in ("png","svg"):
        f=args.out_dir/(name+"."+ext)
        artifacts[f.name]={"bytes":f.stat().st_size,"sha256":hashlib.sha256(f.read_bytes()).hexdigest()}
    result={"schema":"chocho_integrated_butterfly_figure5_v03","status":"SOURCE_LOCKED_FIGURE_RENDERED",
      "values":vals,"sources":PATHS,"digests":digest,"figures":artifacts,
      "limitations":["No causal access or realized colonisation","GBIF indexing not arrival",
         "Host lists pooled 2001-2024 may postdate indexed butterfly response",
         "AIS source modern; no historical vessel itineraries",
         "Held-island and both-endpoint unseen comparisons are distinct population and forecast designs"]}
    receipt=args.out_dir/"FIGURE5_V03_SOURCE_RECEIPT.json"
    receipt.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("FIG5_V03_SOURCE_LOCK_PASS",json.dumps({"value":vals,"artifacts":artifacts},sort_keys=True))
if __name__=="__main__":main()
