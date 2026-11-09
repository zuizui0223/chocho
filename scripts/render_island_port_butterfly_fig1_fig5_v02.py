#!/usr/bin/env python3
"""Render integrated Fig1 + Fig5: ports, island introduced hosts, butterflies.

Source-locked butterfly-centred plots and a strictly noncausal concept chart.
Works on public committed summaries, never private occurrence locations or
unobserved host-transport trajectories. Current marine currents not portrayed
as actual island-to-island movement. Outputs PNG, SVG, PDF and SHA receipt.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.ticker import MaxNLocator
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DEFAULT=ROOT/"figures/integrated_island_port_butterfly_v02"
SOURCE={
    "island":"results/island-mainland-resource-q1/island_mainland_summary.json",
    "global":"benchmarks/exploratory/butterfly_anthropogenic_resource_expansion_result_v0.1.json",
    "port":"provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_validation.json",
    "bridge":"provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_to_port_species_bridge.json",
    "pair":"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_audit.json",
    "early":"provenance/reviewer_defenses/results/original230_dated_host_join_20261031/original230_prior_host_recurrence_audit.json",
    "plant_null":"provenance/reviewer_defenses/results/host_geographic_prevalence_null_20261031/host_geo_frequency_null_results.json",
    "access":"results/public51_butterfly_trait_access_20261109/public51_butterfly_conditional_trait_marine.json",
    "oof":"results/public51_butterfly_trait_access_20261109/public51_trait_summary.json",
    "butterfly":"results/public51_butterfly_species_index_profiles_20261109/butterfly_species_indexed_distance_profiles.csv",
    "migrant49":"results/public51_author_selected_migrant_access_20261109/author_selected_migrant_no9_15_source_null.json",
    "migrant51":"results/public51_author_selected_migrant_access_20261109/author_selected_migrant_two_margin_index_null.json",
    "gbif51":"data/exploratory/public51_butterfly_index_20261109/full51_source_validation.json",
}
def read(key):
    p=ROOT/SOURCE[key]
    if p.suffix==".csv":
        with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
    return json.loads(p.read_text(encoding="utf-8"))

def git_blobs():
    return {name:hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
            for name,p in SOURCE.items() for b in [(ROOT/p).read_bytes()]}

def save(fig,out,stem,caption):
    fig.text(.028,.021,caption,fontsize=8.3,color="#4b5563",ha="left",va="bottom")
    for fmt in ("png","svg","pdf"):
        kw={"dpi":230} if fmt=="png" else {}
        fig.savefig(out/(stem+"."+fmt),bbox_inches="tight",**kw)
    plt.close(fig)

def box(ax,xy,w,h,heading,content,fc="#edf7f6",ec="#4a8691"):
    x,y=xy
    rect=FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.014,rounding_size=0.035",
                        linewidth=1.6,edgecolor=ec,facecolor=fc)
    ax.add_patch(rect)
    ax.text(x+.026,y+h-.055,heading,fontsize=12,fontweight="bold",
            color="#123447",va="top")
    ax.text(x+.026,y+h-.145,content,fontsize=10.2,color="#173247",
            va="top",linespacing=1.36)
    return (x,y,w,h)

def arrow(ax,a,b,label="",dashed=False,rad=0,dy=0):
    p=FancyArrowPatch(a,b,arrowstyle="-|>",mutation_scale=18,lw=2.2,
                      color="#267b83" if not dashed else "#727984",
                      linestyle="-" if not dashed else (0,(5,4)),
                      connectionstyle=f"arc3,rad={rad}")
    ax.add_patch(p)
    if label:
        x=(a[0]+b[0])/2; y=(a[1]+b[1])/2+dy
        ax.text(x,y,label,fontsize=9,ha="center",va="bottom",
                color="#155461" if not dashed else "#616772",
                bbox=dict(facecolor="white",edgecolor="none",alpha=.86,pad=1.8))

def fig1(out):
    ins=read("island")["groups"]["island"]
    resource=read("global")["full_resource_eligible_panel"]
    port=read("port")
    bridge=read("bridge")["aligned_port_candidate_comparison"]
    pair=read("pair")
    assert ins["contemporary_butterfly_region_units"]==10570
    assert resource["species"]==239 and resource["total_introduced_added_species_units"]==14553
    assert port["complete_ports"]==1596 and port["alien_candidate_species"]==687
    assert abs(bridge["spearman_global_added_units_vs_count_of_island_ports_with_candidate_hosts"]-.65242824936)<1e-8
    assert abs(pair["spearman_host_new"]-.29980576483)<1e-8

    fig,ax=plt.subplots(figsize=(14.8,7.8))
    fig.suptitle("Figure 1 | Islands and ports as windows onto butterfly resource geography",
                 x=.025,y=.986,ha="left",fontsize=17,fontweight="bold",color="#14263a")
    ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
    # Boxes are empirical/operational layers, not an observed sequence of organisms transported.
    box(ax,(.025,.56),.27,.23,"KNOWN BUTTERFLY–HOST LINKS",
        "239 butterfly species\nHost list held fixed\nNative vs introduced WCVP ranges",
        fc="#eef1ff",ec="#5f70a8")
    box(ax,(.385,.56),.27,.23,"POTENTIAL RESOURCES",
        "Island +118.9%; mainland +40.9%\nGlobal +54.9% opportunity\n206 / 239 butterfly species gained")
    box(ax,(.745,.56),.23,.23,"BUTTERFLY RECORDS",
        "155-species prediction\n225 × 51 landmass indexing\nNot observed colonisation",
        fc="#eff8e9",ec="#729052")
    arrow(ax,(.305,.68),(.373,.68),"Defined resource\nconstruction",dy=.02)
    arrow(ax,(.665,.68),(.733,.68),"Observed record\ncorrespondence",dashed=True,dy=.02)
    box(ax,(.17,.155),.34,.22,"PORTS AND INTRODUCED FLORA",
        "1,596 surveyed ports / 142 countries\n687 candidate introduced host plants\n98 insular flora models",fc="#fff6eb",ec="#b57736")
    box(ax,(.615,.155),.33,.22,"ISLAND / OCEAN GEOGRAPHY",
        "51 selected landmass groups\nEarlier conspecific record proximity\nMarine resistance: modelled, not routes",
        fc="#f1f4f8",ec="#718096")
    arrow(ax,(.36,.385),(.49,.545),"Same host taxa at\nsome port localities",dashed=True,rad=.17,dy=.01)
    arrow(ax,(.78,.385),(.84,.545),"Association and\nindex prediction",dashed=True,rad=-.05,dy=.005)
    ax.text(.5,.445,"Ports ↔ host geography ↔ butterfly indexes are observational connections",
            ha="center",va="center",fontsize=10,color="#5b6573")
    ax.plot([.04,.96],[.09,.09],color="#a2a7b1",lw=.8)
    ax.text(.045,.067,"Solid = defined data reconstruction",fontsize=10,va="center",color="#155461")
    ax.text(.39,.067,"Dashed = observational association / proposed access comparison",fontsize=10,va="center",color="#606976")
    ax.text(.045,.028,"No historical host shipment, larval feeding, butterfly movement or breeding inferred.",
            fontsize=9.7,va="center",fontweight="bold",color="#5f4b43")
    save(fig,out,"Fig1_island_port_butterfly_framework_v02",
         "Distinct grains: 239 resource butterflies; 98 botanical island regions; 225 candidates/51 landmass groups. Diagrams do not establish a port-to-plant-to-butterfly colonisation path.")

def fig5(out):
    early=read("early")
    plant_null=read("plant_null")
    a=read("access")["periods"]["2018_2022"]["within_geographic_cluster_two_margin"]
    oof=read("oof")["periods"]["2018_2022"]
    profile=read("butterfly")
    migrant49=read("migrant49")
    migrant51=read("migrant51")
    gbif=read("gbif51")

    assert early["rates"]["all_source_eligible"]["eligible_first_record_events"]==113
    assert early["rates"]["all_source_eligible"]["historically_recorded_compatible_host_events"]==34
    obs=early["host_network_structural_null"]["observed_prior_host_index_event_fraction"]
    naive=early["host_network_structural_null"]["null_fraction_mean"]
    geo=plant_null["more_constrained_structural_nulls"]["family_plus_region_signature_and_prevalence_band"]["null_mean_fraction"]
    assert abs(obs-34/113)<1e-12 and abs(naive-.25249356261283923)<1e-12
    assert abs(geo-.302587475836629)<1e-12
    assert oof["cohort"]["prior_source_and_trait_events"]==83
    assert gbif["cells"]==11475
    tests=a["tests"]
    near=tests["geo_access_mean"]; sea=tests["sea_residual_access_mean"]
    assert abs(near["p_two_sided_monte_carlo"]-.0003332222592469177)<1e-12
    assert abs(sea["p_two_sided_monte_carlo"]-.0023325558147284237)<1e-12
    holdout=oof["heldout"]
    sea_gain=np.array([holdout[s]["sea_over_geodesic"]["brier"] for s in ("island","geographic_cluster","species")])
    wing_gain=np.array([holdout[s]["wing_by_geodesic"]["brier"] for s in ("island","geographic_cluster","species")])
    assert abs(sea_gain[1]+6.479351063178629)<1e-8
    oldp=migrant51["periods"]["2018_2022"]["within_original_geographic_region"]["tests"]["residual_uniform_ocean_access"]["p_two_sided_monte_carlo"]
    newp=migrant49["periods"]["2018_2022"]["within_original_geographic_region"]["tests"]["residual_uniform_ocean_access"]["p_two_sided_monte_carlo"]
    assert abs(oldp-.04731756081306231)<1e-12 and abs(newp-.22859046984338555)<1e-12

    fig,axs=plt.subplots(2,2,figsize=(14.6,9.2))
    fig.suptitle("Figure 5 | Butterfly indexing supports geographic association, not a demonstrated island-to-island mechanism",
                 x=.035,y=.989,ha="left",fontsize=14.4,fontweight="bold",color="#172735")
    plt.subplots_adjust(left=.115,right=.975,top=.89,bottom=.12,wspace=.34,hspace=.51)

    # A: Historical host first index, null not independent causal effect.
    ax=axs[0,0]
    vals=[100*obs,100*naive,100*geo]
    bars=ax.bar(np.arange(3),vals,width=.65)
    labels=["Recorded compatible\nhost before butterfly","Family/degree-only\nplant network null","Plant-geography\npreserving null"]
    ax.set_xticks(range(3),labels,fontsize=8.6)
    ax.set_ylabel("Butterfly first-index cells with earlier\ncompatible host indexing (%)")
    ax.set_ylim(0,39)
    for n,b in enumerate(bars):
        ax.annotate(f"{vals[n]:.1f}%",(b.get_x()+b.get_width()/2,b.get_height()+.6),
                    ha="center",fontsize=10,fontweight="bold" if n==0 else "normal")
    ax.text(.02,.965,"2001–2010 host index → 2011–2017 butterfly index\n34 of 113 records; geography-matched null ≈ observed",
            transform=ax.transAxes,ha="left",va="top",fontsize=8.8,
            bbox=dict(facecolor="white",edgecolor="none",alpha=.88))
    ax.set_title("A  Plant-before-butterfly records",loc="left",fontsize=11.1)

    # B: Individual butterfly species, not just island-to-island similarity.
    ax=axs[0,1]
    selected=["Zizula hylax","Nymphalis antiopa","Euchrysops cnejus",
              "Danaus plexippus","Vanessa cardui"]
    subset={p["species"]:p for p in profile if p["period"]=="2018_2022" and p["species"] in selected}
    if len(subset)!=len(selected):raise ValueError("Original butterfly profile examples lost")
    vals=[float(subset[s]["median_nearest_prior_index_distance_km_at_events"]) for s in selected]
    count=[int(subset[s]["new_indexed_islands"]) for s in selected]
    y=np.arange(len(selected))
    ax.hlines(y,0,vals,lw=1.9,color="#b4c4c8")
    ax.scatter(vals,y,s=[45+16*x for x in count],color="#207f98",zorder=3)
    ax.set_yticks(y,selected,fontsize=9.4)
    for t in ax.get_yticklabels(): t.set_fontstyle("italic")
    ax.invert_yaxis()
    ax.set_xlim(0,max(vals)*1.37)
    ax.set_xlabel("Nearest previously indexed conspecific island\nwithin sample (median km; not flight distance)")
    for yy,dist,n in zip(y,vals,count):
        ax.annotate(f"{dist:,.0f} km | {n} new indexes",(dist,yy),xytext=(7,0),
                    textcoords="offset points",va="center",fontsize=8.4)
    ax.set_title("B  Which butterflies? (2018–2022)",loc="left",fontsize=11.1)
    ax.text(.97,.97,"5 examples / 83 trait-complete first-index cases",
            va="top",ha="right",transform=ax.transAxes,fontsize=8.5,color="#435064")

    # C: Geographic association vs actual constrained null
    ax=axs[1,0]
    for i,(name,item) in enumerate([("Near previous index",near),("Lower ocean resistance\nthan distance predicts",sea)]):
        p=item["observed_minus_null"]
        ref=np.asarray(item["null_95"])-item["null_mean"]
        ax.plot(ref,[i,i],lw=7,solid_capstyle="round",color="#b7bfc8")
        ax.scatter([p],[i],s=100,color="#2c7f89",zorder=4)
        ax.annotate(f"p={item['p_two_sided_monte_carlo']:.4g}",(p,i),xytext=(7,6),
                    textcoords="offset points",fontsize=9)
    ax.axvline(0,color="#64748b",ls=":",lw=1)
    ax.set_yticks([0,1],["Geographic\nproximity","Residual\nmarine access"],fontsize=10)
    ax.invert_yaxis()
    ax.set_xlim(-.09,max(near["observed_minus_null"],sea["observed_minus_null"])+.10)
    ax.set_xlabel("Observed minus fixed-margin null mean\n(standardized access; 95% null reference band)")
    ax.set_title("C  Index association after species/island/region controls",loc="left",fontsize=11)
    ax.text(.02,.04,"51 landmass groups | 83 late index cases\nRetrospective constrained Monte Carlo; not observed immigration",
            transform=ax.transAxes,fontsize=8.1,va="bottom")

    # D: Genuine loss of heldout value and heavy-port influence.
    ax=axs[1,1]
    x=np.arange(3)
    width=.34
    aa=ax.bar(x-width/2,sea_gain,width,label="Add marine resistance beyond distance",color="#e9a163")
    bb=ax.bar(x+width/2,wing_gain,width,label="Add wing proxy × distance",color="#2b829b")
    ax.axhline(0,color="#606875",lw=1)
    ax.set_xticks(x,["Island-held-out","Region-held-out","Species-held-out"],fontsize=8.9)
    ax.set_ylabel("Late-index butterfly model Brier\nMSE reduction (%)")
    ax.set_ylim(min(-8,float(sea_gain.min())-1),max(3.5,float(wing_gain.max())+1))
    for group,v in ((aa,sea_gain),(bb,wing_gain)):
        for b,vv in zip(group,v):
            ax.annotate(f"{vv:+.2f}%",(b.get_x()+b.get_width()/2,vv),
                        xytext=(0,4 if vv>=0 else -11),
                        textcoords="offset points",fontsize=8,ha="center")
    ax.legend(loc="lower left",bbox_to_anchor=(.0,.06),fontsize=7.6,frameon=True)
    ax.text(.98,.97,
        f"Independently author-listed migrant candidate subset:\n"
        f"sea-residual p={oldp:.3f} (51 groups) → {newp:.3f} (49)\n"
        "after excluding two port-heavy target/source groups",
        transform=ax.transAxes,ha="right",va="top",fontsize=8.1,
        bbox=dict(facecolor="white",edgecolor="0.9",alpha=.95,pad=4))
    ax.set_title("D  Marine predictive transfer and geographic sensitivity",loc="left",fontsize=11)
    for ax in axs.flat:ax.spines[["top","right"]].set_visible(False)

    save(fig,out,"Fig5_historical_hosts_butterfly_access_falsification_v02",
         "Different data supports: A 113 first-index cells/17 landmasses; B 83 indexed butterflies×landmass cells; C/D 225 candidate taxa×51 sampled landmasses, 199 matched trait taxa. All first INDEX records, not true colonisation.")
    return {
       "historical_observed_pre2011_host_index":float(obs),
       "historical_family_degree_null_mean":float(naive),
       "historical_plant_geography_null_mean":float(geo),
       "late_within_region_geo_p":near["p_two_sided_monte_carlo"],
       "late_within_region_sea_residual_p":sea["p_two_sided_monte_carlo"],
       "late_geographical_holdout_sea_Brier_reduction_pct":float(sea_gain[1]),
       "migrant_late_ocean_p_51":float(oldp),
       "migrant_late_ocean_p_49":float(newp),
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",type=Path,default=DEFAULT)
    args=p.parse_args()
    args.out_dir.mkdir(parents=True,exist_ok=True)
    fig1(args.out_dir)
    tests=fig5(args.out_dir)
    out={p.name:{"size":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()}
         for p in sorted(args.out_dir.iterdir()) if p.name.startswith(("Fig1_","Fig5_"))}
    receipt={
       "schema":"chocho_island_port_butterfly_Figure1_Figure5_v02",
       "status":"RENDERED_AND_VERIFIED_CONCEPT_NO_CAUSAL_ARROW_PLUS_SOURCE_EXACT_FALSIFICATION",
       "canonical_manuscript":"manuscript/integrated_butterfly_resource_access_20261109/ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V02.md",
       "source_git_blobs":git_blobs(),
       "tests":tests,
       "files":out,
       "not_actual_flight_shipping_host_feeding_or_colonisation":True
    }
    r=args.out_dir/"FIGURES_15_SOURCE_RECEIPT_V02.json"
    r.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"status":receipt["status"],"figures":list(out),"tests":tests},indent=2))

if __name__=="__main__":main()
