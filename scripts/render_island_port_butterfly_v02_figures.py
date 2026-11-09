#!/usr/bin/env python3
"""Render original-source-backed Figures 2–4 for ONE island/port butterfly paper.

All plotted values come from previously source-locked GitHub research receipts
and files, not a simulated butterfly distribution or an unobserved shipment.
The user-facing paper is ONE integrated GEB-targeted observational manuscript.
Figures 1 (conceptual) and 5 (access/falsification) remain separate tasks.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
OUT_DEFAULT=ROOT/"figures/integrated_island_port_butterfly_v02"
SOURCES={
 "resource":"benchmarks/exploratory/butterfly_anthropogenic_resource_expansion_result_v0.1.json",
 "island":"results/island-mainland-resource-q1/island_mainland_summary.json",
 "matched":"results/island-mainland-resource-area-v02/area_adjusted_summary.json",
 "ports":"provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_validation.json",
 "allflora":"results/three-gate-port-botanical-effort/botanical_effort_adjusted_port_flora_audit.json",
 "hostconditional":"results/three-gate-port-flora-host-conversion/island_port_flora_host_conversion_audit.json",
 "topplants":"provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_top10_pre2011_port_plant_overlap.json",
 "global_butterflies":"data/frozen/figure_sources/anthropogenic_species_metrics.csv",
 "port_butterflies":"provenance/reviewer_defenses/results/port_connectivity_hosts_20261030/port_connectivity_host_all225.csv",
 "bridge":"provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_to_port_species_bridge.json",
 "pair":"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_pairs.csv",
 "pair_receipt":"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_audit.json",
 "butterfly_host_pairs":"provenance/reviewer_defenses/results/island_shared_host_match_20261030/island_shared_host_match_audit.json",
 "butterfly_host_null":"provenance/reviewer_defenses/results/island_host_null_20261030/island_host_family_null_audit.json",
 "forecast_time":"provenance/reviewer_defenses/results/host_island_synthesis_20261030/forecast_temporal_results.json",
 "forecast_spatial":"provenance/reviewer_defenses/results/host_island_synthesis_20261030/forecast_spatial_results.json",
}
def load(key):
    p=ROOT/SOURCES[key]
    if p.suffix==".csv":
        with p.open(newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))
    return json.loads(p.read_text(encoding="utf8"))

def source_blobs():
    return {k:hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
            for k,path in SOURCES.items()
            for b in [(ROOT/path).read_bytes()]}

def annotate(ax,name,x=-.1,y=1.10):
    ax.text(x,y,name,transform=ax.transAxes,fontsize=13,fontweight="bold",va="top")

def finish(fig,out,name,caption):
    fig.text(.025,.014,caption,fontsize=8,color="0.25",va="bottom")
    fig.savefig(out/(name+".png"),dpi=250,bbox_inches="tight")
    fig.savefig(out/(name+".svg"),bbox_inches="tight")
    plt.close(fig)

def figure2(out):
    a=load("island")
    both=a["groups"]
    match=load("matched")["stratum_test"]
    resource=load("resource")["full_resource_eligible_panel"]
    assert a["panel_total"]["panel_species"]==239
    assert len(match["strata_with_both_classes"])==15
    assert resource["species_expanded_by_introduced_host_ranges"]==206
    ins=both["island"];main=both["mainland"]
    assert (ins["native_butterfly_region_units"],ins["contemporary_butterfly_region_units"])==(4829,10570)
    assert (main["native_butterfly_region_units"],main["contemporary_butterfly_region_units"])==(21215,29901)
    fig,axs=plt.subplots(1,3,figsize=(14.4,4.8),gridspec_kw={"width_ratios":[1.05,1.18,1.05]})
    fig.suptitle("Figure 2 | Islands amplify the potential larval-host resources of butterflies",
                 x=.02,y=.99,ha="left",fontsize=14,fontweight="bold")
    ax=axs[0]
    vals=[ins["native_butterfly_region_units"],ins["contemporary_butterfly_region_units"],
          main["native_butterfly_region_units"],main["contemporary_butterfly_region_units"]]
    labels=["Insular\nnative","Insular\ncurrent","Mainland\nnative","Mainland\ncurrent"]
    bs=ax.bar(np.arange(4),vals,width=.68)
    for q,r in zip(bs,vals):
        ax.annotate(f"{r:,}",(q.get_x()+q.get_width()/2,q.get_height()),
                    xytext=(0,4),textcoords="offset points",ha="center",fontsize=9)
    ax.set_ylim(0,33000)
    ax.set_xticks(np.arange(4),labels,fontsize=8.4)
    ax.set_ylabel("Butterfly × botanical region\npotential host-resource combinations")
    ax.text(.1,.92,"+118.9% insular",transform=ax.transAxes,fontsize=10,fontweight="bold")
    ax.text(.1,.86,"+40.9% mainland",transform=ax.transAxes,fontsize=10)
    ax.set_title("A  Native vs contemporary resources",fontsize=11,loc="left")

    ax=axs[1]
    strs=match["strata_with_both_classes"]
    vals=np.array([100*v["difference"] for v in strs],float)
    wts=np.array([v["weight"] for v in strs],float)
    ix=np.argsort(vals)
    xp=np.arange(len(strs))
    ax.scatter(xp,vals[ix],s=30+10*wts[ix],alpha=.8)
    ax.axhline(0,lw=1,ls=":",color="0.45")
    ax.axhline(100*match["weighted_island_minus_mainland_share"],lw=1.1,
               ls="--",color="0.35")
    ax.set_xticks(xp,[strs[k]["stratum"] for k in ix],rotation=75,fontsize=7)
    ax.set_ylabel("Island minus mainland\nintroduced-added share (pp)")
    ax.set_ylim(min(-3,float(vals.min())-4),max(31,float(vals.max())+4))
    ax.annotate("Weighted +10.7 pp\n15 supported strata",
                xy=(.98,.94),xycoords="axes fraction",ha="right",va="top",fontsize=9)
    ax.set_title("B  Matched biogeographic / area strata",fontsize=11,loc="left")

    ax=axs[2]
    data=resource["host_breadth_strata"]
    names=[("1_family","1 family"),("2_families","2 families"),("3_to_5_families","3–5 families"),("6plus_families","6+ families")]
    r=[data[k]["species_expanded"]/data[k]["species"]*100 for k,_ in names]
    bar=ax.bar(np.arange(4),r,width=.72)
    ax.set_ylim(0,110)
    ax.set_xticks(np.arange(4),[n for k,n in names],rotation=18,ha="right",fontsize=8.5)
    for t,(k,_) in zip(bar,names):
        ax.annotate(f'{data[k]["species_expanded"]}/{data[k]["species"]}',
                    (t.get_x()+t.get_width()/2,t.get_height()),xytext=(0,4),
                    textcoords="offset points",ha="center",fontsize=8)
    ax.set_ylabel("Butterflies with any potential gain (%)")
    ax.text(.03,.08,"Global 239 species: +54.9% resource units\n206 species expand; fixed known host list",
            transform=ax.transAxes,fontsize=8.7,va="bottom")
    ax.set_title("C  Gains across dietary specialization",fontsize=11,loc="left")
    for ax in axs: ax.spines[['top','right']].set_visible(False)
    fig.subplots_adjust(left=.065,right=.985,top=.83,bottom=.27,wspace=.43)
    finish(fig,out,"Fig2_islands_resource_butterflies_v02",
           "Source: frozen HOSTS/WCVP × WGSRPD3; 126 insular, 237 mainland regions. Potential host-resource opportunities; NOT butterfly occurrence or colonisation.")

def figure3(out):
    ports=load("ports")
    models=load("allflora")["models"]
    conditional=load("hostconditional")["model_results"]
    plants=load("topplants")
    assert (ports["registered_ports"],ports["complete_ports"],ports["country_codes"])==(3630,1596,142)
    assert (ports["ports_with_alien_candidates"],ports["alien_candidate_species"],ports["port_alien_host_pairs"])==(1150,687,25742)
    assert plants["top10_species_with_plant_index_at_study_island_ports"]==8
    fig,axs=plt.subplots(1,3,figsize=(16.2,5.6),gridspec_kw={"width_ratios":[.95,1.1,1.4]})
    fig.suptitle("Figure 3 | Port-rich islands contain introduced floras that include butterfly hosts",
                 x=.015,y=.99,ha="left",fontsize=14,fontweight="bold")
    ax=axs[0]
    cats=["Registered WPI ports","Completed plant queries","Ports with candidate\nintroduced larval hosts"]
    vals=[ports["registered_ports"],ports["complete_ports"],ports["ports_with_alien_candidates"]]
    y=np.arange(3)
    bars=ax.barh(y,vals)
    ax.set_yticks(y,cats,fontsize=8.7)
    ax.invert_yaxis()
    for b,x in zip(bars,vals):ax.annotate(f"{x:,}",(x,b.get_y()+b.get_height()/2),
                                           xytext=(5,0),textcoords="offset points",va="center",fontsize=8.7)
    ax.set_xlim(0,4300);ax.set_xlabel("WPI port buffers (counts)")
    ax.text(.02,.95,"142 country codes\n687 host candidate plant taxa\n25,742 plant × port pairs",
            transform=ax.transAxes,va="top",fontsize=8.8,
            bbox={"facecolor":"white","edgecolor":"0.85","alpha":.9})
    ax.set_title("A  Global plant-index inventory",loc="left",fontsize=11)

    ax=axs[1]
    values=[
      models["all_flora_with_gbif_total_effort"]["regression_beta_per_log2_port_unit"],
      models["host_taxa_with_effort"]["regression_beta_per_log2_port_unit"],
      conditional["M1eff_port_host_taxa_given_alien_flora"]["coefficient"],
    ]
    assert max(abs(values[i]-v) for i,v in enumerate([.24447852550694815,.16390595789001258,-.041526655753407464]))<1e-10
    labels=["All introduced\nflora richness","Introduced larval-\nhost taxa","Host taxa given\nwhole alien flora"]
    yy=np.array([2,1,0])
    ax.axvline(0,lw=1,ls=":",color="0.45")
    for y0,v in zip(yy,values):
        ax.hlines(y0,0,v,lw=2)
        ax.scatter([v],[y0],s=60)
        ax.text(v+(0.016 if v>=0 else -0.016),y0+.16,f"β {v:+.3f}",
                fontsize=8.5,ha="left" if v>=0 else "right")
    ax.set_xlim(-.22,.39);ax.set_ylim(-.55,3)
    ax.set_yticks(yy,labels,fontsize=8.8)
    ax.set_xlabel("Port number association (model β)\nper log₂(1 + indexed WPI ports)")
    ax.text(.02,.98,"n=98 effort-covered insular regions\nSeparate models; no causal mediation",
            transform=ax.transAxes,fontsize=8,va="top")
    ax.set_title("B  Port associations with flora and hosts",loc="left",fontsize=11)

    ax=axs[2]
    pp=plants["source_contributor_plants"]
    if len(pp)!=10:raise ValueError("Expected 10 global source contributor taxa")
    pp=sorted(pp,key=lambda z:z["GEB_fraction_of_added_species_region_opportunity_attributed_to_this_plant"],
              reverse=True)
    ys=np.arange(10)
    vals=[int(z["n_port_buffers_with_pre2011_record"]) for z in pp]
    bx=ax.barh(ys,vals,height=.73)
    ax.set_yticks(ys,[z["plant_name"] for z in pp],fontsize=8.3)
    ax.invert_yaxis()
    ax.set_xlabel("Study port buffers with plant GBIF index\nin 2001–2010 (among 38)")
    ax.set_xlim(0,13.3)
    for b,z in zip(bx,pp):
        x=int(z["n_port_buffers_with_pre2011_record"])
        share=100*z["GEB_fraction_of_added_species_region_opportunity_attributed_to_this_plant"]
        ax.text(x+.2,b.get_y()+b.get_height()/2,
                f"{x} buffers  |  {share:.1f}% global resource credit",
                ha="left",va="center",fontsize=7.6)
    ax.text(.98,.98,"8 of 10 taxa indexed early\n17 studied landmasses",
            transform=ax.transAxes,va="top",ha="right",fontsize=8,
            bbox={"facecolor":"white","edgecolor":"0.8","alpha":.9})
    ax.set_title("C  Globally influential hosts already indexed near ports",loc="left",fontsize=11)
    for ax in axs:ax.spines[['top','right']].set_visible(False)
    fig.subplots_adjust(left=.11,right=.99,top=.81,bottom=.22,wspace=.73)
    finish(fig,out,"Fig3_ports_aliens_known_butterfly_hosts_v02",
           "Source: WPI/GBIF port-flora inventory, n=98 effort-covered island botanical regions, 2001–2010 plants at 38 port buffers. Not plant carriage, larval use or island-specific introduction dates.")

def figure4(out):
    butterfly=load("global_butterflies")
    port=load("port_butterflies")
    joined={z["species"]:z for z in butterfly}
    pairs=[]
    for row in port:
        species=row["species"]
        if species not in joined:raise ValueError("Unmatched butterfly across full-225 source "+species)
        pairs.append([float(joined[species]["introduced_added_units"]),
                      float(row["island_ports_with_matching_alien_host"])])
    assert len(pairs)==225 and len({z["species"] for z in port})==225
    xx=np.asarray(pairs,float)
    rho=float(spearmanr(xx[:,0],xx[:,1]).statistic)
    bridge=load("bridge")["aligned_port_candidate_comparison"]
    if abs(rho-bridge["spearman_global_added_units_vs_count_of_island_ports_with_candidate_hosts"])>1e-9:
        raise ValueError(f"Butterfly species port/resource rho differs from source lock: {rho}, expected={bridge['spearman_global_added_units_vs_count_of_island_ports_with_candidate_hosts']}")
    pair=load("pair")
    receipt=load("pair_receipt")
    if len(pair)!=1275:raise ValueError("Expected exactly 1275 dependent landmass pairs")
    host=np.array([float(z["host_jaccard"]) for z in pair])
    new=np.array([float(z["new_butterfly_jaccard"]) for z in pair])
    old=np.array([float(z["earlier_butterfly_jaccard"]) for z in pair])
    rhoh=float(spearmanr(host,new).statistic)
    if abs(rhoh-receipt["spearman_host_new"])>1e-8:
        raise ValueError("Observed Jaccard alignment changed")
    matched=load("butterfly_host_pairs");null=load("butterfly_host_null")
    assert (matched["matched_shared_new_species_pairs"],matched["shared_new_species_pairs"])==(208,1163)
    assert matched["fraction_matched"]>null["null_mean_fraction"]
    def host_improvement(q):
        x=q["results"]["M3_both"]["species_equal_mse"]
        y=q["results"]["M5_host_composition"]["species_equal_mse"]
        return (x-y)/x*100
    a=host_improvement(load("forecast_time"))
    b=host_improvement(load("forecast_spatial"))
    if abs(a-4.12017560793182)>1e-7 or abs(b-2.9952059516988094)>1e-7:
        raise ValueError("Forecast source scores drifted")
    fig,axs=plt.subplots(2,2,figsize=(13.8,9),gridspec_kw={"wspace":.28,"hspace":.34})
    fig.suptitle("Figure 4 | Butterfly records and compatible introduced-host geography are associated",
                 x=.025,y=.99,ha="left",fontsize=14,fontweight="bold")
    ax=axs[0,0]
    ax.scatter(np.log1p(xx[:,1]),np.log1p(xx[:,0]),s=20,alpha=.58)
    ax.set_xlabel("ln(1 + island ports with matching\ncandidate introduced hosts)")
    ax.set_ylabel("ln(1 + global added butterfly × region\npotential host opportunities)")
    ax.text(.03,.97,f"225 butterflies · Spearman ρ={rho:+.3f}\nAdjusted partial rank +0.550",
            transform=ax.transAxes,ha="left",va="top",fontsize=9,
            bbox={"facecolor":"white","edgecolor":"0.85","alpha":.9})
    ax.set_title("A  Butterfly species resource gain vs indexed host ports",loc="left",fontsize=11)

    ax=axs[0,1]
    hb=ax.hexbin(host,new,gridsize=27,mincnt=1,bins="log",extent=(0,1,0,1))
    cb=fig.colorbar(hb,ax=ax,fraction=.045,pad=.02)
    cb.set_label("Number of dependent island pairs (log scale)",fontsize=8)
    ax.set_xlabel("Port candidate introduced-host Jaccard")
    ax.set_ylabel("Later-indexed butterfly Jaccard")
    ax.text(.035,.965,"51 landmasses, 1,275 dependent pairs\nLater ρ=+0.300; adjusted +0.215\nEarlier butterflies ρ=+0.418",
            transform=ax.transAxes,va="top",fontsize=8.5,
            bbox={"facecolor":"white","edgecolor":"0.85","alpha":.9})
    ax.set_title("B  Butterfly communities share host geography",loc="left",fontsize=11)

    ax=axs[1,0]
    xs=["Recorded compatible\nhost in both island ports","Family/degree-only\nnetwork null mean"]
    v=[100*matched["fraction_matched"],100*null["null_mean_fraction"]]
    bar=ax.bar(np.arange(2),v,width=.68)
    for b0,z,an in zip(bar,v,["208 / 1,163","199 randomized networks"]):
        ax.text(b0.get_x()+b0.get_width()/2,b0.get_height()+.55,
                f"{z:.1f}%\n{an}",ha="center",fontsize=9)
    ax.set_xticks([0,1],xs)
    ax.set_ylim(0,26)
    ax.set_ylabel("Compatible candidate host co-index (%)")
    ax.text(.04,.92,"82 shared later-indexed butterflies\nPlant records pooled 2001–2024\nNot actual feeding or an area-matched null",
            transform=ax.transAxes,va="top",fontsize=8.6)
    ax.set_title("C  Species-resolved butterfly–host catalogue overlap",loc="left",fontsize=11)

    ax=axs[1,1]
    bar=ax.bar(["Later record period","Held-out spatial block"],[a,b],width=.64)
    for bp,value in zip(bar,[a,b]):
        ax.text(bp.get_x()+bp.get_width()/2,value+.15,f"{value:+.2f}% MSE",ha="center",fontsize=10)
    ax.set_ylim(0,5.6)
    ax.set_ylabel("Species-equal MSE reduction with host geography (%)")
    ax.text(.05,.92,"155 known butterfly species · 208,482\nrecording evaluation rows\nHost map contemporary, not prior arrival",
            transform=ax.transAxes,va="top",fontsize=8.5,
            bbox={"facecolor":"white","edgecolor":"0.85","alpha":.9})
    ax.set_title("D  Added host geography predicts recorded butterflies",loc="left",fontsize=11)
    for ax in axs.flat:ax.spines[['top','right']].set_visible(False)
    fig.subplots_adjust(left=.095,right=.97,top=.90,bottom=.095,wspace=.40,hspace=.50)
    finish(fig,out,"Fig4_butterflies_and_port_host_geography_v02",
           "Independent original populations: A 225 butterflies; B 51 landmasses/1275 dependent pairs; C 1163 candidate butterfly×pair records; D 155 butterfly model evaluations. No feeding, causal plant shipment or butterfly colonisation inferred.")
    return {"species_rank_rho":rho,"host_pair_rho":rhoh,"forecast_time_pct":a,"forecast_spatial_pct":b}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out-dir",type=Path,default=OUT_DEFAULT)
    opt=ap.parse_args()
    opt.out_dir.mkdir(parents=True,exist_ok=True)
    figure2(opt.out_dir)
    figure3(opt.out_dir)
    score=figure4(opt.out_dir)
    outputs={}
    for p in opt.out_dir.iterdir():
        if p.suffix in (".svg",".png"):outputs[p.name]={"bytes":p.stat().st_size,
            "sha256":hashlib.sha256(p.read_bytes()).hexdigest()}
    receipt={
     "schema":"chocho_integrated_island_port_butterfly_figures_234_v02",
     "status":"FIGURES_2_3_4_RENDERED_FROM_SOURCE_LOCKED_NUMERICAL_RECEIPTS",
     "main_manuscript":"manuscript/integrated_butterfly_resource_access_20261109/ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V02.md",
     "source_git_blobs":source_blobs(),
     "visual_files":outputs,
     "biological_interpretation":"Descriptive potential resource, co-indexing and indexed recording prediction; no true larval feeding, transport, actual marine flight or colonisation",
     "source_crosschecks":score,
     "separate_sample_units":["239 global butterfly resource envelopes and 126/237 WGSRPD3 island/mainland regions",
                              "1596 WPI ports from 3630 registered; 98 island botanical-model regions",
                              "225 butterfly port-resource joins, 51 landmass pair comparisons, 155 butterfly forecast"],
     "remaining_figures":["Figure 1 concept and source-grain diagram","Figure 5 historical host and butterfly-marine-access falsification"]
    }
    (opt.out_dir/"FIGURES_234_SOURCE_RECEIPT_V02.json").write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf8")
    print(json.dumps({"status":receipt["status"],"checks":score,"images":list(outputs)},indent=2))

if __name__=="__main__":
    main()
