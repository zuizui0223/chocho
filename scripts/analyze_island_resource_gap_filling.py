"""Decompose potential host coverage and filling of native-host gaps in the fixed 239-species panel."""
import argparse
import json
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from analyze_wgsrpd3_island_mainland_resource import host_by_butterfly, host_distribution
from analyze_wgsrpd3_island_area_adjusted import area_by_code

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/island_resource_gap_filling"
FIG=ROOT/"figures/integrated_island_port_butterfly_v03/FigS11_native_host_gaps"
REG=ROOT/"results/island-mainland-resource-q1/region_metrics.csv"


def reconstruct(inputs):
    panel=pd.read_csv(ROOT/"data/frozen/figure_sources/anthropogenic_species_metrics.csv")
    regions=pd.read_csv(REG)
    diets=host_by_butterfly(inputs/"original/insect_host_accepted.csv")
    native=host_distribution(inputs/"original/native_extant_nondoubtful_wgsrpd3.csv")
    current=host_distribution(inputs/"original/extant_nondoubtful_wgsrpd3.csv")
    rows=[]
    for sp in panel.itertuples():
        hosts=diets[sp.species]
        n=Counter(code for h in hosts for code in native[h])
        c=set().union(*(current[h] for h in hosts))
        assert set(n)<=c and (len(n),len(c))==(sp.native_resource_units,sp.contemporary_resource_units)
        for code in regions.area_code_l3:
            rows.append((sp.species,code,n[code],int(code in c)))
    d=pd.DataFrame(rows,columns=["species","area_code_l3","native_hosts","contemporary"])
    d.to_csv(OUT/"resource_cells.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    area=area_by_code(inputs/"level3.geojson")
    regions["area_km2"]=regions.area_code_l3.map(area)
    assert regions.area_km2.gt(0).all()
    regions[["area_code_l3","area_km2"]].to_csv(OUT/"areas.csv",index=False)


def metrics(g):
    n=g.native.sum();a=g.added.sum();total=len(g)
    return dict(cells=total,native=n,empty=total-n,added=a,unfilled=total-n-a,
                coverage=n/total,fill=a/(total-n) if total>n else np.nan,
                gain=a/n if n else np.nan,mean_native_hosts=g.native_hosts.mean())


def analyze():
    d=pd.read_csv(OUT/"resource_cells.csv.gz")
    reg=pd.read_csv(REG).set_index("area_code_l3")
    d["native"]=d.native_hosts.gt(0).astype(int);d["added"]=d.contemporary-d.native
    assert d.added.isin([0,1]).all() and len(d)==239*369 and not d.duplicated(["species","area_code_l3"]).any()
    assert d.groupby("area_code_l3")[["native","added"]].sum().eq(reg[["native","added"]]).all().all()
    reg=reg.join(pd.read_csv(OUT/"areas.csv").set_index("area_code_l3"))
    reg=reg[reg["class"].ne("mixed")].copy()
    cuts=np.quantile(np.log(reg.area_km2),[.25,.5,.75])
    reg["block"]=reg.level1.astype(str)+"|"+np.searchsorted(cuts,np.log(reg.area_km2),side="right").astype(str)
    d=d.merge(reg[["class","block"]],left_on="area_code_l3",right_index=True)
    pooled=pd.DataFrame({k:metrics(g) for k,g in d.groupby("class")}).T
    pooled.to_csv(OUT/"pooled.csv",index_label="geography")
    specific=pd.DataFrame([dict(species=sp,group=k,**metrics(g)) for (sp,k),g in d.groupby(["species","class"])])
    specific.to_csv(OUT/"by_species.csv",index=False)
    fills=specific.pivot(index="species",columns="group",values="fill").dropna()
    delta=fills.island-fills.mainland
    blocks=[];weights=[]
    for block,g in d.groupby("block"):
        parts=[part for _,part in g.groupby("class")]
        if len(parts)!=2:continue
        island,main=map(metrics,parts)
        nr=[len(p)//239 for p in parts];weight=min(nr)
        blocks.append(dict(block=block,island_n=nr[0],mainland_n=nr[1],weight=weight,
                           island_fill=island["fill"],mainland_fill=main["fill"],diff_pp=100*(island["fill"]-main["fill"])))
        weights.append(weight)
    strata=pd.DataFrame(blocks);strata.to_csv(OUT/"area_geography_strata.csv",index=False)
    # Randomly reassign whole regions inside each supported block, preserving class sizes.
    rng=np.random.default_rng(20261112)
    perm=np.zeros(4999)
    for row in strata.itertuples():
        group=reg[reg.block.eq(row.block)]
        order=np.argsort(rng.random((len(perm),len(group))),axis=1)
        adds=group.added.to_numpy()[order]
        gaps=(239-group.native).to_numpy()[order]
        left=adds[:,:row.island_n].sum(axis=1)/gaps[:,:row.island_n].sum(axis=1)
        right=adds[:,row.island_n:].sum(axis=1)/gaps[:,row.island_n:].sum(axis=1)
        perm+=row.weight*100*(left-right)/sum(weights)
    # Symmetric decomposition of A/N = [(1-B)/B] F, with B=N/T, F=A/(T-N).
    odds=(1-pooled.coverage)/pooled.coverage
    baseline=(odds.island-odds.mainland)*pooled.fill.mean()
    filling=(pooled.fill.island-pooled.fill.mainland)*odds.mean()
    assert np.isclose(baseline+filling,pooled.gain.island-pooled.gain.mainland)
    out=dict(equal_species_with_both_denominators=len(delta),island_fill_higher=int((delta>0).sum()),tied=int((delta==0).sum()),
             mainland_fill_higher=int((delta<0).sum()),equal_species_diff_pp=100*delta.mean(),
             matched_blocks=len(strata),matched_island_fill_pct=100*np.average(strata.island_fill,weights=weights),
             matched_mainland_fill_pct=100*np.average(strata.mainland_fill,weights=weights),matched_diff_pp=np.average(strata.diff_pp,weights=weights),
             log_area_cutoffs=cuts.tolist(),baseline_component_pp=100*baseline,fill_component_pp=100*filling)
    out["matched_descriptive_permutation_p"]=(1+int(np.sum(abs(perm)>=abs(out["matched_diff_pp"]))))/(1+len(perm))
    (OUT/"summary.json").write_text(json.dumps(out,indent=2)+"\n")
    print(pooled.to_string());print(out)
    fig,axes=plt.subplots(1,3,figsize=(10.6,4.15))
    fig.subplots_adjust(top=.71,bottom=.23,left=.06,right=.985,wspace=.34)
    names=["Native host coverage","Filling of native-host gaps","Gain over native baseline"]
    for ax,var,name,den in zip(axes,["coverage","fill","gain"],names,["All butterfly–region pairs","Pairs without a native host","Pairs with a native host"]):
        vals=100*pooled.loc[["island","mainland"],var].values
        ax.bar([0,1],vals,color=["#167661","#788d98"],width=.62)
        ax.set(xticks=[0,1],xticklabels=["Island","Mainland"],ylim=(0,135 if var=="gain" else 50),ylabel="%")
        ax.set_title(name,loc="left",weight="bold",fontsize=11)
        for i,v in enumerate(vals):ax.text(i,v+(.7 if var!="gain" else 2),f"{v:.1f}%",ha="center",weight="bold")
        ax.spines[["top","right"]].set_visible(False)
        ax.text(.5,-.24,den,ha="center",transform=ax.transAxes,fontsize=8.6)
    fig.text(.045,.95,"Why is proportional host-resource gain larger on islands?",size=16,weight="bold")
    fig.text(.045,.855,"239 butterflies  ·  126 island and 237 mainland botanical regions  ·  native and contemporary ranges of known larval hosts",size=9.5)
    fig.text(.045,.045,f"Gap filling after standardizing area and broad geography: island {out['matched_island_fill_pct']:.1f}%, mainland {out['matched_mainland_fill_pct']:.1f}% ({out['matched_blocks']} shared strata).",size=9.5)
    for ext in ("png","pdf"):fig.savefig(FIG.with_suffix("."+ext),dpi=210)
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--inputs",type=Path)
    inputs=parser.parse_args().inputs
    OUT.mkdir(parents=True,exist_ok=True)
    if inputs:reconstruct(inputs)
    analyze()
