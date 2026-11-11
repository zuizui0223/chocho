"""Route-independent island risk sets and mixed versus native-only host comparisons."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from acquire_route_free_port_flora import OUT,ROOT,GEOD,load_frame
from analyze_native_plus_introduced_27 import classify,estimate,summarize,archived

FIG=ROOT/"figures/integrated_island_port_butterfly_v03/FigS9_route_free_ports"


def construct(workspace):
    ports,annual=load_frame(workspace)
    coverage=pd.read_csv(OUT/"port_coverage.csv")
    assert set(coverage.port_id)==set(ports.loc[~ports.island.isin(ports.loc[ports.country.isin(["DK","DE"]),"island"]),"port_id"])
    dest=np.sort(coverage.loc[coverage.records.gt(0),"island"].unique())
    species=pd.read_csv(workspace/"outputs/global_island_all225_coverage.csv")[["speciesKey","species"]].rename(columns={"speciesKey":"key"})
    assert len(species)==species.key.nunique()==225
    cells=species.merge(pd.DataFrame({"island":dest}),how="cross")
    taxa=pd.read_csv(OUT/"port_taxa.csv.gz").groupby("port_id").species_key.agg(set).to_dict()
    classified=classify(workspace,cells,coverage,taxa,OUT)
    centers={}
    for island,g in ports.groupby("island"):
        lon,lat=np.radians(g[["lon","lat"]].to_numpy()).T
        x,y,z=np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]).mean(axis=1)
        centers[island]=np.degrees([np.arctan2(y,x),np.arctan2(z,np.hypot(x,y))])
    source=np.array(sorted(centers))
    dist=np.array([[GEOD.inv(*centers[i],*centers[j])[2]/1000 if i!=j else np.inf for j in source] for i in dest])
    annual=annual[annual.island.isin(source)]
    first=annual.groupby(["speciesKey","island"]).year.min()
    observed=annual.groupby(["island","year"]).speciesKey.agg(set).to_dict()
    rows=[];no_source=no_source_events=0
    for k in species.key:
        history=(first.loc[k] if k in first.index.get_level_values(0) else pd.Series(dtype=float))
        years=history.reindex(source).fillna(9999).to_numpy()
        here=history.reindex(dest).fillna(9999).to_numpy()
        for y in range(2011,2025):
            earlier=years<y
            if not earlier.any():
                no_source+=int((here>=y).sum());no_source_events+=int((here==y).sum())
                continue
            ix=np.flatnonzero(earlier)
            near=ix[np.argmin(dist[:,ix],axis=1)]
            for j in np.flatnonzero(here>=y):
                i=dest[j];s=near[j]
                assert i!=source[s] and years[s]<y
                rows.append((i,k,y,int(here[j]==y),source[s],dist[j,s],int(bool(observed.get((i,y),set())-{k}))))
    risk=pd.DataFrame(rows,columns=["island","key","year","event","nearest","km","other_record"])
    assert not risk.duplicated(["island","key","year"]).any() and risk.groupby(["island","key"]).event.sum().max()==1
    risk.to_csv(OUT/"risk_years.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    dates=annual.groupby("island").year.agg(["min","max"])
    registry=pd.read_csv(OUT/"island_registry.csv")
    old=archived("species_islands.csv")
    registry["both_periods"]=registry.island.isin(dates[(dates["min"]<=2010)&(dates["max"]>=2011)].index)
    registry["old_destination"]=registry.island.isin(old.island)
    registry.to_csv(OUT/"island_registry.csv",index=False)
    pd.Series(sorted(old.key.unique()),name="key").to_csv(OUT/"previous_species.csv",index=False)
    (OUT/"screen.json").write_text(json.dumps(dict(source_islands=len(source),outside_islands=coverage.island.nunique(),ports=len(coverage),
        positive_ports=int(coverage.records.gt(0).sum()),flora_and_butterfly_islands=len(dest),candidates=len(species),no_source_years=no_source,
        no_source_events=no_source_events),indent=2)+"\n")


def analyze():
    risk=pd.read_csv(OUT/"risk_years.csv.gz")
    cells=pd.read_csv(OUT/"species_island_hosts.csv")
    registry=pd.read_csv(OUT/"island_registry.csv")
    d=risk.merge(cells,on=["key","island"],validate="many_to_one").merge(registry,on="island",validate="many_to_one")
    assert len(d)==len(risk)
    d["both"]=d.introduced.gt(0).astype(int)
    counts=d.groupby([d.native.gt(0).rename("native_recorded"),"both"]).event.agg(["size","sum"])
    counts.to_csv(OUT/"all_risk_counts.csv")
    comp=d[d.native.gt(0)].copy()
    shared_i=comp.groupby("island").both.nunique();shared_s=comp.groupby("key").both.nunique()
    overlap=comp[comp.island.isin(shared_i[shared_i.eq(2)].index)&comp.key.isin(shared_s[shared_s.eq(2)].index)]
    # Retain both statuses within each island and each species of this subset.
    while len(overlap):
        keep=overlap.groupby("island").both.transform("nunique").eq(2) & overlap.groupby("key").both.transform("nunique").eq(2)
        if keep.all():break
        overlap=overlap[keep]
    groups=pd.DataFrame([dict(both=int(k),**summarize(g)) for k,g in comp.groupby("both")])
    groups.to_csv(OUT/"group_rates.csv",index=False)
    geo=comp.groupby("country").agg(islands=("island","nunique"),species=("key","nunique"),risk_years=("event","size"),events=("event","sum"))
    geo.to_csv(OUT/"regional_support.csv")
    comp.groupby(["island","country","both"]).event.agg(["size","sum"]).to_csv(OUT/"island_support.csv")
    old=pd.read_csv(OUT/"previous_species.csv")
    base="bs(np.log(km), df=4) + C(year)";full=base+" + C(island) + C(key)"
    fits=[]
    for label,g,model in [("Distance + year",comp,base),("Island + species effects",comp,full),
                          ("Shared islands and species",overlap,full),("Earlier 132 candidates",comp[comp.key.isin(old.key)],full),
                          ("Additional islands only",comp[~comp.old_destination],full),
                          ("Without Britain and Ireland",comp[~comp.country.str.contains(r"(?:^|\|)(?:GB|IE)(?:\||$)")],full),
                          ("Both observation periods",comp[comp.both_periods],full),
                          ("Other butterfly recorded",comp[comp.other_record.eq(1)],full)]:
        fits.append(dict(model=label,**estimate(g,model)))
    fits=pd.DataFrame(fits);fits.to_csv(OUT/"rate_differences.csv",index=False)
    support=dict(**summarize(comp),shared_islands=int(shared_i.eq(2).sum()),shared_species=int(shared_s.eq(2).sum()),
                 additional_comparison_islands=comp.loc[~comp.old_destination,"island"].nunique(),registry_groups=len(geo))
    (OUT/"support.json").write_text(json.dumps(support,indent=2)+"\n")
    print(support);print(groups.to_string(index=False));print(fits.round(3).to_string(index=False));print(geo.sort_values("events",ascending=False).head(8))
    fig,(a,b)=plt.subplots(1,2,figsize=(10.8,5.2),gridspec_kw={"width_ratios":[.82,1.40]})
    fig.subplots_adjust(left=.07,right=.98,top=.78,bottom=.27,wspace=.78)
    a.bar([0,1],groups.rate_pct,color=["#67858b","#007866"],width=.56)
    for i,r in enumerate(groups.itertuples()):
        a.text(i,r.rate_pct,f"{r.rate_pct:.2f}%\n{r.events:,} / {r.rows:,}",ha="center",va="bottom",fontsize=9)
    a.set(xticks=[0,1],xticklabels=["Native only","Native +\nintroduced"],ylim=(0,groups.rate_pct.max()*1.5),ylabel="First records per 100 at-risk years")
    a.set_title("A   Recorded annual rates",loc="left",weight="bold",fontsize=12)
    indices=[0,1,2,3,4,5,6,7]
    b.axvline(0,ls="--",color="#94a09f",lw=1)
    for i in indices:
        r=fits.iloc[i];b.errorbar(r.diff_pp,i,xerr=[[r.diff_pp-r.lo],[r.hi-r.diff_pp]],color="#007866" if i==1 else "#597783",fmt="o",capsize=3)
    b.set(yticks=indices,yticklabels=fits.model,ylim=(7.65,-.65),xlabel="Both minus native only (percentage points)")
    b.tick_params(axis="y",labelsize=8)
    b.set_title("B   Adjusted differences and 95% intervals",loc="left",weight="bold",fontsize=11)
    for ax in [a,b]:ax.spines[["top","right"]].set_visible(False)
    fig.text(.035,.955,"Island butterfly first records outside Denmark and Germany",size=16,weight="bold")
    fig.text(.035,.88,f"No route requirement  ·  {support['islands']} islands  ·  {support['species']} butterfly species  ·  {support['events']} first records",size=11)
    fig.text(.035,.117,f"Both groups occurred within {support['shared_islands']} islands and, separately, within {support['shared_species']} species. Plant inventories: 2001–2010; butterfly records: 2011–2024.",size=9)
    fig.text(.035,.06,"Distance is updated annually to the nearest earlier conspecific port island. All rows below the first also adjust for island and species.",size=9)
    for ext in ("png","pdf"):fig.savefig(FIG.with_suffix("."+ext),dpi=200)
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--workspace",type=Path)
    args=parser.parse_args()
    with threadpool_limits(limits=2):
        if args.workspace:construct(args.workspace)
        analyze()
