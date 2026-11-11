"""Compare annual first records with a prior native host, by prior introduced host status."""
import argparse
import json
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from patsy import dmatrix
from pyproj import Geod
from scipy.stats import t
from shapely.geometry import Polygon
from trace_top_hosts_to_island_records import archived

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/native_plus_introduced_27"
FIG=ROOT/"figures/integrated_island_port_butterfly_v03/FigS8_native_plus_introduced"


def classify(workspace):
    w=workspace/"work/barriers"
    s=workspace/"work/remote_sensing/host_source/extracted/inputs"
    previous=archived("species_islands.csv")
    coverage=archived("plant_coverage.csv")
    coverage=coverage[coverage.island.isin(previous.island.unique())]
    edges=pd.read_parquet(w/"all225_named_host_edges_local.parquet").drop_duplicates()
    lookup=pd.read_parquet(w/"host_flora_name_join_audit.parquet").query("join_status=='SPECIES_FACET_COMPATIBLE'")
    names=lookup.set_index("accepted_gbif_species_key").host_name.to_dict()
    reg=json.loads((w/"ports_2019.geojson").read_text())["features"]
    countries={p["properties"]["OBJECTID"]:p["properties"]["COUNTRY"] for p in reg}
    registry=coverage.assign(country=coverage.port_id.map(countries)).groupby("island").country.agg(lambda x:"|".join(sorted(set(x))))
    registry.to_csv(OUT/"island_registry.csv")
    xy={p["properties"]["OBJECTID"]:p["geometry"]["coordinates"][:2] for p in reg}
    units=gpd.read_file(s/"level3.geojson").to_crs(4326)
    geod=Geod(ellps="WGS84")
    ids=pd.read_csv(s/"original/insect_host_accepted.csv",dtype=str).set_index("accepted_plant_name_id").accepted_name.to_dict()
    for butterfly in json.loads((s/"bce/accepted_host_edges.json").read_text())["per_butterfly"]:
        for h in butterfly["details"]:
            if h.get("match_status")=="ACCEPTED_ID_ABSENT_FROM_FROZEN_HOSTS_FOR_BUTTERFLY":
                ids[str(h["accepted_plant_name_id"])]=h["accepted_name"]
    def catalogue(original,extra):
        d=pd.concat([pd.read_csv(s/original,dtype=str),pd.read_csv(s/extra,dtype=str)])
        return set(zip(d.accepted_plant_name_id.map(ids),d.area_code_l3))
    native=catalogue("original/native_extant_nondoubtful_wgsrpd3.csv","bce/additional_native_wgsrpd3.csv")
    present=catalogue("original/extant_nondoubtful_wgsrpd3.csv","bce/additional_contemporary_wgsrpd3.csv")
    rows=[]
    for r in coverage.itertuples():
        if r.method=="prior_complete":
            d=json.loads((w/f"port_flora_2001_2010/{r.port_id}.json").read_text())
            keys={int(x["name"]) for x in d["species_counts"]}
        else:
            d=json.loads((w/f"expanded_port_flora_2001_2010/{r.port_id}.json").read_text())
            keys=set(d["keys"])
        assert int(r.plant_records)==d["count"]
        polygon=Polygon([geod.fwd(*xy[r.port_id],-angle,5000)[:2] for angle in range(0,360,15)])
        codes=set(units.iloc[units.sindex.query(polygon,predicate="intersects")].LEVEL3_COD.dropna())
        for key in sorted(keys & names.keys()):
            host=names[key]
            status="unresolved"
            if codes and all((host,c) in native for c in codes):
                status="native"
            elif codes and all((host,c) in present and (host,c) not in native for c in codes):
                status="introduced"
            rows.append((r.island,r.port_id,host,status))
    plants=pd.DataFrame(rows,columns=["island","port_id","host_name","status"])
    plants.to_csv(OUT/"port_host_status.csv",index=False)
    species=previous[["key","species"]].drop_duplicates()
    cells=species.merge(previous[["island"]].drop_duplicates(),how="cross")
    links=plants.merge(edges,on="host_name").merge(species,on="species")
    for status in ("native","introduced","unresolved"):
        counts=links[links.status.eq(status)].groupby(["island","key"]).host_name.nunique().rename(status)
        cells=cells.join(counts,on=["island","key"])
        cells[status]=cells[status].fillna(0).astype(int)
    cells.to_csv(OUT/"species_island_hosts.csv",index=False)
    assert len(cells)==27*132
    return cells


def estimate(frame,model):
    X=np.asarray(dmatrix(model,frame))
    Y=frame[["both","event"]].to_numpy()
    coef,_,rank,_=np.linalg.lstsq(X,Y,rcond=None)
    h,y=(Y-X@coef).T
    assert h@h>1e-6
    b=h@y/(h@h)
    score=pd.Series(h*(y-b*h)/(h@h))
    def meat(cluster):
        sums=score.groupby(np.asarray(cluster)).sum()
        return len(sums)/(len(sums)-1)*(sums@sums)
    variance=(meat(frame.island)+meat(frame.key)-meat(frame.groupby(["island","key"]).ngroup()))*(len(frame)-1)/(len(frame)-rank-1)
    assert variance>=0
    margin=t.ppf(.975,min(frame.island.nunique(),frame.key.nunique())-1)*np.sqrt(variance)
    return dict(n=len(frame),events=int(frame.event.sum()),islands=frame.island.nunique(),species=frame.key.nunique(),
                diff_pp=100*b,lo=100*(b-margin),hi=100*(b+margin))


def summarize(d):
    return dict(rows=len(d),pairs=len(d[["island","key"]].drop_duplicates()),islands=d.island.nunique(),species=d.key.nunique(),
                events=int(d.event.sum()),event_islands=d.loc[d.event.eq(1),"island"].nunique(),rate_pct=100*d.event.mean())


def analysis(cells):
    risk=pd.read_csv(ROOT/"results/annual_port_host_27/risk_years.csv.gz")
    risk=risk.merge(cells,on=["island","key"],validate="many_to_one")
    assert len(risk)==42228 and risk.event.sum()==266
    risk["both"]=risk.introduced.gt(0).astype(int)
    risk["group"]=np.select([risk.native.gt(0)&risk.both.eq(1),risk.native.gt(0),risk.both.eq(1)],
                            ["Both recorded","Native only recorded","Introduced only recorded"],default="Neither confirmed")
    counts=pd.DataFrame([dict(group=k,**summarize(g)) for k,g in risk.groupby("group")])
    counts.to_csv(OUT/"risk_groups.csv",index=False)
    comp=risk[risk.native.gt(0)].copy()
    islands=comp.groupby("island").both.nunique()
    species=comp.groupby("key").both.nunique()
    shared=comp[comp.island.isin(islands[islands.eq(2)].index)&comp.key.isin(species[species.eq(2)].index)]
    full="bs(np.log(km), df=4) + C(year) + C(island) + C(key)"
    fits=[]
    for name,frame,model in [("Distance + year",comp,"bs(np.log(km), df=4) + C(year)"),
                             ("+ island + species",comp,full),("Shared islands and species",shared,full),
                             ("Other butterfly recorded",comp[comp.other_record.eq(1)],full)]:
        fits.append(dict(model=name,**estimate(frame,model)))
    fits=pd.DataFrame(fits)
    fits.to_csv(OUT/"rate_differences.csv",index=False)
    audit=dict(shared_islands=int(islands.eq(2).sum()),shared_species=int(species.eq(2).sum()),
               legacy_introduced_flag_differences=int(risk.host.ne(risk.both).sum()),comparison=summarize(comp))
    (OUT/"support.json").write_text(json.dumps(audit,indent=2)+"\n")
    by_island=comp.groupby(["island","group"]).event.agg(["size","sum"])
    by_island.to_csv(OUT/"island_support.csv")
    geographic=by_island.reset_index().merge(pd.read_csv(OUT/"island_registry.csv"),on="island",validate="many_to_one")
    geographic.groupby("country").agg(islands=("island","nunique"),events=("sum","sum"),at_risk=("size","sum")).to_csv(OUT/"regional_support.csv")
    print(counts.to_string(index=False));print(fits.to_string(index=False));print(audit)
    return counts,fits,audit


def draw(counts,fits,audit):
    fig,(a,b)=plt.subplots(1,2,figsize=(10.8,4.9),gridspec_kw={"width_ratios":[.85,1.42]})
    fig.subplots_adjust(left=.07,right=.97,bottom=.31,top=.75,wspace=.66)
    focus=counts.set_index("group").loc[["Native only recorded","Both recorded"]]
    colors=["#789098","#087260"]
    a.bar([0,1],focus.rate_pct,color=colors,width=.54)
    for i,r in enumerate(focus.itertuples()):
        a.text(i,r.rate_pct,f"{r.rate_pct:.2f}%\n{r.events:,} / {r.rows:,}",ha="center",va="bottom",fontsize=10)
    a.set(xticks=[0,1],xticklabels=["Native only","Native +\nintroduced"],ylabel="First records per 100 at-risk years",ylim=(0,focus.rate_pct.max()*1.50))
    a.set_title("A   Observed annual rates",loc="left",weight="bold",fontsize=12)
    b.axvline(0,color="#879493",ls="--",lw=1)
    for i,r in enumerate(fits.itertuples()):
        b.errorbar(r.diff_pp,i,xerr=[[r.diff_pp-r.lo],[r.hi-r.diff_pp]],fmt="o",color="#087260" if i==1 else "#567483",capsize=4)
    b.set(yticks=range(len(fits)),yticklabels=["Distance + year","+ island + species","Shared islands\nand species","Other butterfly\nrecorded"],xlabel="Annual rate difference (percentage points)")
    b.set_ylim(3.6,-.65)
    b.set_title("B   Both minus native only (95% intervals)",loc="left",weight="bold",fontsize=11)
    for ax in (a,b):ax.spines[["top","right"]].set_visible(False)
    fig.text(.035,.95,"Does an additional introduced host predict island butterfly first records?",size=15,weight="bold")
    q=audit["comparison"]
    fig.text(.035,.86,f"27 islands × 132 species screened; native-host comparison: {q['islands']} islands, {q['species']} species, {q['events']} first records.",size=10)
    fig.text(.035,.13,f"Both exposure groups contributed within {audit['shared_islands']} islands and, separately, within {audit['shared_species']} butterfly species.",size=10)
    fig.text(.035,.065,"Hosts recorded at island ports, 2001–2010; first butterfly records, 2011–2024. Unconfirmed hosts can still be present.",size=9)
    for ext in ("png","pdf"):fig.savefig(FIG.with_suffix("."+ext),dpi=200)
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--workspace",type=Path,help="Reconstruct classified plants from the archived port facets")
    args=parser.parse_args()
    OUT.mkdir(exist_ok=True,parents=True)
    cells=classify(args.workspace) if args.workspace else pd.read_csv(OUT/"species_island_hosts.csv")
    draw(*analysis(cells))
