"""Rebuild annual first-index risk sets for the fixed 27 destinations and 132 species."""
import argparse
import io
import json
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from patsy import dmatrix
from pyproj import Geod
from scipy.stats import t

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("workspace", type=Path, help="choch workspace containing outputs and work/barriers")
base = parser.parse_args().workspace
w = base / "work/barriers"
out = root / "results/annual_port_host_27"
out.mkdir(parents=True, exist_ok=True)
figdir = root / "figures/integrated_island_port_butterfly_v03"
def archived(name):
    ref = "ff642c8:provenance/reviewer_defenses/results/expanded_port_islands/" + name
    return pd.read_csv(io.BytesIO(subprocess.check_output(["git","show",ref],cwd=root)))
previous = archived("species_islands.csv")
species = previous.set_index("key").species.to_dict()
destinations = sorted(previous.island.unique())
assert (len(destinations),len(species)) == (27,132)
features = json.loads((w/"ports_2019.geojson").read_text(encoding="utf8"))["features"]
registry = {f["properties"]["OBJECTID"]: f for f in features}
ports = pd.read_csv(base/"outputs/port_buffer_landmass_audit.csv").query("status=='SINGLE_LANDMASS' and islands_intersecting_5km_buffer==1").copy()
ports["island"] = ports.single_landmass_id.astype(int)
countries = ports.port_id.map(lambda i:registry[i]["properties"]["COUNTRY"])
ports = ports[~ports.island.isin(ports.loc[countries.isin(["GB","IE"]),"island"])]
geod = Geod(ellps="WGS84")
centers = {}
for i,g in ports.groupby("island"):
    lon,lat=np.radians([registry[k]["geometry"]["coordinates"][:2] for k in g.port_id]).T
    x,y,z=np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]).mean(axis=1)
    centers[i]=np.degrees([np.arctan2(y,x),np.arctan2(z,np.hypot(x,y))])
assert len(centers)==200
km={(i,j):geod.inv(*centers[i],*centers[j])[2]/1000 for i in destinations for j in centers if i!=j}
annual=pd.read_parquet(w/"global_original_species_island_year.parquet")
annual["island"]=annual.landmass_id.astype(int)
annual=annual[annual.island.isin(centers) & annual.records.gt(0) & annual.year.between(2001,2024)]
first=annual[annual.speciesKey.isin(species)].groupby(["speciesKey","island"]).year.min()
observed=annual.groupby(["island","year"]).speciesKey.agg(set).to_dict()
flora=archived("prior_hosts.csv").groupby("island").host.agg(set).to_dict()
menu=set(pd.read_csv(base/"outputs/port_alien_host_broad_species.csv").host_name)
assert len(menu)==687
diets=pd.read_parquet(w/"all225_named_host_edges_local.parquet")
diets=diets[diets.host_name.isin(menu)].groupby("species").host_name.agg(set).to_dict()
rows=[]
excluded=[]
for key,name in sorted(species.items()):
    history=first.loc[key].to_dict()
    for island in destinations:
        first_here=history.get(island,9999)
        host=int(bool(flora.get(island,set()) & diets[name]))
        for year in range(2011,min(2024,first_here)+1):
            sources=[j for j,y in history.items() if y<year and j!=island]
            if not sources:
                excluded.append((year,key,island,int(first_here==year)))
                continue
            near=min(sources,key=lambda j:(km[island,j],j))
            other=bool(observed.get((island,year),set())-{key})
            rows.append((year,key,island,host,int(first_here==year),near,history[near],km[island,near],int(other)))
d=pd.DataFrame(rows,columns=["year","key","island","host","event","nearest","nearest_first","km","other_record"])
assert (d.nearest_first<d.year).all() and not d.duplicated(["year","key","island"]).any()
assert d.groupby(["key","island"]).event.sum().max()==1
overlap=previous.merge(d.query("year==2024"),on=["key","island"],validate="one_to_one")
assert len(overlap)==len(previous) and (overlap.host_x==overlap.host_y).all() and (overlap.event_x==overlap.event_y).all()
assert (overlap.km_y<=overlap.km_x+1e-7).all()
pd.DataFrame(excluded,columns=["year","key","island","event"]).to_csv(out/"no_prior_source.csv",index=False)
d.to_csv(out/"risk_years.csv.gz",index=False,compression={"method":"gzip","mtime":0})
counts=d.groupby(["year","host"]).event.agg(["size","sum"]).rename(columns={"size":"at_risk","sum":"first_indices"})
counts.to_csv(out/"annual_counts.csv")
d["logkm"]=np.log(d.km)
def estimate(frame,adjustment):
    # OLS score for the host coefficient after partialling out nuisance regressors.
    X=np.asarray(dmatrix(adjustment,frame))
    h=frame.host.to_numpy(); y=frame.event.to_numpy()
    coef,_,rank,_=np.linalg.lstsq(X,np.column_stack([h,y]),rcond=None)
    v,u=(np.column_stack([h,y])-X@coef).T
    b=v@u/(v@v)
    score=v*(u-v*b)/(v@v)
    def meat(group):
        sums=pd.Series(score).groupby(np.asarray(group)).sum()
        n=len(sums)
        return n/(n-1)*(sums@sums)
    variance=(meat(frame.island)+meat(frame.key)-meat(frame.groupby(["island","key"]).ngroup()))*(len(h)-1)/(len(h)-rank-1)
    assert variance>=0
    margin=t.ppf(.975,min(frame.island.nunique(),frame.key.nunique())-1)*np.sqrt(variance)
    return dict(n=len(h),events=int(y.sum()),diff_pp=100*b,lower_pp=100*(b-margin),upper_pp=100*(b+margin))
spline="bs(logkm, df=4)"
full=spline+" + C(year) + C(island) + C(key)"
fits=[]
for year,g in d.groupby("year"):
    fits.append(dict(label=str(year),**estimate(g,spline)))
for label,frame,spec in [("Pooled: distance + year",d,spline+" + C(year)"),
                         ("Pooled: + island + species",d,full),
                         ("Other butterfly recorded",d[d.other_record.eq(1)],full),
                         ("Without Sakhalin",d[d.island.ne(28)],full)]:
    fits.append(dict(label=label,**estimate(frame,spec)))
res=pd.DataFrame(fits)
res.to_csv(out/"rate_differences.csv",index=False)
fig,(a,b)=plt.subplots(1,2,figsize=(11,4.6),gridspec_kw={"width_ratios":[1.15,1]})
fig.subplots_adjust(left=.075,right=.98,top=.85,bottom=.17,wspace=.32)
for h,label,color in [(1,"Prior host record","#16725f"),(0,"Unconfirmed","#775b83")]:
    z=counts.xs(h,level="host"); a.plot(z.index,100*z.first_indices/z.at_risk,"o-",label=label,color=color)
a.legend(frameon=False);a.set(xlabel="First-index year",ylabel="New records per 100 at-risk combinations",xticks=[2011,2014,2017,2020,2024],ylim=(0,None))
a.set_title("A   Annual observed rates",loc="left",weight="bold")
ann=res.iloc[:14]
b.axhline(0,color="#8d9598",ls="--",lw=1)
def interval(ax,x,r,color):
    ax.errorbar(x,r.diff_pp,yerr=[r.diff_pp-r.lower_pp,r.upper_pp-r.diff_pp],fmt="o",color=color,capsize=2.3,ms=4)
interval(b,np.arange(14),ann,"#3d657a")
interval(b,[15,17],res.iloc[14:16],"#a45328")
b.axvline(14,color="#ccd1d2",lw=1)
b.set(xticks=[0,4,8,13,15,17],xticklabels=["2011","2015","2019","2024","All¹","All²"],ylabel="Rate difference (percentage points)")
b.set_title("B   Distance-adjusted contrasts",loc="left",weight="bold")
for ax in (a,b):
    ax.spines[["top","right"]].set_visible(False)
fig.text(.09,.035,"¹ Distance and year.   ² Distance, year, destination island and butterfly species.   Bars show approximate 95% intervals.",size=9)
for ext in ["png","pdf"]:
    fig.savefig(figdir/("FigS3_annual_port_host."+ext),dpi=200)
print(counts.to_string());print(res.round(3).to_string(index=False))
print("Totals",d.groupby("host").event.agg(["size","sum"]).to_dict(),"missing source",len(excluded),"pairs",d.groupby(["island","key"]).ngroups)
