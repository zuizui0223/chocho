"""Run with the original choch workspace as the working directory."""
from pathlib import Path
from itertools import combinations
from io import BytesIO
import subprocess
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = Path(__file__).resolve().parent
repo = Path("work/remote-sensing-runner")
o = Path("outputs")
source = "bb9f7d08e09416635725359eb4c7c3c48b866a33:provenance/reviewer_defenses/results/island_observational_robustness_20261031/ais_2020_wpi51_island_5km."
def archived(suffix):
    return pd.read_csv(BytesIO(subprocess.check_output(["git","-C",str(repo),"show",source+suffix])))

ports = archived("WPI_port_AIS_source_coverage.csv").query("ais2020_source_endpoint_rows > 0").set_index("wpi_port_id")
edges = archived("observed_port_links.csv")
voyages = {tuple(sorted((r.wpi_port_id_a,r.wpi_port_id_b))):r.ais_2020_sampled_voyage_rows for r in edges.itertuples()}
inventory = pd.read_csv(o/"port_alien_host_broad_all_ports.csv")
assert set(ports.index) <= set(inventory.loc[inventory.query_complete,"port_id"])
plants = pd.read_parquet("work/barriers/port_named_host_regional_status.parquet")
top = pd.read_csv(o/"top_shared_host_island_record_summary.csv").host_name.tolist()
links = pd.read_csv(o/"top_shared_host_link_provenance.csv")
assert len(top)==6 and set(top)==set(links.host_name)
diet = links.groupby("host_name").species.agg(set).to_dict()
users = set(links.species)
local = plants[plants.host_name.isin(top) & plants.regional_catalogue_status.eq("contemporary_not_native_in_catalogue")]
host = local.groupby("port_id").host_name.agg(set).to_dict()
history = pd.read_csv(o/"global_original_species_island_history.csv").merge(pd.read_csv(o/"forecast_candidate_disposition.csv")[["species","taxon_key"]],left_on="speciesKey",right_on="taxon_key",validate="many_to_one")
def assemblage(lo,hi):
    return history[history.species.isin(users) & history.first_record_year.between(lo,hi)].groupby("landmass_id").species.agg(set).to_dict()
new = assemblage(2011,2022)
post = assemblage(2021,2022)
old = assemblage(2001,2010)
def jac(a,b):
    return len(a&b)/len(a|b) if a|b else np.nan
coords = {f["properties"]["OBJECTID"]:np.radians(f["geometry"]["coordinates"]) for f in json.loads(Path("work/barriers/ports_2019.geojson").read_text(encoding="utf-8"))["features"]}
def km(i,j):
    x,y=coords[i];u,v=coords[j]
    return 12742.0176*np.arcsin(np.sqrt(np.clip(np.sin((y-v)/2)**2+np.cos(y)*np.cos(v)*np.sin((x-u)/2)**2,0,1)))
rows=[]
for i,j in combinations(sorted(ports.index),2):
    a,b=(int(ports.at[k,"landmass_id"]) for k in (i,j))
    if a==b:continue
    h,g=(host.get(k,set()) for k in (i,j))
    n,m=(new.get(k,set()) for k in (a,b))
    share=h&g
    match=sorted(set().union(*(diet[t] for t in share)) & n & m)
    v=voyages.get((i,j),0)
    rows.append(dict(port_a=i,port_b=j,island_a=a,island_b=b,distance_km=km(i,j),
                     observed=int(v>0),sampled_voyages=v,shared_count=len(share),shared_hosts="|".join(sorted(share)),
                     host_j=jac(h,g),new_j=jac(n,m),post_j=jac(post.get(a,set()),post.get(b,set())),
                     early_j=jac(old.get(a,set()),old.get(b,set())),matched_species="|".join(match)))
d=pd.DataFrame(rows)
assert d.sampled_voyages.sum()==180
d.to_csv(out/"port_pairs.csv",index=False)
# Pool catalogues across AIS-covered ports for a single observation per island pair.
islandhosts={i:set().union(*(host.get(k,set()) for k in g.index)) for i,g in ports.groupby("landmass_id")}
d["pair"]=[tuple(sorted((a,b))) for a,b in zip(d.island_a,d.island_b)]
ir=[]
for (a,b),q in d.groupby("pair"):
    h,g=islandhosts[a],islandhosts[b]
    ir.append(dict(island_a=a,island_b=b,observed=q.observed.max(),sampled_voyages=q.sampled_voyages.sum(),
                   shared_count=len(h&g),host_j=jac(h,g),new_j=q.new_j.iloc[0],post_j=q.post_j.iloc[0],
                   early_j=q.early_j.iloc[0],distance_km=q.distance_km.min()))
isl=pd.DataFrame(ir)
isl.to_csv(out/"island_pairs.csv",index=False)
def rho(x,y):
    if len(x)<3 or x.nunique()<2 or y.nunique()<2:return np.nan
    return x.rank().corr(y.rank())
stats=[]
for grain,table in [("port",d),("island",isl)]:
    for x,y in [("observed","shared_count"),("observed","host_j"),("sampled_voyages","host_j"),("host_j","new_j"),("host_j","post_j"),("observed","new_j")]:
        q=table.dropna(subset=[x,y])
        if x=="sampled_voyages":q=q[q.observed.eq(1)]
        islands=set(q.island_a)|set(q.island_b)
        leave=[rho(z[x],z[y]) for i in islands for z in [q[(q.island_a!=i)&(q.island_b!=i)]]]
        stats.append(dict(grain=grain,x=x,y=y,n=len(q),documented=int(q.observed.sum()),rho=rho(q[x],q[y]),
                          omitted_island_min=pd.Series(leave).min(),omitted_island_max=pd.Series(leave).max()))
pd.DataFrame(stats).to_csv(out/"correlations.csv",index=False)
plt.rcParams["font.family"]="Yu Gothic"
fig,(ax,cx,bx)=plt.subplots(1,3,figsize=(12,3.9),layout="constrained")
for pos,v in enumerate([0,1]):
    z=d[d.observed.eq(v)].shared_count.value_counts().sort_index()
    ax.scatter(np.full(len(z),pos),z.index,s=18+3*np.sqrt(z)*6,color="#286583",alpha=.65)
    ax.plot([pos-.2,pos+.2],[d[d.observed.eq(v)].shared_count.mean()]*2,color="black")
ax.set(xticks=[0,1],xticklabels=["直航記録なし","直航記録あり"],ylabel="両港で記録された上位食草数",title="2020年AISで照合できた港ペア",xlim=(-.5,1.5))
positive=d[d.observed.eq(1)].dropna(subset=["host_j"]).copy()
gb=positive.pair.map(lambda pair: pair==(15,26))
for mask,col,lab in [(~gb,"#286583","その他の島間"),(gb,"#bd5528","英国・アイルランド間")]:
    z=positive[mask].groupby(["sampled_voyages","host_j"]).size().reset_index(name="n")
    cx.scatter(z.sampled_voyages,z.host_j,s=22+14*np.sqrt(z.n),alpha=.72,label=lab,color=col)
cx.set(xscale="log",xlabel="記録航海数（対数目盛）",ylabel="上位食草の共有率",title="直航記録のある港ペア")
cx.legend(fontsize=8)
print("GB-Ireland only:",len(positive[gb]),rho(positive[gb].sampled_voyages,positive[gb].host_j))
print("GB-Ireland excluded:",len(positive[~gb]),rho(positive[~gb].sampled_voyages,positive[~gb].host_j))
for key,group in [("pooled",positive),("nonGBIE",positive[~gb])]:
    z=group[["sampled_voyages","host_j","distance_km"]].rank()
    if z.host_j.nunique()<2:continue
    design=np.column_stack([np.ones(len(z)), z.distance_km])
    res=z.iloc[:,:2].to_numpy()-design@np.linalg.lstsq(design,z.iloc[:,:2].to_numpy(),rcond=None)[0]
    print(key,"distance adjusted",np.corrcoef(res.T)[0,1])
q=isl.dropna(subset=["host_j","new_j"])
for v,col,lab in [(0,"#a9a9a9","直航記録なし"),(1,"#bd5528","直航記録あり")]:
    z=q[q.observed.eq(v)].groupby(["host_j","new_j"]).size().reset_index(name="n")
    bx.scatter(z.host_j,z.new_j,s=18+12*np.sqrt(z.n),color=col,alpha=.7,label=lab)
bx.set(xlabel="上位食草の共有率",ylabel="利用チョウの新規記録種の共有率",title="島ペアごとの種組成（Jaccard指数）")
bx.legend(fontsize=8)
fig.savefig(out/"correlations.png",dpi=170)
print(pd.DataFrame(stats).round(3).to_string(index=False))
print(d.groupby("observed").agg(n=("observed","size"),mean_shared=("shared_count","mean"),any_shared=("shared_count",lambda s:int((s>0).sum()))))
print("ports/islands/users",len(ports),ports.landmass_id.nunique(),len(users))
print(d[d.observed.eq(1) & d.shared_count.gt(0)][["shared_hosts","matched_species"]].value_counts().to_string())
