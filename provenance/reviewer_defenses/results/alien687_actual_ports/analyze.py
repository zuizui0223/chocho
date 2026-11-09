"""Run in the choch workspace with the original botanical and butterfly inventories."""
from pathlib import Path
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out=Path(__file__).resolve().parent
o=Path("outputs")
w=Path("work/barriers")
base=pd.read_csv(out.parent/"top6_actual_ports_20261109/port_pairs.csv")
base=base[["port_a","port_b","island_a","island_b","sampled_voyages","observed","distance_km"]]
assert len(base)==1009 and base.observed.sum()==29 and base.sampled_voyages.sum()==180
candidate=set(pd.read_csv(o/"port_alien_host_broad_species.csv").host_name)
assert len(candidate)==687
plants=pd.read_parquet(w/"port_named_host_regional_status.parquet")
plants=plants[plants.regional_catalogue_status.eq("contemporary_not_native_in_catalogue")]
assert set(plants.host_name)==candidate
hosts=plants.groupby("port_id").host_name.agg(set).to_dict()
links=pd.read_parquet(w/"all225_named_host_edges_local.parquet")
links=links[links.host_name.isin(candidate)]
users=set(links.species)
assert len(users)==192
diet=links.groupby("host_name").species.agg(set).to_dict()
names=pd.read_csv(o/"forecast_candidate_disposition.csv")
history=pd.read_csv(o/"global_original_species_island_history.csv").merge(names[["species","taxon_key"]],left_on="speciesKey",right_on="taxon_key",validate="many_to_one")
history=history[history.species.isin(users)]
years=history.set_index(["landmass_id","species"]).first_record_year.to_dict()
def new_set(lo):
    return history[history.first_record_year.between(lo,2022)].groupby("landmass_id").species.agg(set).to_dict()
new=new_set(2011)
post=new_set(2021)
portnames={f["properties"]["OBJECTID"]:f["properties"]["PORT_NAME"] for f in json.loads((w/"ports_2019.geojson").read_text(encoding="utf-8"))["features"]}
def jaccard(a,b):
    return len(a&b)/len(a|b) if a|b else np.nan
def metrics(a,b,h,g):
    n,m=new.get(a,set()),new.get(b,set())
    return dict(shared=len(h&g),host_j=jaccard(h,g),new_j=jaccard(n,m),
                post_j=jaccard(post.get(a,set()),post.get(b,set())),
                gb_pair={a,b}=={15,26},touches_gb=bool({a,b}&{15,26}))
portrows=[]; matches=[]; islandhosts={}
for r in base.itertuples(index=False):
    i,j=r.port_a,r.port_b
    a,b=r.island_a,r.island_b
    h,g=hosts.get(i,set()),hosts.get(j,set())
    islandhosts.setdefault(a,set()).update(h)
    islandhosts.setdefault(b,set()).update(g)
    portrows.append(dict(**r._asdict(),**metrics(a,b,h,g),name_a=portnames[i],name_b=portnames[j]))
    if r.observed:
        for plant in sorted(h&g):
            spp=diet.get(plant,set()) & new.get(a,set()) & new.get(b,set())
            for sp in sorted(spp) or [""]:
                matches.append(dict(port_a=i,port_b=j,name_a=portnames[i],name_b=portnames[j],island_a=a,island_b=b,host=plant,butterfly=sp,year_a=years.get((a,sp)),year_b=years.get((b,sp))))
p=pd.DataFrame(portrows)
p["islands"]=[tuple(sorted(x)) for x in zip(p.island_a,p.island_b)]
islands=[]
for (a,b),group in p.groupby("islands"):
    islands.append(dict(island_a=a,island_b=b,observed=group.observed.max(),sampled_voyages=group.sampled_voyages.sum(),**metrics(a,b,islandhosts[a],islandhosts[b])))
s=pd.DataFrame(islands)
t=pd.DataFrame(matches)
p.to_csv(out/"port_pairs.csv",index=False)
s.to_csv(out/"island_pairs.csv",index=False)
t.to_csv(out/"shared_hosts_and_users.csv",index=False)
stats=[]
for grain,table in [("port",p),("island",s)]:
    for scope,subset in [("all",table),("omit_GB_IE_link",table[~table.gb_pair]),("omit_both_islands",table[~table.touches_gb])]:
        for edges,data in [("all",subset),("direct",subset[subset.observed.eq(1)])]:
            for x,y in [("observed","shared"),("sampled_voyages","host_j"),("host_j","new_j"),("host_j","post_j")]:
                if (x=="observed" and edges=="direct") or (x=="sampled_voyages" and edges=="all"):continue
                q=data[[x,y]].dropna()
                rho=q[x].rank().corr(q[y].rank()) if len(q)>=3 and q.nunique().min()>1 else np.nan
                stats.append(dict(grain=grain,scope=scope,edges=edges,x=x,y=y,n=len(q),rho=rho))
stats=pd.DataFrame(stats)
stats.to_csv(out/"correlations.csv",index=False)
others=p[p.observed.eq(1)&~p.gb_pair]
details=t.merge(others[["port_a","port_b","islands"]],on=["port_a","port_b"])
rank=details.drop_duplicates(["port_a","port_b","host"]).groupby("host").agg(port_pairs=("port_a","size"),island_pairs=("islands","nunique"))
rank["both_new_users"]=details[details.butterfly.ne("")].groupby("host").butterfly.nunique()
rank=rank.fillna(0).astype(int).sort_values(["port_pairs","island_pairs"],ascending=False)
rank.to_csv(out/"shared_hosts_other_routes.csv")
print(rank.head(12));print(others.groupby(["name_a","name_b"]).shared.max().to_string())
print(details[details.butterfly.ne("")][["host","butterfly","name_a","name_b","year_a","year_b"]].drop_duplicates().to_string(index=False))
for group in [p,p[~p.gb_pair],p[~p.touches_gb]]:
    print(group.groupby("observed").shared.agg(["size","mean",lambda x:(x>0).sum()]))
print(stats.query("scope=='omit_GB_IE_link'").round(4).to_string(index=False))
plt.rcParams["font.family"]="Yu Gothic"
fig,axes=plt.subplots(1,3,figsize=(13,4.2),layout="constrained")
ax=axes[0]; display=rank.head(8).iloc[::-1]
ax.barh(display.index,display.port_pairs,color="#347a72")
ax.set(xlabel="共通の記録がある直航港ペア数",title="共通食草（上位8種まで）")
for ax,df,x,y,xlab,ylab,title in [
    (axes[1],others,"sampled_voyages","host_j","2020年の記録航海数","港の食草共有率","直航港ペア"),
    (axes[2],s[s.observed.eq(1)&~s.gb_pair],"host_j","new_j","島の食草共有率","利用チョウの新規記録共有率","直航のある島ペア")]:
    counts=df.groupby([x,y]).size().reset_index(name="n")
    ax.scatter(counts[x],counts[y],s=18+18*counts.n,color="#286583",alpha=.70)
    ax.set(xlabel=xlab,ylabel=ylab,title=title)
fig.suptitle("グレートブリテン島―アイルランド島間を除いた比較")
fig.savefig(out/"comparisons.png",dpi=165)
