"""Run from the original choch data workspace."""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out=Path(__file__).resolve().parent
o=Path("outputs")
w=Path("work/barriers")
ports=pd.read_csv(o/"island_port_early_inventory_status.csv")
assert len(ports)==194 and ports.complete.all()
ids=set(ports.landmass_id)
assert len(ids)==51
ships=pd.read_csv(out.parent/"alien687_actual_ports/port_pairs.csv")
covered=set(ships.island_a)|set(ships.island_b)
edges=ships[ships.observed.eq(1)]
assert len(covered)==21 and len(edges)==29 and edges.sampled_voyages.sum()==180
distances=pd.read_csv(o/"island_host_network_pairs.csv")
dist={(int(r.island_a),int(r.island_b)):r.port_distance_km for r in distances.itertuples()}
dist.update({(j,i):v for (i,j),v in list(dist.items())})
assert len(dist)==51*50
plants=pd.read_parquet(w/"port_named_hosts_2001_2010.parquet")
assert plants.year_start.eq(2001).all() and plants.year_end.eq(2010).all()
hosts=set(pd.read_csv(o/"port_alien_host_broad_species.csv").host_name)
assert len(hosts)==687
plants=plants[plants.regional_catalogue_status.eq("contemporary_not_native_in_catalogue") & plants.host_name.isin(hosts)].copy()
plants["landmass_id"]=plants.landmass_id.astype(int)
catalog=plants.groupby("landmass_id").host_name.agg(set).to_dict()
portcatalog=plants.groupby("port_id").host_name.agg(set).to_dict()
plant_effort=ports.groupby("landmass_id").early_plant_records.sum()
links=pd.read_parquet(w/"all225_named_host_edges_local.parquet")
diet=links[links.host_name.isin(hosts)].groupby("species").host_name.agg(set).to_dict()
species=pd.read_csv(o/"forecast_candidate_disposition.csv").set_index("taxon_key").species
assert len(species)==225 and len(diet)==192
annual=pd.read_parquet(w/"global_original_species_island_year.parquet")
annual["landmass_id"]=pd.to_numeric(annual.landmass_id,errors="raise").astype(int)
annual=annual[annual.landmass_id.isin(ids)&annual.speciesKey.isin(species.index)]
assert len(set(annual.landmass_id))==51
first=annual.groupby(["speciesKey","landmass_id"]).year.min().to_dict()
connected={}
for r in edges.itertuples():
    for a,b,p in [(r.island_a,r.island_b,r.port_a),(r.island_b,r.island_a,r.port_b)]:
        connected.setdefault((a,b),set()).add(p)
rows=[]
for k,name in species.items():
    if name not in diet:continue
    sources={j for j in ids if first.get((k,j),9999)<=2019}
    if not sources&covered:continue
    for i in sorted(covered):
        yr=first.get((k,i),9999)
        if yr<=2020 or plant_effort[i]==0:continue
        nearest=min(dist[i,j] for j in sources)
        targetports=set().union(*(connected.get((i,j),set()) for j in sources))
        prior=diet[name]&catalog.get(i,set())
        exact=prior&set().union(*(portcatalog.get(p,set()) for p in targetports))
        rows.append(dict(key=k,species=name,island=i,nearest_km=nearest,route=int(bool(targetports)),
                         prior_host=int(bool(prior)),prior_names="|".join(sorted(prior)),
                         route_ports="|".join(map(str,sorted(targetports))),exact_host=int(bool(exact)),
                         source_islands="|".join(map(str,sorted(sources))),
                         linked_sources="|".join(map(str,sorted(j for j in sources if (i,j) in connected))),
                         event_year=yr if yr<9999 else np.nan))
panel=pd.DataFrame(rows)
panel.to_csv(out/"eligible_species_islands.csv",index=False)
rates=[];effects=[]
active=set(edges.island_a)|set(edges.island_b)
for end in [2022,2024]:
    observable=set(annual.loc[annual.year.between(2021,end),"landmass_id"])
    prior_observable=set(annual.loc[annual.year.le(2019),"landmass_id"])
    far=panel[panel.nearest_km.ge(500)&panel.island.isin(observable&prior_observable)].copy()
    far["event"]=far.event_year.between(2021,end).astype(int)
    for scope,frame in [("all",far),("no_GB_IE_targets",far[~far.island.isin([15,26])]),("known_link_targets",far[far.island.isin(active)])]:
        for r in [0,1]:
            for h in [0,1]:
                g=frame[frame.route.eq(r)&frame.prior_host.eq(h)]
                rates.append(dict(end=end,scope=scope,route=r,host=h,n=len(g),events=g.event.sum(),
                                  rate=g.event.mean(),islands=g.island.nunique(),species=g.key.nunique()))
        def diff(frame):
            z=frame.groupby(["route","prior_host"]).event.mean()
            if len(z)!=4:return np.nan
            return (z[1,1]-z[1,0])-(z[0,1]-z[0,0])
        deletes=[diff(frame[frame.island.ne(i)]) for i in frame.island.unique()]
        effects.append(dict(end=end,scope=scope,interaction_pp=100*diff(frame),drop_one_min_pp=100*np.nanmin(deletes),drop_one_max_pp=100*np.nanmax(deletes)))
    print(end,"targets",far.island.nunique(),"species",far.key.nunique(),"same endpoint n/events",far.query("exact_host==1").event.agg(["size","sum"]).to_dict())
rates=pd.DataFrame(rates)
rates.to_csv(out/"rates.csv",index=False)
pd.DataFrame(effects).to_csv(out/"contrasts.csv",index=False)
print(rates.query("scope=='all'").to_string(index=False))
print(pd.DataFrame(effects).to_string(index=False))
plt.rcParams["font.family"]="Yu Gothic"
fig,axes=plt.subplots(1,2,figsize=(9.3,4.5),sharey=True,layout="constrained")
for ax,(end,group) in zip(axes,rates.query("scope=='all'").groupby("end")):
    for h,color,dx,label in [(0,"#9aa4ae",-.19,"外来食草の先行記録なし"),(1,"#247b69",.19,"外来食草の先行記録あり")]:
        z=group[group.host.eq(h)].sort_values("route")
        bars=ax.bar(np.arange(2)+dx,z.rate*100,width=.36,color=color,label=label)
        ax.bar_label(bars,labels=[f"{v.rate:.1%}\n{v.events}/{v.n}" for v in z.itertuples()],padding=4)
    ax.set_xticks([0,1],["直航未確認","直航記録あり"])
    ax.set_title(f"2021–{end}年の新規記録")
    ax.spines[["right","top"]].set_visible(False)
axes[0].set_ylim(0,max(1,rates.query("scope=='all'").rate.max()*150))
axes[0].set_ylabel("島での新規記録率（％）")
fig.suptitle("同種の既往記録島から500 km以上離れた島")
fig.legend(*axes[0].get_legend_handles_labels(),loc="outside lower center",ncol=2,frameon=False)
fig.savefig(out/"source_route_host.png",dpi=165)
