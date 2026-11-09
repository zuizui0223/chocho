"""Species x island risks, observed directed arrivals 2020-21, follow-up 2022-24."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out=Path(__file__).resolve().parent
w=Path('work/barriers');o=Path('outputs')
edges=pd.read_csv(out/'verified_links.csv')
assert set(edges.year)=={2020,2021} and not (set(edges.a)|set(edges.b))&{15,26,456}
ports=pd.read_csv(o/'island_port_early_inventory_status.csv')
assert ports.complete.all() and len(ports)==194
ids=set(ports.landmass_id)-{15,26,456}
eff=ports.groupby('landmass_id').early_plant_records.sum()
plants=pd.read_parquet(w/'port_named_hosts_2001_2010.parquet')
assert plants.year_end.eq(2010).all()
plants['landmass_id']=plants.landmass_id.astype(int)
hosts=set(pd.read_csv(o/'port_alien_host_broad_species.csv').host_name)
menu=pd.read_parquet(w/'all225_named_host_edges_local.parquet')
menu=menu[menu.host_name.isin(hosts)].groupby('species').host_name.agg(set).to_dict()
assert len(hosts)==687 and len(menu)==192
alien=plants[plants.regional_catalogue_status.eq('contemporary_not_native_in_catalogue') & plants.host_name.isin(hosts)].groupby('landmass_id').host_name.agg(set).to_dict()
sp=pd.read_csv(o/'forecast_candidate_disposition.csv').set_index('taxon_key').species
assert len(sp)==225
annual=pd.read_parquet(w/'global_original_species_island_year.parquet')
annual['landmass_id']=annual.landmass_id.astype(int)
annual=annual[annual.landmass_id.isin(ids)&annual.speciesKey.isin(sp.index)]
first=annual.groupby(['speciesKey','landmass_id']).year.min().to_dict()
eligible_islands=set(annual[annual.year.le(2019)].landmass_id)&set(annual[annual.year.between(2022,2024)].landmass_id)&set(eff[eff.gt(0)].index)
d=pd.read_csv(o/'island_host_network_pairs.csv')
km={(r.island_a,r.island_b):r.port_distance_km for r in d.itertuples()}
km.update({(b,a):v for (a,b),v in list(km.items())})
panel=[]
network_nodes=set(edges.a)|set(edges.b)
for dataset,e in edges.groupby('source'):
    available=set(e.a)
    for island,g in e.groupby('b'):
        if island not in eligible_islands:continue
        ny=g.groupby('a').year.nunique()
        for key,name in sp.items():
            if name not in menu:continue
            orig={j for j in ids if first.get((key,j),9999)<=2019}
            yr=first.get((key,island),9999)
            if yr<=2021 or not orig&network_nodes:continue
            nearest=min(km[island,j] for j in orig)
            host=sorted(menu[name]&alien.get(island,set()))
            linked=orig&set(ny.index)
            panel.append(dict(dataset=dataset,island=island,key=key,species=name,km=nearest,
                              host=int(bool(host)),hosts='|'.join(host),route=int(bool(linked)),
                              both_years=int(bool(orig&set(ny[ny.eq(2)].index))),same_source_domain=int(bool(orig&available)),sources='|'.join(map(str,sorted(linked))),
                              event=int(2022<=yr<=2024),first_year=yr if yr<9999 else np.nan))
p=pd.DataFrame(panel)
assert not p.duplicated(['island','key']).any()
p.to_csv(out/'species_islands.csv',index=False)
rates=[]
for scope,q in [('far',p[p.km.ge(500)]),('far_source_matched',p[p.km.ge(500)&p.same_source_domain.eq(1)]),('all_distances',p)]+[(s,p[p.km.ge(500)&p.dataset.eq(s)]) for s in sorted(p.dataset.unique())]:
    for route in ['route','both_years']:
        for r in [0,1]:
            for h in [0,1]:
                g=q[q[route].eq(r)&q.host.eq(h)]
                rates.append(dict(scope=scope,definition=route,route=r,host=h,n=len(g),events=int(g.event.sum()),rate=g.event.mean(),islands=g.island.nunique(),species=g.key.nunique()))
t=pd.DataFrame(rates)
t.to_csv(out/'rates.csv',index=False)
far=p[p.km.ge(500)]
print('far',len(far),far.island.nunique(),far.key.nunique(),far.groupby(['dataset','island']).event.agg(['sum','size']).to_string())
print(t[t.scope.eq('far')].to_string(index=False))
print(far.query('host==1 and route==1 and event==1').to_string(index=False))
plt.rcParams['font.family']='Yu Gothic'
fig,axes=plt.subplots(1,2,figsize=(9,4.9),layout='constrained',sharey=True)
for ax,(col,title) in zip(axes,[('route','いずれかの年で直航を確認'),('both_years','同じ島間の直航を両年で確認')]):
    q=t[(t.scope=='far')&(t.definition==col)]
    for h,dx,label,color in [(0,-.19,'先行食草記録なし','#9ba3ad'),(1,.19,'先行食草記録あり','#207868')]:
        a=q[q.host.eq(h)].sort_values('route')
        bars=ax.bar(np.arange(2)+dx,a.rate.fillna(0)*100,.36,color=color,label=label)
        ax.bar_label(bars,labels=[f'{v.rate:.1%}\n{v.events}/{v.n}' if v.n else '対象なし' for v in a.itertuples()],padding=3)
    ax.set_xticks([0,1],['未確認','確認']);ax.set_title(title)
    ax.set_xlabel('以前の同種記録島からの入港')
    ax.spines[['right','top']].set_visible(False)
axes[0].set_ylabel('2022–2024年の新規記録率（%）')
axes[0].set_ylim(0,max(4,150*t[t.scope.eq('far')].rate.max()))
fig.suptitle('外来食草の先行記録と2020–2021年の直航（500 km以上）')
fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=2,frameon=False)
fig.savefig(out/'multiyear_comparison.png',dpi=175)
