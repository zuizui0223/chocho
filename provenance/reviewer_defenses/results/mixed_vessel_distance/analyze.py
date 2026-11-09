"""Species x island risks, observed directed arrivals 2020-21, follow-up 2022-24."""
from pathlib import Path
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out=Path(__file__).resolve().parent
w=Path('work/barriers');o=Path('outputs')
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
edges=pd.read_csv(out/'links.csv')
nodes=pd.read_csv(out/'endpoint_coverage.csv')
assert set(edges.kind)=={'container','ferry','LNG'} and not (set(edges.a)|set(edges.b))-{*ids}
dist=pd.read_csv(o/'island_host_network_pairs.csv')
km={(r.island_a,r.island_b):r.port_distance_km for r in dist.itertuples()}
km.update({(j,i):v for (i,j),v in list(km.items())})
rows=[]
for last in [2021,2023]:
    early=set(annual[annual.year.le(2019)].landmass_id)
    late=set(annual[annual.year.between(last+1,2024)].landmass_id)
    for include in ['all','container_ferry']:
        n=nodes[nodes.year.le(last)]
        e=edges[edges.year.le(last)]
        if include=='container_ferry':
            n=n[n.kind.ne('LNG')];e=e[e.kind.ne('LNG')]
        send=set(n[n.direction.eq('out')].island)
        for i,g in n[n.direction.eq('in')].groupby('island'):
            if g.year.nunique()<2 or i not in early&late or eff[i]<=0:continue
            incoming=e[e.b.eq(i)]
            same_registry=set(n[n.direction.eq('out') & n.source.isin(g.source)].island)
            for k,name in sp.items():
                yr=first.get((k,i),9999)
                if yr<=last or name not in menu:continue
                orig={j for j in ids if first.get((k,j),9999)<=2019}
                if not orig&send:continue
                links=incoming[incoming.a.isin(orig)]
                host=menu[name]&alien.get(i,set())
                rows.append(dict(last=last,network=include,island=i,species=name,km=min(km[i,j] for j in orig),
                                 route=int(not links.empty),host=int(bool(host)),hosts='|'.join(sorted(host)),
                                 kinds='|'.join(sorted(set(links.kind))),linked_islands='|'.join(map(str,sorted(set(links.a.astype(int))))),
                                 same_registry=bool(orig&same_registry),event=int(yr<=2024),year=yr if yr<9999 else np.nan))
p=pd.DataFrame(rows)
summ=[];matched=[]
def overlap(q,by_island):
    groups=q.groupby('island') if by_island else [(0,q)]
    keep=[]
    for _,g in groups:
        if g.route.nunique()!=2:continue
        by=g.groupby('route').km
        z=g[g.km.between(by.min().max(),by.max().min())]
        if z.route.nunique()==2:keep.extend(z.index)
    return keep
for (last,net,host),g in p.groupby(['last','network','host']):
    for scope,q in [('all',g),('common_distance',g.loc[overlap(g,False)]),('same_island_distance',g.loc[overlap(g,True)]),
                    ('same_registry_island_distance',g.loc[overlap(g[g.same_registry],True)])]:
        for r in [0,1]:
            sub=q[q.route.eq(r)]
            summ.append(dict(last=last,network=net,host=host,scope=scope,route=r,n=len(sub),events=int(sub.event.sum()),
                             islands=sub.island.nunique(),rate=sub.event.mean(),min_km=sub.km.min(),max_km=sub.km.max()))
        if scope=='same_island_distance' and len(q):
            for i,z in q.groupby('island'):
                left=z[z.route.eq(1)];right=z[z.route.eq(0)]
                for t in left.itertuples():
                    # Same distance on the same island can have several valid comparison species.
                    err=np.abs(np.log(right.km/t.km))
                    controls=right[np.isclose(err,err.min())]
                    matched.append(dict(last=last,network=net,host=host,island=i,species=t.species,km=t.km,
                                        control_km=controls.km.mean(),control_species='|'.join(controls.species),
                                        event=t.event,control_rate=controls.event.mean(),distance_ratio=np.exp(err.min())))
t=pd.DataFrame(summ);t.to_csv(out/'comparisons.csv',index=False)
pd.DataFrame(matched).to_csv(out/'nearest_controls.csv',index=False)
p.to_csv(out/'species_islands.csv',index=False)
primary=p[p['last'].eq(2023)&p.network.eq('all')];h=primary[primary.host.eq(1)].copy()
local=set(overlap(h,True)); pooled=set(overlap(h,False))
near=[]
for (last,network,host,i),g in p.groupby(['last','network','host','island']):
    controls=g[g.route.eq(0)]
    if controls.empty:continue
    for z in g[g.route.eq(1)].itertuples():
        dx=np.abs(np.log(controls.km/z.km))
        c=controls[np.isclose(dx,dx.min())]
        near.append(dict(last=last,network=network,host=host,island=i,species=z.species,
                         distance_ratio=np.exp(dx.min()),event=z.event,control_rate=c.event.mean(),control_n=len(c)))
pd.DataFrame(near).to_csv(out/'within_island_nearest.csv',index=False)
print('all',len(primary),primary.island.nunique(),primary.species.nunique(),'host',len(h),h.island.nunique())
print(t.query("last==2023 and network=='all' and host==1").to_string(index=False))
print(h[['island','species','km','route','kinds','event']].to_string(index=False))
plt.rcParams['font.family']='Yu Gothic'
fig,ax=plt.subplots(figsize=(10,5.7),layout='constrained')
labels={9:'ボルネオ',13:'本州',21:'ニューファンドランド',27:'北海道',28:'サハリン',43:'九州',44:'台湾',79:'ケープブレトン',152:'マヨルカ',169:'フュン',715:'アルス',1577:'エーロ'}
islands=sorted(h.island.unique());rng=np.random.default_rng(28)
if pooled:
    by=h.groupby('route').km
    ax.axvspan(by.min().max(),by.max().min(),color='#edc574',alpha=.18,zorder=0)
for j,i in enumerate(islands):
    g=h[h.island.eq(i)]
    sub=g.loc[g.index.intersection(local)]
    if len(sub):
        ax.plot([sub.km.min(),sub.km.max()],[j,j],color='#4e836a',lw=5,alpha=.35,zorder=0)
    for (r,dist),z in g.groupby(['route','km']):
        offsets=np.linspace(-.10,.10,len(z)) if len(z)>1 else np.zeros(1)
        ax.scatter([dist]*len(z),j+(.16 if r else -.16)+offsets,marker='s' if r else 'o',
                   c=np.where(z.event,'#bf4633','#9ba9b5'),s=49,edgecolors='white',linewidths=.5)
ax.set_xscale('log')
ax.set_yticks(range(len(islands)),[labels.get(i,str(i)) for i in islands])
ax.set_xlabel('以前の同種記録がある最寄りの島までの距離（km、対数軸）')
from matplotlib.ticker import FuncFormatter
ticks=[k for k in [20,50,100,200,400,600,1000,2000,4000,8000,16000] if h.km.min()*.85<=k<=h.km.max()*1.2]
ax.set_xticks(ticks,labels=[f'{k:,}' for k in ticks]);ax.set_xlim(h.km.min()*.85,h.km.max()*1.2)
ax.grid(axis='x',alpha=.18);ax.spines[['top','right']].set_visible(False)
ax.set_title('外来食草の先行記録がある種×島：2020–2023年の航路と2024年の初記録')
from matplotlib.lines import Line2D
items=[Line2D([],[],ls='',marker=m,color=c,label=t) for m,c,t in [('s','#606e7a','航路確認'),('o','#606e7a','航路未確認'),('o','#bf4633','赤＝新規記録あり'),('o','#9ba9b5','灰＝新規記録なし')]]
if local:items.append(Line2D([],[],color='#9bbb9e',lw=4,label='同じ島の共通距離域（採用点）'))
if pooled:
    from matplotlib.patches import Patch
    items.append(Patch(facecolor='#faecd0',label='全島をまとめた共通距離域'))
fig.legend(handles=items,loc='outside lower center',ncol=3,frameon=False,fontsize=9)
fig.savefig(out/'overlap.png',dpi=170)
