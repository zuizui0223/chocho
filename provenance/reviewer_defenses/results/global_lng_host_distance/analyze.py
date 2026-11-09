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
last=int(sys.argv[1]) if len(sys.argv)>1 else 2021
assert last in [2021,2022,2023]
if last!=2021:out=out/f"through_{last}"
w=Path('work/barriers');o=Path('outputs')
edges=pd.read_csv(out/'links.csv')
endpoint=pd.read_csv(out/'endpoint_coverage.csv')
assert set(edges.year)==set(range(2020,last+1)) and not (set(edges.a)|set(edges.b))&{15,26,456}
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
eligible_islands=set(annual[annual.year.le(2019)].landmass_id)&set(annual[annual.year.between(last+1,2024)].landmass_id)&set(eff[eff.gt(0)].index)
d=pd.read_csv(o/'island_host_network_pairs.csv')
km={(r.island_a,r.island_b):r.port_distance_km for r in d.itertuples()}
km.update({(b,a):v for (a,b),v in list(km.items())})
send=set(endpoint[endpoint.direction=='out'].island)
inbound=endpoint[endpoint.direction=='in'].groupby('island').year.nunique()
rows=[]
for i,ny in inbound.items():
    if i not in eligible_islands:continue
    receives=set(edges.loc[edges.b.eq(i),'a'])
    recurrent=set(edges.loc[edges.b.eq(i)].groupby('a').year.nunique().loc[lambda x:x.ge(2)].index)
    for k,name in sp.items():
        if name not in menu or first.get((k,i),9999)<=last:continue
        orig={j for j in ids if first.get((k,j),9999)<=2019}
        if not orig&send:continue
        rows.append(dict(island=i,species=name,km=min(km[i,j] for j in orig),
                         route=int(bool(orig&receives)),host=int(bool(menu[name]&alien.get(i,set()))),
                         repeated_route=int(bool(orig&recurrent)),linked_prior='|'.join(map(str,sorted(int(j) for j in orig&receives))),
                         host_names='|'.join(sorted(menu[name]&alien.get(i,set()))),
                         first_year=first.get((k,i),np.nan),sources=len(orig&send),years=ny,event=int(first.get((k,i),9999)<=2024)))
p=pd.DataFrame(rows)
p.to_csv(out/'species_islands.csv',index=False)
print('islands',p.groupby(['island','years']).event.agg(['size','sum']))
print(p.groupby(['route','host']).event.agg(['size','sum']))
assert not p.duplicated(['island','species']).any()
def counts(exposure):
    return p.groupby(['host',exposure]).agg(n=('event','size'),events=('event','sum'),
        islands=('island','nunique'),min_km=('km','min'),max_km=('km','max')).reset_index().rename(columns={exposure:'connected'}).assign(definition=exposure)
pd.concat([counts('route'),counts('repeated_route')]).to_csv(out/'counts.csv',index=False)
# Compare within-island route contrasts on the common continuous-distance range only.
support=[]
for (host,i),g in p.groupby(['host','island']):
    halves=[g[g.route.eq(x)] for x in (0,1)]
    left=max(z.km.min() for z in halves) if all(len(z) for z in halves) else np.nan
    right=min(z.km.max() for z in halves) if all(len(z) for z in halves) else np.nan
    supported=g[g.km.between(left,right)]
    if supported.route.nunique()!=2:supported=g.iloc[:0]
    support.append(dict(host=host,island=i,n=len(g),common_n=len(supported),common_events=int(supported.event.sum())))
pd.DataFrame(support).to_csv(out/'overlap.csv',index=False)
# A descriptive linear probability difference in percentage points; distance enters
# continuously. Estimate only strata with within-island common support on both routes.
import statsmodels.formula.api as smf
checks=[]
for host,g in p.groupby('host'):
    pieces=[]
    for i,z in g.groupby('island'):
        if z.route.nunique()<2:continue
        lo=z.groupby('route').km.min().max();hi=z.groupby('route').km.max().min()
        z=z[z.km.between(lo,hi)]
        if z.route.nunique()==2:pieces.append(z)
    if not pieces:continue
    z=pd.concat(pieces)
    fit=smf.ols('event ~ route + np.log2(km) + np.log1p(sources) + C(island)',data=z).fit()
    assert np.linalg.matrix_rank(fit.model.exog)==fit.model.exog.shape[1]
    checks.append(dict(host=host,n=len(z),islands=z.island.nunique(),events=int(z.event.sum()),
                       route_difference_pp=100*fit.params['route'],distance_doubling_pp=100*fit.params['np.log2(km)']))
pd.DataFrame(checks).to_csv(out/'distance_adjusted.csv',index=False)
print(checks)
plt.rcParams['font.family']='Yu Gothic'
fig,axes=plt.subplots(2,1,figsize=(9,5.8),sharex=True,layout='constrained')
for ax,h in zip(axes,[0,1]):
    for route,y in [(0,0),(1,1)]:
        q=p[p.host.eq(h)&p.route.eq(route)].sort_values(['km','event','species'])
        for dist,group in q.groupby('km'):
            yy=y+np.linspace(-.21,.21,len(group)) if len(group)>1 else [y]
            ax.scatter(np.full(len(group),dist),yy,c=np.where(group.event,'#bb442d','#bac6cf'),s=24,alpha=.9,edgecolors='white',linewidths=.25)
        ax.text(1.01,(y+.35)/1.70,f'{q.event.sum()}/{len(q)}',transform=ax.transAxes,ha='left',va='center',fontsize=10)
    ax.set(ylim=(-.35,1.35),yticks=[0,1],yticklabels=['航路未確認','航路確認'])
    ax.set_title('港で外来食草の先行記録'+('あり' if h else 'なし'),loc='left',fontsize=12)
    ax.spines[['top','right']].set_visible(False)
axes[1].set_xscale('log')
lo,hi=p.km.min()/1.15,p.km.max()*1.15
ticks=[x for x in [100,200,400,600,1000,2000,4000,6000,10000,20000] if lo<x<hi]
axes[1].set(xticks=ticks,xticklabels=[f'{x:,}' for x in ticks],xlim=(lo,hi),xlabel='以前の同種記録がある最寄りの島までの距離（km、対数軸）')
window=f'{last+1}–2024' if last<2023 else '2024'
fig.suptitle(f'世界のLNG船航路と島での新規記録（{window}年）')
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([],[],marker='o',linestyle='',color=c,label=t) for c,t in [('#bb442d','新規記録あり'),('#bac6cf','新規記録なし')]],loc='outside lower center',ncol=2,frameon=False)
fig.savefig(out/'continuous_distance.png',dpi=170)
panel=[]
