from frame import *
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyproj import Geod

flora=pd.read_csv(out/'plant_coverage.csv')
assert set(flora.port_id)==set(ports[ports.island.isin(targets)].port_id)
plants=pd.read_csv(out/'prior_hosts.csv').groupby('island').host.agg(set).to_dict()
eff=flora.groupby('island').plant_records.sum()
old=set(pd.read_csv(o/'island_port_early_inventory_status.csv').landmass_id)
menu=set(pd.read_csv(o/'port_alien_host_broad_species.csv').host_name)
diets=pd.read_parquet(w/'all225_named_host_edges_local.parquet')
diets=diets[diets.host_name.isin(menu)].groupby('species').host_name.agg(set).to_dict()
assert len(diets)==192
# Same port-coordinate method on every island, including potential sources without an observed route.
centers={}
for i,g in ports.groupby('island'):
    ll=np.radians(np.array([registered[k]['geometry']['coordinates'][:2] for k in g.port_id]))
    lon,lat=ll.T
    x,y,z=np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]).mean(axis=1)
    centers[i]=np.degrees([np.arctan2(y,x),np.arctan2(z,np.hypot(x,y))])
geod=Geod(ellps='WGS84')
km={(i,j):geod.inv(*centers[i],*centers[j])[2]/1000 for i in targets for j in islands if i!=j}
first=annual.groupby(['speciesKey','island']).year.min().to_dict();send=set(nodes[nodes.direction.eq('out')].island)
cases=[]
for i in sorted(targets):
    if eff[i]<=0:continue
    e=edges[edges.b.eq(i)]
    src_tables=set(nodes.loc[nodes.island.eq(i)&nodes.direction.eq('in'),'source'])
    source_domain=set(nodes.loc[nodes.source.isin(src_tables)&nodes.direction.eq('out'),'island'])
    for k,name in candidate.items():
        yr=first.get((k,i),9999)
        orig={j for j in islands if first.get((k,j),9999)<=2018}
        if name not in diets or yr<=2023 or not orig&send:continue
        recorded=diets[name]&plants.get(i,set())
        hits=e[e.a.isin(orig)]
        nearest=min(orig,key=lambda j:(km[i,j],j))
        cases.append(dict(island=i,original_51=i in old,key=k,species=name,nearest_island=nearest,km=km[i,nearest],
                          host=int(bool(recorded)),hosts='|'.join(sorted(recorded)),route=int(not hits.empty),
                          linked_sources='|'.join(map(str,sorted(set(hits.a.astype(int))))),kind='|'.join(sorted(set(hits.kind))),
                          same_registry=bool(orig&source_domain),event=int(yr==2024)))
p=pd.DataFrame(cases);assert not p.duplicated(['key','island']).any()
p.to_csv(out/'species_islands.csv',index=False)
h=p[p.host.eq(1)].copy()
def pairs(data,ratio,group='island',tie_side=None,sign=0):
    rows=[];pool=[]
    for _,g in ([(None,data)] if group=='any' else data.groupby(group)):
        yes=g[g.route.eq(1)].sort_values('key');no=g[g.route.eq(0)].sort_values('key')
        if yes.empty or no.empty:continue
        dif=np.abs(np.log(yes.km.to_numpy()[:,None]/no.km.to_numpy()))
        allowed=dif<=np.log(ratio)+1e-10
        cost=np.where(allowed,dif,1e6)
        a,b=linear_sum_assignment(cost)
        if tie_side:
            score=yes.event.to_numpy()[:,None] if tie_side=='yes' else no.event.to_numpy()[None,:]
            aa,bb=linear_sum_assignment(cost+np.where(allowed,sign*1e-8*score,0))
            valid=allowed[a,b];alt=allowed[aa,bb]
            assert alt.sum()==valid.sum() and abs(dif[a[valid],b[valid]].sum()-dif[aa[alt],bb[alt]].sum())<1e-11
            a,b=aa,bb
        rows.extend((yes.index[x],no.index[y]) for x,y in zip(a,b) if allowed[x,y])
        pool.extend(yes.index[allowed.any(axis=1)]);pool.extend(no.index[allowed.any(axis=0)])
    return rows,pool
summ=[];matches=[]
for label,data in [('all',h),('near_250km',h[h.km.le(250)]),('same_registry',h[h.same_registry]),('outside_old51',h[~h.original_51])]:
    for cutoff,group in [(r,g) for r in [1.,1.10,1.25,1.50] for g in ['any','species','island']]:
        restriction=group
        chosen,pool=pairs(data,cutoff,group)
        for group,indices in [('all_eligible',data.index),('has_counterpart',pool),('paired',[j for pair in chosen for j in pair])]:
            q=data.loc[list(indices)]
            for r in [0,1]:
                z=q[q.route.eq(r)]
                summ.append(dict(sample=label,ratio=cutoff,selection=group,restriction=restriction,route=r,n=len(z),events=int(z.event.sum()),rate=z.event.mean(),islands=z.island.nunique()))
        for a,b in chosen:
            matches.append(dict(sample=label,ratio=cutoff,restriction=restriction,yes_island=data.at[a,'island'],no_island=data.at[b,'island'],yes_species=data.at[a,'species'],no_species=data.at[b,'species'],
                                yes_km=data.at[a,'km'],no_km=data.at[b,'km'],yes_event=data.at[a,'event'],no_event=data.at[b,'event']))
s=pd.DataFrame(summ);s.to_csv(out/'matched_rates.csv',index=False)
pd.DataFrame(matches).to_csv(out/'pairs.csv',index=False)
coverage=[]
for i in sorted(islands):
    q=p[p.island.eq(i)]
    coverage.append(dict(island=i,old51=i in old,ports=ports.island.eq(i).sum(),two_year_inbound=i in dest,
                         early_butterfly=i in prior,late_butterfly=i in later,early_plants=eff.get(i,np.nan),
                         n=len(q),host_positive=int(q.host.sum()),events=int(q.event.sum())))
pd.DataFrame(coverage).to_csv(out/'island_coverage.csv',index=False)
panel=s[(s.ratio==1.25)&s['sample'].isin(['all','near_250km'])&s.selection.eq('paired')]
ties=[]
for label,g in [('all',h),('near_250km',h[h.km.le(250)])]:
    for group in ['any','species','island']:
        for route,side in [(0,'no'),(1,'yes')]:
            counts=[]
            for sign in [1,-1]:
                pairs0,_=pairs(g,1.25,group,side,sign)
                counts.append(sum(int(g.at[a if route else b,'event']) for a,b in pairs0))
            ties.append(dict(sample=label,restriction=group,route=route,n=len(pairs0),min_events=counts[0],max_events=counts[1]))
tie=pd.DataFrame(ties);tie.to_csv(out/'tie_sensitivity.csv',index=False)
print('frame',len(islands),len(dest),len(targets),p.island.nunique(),p.species.nunique(),len(p),len(h),h.island.nunique())
print(panel.to_string(index=False))
plt.rcParams['font.family']='Yu Gothic'
fig,axes=plt.subplots(1,2,sharey=True,figsize=(11.5,5),layout='constrained')
for ax,scope,title in zip(axes,['all','near_250km'],['全距離から比較可能なペアを抽出','既往記録島が250 km以内の組に限定']):
    for r,c,dx,name in [(0,'#91a1ac',-.19,'航路未確認'),(1,'#286e77',.19,'航路確認')]:
        z=panel[panel.route.eq(r)&panel['sample'].eq(scope)].set_index('restriction').loc[['any','species','island']]
        for x,row in enumerate(z.itertuples()):
            if row.n:
                ax.bar(x+dx,100*row.rate,.36,color=c)
                rng=tie[(tie['sample']==scope)&(tie.restriction==row.Index)&(tie.route==r)].iloc[0]
                lo,hi=100*rng.min_events/row.n,100*rng.max_events/row.n
                label=f'{row.rate:.1%}\n{row.events}/{row.n}'
                if lo!=hi:
                    ax.errorbar(x+dx,(lo+hi)/2,(hi-lo)/2,fmt='none',color='#43474a',capsize=5)
                    label=f'{lo:.0f}–{hi:.0f}%\n{rng.min_events}–{rng.max_events}/{row.n}'
                ax.annotate(label,(x+dx,hi),xytext=(0,4),textcoords='offset points',ha='center',fontsize=9)
            elif r==1:ax.annotate('ペアなし',(x,0),xytext=(0,4),textcoords='offset points',ha='center',fontsize=9)
        ax.bar([],[],color=c,label=name)
    ax.set_xticks(range(3),['距離のみ','同じ種＋距離','同じ島＋距離'])
    ax.set(xlim=(-.6,2.6),ylim=(0,max(10,1.5*panel.rate.max()*100)),title=title)
    ax.spines[['top','right']].set_visible(False)
axes[0].set_ylabel('2024年のチョウ新規記録率（%）')
fig.suptitle('外来食草の先行記録がある種×島：距離比1.25倍以内で1対1照合')
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=c,label=n) for c,n in [('#91a1ac','航路未確認'),('#286e77','航路確認')]],loc='outside lower center',ncol=2,frameon=False)
fig.savefig(out/'matched_record_rates.png',dpi=175)
