"""Retrospective island-pair composition association; no inferred routes."""
from pathlib import Path
from itertools import combinations
import hashlib,json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr,rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers'
files=[O/'island_host_network_lock.json',O/'port_buffer_landmass_audit.csv',O/'port_alien_host_broad_all_ports.csv',
       W/'port_named_host_regional_status.parquet',O/'global_original_species_island_history.csv',
       O/'forecast_candidate_disposition.csv',W/'ports_2019.geojson']
lock=json.loads(files[0].read_text());assert lock['candidate_population']==225
geo=pd.read_csv(files[1],dtype={'single_landmass_id':str});ports=pd.read_csv(files[2],keep_default_na=False)
complete=set(ports.loc[ports.query_complete.astype(str).eq('True'),'port_id'])
geo=geo[geo.status.eq('SINGLE_LANDMASS') & geo.islands_intersecting_5km_buffer.eq(1) & geo.port_id.isin(complete)]
mapping=geo.set_index('port_id').single_landmass_id
h=pd.read_parquet(files[3]);h=h[h.port_id.isin(mapping.index)].copy();h['island']=h.port_id.map(mapping)
hostsets={}
for lid,g in h.groupby('island'):
    ss=g.groupby('host_name').regional_catalogue_status.agg(set)
    hostsets[lid]=set(ss[ss.map(lambda v:v=={'contemporary_not_native_in_catalogue'})].index)
b=pd.read_csv(files[4],dtype={'landmass_id':str})
c=pd.read_csv(files[5])[['species','taxon_key']];assert len(c)==225
assert set(b.speciesKey)<=set(c.taxon_key)
early={lid:set(g.speciesKey) for lid,g in b[b.first_record_year.le(2010)].groupby('landmass_id')}
new={lid:set(g.speciesKey) for lid,g in b[b.first_record_year.between(2011,2022)].groupby('landmass_id')}
cat=json.loads(files[6].read_text(encoding='utf-8'))
coords={f['properties']['OBJECTID']:f['geometry']['coordinates'] for f in cat['features']}
representative={}
for lid,g in geo.groupby('single_landmass_id'):
    ll=np.radians(np.array([coords[i] for i in g.port_id]));lon,lat=ll.T
    xyz=np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]).mean(axis=1)
    xyz/=np.linalg.norm(xyz);representative[lid]=xyz
rows=[];islands=[]
for lid in sorted(set(mapping)):
    eligible=bool(hostsets.get(lid)) and bool(early.get(lid))
    islands.append(dict(landmass_id=lid,ports=int((mapping==lid).sum()),alien_host_species=len(hostsets.get(lid,set())),
                        earlier_butterfly_species=len(early.get(lid,set())),new_butterfly_species=len(new.get(lid,set())),
                        eligible=eligible))
islands=pd.DataFrame(islands);ids=islands.loc[islands.eligible,'landmass_id'].tolist()
def jac(a,b):return len(a&b)/len(a|b) if a|b else np.nan
for a,z in combinations(ids,2):
    aa=new.get(a,set());zz=new.get(z,set())
    rows.append(dict(island_a=a,island_b=z,host_jaccard=jac(hostsets[a],hostsets[z]),
                     earlier_butterfly_jaccard=jac(early[a],early[z]),new_butterfly_jaccard=jac(aa,zz),
                     shared_new_species=len(aa&zz),union_new_species=len(aa|zz),
                     port_distance_km=6371.0088*np.arccos(np.clip(representative[a]@representative[z],-1,1))))
d=pd.DataFrame(rows);q=d.dropna(subset=['new_butterfly_jaccard']).copy()
assert len(q)>3
x=rankdata(q.host_jaccard);y=rankdata(q.new_butterfly_jaccard)
design=np.column_stack([np.ones(len(q)),rankdata(np.log1p(q.port_distance_km)),rankdata(q.earlier_butterfly_jaccard)])
xr=x-design@np.linalg.lstsq(design,x,rcond=None)[0];yr=y-design@np.linalg.lstsq(design,y,rcond=None)[0]
q['host_rank_residual']=xr;q['new_butterfly_rank_residual']=yr
degree=pd.concat([q.island_a,q.island_b]).value_counts()
weight=.5*(q.island_a.map(1/degree)+q.island_b.map(1/degree))
def wcorr(a,b,w):
    a=a-np.average(a,weights=w);b=b-np.average(b,weights=w)
    return float(np.sum(w*a*b)/np.sqrt(np.sum(w*a*a)*np.sum(w*b*b)))
stats=dict(status='RETROSPECTIVE_ISLAND_PAIR_ASSOCIATION_NOT_TRANSPORT',eligible_islands=len(ids),all_pairs=len(d),
           usable_pairs=len(q),empty_new_union_pairs=int(d.new_butterfly_jaccard.isna().sum()),
           spearman_host_new=float(spearmanr(q.host_jaccard,q.new_butterfly_jaccard).statistic),
           spearman_host_earlier=float(spearmanr(q.host_jaccard,q.earlier_butterfly_jaccard).statistic),
           partial_rank_host_new=float(np.corrcoef(xr,yr)[0,1]),
           endpoint_equal_rank_correlation=wcorr(x,y,weight),
           limits=lock['inference'],input_sha256={str(f.relative_to(R)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
islands.to_csv(O/'island_host_network_islands.csv',index=False)
d.to_csv(O/'island_host_network_pairs.csv',index=False)
q.to_csv(O/'island_host_network_rank_residuals.csv',index=False)
c['earlier_records_on_eligible_islands']=c.taxon_key.map(b[b.landmass_id.isin(ids)&b.first_record_year.le(2010)].groupby('speciesKey').size()).fillna(0).astype(int)
c['new_records_on_eligible_islands']=c.taxon_key.map(b[b.landmass_id.isin(ids)&b.first_record_year.between(2011,2022)].groupby('speciesKey').size()).fillna(0).astype(int)
c.to_csv(O/'island_host_network_all225.csv',index=False)
(O/'island_host_network_audit.json').write_text(json.dumps(stats,indent=2),encoding='utf-8')
fig,axes=plt.subplots(1,2,figsize=(10,4.7))
for ax,xx,yy,title in [(axes[0],q.host_jaccard,q.new_butterfly_jaccard,f'Raw association: rho={stats["spearman_host_new"]:.3f}'),
                       (axes[1],xr,yr,f'After distance and earlier composition: r={stats["partial_rank_host_new"]:.3f}')]:
    ax.scatter(xx,yy,s=9,alpha=.25,color='#177e89',edgecolors='none');ax.set_title(title,fontsize=10)
    ax.grid(alpha=.15)
axes[0].set(xlabel='Alien host Jaccard similarity',ylabel='New butterfly Jaccard similarity',xlim=(-.02,1.02),ylim=(-.02,1.02))
axes[1].set(xlabel='Host similarity: rank residual',ylabel='New butterfly similarity: rank residual')
fig.suptitle(f'{len(ids)} islands; {len(q):,} island pairs (pairs are not independent)',fontsize=12)
fig.text(.02,.015,'New records: 2011–2022; earlier records: 2001–2010. Plants pooled 2001–2024. Not a dispersal network or a forecast.',fontsize=8)
fig.tight_layout(rect=[0,.04,1,.94]);fig.savefig(O/'island_host_network.png',dpi=180);plt.close(fig)
print(json.dumps({k:v for k,v in stats.items() if k!='input_sha256'}))
