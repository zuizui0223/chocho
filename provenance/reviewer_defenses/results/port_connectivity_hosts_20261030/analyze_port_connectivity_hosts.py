from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers'
files=[O/'port_connectivity_host_lock.json',O/'port_connectivity_exact_matches.csv',W/'unctad_port_connectivity_2006_2021.parquet',
       O/'port_alien_host_broad_all_ports.csv',W/'port_named_host_regional_status.parquet',O/'port_buffer_landmass_audit.csv',
       W/'all225_named_host_edges_local.parquet',O/'forecast_candidate_disposition.csv']
lock=json.loads(files[0].read_text());m=pd.read_csv(files[1],keep_default_na=False)
s=pd.read_parquet(files[2]);s=s[s.year.between(2006,2010)]
index=s.groupby('series_code').agg(quarters=('value','count'),connectivity=('value','mean')).reset_index()
m=m.merge(index,on='series_code',validate='one_to_one');m=m[m.quarters.ge(16)]
inventory=pd.read_csv(files[3],keep_default_na=False);complete=set(inventory.loc[inventory.query_complete.astype(str).eq('True'),'port_id'])
m=m[m.port_id.isin(complete)].copy();h=pd.read_parquet(files[4]);h=h[h.port_id.isin(m.port_id)]
geo=pd.read_csv(files[5],dtype={'single_landmass_id':str})
geo=geo[geo.status.eq('SINGLE_LANDMASS')&geo.islands_intersecting_5km_buffer.eq(1)]
m=m.merge(geo[['port_id','single_landmass_id']],on='port_id',how='left',validate='one_to_one')
def counts(g):
    states=g.groupby('host_name').regional_catalogue_status.agg(set)
    alien=sum(x=={'contemporary_not_native_in_catalogue'} for x in states)
    native=sum(x=={'native_in_catalogue'} for x in states)
    return dict(alien_hosts=alien,native_hosts=native,alien_share=alien/(alien+native) if alien+native else np.nan)
byport={pid:counts(g) for pid,g in h.groupby('port_id')}
for field in ['alien_hosts','native_hosts','alien_share']:
    m[field]=m.port_id.map(lambda pid:byport.get(pid,dict(alien_hosts=0,native_hosts=0,alien_share=np.nan))[field])
m.to_csv(O/'port_connectivity_host_port_rows.csv',index=False)
island=[]
for lid,g in m.dropna(subset=['single_landmass_id']).groupby('single_landmass_id'):
    row=dict(landmass_id=lid,ports=len(g),connectivity=g.connectivity.max());row.update(counts(h[h.port_id.isin(g.port_id)]));island.append(row)
i=pd.DataFrame(island);i.to_csv(O/'port_connectivity_host_island_rows.csv',index=False)
rows=[]
for label,d in [('all_completed_ports',m),('single_island_ports',m[m.single_landmass_id.notna()]),('island_union_max_connectivity',i)]:
    for outcome in ['alien_hosts','alien_share']:
        q=d.dropna(subset=['connectivity',outcome]);rho=float(spearmanr(q.connectivity,q[outcome]).statistic) if len(q)>2 and q[outcome].nunique()>1 else np.nan
        rows.append(dict(cohort=label,outcome=outcome,n=len(q),spearman_rho=rho))
result=pd.DataFrame(rows);result.to_csv(O/'port_connectivity_host_summary.csv',index=False)
edges=pd.read_parquet(files[6]);candidate=pd.read_csv(files[7])[['species','taxon_key']];assert len(candidate)==225
alien=h[h.regional_catalogue_status.eq('contemporary_not_native_in_catalogue')].merge(edges,on='host_name',validate='many_to_many')
for label,ids in [('all_ports',set(m.port_id)),('island_ports',set(m.loc[m.single_landmass_id.notna(),'port_id']))]:
    candidate[label+'_with_matching_alien_host']=candidate.species.map(alien[alien.port_id.isin(ids)].groupby('species').port_id.nunique()).fillna(0).astype(int)
candidate.to_csv(O/'port_connectivity_host_all225.csv',index=False)
audit=dict(status='EXPLORATORY_ARCHIVE_NAME_MATCH_ASSOCIATION_NOT_TRANSPORT_EFFECT',all_ports=len(m),island_ports=int(m.single_landmass_id.notna().sum()),
           islands=len(i),species_linked_island_hosts=int(candidate.island_ports_with_matching_alien_host.gt(0).sum()),
           country_counts=m.groupby('country_code').size().to_dict(),limits=lock['limits'],
           input_sha256={str(f.relative_to(R)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
(O/'port_connectivity_host_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
fig,axes=plt.subplots(1,2,figsize=(9,4.6));z=m[m.single_landmass_id.notna()]
for ax,col,label in [(axes[0],'alien_hosts','Recorded alien candidate host richness'),(axes[1],'alien_share','Alien share of classified host species')]:
    ax.scatter(z.connectivity,z[col],s=35,color='#177e89',alpha=.75);ax.set_xlabel('Mean connectivity index, 2006–2010');ax.set_ylabel(label);ax.grid(alpha=.15)
fig.suptitle(f'{len(z)} island ports: archived connectivity vs recorded host plants')
fig.text(.02,.015,'Plants pooled2001–2024; name-matched ports. No causal interpretation; no adjustment for recording or geography.',fontsize=8)
fig.tight_layout(rect=[0,.05,1,.93]);fig.savefig(O/'port_connectivity_host.png',dpi=180);plt.close(fig)
print(result.to_string(index=False));print(json.dumps({k:v for k,v in audit.items() if k!='input_sha256'}))
