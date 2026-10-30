from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(__file__).resolve().parents[2]; O=R/'outputs'; W=R/'work/barriers'
inputs=[O/'island_host_continuity_lock.json',O/'global_original_sequence_observation_opportunity.csv',
        O/'port_buffer_landmass_audit.csv',O/'port_alien_host_broad_all_ports.csv',
        W/'port_named_host_regional_status.parquet',W/'all225_named_host_edges_local.parquet',
        W/'global_original_species_island_year.parquet',O/'forecast_candidate_disposition.csv']
lock=json.loads(inputs[0].read_text())
assert lock['index_years']==[2011,2017]
events=pd.read_csv(inputs[1],dtype={'landmass_id':str})
events=events[events.first_record_year.between(*lock['index_years'])].copy()
assert not events.duplicated(['speciesKey','landmass_id']).any()
geo=pd.read_csv(inputs[2],dtype={'single_landmass_id':str})
ports=pd.read_csv(inputs[3],keep_default_na=False)
complete=set(ports.loc[ports.query_complete.astype(str).eq('True'),'port_id'])
geo=geo[geo.status.eq('SINGLE_LANDMASS') & geo.islands_intersecting_5km_buffer.eq(1) & geo.port_id.isin(complete)]
mapping=geo.set_index('port_id').single_landmass_id
hosts=pd.read_parquet(inputs[4]); hosts=hosts[hosts.port_id.isin(mapping.index)].copy()
hosts['landmass_id']=hosts.port_id.map(mapping)
edges=pd.read_parquet(inputs[5])[['species','host_name']].drop_duplicates()
linked=hosts.merge(edges,on='host_name',validate='many_to_many')
resources={}
for key,g in linked.groupby(['species','landmass_id']):
    states=g.groupby('host_name').regional_catalogue_status.agg(set)
    resources[key]=(len(states),sum(x=={'contemporary_not_native_in_catalogue'} for x in states),
                    sum(x=={'native_in_catalogue'} for x in states))
years=pd.read_parquet(inputs[6]); years.landmass_id=years.landmass_id.astype(str)
history=years.groupby(['speciesKey','landmass_id']).year.agg(set).to_dict()
known_islands=set(mapping)
rows=[]
for e in events.itertuples():
    hs=resources.get((e.species,e.landmass_id),(0,0,0)) if e.landmass_id in known_islands else None
    group=('NOT_EVALUABLE' if hs is None else 'ALIEN_CANDIDATE_RECORDED' if hs[1] else
           'NATIVE_ONLY_CLASSIFIED' if hs[2] else 'UNRESOLVED_ONLY' if hs[0] else 'NO_MATCHING_HOST_RECORD')
    yy=history[(e.speciesKey,e.landmass_id)]
    assert min(yy)==e.first_record_year
    follow=sorted(y for y in yy if e.first_record_year<y<=e.first_record_year+5)
    rows.append(dict(species=e.species,speciesKey=e.speciesKey,landmass_id=e.landmass_id,
                     first_record_year=e.first_record_year,host_group=group,
                     matching_hosts=np.nan if hs is None else hs[0],alien_hosts=np.nan if hs is None else hs[1],
                     followup_record_years=len(follow),any_followup_record=bool(follow),
                     record_in_year5=e.first_record_year+5 in yy,
                     followup_year_list=';'.join(map(str,follow))))
d=pd.DataFrame(rows); d.to_csv(O/'island_host_continuity_events.csv',index=False)
metrics=lock['outcomes']
summary=[]
for group,g in d.groupby('host_group'):
    row=dict(host_group=group,events=len(g),species=g.species.nunique(),islands=g.landmass_id.nunique())
    for m in metrics:
        row[m+'_event_mean']=g[m].mean()
        row[m+'_species_equal_mean']=g.groupby('species')[m].mean().mean()
    summary.append(row)
s=pd.DataFrame(summary);s.to_csv(O/'island_host_continuity_summary.csv',index=False)
paired=[]
for sp,g in d.groupby('species'):
    a=g[g.host_group.eq('ALIEN_CANDIDATE_RECORDED')]; b=g[g.host_group.eq('NATIVE_ONLY_CLASSIFIED')]
    if len(a) and len(b):
        row=dict(species=sp,alien_events=len(a),native_events=len(b))
        for m in metrics: row[m+'_difference']=a[m].mean()-b[m].mean()
        paired.append(row)
p=pd.DataFrame(paired);p.to_csv(O/'island_host_continuity_within_species.csv',index=False)
c=pd.read_csv(inputs[7])[['species','taxon_key']]; assert len(c)==225
counts=d.groupby('species').size(); valid=d[d.host_group.ne('NOT_EVALUABLE')].groupby('species').size()
c['index_events']=c.species.map(counts).fillna(0).astype(int)
c['host_context_events']=c.species.map(valid).fillna(0).astype(int)
c['paired_comparison']=c.species.isin(p.species if len(p) else [])
c['disposition']=np.where(c.index_events.eq(0),'NO_INDEX_EVENT',np.where(c.host_context_events.eq(0),'NO_RESOLVED_HOST_CONTEXT','HOST_CONTEXT_AVAILABLE'))
c.to_csv(O/'island_host_continuity_all225.csv',index=False)
audit=dict(status='RETROSPECTIVE_RECORD_CONTINUITY_NOT_ESTABLISHMENT',events=len(d),species=d.species.nunique(),
           evaluated_events=int(d.host_group.ne('NOT_EVALUABLE').sum()),paired_species=len(p),
           mean_within_species_difference={m:float(p[m+'_difference'].mean()) for m in metrics} if len(p) else {},
           input_sha256={str(x.relative_to(R)).replace('\\','/'):hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs},
           limits=lock['limits'])
(O/'island_host_continuity_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
groups=['ALIEN_CANDIDATE_RECORDED','NATIVE_ONLY_CLASSIFIED','UNRESOLVED_ONLY','NO_MATCHING_HOST_RECORD','NOT_EVALUABLE']
labels=['Alien host candidate recorded','Native classified host; no alien record','Unresolved host records only','No matching host record','Host context unavailable']
fig,(ax,ax2)=plt.subplots(1,2,figsize=(13,6),gridspec_kw={'width_ratios':[1.5,1]}); q=s.set_index('host_group').reindex(groups)
for i,g in enumerate(groups):
    if pd.notna(q.loc[g,'events']):
        ax.scatter(q.loc[g,'any_followup_record_event_mean']*100,i,color='#177e89',s=65,label='Events weighted equally' if i==0 else None)
        ax.scatter(q.loc[g,'any_followup_record_species_equal_mean']*100,i,marker='D',color='#b85c38',s=45,label='Species weighted equally' if i==0 else None)
ax.set_yticks(range(5),[f'{lab}\n(n={int(q.loc[g,"events"]) if pd.notna(q.loc[g,"events"]) else 0} events)' for g,lab in zip(groups,labels)])
ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('Any subsequent record in years 1–5 (%)')
ax.set_title('All events: different species mixtures',loc='left',weight='bold',fontsize=11)
ax.grid(axis='x',alpha=.2);ax.legend(loc='lower right',fontsize=8)
for side in ['top','right']:ax.spines[side].set_visible(False)
if len(p):
    differences=p.any_followup_record_difference*100
    ax2.scatter(differences,np.arange(len(p)),s=28,color='#666666')
    ax2.axvline(0,color='#999999',lw=1)
    ax2.axvline(differences.mean(),color='#b85c38',lw=2,label=f'Mean: {differences.mean():+.1f} pp')
    ax2.set_yticks(np.arange(len(p)),[f'{r.species} ({r.alien_events}/{r.native_events})' for r in p.itertuples()],fontsize=7)
    ax2.set_xlim(-105,105);ax2.set_xlabel('Alien minus native context (percentage points)')
    ax2.set_title(f'Within-species contrasts (n={len(p)})\nLabels: alien/native event counts',loc='left',weight='bold',fontsize=10)
    ax2.legend(fontsize=8,loc='lower right')
    ax2.grid(axis='x',alpha=.2)
    for side in ['top','right']:ax2.spines[side].set_visible(False)
fig.text(.02,.015,'First records 2011–2017; pooled plant records 2001–2024. Descriptive association, not establishment or a host effect.',fontsize=8)
fig.tight_layout(rect=[0,.04,1,1]);fig.savefig(O/'island_host_continuity.png',dpi=180);plt.close(fig)
print(s.to_string(index=False));print(json.dumps({k:v for k,v in audit.items() if k!='input_sha256'}))
