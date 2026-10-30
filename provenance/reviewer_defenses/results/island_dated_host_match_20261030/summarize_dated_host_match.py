from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers'
inventory=json.loads((O/'island_port_early_inventory_audit.json').read_text())
assert inventory['status']=='COMPLETE_DATED_INVENTORY'
dated=json.loads((O/'island_host_family_null_dated_audit.json').read_text())
pooled=json.loads((O/'island_host_family_null_audit.json').read_text())
observed=pd.read_csv(O/'island_host_family_null_dated_observed_pairs.csv',dtype={'island_a':str,'island_b':str})
h=pd.read_parquet(W/'port_named_hosts_2001_2010.parquet')
assert h.year_start.eq(2001).all() and h.year_end.eq(2010).all()
edge=pd.read_parquet(W/'all225_named_host_edges_local.parquet')[['species','host_name']].drop_duplicates()
diet={sp:set(g.host_name) for sp,g in edge.groupby('species')}
resources={}
for lid,g in h.groupby('landmass_id'):
    states=g.groupby('host_name').regional_catalogue_status.agg(set)
    resources[str(lid)]=set(states[states.map(lambda x:x=={'contemporary_not_native_in_catalogue'})].index)
matched=[];links=[]
for r in observed.itertuples():
    names=diet.get(r.species,set())&resources.get(r.island_a,set())&resources.get(r.island_b,set())
    matched.append(bool(names))
    for host in names:links.append(dict(species=r.species,island_a=r.island_a,island_b=r.island_b,host_name=host))
np.testing.assert_array_equal(matched,observed.matching_shared_early_host)
assert len(observed)==1163
assert not observed.duplicated(['species','island_a','island_b']).any()
plants=pd.DataFrame(links,columns=['species','island_a','island_b','host_name'])
plants.to_csv(O/'island_dated_host_match_links.csv',index=False)
top=plants.groupby('host_name').agg(butterfly_species=('species','nunique'),species_pair_links=('species','size')).reset_index()
top=top.sort_values(['butterfly_species','species_pair_links','host_name'],ascending=[False,False,True])
top.to_csv(O/'island_dated_host_match_top_plants.csv',index=False)
c=pd.read_csv(O/'forecast_candidate_disposition.csv')[['species','taxon_key']]
summary=observed.groupby('species').agg(both_new_pairs=('species','size'),matching_early_pairs=('matching_shared_early_host','sum'),fraction_matching_early=('matching_shared_early_host','mean'))
c=c.merge(summary,left_on='species',right_index=True,how='left',validate='one_to_one')
c['comparison_evaluable']=c.both_new_pairs.notna();assert len(c)==225
c.to_csv(O/'island_dated_host_match_all225.csv',index=False)
rows=[]
for label,a in [('Plants2001-2024',pooled),('Plants2001-2010',dated)]:
    rows.append(dict(plant_window=label,observed=a['observed_fraction'],null_mean=a['null_mean_fraction'],
                     observed_species_equal=a['observed_species_equal_fraction'],null_species_equal=a['null_mean_species_equal'],
                     draws_at_least_observed=a['draws_at_least_observed'],draws_at_least_species_equal=a['draws_at_least_observed_species_equal']))
s=pd.DataFrame(rows);s.to_csv(O/'island_dated_host_match_comparison.csv',index=False)
v=dict(status='PASS_DATED_MATCHES_RECOMPUTED',records_checked=len(observed),matched_pairs=int(sum(matched)),
       matched_species=int(observed.loc[observed.matching_shared_early_host,'species'].nunique()),
       all225_accounted=True,changed_count_ports=inventory['changed_count_ports'],
       early_host_port_pairs=len(h),early_host_species=h.host_name.nunique(),
       early_alien_candidate_species=h.loc[h.regional_catalogue_status.eq('contemporary_not_native_in_catalogue'),'host_name'].nunique(),
       matched_pairs_with_no_pooled_match=int((observed.matching_shared_early_host&~observed.shared_matching_host).sum()),
       checks=['all1163 matches reconstructed via set intersections','host years2001-2010','no duplicated species-island pairs','all225 candidate accounting'],
       limits='Recorded temporal precedence only. Missing plant records are not absence. Static catalogue and2019 port locations are retrospective. Does not validate forecasting or transport.')
(O/'island_dated_host_match_validation.json').write_text(json.dumps(v,indent=2),encoding='utf-8')
fig,axs=plt.subplots(1,2,figsize=(9,4.6))
for ax,metric,null,title in [(axs[0],'observed','null_mean','Species-pair weighted'),(axs[1],'observed_species_equal','null_species_equal','Species weighted equally')]:
    x=np.arange(2);ax.bar(x-.17,s[metric]*100,width=.34,label='Observed',color='#177e89');ax.bar(x+.17,s[null]*100,width=.34,label='Family-null mean',color='#b8cece')
    ax.set_xticks(x,['Plants2001-2024','Plants2001-2010']);ax.set_ylabel('Matching shared host (%)');ax.set_title(title);ax.legend(fontsize=8)
fig.text(.02,.015,'Same 51 islands, 82 butterfly species, 1,163 shared-new species-pair records. No arrival or transport inference.',fontsize=8)
fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(O/'island_dated_host_match.png',dpi=180);plt.close(fig)
print(s.to_string(index=False));print(json.dumps(v));print(top.head(8).to_string(index=False))
