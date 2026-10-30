from pathlib import Path
import json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers'
d=pd.read_csv(O/'island_shared_host_match_species_pairs.csv',dtype={'island_a':str,'island_b':str},keep_default_na=False)
p=pd.read_csv(O/'island_host_network_pairs.csv',dtype={'island_a':str,'island_b':str})
h=pd.read_csv(O/'global_original_species_island_history.csv',dtype={'landmass_id':str})
h=h[h.first_record_year.between(2011,2022)]
a=p[['island_a','island_b']].merge(h[['landmass_id','speciesKey']],left_on='island_a',right_on='landmass_id').drop(columns='landmass_id')
b=p[['island_a','island_b']].merge(h[['landmass_id','speciesKey']],left_on='island_b',right_on='landmass_id').drop(columns='landmass_id')
joined=a.merge(b,on=['island_a','island_b','speciesKey'],how='outer',indicator=True,validate='one_to_one').rename(columns={'speciesKey':'taxon_key'})
check=d.merge(joined,on=['island_a','island_b','taxon_key'],how='outer',validate='one_to_one')
assert len(check)==len(d) and check.species.notna().all()
np.testing.assert_array_equal(check.new_in_both,check['_merge'].eq('both'))
edges=pd.read_parquet(W/'all225_named_host_edges_local.parquet')[['species','host_name']].drop_duplicates()
exploded=d.assign(host_name=d.matching_hosts.str.split(';')).explode('host_name')
exploded=exploded[exploded.host_name.ne('')]
tested=exploded.merge(edges,on=['species','host_name'],how='left',indicator=True,validate='many_to_one')
assert tested['_merge'].eq('both').all()
np.testing.assert_array_equal(d.matching_host_count,d.matching_hosts.map(lambda s:len(s.split(';')) if s else 0))
np.testing.assert_array_equal(d.shared_matching_host,d.matching_host_count.gt(0))
source=pd.read_parquet(W/'port_named_host_regional_status.parquet')
geo=pd.read_csv(O/'port_buffer_landmass_audit.csv',dtype={'single_landmass_id':str})
ports=pd.read_csv(O/'port_alien_host_broad_all_ports.csv',keep_default_na=False)
complete=set(ports.loc[ports.query_complete.astype(str).eq('True'),'port_id'])
geo=geo[geo.status.eq('SINGLE_LANDMASS')&geo.islands_intersecting_5km_buffer.eq(1)&geo.port_id.isin(complete)]
source=source.merge(geo[['port_id','single_landmass_id']],on='port_id',validate='many_to_one')
states=source.groupby(['single_landmass_id','host_name']).regional_catalogue_status.agg(set).to_dict()
for row in exploded.itertuples():
    assert states[(row.island_a,row.host_name)]=={'contemporary_not_native_in_catalogue'}
    assert states[(row.island_b,row.host_name)]=={'contemporary_not_native_in_catalogue'}
top=pd.read_csv(O/'island_shared_host_match_top_plants.csv')
actual=exploded[exploded.new_in_both].groupby('host_name').species.nunique()
np.testing.assert_array_equal(top.butterfly_species,actual.reindex(top.host_name))
c=pd.read_csv(O/'island_shared_host_match_all225.csv');assert len(c)==225 and c.species.nunique()==225
v=dict(status='PASS',species_pairs_checked=len(d),positive_host_links_checked=len(exploded),
       checks=['all union/both outcomes reconstructed by relational joins','each reported host edge exists','each reported host is alien candidate on both islands','counts and top-plant species richness','225 candidate accounting'],
       limitations='Does not independently establish completeness of catalogue edges, plant observations or transport; checks positive matches, not exhaustive absence of omitted matches.')
(O/'island_shared_host_match_validation.json').write_text(json.dumps(v,indent=2),encoding='utf-8');print(json.dumps(v))
