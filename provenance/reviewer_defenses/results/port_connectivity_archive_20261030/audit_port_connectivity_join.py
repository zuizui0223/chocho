from pathlib import Path
import sys,json,unicodedata,re
import pandas as pd
R=Path(__file__).resolve().parents[2];sys.path.insert(0,str(R/'work/deps'))
import pycountry
O=R/'outputs';W=R/'work/barriers'
d=pd.read_csv(O/'port_connectivity_archive_catalogue.csv',keep_default_na=False)
def norm(s):return re.sub('[^A-Z0-9]','',unicodedata.normalize('NFKD',s).encode('ascii','ignore').decode().upper())
def country(s):
    try:return pycountry.countries.lookup(s).alpha_2
    except LookupError:return None
d['country_code']=d.country_label.map(country);d['port_normalized']=d.port_label.map(norm)
features=json.loads((W/'ports_2019.geojson').read_text())['features']
p=pd.DataFrame([dict(port_id=f['properties']['OBJECTID'],country_code=f['properties']['COUNTRY'],port_name=f['properties']['PORT_NAME']) for f in features])
p['port_normalized']=p.port_name.map(norm)
pc=p.groupby(['country_code','port_normalized']).size();dc=d.dropna(subset=['country_code']).groupby(['country_code','port_normalized']).size()
unique_p=p[[pc.get((r.country_code,r.port_normalized),0)==1 for r in p.itertuples()]]
unique_d=d[[pd.notna(r.country_code) and dc.get((r.country_code,r.port_normalized),0)==1 for r in d.itertuples()]]
matches=unique_p.merge(unique_d,on=['country_code','port_normalized'],validate='one_to_one')
assert matches.port_id.nunique()==len(matches) and matches.series_code.nunique()==len(matches)
matches.to_csv(O/'port_connectivity_exact_matches.csv',index=False)
q=p.merge(matches[['port_id','series_code']],on='port_id',how='left',validate='one_to_one')
q['status']=q.series_code.notna().map({True:'UNIQUE_COUNTRY_NAME_MATCH_REQUIRES_SPATIAL_CONFIRMATION',False:'UNMATCHED_OR_AMBIGUOUS_NOT_ZERO'})
q.to_csv(O/'port_connectivity_all3630_disposition.csv',index=False)
series=pd.read_parquet(W/'unctad_port_connectivity_2006_2021.parquet')
historical=series[series.year.between(2006,2010)]
stats=historical.groupby('series_code').agg(available_quarters=('value','count'),mean_2006_2010=('value','mean')).reset_index()
matches=matches.merge(stats,on='series_code',validate='one_to_one')
geo=pd.read_csv(O/'port_buffer_landmass_audit.csv',dtype={'single_landmass_id':str})
single=geo[geo.status.eq('SINGLE_LANDMASS')&geo.islands_intersecting_5km_buffer.eq(1)]
islandmatches=matches.merge(single[['port_id','single_landmass_id']],on='port_id',validate='one_to_one')
inventory=pd.read_csv(O/'port_alien_host_broad_all_ports.csv',keep_default_na=False)
completed=set(inventory.loc[inventory.query_complete.astype(str).eq('True'),'port_id'])
islandmatches['plant_inventory_complete']=islandmatches.port_id.isin(completed)
target=pd.read_csv(O/'island_host_network_islands.csv',dtype={'landmass_id':str});target=set(target.loc[target.eligible,'landmass_id'])
islandmatches['in_previous51']=islandmatches.single_landmass_id.isin(target)
islandmatches.to_csv(O/'port_connectivity_island_support.csv',index=False)
complete=islandmatches[islandmatches.plant_inventory_complete&islandmatches.available_quarters.ge(16)]
audit=dict(status='NAME_JOIN_CANDIDATES_NOT_VALIDATED_PORT_IDENTITIES_OR_EFFECT',
           provider_ports=len(d),wpi_ports=len(p),exact_country_name_matches=len(matches),
           unmapped_provider_country_labels=sorted(set(d.loc[d.country_code.isna(),'country_label'])),
           single_island_port_candidates=len(islandmatches),single_island_candidates=islandmatches.single_landmass_id.nunique(),
           complete_plant_ports_with16of20_early_quarters=len(complete),
           complete_plant_islands_with16of20_early_quarters=complete.single_landmass_id.nunique(),
           previous51_islands_with16of20_early_quarters=complete.loc[complete.in_previous51,'single_landmass_id'].nunique(),
           limits='Exact normalized country+name only; no geographic identity confirmation in PLSCI archive. Port complexes and near-duplicate names need review. Missing values not zero; 16/20 is coverage description, not fitted-model selection.')
(O/'port_connectivity_join_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps(audit))
