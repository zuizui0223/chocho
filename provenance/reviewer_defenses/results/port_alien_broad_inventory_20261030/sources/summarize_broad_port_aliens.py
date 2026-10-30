"""Unadjusted broad port host inventory; pending queries are not absences."""
from pathlib import Path
import json
import pandas as pd
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs'
p=pd.read_parquet(W/'port_named_host_regional_status.parquet')
complete=pd.read_csv(O/'port_flora_regional_support.csv')
ports=pd.DataFrame([dict(port_id=x['properties']['OBJECTID'],country=x['properties']['COUNTRY']) for x in json.loads((W/'ports_2019.geojson').read_text())['features']])
a=p[p.regional_catalogue_status.eq('contemporary_not_native_in_catalogue')]
q=complete[['port_id']].merge(ports,on='port_id',validate='one_to_one')
q['has_alien_candidate']=q.port_id.isin(a.port_id)
q['has_classified_host']=q.port_id.isin(p[p.regional_catalogue_status.isin(['native_in_catalogue','contemporary_not_native_in_catalogue'])].port_id)
s=q.groupby('country',dropna=False).agg(complete_ports=('port_id','size'),ports_with_classified_hosts=('has_classified_host','sum'),ports_with_alien_candidates=('has_alien_candidate','sum')).sort_values('complete_ports',ascending=False)
s.to_csv(O/'port_alien_host_broad_country_coverage.csv')
t=a.groupby('host_name').agg(ports=('port_id','nunique')).sort_values('ports',ascending=False)
t.to_csv(O/'port_alien_host_broad_species.csv')
counts=pd.crosstab(p.port_id,p.regional_catalogue_status)
counts=counts.reindex(complete.port_id,fill_value=0).reset_index()
allports=ports.merge(counts,on='port_id',how='left',validate='one_to_one')
allports['query_complete']=allports.port_id.isin(complete.port_id)
# Missing cells outside complete queries remain missing; zeros mean no records.
allports.to_csv(O/'port_alien_host_broad_all_ports.csv',index=False)
assert len(allports)==3630 and allports.port_id.is_unique
assert allports.loc[~allports.query_complete,counts.columns[1:]].isna().all().all()
assert s.complete_ports.sum()==len(complete)
assert s.ports_with_alien_candidates.sum()==a.port_id.nunique()
assert len(t)==a.host_name.nunique() and t.ports.sum()==len(a)
receipt=dict(status='PASS_BROAD_INVENTORY_AGGREGATES',registered_ports=len(allports),complete_ports=len(complete),country_codes=len(s),ports_with_alien_candidates=int(a.port_id.nunique()),alien_candidate_species=len(t),port_alien_host_pairs=len(a),unknown_origin_port_host_pairs=int(p.regional_catalogue_status.isin(['unknown_catalogue_status','no_region_intersection','ambiguous_across_regions']).sum()),checks=['country totals reconcile','species totals reconcile','all registered ports retained','pending queries retain missing values'],limits=['record inventory not introduction mechanism','no nonport comparison or effort adjustment','port buffers overlap','regional contemporary minus native is candidate status only'])
(O/'port_alien_host_broad_validation.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt))
