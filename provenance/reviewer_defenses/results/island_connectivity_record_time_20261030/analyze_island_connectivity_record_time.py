from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers'
files=[O/'island_connectivity_record_time_lock.json',O/'port_connectivity_island_support.csv',O/'global_original_species_island_history.csv',O/'forecast_candidate_disposition.csv']
lock=json.loads(files[0].read_text());ports=pd.read_csv(files[1],dtype={'single_landmass_id':str})
ports=ports[ports.available_quarters.ge(16)]
exposure=ports.groupby('single_landmass_id').agg(connectivity=('mean_2006_2010','max'),matched_ports=('port_id','size')).rename_axis('landmass_id').reset_index()
h=pd.read_csv(files[2],dtype={'landmass_id':str});c=pd.read_csv(files[3])[['species','taxon_key']];assert len(c)==225
prior=set(h.loc[h.first_record_year.le(2010),'landmass_id'])
exposure['prior_butterfly_record']=exposure.landmass_id.isin(prior)
exposure.to_csv(O/'island_connectivity_record_time_support.csv',index=False)
d=h[h.first_record_year.le(2022)].merge(exposure[exposure.prior_butterfly_record],on='landmass_id',validate='many_to_one')
d=d.merge(c,left_on='speciesKey',right_on='taxon_key',validate='many_to_one');d['new_record']=d.first_record_year.ge(2011)
assert not d.duplicated(['speciesKey','landmass_id']).any()
d.to_csv(O/'island_connectivity_record_time_rows.csv',index=False)
islands=d.groupby('landmass_id').agg(connectivity=('connectivity','first'),recorded_species=('species','nunique'),new_species=('new_record','sum'),new_share=('new_record','mean')).reset_index()
islands.to_csv(O/'island_connectivity_record_time_islands.csv',index=False)
rows=[]
for r in c.itertuples():
    q=d[d.species.eq(r.species)];valid=len(q)>=5 and q.new_record.nunique()>1 and q.connectivity.nunique()>1
    rows.append(dict(species=r.species,taxon_key=r.taxon_key,islands=len(q),new_islands=int(q.new_record.sum()),
                     rho=float(spearmanr(q.connectivity,q.new_record.astype(int)).statistic) if valid else np.nan,
                     disposition='EVALUABLE' if valid else 'NO_RECORDED_SUPPORT' if not len(q) else 'INSUFFICIENT_ISLANDS_OR_VARIATION'))
species=pd.DataFrame(rows);species.to_csv(O/'island_connectivity_record_time_all225.csv',index=False)
valid=species.dropna(subset=['rho'])
audit=dict(status='RECORD_CHRONOLOGY_ASSOCIATION_NOT_COLONIZATION',index_supported_islands=len(exposure),
           islands=len(islands),species=d.species.nunique(),species_island_records=len(d),new_records=int(d.new_record.sum()),
           island_new_share_rho=float(spearmanr(islands.connectivity,islands.new_share).statistic),
           within_species_evaluable=len(valid),mean_species_rho=float(valid.rho.mean()),median_species_rho=float(valid.rho.median()),
           positive_species=int(valid.rho.gt(0).sum()),limits=lock['limits'],
           input_sha256={str(f.relative_to(R)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
(O/'island_connectivity_record_time_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
# Independent summary checks against stored row-level values.
assert len(species)==225 and species.islands.sum()==len(d)
for r in valid.itertuples():
    q=d[d.species.eq(r.species)];actual=q.assign(new_numeric=q.new_record.astype(int))[['connectivity','new_numeric']].rank().corr().iloc[0,1]
    np.testing.assert_allclose(actual,r.rho,atol=1e-12)
audit['rank_recalculation_validation']='PASS'
(O/'island_connectivity_record_time_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in audit.items() if k!='input_sha256'}))
