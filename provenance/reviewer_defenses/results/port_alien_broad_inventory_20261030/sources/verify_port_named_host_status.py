from pathlib import Path
import json
import pandas as pd
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs'
p=pd.read_parquet(W/'port_named_host_regional_status.parquet')
names=pd.read_parquet(W/'host_flora_name_join_audit.parquet').set_index('host_name')
checked=0
for pid,q in p.groupby('port_id'):
    a=json.loads((W/'port_flora_inventory'/f'{pid}.json').read_text())
    assert a['status']=='COMPLETE_FACETS'
    counts={int(x['name']):x['count'] for x in a['facets']['SPECIES_KEY']}
    for r in q.itertuples():
        assert counts[r.gbif_species_key]==r.plant_records
        assert names.loc[r.host_name,'accepted_gbif_species_key']==r.gbif_species_key
        assert r.source_response_sha256==a['response_sha256']
        checked+=1
s=pd.read_csv(O/'port_named_host_status_summary.csv').set_index('regional_catalogue_status')
for status,q in p.groupby('regional_catalogue_status'):
    assert s.loc[status,'port_host_pairs']==len(q)
    assert s.loc[status,'ports']==q.port_id.nunique()
    assert s.loc[status,'host_species']==q.gbif_species_key.nunique()
    assert s.loc[status,'plant_records']==q.plant_records.sum()
receipt=dict(status='PASS',checked_port_host_pairs=checked,
             checks=['all counts and source hashes match cached complete GBIF facets',
                     'all names match accepted species key audit', 'all published summary counts recomputed'],
             not_verified=['independent introduced status', 'true plant abundance or absence',
                           'host arrival chronology','independent accuracy of source regional catalogue'])
(O/'port_named_host_status_validation.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt))
