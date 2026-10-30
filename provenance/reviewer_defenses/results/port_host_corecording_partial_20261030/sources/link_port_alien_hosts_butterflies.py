"""Same queried port neighborhood co-recording; not interaction or arrival evidence."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import shapely

R=Path(__file__).resolve().parents[2]; W=R/'work/barriers'; O=R/'outputs'
audit=json.loads((O/'port_host_distance_join_audit.json').read_text())
sources=json.loads((W/'port_host_distance_source_pages.json').read_text())
portids=sorted({x['port_id'] for x in sources})
assert len(portids)==audit['complete_ports']
geoms=[]
for pid in portids:
    status=json.loads((W/f'port_plant_coordinates/{pid}/status.json').read_text())
    assert status['status']=='COMPLETE'
    geoms.append(shapely.from_wkt(status['query']['geometry']))
tree=shapely.STRtree(geoms)
parts=sorted((R/'work/temporal/occurrence_parts').glob('part_*.parquet'))
rows=[]; scanned=0
for i,f in enumerate(parts):
    q=pd.read_parquet(f,columns=['gbifID','speciesKey','decimalLongitude','decimalLatitude','year'])
    scanned+=len(q)
    q=q[q.year.between(2001,2024)&q.decimalLongitude.between(-180,180)&q.decimalLatitude.between(-90,90)].reset_index(drop=True)
    point_idx,poly_idx=tree.query(shapely.points(q.decimalLongitude,q.decimalLatitude),predicate='intersects')
    if len(point_idx):
        x=q.iloc[point_idx].copy(); x['port_id']=np.array(portids)[poly_idx]; rows.append(x)
    if (i+1)%10==0:print(f'scanned {i+1}/{len(parts)} parts',flush=True)
assert rows
b=pd.concat(rows,ignore_index=True).drop_duplicates(['gbifID','port_id'])
b.to_parquet(W/'butterfly_records_in_complete_plant_port_queries.parquet',index=False)
candidates=pd.read_csv(O/'forecast_candidate_disposition.csv')[['species','taxon_key']]
assert len(candidates)==225 and candidates.taxon_key.is_unique
b=b.merge(candidates.rename(columns={'taxon_key':'speciesKey'}),on='speciesKey',validate='many_to_one')
plants=pd.read_parquet(W/'port_host_occurrence_distance_joined.parquet')
members=pd.read_parquet(W/'port_plant_occurrence_membership.parquet').rename(columns={'query_port_id':'port_id'})
alien=plants[plants.host_origin.eq('NONNATIVE_CANDIDATE')][['key','host_name']]
atports=alien.merge(members,on='key',validate='one_to_many')[['port_id','host_name']].drop_duplicates()
edges=pd.read_parquet(W/'all225_named_host_edges_local.parquet')[['species','host_name']].drop_duplicates()
linked=atports.merge(edges,on='host_name',validate='many_to_many')
resources=linked.groupby(['species','port_id']).host_name.nunique().rename('recorded_alien_candidate_host_species').reset_index()
observed=b.groupby(['species','port_id']).agg(butterfly_records=('gbifID','nunique'),first_record_year=('year','min'),last_record_year=('year','max')).reset_index()
grid=pd.MultiIndex.from_product([candidates.species,portids],names=['species','port_id']).to_frame(index=False)
out=grid.merge(resources,how='left',on=['species','port_id'],validate='one_to_one').merge(observed,how='left',on=['species','port_id'],validate='one_to_one')
for col in ['recorded_alien_candidate_host_species','butterfly_records']:out[col]=out[col].fillna(0).astype(int)
out['both_recorded']=out.recorded_alien_candidate_host_species.gt(0)&out.butterfly_records.gt(0)
out.to_csv(O/'port_alien_host_butterfly_corecording.csv',index=False)
summary=out.groupby('species').agg(ports_with_alien_candidate_hosts=('recorded_alien_candidate_host_species',lambda x:int(x.gt(0).sum())),ports_with_butterfly_records=('butterfly_records',lambda x:int(x.gt(0).sum())),ports_with_both=('both_recorded','sum')).reset_index()
summary=candidates.merge(summary,on='species',validate='one_to_one')
assert len(out)==225*len(portids) and len(summary)==225
assert summary.ports_with_both.le(summary[['ports_with_alien_candidate_hosts','ports_with_butterfly_records']].min(axis=1)).all()
summary.to_csv(O/'port_alien_host_butterfly_all225_summary.csv',index=False)
receipt=dict(status='PARTIAL_PORT_SAME_NEIGHBORHOOD_CORECORDING_NOT_INTERACTION',plant_snapshot_unique_records=audit['unique_records'],candidate_species=225,complete_plant_ports=len(portids),butterfly_source_parts=len(parts),butterfly_records_scanned=scanned,butterfly_unique_records_in_queries=int(b.gbifID.nunique()),species_with_alien_candidate_host_records=int(summary.ports_with_alien_candidate_hosts.gt(0).sum()),species_with_butterfly_records=int(summary.ports_with_butterfly_records.gt(0).sum()),species_with_both_at_same_port=int(summary.ports_with_both.gt(0).sum()),species_port_pairs_with_both=int(out.both_recorded.sum()),limits=['original coordinates intersect the exact same 5km query polygon as plants','pooled2001-2024 co-recording does not establish simultaneous occurrence or host use','zero means no record in this snapshot, never ecological absence','overlapping port memberships are retained; unique butterfly count is deduplicated','partial17port cohort is not a global estimate','no matched controls or observation effort adjustment'])
(O/'port_alien_host_butterfly_corecording_audit.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt),flush=True)
