from pathlib import Path
import pandas as pd,json,shapely
r=Path(__file__).resolve().parents[2];w=r/'work/barriers';o=r/'outputs'
b=pd.read_parquet(w/'butterfly_records_in_complete_plant_port_queries.parquet')
for pid,q in b.groupby('port_id'):
    a=json.loads((w/f'port_plant_coordinates/{pid}/status.json').read_text())
    assert shapely.covers(shapely.from_wkt(a['query']['geometry']),shapely.points(q.decimalLongitude,q.decimalLatitude)).all()
assert not b.duplicated(['gbifID','port_id']).any()
s=pd.read_csv(o/'port_alien_host_butterfly_all225_summary.csv')
x=pd.read_csv(o/'port_alien_host_butterfly_corecording.csv')
edges=pd.read_parquet(w/'all225_named_host_edges_local.parquet')
p=pd.read_parquet(w/'port_host_occurrence_distance_joined.parquet')
m=pd.read_parquet(w/'port_plant_occurrence_membership.parquet')
alien=set(p.loc[p.host_origin.eq('NONNATIVE_CANDIDATE'),'key'])
hostname=p.set_index('key').host_name.to_dict()
resources={pid:{hostname[k] for k in q.key if k in alien} for pid,q in m.groupby('query_port_id')}
network={name:set(q.host_name) for name,q in edges.groupby('species')}
keys=s.set_index('species').taxon_key.to_dict()
obs=b.groupby(['speciesKey','port_id']).gbifID.nunique().to_dict()
for z in x.itertuples():
    assert len(network.get(z.species,set()) & resources.get(z.port_id,set()))==z.recorded_alien_candidate_host_species
    assert obs.get((keys[z.species],z.port_id),0)==z.butterfly_records
    assert bool(z.both_recorded)==(z.recorded_alien_candidate_host_species>0 and z.butterfly_records>0)
for z in s.itertuples():
    q=x[x.species.eq(z.species)]
    assert q.both_recorded.sum()==z.ports_with_both
    assert q.butterfly_records.gt(0).sum()==z.ports_with_butterfly_records
    assert q.recorded_alien_candidate_host_species.gt(0).sum()==z.ports_with_alien_candidate_hosts
a={'status':'PASS','all_species_port_rows_checked':len(x),'coordinate_memberships_checked':len(b),'checks':['all retained butterfly coordinates covered by query polygon','unique record-port memberships','all host counts recomputed using set intersections','all butterfly counts and joint-recording flags','all225 species summaries']}
(o/'port_alien_host_butterfly_corecording_validation.json').write_text(json.dumps(a,indent=2))
print(json.dumps(a))
