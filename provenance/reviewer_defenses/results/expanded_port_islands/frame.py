from pathlib import Path
import json
import pandas as pd
import numpy as np
out=Path(__file__).resolve().parent
w=Path('work/barriers')
o=Path('outputs')
features=json.loads((w/'ports_2019.geojson').read_text(encoding='utf-8'))['features']
registered={f['properties']['OBJECTID']:f for f in features}
ports=pd.read_csv(o/'port_buffer_landmass_audit.csv').query("status=='SINGLE_LANDMASS' and islands_intersecting_5km_buffer==1").copy()
ports['island']=ports.single_landmass_id.astype(int)
excluded=set(ports.loc[ports.port_id.map(lambda i:registered[i]['properties']['COUNTRY']).isin(['GB','IE']),'island'])
islands=set(ports.island)-excluded
ports=ports[ports.island.isin(islands)]
annual=pd.read_parquet(w/'global_original_species_island_year.parquet')
annual['island']=annual.landmass_id.astype(int)
candidate=pd.read_csv(o/'forecast_candidate_disposition.csv').set_index('taxon_key').species
assert len(candidate)==225
annual=annual[annual.island.isin(islands)&annual.speciesKey.isin(candidate.index)]
nodes=pd.read_csv(out/'endpoint_coverage.csv')
edges=pd.read_csv(out/'links.csv')
ne=[];nn=[]
for r in pd.read_csv(out/'ncfo_segments.csv').itertuples():
    for a,b in [(r.a,r.b),(r.b,r.a)]:
        if a in islands:nn.append((r.year,'out',int(a)))
        if b in islands:nn.append((r.year,'in',int(b)))
        if a in islands and b in islands:ne.append((r.year,int(a),int(b)))
def bts(rows,cols):
    return pd.DataFrame(rows,columns=cols).drop_duplicates().assign(source='US BTS NCFO',kind='ferry')
nodes=pd.concat([nodes,bts(nn,['year','direction','island'])]).drop_duplicates()
edges=pd.concat([edges,bts(ne,['year','a','b'])]).drop_duplicates()
nodes=nodes[nodes.island.isin(islands)]
edges=edges[edges.a.isin(islands)&edges.b.isin(islands)]
dest=set(nodes[nodes.direction.eq('in')].groupby('island').year.nunique().loc[lambda x:x>=2].index)
prior=set(annual[annual.year.le(2018)].island);later=set(annual[annual.year.eq(2024)].island)
targets=dest&prior&later
if __name__=='__main__':
    print('frame',len(islands),len(ports),'multi-year',len(dest),'with prior and late butterfly',len(targets),'ports',ports.island.isin(targets).sum())
    print(sorted(targets))
