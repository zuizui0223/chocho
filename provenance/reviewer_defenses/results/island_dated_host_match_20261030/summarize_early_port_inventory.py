"""Audit dated host inventory, keeping unqueried ports distinct from zero records."""
from pathlib import Path
import json,hashlib
import pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers';S=R/'work/remote_sensing/host_source/extracted/inputs'
plan=pd.read_csv(O/'island_port_early_plant_availability.csv',dtype={'landmass_id':str})
regions=pd.read_csv(O/'port_flora_regional_support.csv',keep_default_na=False).set_index('port_id').candidate_regions.to_dict()
names=pd.read_csv(S/'original/insect_host_accepted.csv',dtype=str)[['accepted_plant_name_id','accepted_name']].drop_duplicates()
mapping=dict(names.itertuples(index=False,name=None))
for butterfly in json.loads((S/'bce/accepted_host_edges.json').read_text())['per_butterfly']:
    for item in butterfly['details']:
        if item.get('match_status')=='ACCEPTED_ID_ABSENT_FROM_FROZEN_HOSTS_FOR_BUTTERFLY':mapping[str(item['accepted_plant_name_id'])]=item['accepted_name']
def regional(a,b):
    x=pd.concat([pd.read_csv(S/a,dtype=str),pd.read_csv(S/b,dtype=str)])
    return set(zip(x.accepted_plant_name_id.map(mapping),x.area_code_l3))
native=regional('original/native_extant_nondoubtful_wgsrpd3.csv','bce/additional_native_wgsrpd3.csv')
present=regional('original/extant_nondoubtful_wgsrpd3.csv','bce/additional_contemporary_wgsrpd3.csv')
lookup=pd.read_parquet(W/'host_flora_name_join_audit.parquet');lookup=lookup[lookup.join_status.eq('SPECIES_FACET_COMPATIBLE')]
assert not lookup.accepted_gbif_species_key.duplicated().any()
keynames=dict(zip(lookup.accepted_gbif_species_key.astype(int),lookup.host_name))
ports=[];hosts=[];hashes={}
for r in plan.itertuples():
    file=W/'port_flora_2001_2010'/f'{r.port_id}.json';record=None
    if file.exists():
        try:record=json.loads(file.read_text(encoding='utf-8'))
        except json.JSONDecodeError:pass
    if r.early_plant_records==0:
        status='CACHED_ZERO_EARLY_RECORDS';count=0;delta=0;values=[]
        original=W/'port_flora_inventory'/f'{r.port_id}.json'
        hashes[str(original.relative_to(R)).replace('\\','/')]=hashlib.sha256(original.read_bytes()).hexdigest()
        source=json.loads(original.read_text(encoding='utf-8'))
        assert sum(x['count'] for x in source['facets']['YEAR'] if 2001<=int(x['name'])<=2010)==0
    elif record and record.get('status')=='COMPLETE_EARLY_FACETS':
        status=record['status'];count=record['count'];delta=count-int(r.early_plant_records);values=record['species_counts']
        assert ['year','2001,2010'] in record['query']
        assert record['pages'][-1]['rows']<1000
        assert len(values)==len({x['name'] for x in values})
        assert sum(x['count'] for x in values)<=count
        hashes[str(file.relative_to(R)).replace('\\','/')]=hashlib.sha256(file.read_bytes()).hexdigest()
    else:
        status=record.get('status','PENDING') if record else 'PENDING';count=None;delta=None;values=[]
    complete=status in ['CACHED_ZERO_EARLY_RECORDS','COMPLETE_EARLY_FACETS']
    ports.append(dict(port_id=r.port_id,landmass_id=r.landmass_id,status=status,complete=complete,
                      original_cached_count=int(r.early_plant_records),early_plant_records=count,count_change=delta))
    codes=regions[r.port_id].split(';') if regions[r.port_id] else []
    for item in values:
        key=int(item['name'])
        if key not in keynames:continue
        name=keynames[key]
        states=['native_in_catalogue' if (name,c) in native else 'contemporary_not_native_in_catalogue' if (name,c) in present else 'unknown_catalogue_status' for c in codes]
        st=states[0] if states and len(set(states))==1 else 'ambiguous_across_regions' if states else 'no_region_intersection'
        hosts.append(dict(port_id=r.port_id,landmass_id=r.landmass_id,host_name=name,gbif_species_key=key,
                          plant_records=item['count'],regional_catalogue_status=st,year_start=2001,year_end=2010))
p=pd.DataFrame(ports);p.to_csv(O/'island_port_early_inventory_status.csv',index=False)
h=pd.DataFrame(hosts,columns=['port_id','landmass_id','host_name','gbif_species_key','plant_records','regional_catalogue_status','year_start','year_end'])
assert not h.duplicated(['port_id','gbif_species_key']).any()
h.to_parquet(W/'port_named_hosts_2001_2010.parquet',index=False)
islands=p.groupby('landmass_id').agg(planned_ports=('port_id','size'),complete_ports=('complete','sum'),all_ports_complete=('complete','all'))
islands.to_csv(O/'island_port_early_inventory_islands.csv')
audit=dict(status='COMPLETE_DATED_INVENTORY' if p.complete.all() else 'PARTIAL_DATED_INVENTORY_NO_ECOLOGICAL_COMPARISON',
           planned_ports=len(p),complete_ports=int(p.complete.sum()),live_query_completed_ports=int(p.status.eq('COMPLETE_EARLY_FACETS').sum()),
           cached_zero_ports=int(p.status.eq('CACHED_ZERO_EARLY_RECORDS').sum()),
           completely_resolved_islands=int(islands.all_ports_complete.sum()),
           changed_count_ports=int(p.count_change.fillna(0).ne(0).sum()),host_port_pairs=len(h),
           status_counts=p.status.value_counts().to_dict(),snapshot_sha256=hashes,
           limits='No ecological absence from zero records. Geographic acquisition order is nonrandom; do not compare incomplete subset as global result. Contemporary alien status is not historical introduction status. Counts across ports may overlap.')
(O/'island_port_early_inventory_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in audit.items() if k!='snapshot_sha256'}))
