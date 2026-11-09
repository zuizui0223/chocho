from frame import w,o,out,ports,targets,registered
import json,hashlib,time
import requests,pandas as pd, geopandas as gpd,shapely
from pyproj import Geod
from pathlib import Path
S=Path('work/remote_sensing/host_source/extracted/inputs')
cache=w/'expanded_port_flora_2001_2010';cache.mkdir(exist_ok=True)
menu=set(pd.read_csv(o/'port_alien_host_broad_species.csv').host_name);assert len(menu)==687
lookup=pd.read_parquet(w/'host_flora_name_join_audit.parquet').query("join_status=='SPECIES_FACET_COMPATIBLE'")
keynames=lookup[lookup.host_name.isin(menu)].set_index('accepted_gbif_species_key').host_name.to_dict()
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
regions=gpd.read_file(S/'level3.geojson').to_crs(4326)
hist=pd.read_csv(o/'island_port_early_inventory_status.csv').set_index('port_id')
orig=pd.read_parquet(w/'port_named_hosts_2001_2010.parquet').query("regional_catalogue_status=='contemporary_not_native_in_catalogue'")
geod=Geod(ellps='WGS84')
session=requests.Session()
def collect(i):
    f=cache/f'{i}.json'
    lon,lat=registered[i]['geometry']['coordinates'][:2]
    ring=[geod.fwd(lon,lat,-a,5000)[:2] for a in range(0,360,15)]
    geom=shapely.Polygon(ring)
    assert geom.bounds[2]-geom.bounds[0]<180
    if f.exists():s=json.loads(f.read_text())
    else:
        whole=w/'port_flora_inventory'/f'{i}.json'
        old=json.loads(whole.read_text()) if whole.exists() else {}
        if old.get('status')=='COMPLETE_FACETS' and not any(2001<=int(x['name'])<=2010 for x in old['facets']['YEAR']):
            s=dict(count=0,keys=[],hashes=[old['response_sha256']],method='cached_year_zero')
        else:
            params=dict(taxonKey=7707728,hasCoordinate='true',hasGeospatialIssue='false',occurrenceStatus='PRESENT',
                        year='2001,2010',geometry=shapely.to_wkt(geom,rounding_precision=6),limit=0,facet='speciesKey',facetLimit=1000)
            keys=[];hashes=[];total=None
            while True:
                time.sleep(.75)
                for attempt in range(3):
                    resp=session.get('https://api.gbif.org/v1/occurrence/search',params=dict(params,facetOffset=len(keys)),timeout=45)
                    if resp.status_code==429:time.sleep(max(10,int(resp.headers.get('Retry-After',30))));continue
                    if resp.ok:break
                    time.sleep(3)
                resp.raise_for_status();j=resp.json()
                assert total is None or j['count']==total
                total=j['count'];batch=[int(x['name']) for fac in j.get('facets',[]) if fac['field']=='SPECIES_KEY' for x in fac['counts']]
                assert not set(keys).intersection(batch)
                keys.extend(batch);hashes.append(hashlib.sha256(resp.content).hexdigest())
                if len(batch)<1000:break
            s=dict(count=total,keys=keys,hashes=hashes,method='dated_facets')
        f.write_text(json.dumps(s))
    codes=set(regions.iloc[regions.sindex.query(geom,predicate='intersects')].LEVEL3_COD.dropna())
    alien={keynames[k] for k in s['keys'] if k in keynames and codes and all((keynames[k],c) in present and (keynames[k],c) not in native for c in codes)}
    return int(s['count']),alien,s['method'],'|'.join(s['hashes'])
plan=ports[ports.island.isin(targets)]
print('planned',len(plan),'new or cached outside earlier inventory',sum(~plan.port_id.isin(hist.index)),flush=True)
rows=[];hosts=[]
for r in plan.itertuples():
    i=r.port_id
    if i in hist.index:
        assert int(hist.at[i,'landmass_id'])==r.island
        count=hist.at[i,'early_plant_records'];found=set(orig.loc[orig.port_id.eq(i),'host_name'])&menu;method='prior_complete';digest=''
    else:
        count,found,method,digest=collect(i)
    rows.append(dict(port_id=i,island=r.island,plant_records=count,method=method,response_hashes=digest))
    hosts.extend(dict(port_id=i,island=r.island,host=h) for h in sorted(found))
    if len(rows)%20==0:print(len(rows),flush=True)
pd.DataFrame(rows).to_csv(out/'plant_coverage.csv',index=False)
pd.DataFrame(hosts).to_csv(out/'prior_hosts.csv',index=False)
