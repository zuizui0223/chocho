from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,time
import pandas as pd,requests
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers';D=W/'port_flora_2001_2010';D.mkdir(exist_ok=True)
plan=pd.read_csv(O/'island_port_early_plant_availability.csv',dtype={'landmass_id':str})
assert len(plan)==194 and plan.port_id.nunique()==194
lock=json.loads((O/'island_port_early_inventory_lock.json').read_text())
for idx,r in enumerate(plan.itertuples(),1):
    path=D/f'{r.port_id}.json'
    if path.exists():
        saved=json.loads(path.read_text(encoding='utf-8'))
        if saved['status'] in ['COMPLETE_EARLY_FACETS','CACHED_ZERO_EARLY_RECORDS']:continue
    source=json.loads((W/'port_flora_inventory'/f'{r.port_id}.json').read_text(encoding='utf-8'))
    assert source['status']=='COMPLETE_FACETS'
    out=dict(port_id=int(r.port_id),landmass_id=r.landmass_id,year_range=[2001,2010],
             original_cached_early_count=int(r.early_plant_records),original_response_sha256=source['response_sha256'],
             checked_at=datetime.now(timezone.utc).isoformat(),pages=[],species_counts=[])
    if r.early_plant_records==0:
        out.update(status='CACHED_ZERO_EARLY_RECORDS',count=0)
    else:
        params=[(k,v) for k,v in source['query'] if k not in ['year','facet','facetOffset','facetLimit']]
        params.extend([('year','2001,2010'),('facet','speciesKey'),('facetLimit',1000)])
        out['query']=params;offset=0;seen=set();total=None
        try:
            while True:
                time.sleep(3)
                res=requests.get('https://api.gbif.org/v1/occurrence/search',params=params+[('facetOffset',offset)],timeout=45)
                if res.status_code==429:
                    out.update(status='RATE_LIMITED',retry_after=res.headers.get('Retry-After'));break
                res.raise_for_status();data=res.json()
                if total is not None and data['count']!=total:
                    out.update(status='SOURCE_CHANGED_DURING_PAGING');break
                total=data['count'];out['count']=total
                values=next(x['counts'] for x in data['facets'] if x['field']=='SPECIES_KEY') if data.get('facets') else []
                assert not seen.intersection(x['name'] for x in values)
                seen.update(x['name'] for x in values);out['species_counts'].extend(values)
                out['pages'].append(dict(offset=offset,rows=len(values),response_sha256=hashlib.sha256(res.content).hexdigest()))
                if len(values)<1000:
                    out.update(status='COMPLETE_EARLY_FACETS',count_change_from_cache=total-int(r.early_plant_records));break
                offset+=1000
            assert sum(x['count'] for x in out['species_counts'])<=out.get('count',float('inf'))
        except Exception as ex:out.update(status='ERROR',error=str(ex)[:400])
    path.write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(dict(index=idx,total=len(plan),port_id=int(r.port_id),status=out['status'],count=out.get('count'))),flush=True)
    if out['status']=='RATE_LIMITED':raise SystemExit(2)
