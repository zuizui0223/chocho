#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,time,urllib.parse,urllib.request
from pathlib import Path
from shapely.geometry import shape
from shapely import to_wkt

GBIF='https://api.gbif.org/v1'

def get_json(path,params,timeout=45,retries=3):
    url=f"{GBIF}/{path}?"+urllib.parse.urlencode(params)
    last=None
    for attempt in range(retries):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'chocho-temporal-lag-audit/0.1'})
            with urllib.request.urlopen(req,timeout=timeout) as r:
                return json.load(r), url
        except Exception as e:
            last=e
            if attempt+1<retries: time.sleep(2**attempt)
    raise RuntimeError(f'GBIF request failed: {last}; URL length={len(url)}')

def resolve(name):
    x,_=get_json('species/match',{'name':name,'rank':'SPECIES'},timeout=30)
    key=x.get('usageKey'); canonical=str(x.get('canonicalName') or '')
    rank=str(x.get('rank') or '').upper(); mt=str(x.get('matchType') or '').upper()
    confidence=x.get('confidence')
    ok=bool(key) and rank=='SPECIES' and canonical
    return {'usage_key':key,'canonical_name':canonical,'rank':rank,'match_type':mt,'confidence':confidence,'accepted':ok}

def geometry_by_code(path):
    x=json.loads(path.read_text())
    out={}
    for f in x.get('features',[]):
        code=str((f.get('properties') or {}).get('LEVEL3_COD') or '').strip()
        if code and f.get('geometry'): out[code]=shape(f['geometry'])
    return out

def query_years(taxon_key,geom):
    last=None
    for tol in (0.01,0.03,0.05,0.1,0.2):
        g=geom.simplify(tol,preserve_topology=True)
        wkt=to_wkt(g,rounding_precision=4,trim=True)
        params={'taxonKey':int(taxon_key),'geometry':wkt,'hasCoordinate':'true','hasGeospatialIssue':'false','occurrenceStatus':'PRESENT','year':'1750,2026','facet':'year','facetLimit':300,'facetMincount':1,'limit':0}
        url=f"{GBIF}/occurrence/search?"+urllib.parse.urlencode(params)
        if len(url)>12000:
            last=f'URL too long at tol={tol}: {len(url)}'; continue
        try:
            x,_=get_json('occurrence/search',params,timeout=60,retries=3)
        except Exception as e:
            last=str(e); continue
        years=[]
        for facet in x.get('facets') or []:
            if str(facet.get('field') or '').upper()=='YEAR':
                for c in facet.get('counts') or []:
                    try: years.append((int(c['name']),int(c['count'])))
                    except Exception: pass
        years=sorted(years)
        return {'count':int(x.get('count') or 0),'years':years,'geometry_simplify_tolerance':tol,'request_url_length':len(url)}
    raise RuntimeError(last or 'unable to query GBIF years')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--event-host-csv',type=Path,required=True)
    ap.add_argument('--level3-geojson',type=Path,required=True)
    ap.add_argument('--output-csv',type=Path,required=True)
    ap.add_argument('--output-json',type=Path,required=True)
    a=ap.parse_args()
    geom=geometry_by_code(a.level3_geojson)
    rows=list(csv.DictReader(a.event_host_csv.open(newline='',encoding='utf-8')))
    resolved={}; out=[]
    for i,r in enumerate(rows,1):
        name=r['host_name'].strip(); code=r['wgsrpd3_code'].strip()
        if name not in resolved: resolved[name]=resolve(name)
        meta=resolved[name]
        q={'count':0,'years':[],'geometry_simplify_tolerance':None,'request_url_length':None}; error=''
        if not meta['accepted']: error='REJECTED_TAXON_MATCH'
        elif code not in geom: error='MISSING_WGSRPD3_GEOMETRY'
        else:
            try:q=query_years(meta['usage_key'],geom[code])
            except Exception as e:error=str(e)
        years=[y for y,n in q['years'] if 1750<=y<=2026 and n>0]
        earliest=min(years) if years else None
        butterfly=int(r['butterfly_first_record_year'])
        out.append({**r,'gbif_host_usage_key':meta.get('usage_key'),'gbif_host_canonical_name':meta.get('canonical_name'),'gbif_host_match_type':meta.get('match_type'),'gbif_host_match_confidence':meta.get('confidence'),'gbif_coordinate_records_in_region':q['count'],'gbif_earliest_host_record_year':earliest,'recorded_host_to_butterfly_lag_years':None if earliest is None else butterfly-earliest,'host_record_precedes_butterfly':None if earliest is None else int(earliest<butterfly),'host_record_by_2017':None if earliest is None else int(earliest<=2017),'geometry_simplify_tolerance':q['geometry_simplify_tolerance'],'request_url_length':q['request_url_length'],'error':error})
        print(i,len(rows),name,code,'earliest',earliest,'butterfly',butterfly,'error',error,flush=True)
    fields=list(out[0])
    a.output_csv.parent.mkdir(parents=True,exist_ok=True)
    with a.output_csv.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
    cells={}
    for r in out:
        key=(r['species'],r['wgsrpd3_code'],int(r['butterfly_first_record_year']))
        y=r['gbif_earliest_host_record_year']
        if y not in ('',None): cells[key]=min(cells.get(key,9999),int(y))
    cell_rows=[]
    for (sp,code,bfy),hy in sorted(cells.items()):
        cell_rows.append({'species':sp,'wgsrpd3_code':code,'butterfly_first_record_year':bfy,'earliest_recorded_introduced_known_host_year':hy,'recorded_resource_lead_years':bfy-hy,'host_record_by_2017':hy<=2017})
    payload={'schema':'chocho_butterfly_event_host_first_record_feasibility_v0.1','status':'POSTHOC_STAGE_B_FEASIBILITY_AUDIT','host_region_queries':len(out),'queries_with_dated_host_record':sum(r['gbif_earliest_host_record_year'] not in ('',None) for r in out),'event_cells':len(set((r['species'],r['wgsrpd3_code']) for r in out)),'event_cells_with_dated_host_record':len(cell_rows),'host_first_event_cells':sum(r['recorded_resource_lead_years']>0 for r in cell_rows),'same_year_event_cells':sum(r['recorded_resource_lead_years']==0 for r in cell_rows),'butterfly_first_event_cells':sum(r['recorded_resource_lead_years']<0 for r in cell_rows),'event_cells_host_record_by_2017':sum(r['host_record_by_2017'] for r in cell_rows),'cells':cell_rows,'claim_boundary':['GBIF first record is a detection date, not an establishment date.','This audit is conditioned on the nine Stage-A introduced-resource cells with subsequent butterfly detections and is descriptive, not an effect estimate.','WCVP classifies these host-region combinations as introduced; GBIF occurrence records are used only to date recorded presence.']}
    a.output_json.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(payload,indent=2))
if __name__=='__main__':main()
