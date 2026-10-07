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
    resolved={}; butterfly_resolved={}; butterfly_region_first={}; out=[]
    for i,r in enumerate(rows,1):
        name=r['host_name'].strip(); code=r['wgsrpd3_code'].strip(); butterfly_name=r['species'].strip()
        bkey=(butterfly_name,code)
        if bkey not in butterfly_region_first:
            if butterfly_name not in butterfly_resolved: butterfly_resolved[butterfly_name]=resolve(butterfly_name)
            bmeta=butterfly_resolved[butterfly_name]
            if bmeta['accepted'] and code in geom:
                try:
                    bq=query_years(bmeta['usage_key'],geom[code])
                    byears=[y for y,n in bq['years'] if 1750<=y<=2026 and n>0]
                    butterfly_region_first[bkey]={'earliest':min(byears) if byears else None,'count':bq['count'],'meta':bmeta}
                except Exception as e:
                    butterfly_region_first[bkey]={'earliest':None,'count':0,'meta':bmeta,'error':str(e)}
            else:
                butterfly_region_first[bkey]={'earliest':None,'count':0,'meta':bmeta}
        binfo=butterfly_region_first[bkey]
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
        butterfly_2010plus=int(r['butterfly_first_record_year'])
        butterfly_full=binfo.get('earliest')
        out.append({**r,'gbif_butterfly_usage_key':binfo.get('meta',{}).get('usage_key'),'gbif_butterfly_coordinate_records_in_region':binfo.get('count'),'gbif_earliest_butterfly_record_year_full':butterfly_full,'butterfly_has_pre2010_record':None if butterfly_full is None else int(int(butterfly_full)<2010),'gbif_host_usage_key':meta.get('usage_key'),'gbif_host_canonical_name':meta.get('canonical_name'),'gbif_host_match_type':meta.get('match_type'),'gbif_host_match_confidence':meta.get('confidence'),'gbif_coordinate_records_in_region':q['count'],'gbif_earliest_host_record_year':earliest,'recorded_host_to_butterfly_lag_years':None if earliest is None or butterfly_full is None else int(butterfly_full)-earliest,'host_record_precedes_butterfly':None if earliest is None or butterfly_full is None else int(earliest<int(butterfly_full)),'host_record_by_2017':None if earliest is None else int(earliest<=2017),'geometry_simplify_tolerance':q['geometry_simplify_tolerance'],'request_url_length':q['request_url_length'],'error':error})
        print(i,len(rows),name,code,'host earliest',earliest,'butterfly full earliest',butterfly_full,'2010plus',butterfly_2010plus,'error',error,flush=True)
    fields=list(out[0])
    a.output_csv.parent.mkdir(parents=True,exist_ok=True)
    with a.output_csv.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
    cells={}
    for r in out:
        bfull=r['gbif_earliest_butterfly_record_year_full']
        key=(r['species'],r['wgsrpd3_code'],None if bfull in ('',None) else int(bfull),int(r['butterfly_first_record_year']))
        y=r['gbif_earliest_host_record_year']
        if y not in ('',None) and bfull not in ('',None): cells[key]=min(cells.get(key,9999),int(y))
    cell_rows=[]
    for (sp,code,bfy_full,bfy_2010plus),hy in sorted(cells.items()):
        cell_rows.append({'species':sp,'wgsrpd3_code':code,'butterfly_first_record_year_full_gbif':bfy_full,'butterfly_first_record_year_in_2010_2026_snapshot':bfy_2010plus,'butterfly_has_pre2010_record':bfy_full<2010,'earliest_recorded_introduced_known_host_year':hy,'recorded_resource_lead_years':bfy_full-hy,'host_record_by_2017':hy<=2017})
    payload={'schema':'chocho_butterfly_event_host_first_record_feasibility_v0.1','status':'POSTHOC_STAGE_B_FEASIBILITY_AUDIT','host_region_queries':len(out),'queries_with_dated_host_record':sum(r['gbif_earliest_host_record_year'] not in ('',None) for r in out),'event_cells':len(set((r['species'],r['wgsrpd3_code']) for r in out)),'event_cells_with_dated_host_record':len(cell_rows),'host_first_event_cells':sum(r['recorded_resource_lead_years']>0 for r in cell_rows),'same_year_event_cells':sum(r['recorded_resource_lead_years']==0 for r in cell_rows),'butterfly_first_event_cells':sum(r['recorded_resource_lead_years']<0 for r in cell_rows),'event_cells_host_record_by_2017':sum(r['host_record_by_2017'] for r in cell_rows),'event_cells_with_pre2010_butterfly_record':sum(r['butterfly_has_pre2010_record'] for r in cell_rows),'cells':cell_rows,'claim_boundary':['GBIF first record is a detection date, not an establishment date.','This audit is conditioned on the nine Stage-A introduced-resource cells with subsequent butterfly detections and is descriptive, not an effect estimate.','The full-history butterfly query audits whether a 2018-2025 snapshot detection was actually absent from GBIF before 2010.','WCVP classifies these host-region combinations as introduced; GBIF occurrence records are used only to date recorded presence.']}
    a.output_json.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(payload,indent=2))
if __name__=='__main__':main()
