#!/usr/bin/env python3
"""Bounded, post-hoc iNaturalist larval-presence audit of the fixed 32-species panel.

An annotated larva is NOT proof of feeding on a mapped host, successful breeding,
resource-driven expansion or persistent establishment. Unknown/obscured locations
are not upgraded. Taxa and the five-page-per-species limit are frozen in advance.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, math, time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from shapely.geometry import Point, shape
from shapely.strtree import STRtree

API='https://api.inaturalist.org/v1/'
START='2010-01-01'
END='2026-10-11'

def resolve(name, candidates):
    exact=[t for t in candidates if t.get('rank')=='species' and t.get('name','').casefold()==name.casefold()]
    if len(exact)==1: return exact[0]
    synonyms=[t for t in candidates if t.get('rank')=='species' and t.get('matched_term','').casefold()==name.casefold()]
    if len(synonyms)==1: return synonyms[0]
    raise ValueError('No unique exact scientific-name or explicitly matched synonym: '+name)

def eligible(o, taxon):
    if o.get('quality_grade')!='research': return False,'not_research_grade'
    if o.get('captive') is not False: return False,'captive_or_unknown'
    if o.get('obscured') or o.get('geoprivacy') in ('obscured','private') or o.get('taxon_geoprivacy') in ('obscured','private'):
        return False,'obscured_or_private'
    t=o.get('taxon') or {}
    if t.get('id')!=taxon and taxon not in (t.get('ancestor_ids') or []): return False,'taxon_mismatch'
    ann=o.get('annotations') or []
    if not any(a.get('controlled_attribute_id')==1 and a.get('controlled_value_id')==6 and (a.get('vote_score') or 0)>=0 for a in ann):
        return False,'no_nonnegative_larval_annotation'
    date=o.get('observed_on') or ''
    if not START<=date[:10]<=END: return False,'date_outside_frozen_window'
    xy=(o.get('geojson') or {}).get('coordinates')
    if not xy or len(xy)!=2 or not all(isinstance(v,(int,float)) and math.isfinite(v) for v in xy): return False,'no_public_coordinates'
    if not (-180<=xy[0]<=180 and -90<=xy[1]<=90): return False,'invalid_coordinates'
    return True,'eligible'

class RegionMapper:
    def __init__(self,codes,geometries):
        self.codes=codes;self.geometries=geometries;self.tree=STRtree(geometries)
        self.by_id={id(g):i for i,g in enumerate(geometries)}
    def map_point(self,lon,lat):
        p=Point(lon,lat); hits=[]
        for item in self.tree.query(p):
            i=self.by_id[id(item)] if hasattr(item,'geom_type') else int(item)
            if self.geometries[i].covers(p):hits.append(self.codes[i])
        hits=sorted(set(hits))
        return (hits[0],'unique_polygon') if len(hits)==1 else (None,'no_polygon' if not hits else 'multiple_polygons')

def summarize_record(o,species,taxon,mapper,native,contemporary):
    ok,why=eligible(o,taxon)
    if not ok:return {'observation_id':o.get('id'),'focal_species':species,'excluded_reason':why}
    lon,lat=o['geojson']['coordinates'];region,mapping=mapper.map_point(lon,lat)
    cls=('native_host_mapped' if region in native else 'introduced_only_in_frozen_hosts' if region in contemporary else 'no_known_host_mapped') if region else 'unmapped'
    accuracy=o.get('positional_accuracy')
    precision='unknown' if accuracy is None else 'known_le_10000m' if 0<=accuracy<=10000 else 'over_10000m'
    fields=[{'name':f.get('name'),'field_id':f.get('field_id'),'value':f.get('value'),'taxon_id':f.get('taxon_id')} for f in (o.get('ofvs') or [])]
    photos=[{'url':p.get('url'),'license_code':p.get('license_code')} for p in (o.get('photos') or [])]
    return {'observation_id':o['id'],'focal_species':species,'observed_taxon':o['taxon'].get('name'),
            'observed_on':o['observed_on'],'url':o.get('uri') or f"https://www.inaturalist.org/observations/{o['id']}",
            'observer_id':(o.get('user') or {}).get('id'),'longitude':lon,'latitude':lat,'positional_accuracy_m':accuracy,
            'coordinate_precision':precision,'region':region,'mapping':mapping,'resource_class':cls,
            'evidence_level':'L1_LARVAL_PRESENCE_NOT_HOST_USE','description':o.get('description'),
            'observation_fields':fields,'photos':photos,'license_code':o.get('license_code'),
            'annotations':o.get('annotations'),'quality_grade':o.get('quality_grade')}

def rows(path):
    with Path(path).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))

def load_resources(pair_csv,native_csv,contemporary_csv):
    hosts=defaultdict(set);native=defaultdict(set);current=defaultdict(set)
    for r in rows(pair_csv):hosts[r['insect_species']].add(r['accepted_plant_name_id'])
    for r in rows(native_csv):native[r['accepted_plant_name_id']].add(r['area_code_l3'])
    for r in rows(contemporary_csv):current[r['accepted_plant_name_id']].add(r['area_code_l3'])
    return hosts,native,current

class Client:
    def __init__(self):self.calls=[];self.last=0
    def get(self,endpoint,params):
        time.sleep(max(0,1.3-(time.monotonic()-self.last)))
        url=API+endpoint+'?'+urlencode(params)
        req=Request(url,headers={'User-Agent':'chocho-larval-presence-audit/20261011 (public ecological research)','Accept':'application/json'})
        self.last=time.monotonic()
        try:
            with urlopen(req,timeout=35) as f: raw=f.read()
        except Exception as ex:
            self.calls.append({'url':url,'retrieved_utc':datetime.now(timezone.utc).isoformat(),'error':str(ex)})
            raise
        payload=json.loads(raw)
        self.calls.append({'url':url,'retrieved_utc':datetime.now(timezone.utc).isoformat(),'sha256':hashlib.sha256(raw).hexdigest(),'total_results':payload.get('total_results'),'returned':len(payload.get('results',[]))})
        return payload

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sources',type=Path,required=True);ap.add_argument('--geometry',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--panel-csv',type=Path,default=Path('data/frozen/figure_sources/occurrence_resource_validation_species.csv'))
    ap.add_argument('--max-pages',type=int,default=5);args=ap.parse_args()
    if args.max_pages!=5:raise ValueError('Frozen cap is five 200-record pages per species')
    out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    panel=[r['species'] for r in rows(args.panel_csv)]
    if len(panel)!=32 or len(set(panel))!=32:raise ValueError('Fixed prior panel drift')
    hp,np,cp=load_resources(args.sources/'native/insect_host_accepted.csv',args.sources/'native/native_extant_nondoubtful_wgsrpd3.csv',args.sources/'contemporary/extant_nondoubtful_wgsrpd3.csv')
    geo=json.loads(args.geometry.read_text());features=geo['features']
    mapper=RegionMapper([f['properties']['LEVEL3_COD'] for f in features],[shape(f['geometry']) for f in features])
    client=Client();results=[];catalog=[]
    for sp in panel:
        audit={'species':sp}; normalized=[]; fetched=0; seen=set()
        try:
            tax=resolve(sp,client.get('taxa',{'q':sp,'rank':'species','per_page':30})['results'])
            audit.update(taxon_id=tax['id'],resolved_name=tax['name'],matched_term=tax.get('matched_term'))
            n=set().union(*(np[h] for h in hp[sp]));c=set().union(*(cp[h] for h in hp[sp]))
            total=None
            for page in range(1,args.max_pages+1):
                data=client.get('observations',{'taxon_id':tax['id'],'term_id':1,'term_value_id':6,'quality_grade':'research','captive':'false','geo':'true','photos':'true','d1':START,'d2':END,'per_page':200,'page':page,'order_by':'id','order':'asc'})
                if total is None:total=data['total_results']
                for o in data['results']:
                    if o['id'] in seen:raise ValueError('Duplicate observation across pages')
                    seen.add(o['id']);normalized.append(summarize_record(o,sp,tax['id'],mapper,n,c))
                fetched+=len(data['results'])
                if fetched>=total or len(data['results'])<200:break
            audit.update(provider_total=total,fetched=fetched,complete=(fetched>=total),status='COMPLETE' if fetched>=total else 'CAPPED_OLDEST_1000_NOT_RANDOM')
        except Exception as ex:
            audit.update(status='SOURCE_OR_RESOLUTION_ERROR',error=str(ex),fetched=fetched,complete=False)
        results.extend(normalized);catalog.append(audit)
        (out/'source_requests.json').write_text(json.dumps(client.calls,indent=2)+'\n')
        (out/'species_retrieval_audit.json').write_text(json.dumps(catalog,indent=2)+'\n')
        print(json.dumps(audit),flush=True)
    with gzip.open(out/'larval_observations_normalized.json.gz','wt',encoding='utf-8') as f:json.dump(results,f)
    valid=[r for r in results if 'excluded_reason' not in r]
    target=[r for r in valid if r['resource_class']=='introduced_only_in_frozen_hosts']
    precise=[r for r in target if r['coordinate_precision']=='known_le_10000m']
    summary={'status':'POSTHOC_BOUNDED_LARVAL_PRESENCE_SCREEN_NOT_HOST_USE_VALIDATION','panel_species':32,'date_window':[START,END],
             'complete_species':sum(x['complete'] for x in catalog),'capped_species':sum(x['status'].startswith('CAPPED') for x in catalog),
             'source_error_species':sum(x['status']=='SOURCE_OR_RESOLUTION_ERROR' for x in catalog),'downloaded_observations':len(results),
             'eligible_larval_observations':len(valid),'excluded':dict(Counter(r['excluded_reason'] for r in results if 'excluded_reason' in r)),
             'resource_classes':dict(Counter(r['resource_class'] for r in valid)),'target_larval_records':len(target),
             'target_species':sorted({r['focal_species'] for r in target}),'target_species_regions':len({(r['focal_species'],r['region']) for r in target}),
             'target_known_accuracy_le_10km_records':len(precise),'target_known_accuracy_le_10km_species_regions':len({(r['focal_species'],r['region']) for r in precise}),
             'limitation':'Descriptive first-stage screen. Larval annotation and point mapping do not identify the feeding host, completion of development or established populations. No new-area absence or detection probability is inferred. Oldest-first caps are not random samples. Minimum taxonomy consistency does not imply complete host records. Positional accuracy is not a boundary-buffer proof.'}
    (out/'larval_presence_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'introduced_only_larval_candidates.json').write_text(json.dumps(target,indent=2)+'\n')
    print('SUMMARY '+json.dumps(summary),flush=True)
    if summary['source_error_species']==32:raise RuntimeError('All source queries failed')

if __name__=='__main__': main()
