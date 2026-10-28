#!/usr/bin/env python3
"""UKBMS GBIF historical dataset source fields, no raw occurrence row export."""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
from hashlib import sha256
import json,time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen

DATASET='1e266c3d-92ef-4d5a-8e4a-c04742c772c3'
BASE='https://api.gbif.org/v1/'
YEARS={'1970_1990':'1970,1990','2010_2024':'2010,2024'}
FIELDS=('eventID','samplingProtocol','locationID','individualCount',
        'occurrenceStatus','eventDate','basisOfRecord')
SCHEMA='chocho_UKBMS_historical_event_and_absence_GBIF_source_audit_v01'


def validate_year_facets(d):
    if d.get('limit')!=0 or d.get('results')!=[] or type(d.get('count')) is not int:
        raise ValueError('Source should be limit0 year facets')
    facets=d.get('facets')
    if not isinstance(facets,list):raise ValueError('Missing YEAR facet')
    if not facets and d['count']==0:return {'records':0,'years':{}}
    if len(facets)!=1 or facets[0].get('field')!='YEAR':
        raise ValueError('YEAR facet absent or malformed')
    items=facets[0].get('counts')
    if not isinstance(items,list) or len(items)>=100:
        raise ValueError('YEAR facet incomplete')
    years={}
    for item in items:
        y=item.get('name')
        n=item.get('count')
        if not isinstance(y,str) or not y.isdigit() or not isinstance(n,int) or n<0 or y in years:
            raise ValueError('YEAR source facet invalid')
        years[y]=n
    if sum(years.values())>d['count']:
        raise ValueError('YEAR sum exceeds source total')
    return {'records':d['count'],'years':dict(sorted(years.items()))}


def field_presence(d):
    if not isinstance(d,dict) or not isinstance(d.get('results'),list):
        raise ValueError('Invalid source occurrence schema')
    rows=d['results']
    if len(rows)>12:raise ValueError('Sample exceeds expected bound')
    return {'returned_sample_rows':len(rows),
            'number_with_field_nonblank':{f:sum(x.get(f) not in (None,'') for x in rows) for f in FIELDS},
            'recording_source_samples_are_not_a_site_visit_design':True}


def fetch(path,queries):
    url=BASE+path+'?'+urlencode(queries) if queries else BASE+path
    req=Request(url,headers={'User-Agent':'chocho-ukbms-historical-field-audit/0.1','Accept':'application/json'})
    with urlopen(req,timeout=45) as r:
        if r.status!=200:raise ValueError('GBIF source HTTP non200')
        raw=r.read(900_001)
    if len(raw)>900_000:raise ValueError('GBIF source record body exceeded bounded query')
    parsed=json.loads(raw)
    if not isinstance(parsed,dict):raise ValueError('Unexpected GBIF response')
    return parsed,{'query_url':url,'source_response_sha256':sha256(raw).hexdigest(),
                   'source_response_bytes':len(raw)}


def main(path):
    dataset,receipt=fetch('dataset/'+DATASET,[])
    if dataset.get('key')!=DATASET:
        raise ValueError('Fixed UKBMS GBIF dataset key does not match original registration')
    result={'schema':SCHEMA,'status':'SOURCE_TRANSPORT_BEGIN',
            'dataset_registration':DATASET,
            'dataset_metadata':{'title':dataset.get('title'),
                                'license':dataset.get('license'),
                                'publishingOrganizationKey':dataset.get('publishingOrganizationKey'),
                                'source_receipt':receipt},
            'original_occurrence_identifiers_locations_photographs_exported':False,
            'windows':{}}
    for label,window in YEARS.items():
        base={'datasetKey':DATASET,'country':'GB','hasCoordinate':'true','year':window}
        data,cr=fetch('occurrence/search',list(base.items())+[
            ('occurrenceStatus','PRESENT'),('limit','0'),('facet','year'),('facetLimit','100')])
        present=validate_year_facets(data)
        missing_data,ar=fetch('occurrence/search',list(base.items())+[
            ('occurrenceStatus','ABSENT'),('limit','0')])
        if missing_data.get('limit')!=0 or type(missing_data.get('count')) is not int:
            raise ValueError('UKBMS source absence count invalid')
        sample, sr=fetch('occurrence/search',list(base.items())+[
            ('occurrenceStatus','PRESENT'),('limit','12')])
        sample_data=field_presence(sample)
        result['windows'][label]={
            'presence_records_in_GBIF_index':present['records'],
            'year_index_counts':present['years'],
            'explicit_GBIF_ABSENT_occurrence_records':missing_data['count'],
            'source_field_availability_from_fixed_first_12_rows':sample_data,
            'source_checksums':{'presence':cr,'absent':ar,'field_sample':sr},
            'historical_zero_detections_identified':False}
        print(json.dumps({'window':label,'indexed_present':present['records'],
                          'indexed_ABSENT':missing_data['count'],
                          'sample_fields':sample_data['number_with_field_nonblank']},sort_keys=True),flush=True)
        time.sleep(.4)
    result['status']='UKBMS_GBIF_DATASET_HISTORICAL_SOURCE_AUDITED'
    result['inference']={
        'explicit_same_site_survey_visits_with_non_detections_verified':False,
        'GBIF_absent_occurrence_status_is_not_complete_survey_non_detections':True,
        'site_year_occupancy_inferred':False,
        'host_climate_hindcast_estimated':False,
        'raw_individual_visitor_or_site_records_saved':False,
        'note':'Original source schema only. UKBMS monitored transect visit records may require a separate licence.'}
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    main(a.output)
