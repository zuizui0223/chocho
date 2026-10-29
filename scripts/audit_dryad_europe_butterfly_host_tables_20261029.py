#!/usr/bin/env python3
"""Inspect public original Clarke European butterfly larval host CSV files.

No host relationship rows released; only provider original CSV hashes, headers
and source number of records. A blocked file endpoint is NOT zero host use.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from urllib.parse import quote

DOI='doi:10.5061/dryad.1vhhmgr12'
VERSION=274312
FILES={
'larval_host_plant.csv':688783,
'larval_plants.csv':2381966,
'larval_host_plant_status.csv':679,
'butterfly_synonym.csv':807480,
'README.md':18755,
}


def get(url,limit):
    if not url.startswith('https://datadryad.org/'):
        raise ValueError('Expected public Dryad original source URL')
    req=Request(url,headers={'User-Agent':'chocho-source-verified-europe-butterfly-hosts/1.0',
                             'Accept':'text/csv,application/json,*/*'})
    with urlopen(req,timeout=75) as response:
        data=response.read(limit+1)
        status=response.status
        final_host=response.geturl().split('/')[2].split(':')[0]
    if status!=200 or not data or len(data)>limit:
        raise ValueError('Dryad source download oversized or invalid')
    return data,final_host


def discover_file_records(data):
    embedded=data.get('_embedded') or {}
    fs=embedded.get('stash:files')
    if not isinstance(fs,list):
        raise ValueError('No official Dryad published file catalogue')
    candidates={}
    for row in fs:
        name=str(row.get('path') or row.get('filename') or '').split('/')[-1]
        if name in FILES:
            if name in candidates:raise ValueError('Duplicated exact published Dryad file '+name)
            candidates[name]=row
    if len(candidates)!=len(FILES):
        raise ValueError('Original expected Dryad file set incomplete')
    return candidates


def csv_metadata(raw,expected_name):
    if not raw or raw[:1]==b'<' or b'<html' in raw[:500].lower():
        raise ValueError('Dryad returned HTML not original CSV')
    if b'\x00' in raw[:300]:raise ValueError('Dryad returned binary not original tabular source')
    text=raw.decode('utf-8-sig')
    reader=csv.reader(io.StringIO(text,newline=''))
    header=next(reader,[])
    if len(header)<2:raise ValueError('Dryad original file lacks tabular columns')
    rows=sum(1 for _ in reader)
    if rows==0:raise ValueError('Original host file contains no source records')
    return {'file':expected_name,'header':[str(x).strip()[:120] for x in header][:120],
            'original_rows':rows,'source_sha256':hashlib.sha256(raw).hexdigest(),
            'source_bytes':len(raw),'original_data_rows_not_exported':True}


def inspect(path):
    result={'schema':'chocho_Dryad_original_Clarke_2024_Europe_butterfly_larval_host_source_20261029_v01',
            'source_doi':DOI,'publisher':'Dryad original published dataset',
            'original_version_id':VERSION,'original_host_plant_rows_exported':False}
    catalogue_url='https://datadryad.org/api/v2/versions/'+str(VERSION)+'/files'
    try:
        payload,_=get(catalogue_url,1_500_000)
        catalog=discover_file_records(json.loads(payload))
        entries=[]
        for name,size in FILES.items():
            row=catalog[name]
            fid=row.get('id')
            if not isinstance(fid,int) or fid<=0:
                raise ValueError('Published Dryad file id absent')
            md={'original_file':name,'original_file_id':fid,'published_size_bytes':row.get('size'),
                'source_download_status':'NOT_ATTEMPTED'}
            try:
                endpoint='https://datadryad.org/api/v2/files/'+str(fid)+'/download'
                b,host=get(endpoint,min(size*3+20000,6_000_000))
                if name.lower().endswith('.csv'):
                    md.update(csv_metadata(b,name))
                else:
                    if not b'README' in b[:1000].upper() and len(b)<100:
                        raise ValueError('Original README not valid')
                    md.update({'source_sha256':hashlib.sha256(b).hexdigest(), 'source_bytes':len(b)})
                md['source_download_status']='ORIGINAL_DRYAD_FILE_VERIFIED'
                md['final_source_data_host']=host
            except (OSError,ValueError,UnicodeError,HTTPError,URLError) as e:
                md['source_download_status']='ORIGINAL_FILE_ACCESS_BLOCKED'
                md['blocker_type']=type(e).__name__
                md['blocker_reason']=str(e)[:200]
            entries.append(md)
        result['files']=entries
        result['available_original_host_tables']=sum(x['source_download_status']=='ORIGINAL_DRYAD_FILE_VERIFIED' and x['original_file'].endswith('.csv') for x in entries)
        result['status']='ORIGINAL_DRYAD_HOST_CATALOGUE_WITH_VERIFIED_CSVS' if result['available_original_host_tables']>=2 else 'ORIGINAL_DRYAD_CATALOGUE_ONLY_OR_BLOCKED_DOWNLOADS'
    except (OSError,ValueError,TypeError,HTTPError,URLError) as e:
        result['status']='DRYAD_ORIGINAL_FILE_CATALOGUE_UNAVAILABLE'
        result['failure_type']=type(e).__name__
        result['failure_reason']=str(e)[:200]
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'files':[{k:v.get(k) for k in ('original_file','source_download_status','original_rows','header','blocker_reason')} for v in result.get('files',[])]},ensure_ascii=False),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    r=inspect(p.parse_args().output)
    if r['status']!='ORIGINAL_DRYAD_HOST_CATALOGUE_WITH_VERIFIED_CSVS':raise SystemExit(2)
