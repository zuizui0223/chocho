#!/usr/bin/env python3
"""Inspect the *inner* published BSBI hectad distribution ZIP and Dryad files.

Only source schema and field-name metadata are released; never raw botanical
10km presences, taxon-level outcomes, site coordinates or private observations.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,zipfile
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen

from audit_bsbi_atlas_and_clarke_source_20261029 import fetch,SOURCES

OUTER_MD5='66a8bb4305df484839587b5fe8acf1d2'
OUTER_SHA='fdd9db2166ba93fda66d24af9ba5b561f6aa4108a9d71c28f0e453c80a9b9ddc'
EXPECTED_NESTED='hectadDistributions.zip'
MAX_INNER_UNCOMPRESSED=400_000_000


def nested_schema(archive:bytes):
    if hashlib.md5(archive).hexdigest()!=OUTER_MD5 or hashlib.sha256(archive).hexdigest()!=OUTER_SHA:
        raise ValueError('Original outer BSBI atlas archive source drift')
    with zipfile.ZipFile(io.BytesIO(archive)) as outer:
        nested=[x for x in outer.namelist() if Path(x).name==EXPECTED_NESTED]
        if len(nested)!=1:raise ValueError('Original published 10km nested archive missing')
        raw=outer.read(nested[0])
        inner_sha=hashlib.sha256(raw).hexdigest()
    if raw[:4]!=b'PK\x03\x04':raise ValueError('Original nested distribution is not ZIP')
    with zipfile.ZipFile(io.BytesIO(raw)) as inner:
        if inner.testzip() is not None:raise ValueError('Nested source ZIP checksum error')
        infos=[x for x in inner.infolist() if not x.is_dir()]
        if sum(x.file_size for x in infos)>MAX_INNER_UNCOMPRESSED:
            raise ValueError('Original published hectad data exceeds source-safe bound')
        out=[]
        for i,f in enumerate(infos):
            p=Path(f.filename)
            if p.is_absolute() or '..' in p.parts:raise ValueError('Nested ZIP unsafe source path')
            meta={'entry_number':i+1,'filename':p.name[:180],
                  'suffix':p.suffix.lower(),'unpacked_bytes':f.file_size}
            if p.suffix.lower() in ('.csv','.tsv'):
                with inner.open(f) as z:
                    reader=csv.reader(io.TextIOWrapper(z,encoding='utf-8-sig',errors='replace',newline=''),delimiter='\t' if p.suffix.lower()=='.tsv' else ',')
                    meta['header_fields']=[str(x).strip()[:100] for x in next(reader,[])][:100]
                    # Do not inspect actual plant-by-hectad survey outcomes.
                    meta['original_rows_read']=0
            out.append(meta)
    return {'original_nested_zip_sha256':inner_sha,'nested_zip_bytes':len(raw),
            'distribution_file_count':len(out),'distribution_file_schemas':out[:100]}


def dryad_file_metadata():
    r={'version_id':274312,'source':'https://doi.org/10.5061/dryad.1vhhmgr12',
       'original_source_biology_content_read':False}
    try:
        uri='https://datadryad.org/api/v2/versions/274312/files'
        body=fetch(uri,1_500_000)
        j=json.loads(body)
        embedded=j.get('_embedded',{})
        raw=embedded.get('stash:files',[])
        if not isinstance(raw,list):raise ValueError('Dryad files metadata missing')
        r['public_file_rows_returned']=len(raw)
        r['public_download_names_sizes_only']=[{'file_name':str(x.get('path') or x.get('filename') or '').split('/')[-1][:150],
                                                 'size_bytes':x.get('size')} for x in raw][:100]
        wanted={'larval_host_plant.csv','larval_plants.csv','larval_host_plant_status.csv'}
        r['required_larval_source_file_names_present']=sorted(wanted & {x['file_name'] for x in r['public_download_names_sizes_only']})
        r['status']='DRYAD_ORIGINAL_FILE_CATALOGUE_VERIFIED' if r['required_larval_source_file_names_present'] else 'DRYAD_NO_EXPECTED_FILE_IN_RETURNED_PAGE'
    except (OSError,ValueError,TypeError,HTTPError,URLError) as e:
        r['status']='DRYAD_FILE_METADATA_UNAVAILABLE'
        r['error_class']=type(e).__name__
        r['failure_reason']=str(e)[:240]
    return r


def run(output):
    result={'schema':'chocho_BSBI_10km_historical_dateclass_inner_schema_20261029_v01',
            'source_atlas_zenodo_record':10998003,
            'original_site_record_or_botanical_presence_exported':False}
    try:
        uri='https://zenodo.org/records/10998003/files/PA2020_10by10kmSquareData.zip?download=1'
        original=fetch(uri,SOURCES[0][-1])
        result.update(nested_schema(original))
        result['status']='ORIGINAL_BSBI_INNER_DISTRIBUTION_SCHEMA_VERIFIED'
    except (OSError,ValueError,zipfile.BadZipFile,HTTPError,URLError) as e:
        result['status']='ORIGINAL_BSBI_INNER_SOURCE_FAILED'
        result['failure_class']=type(e).__name__
        result['failure_reason']=str(e)[:230]
    result['dryad_files']=dryad_file_metadata()
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'nested':result.get('distribution_file_schemas'),
                       'dryad':result['dryad_files']['status']},sort_keys=True))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    x=run(p.parse_args().output)
    if x['status']!='ORIGINAL_BSBI_INNER_DISTRIBUTION_SCHEMA_VERIFIED':raise SystemExit(2)
