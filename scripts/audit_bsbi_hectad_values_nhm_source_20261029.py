#!/usr/bin/env python3
"""Original BSBI 10km era coding, BNG LTMN overlap, and NHM HOSTS source gate.

Source-locked metadata/aggregate only. Original site/grid IDs, botanical
source presences, recorder/observer fields and larval host rows are not saved.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,re,statistics,zipfile
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urlencode,urlparse
from urllib.error import HTTPError,URLError
from audit_bsbi_atlas_and_clarke_source_20261029 import fetch
from audit_bsbi_nested_era_host_schema_20261029 import OUTER_MD5,OUTER_SHA,EXPECTED_NESTED

HECTAD_COLS=('hectad','to 1929','1930 - 1969','1970 - 1986','1987 - 1999','2000 - 2009','2010 - 2019','hectadstatus','distinct tetrads')
NHM_RESOURCE='877f387a-36a3-486c-a0c1-b8d5fb69f85a'


def grid10(e,n):
    """Official British National Grid two-letter 100km+single-digit 10km square."""
    if not isinstance(e,(float,int)) or not isinstance(n,(float,int)) or not(0<=e<700000 and 0<=n<1300000):
        return None
    e100,n100=int(e//100000),int(n//100000)
    l1=(19-n100)-((19-n100)%5)+((e100+10)//5)
    l2=((19-n100)*5%25)+(e100%5)
    if l1>7:l1+=1
    if l2>7:l2+=1
    if not(0<=l1<26 and 0<=l2<26):return None
    return chr(l1+65)+chr(l2+65)+str(int(e//10000)%10)+str(int(n//10000)%10)


def nested_distributions_scan(original:bytes, original_uk_public_grid_codes:set[str]|None=None):
    if hashlib.md5(original).hexdigest()!=OUTER_MD5 or hashlib.sha256(original).hexdigest()!=OUTER_SHA:
        raise ValueError('BSBI original MD5/SHA256 mismatch')
    with zipfile.ZipFile(io.BytesIO(original)) as outer:
        raw=outer.read(EXPECTED_NESTED)
        lookup=outer.read('nameLookup.csv')
    meta={'nested_sha256':hashlib.sha256(raw).hexdigest(),'nameLookup_sha256':hashlib.sha256(lookup).hexdigest()}
    lo=list(csv.DictReader(io.StringIO(lookup.decode('utf-8-sig',errors='replace'))))
    lookup_headers=list(lo[0]) if lo else []
    names=[str(x.get('ddbidLookup') or '') for x in lo]
    all_names=set(names)
    output=collections.Counter()
    codes={col:collections.Counter() for col in HECTAD_COLS[1:7]}
    types=collections.Counter();matched=collections.Counter()
    nfiles=0;unexpected=0;file_keys=[]
    with zipfile.ZipFile(io.BytesIO(raw)) as inner:
        for info in inner.infolist():
            if info.is_dir() or not info.filename.lower().endswith('.csv'):continue
            nfiles+=1
            name=Path(info.filename).stem
            file_keys.append(name)
            with inner.open(info) as bf:
                data=csv.DictReader(io.TextIOWrapper(bf,encoding='utf-8-sig',errors='replace',newline=''))
                if tuple(data.fieldnames or ())!=HECTAD_COLS:
                    raise ValueError('Published per-taxon hectad era column schema drift')
                # Only source *coding vocabulary* and counts, not individual hectads exported.
                for row in data:
                    output['source_species_by_hectad_rows']+=1
                    hectad=str(row.get('hectad') or '').strip().upper()
                    if re.fullmatch(r'[A-Z]{2}[0-9]{2}',hectad):types['UK_BNG_hectad_code']+=1
                    else:types['other_hectad_code_format']+=1
                    if original_uk_public_grid_codes is not None and hectad in original_uk_public_grid_codes:
                        matched['source_taxon_era_rows_in_original_uk_transect_10km_grids']+=1
                    for col in codes:
                        v=str(row.get(col) or '').strip().upper()
                        codes[col][v or '[empty]']+=1
    lookup_match=sum(file in all_names for file in file_keys)
    return {**meta,'source_plant_names_lookup_rows':len(lo),
            'source_lookup_headers':lookup_headers,
            'source_nested_plant_csv_files':nfiles,
            'file_stems_matching_ddbidLookup':lookup_match,
            'file_stems_not_matching_ddbidLookup':nfiles-lookup_match,
            'hectad_source_rows_total':output['source_species_by_hectad_rows'],
            'hectad_format_counts':dict(types),
            'original_source_period_value_vocabulary':{k:dict(v.most_common(18)) for k,v in codes.items()},
            'original_UKBMS_transect_hectad_count':len(original_uk_public_grid_codes or []),
            'uk_transect_grid_intersections_source_taxon_rows':matched['source_taxon_era_rows_in_original_uk_transect_10km_grids'],
            'no_site_level_record_or_original_botanical_taxon_observation_exported':True}


def original_geocoded_19_site_grid_codes():
    # Use prior original Natural England SHA-pinned crosswalk/UKBMS verified
    # survey gate; no site or grid identifiers serialized.
    from audit_uk_ltmn_site_geography_20261028 import read_crosswalk
    from audit_uk_ltmn_butterfly_detection_panel_20261028 import site_effort
    from audit_uk_ltmn_historical_event_join_20261028 import get_source,excel_rows,FIELDS_VISITS
    cross,qc=read_crosswalk()
    v,sites,period=site_effort(excel_rows(get_source('visits'),FIELDS_VISITS),cross)
    if len(sites)!=19:raise ValueError('Original immutable UKBMS 19 source transect panel drift')
    grid={grid10(cross[site][0],cross[site][1]) for site in sites}
    if None in grid:raise ValueError('Invalid original British National Grid projection')
    return grid,{'original_verified_transects':len(sites),'original_unique_10km_hectads':len(grid),
                 'UKBMS_source_site_coordinates_serialized':False}


def nhm_original_hosts_metadata():
    result={'original_nhm_hosts_resource_key':NHM_RESOURCE,
            'source_page':'https://data.nhm.ac.uk/dataset/hosts/resource/'+NHM_RESOURCE,
            'copyright':'NHM HOSTS CCZero Robinson et al. 2023 (DOI:10.5519/havt50xw)',
            'raw_larval_relation_rows_or_locations_exported':False}
    try:
        api='https://data.nhm.ac.uk/api/3/action/resource_show?'+urlencode({'id':NHM_RESOURCE})
        with urlopen(Request(api,headers={'User-Agent':'chocho-UK-BSBI-hostatlas/1.0'}),timeout=45) as r:
            payload=r.read(1_500_001)
        if len(payload)>1_500_000:raise ValueError('Original NHM metadata too large')
        j=json.loads(payload)
        info=j.get('result',{})
        if not j.get('success') or info.get('id')!=NHM_RESOURCE:
            raise ValueError('Original NHM HOSTS resource identity changed')
        dataurl=info.get('url','')
        if not dataurl.startswith(('https://data.nhm.ac.uk/','https://media.nhm.ac.uk/')):
            raise ValueError('NHM HOSTS CSV redirected to unapproved external publisher')
        with urlopen(Request(dataurl,headers={'User-Agent':'chocho-UK-BSBI-hostatlas/1.0'}),timeout=80) as r:
            data=r.read(25_000_001)
        if not data or len(data)>25_000_000:raise ValueError('NHM original HOSTS CSV exceeded expected size')
        reader=csv.reader(io.StringIO(data.decode('utf-8-sig',errors='replace')))
        header=next(reader,[])
        rows=sum(1 for _ in reader)
        result.update({'status':'ORIGINAL_NHM_HOSTS_CSV_AVAILABLE',
                       'source_format':str(info.get('format'))[:40],
                       'source_bytes':len(data),
                       'source_sha256':hashlib.sha256(data).hexdigest(),
                       'original_header_fields':[str(a).strip()[:90] for a in header],
                       'original_relation_rows':rows})
    except (OSError,ValueError,TypeError,KeyError,HTTPError,URLError) as e:
        result.update({'status':'ORIGINAL_NHM_HOSTS_UNAVAILABLE',
                       'error_type':type(e).__name__,
                       'error_reason':str(e)[:250]})
    return result


def run(output):
    grid,source=original_geocoded_19_site_grid_codes()
    archive=fetch('https://zenodo.org/records/10998003/files/PA2020_10by10kmSquareData.zip?download=1',25_000_000)
    atlas=nested_distributions_scan(archive,grid)
    nhm=nhm_original_hosts_metadata()
    result={'schema':'chocho_BSBI_hectad_period_codes_original_UKBMS_grid_and_NHM_hosts_source_20261029_v01',
            'status':'BSBI_3497_PLANT_10KM_SOURCE_CODES_AND_UKBMS_GRID_VERIFIED',
            'original_BSBI_source':atlas,
            'original_UKBMS_1976_2021_transect_source':source,
            'NHM_host_caterpillar_original_source':nhm,
            'no_original_plant_hectad_species_status_or_butterfly_site_ids_exported':True}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+'\n')
    print(json.dumps({'status':result['status'],
                     'original_plant_files':atlas['source_nested_plant_csv_files'],
                     'lookup_file_match':atlas['file_stems_matching_ddbidLookup'],
                     'hectad_rows':atlas['hectad_source_rows_total'],
                     'plant_status_code_vocab':atlas['original_source_period_value_vocabulary'],
                     'n_uk_transects':source['original_verified_transects'],
                     'n_original_uk_10km_grids':source['original_unique_10km_hectads'],
                     'matched_taxon_grid_rows':atlas['uk_transect_grid_intersections_source_taxon_rows'],
                     'nhm_source_status':nhm['status'],'NHM_headers':nhm.get('original_header_fields')},ensure_ascii=False),flush=True)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
