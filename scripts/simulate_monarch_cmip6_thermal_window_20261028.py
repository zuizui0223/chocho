#!/usr/bin/env python3
"""Physical CMIP6 thermal-window scenarios for public monarch host-record sites.

NOT a forecast of butterfly counts, colonization, winter breeding or plant survival.
Public larval records identify source sampling locations. Results are hex-equal
thermal opportunities conditional on hypothesized host leaf persistence, with
observer and exact locations never serialized.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import numpy as np
from rasterio.io import MemoryFile
from audit_monarch_annotated_host_source_20261028 import source_download, SHARE_URL
from audit_monarch_host_label_categories_20261028 import ORIGINAL_ZIP_SHA256, CORE, THINNED
from analyze_monarch_cultivation_season_20261028 import canonical_host

SCHEMA = 'chocho_monarch_cmip6_recorded_locality_thermal_opportunity_v01'
BASE_URL = 'https://geodata.ucdavis.edu/climate/worldclim/2_1/base/'
FUT_URL = 'https://geodata.ucdavis.edu/climate/worldclim/2_1/fut/10m/'
MODELS = ('MIROC6', 'IPSL-CM6A-LR')
SSPS = ('ssp126', 'ssp585')
PERIODS = ('2041-2060', '2081-2100')
VARS = ('tmin', 'tmax')
WINTER = (11, 12, 1, 2)
SUMMER = (6, 7, 8)
DAYS = {1:31,2:28,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
LOWER_C = 11.5
UPPER_C = 33.0
ONE_GENERATION_GDD = 323.0
HEAT_SCREEN_TMAX_C = 38.0
QUALITY_COORD_LIMIT_METERS = 20000.0
MIN_QUALIFIED_SOURCE_HEX = 30


def archive_name(var, model=None, ssp=None, period=None):
    if var not in VARS: raise ValueError('Unsupported climate variable')
    if model is None:
        if ssp is not None or period is not None: raise ValueError('Inconsistent baseline spec')
        return f'wc2.1_10m_{var}.zip'
    if model not in MODELS or ssp not in SSPS or period not in PERIODS:
        raise ValueError('Climate model/SSP/period not in frozen scenario set')
    return f'wc2.1_10m_{var}_{model}_{ssp}_{period}.zip'


def sha256_file(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def stable_core_coords(archive):
    coords={}
    conflicting=set()
    qc=Counter()
    with archive.open(CORE) as stream:
        reader=csv.DictReader(io.TextIOWrapper(stream,encoding='utf-8-sig',newline=''))
        required={'uri','latitude','longitude','positional_accuracy'}
        if not required<=set(reader.fieldnames or []):raise ValueError('Original core coordinates schema drift')
        for row in reader:
            qc['core_rows']+=1
            key=str(row.get('uri') or '').strip()
            if not key:continue
            try:
                lat=float(row.get('latitude') or 'nan');lon=float(row.get('longitude') or 'nan')
                acc=float(row.get('positional_accuracy') or 'nan')
            except (ValueError,TypeError):continue
            if not(-90<lat<90 and -180<lon<180):
                qc['unscorable_lat_lon']+=1;continue
            if math.isfinite(acc) and (acc < 0 or acc > QUALITY_COORD_LIMIT_METERS):
                qc['source_uncertainty_exceeds_20km']+=1;continue
            value=(lat,lon)
            if key in coords and coords[key]!=value:conflicting.add(key)
            coords[key]=value
    for key in conflicting:coords.pop(key,None)
    qc['ambiguous_location_keys']=len(conflicting)
    qc['georeferenced_core_uris']=len(coords)
    if qc['core_rows']!=30320:
        raise ValueError('Original published source core 30320 rows changed')
    return coords,dict(qc)


def source_hex_centroids(archive, coords):
    hexes=defaultdict(list)
    groups=defaultdict(set)
    qc=Counter()
    thinned_uri_seen=set()
    with archive.open(THINNED) as stream:
        reader=csv.DictReader(io.TextIOWrapper(stream,encoding='utf-8-sig',newline=''))
        required={'uri','hex_id','host','rearing','setting','confidence','EW','native_status','cultivated','month'}
        if not required<=set(reader.fieldnames or []):raise ValueError('Original source thinned schema drift')
        for row in reader:
            qc['source_rows']+=1
            key=str(row.get('uri') or '').strip()
            h=str(row.get('hex_id') or '').strip()
            host=canonical_host(str(row.get('host') or ''))
            if not key or not h or not host:
                qc['missing_uri_hex_or_species_host']+=1;continue
            if key in thinned_uri_seen:
                qc['duplicate_original_observation_uri']+=1;continue
            thinned_uri_seen.add(key)
            if str(row['rearing']).strip().lower()!='no' or str(row['setting']).strip().lower()=='obscured':
                qc['reared_or_obscured']+=1;continue
            if str(row['confidence']).strip().lower() not in ('certain','high','medium') or str(row['EW']).strip().lower()=='na':
                qc['unknown_host_or_biogeographic_region']+=1;continue
            native=str(row['native_status']).strip().lower()
            cultivated=str(row['cultivated']).strip().lower()
            if native not in ('native','non-native') or cultivated not in ('yes','no','likely','unlikely','maybe'):
                qc['unknown_origin_or_cultivation']+=1;continue
            location=coords.get(key)
            if location is None:
                qc['no_original_core_coordinate_join']+=1;continue
            try: month=int(str(row['month']).strip())
            except ValueError:
                qc['invalid_original_month']+=1;continue
            if not 1<=month<=12:
                qc['invalid_original_month']+=1;continue
            qc['accepted_thinned_source_records']+=1
            hexes[h].append(location)
            groups['all'].add(h)
            groups[f'host_origin_{native}'].add(h)
            groups[f'cultivated_{cultivated}'].add(h)
            if month in WINTER:groups['source_winter_larval_record'].add(h)
    if qc['source_rows']!=12942:
        raise ValueError('Original thinned cohort changed')
    # Approximate one study-defined 100 km locality with the median of source photo
    # coordinates in that cell. No centroid, coordinates or hex ID is exported.
    centers={h:(statistics.median(lat for lat,_ in points),
                statistics.median(lon for _,lon in points)) for h,points in hexes.items()}
    if len(centers)<MIN_QUALIFIED_SOURCE_HEX:
        raise ValueError('Original source-core geo join too sparse for frozen hex gate')
    # Never return observer, URIs, original source rows or raw source coordinates to report.
    return centers, groups, dict(qc)


def grid_cell(lat,lon):
    if not (-90<lat<90 and -180<lon<180):return None
    rr=int((90.-lat)*6.)
    cc=int((lon+180.)*6.)
    if 0<=rr<1080 and 0<=cc<2160:return (rr,cc)
    return None


def month_from_file(filename):
    name=Path(filename).name
    m=re.search(r'_(\d{1,2})\.tif$',name,flags=re.I)
    return int(m.group(1)) if m and 1<=int(m.group(1))<=12 else None


def _validate_temperature_grid(ds):
    if ds.crs.to_epsg()!=4326 or (ds.height,ds.width)!=(1080,2160):
        raise ValueError('Unexpected WorldClim WGS84 10-minute grid')
    tf=ds.transform
    if not(math.isclose(tf.a,1/6,abs_tol=1e-7) and
           math.isclose(-tf.e,1/6,abs_tol=1e-7) and
           math.isclose(tf.c,-180,abs_tol=1e-6) and
           math.isclose(tf.f,90,abs_tol=1e-6)):
        raise ValueError('WorldClim grid origin mismatch')


def _sample_tif_band(ds, band, pairlist):
    _validate_temperature_grid(ds)
    rows=[r for r,c in pairlist];cols=[c for r,c in pairlist]
    data=ds.read(band,masked=True)
    if np.any(np.ma.getmaskarray(data)[rows,cols]):
        raise ValueError('Missing climate raster pixel at fixed source locale')
    raw=np.asarray(data[rows,cols],dtype=float)
    if not np.all(np.isfinite(raw)):
        raise ValueError('Non-finite WorldClim source climate pixel')
    return raw,float(ds.scales[band-1])


def read_sampled_archive(path, pairs, scale=None):
    """WorldClim historical: 12 single-band TIFs; CMIP6: one 12-band TIF.

    Both are valid source representations of 12 months. Fail closed on mixed,
    incomplete or ambiguous archives, without substituting climate variables.
    """
    path=Path(path)
    if path.stat().st_size<100000 or path.stat().st_size>200000000:
        raise ValueError('Unexpected WorldClim archive size')
    pairlist=sorted(set(pairs))
    vals={};scales={}
    with ZipFile(path) as z:
        members=[x for x in z.namelist() if x.lower().endswith('.tif')]
        if len(members)==12:
            names={month_from_file(m):m for m in members}
            if set(names)!=set(range(1,13)):
                raise ValueError('12-file WorldClim monthly archive has ambiguous months')
            for month in range(1,13):
                with MemoryFile(z.read(names[month])) as mem, mem.open() as ds:
                    if ds.count!=1:
                        raise ValueError('Expected single-band monthly original GeoTIFF')
                    vals[month],scales[month]=_sample_tif_band(ds,1,pairlist)
        elif len(members)==1:
            with MemoryFile(z.read(members[0])) as mem, mem.open() as ds:
                if ds.count!=12:
                    raise ValueError(f'Expected 12-band CMIP6 source TIF; found {ds.count} bands')
                for month in range(1,13):
                    vals[month],scales[month]=_sample_tif_band(ds,month,pairlist)
        else:
            raise ValueError(f'WorldClim source invalid TIFF count {len(members)}, expected 12 single files or one 12-band TIF')
    return vals,scales,sha256_file(path),pairlist


def choose_temperature_scale(summer_max_samples,declared):
    """Storage encoding QC (not a temperature/ecological threshold fitted to outcomes)."""
    values=np.asarray(summer_max_samples,dtype=float)
    if len(values)==0 or not np.all(np.isfinite(values)):
        raise ValueError('Empty temperature samples')
    if any(math.isfinite(float(x)) and abs(float(x)-1)>1e-7 for x in declared):
        d={round(float(x),6) for x in declared}
        if len(d)!=1 or next(iter(d)) not in (0.1,0.01,1.):
            raise ValueError('Inconsistent explicit TIF scale encoding')
        return next(iter(d))
    m=float(np.median(values))
    if 8.<=m<=60.:return 1.0
    if 80.<=m<=600.:return .1
    raise ValueError(f'Cannot identify WorldClim temperature storage units from median summer Tmax {m:.3f}')


def validate_climate(minima,maxima):
    if set(minima)!=set(range(1,13)) or set(maxima)!=set(range(1,13)):
        raise ValueError('Incomplete annual climate source')
    for month in range(1,13):
        x=np.asarray(minima[month],dtype=float)
        y=np.asarray(maxima[month],dtype=float)
        if x.shape != y.shape or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
            raise ValueError('Climate paired source missing')
        if np.any(x > y+1e-3) or np.any(x < -90) or np.any(y > 75):
            raise ValueError('Physically implausible monthly source climate')


def scenario_thermal(minima,maxima,development_base=LOWER_C,development_cap=UPPER_C):
    validate_climate(minima,maxima)
    if not (0 < development_base < development_cap <= 40):
        raise ValueError('Invalid biologically fixed degree-day limits')
    gdd_winter=np.zeros_like(minima[1],dtype=float)
    winter_early=np.zeros_like(gdd_winter)
    winter_late=np.zeros_like(gdd_winter)
    for month in WINTER:
        avg=.5*(minima[month]+maxima[month])
        gdd=np.maximum(0.,np.minimum(avg,development_cap)-development_base)*DAYS[month]
        gdd_winter+=gdd
        if month in (11,12):winter_early+=gdd
        else:winter_late+=gdd
    heat_flag=np.zeros_like(gdd_winter,dtype=bool)
    for month in SUMMER:
        heat_flag |= maximas_bool(maxima=maxima[month],threshold=HEAT_SCREEN_TMAX_C)
    return {'winter_gdd':gdd_winter,'winter_novdec_gdd':winter_early,
            'winter_janfeb_gdd':winter_late,'summer_monthly_mean_max_above38':heat_flag}


def maximas_bool(maxima,threshold):return np.asarray(maxima)>threshold


def summarize_group(metrics, indices):
    # Retain one count per original 100-km hex. Distinct hexes on the same
    # coarse 10-minute climate pixel remain correlated physical samples.
    take=np.asarray(sorted(indices),dtype=int)
    if len(take)==0: return {'n_100km_hex_localities':0, 'thermal_fraction_ge_323':None}
    winter=metrics['winter_gdd'][take]
    eligible=winter>=ONE_GENERATION_GDD
    return {'n_100km_hex_localities':int(len(take)),
            'winter_gdd_mean':float(np.mean(winter)),
            'winter_gdd_median':float(np.median(winter)),
            'winter_thermal_323_count':int(np.sum(eligible)),
            'thermal_fraction_ge_323':float(np.mean(eligible)),
            'novdec_only_fraction_ge_323':float(np.mean(metrics['winter_novdec_gdd'][take]>=ONE_GENERATION_GDD)),
            'janfeb_only_fraction_ge_323':float(np.mean(metrics['winter_janfeb_gdd'][take]>=ONE_GENERATION_GDD)),
            'winter_required_leaf_fraction_median_among_thermally_feasible':(
                float(np.median(ONE_GENERATION_GDD/winter[eligible])) if np.any(eligible) else None),
            'winter_available_50pct_323_fraction':float(np.mean(winter*.5>=ONE_GENERATION_GDD)),
            'summer_mean_tmax_over38_screen_fraction':float(np.mean(metrics['summer_monthly_mean_max_above38'][take])),
           }


def compute(args):
    # Climate archive verification and source join always use exact original files.
    with TemporaryDirectory(prefix='chocho-monarch-cmip6-') as tmp:
        pub=Path(tmp)/'author-original-public-source.zip'
        source_download(SHARE_URL,pub)
        pubsha=sha256_file(pub)
        if pubsha!=ORIGINAL_ZIP_SHA256:
            raise ValueError('Author public monarch archive source SHA changed')
        with ZipFile(pub) as z:
            coords,core_qc=stable_core_coords(z)
            centers,groups,qc=source_hex_centroids(z,coords)
    pixel_for_hex={h:grid_cell(*center) for h,center in centers.items()}
    if any(v is None for v in pixel_for_hex.values()):
        raise ValueError('Source hex centroid outside WGS84')
    unique_pixels=sorted(set(pixel_for_hex.values()))
    pixel_id={v:i for i,v in enumerate(unique_pixels)}
    group_pixels={g:[pixel_id[pixel_for_hex[h]] for h in hexes]
                  for g,hexes in groups.items()}
    # Reusing exact 100km hex key as unit avoids photo-record pseudoreplication.
    # Multiple hexes can still share one WorldClim pixel; both are reported as
    # non-independent climate support in output.
    report={'schema':SCHEMA,'status':'SOURCE_LOCKED_CMIP6_THERMAL_SCENARIO_COMPLETE',
            'source':{'monarch_author_zip_sha256':pubsha,
                      'climate_base_url':BASE_URL,'climate_cmip6_url':FUT_URL,
                      'historical_period':'1970-2000','source_core_qc':core_qc,
                      'source_thinned_qc':qc,
                      'author_spatially_thinned_hexes':len(centers),
                      'distinct_sampled_worldclim_pixels':len(unique_pixels),
                      'georeferenced_source_coordinates_serialized':False},
            'model':{'gdd_base_c':LOWER_C,'gdd_cap_c':UPPER_C,
                     'thermal_egg_to_adult_target_degree_days':ONE_GENERATION_GDD,
                     'screen_monthly_mean_tmax_c':HEAT_SCREEN_TMAX_C,
                     'season_months_winter':[11,12,1,2],
                     'assumed_month_lengths':DAYS,
                     'seasonal_capacity_is_monthly_climate_proxy_not_daily_development':True,
                     'leaf_persistence_assumptions':{
                       'none':0.0,'half_of_thermal_opportunity':.5,'all_winter':1.0,
                       'novdec_only_months':[11,12],'janfeb_only_months':[1,2]},
                     'no_host_occupancy_or_future_monarch_distribution_modeled':True},
            'archives':{},'scenarios':[]}
    # Derive one fixed WorldClim stored-unit scale from BASE Tmax July.
    base_raw={}
    for var in VARS:
        path=args.climate_dir/archive_name(var)
        raw,tag,digest,pairs=read_sampled_archive(path,unique_pixels)
        if pairs!=unique_pixels:raise ValueError('Misaligned climate sampling grid')
        base_raw[var]=raw
        report['archives'][path.name]={'sha256':digest,'bytes':path.stat().st_size}
        if var=='tmax':
            fixed_scale=choose_temperature_scale(raw[7],tag.values())
    report['source']['temperature_encoding_factor_from_original_july_tmax']=fixed_scale
    base_min={m:base_raw['tmin'][m]*fixed_scale for m in range(1,13)}
    base_max={m:base_raw['tmax'][m]*fixed_scale for m in range(1,13)}
    baseline=scenario_thermal(base_min,base_max)
    base_summary={g:summarize_group(baseline,indices) for g,indices in sorted(group_pixels.items())}
    if base_summary['all']['n_100km_hex_localities']!=len(centers):
        raise ValueError('Original 100-km hex source sampling unit lost in raster join')
    report['baseline']=base_summary
    for model in MODELS:
        for ssp in SSPS:
            for period in PERIODS:
                scenario={}
                for var in VARS:
                    path=args.climate_dir/archive_name(var,model,ssp,period)
                    raw,tag,digest,pairs=read_sampled_archive(path,unique_pixels)
                    if pairs!=unique_pixels:raise ValueError('Future climate grid differs from original')
                    scale=choose_temperature_scale(raw[7] if var=='tmax' else base_raw['tmax'][7],tag.values())
                    if not math.isclose(scale,fixed_scale,abs_tol=1e-8):
                        raise ValueError('Historical/future temperature encoding mismatch')
                    scenario[var]={m:raw[m]*scale for m in range(1,13)}
                    report['archives'][path.name]={'sha256':digest,'bytes':path.stat().st_size}
                future=scenario_thermal(scenario['tmin'],scenario['tmax'])
                bygroup={g:summarize_group(future,indices) for g,indices in sorted(group_pixels.items())}
                for key in bygroup:
                    a=base_summary[key];b=bygroup[key]
                    b['delta_thermal_fraction_ge_323_vs_1970_2000']=(
                        b['thermal_fraction_ge_323']-a['thermal_fraction_ge_323'])
                    b['delta_winter_gdd_mean_vs_1970_2000']=b['winter_gdd_mean']-a['winter_gdd_mean']
                report['scenarios'].append({'gcm':model,'ssp':ssp,'period':period,'groups':bygroup})
    if len(report['scenarios'])!=8 or any(s['groups']['all']['n_100km_hex_localities']!=len(centers)
                                         for s in report['scenarios']):
        raise ValueError('Incomplete ensemble or changing study-source support')
    report['inference_gate']={
        'thermal_projections_from_real_downscaled_climate':True,
        'host_leaf_phenology_measured':False,
        'butterfly_stage_survival_measured':False,
        'GBIF_sampling_effort_standardized':False,
        'independent_plant_stem_absence_denominator':False,
        'absolute_population_predictions_authorized':False,
        'claimed_scope':'Fixed geographies with observed monarch larval host photos, thermal exposure and explicit leaf-continuity scenarios only',
        'data_leak_prevention':'Only group aggregate metrics, no hex IDs, photo URIs, observation rows, coordinates, or observer identities',
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':report['status'],'n_hexes':len(centers),
                     'n_climate_pixels':len(unique_pixels),
                     'base':base_summary.get('all'),
                     'future':[{k:s[k] for k in ('gcm','ssp','period')}|
                               {'all':s['groups']['all']} for s in report['scenarios']]},sort_keys=True))
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--climate-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    compute(args)
