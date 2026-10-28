#!/usr/bin/env python3
"""CMIP6 thermal screening at independently USA-NPN-monitored milkweed sites.

Observed month-specific leaf yes/no is held FIXED as a counterfactual.
Neither winter leaf change under warming nor monarch population outcomes are modeled.
No raw site coordinates, individuals, observation IDs or exact locations exported.
"""
from __future__ import annotations
import argparse, hashlib, json, math, statistics, time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
from audit_usanpn_milkweed_registry_20261028 import YEARS
from audit_usanpn_milkweed_leaf_status_20261028 import frozen_registry,get_records,parse_valid
from analyze_usanpn_winter_leaf_continuity_20261028 import site_winter_months,WINTER
from simulate_monarch_cmip6_thermal_window_20261028 import (
    VARS,MODELS,SSPS,PERIODS,DAYS,LOWER_C,UPPER_C,ONE_GENERATION_GDD,
    archive_name,sha256_file,grid_cell,read_sampled_archive,
    choose_temperature_scale,climate_validity_mask,maximas_bool
)

TARGET=('Asclepias syriaca','Asclepias tuberosa')
SCHEMA='chocho_usanpn_same_site_cmip6_fixed_monthly_leaf_scenarios_20261028_v01'
EXPECTED_CASES={'Asclepias syriaca':30,'Asclepias tuberosa':21}


def acquire_source(registry, original_receipt):
    selected,regsha=frozen_registry(registry)
    old=json.loads(original_receipt.read_text())
    if old['status']!='COMPLETE_PHENOPHASE_SOURCE_WITH_INFORMATION_GATES' or old['registry_sha256']!=regsha:
        raise ValueError('Previous complete source receipt missing')
    by_year={int(x['year']):x['response_sha256'] for x in old['requests']}
    accepted=[];coords=defaultdict(list);qc=Counter()
    for year in YEARS:
        rows,receipt=get_records(year,selected)
        if receipt['status']!='COMPLETE_SOURCE_HTTP_200' or receipt['response_sha256']!=by_year[year]:
            raise ValueError(f'Original USA-NPN source response drifted at {year}, no automatic relabel')
        valid,checks=parse_valid(rows,year,selected)
        accepted.extend(valid)
        qc.update(checks)
        for r in rows:
            try:
                spid=int(r.get('species_id'));site=int(r.get('site_id'))
                lat=float(r.get('latitude'));lon=float(r.get('longitude'))
            except (ValueError,TypeError):continue
            if site<=0 or not(-90<lat<90 and -180<lon<180):continue
            if spid not in {selected[z]['species_id'] for z in TARGET}:continue
            # NPN site ID is a repeated geographical station, not an individual.
            coords[site].append((lat,lon))
        print(json.dumps({'year':year,'source_rows':len(rows),'accepted_unique_individual_dates':len(valid)}),flush=True)
        time.sleep(.6)
    if len(accepted)!=27057 or qc['source_rows']!=30476:
        raise ValueError('NPN original eight-year accepted source denominator changed')
    sites={};bad=set()
    for site,xy in coords.items():
        lat=statistics.median(a for a,b in xy)
        lon=statistics.median(b for a,b in xy)
        if max(abs(a-lat) for a,b in xy)>.03 or max(abs(b-lon) for a,b in xy)>.03:
            bad.add(site);continue
        sites[site]=(lat,lon)
    return accepted,sites,{'registered_taxa':list(TARGET),'original_receipt_sha256':hashlib.sha256(original_receipt.read_bytes()).hexdigest(),
                           'source_year_replies_exact_sha_verified':len(YEARS),
                           'unique_source_records':len(accepted),
                           'source_site_geo_conflict_count':len(bad),
                           'source_sites_with_stable_coordinates':len(sites)}


def masks_from_observed_site_month(records,sites):
    monthly=site_winter_months(records)
    by_case=defaultdict(dict)
    for (sp,state,site,win,month),data in monthly.items():
        if sp not in TARGET:continue
        by_case[(sp,state,site,win)][month]=data
    cases=[];qc=Counter()
    for (sp,state,site,win),values in sorted(by_case.items()):
        if not all(m in values and len(values[m]['dates'])>=2 for m in WINTER):continue
        qc[sp]+=1
        loc=sites.get(site)
        if loc is None:
            qc[f'{sp}_no_stable_geocoordinate']+=1;continue
        pix=grid_cell(*loc)
        if pix is None:
            qc[f'{sp}_invalid_grid']+=1;continue
        optimistic={m:int(values[m]['yes']>0) for m in WINTER}
        conservative={m:int(values[m]['yes']>0 and values[m]['no']==0) for m in WINTER}
        if any(conservative[m]>optimistic[m] for m in WINTER):raise ValueError('Leaf source mask order invalid')
        # Site ID used only internally to avoid repeated winter pseudoreplication.
        cases.append({'species':sp,'site':site,'state':state,'winter':win,'pixel':pix,
                      'observed_any_yes_month':optimistic,
                      'observed_all_yes_only_month':conservative})
    for sp,n in EXPECTED_CASES.items():
        if qc[sp]!=n: raise ValueError(f'Frozen 4-month monitored site-year count changed for {sp}')
    if sum(qc[sp] for sp in TARGET)!=51:raise ValueError('Original full winter cohort changed')
    if len(cases)<30:raise ValueError('Too many original full-winter site-years have invalid locations')
    return cases,dict(qc)


def seasonal_degree_days(minima,maxima):
    out={}
    for m in WINTER:
        avg=.5*(np.asarray(minima[m],dtype=float)+np.asarray(maxima[m],dtype=float))
        out[m]=np.maximum(0.,np.minimum(avg,UPPER_C)-LOWER_C)*DAYS[m]
    return out


def evaluate(cases,gdd,pixel_index,valid_cases=None):
    """Site-equal average of winter-year binary screens, not 51 independent populations."""
    per_site=defaultdict(list)
    per_species=defaultdict(set)
    for case in cases:
        sp=case['species']
        pix=pixel_index[case['pixel']]
        total=sum(float(gdd[m][pix]) for m in WINTER)
        optimistic=sum(float(gdd[m][pix])*case['observed_any_yes_month'][m] for m in WINTER)
        conservative=sum(float(gdd[m][pix])*case['observed_all_yes_only_month'][m] for m in WINTER)
        if not 0<=conservative<=optimistic+1e-8 or optimistic>total+1e-8:
            raise ValueError('Observed-leaf scenario requires no more heat than thermal-only')
        per_site[(sp,case['site'])].append({
            'thermal_only':int(total>=ONE_GENERATION_GDD),
            'optimistic_observed_months':int(optimistic>=ONE_GENERATION_GDD),
            'conservative_observed_months':int(conservative>=ONE_GENERATION_GDD),
            'thermal_gdd':total,'optimistic_leaf_gdd':optimistic,'conservative_leaf_gdd':conservative})
        per_species[sp].add(case['pixel'])
    out={}
    for sp in TARGET:
        sites=[ls for (species,site),ls in per_site.items() if species==sp]
        if not sites:
            out[sp]={'n_sites':0,'information':'no retained sites'}
            continue
        estimands=('thermal_only','optimistic_observed_months','conservative_observed_months')
        site_scores={key:[sum(x[key] for x in windows)/len(windows) for windows in sites] for key in estimands}
        out[sp]={
            'n_independent_sites':len(sites),
            'n_site_winter_years':sum(len(x) for x in sites),
            'n_distinct_worldclim_10min_pixels':len(per_species[sp]),
            'site_equal_proportion_threshold_thermal_only':sum(site_scores['thermal_only'])/len(sites),
            'site_equal_proportion_threshold_optimistic_monthly_leaf_mask':sum(site_scores['optimistic_observed_months'])/len(sites),
            'site_equal_proportion_threshold_conservative_monthly_leaf_mask':sum(site_scores['conservative_observed_months'])/len(sites),
            'site_equal_mean_thermal_gdd':sum(sum(y['thermal_gdd'] for y in x)/len(x) for x in sites)/len(sites),
            'per_site_years_not_independent':True}
        d=out[sp]
        if not d['site_equal_proportion_threshold_conservative_monthly_leaf_mask']<=d['site_equal_proportion_threshold_optimistic_monthly_leaf_mask']<=d['site_equal_proportion_threshold_thermal_only']+1e-10:
            raise ValueError('Source leaf GDD physics consistency lost')
    return out


def run(args):
    accepted,coords,source_qc=acquire_source(args.registry,args.npn_source)
    cases,case_qc=masks_from_observed_site_month(accepted,coords)
    original_json=json.loads(args.cmip6_source.read_text(encoding='utf-8'))
    if original_json['status']!='SOURCE_LOCKED_CMIP6_THERMAL_SCENARIO_COMPLETE':
        raise ValueError('Previous source-validated CMIP6 source receipt missing')
    expected_archives={k:v['sha256'] for k,v in original_json['archives'].items()}
    pairs=sorted({c['pixel'] for c in cases})
    all_panels=[];archives={}
    for model,ssp,period in [('1970-2000','baseline','historical')]+[(m,s,p) for m in MODELS for s in SSPS for p in PERIODS]:
        sampled={}
        for var in VARS:
            name=archive_name(var) if model=='1970-2000' else archive_name(var,model,ssp,period)
            path=args.climate_dir/name
            if sha256_file(path)!=expected_archives[name]:raise ValueError('Climate original archive source hash changed '+name)
            raw,scales,digest,sorted_pairs=read_sampled_archive(path,pairs)
            if sorted_pairs!=pairs:raise ValueError('Climate source grid does not match fixed NPN site pixels')
            archives[name]={'sha256':digest,'bytes':path.stat().st_size}
            sampled[var]=(raw,scales)
        scale=choose_temperature_scale(sampled['tmax'][0][7],sampled['tmax'][1].values())
        if model=='1970-2000':base_scale=scale
        elif not math.isclose(scale,base_scale,abs_tol=1e-8):
            raise ValueError('Baseline/future WorldClim source scale mismatch')
        tmin={m:sampled['tmin'][0][m]*scale for m in range(1,13)}
        tmax={m:sampled['tmax'][0][m]*scale for m in range(1,13)}
        good,bad=climate_validity_mask(tmin,tmax)
        all_panels.append({'model':model,'ssp':ssp,'period':period,'tmin':tmin,'tmax':tmax,'good':good,'bad':bad})
    common=np.ones(len(pairs),dtype=bool)
    for p in all_panels:common&=p['good']
    retained=int(np.sum(common))
    if retained<.9*len(pairs):raise ValueError('More than 10% of NPN-site climate pixels failed common physical QA')
    final_pairs=[p for p,good in zip(pairs,common) if good]
    pixel_index={p:i for i,p in enumerate(final_pairs)}
    cases=[c for c in cases if c['pixel'] in pixel_index]
    if len(cases)<30:raise ValueError('Insufficient source-matched climates after common physical QA')
    projections=[]
    for panel in all_panels:
        tmin={m:v[common] for m,v in panel['tmin'].items()}
        tmax={m:v[common] for m,v in panel['tmax'].items()}
        gdd=seasonal_degree_days(tmin,tmax)
        projected=evaluate(cases,gdd,pixel_index)
        projections.append({'model':panel['model'],'ssp':panel['ssp'],'period':panel['period'],'by_milkweed':projected})
    report={
        'schema':SCHEMA,'status':'EXPLORATORY_REAL_CMIP6_FIXED_OBSERVED_LEAF_MONTH_SITE_OVERLAY',
        'source':source_qc,
        'source_site_winter_gate':case_qc,
        'physical_common_site_grid_quality':{
            'NPN_unique_source_climate_pixels':len(pairs),'retained_valid_grid_pixels':len(final_pairs),
            'excluded_invalid_grid_pixels':len(pairs)-len(final_pairs),
            'retained_observed_winter_site_years':len(cases),
            'all_nine_climate_panels_share_exact_same_source_geography':True,
        },
        'CMIP6_original_source_archive_sha256':archives,
        'baseline_and_eight_scenarios':projections,
        'inference':{
            'original_monarch_840_hexes_directly_joined':False,
            'leaf_availability_future_projected':False,
            'observed_NPN_month_leaf_status_is_FIXED_counterfactual':True,
            'optimistic_leaf_status_assumes_whole_month_if_any_sampled_date_yes':True,
            'conservative_leaf_status_assumes_whole_month_only_if_all_sampled_dates_yes':True,
            'neither_monthly_mask_observes_unsampled_days_or_plant_biomass':True,
            'thermal_condition_is_not_monarch_development_or_population_growth':True,
            'source_years_2018to2025_vs_WorldClim_baseline_1970to2000_epoch_mismatch':True,
            'site_sample_is_nonrandom_and_small':True,
            'no_individual_observer_site_identifiers_or_coordinates_exported':True,
        }
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'status':report['status'],'source':source_qc,
        'common_quality':report['physical_common_site_grid_quality'],
        'baseline':projections[0]['by_milkweed'],
        'late_ssp585':[(x['model'],x['by_milkweed']) for x in projections if x['ssp']=='ssp585' and x['period']=='2081-2100']},ensure_ascii=False),flush=True)
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--registry',type=Path,required=True)
    p.add_argument('--npn-source',type=Path,required=True)
    p.add_argument('--cmip6-source',type=Path,required=True)
    p.add_argument('--climate-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a)
