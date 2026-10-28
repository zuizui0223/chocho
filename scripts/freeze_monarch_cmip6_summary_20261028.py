#!/usr/bin/env python3
"""Validate exact source-only CMIP6 scenario receipt and freeze aggregate summary.

No observation or model reruns; no plot-level/site-level/observer data exported.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path

SOURCE_SCHEMA='chocho_monarch_cmip6_recorded_locality_thermal_opportunity_v01'


def verify(payload):
    if payload.get('schema')!=SOURCE_SCHEMA or len(payload.get('scenarios',[]))!=8:
        raise ValueError('Incomplete original source ensemble')
    if len(payload.get('archives',{}))!=18:raise ValueError('Incomplete source archive checksum ledger')
    source=payload['source']
    if source['monarch_author_zip_sha256']!='4e4def63669b34738bcbfb8ba32ccf6c61d02cc6992b8f1ebfa8bc350a3987db':
        raise ValueError('Monarch source digest drift')
    if source['source_thinned_qc']['source_rows']!=12942 or source['source_core_qc']['core_rows']!=30320:
        raise ValueError('Original author source population drift')
    if source['source_hexes_after_common_climate_qa']!=840 or source['sampled_worldclim_pixels_after_common_qa']!=835:
        raise ValueError('Physical common-source QA unexpectedly changed')
    assert source['common_source_climate_quality_audit']['excluded_pixel_count']==1
    seen=set()
    assert payload['model']['thermal_egg_to_adult_target_degree_days']==323.
    for s in payload['scenarios']:
        key=(s['gcm'],s['ssp'],s['period'])
        if key in seen:raise ValueError('Duplicate scenario')
        seen.add(key)
        a=s['groups']['all']
        b=payload['baseline']['all']
        if a['n_100km_hex_localities']!=840 or b['n_100km_hex_localities']!=840:
            raise ValueError('Changing source sample across scenario')
        for n in ('thermal_fraction_ge_323','winter_available_50pct_323_fraction','summer_mean_tmax_over38_screen_fraction'):
            if not 0<=a[n]<=1.:raise ValueError('Invalid source-derived fraction')
        if a['winter_available_50pct_323_fraction']>a['thermal_fraction_ge_323']:
            raise ValueError('Impossible leaf constraint improves potential thermal accessibility')
        if a['winter_thermal_323_count']!=round(a['thermal_fraction_ge_323']*840):
            raise ValueError('Numerator and denominator mismatch')
        if abs(a['delta_thermal_fraction_ge_323_vs_1970_2000'] -
                (a['thermal_fraction_ge_323']-b['thermal_fraction_ge_323']))>1e-10:
            raise ValueError('Future-baseline derived delta wrong')
    if len(seen)!=8:
        raise ValueError('Future GCM/SSP/period product incomplete')
    for key in ('absolute_population_predictions_authorized','host_leaf_phenology_measured',
                'butterfly_stage_survival_measured'):
        if payload['inference_gate'][key] is not False:
            raise ValueError('Unsupported biology promoted')
    forbidden=('"uri":','"observer":','"latitude":','"longitude":','"hex_id":')
    data=json.dumps(payload,sort_keys=True).lower()
    if any(word in data for word in forbidden):raise ValueError('Source private geography or observer leak')
    return True


def rows_from_source(payload):
    base=payload['baseline']['all']
    rows=[]
    def to_row(gcm,ssp,period,r):
        return {
            'gcm':gcm,'ssp':ssp,'period':period,
            'n_source_100km_hexes':r['n_100km_hex_localities'],
            'mean_winter_GDD':round(r['winter_gdd_mean'],3),
            'winter_323GDD_screen_pct':round(100*r['thermal_fraction_ge_323'],3),
            'winter_half_leaf_323GDD_screen_pct':round(100*r['winter_available_50pct_323_fraction'],3),
            'winter_NovDec_only_323GDD_screen_pct':round(100*r['novdec_only_fraction_ge_323'],3),
            'winter_JanFeb_only_323GDD_screen_pct':round(100*r['janfeb_only_fraction_ge_323'],3),
            'summer_monthly_mean_Tmax_above_38C_pct':round(100*r['summer_mean_tmax_over38_screen_fraction'],3),
            'pct_point_change_winter_323GDD_vs_baseline':round(100*(r['thermal_fraction_ge_323']-base['thermal_fraction_ge_323']),3),
        }
    rows.append(to_row('historical','none','1970-2000',base))
    for s in payload['scenarios']:
        rows.append(to_row(s['gcm'],s['ssp'],s['period'],s['groups']['all']))
    return rows


def freeze(input_path,output_json,output_csv):
    raw=Path(input_path).read_bytes();d=json.loads(raw)
    verify(d)
    output_json.parent.mkdir(parents=True,exist_ok=True)
    output_csv.parent.mkdir(parents=True,exist_ok=True)
    output_json.write_bytes(raw)
    rows=rows_from_source(d)
    with output_csv.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    print(json.dumps({'status':'EXACT_SOURCE_SOURCE_ONLY_CMIP6_SCENARIOS_FROZEN',
                      'author_source_hexes':d['source']['author_spatially_thinned_hexes'],
                      'modelled_hexes':d['baseline']['all']['n_100km_hex_localities'],
                      'scenarios':len(d['scenarios']),
                      'result_sha256':hashlib.sha256(raw).hexdigest(),
                      'csv_sha256':hashlib.sha256(output_csv.read_bytes()).hexdigest()}))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--json',type=Path,required=True)
    p.add_argument('--csv',type=Path,required=True)
    a=p.parse_args();freeze(a.input,a.json,a.csv)
