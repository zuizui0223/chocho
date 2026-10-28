#!/usr/bin/env python3
"""External Daymet weather and USA-NPN live leaves: held-out SITE forecasting.

Source-blind, effort-only station selection, 5-fold station held-out model.
No original site coordinates, person identifiers, observations or weather rows exported.
Not a future milkweed leaf, monarch population or climate causal model.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import time
from collections import Counter, defaultdict
from datetime import date, timedelta
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from audit_usanpn_milkweed_registry_20261028 import SPECIES, YEARS
from audit_usanpn_milkweed_leaf_status_20261028 import frozen_registry, get_records, parse_valid

SCHEMA = 'chocho_usanpn_daymet_site_holdout_leaf_20261028_v01'
DAYMET_URL = 'https://daymet.ornl.gov/single-pixel/api/data'
DAYMET_DOI = 'https://doi.org/10.3334/ORNLDAAC/2361'
DAYMET_VARS = ('tmax', 'tmin', 'prcp')
RETRIEVAL_YEARS = tuple(range(2018, 2026))
MAX_SITES = 72
PRIORITY_PER_SPECIES = 12
MIN_DISTINCT_DATES = 20
MIN_YEARS = 2
MIN_MONTHS = 4
N_FOLDS = 5
WEATHER_DAYS = 30
C_REG = 0.1
MAX_BYTES_PER_SITE = 2_000_000
MAX_NETWORK_FAILURES = 7
RANDOM_SEED = 20261028


def source_day_records(registry: Path, receipt: Path):
    selected, regsha = frozen_registry(registry)
    original = json.loads(receipt.read_text(encoding='utf-8'))
    if (original.get('status') != 'COMPLETE_PHENOPHASE_SOURCE_WITH_INFORMATION_GATES'
            or original.get('registry_sha256') != regsha):
        raise ValueError('Original eight-year NPN source provenance missing')
    old = {int(x['year']): x['response_sha256'] for x in original['requests']}
    if set(old) != set(RETRIEVAL_YEARS):
        raise ValueError('NPN source year coverage changed')
    site_xy = defaultdict(list)
    site_species_dates = defaultdict(set)
    records = []
    qc = Counter()
    accepted_ids = {v['species_id'] for v in selected.values()}
    for year in YEARS:
        raw, source = get_records(year, selected)
        if (source['status'] != 'COMPLETE_SOURCE_HTTP_200'
                or source['response_sha256'] != old[year]):
            raise ValueError('NPN original year source body changed: fail closed')
        eligible, checks = parse_valid(raw, year, selected)
        qc.update(checks)
        records.extend(eligible)
        for (sp, state, site, plant, observed_date), _ in eligible:
            site_species_dates[(site, sp)].add(observed_date)
        for item in raw:
            try:
                ident = int(item['site_id'])
                species_id = int(item['species_id'])
                lat = float(item['latitude'])
                lon = float(item['longitude'])
            except (ValueError, TypeError, KeyError):
                continue
            if species_id not in accepted_ids or ident <= 0:
                continue
            if not (-90 < lat < 90 and -180 < lon < 180):
                continue
            site_xy[ident].append((lat, lon))
        print(json.dumps({'original_usanpn_complete_year': int(year),
                          'raw_status_entries': len(raw),
                          'accepted_plant_days': len(eligible)}), flush=True)
        time.sleep(0.25)
    if len(records) != 27057 or qc['source_rows'] != 30476:
        raise ValueError('NPN original 27057 leaf statuses changed')
    coords = {}
    n_bad = 0
    for site, points in site_xy.items():
        y = sorted(v[0] for v in points)[len(points)//2]
        x = sorted(v[1] for v in points)[len(points)//2]
        if max(abs(lat-y) for lat, lon in points) > .03 or max(abs(lon-x) for lat, lon in points) > .03:
            n_bad += 1
            continue
        if 14.5 <= y <= 52.0 and -131 <= x <= -53:
            coords[site] = (y, x)
    return records, coords, site_species_dates, {
        'original_accepted_plant_date_records': len(records),
        'original_observations': qc['source_rows'],
        'source_complete_years': len(RETRIEVAL_YEARS),
        'source_year_receipt_sha256': hashlib.sha256(receipt.read_bytes()).hexdigest(),
        'source_registry_sha256': regsha,
        'excluded_nonstable_site_geometries': n_bad,
        'registered_species_count': len(selected),
        'no_site_id_observer_or_coordinates_exported': True,
    }


def choose_sites(source_dates, coords, max_sites=MAX_SITES):
    """Site selection is strictly sample availability; never check YES/NO flags."""
    by_species = defaultdict(list)
    effort = defaultdict(set)
    for (site, species), dates in source_dates.items():
        if site not in coords:
            continue
        if not (len(dates) >= MIN_DISTINCT_DATES
                and len({d.year for d in dates}) >= MIN_YEARS
                and len({d.month for d in dates}) >= MIN_MONTHS):
            continue
        by_species[species].append((site, len(dates)))
        effort[site].update(dates)
    for species in by_species:
        by_species[species].sort(key=lambda item: (-item[1], item[0]))
    selected = set()
    for species in SPECIES:
        for site, count in by_species[species][:PRIORITY_PER_SPECIES]:
            if len(selected) < max_sites:
                selected.add(site)
    remainder = sorted(effort, key=lambda site: (-len(effort[site]), site))
    for site in remainder:
        if len(selected) >= max_sites:
            break
        selected.add(site)
    if len(selected) > max_sites:
        raise ValueError('Selected too many original monitoring sites')
    return selected, {
        'eligible_site_species_groups_by_species': {s: len(by_species[s]) for s in SPECIES},
        'eligible_distinct_site_count': len(effort),
        'source_selected_site_count': len(selected),
        'source_selected_state_or_weather_not_used_to_rank': True,
    }


def _safe_daymet_decode(body):
    if not body or len(body) > MAX_BYTES_PER_SITE:
        raise ValueError('Daymet response empty/oversized')
    text = body.decode('utf-8-sig')
    if '<html' in text[:1000].lower() or 'login' in text[:1000].lower():
        raise ValueError('Daymet returned login or non-CSV response')
    lines = text.splitlines()
    headers = [i for i, line in enumerate(lines)
               if 'year' in line.lower() and 'yday' in line.lower()
               and 'tmax' in line.lower() and 'tmin' in line.lower()
               and 'prcp' in line.lower()]
    if len(headers) != 1:
        raise ValueError('Daymet CSV header is missing or ambiguous')
    it = csv.DictReader(io.StringIO('\n'.join(lines[headers[0]:])))
    keys = list(it.fieldnames or [])
    if len(keys) != 5:
        raise ValueError('Daymet source only five columns expected')
    mapping = {}
    for canonical in ('year', 'yday', 'tmax', 'tmin', 'prcp'):
        options = [key for key in keys if key.split(' (')[0].strip().lower() == canonical]
        if len(options) != 1:
            raise ValueError('Daymet missing/ambiguous source column ' + canonical)
        mapping[canonical] = options[0]
    result = {}
    byyear = Counter()
    for row in it:
        try:
            year = int(row[mapping['year']])
            yday = int(row[mapping['yday']])
            tmax = float(row[mapping['tmax']])
            tmin = float(row[mapping['tmin']])
            precip = float(row[mapping['prcp']])
        except (TypeError, ValueError) as exc:
            raise ValueError('Daymet invalid field') from exc
        if (year not in RETRIEVAL_YEARS or not 1 <= yday <= 365
                or not -90 < tmin <= tmax < 75
                or not 0 <= precip < 1000 or not all(map(math.isfinite, [tmin, tmax, precip]))):
            raise ValueError('Daymet years/date/physical values invalid')
        dt = date(year, 1, 1) + timedelta(days=yday-1)
        if dt in result:
            raise ValueError('Daymet duplicate calendar date')
        result[dt] = (tmin, tmax, precip)
        byyear[year] += 1
    if byyear != {year:365 for year in RETRIEVAL_YEARS}:
        raise ValueError('Daymet source did not return complete 8 years×365 calendar days')
    return result


def weather_for_site(lat, lon, fetcher=None, pause=0.45):
    if not (14.5 <= lat <= 52 and -131 <= lon <= -53):
        raise ValueError('NPN monitoring site outside Daymet North American coverage')
    data = {'lat': format(lat, '.6f'), 'lon': format(lon, '.6f'),
            'vars': ','.join(DAYMET_VARS),
            'years': ','.join(map(str, RETRIEVAL_YEARS))}
    if fetcher is None:
        def fetcher():
            # Exact coordinates stay in the encrypted GitHub runner process only.
            req = Request(DAYMET_URL + '?' + urlencode(data),
                          headers={'User-Agent': 'chocho-milkweed-daymet-validation/1.0',
                                   'Accept': 'text/csv,text/plain,*/*'})
            with urlopen(req, timeout=80) as resp:
                if resp.status != 200:
                    raise ValueError('Daymet HTTP response not 200')
                raw = resp.read(MAX_BYTES_PER_SITE+1)
            return raw
    err = None
    for attempt in range(3):
        try:
            raw = fetcher()
            values = _safe_daymet_decode(raw)
            return values, hashlib.sha256(raw).hexdigest()
        except (OSError, UnicodeError, ValueError, HTTPError, URLError) as exc:
            # Never serialize an exception or source URL containing a source station coordinate.
            err = type(exc).__name__
            time.sleep(pause*(attempt+1))
    raise ValueError('Daymet source cannot be validated (' + str(err) + ')')


def rolling_weather(daymet_daily, observed_date):
    days = [observed_date - timedelta(days=k) for k in range(WEATHER_DAYS)]
    if any(day not in daymet_daily for day in days):
        return None
    z = [daymet_daily[day] for day in days]
    mean_t = sum(.5*(a+b) for a,b,p in z)/WEATHER_DAYS
    precip = sum(p for a,b,p in z)
    frost = sum(a < 0 for a,b,p in z)
    return (mean_t, math.log1p(precip), frost)


def design_rows(observations, selected_sites, coordinates, daily_source):
    source = []
    qc = Counter()
    for (species, state, site, plant, observed_date), flag in observations:
        if site not in selected_sites or site not in daily_source:
            continue
        climate = rolling_weather(daily_source[site], observed_date)
        if climate is None:
            qc['missing_previous_30_day_weather_or_Daymet_leap_calendar']+=1
            continue
        if flag not in (0,1):
            raise ValueError('NPN nonbinary leaf status reached analysis')
        lat,lon = coordinates[site]
        day = observed_date.timetuple().tm_yday
        phase = 2*math.pi*(day-1)/365.25
        baseline = [float(species == name) for name in SPECIES] + [
            math.sin(phase),math.cos(phase),math.sin(2*phase),math.cos(2*phase),
            lat,lon,observed_date.year-2021.5,
        ]
        source.append((site,state,species,int(flag),baseline,baseline+list(climate)))
        qc['accepted_plant_date_outcome_and_weather']+=1
    return source,dict(qc)


def fold(site):
    return int(hashlib.sha256(str(site).encode('utf-8')).hexdigest()[:8],16) % N_FOLDS


def score_sites(probs, labels, site_ids):
    by_site=defaultdict(list)
    for prob,label,site in zip(probs,labels,site_ids):
        p=min(1-1e-8,max(1e-8,float(prob)))
        ll= -(int(label)*math.log(p)+(1-int(label))*math.log(1-p))
        br=(float(p)-int(label))**2
        by_site[site].append((ll,br))
    return (sum(sum(x[0] for x in values)/len(values) for values in by_site.values())/len(by_site),
            sum(sum(x[1] for x in values)/len(values) for values in by_site.values())/len(by_site))


def model_and_predict(train_x,train_y,train_site,test_x,*,lib='sklearn'):
    if len(set(train_y)) < 2:
        raise ValueError('Source model training has only one leaf class')
    scaler=StandardScaler()
    x=scaler.fit_transform(train_x)
    xt=scaler.transform(test_x)
    by=Counter(train_site)
    # Give each independent monitored station the same total weight.
    per_site_mass=len(train_x)/len(by)
    weights=np.asarray([per_site_mass/by[site] for site in train_site],dtype=float)
    fit=LogisticRegression(C=C_REG,max_iter=1000,solver='lbfgs',random_state=RANDOM_SEED)
    fit.fit(x,train_y,sample_weight=weights)
    return fit.predict_proba(xt)[:,1]


def evaluate_prediction(rows):
    sites={r[0] for r in rows}
    states={r[1] for r in rows}
    species={r[2] for r in rows}
    dates={}
    if len(sites)<30 or len(states)<5 or len(species)<3:
        return {'status':'INSUFFICIENT_INDEPENDENT_SITE_TAXON_COVERAGE',
                'n_sites':len(sites),'n_states':len(states),'n_species':len(species)}
    out=[]
    for f in range(N_FOLDS):
        train=[r for r in rows if fold(r[0]) != f]
        test=[r for r in rows if fold(r[0]) == f]
        test_sites={r[0] for r in test}
        labels=[r[3] for r in test]
        if len(test_sites)<5 or len(set(labels))<2 or len({r[0] for r in train})<20:
            return {'status':'INSUFFICIENT_FIVE_FOLD_HELD_OUT_SITE_INFORMATION',
                    'n_sites':len(sites),'n_states':len(states),'n_species':len(species),
                    'failed_fold':f,'n_test_sites':len(test_sites)}
        y=[r[3] for r in train]
        site=[r[0] for r in train]
        baseline=model_and_predict(np.array([r[4] for r in train]), y, site,
                                   np.array([r[4] for r in test]))
        weather=model_and_predict(np.array([r[5] for r in train]), y, site,
                                  np.array([r[5] for r in test]))
        b_ll,b_br=score_sites(baseline,labels,[r[0] for r in test])
        w_ll,w_br=score_sites(weather,labels,[r[0] for r in test])
        out.append({'fold':f,'test_sites':len(test_sites),'test_plant_visits':len(test),
                    'test_yes_visits':sum(labels),'test_no_visits':len(labels)-sum(labels),
                    'seasonal_geography_log_loss':b_ll,'seasonal_weather_log_loss':w_ll,
                    'heldout_log_loss_lift':b_ll-w_ll,
                    'seasonal_geography_brier':b_br,'seasonal_weather_brier':w_br,
                    'heldout_brier_lift':b_br-w_br})
    nsites=[a['test_sites'] for a in out]
    def combined(key):return sum(a[key]*a['test_sites'] for a in out)/sum(nsites)
    primary=combined('heldout_log_loss_lift')
    positive=sum(a['heldout_log_loss_lift']>0 for a in out)
    return {'status':('HELD_OUT_WEATHER_PREDICTIVE_INCREMENT_REPLICATED'
                      if primary>0 and positive>=4
                      else 'NO_REPLICATED_HELD_OUT_WEATHER_INCREMENT'),
            'n_sites':len(sites),'n_states':len(states),'n_species':len(species),
            'source_independent_site_folds':out,
            'site_equal_baseline_log_loss':combined('seasonal_geography_log_loss'),
            'site_equal_daymet_log_loss':combined('seasonal_weather_log_loss'),
            'site_equal_delta_log_loss_baseline_minus_daymet':primary,
            'site_equal_delta_brier_baseline_minus_daymet':combined('heldout_brier_lift'),
            'positive_weather_increment_folds':positive,
            'required_positive_folds':4,
            'no_statistical_p_value_or_demographic_effect':True}


def run(registry,receipt,output,max_sites=MAX_SITES,prechecked_weather=None):
    if not (30<=max_sites<=MAX_SITES):
        raise ValueError('Site cap cannot exceed pre-specified 72 sites')
    observations,coords,source_dates,qc=source_day_records(registry,receipt)
    selected,selection=choose_sites(source_dates,coords,max_sites)
    source_by_site=Counter(s for s,sp in source_dates if s in selected)
    state_by_site=defaultdict(set)
    for (sp,state,site,plant,dt),flag in observations:
        if site in selected:
            state_by_site[site].add(state)
    result={'schema':SCHEMA,'study_source':{'usanpn':qc,
            'daymet_public_service':DAYMET_URL,'daymet_citation':DAYMET_DOI,
            'daymet_years':list(RETRIEVAL_YEARS),'weather_day_lookback':WEATHER_DAYS,
            'source_site_effort_selection':selection},
            'model_specification':{'site_folds':N_FOLDS,'ridge_C':C_REG,
                'seasonal_harmonics':2,'source_season_year_and_latlon_baseline':True,
                'weather_additions':['mean30_temperature_C','log1p_sum30_precip_mm',
                                     'n30days_tmin_below_0C'],
                'training_fold_only_scaling':True,
                'site_equal_training_and_heldout_scoring':True},
            'daymet_data_not_exported':True,
            'original_plant_site_individual_identifiers_and_locations_not_exported':True,
            'daymet_original_full_site_request_count':0}
    if len(selected)<30:
        result['decision']='INSUFFICIENT_SOURCE_EFFORT_SELECTED_SITES'
    else:
        weather={}
        digests=[]
        failures=0
        for i,site in enumerate(sorted(selected),1):
            lat,lon=coords[site]
            try:
                series,digest=weather_for_site(lat,lon)
                weather[site]=series
                digests.append(digest)
            except (ValueError,OSError):
                failures+=1
            if i%12 == 0 or i == len(selected):
                print(json.dumps({'source_weather_sites_processed':i,
                                  'source_weather_sites_retrieved':len(weather),
                                  'source_weather_site_rejections':failures}),flush=True)
            if failures > MAX_NETWORK_FAILURES:
                break
            time.sleep(.35)
        result['daymet_original_full_site_request_count']=min(len(selected),len(weather)+failures)
        result['daymet_complete_sites']=len(weather)
        result['daymet_rejected_sites']=failures
        result['daymet_aggregate_sha256']=(hashlib.sha256(('\n'.join(sorted(digests))).encode()).hexdigest()
                                          if digests else None)
        if len(weather) < .8*len(selected) or len(weather) < 30 or failures > MAX_NETWORK_FAILURES:
            result['decision']='DAYMET_SOURCE_COVERAGE_OR_ACCESS_GATE_FAILED'
        else:
            rows,counts=design_rows(observations,selected,coords,weather)
            result['eligible_response_weather']=counts
            result['quality_after_complete_weather']={
                'response_days':len(rows),
                'sites':len({r[0] for r in rows}),
                'states':len({r[1] for r in rows}),
                'host_species':len({r[2] for r in rows}),
                'yes_statuses':sum(r[3] for r in rows),
                'no_statuses':sum(not r[3] for r in rows),
                'source_species_counts':{sp:sum(r[2]==sp for r in rows) for sp in SPECIES},
                'site_identifiers_not_exported':True,
            }
            result['cv']=evaluate_prediction(rows)
            result['decision']=result['cv']['status']
    result['phenology_future_climate_projection_authorized']=(
        result['decision']=='HELD_OUT_WEATHER_PREDICTIVE_INCREMENT_REPLICATED')
    result['interpretation_limits']=[
      'Held-out-site weather prediction is observational predictive evidence, not temperature causation',
      'Nature\'s Notebook monitored plants are a nonrandom sample; missing leaves visits never encoded no',
      'Weather changes may be confounded with observation dates, horticultural management, irrigation, state or elevation',
      'This result does not predict actual future leaf senescence from 20-year monthly CMIP6',
      'No milkweed edible biomass, monarch egg/larval survival or monarch population model is estimated',
      'Source coordinates, individual observations and Daymet daily climate are not redistributed',
    ]
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'decision':result['decision'],
                      'daymet_complete_sites':result.get('daymet_complete_sites'),
                      'source_cohort':result.get('quality_after_complete_weather'),
                      'cv_gate':result.get('cv',{}).get('status'),
                      'heldout_delta_logloss':result.get('cv',{}).get('site_equal_delta_log_loss_baseline_minus_daymet')},
                     sort_keys=True),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--registry',type=Path,required=True)
    p.add_argument('--npn-source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    result=run(args.registry,args.npn_source,args.output)
    if result['decision'].startswith('DAYMET_SOURCE') or result['decision']=='INSUFFICIENT_SOURCE_EFFORT_SELECTED_SITES':
        raise SystemExit(2)
