#!/usr/bin/env python3
"""USA-NPN exact live Leaves status among six preregistered Asclepias species.

No raw observations, individual identifiers, site coordinates, or observer data
are written. A no-record month is UNSAMPLED, never assigned leaf absence.
"""
from __future__ import annotations
import argparse, hashlib, json, re, time
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

from audit_usanpn_milkweed_registry_20261028 import SPECIES, YEARS, BASE

API = BASE + 'observations/getObservations.json'
MAX_BYTES = 55_000_000
TIMEOUT = 110
REQ_SOURCE = 'chocho NPN leaf phenology audit'


def frozen_registry(path):
    d=json.loads(path.read_text(encoding='utf-8'))
    if d.get('schema')!='chocho_usanpn_milkweed_live_leaves_registry_v01' or d.get('status')!='REGISTRY_AND_EXACT_LIVE_LEAVES_PROTOCOL_RESOLVED':
        raise ValueError('USA-NPN registry provenance is not verified')
    if tuple(d.get('calendar_years',()))!=YEARS or set(d.get('taxa',{}))!=set(SPECIES):
        raise ValueError('Frozen USA-NPN six-species/2018–2025 roster modified')
    accepted={}
    for name,obj in d['taxa'].items():
        spid=obj.get('species_id')
        if spid is None:continue
        if not isinstance(spid,int) or spid<=0:raise ValueError('Invalid registry species id')
        ids=obj['leaves_pheno_id_by_year']
        if set(ids)!={str(x) for x in YEARS}:
            raise ValueError('NPN leaves phase mapping incomplete')
        accepted[name]={'species_id':spid,'year_to_leaf_id':ids}
    if not accepted:raise ValueError('No accepted milkweed species for querying')
    return accepted,hashlib.sha256(path.read_bytes()).hexdigest()


def get_records(year, selected, request_fn=None):
    """Query only registry-confirmed taxa and their year-exact Leaves phenophase."""
    ids={int(obj['year_to_leaf_id'][str(year)]) for obj in selected.values() if obj['year_to_leaf_id'].get(str(year)) is not None}
    if not ids:
        return [],{'status':'NO_APPLICABLE_LEAVES_PROTOCOL','year':year}
    keys=sorted(obj['species_id'] for obj in selected.values() if obj['year_to_leaf_id'].get(str(year)) is not None)
    params={'request_src':REQ_SOURCE,'climate_data':'0','start_date':f'{year}-01-01','end_date':f'{year}-12-31'}
    for i,k in enumerate(keys,1):params[f'species_id[{i}]']=str(k)
    for i,k in enumerate(sorted(ids),1):params[f'phenophase_id[{i}]']=str(k)
    if request_fn is None:
        def request_fn():
            data=urlencode(params).encode('utf-8')
            req=Request(API, data=data, headers={'User-Agent':'chocho-usanpn-leaf-status/1.0',
                'Accept':'application/json','Content-Type':'application/x-www-form-urlencoded'})
            with urlopen(req,timeout=TIMEOUT) as response:
                if response.status!=200:raise ValueError('HTTP '+str(response.status))
                if 'json' not in response.headers.get('Content-Type','').lower():
                    raise ValueError('USA-NPN returned non-JSON source')
                raw=response.read(MAX_BYTES+1)
            return raw
    raw=request_fn()
    if len(raw)>MAX_BYTES:raise ValueError('USA-NPN response too large for bounded source')
    try:records=json.loads(raw)
    except json.JSONDecodeError as exc: raise ValueError('Malformed JSON source') from exc
    if not isinstance(records,list): raise ValueError('NPN observations must be a JSON list')
    # Source query receipt only, no raw body redistribution.
    return records, {'status':'COMPLETE_SOURCE_HTTP_200','year':year,
        'response_sha256':hashlib.sha256(raw).hexdigest(),'source_records':len(records),
        'requested_species':keys,'requested_phenophases':sorted(ids),'response_bytes':len(raw)}


def parse_valid(records, year, selected):
    """Identity-based QC: status 0 is explicitly observed no, not missing."""
    sid_to_name={obj['species_id']:name for name,obj in selected.items()}
    unique={}
    qc=Counter()
    for r in records:
        qc['source_rows']+=1
        if not isinstance(r,dict):raise ValueError('NPN observation not object')
        try:
            spid=int(r['species_id'])
            pid=int(r['phenophase_id'])
        except (ValueError,TypeError,KeyError) as exc:
            raise ValueError('NPN source missing exact species/phase id') from exc
        if spid not in sid_to_name: raise ValueError('Unexpected NPN taxon identity')
        name=sid_to_name[spid]
        if pid != int(selected[name]['year_to_leaf_id'][str(year)]):
            raise ValueError('Unexpected phenophase in leaves source')
        d=str(r.get('observation_date',''))
        try:dt=date.fromisoformat(d)
        except ValueError:qc['invalid_dates']+=1;continue
        if dt.year!=year: raise ValueError('Response includes nonrequested year')
        state=str(r.get('state','') or '').strip()
        if not re.fullmatch(r'[A-Z]{2}',state):qc['invalid_state']+=1;continue
        site=r.get('site_id'); individual=r.get('individual_id')
        try:site=int(site);individual=int(individual)
        except (ValueError,TypeError):qc['invalid_individual_or_site_id']+=1;continue
        if site<=0 or individual<=0:qc['invalid_individual_or_site_id']+=1;continue
        status=str(r.get('phenophase_status',''))
        if status not in ('0','1'):
            qc['unknown_or_uncertain_status']+=1
            continue
        key=(name,state,site,individual,dt)
        val=int(status)
        if key in unique and unique[key] != val:
            unique[key]=None  # conflict remains unresolved; never impute a no
            qc['conflicting_same_individual_day']+=1
        elif key not in unique:
            unique[key]=val
        else:qc['duplicate_same_individual_day']+=1
    accepted=[(key,flag) for key,flag in unique.items() if flag in (0,1)]
    qc['accepted_unique_plant_dates']=len(accepted)
    qc['conflict_or_missing_after_unique']=len(unique)-len(accepted)
    return accepted,dict(qc)


def aggregate(all_entries, qc):
    """Build observed phenology only; visit or missing month != whole-month leaf presence."""
    by_individual_month=defaultdict(lambda:{'dates':set(),'yes':0,'no':0})
    by_site_month=defaultdict(lambda:{'dates':set(),'yes':0,'no':0})
    by_winter=defaultdict(lambda:{'dates':set(),'yes':0,'no':0})
    states_by_species=defaultdict(set)
    visit_totals=Counter()
    for (species,state,site,ind,dt),flag in all_entries:
        m=dt.month; year=dt.year; win=year+1 if m>=11 else year
        k=(species,state,site,ind,year,m)
        x=by_individual_month[k]; x['dates'].add(dt);x['yes' if flag==1 else 'no']+=1
        k2=(species,state,site,year,m)
        y=by_site_month[k2];y['dates'].add(dt);y['yes' if flag==1 else 'no']+=1
        if m in (11,12,1,2) and win in range(2019,2026):
            k3=(species,state,site,win)
            z=by_winter[k3];z['dates'].add(dt);z['yes' if flag==1 else 'no']+=1
        states_by_species[species].add(state)
        visit_totals[species,'yes' if flag==1 else 'no']+=1
    month=defaultdict(lambda:Counter())
    for (species,state,site,year,m),v in by_site_month.items():
        key=(species,year,m)
        z=month[key]
        z['sites_any_checked']+=1
        if len(v['dates'])>=2:
            z['sites_at_least_two_dates']+=1
            if v['yes'] and v['no']:z['mixed_status_site_months']+=1
            elif v['yes']:z['all_yes_site_months']+=1
            elif v['no']:z['all_no_site_months']+=1
    winters=[];gate={}
    for species in SPECIES:
        w=[(s,v) for (sp,state,site,win),v in by_winter.items()
           if sp==species and len(v['dates'])>=2 for s in [state]]
        states={st for st,v in w}
        yes=sum(v['yes'] for _,v in w)
        no=sum(v['no'] for _,v in w)
        n=len(w)
        status=('PASS_NPN_DIRECT_WINTER_LEAF_SUPPORT' if n>=15 and len(states)>=3 and yes>0 and no>0
                else 'FAIL_NPN_WINTER_SITE_YEAR_INFORMATION_GATE')
        winters.append({'species':species,'two_date_winter_site_years':n,
             'states_two_date_winter_coverage':len(states),
             'explicit_winter_yes_visits':yes,'explicit_winter_no_visits':no,
             'status':status,'raw_unique_leaf_yes_observation_dates':visit_totals[species,'yes'],
             'raw_unique_leaf_no_observation_dates':visit_totals[species,'no'],
             'any_recorded_states':len(states_by_species[species])})
        gate[species]=status
    monthly=[{'species':sp,'year':year,'month':m,**dict(c)} for (sp,year,m),c in sorted(month.items())]
    return {'monthly_source':monthly,'winter_data_gate':winters,'species_winter_gates':gate,
            'unique_site_months':len(by_site_month),'unique_individual_months':len(by_individual_month),
            'n_accepted_observed_individual_days':len(all_entries)}


def run(registry,output):
    selected,digest=frozen_registry(registry)
    all_entries=[];receipts=[];quality=Counter();failures=[]
    for year in YEARS:
        try:
            records,receipt=get_records(year,selected)
            if receipt['status']=='NO_APPLICABLE_LEAVES_PROTOCOL':
                receipts.append(receipt);continue
            valid,qc=parse_valid(records,year,selected)
            all_entries.extend(valid)
            quality.update(qc)
            receipts.append(receipt)
            print(json.dumps({'year':year,'retrieved_rows':len(records),
                'accepted_plant_dates':len(valid),'source_sha256':receipt['response_sha256']},sort_keys=True),flush=True)
            time.sleep(.8)
        except (OSError,ValueError,TimeoutError,HTTPError,URLError) as exc:
            failures.append({'year':year,'failure_class':type(exc).__name__,'message':str(exc)[:230]})
            break # no silent partial result
    result={'schema':'chocho_usanpn_six_milkweeds_live_leaf_status_20261028_v01',
         'source_usa_npn':'https://services.usanpn.org/npn_portal/',
         'registered_species':list(SPECIES),'years':list(YEARS),
         'registry_sha256':digest,'requests':receipts,
         'source_failure':failures,'qc':dict(quality),
         'no_raw_locations_observer_names_or_plant_id_exported':True,
         'disclaimer':'Direct observational live-leaf yes/no only, not biomass, host use or future senescence.'}
    if not failures and sum(x['status']=='COMPLETE_SOURCE_HTTP_200' for x in receipts)>=6:
        result.update(aggregate(all_entries,quality))
        result['status']='COMPLETE_PHENOPHASE_SOURCE_WITH_INFORMATION_GATES'
    else:
        result['status']='SOURCE_UNAVAILABLE_OR_INCOMPLETE'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'source_years_completed':len([r for r in receipts if r['status']=='COMPLETE_SOURCE_HTTP_200']),
             'usable_observation_days':len(all_entries),
             'winter_gates':result.get('species_winter_gates',{}),'failed':failures},sort_keys=True),flush=True)
    if result['status']!='COMPLETE_PHENOPHASE_SOURCE_WITH_INFORMATION_GATES':raise SystemExit(2)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--registry',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.registry,a.output)
