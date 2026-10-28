#!/usr/bin/env python3
"""Posthoc observation coverage check before any leaf-constrained CMIP6 reforecast.

Ask whether each source site was monitored throughout ALL FOUR winter months.
Records from weeks/months with no visits are unknown, NOT leaf absences.
"""
from __future__ import annotations
import argparse, hashlib, json
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

from audit_usanpn_milkweed_registry_20261028 import SPECIES,YEARS
from audit_usanpn_milkweed_leaf_status_20261028 import frozen_registry,get_records,parse_valid

WINTER=(11,12,1,2)

def site_winter_months(records):
    # Site-month is a conservative geographic observation unit. Raw keys
    # remain in memory only and are discarded before writing summaries.
    agg=defaultdict(lambda:{'dates':set(),'yes':0,'no':0})
    for (species,state,site,individual,day),flag in records:
        if day.month not in WINTER:continue
        wy=day.year+1 if day.month in (11,12) else day.year
        if wy<2019 or wy>2025:continue
        obj=agg[(species,state,site,wy,day.month)]
        obj['dates'].add(day)
        obj['yes' if flag==1 else 'no']+=1
    return agg

def complete_winter_gate(site_month):
    by_winter=defaultdict(dict)
    by_state_month=defaultdict(Counter)
    for (sp,state,site,wy,month),x in site_month.items():
        by_winter[(sp,state,site,wy)][month]=x
        if len(x['dates'])>=2:
            e=by_state_month[(sp,state,wy,month)]
            e['sampled_site_months']+=1
            if x['yes'] and x['no']:
                e['mixed_sampled_status']+=1
            elif x['yes']:
                e['all_yes_sampled_dates']+=1
            elif x['no']:
                e['all_no_sampled_dates']+=1
    results=[]
    for species in SPECIES:
        selected=[(state,winter) for (sp,state,site,wy),winter in by_winter.items() if sp==species]
        all4=[(state,m) for state,m in selected if all(month in m for month in WINTER)]
        strict=[(state,m) for state,m in all4 if all(len(m[month]['dates'])>=2 for month in WINTER)]
        full_yes=sum(all(m[k]['yes']>0 for k in WINTER) for _,m in strict)
        all_yes_dates=sum(all(m[k]['yes']>0 and not m[k]['no'] for k in WINTER) for _,m in strict)
        any_no=sum(any(m[k]['no']>0 and not m[k]['yes'] for k in WINTER) for _,m in strict)
        nstates=len({state for state,_ in strict})
        gate='PASS_FULL_FOUR_MONTH_WINTER_OBSERVATION_GATE' if len(strict)>=10 and nstates>=2 else 'FAIL_FULL_WINTER_SAMPLING_GATE'
        results.append({'species':species,'site_winter_years_with_any_checked_month':len(selected),
                        'site_winter_years_all_four_months_observed':len(all4),
                        'site_winter_years_all_four_months_two_or_more_dates_each':len(strict),
                        'states_with_strict_four_month_winter_sampling':nstates,
                        'strict_winter_siteyears_all_four_months_any_yes':full_yes,
                        'strict_winter_siteyears_all_four_months_only_yes_on_sampled_dates':all_yes_dates,
                        'strict_winter_siteyears_at_least_one_month_only_no':any_no,
                        'information_gate':gate})
    # A same-State×winter-year×calendar-month comparison of A. curassavica vs
    # A. syriaca still cannot control site or taxonomic physiology.
    a,b='Asclepias curassavica','Asclepias syriaca'
    matched=[]
    keys={(state,wy,m) for sp,state,wy,m in by_state_month if sp==a}
    keys &= {(state,wy,m) for sp,state,wy,m in by_state_month if sp==b}
    for key in sorted(keys):
        c=by_state_month[(a,*key)];s=by_state_month[(b,*key)]
        if min(c['sampled_site_months'],s['sampled_site_months'])>=2:
            matched.append((key,c,s))
    species_site_months={name:sum(v['sampled_site_months'] for (sp,_,_,_),v in by_state_month.items() if sp==name) for name in (a,b)}
    enough=len(matched)>=10 and len({x[0][0] for x in matched})>=2
    return {'species_fully_monitored_winter_gate':results,
            'strict_repeated_site_months':species_site_months,
            'within_state_winter_month_2018_2025':{
                'candidate_paired_state_winter_month_cells':len(keys),
                'qualifying_cells_min2_repeated_sites_each_species':len(matched),
                'distinct_qualifying_states':len({x[0][0] for x in matched}),
                'curassavica_all_yes_site_months':sum(c['all_yes_sampled_dates'] for _,c,s in matched),
                'curassavica_all_no_site_months':sum(c['all_no_sampled_dates'] for _,c,s in matched),
                'syriaca_all_yes_site_months':sum(s['all_yes_sampled_dates'] for _,c,s in matched),
                'syriaca_all_no_site_months':sum(s['all_no_sampled_dates'] for _,c,s in matched),
                'gate':'PASS_EXPLORATORY_SHARED_STATE_COMPARABILITY' if enough else 'INSUFFICIENT_SHARED_STATE_COMPARABILITY',
                'interpretation':'Descriptive observed site-months, not cultivated-vs-wild causal contrast.'
            }}


def run(registry,output):
    selected,sha=frozen_registry(registry)
    records=[];receipts=[];quality=Counter()
    for year in YEARS:
        data,receipt=get_records(year,selected)
        if receipt['status']!='COMPLETE_SOURCE_HTTP_200':
            raise ValueError('Unknown year missing exact phenophase source')
        parsed,qc=parse_valid(data,year,selected)
        records.extend(parsed);quality.update(qc);receipts.append(receipt)
        print(json.dumps({'year':year,'source_rows':len(data),'accepted':len(parsed)}),flush=True)
    result=complete_winter_gate(site_winter_months(records))
    result.update({'schema':'chocho_usanpn_winter_full_leaf_sampling_gate_v01',
        'status':'COMPLETE_INDEPENDENT_SOURCE_WINTER_COVERAGE_AUDIT',
        'source':'USA National Phenology Network Nature’s Notebook (2018–2025)',
        'primary_pre_outcome_protocol':'USANPN_LIVE_LEAF_OBSERVATION_PROTOCOL_20261028.md',
        'scope':'POSTHOC_FULL_FOUR_MONTH_SOURCE_COVERAGE_SENSITIVITY',
        'source_registry_sha256':sha,'source_year_sha256':[x['response_sha256'] for x in receipts],
        'retained_source_individual_days':len(records),
        'no_raw_plant_site_individual_or_observer_exported':True,
        'no_future_leaf_pheno_response_assumed':True,
        'not_directly_merged_with_monarch_840_hexes':True})
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],
        'full_winter_gate':[(r['species'],r['site_winter_years_all_four_months_two_or_more_dates_each']) for r in result['species_fully_monitored_winter_gate']],
        'matched':result['within_state_winter_month_2018_2025']},ensure_ascii=False),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--registry',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.registry,a.output)
