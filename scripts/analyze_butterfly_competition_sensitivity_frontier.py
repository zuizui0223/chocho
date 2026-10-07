#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,math
from collections import defaultdict
from pathlib import Path
from statistics import median

def load_desc(p):
    out={}
    with open(p,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            out[r['species']]={'b':float(r['host_family_count']),'n':int(r['host_wgsrpd3_unit_count'])}
    return out

def load_pairs(p):
    by=defaultdict(set)
    with open(p,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            s=r['insect_species'].strip(); h=r['accepted_plant_name_id'].strip()
            if s and h: by[s].add(h)
    return by

def load_units(p):
    out=defaultdict(set)
    with open(p,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            h=r['accepted_plant_name_id'].strip(); u=r['area_code_l3'].strip()
            if h and u: out[h].add(u)
    return out

def q(xs,p):
    xs=sorted(xs)
    if not xs:return None
    z=p*(len(xs)-1); i=int(math.floor(z)); j=int(math.ceil(z))
    if i==j:return xs[i]
    a=z-i; return xs[i]*(1-a)+xs[j]*a

def summarize(vals):
    vals=list(vals)
    return {'n':len(vals),'mean':sum(vals)/len(vals) if vals else None,'median':median(vals) if vals else None,
            'q25':q(vals,.25),'q75':q(vals,.75)}

def main():
    ap=argparse.ArgumentParser()
    for x in ['descriptors_csv','insect_host_csv','native_distribution_csv','contemporary_distribution_csv','climate_crossfit_csv']:
        ap.add_argument('--'+x.replace('_','-'),type=Path,required=True)
    ap.add_argument('--output-json',type=Path,required=True)
    ap.add_argument('--output-frontier-csv',type=Path,required=True)
    ap.add_argument('--output-region-csv',type=Path,required=True)
    args=ap.parse_args()

    desc=load_desc(args.descriptors_csv); pairs=load_pairs(args.insect_host_csv)
    native=load_units(args.native_distribution_csv); cont=load_units(args.contemporary_distribution_csv)
    spp=[s for s,d in desc.items() if d['b']>0 and d['n']>0]
    assert len(spp)==239
    focal=set(spp)

    n_union={}; c_union={}; a_union={}
    host_consumers=defaultdict(set)
    for s in spp:
        for h in pairs[s]: host_consumers[h].add(s)
    for s in spp:
        n=set(); c=set()
        for h in pairs[s]:
            n |= native.get(h,set()); c |= cont.get(h,set())
        assert len(n)==desc[s]['n']
        n_union[s]=n; c_union[s]=c; a_union[s]=c-n
    assert sum(map(len,a_union.values()))==14553

    ncomp=defaultdict(set); ccomp=defaultdict(set)
    for h, consumers0 in host_consumers.items():
        consumers=sorted(consumers0 & focal)
        if len(consumers)<2: continue
        for u in native.get(h,set()):
            for s in consumers: ncomp[(s,u)].update(x for x in consumers if x!=s)
        for u in cont.get(h,set()):
            for s in consumers: ccomp[(s,u)].update(x for x in consumers if x!=s)

    native_k=[len(ncomp[(s,u)]) for s in spp for u in n_union[s]]
    cont_k=[len(ccomp[(s,u)]) for s in spp for u in c_union[s]]
    added_k=[len(ccomp[(s,u)]) for s in spp for u in a_union[s]]

    alphas=[0,0.02,0.05,0.1,0.25,0.5,1,2,5]
    comp_sensitivity={}
    for form in ['reciprocal','exponential']:
        arr=[]
        for alpha in alphas:
            def w(k):
                return 1/(1+alpha*k) if form=='reciprocal' else math.exp(-alpha*k)
            en=[]; ec=[]
            for s in spp:
                a=sum(w(len(ncomp[(s,u)])) for u in n_union[s])
                b=sum(w(len(ccomp[(s,u)])) for u in c_union[s])
                en.append(a); ec.append(b)
            agg_n=sum(en); agg_c=sum(ec)
            ratios=[b/a if a>0 else None for a,b in zip(en,ec)]
            valid=[x for x in ratios if x is not None]
            arr.append({'alpha':alpha,'effective_native_opportunity':agg_n,'effective_contemporary_opportunity':agg_c,
                        'aggregate_effective_gain_fraction':agg_c/agg_n-1,
                        'species_with_positive_effective_gain':sum(b>a for a,b in zip(en,ec)),
                        'median_species_effective_ratio':median(valid)})
        comp_sensitivity[form]=arr

    expanded=[s for s in spp if a_union[s]]
    dilution=[len(n_union[s])/len(c_union[s]) for s in expanded]
    intra={'assumption':'If total conspecific abundance is fixed and spreads uniformly across resource regions, per-region conspecific density scales as 1/resource breadth.',
           'expanded_species':len(expanded),'median_contemporary_to_native_per_region_density_ratio':median(dilution),
           'median_density_reduction_fraction':1-median(dilution),'density_ratio_q25':q(dilution,.25),'density_ratio_q75':q(dilution,.75)}

    all_regions=sorted(set().union(*c_union.values()))
    region_rows=[]
    for u in all_regions:
        n=sum(u in n_union[s] for s in spp); a=sum(u in a_union[s] for s in spp); c=sum(u in c_union[s] for s in spp)
        ks=[len(ccomp[(s,u)]) for s in spp if u in c_union[s]]
        region_rows.append({'wgsrpd3_code':u,'native_resource_butterflies':n,'introduced_added_butterflies':a,
                            'contemporary_resource_butterflies':c,'introduced_dependency_share':a/c if c else 0,
                            'mean_exact_host_cousers':sum(ks)/len(ks) if ks else 0})
    region_rows.sort(key=lambda r:(-r['introduced_dependency_share'],-r['introduced_added_butterflies'],r['wgsrpd3_code']))

    climate=[]
    with open(args.climate_crossfit_csv,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f): climate.append(r)
    eval_by=defaultdict(list)
    for r in climate:
        try:m=float(r['climate_mismatch_to_train_niche'])
        except:continue
        if r['observed_unit_split']=='eval_observed' and str(r['effort_supported']) in {'1','True','TRUE','true'}:
            eval_by[r['species']].append(m)
    frontier=[]
    for r in climate:
        s=r['species']; u=r['wgsrpd3_code']
        if s not in focal or u not in a_union.get(s,set()): continue
        if r['observed_unit_split']!='never_observed': continue
        if str(r['effort_supported']) not in {'1','True','TRUE','true'}: continue
        try:m=float(r['climate_mismatch_to_train_niche'])
        except:continue
        if not math.isfinite(m) or not eval_by.get(s): continue
        ev=eval_by[s]
        compat=sum(x>=m for x in ev)/len(ev)
        frontier.append({'species':s,'wgsrpd3_code':u,'wgsrpd3_name':r['wgsrpd3_name'],
                         'climate_mismatch':m,'climate_compatibility_rank':compat,
                         'other_pilot_record_effort':int(float(r['other_pilot_record_effort'])),
                         'eval_observed_units_for_calibration':len(ev)})
    frontier.sort(key=lambda r:(-r['climate_compatibility_rank'],-r['other_pilot_record_effort'],r['species'],r['wgsrpd3_code']))

    payload={'schema':'chocho_butterfly_competition_sensitivity_frontier_v0.1',
      'status':'SUCCESS_POSTHOC_COMPETITION_SENSITIVITY_AND_TESTABLE_FRONTIER',
      'resource_crowding':{'native':summarize(native_k),'contemporary':summarize(cont_k),'introduced_added':summarize(added_k),
                           'native_zero_competitor_fraction':sum(k==0 for k in native_k)/len(native_k),
                           'contemporary_zero_competitor_fraction':sum(k==0 for k in cont_k)/len(cont_k),
                           'introduced_added_zero_competitor_fraction':sum(k==0 for k in added_k)/len(added_k)},
      'competition_discount_sensitivity':{'definition':'Each butterfly×region resource opportunity is down-weighted by the number k of other focal butterflies sharing at least one exact host species in that region. Reciprocal weight=1/(1+alpha*k); exponential weight=exp(-alpha*k). Alpha is an assumed competition-strength sensitivity parameter, not estimated from data.',
                                           'scenarios':comp_sensitivity},
      'intraspecific_dilution_scenario':intra,
      'regional_introduced_resource_dependency':{'regions':len(region_rows),
          'top_dependency_regions_min10_contemporary':[r for r in region_rows if r['contemporary_resource_butterflies']>=10][:30],
          'top_added_count_regions':sorted(region_rows,key=lambda r:(-r['introduced_added_butterflies'],-r['introduced_dependency_share']))[:30]},
      'testable_colonization_frontier':{'definition':'Introduced-only host-resource regions in the 24-species climate panel that remain never observed despite effort support. Climate compatibility is the fraction of effort-supported held-out observed units with mismatch >= the candidate mismatch; larger values mean the candidate climate is at least as close to the training niche as more observed units.',
          'candidate_units':len(frontier),'species':len(set(r['species'] for r in frontier)),
          'compatibility_ge_0_5':sum(r['climate_compatibility_rank']>=.5 for r in frontier),
          'compatibility_ge_0_8':sum(r['climate_compatibility_rank']>=.8 for r in frontier),
          'top_candidates':frontier[:50]},
      'claim_boundary':'Competition scenarios are sensitivity analyses, not estimates of realized competition. The frontier is a falsifiable priority list for future occurrence/field validation, not a forecast of guaranteed colonization; dispersal, host quality, abundance, phenology, enemies and unmeasured habitat can prevent realization.'}
    args.output_json.parent.mkdir(parents=True,exist_ok=True);args.output_json.write_text(json.dumps(payload,indent=2)+'\n')
    with open(args.output_frontier_csv,'w',newline='',encoding='utf-8') as f:
        if frontier:
            w=csv.DictWriter(f,fieldnames=list(frontier[0]));w.writeheader();w.writerows(frontier)
    with open(args.output_region_csv,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(region_rows[0]));w.writeheader();w.writerows(region_rows)
    print(json.dumps(payload,indent=2))
if __name__=='__main__':main()
