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
    out=defaultdict(set)
    with open(p,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            s=r['insect_species'].strip(); h=r['accepted_plant_name_id'].strip()
            if s and h: out[s].add(h)
    return out

def load_units(p):
    out=defaultdict(set)
    with open(p,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            h=r['accepted_plant_name_id'].strip(); u=r['area_code_l3'].strip()
            if h and u: out[h].add(u)
    return out

def weight(k,a,mode):
    if mode=='reciprocal': return 1/(1+a*k)
    if mode=='exponential': return math.exp(-a*k)
    raise ValueError(mode)

def score(counts,a,mode):
    return sum(weight(k,a,mode) for k in counts)

def root(fn,max_alpha=1e6):
    lo=0.0; flo=fn(lo)
    if flo<=0: return 0.0
    hi=1e-6
    fhi=fn(hi)
    while fhi>0 and hi<max_alpha:
        hi*=2
        fhi=fn(hi)
    if fhi>0: return None
    for _ in range(80):
        mid=(lo+hi)/2
        fm=fn(mid)
        if fm>0: lo=mid
        else: hi=mid
    return hi

def strata(b):
    if b<=1:return '1_family'
    if b<=2:return '2_families'
    if b<=5:return '3_to_5_families'
    return '6plus_families'

def main():
    ap=argparse.ArgumentParser()
    for x in ['descriptors_csv','insect_host_csv','native_distribution_csv','contemporary_distribution_csv']:
        ap.add_argument('--'+x.replace('_','-'),type=Path,required=True)
    ap.add_argument('--output-json',type=Path,required=True)
    ap.add_argument('--output-species-csv',type=Path,required=True)
    args=ap.parse_args()
    d=load_desc(args.descriptors_csv); pairs=load_pairs(args.insect_host_csv)
    native=load_units(args.native_distribution_csv); cont=load_units(args.contemporary_distribution_csv)
    spp=[s for s,x in d.items() if x['b']>0 and x['n']>0]; assert len(spp)==239
    consumers=defaultdict(set)
    for s in spp:
        for h in pairs[s]: consumers[h].add(s)
    nc={};cc={};nr={};cr={}
    for s in spp:
        nset=set();cset=set();ncomp=defaultdict(set);ccomp=defaultdict(set)
        for h in pairs[s]:
            oth=consumers[h]-{s}
            for u in native.get(h,set()): nset.add(u);ncomp[u].update(oth)
            for u in cont.get(h,set()): cset.add(u);ccomp[u].update(oth)
        assert len(nset)==d[s]['n']
        nr[s]=nset;cr[s]=cset;nc[s]=[len(ncomp[u]) for u in nset];cc[s]=[len(ccomp[u]) for u in cset]
    expanded=[s for s in spp if len(cr[s])>len(nr[s])];assert len(expanded)==206
    out={'schema':'chocho_butterfly_competition_break_even_v0.1',
         'status':'SUCCESS_HYPOTHETICAL_COMPETITION_BREAK_EVEN_SENSITIVITY',
         'definition':'Break-even alpha is the smallest hypothetical crowding-penalty strength at which competition-discounted contemporary resource opportunity no longer exceeds native opportunity. Alpha is not estimated from observations.'}
    rows=[]
    for mode in ['reciprocal','exponential']:
        def agg(a):
            return sum(score(cc[s],a,mode)-score(nc[s],a,mode) for s in spp)
        ar=root(agg)
        vals=[];inf=0
        for s in expanded:
            rr=root(lambda a,s=s: score(cc[s],a,mode)-score(nc[s],a,mode))
            if rr is None: inf+=1
            else: vals.append(rr)
            rows.append({'species':s,'host_family_count':d[s]['b'],'breadth_class':strata(d[s]['b']),
                         'model':mode,'break_even_alpha':rr,'raw_added_regions':len(cr[s])-len(nr[s]),
                         'native_zero_competitor_regions':sum(k==0 for k in nc[s]),
                         'contemporary_zero_competitor_regions':sum(k==0 for k in cc[s])})
        groups={}
        for g in ['1_family','2_families','3_to_5_families','6plus_families']:
            vv=[r['break_even_alpha'] for r in rows if r['model']==mode and r['breadth_class']==g and r['break_even_alpha'] is not None]
            total=sum(r['model']==mode and r['breadth_class']==g for r in rows)
            groups[g]={'expanded_species':total,'finite_break_even_species':len(vv),
                       'no_break_even_within_search':total-len(vv),
                       'median_finite_break_even_alpha':median(vv) if vv else None}
        out[mode]={
          'aggregate_break_even_alpha':ar,
          'expanded_species':len(expanded),
          'finite_break_even_species':len(vals),
          'no_break_even_within_search':inf,
          'finite_threshold_summary':{
             'median':median(vals) if vals else None,
             'q25':sorted(vals)[int(.25*(len(vals)-1))] if vals else None,
             'q75':sorted(vals)[int(.75*(len(vals)-1))] if vals else None,
             'alpha_le_0_05':sum(v<=.05 for v in vals),
             'alpha_le_0_1':sum(v<=.1 for v in vals),
             'alpha_le_0_25':sum(v<=.25 for v in vals),
             'alpha_le_0_5':sum(v<=.5 for v in vals),
             'alpha_le_1':sum(v<=1 for v in vals),
             'alpha_le_2':sum(v<=2 for v in vals),
          },
          'by_breadth_class':groups
        }
    finite_exp=[r for r in rows if r['model']=='exponential' and r['break_even_alpha'] is not None]
    finite_exp.sort(key=lambda r:(r['break_even_alpha'],r['species']))
    out['most_fragile_exponential']=finite_exp[:25]
    out['most_robust_finite_exponential']=finite_exp[-25:]
    out['claim_boundary']='These thresholds are dimensionless sensitivity break-points under stylized penalties, not measured competition coefficients. They show how strong resource-sharing costs would need to be to erase reconstructed geographic gains; abundance, carrying capacity, enemy sharing and demographic effects remain unobserved.'
    args.output_json.parent.mkdir(parents=True,exist_ok=True);args.output_json.write_text(json.dumps(out,indent=2)+'\n')
    with open(args.output_species_csv,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
